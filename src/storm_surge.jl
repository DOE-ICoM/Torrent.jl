
#     DISCLAIMER
#
# This material was prepared as an account of work sponsored by an agency
# of the United States Government.  Neither the United States Government
# nor the United States Department of Energy, nor Battelle, nor any of
# their employees, nor any jurisdiction or organization that has cooperated
# in the development of these materials, makes any warranty, express or
# implied, or assumes any legal liability or responsibility for the accuracy,
# completeness, or usefulness or any information, apparatus, product, software,
# or process disclosed, or represents that its use would not infringe privately
# owned rights. Reference herein to any specific commercial product, process,
# or service by trade name, trademark, manufacturer, or otherwise does not
# necessarily constitute or imply its endorsement, recommendation, or favoring
# by the United States Government or any agency thereof, or Battelle Memorial
# Institute. The views and opinions of authors expressed herein do not
# necessarily state or reflect those of the United States Government or any
# agency thereof.
#
#    PACIFIC NORTHWEST NATIONAL LABORATORY
#        operated by
#    BATTELLE
#        for the
#    UNITED STATES DEPARTMENT OF ENERGY
#         under Contract DE-AC05-76RL01830

"""
    binarize(g::Grid; level::Float64=0.01, invert::Bool=false)

Thresholds a grid at `level`, replacing all values at or above `level` with
1.0 and all values below `level` with 0.0. By passing the keyword argument
`invert` the resulting mask can be inverted.
"""
function binarize(g::Grid; level::Float64=0.01, invert::Bool=false)
  broadcast(x -> ((invert ? x>=level : x<level) ? 0.0 : 1.0), g)
end


"""
    surge_source_field(
      dem::Grid,
      smoothing_scale::Float64,
      inset_from_smoothed_boundary::Float64,
      width_of_source_line::Float64,
      backstop_distance::Float64
    ) :: Tuple{Grid,Grid,Grid}

Takes a DEM as input and creates a smooth source line a distance
`smoothing_scale - inset_from_smoothed_boundary` off the coast.
Returns three `Grid` instances, each the same size as the DEM.

- The source field, with cells in a source line set to 1.0 and
all other cells set to 0.0.
- The backstop wall mask, some multiple of which can be added to
the DEM to prevent surge-based rivulets from just flowing back
out into the open ocean, but which could be lowered over time to
allow rivulets to flow over the wall and out of the domain.
- And a filled-in backstop mask that can be used for computing
a depth-volume curve.
"""
function surge_source_field(
  dem::Grid,
  smoothing_scale::Float64,
  inset_from_smoothed_boundary::Float64,
  width_of_source_line::Float64,
  backstop_distance::Float64
) :: Tuple{Grid,Grid,Grid}

  # create masks we'll need to generate the source line
  sea_mask = binarize(smooth(dem, smoothing_scale), invert=true)
  backstop = binarize(smooth(sea_mask, inset_from_smoothed_boundary - backstop_distance))
  outset_backstop =  binarize(smooth(sea_mask, inset_from_smoothed_boundary - backstop_distance - width_of_source_line))
  base_source_line = binarize(smooth(sea_mask, inset_from_smoothed_boundary))
  inset_source_line = binarize(smooth(sea_mask, inset_from_smoothed_boundary + width_of_source_line))

  backstop_wall = Grid(dem.registration, similar(backstop.data))

  # then basically do a binary subtraction to create the source line
  for row in 1:dem.registration.nrows
    for col in 1:dem.registration.ncols
      inset_source_line[row,col] = inset_source_line[row,col] > 0.5 && base_source_line[row,col] < 0.5 ? 1.0 : 0.0
      backstop_wall[row,col] = backstop[row,col] > 0.5 && outset_backstop[row,col] < 0.5 ? 1.0 : 0.0
    end
  end

  (inset_source_line, backstop_wall, backstop)
end



"""
    coord_translation(x::Real, y::Real, x0::Real, y0::Real, θ::Real)::Tuple{Float64,Float64}

Transpose to a translated and rotated coordinate system. `x` and `y` are the coordinates of
the point you want to convert. `x0` and `y0` indicate the point of origin of the new reference
frame and `θ` the rotation of the coordinate system about that new point of origin.
"""
function coord_translation(x::Real, y::Real, x0::Real, y0::Real, θ::Real)::Tuple{Float64,Float64}
  x1 = (x-x0) * cos(θ * π/180.0) + (y-y0) * sin(θ * π/180.0)
  y1 = -(x-x0) * sin(θ * π/180.0) + (y-y0) * cos(θ * π/180.0)
  (x1,y1)
end


"""
    coord_translation(index::Index, index0::Index, θ::Real)::Index

Transpose to a translated and rotated coordinate system. `index` are the coordinates of
the point you want to convert. `index0` indicates the point of origin of the new reference
frame, and `θ` the rotation of the coordinate system about that new point of origin.
"""
function coord_translation(index::Index, index0::Index, θ::Real)::Index
  (x1,y1) = coord_translation(index.col, index.row, index0.col, index0.row, θ)
  Index(round(Int,y1), round(Int, x1))
end



"""
    surge_depth_volume_curve(dem::Grid, max_surge_depth::Float64, depth_increment::Float64)

Computes the volume of fluid that would be required to fill the `dem` to different depths
ranging from `depth_increment` to `max_surge_depth` in steps of `depth_increment`.
"""
function surge_depth_volume_curve(dem::Grid, max_surge_depth::Float64, depth_increment::Float64)

  dV = [begin
    # if topographic elevation is lower than surge depth, mark cell as filled
    broadcast(elev -> ( (elev<depth) ? 1.0 : 0.0), dem) |> 
      # count number of filled cells
      (g -> sum(g.data)) |>
      # and compute their total volume
      (x -> x * depth_increment * dem.registration.cell_size_meters^2)
  end for depth in (0.0:depth_increment:max_surge_depth)]

  V = accumulate((+), dV)

  [[depth_increment*(i-1), V[i]] for i in 1:length(V)]

end


"""
    temporal_surge_flux_old(
      backstopped_dem::Grid,
      peak_depth::Float64,
      peak_time_step::Float64,
      surge_period_in_steps::Float64,
      time_step::Float64,
      max_steps::Int
    )

Returns an array of time step versus volume flux that approximates incoming
surge with `peak_depth` center at `peak_time_step` and with a width of
`surge_period_in_steps` approximated as a normal distribution.

NOTE: At the moment it only returns positive fluxes. Negative fluxes are
set to 0.0.
"""
function temporal_surge_flux_old(
  backstopped_dem::Grid,
  peak_depth::Float64,
  peak_time_step::Float64,
  surge_period_in_steps::Float64,
  time_step::Float64,
  max_steps::Int
)::Vector{Vector{Float64}}

  points = surge_depth_volume_curve(backstopped_dem, peak_depth, 0.25)
  f_depth_volume = interpolate_points(points)

  # compute temporal array of volumes necessary to fill DEM proposed level
  # this level is given by a gaussian curve centered at `peak_time_step` and
  # with width `surge_period_in_steps`
  surge_volumes = [ f_depth_volume(peak_depth*exp(-(i-peak_time_step)^2/(2.0*surge_period_in_steps^2))) for i in 1:max_steps]

  # return an array of time step versus requisite flux to achieve proposed volumes
  surge_fluxes = [ [Float64(i), i == 1 ? 0.0 : max((surge_volumes[i]-surge_volumes[i-1])/time_step, 0.0)] for i in 1:length(surge_volumes)]
  surge_fluxes
end


"""
    temporal_surge_flux(
      backstop_wall_mask::Grid,
      surge_line_width::Float64,
      manning_coef::Float64,
      peak_depth::Float64,
      peak_time_step::Float64,
      surge_period::Float64,
      time_step_seconds::Float64,
      max_steps::Int
    )::Vector{Vector{Float64}}

TBW
"""
function temporal_surge_flux(
  backstop_wall_mask::Grid,
  surge_line_width::Float64,
  manning_coef::Float64,
  peak_depth::Float64,
  peak_time_step::Float64,
  surge_period::Float64,
  time_step_seconds::Float64,
  max_steps::Int
)::Vector{Vector{Float64}}

  # note that in the calculation of the length of the surge source line, we
  # divide by the line width plus 2. the plus 2 is for he effects of anti-aliasing
  # when generating the source line which makes the line wider than the proposed
  # value
  cell_size_meters = backstop_wall_mask.registration.cell_size_meters
  length_of_surge_source_line = sum(backstop_wall_mask.data) / (surge_line_width + 2.0)

  # where t is in time steps rather than seconds
  function phi(t::Float64)
    exp_part = peak_depth * exp(-(t-peak_time_step)^2/(2*surge_period^2))
    max(
      (cell_size_meters^2 + cell_size_meters/manning_coef^2 * exp_part^(7/3)) * (peak_time_step-t)/surge_period^2 * exp_part,
      0.0
    ) / time_step_seconds
  end

  surge_fluxes = [ [Float64(i), phi(Float64(i)) * length_of_surge_source_line] for i in 1:max_steps]
  surge_fluxes

end



function temporal_surge_flux_b(
  backstop_wall_mask::Grid,
  surge_line_width::Float64,
  manning_coef::Float64,
  peak_depth::Float64,
  peak_time_step::Float64,
  surge_period::Float64,
  time_step_seconds::Float64,
  max_steps::Int
)::Vector{Vector{Float64}}

  # note that in the calculation of the length of the surge source line, we
  # divide by the line width plus 2. the plus 2 is for he effects of anti-aliasing
  # when generating the source line which makes the line wider than the proposed
  # value
  cell_size_meters = backstop_wall_mask.registration.cell_size_meters
  length_of_surge_source_line = sum(backstop_wall_mask.data) / (surge_line_width + 2.0)

  # where t is in time steps rather than seconds
  function phi(t::Float64)
    exp_part = peak_depth * exp(-(t-peak_time_step)^2/(2*surge_period^2))
    t_part = max( (peak_time_step-t)/surge_period^2, 0)
    max(
      cell_size_meters^2 * t_part * exp_part + 
        cell_size_meters/manning_coef^(2/3) * t_part^(1/3) * exp_part^(16/9),
      0.0
    ) / time_step_seconds
  end

  surge_fluxes = [ [Float64(i), phi(Float64(i)) * length_of_surge_source_line] for i in 1:max_steps]
  surge_fluxes

end




