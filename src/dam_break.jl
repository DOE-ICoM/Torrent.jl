
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
    reservoir_volume(;
      breach_width_top::Float64,
      breach_width_bottom::Float64,
      initial_reservoir_volume::Float64,
      reservoir_depth_curve::Function,
      initial_dam_height::Float64,
      final_dam_height::Float64,
      failure_period::Float64,  # seconds
      time_step::Float64,  # seconds
      num_steps::Int
    )

Produces time traces of the reservoir volume and dam height (returned
as a tuple in that order) based on a dam breach scenario.

# Arguments
- `breach_width_top`: Width of the top of the breach (in meters).
- `breach_width_bottom`: Width of the bottom of the breach (in meters).
- `initial_reservoir_volume`: Initial reservoir volume at start of the breach
  (in cubic meters).
- `reservoir_depth_curve`: A function representing the reservoir depth curve.
  This function should take a reservoir volume (in cubic meters) and return
  a depth (in meters).
- `initial_dam_height`: Dam height at start of simulation
- `final_dam_height`: Dam height at end of breach
- `failure_period`: Time period (seconds) over which breach occurs, resulting
  linear reduction in the dam height.
- `time_step`: Length of a time step in seconds.
- `num_steps`: The number of time steps to take during the simulation.

"""
function reservoir_volume(;
  breach_width_top::Float64,
  breach_width_bottom::Float64,
  initial_reservoir_volume::Float64,
  reservoir_depth_curve::Function,
  initial_dam_height::Float64,
  final_dam_height::Float64,
  failure_period::Float64,  # seconds
  time_step::Float64,  # seconds
  num_steps::Int
)

  # ensure realistic breach geometry
  if breach_width_bottom > breach_width_top
    error("Bottom of breach must be narrower than top.")
  end

  # weir coefficient
  Cw = 1.84  # from https://learn.hydrologystudio.com/studio-express/knowledge-base/weirs/
  
  # breach slope (horizontal over vertical)
  Z = 0.5 * (breach_width_top - breach_width_bottom) / (initial_dam_height - final_dam_height)

  # vectors to store reservoir volume and dam height traces
  V = Vector{Float64}(undef, num_steps)
  dam_height = Vector{Float64}(undef, num_steps)

  # compute the initial reservoir volume
  V[1] = initial_reservoir_volume

  # set the initial dam height as well as computing the linear change in the
  # dam height during the period of failure
  dam_height[1] = initial_dam_height
  dz0 = time_step * (final_dam_height - initial_dam_height)/failure_period

  for i in 2:num_steps

    # Compute the flux:
    # if the reservoir depth is greater than the dam height we need to compute
    # it using the weir equation; if it's less than the dam height, the flux
    # is just zero.
    flux = if reservoir_depth_curve(V[i-1]) > dam_height[i-1]
      # flux through rectangular portion
      (Cw * breach_width_bottom * (reservoir_depth_curve(V[i-1]) - dam_height[i-1])^(1.5)) + 
      # flux through angled edges
      2.0 * Cw * Z * (reservoir_depth_curve(V[i-1]) - dam_height[i-1])^2.5
    else
      0.0
    end

    # and update the reservoir volume
    V[i] = V[i-1] - time_step * flux

    # also update the current dam height
    if i * time_step < failure_period
      dam_height[i] = dam_height[i-1] + dz0
    else
      dam_height[i] = dam_height[i-1]
    end
 
  end
  return (V, dam_height)
end



"""
    reservoir_depth_curve_rectangular(area::Float64)::Function

Returns a reservoir depth curve that approximates the reservoir as
a rectangular volume.
"""
function reservoir_depth_curve(area::Float64)::Function
  (volume) -> volume/area
end


"""
    reservoir_depth_curve(curve_points::Vector{Vector{Float64}})

User should specify a set of volume/depth curves with each row containing
a `[volume, depth]` pair. Returns a function that when supplied with a
volume will return a linearly-interpolated depth based on the provided curve.
"""
function reservoir_depth_curve(curve_points::Vector{Vector{Float64}})

  (volume) -> begin
    upper_idx = findfirst(row -> volume<row[1], curve_points)
    if upper_idx == 1
      curve_points[1][2]
    elseif isnothing(upper_idx)
      curve_points[lastindex(curve_points)][2]
    else
      frac = (volume - curve_points[upper_idx-1][1]) / (curve_points[upper_idx][1] - curve_points[upper_idx-1][1])
      curve_points[upper_idx-1][2] + frac * (curve_points[upper_idx][2] - curve_points[upper_idx-1][2])
  
    end
  end
end


"""
    hydrograph(;
      breach_width_top::Float64,
      breach_width_bottom::Float64,
      initial_reservoir_volume::Float64,
      reservoir_depth_curve::Function,
      initial_dam_height::Float64,
      final_dam_height::Float64,
      failure_period::Float64,  # seconds
      time_step::Float64,  # seconds
      num_steps::Int
    )

Produces time traces of the flux through the breach and dam height (returned
as a tuple in that order) based on a dam breach scenario.

# Arguments
- `breach_width_top`: Width of the top of the breach (in meters).
- `breach_width_bottom`: Width of the bottom of the breach (in meters).
- `initial_reservoir_volume`: Initial reservoir volume at start of the breach
  (in cubic meters).
- `reservoir_depth_curve`: A function representing the reservoir depth curve.
  This function should take a reservoir volume (in cubic meters) and return
  a depth (in meters).
- `initial_dam_height`: Dam height at start of simulation
- `final_dam_height`: Dam height at end of breach
- `failure_period`: Time period (seconds) over which breach occurs, resulting
  linear reduction in the dam height.
- `time_step`: Length of a time step in seconds.
- `num_steps`: The number of time steps to take during the simulation.
"""
function hydrograph(;
  breach_width_top::Float64,
  breach_width_bottom::Float64,
  initial_reservoir_volume::Float64,
  reservoir_depth_curve::Function,
  initial_dam_height::Float64,
  final_dam_height::Float64,
  failure_period::Float64,  # seconds
  time_step::Float64,  # seconds
  num_steps::Int
)

  # compute the time trace of the reservoir volume
  (V, dam_height) = reservoir_volume(
    breach_width_top = breach_width_top,
    breach_width_bottom = breach_width_bottom,
    initial_reservoir_volume = initial_reservoir_volume,
    reservoir_depth_curve = reservoir_depth_curve,
    initial_dam_height = initial_dam_height,
    final_dam_height = final_dam_height,
    failure_period = failure_period,
    time_step = time_step,
    num_steps = num_steps
  )

  # then differentiate to get the flux
  ([-(V[i]-V[i-1])/time_step for i in 2:length(V)], dam_height)  # negative sign because we're treating outflow as positive flux

end

