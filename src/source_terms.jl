
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


"Conversion factor for input precipitation fluxes."
mmPerHourToMetersPerSecond = Float32(1.0 / 1000.0 / 3600.0)

"""
A precipitation event. Each event is valid for a certain period of time
defined by the function, is_valid_when.
"""
struct PrecipitationPeriod

  "The spatial distribution of source locations representing this precipitation period."
  sources::Array{Index}

  "The total volumetric flux (m^3/s) of rainfall during this period."
  total_flux::Float32

  """
  This function should take the current time step as argument and
  return a boolean value reflecting whether the precipitation event
  is valid at that time.
  """
  is_valid_when::Function

end


"A set of precipitation periods ordered in time."
struct PrecipitationTimeSeries
  periods::DataStructures.Queue{PrecipitationPeriod}
end


"""
    get_current_period(pts::PrecipitationTimeSeries, time_step::Int)::PrecipitationPeriod

    Returns the precipitation period corresponding to the specified
    simulation time step.
"""
function get_current_period(pts::PrecipitationTimeSeries, time_step::Int)::PrecipitationPeriod

  # check to see whether the first period on the queue is valid
  if first(pts.periods).is_valid_when(time_step)
    first(pts.periods)

  # otherwise, move through the queue until we find the first valid
  # precipitation period. precipitation periods are assumed to be
  # ordered in the queue and comprehensive, so generally this should
  # just step forward to the next period on the queue.
  else
    while !first(pts.periods).is_valid_when(time_step)
      DataStructures.dequeue!(pts.periods)
    end
  end

  # the applicable precipitation period should now be the first
  first(pts.periods)

end


"""
    create_period(
      num_sources::Int32,
      precipitation_grid::Vector{Vector{Float32}},
      registration::GeoRegistration,
      when_to_use::Function
    )

    Create a new instance of a precipitation period.
"""
function create_period(
  num_sources::Int,
  precipitation_field::Array{Float32, 2},
  from_registration::GeoRegistration,
  to_registration::GeoRegistration,
  is_horizontal::Bool,
  post_filter_source_locations::Bool,
  is_valid_when::Function
)

  # create a spatially distributed set of source locations whose spatial
  # density is correlate with the rainfall rate.
  (source_locations, total_flux) = if post_filter_source_locations
    stochastic_source_distribution_c(
      precipitation_field, num_sources, from_registration, to_registration, is_horizontal
    )
  else
    stochastic_source_distribution_c(
      precipitation_field, num_sources, from_registration, to_registration, is_horizontal
    )
  end


  # create a new precipitation period instance from these source locations
  PrecipitationPeriod(source_locations, total_flux, is_valid_when)

end


"""
    stochastic_source_distribution_old(
      precipitation_field::Array{Float32, 2},
      num_sources::Int,
      from_registration::GeoRegistration,
      to_registration::GeoRegistration
    ) 

    Efficiently creates a stochastic spatial distribution of num_sources
    sources based on the provided precipitation field. The spatial density
    of the source locations is correlated with the rainfall rate distribution.
    The returned indices are in the reference frame of to_registration.
    This is to allow a precipitation field to be used to generate a source
    distribution for a use with a DEM that may have different cell size and/or
    potentially slightly different bounding area (the precipitation field might
    be coarser than the DEM, for example).
"""
function stochastic_source_distribution_a(
  precipitation_field::Array{Float32, 2},
  num_sources::Int,
  from_registration::GeoRegistration,
  to_registration::GeoRegistration
)

  # grab the width of the precipitation raster
  width = size(precipitation_field)[2]

  # Create a normalized, flattened, cumulative distribution of
  # relative rainfall rates
  dist = normalized_cumulative(precipitation_field)

  # along with an ordered sequence of random numbers
  rs = sort(rand(Float32, num_sources))

  # we now determine at what index in the normalized cumulative
  # rainfall distribution each of our random numbers falls. these
  # indices will be distributed in proportion to the rainfall rate.
  indices = DataStructures.MutableLinkedList{Int}()
  rs_idx = 1
  dist_idx = 1
  while rs_idx < num_sources
    while dist_idx < length(dist) && dist[dist_idx] < rs[rs_idx]
      dist_idx += 1
    end
    append!(indices, dist_idx)
    rs_idx += 1
  end

  # next we need to turn these linear indices into 2D grid locations.
  # additionally, we want to randomly shuffle them so that we can simply
  # step through them in order, but have new source locations generated
  # pseudo-randomly with likelihood proportional to the local rainfall.
  # this prevents us from having to generate a random number many times
  # within the simulation, which is relatively costly.
  indices_2d = Random.shuffle!(collect(Index, map(to_table_indices(width), indices)))

  # finally, let's compute the total precipitation flux across the field
  total_flux = sum(precipitation_field) * from_registration.cell_size_meters^2

  # convert indices to locations in target georegistration
  target_indices = convert_source_locations(indices_2d, from_registration, to_registration, true)

  # we return the set of spatial source locations and the total flux
  # for this period.
  (target_indices, total_flux)
end



"""
    stochastic_source_distribution(
      precipitation_field::Array{Float32, 2},
      num_sources::Int,
      from_registration::GeoRegistration,
      to_registration::GeoRegistration
    )

    The old method creates num_sources across the whole from_registration area,
    then filters these down to only those within to_registration. In some cases,
    this can be a tiny fraction (or none) of num_sources. 
"""
function stochastic_source_distribution_b(
  flux_field::Array{Float32, 2},
  num_sources::Int,
  from_registration::GeoRegistration,
  to_registration::GeoRegistration,
  is_horizontal::Bool = false
)

  # rerasterize the precipitation field to the same grid as to_registration
  to_flux_field = rerasterize(flux_field, from_registration, to_registration)

  # grab the width of the precipitation raster
  width = size(to_flux_field)[2]

  # Create a normalized, flattened, cumulative distribution of
  # relative rainfall rates
  dist = normalized_cumulative(to_flux_field)

  # along with an ordered sequence of random numbers
  rs = sort(rand(Float32, num_sources))

  # we now determine at what index in the normalized cumulative
  # rainfall distribution each of our random numbers falls. these
  # indices will be distributed in proportion to the rainfall rate.
  indices = DataStructures.MutableLinkedList{Int}()
  rs_idx = 1
  dist_idx = 1
  while rs_idx < num_sources
    while dist_idx < length(dist) && dist[dist_idx] < rs[rs_idx]
      dist_idx += 1
    end
    append!(indices, dist_idx)
    rs_idx += 1
  end

  # next we need to turn these linear indices into 2D grid locations.
  # additionally, we want to randomly shuffle them so that we can simply
  # step through them in order, but have new source locations generated
  # pseudo-randomly with likelihood proportional to the local rainfall.
  # this prevents us from having to generate a random number many times
  # within the simulation, which is relatively costly.
  indices_2d = Random.shuffle!(collect(Index, map(to_table_indices(width), indices)))

  # finally, let's compute the total precipitation flux across the field
  scale_factor = is_horizontal ? to_registration.cell_size_meters : to_registration.cell_size_meters^2
  total_flux = sum(to_flux_field) * scale_factor

  # we return the set of spatial source locations and the total flux
  # for this period.
  (indices_2d, total_flux)
end



"""
  stochastic_source_distribution_c(
    precipitation_field::Array{Float32, 2},
    num_sources::Int,
    from_registration::GeoRegistration,
    to_registration::GeoRegistration
  )

  TBW
"""
function stochastic_source_distribution_c(
  flux_field::Array{Float32, 2},
  num_sources::Int,
  from_registration::GeoRegistration,
  to_registration::GeoRegistration,
  is_horizontal::Bool = false
)

  # figure out the indices of the upper left and lower right of the possibly smaller,
  # to_registration grid in the coordinates of the larger, from_registration
  ul_bound, lr_bound = convert_source_locations(
    [Index(1,1), Index(to_registration.nrows, to_registration.ncols)],
    to_registration,
    from_registration,
    false;
    filter_in_bounds = false
  )

  upper_left = Index(ul_bound.row < 1 ? 1 : ul_bound.row, ul_bound.col < 1 ? 1 : ul_bound.col)
  lower_right = Index(
    lr_bound.row > from_registration.nrows ? from_registration.nrows : lr_bound.row,
    lr_bound.col > from_registration.ncols ? from_registration.ncols : lr_bound.col
  )

  # then crop the original flux field down to these bounds
  cropped_flux_field = flux_field[upper_left.row:lower_right.row, upper_left.col:lower_right.col]

  # and create a new registration for it since we'll need this step to convert indices at the end
  cropped_registration = GeoRegistration(
    size(cropped_flux_field)[2],
    size(cropped_flux_field)[1],
    from_registration.xll + upper_left.col * from_registration.cell_size,
    from_registration.yll + (from_registration.nrows - lower_right.row) * from_registration.cell_size,
    from_registration.cell_size,
    from_registration.no_data_value,
    from_registration.cell_size_meters,
    from_registration.proj_string,
    from_registration.units
  )

  # grab the width of the precipitation raster
  width = size(cropped_flux_field, 2)

  # Create a normalized, flattened, cumulative distribution of
  # relative rainfall rates
  dist :: Union{Array{Float32,1},Nothing} = nothing
  try
    dist = normalized_cumulative(cropped_flux_field)
  catch ex
    printstyled("ERROR: You might check whether the DEM and precipitation\n"; color=198)
    printstyled("data are in the same CRS (both should be in 4269), that\n"; color=198)
    printstyled("the precipitation grid resolution is symmetric, and that\n"; color=198)
    printstyled("the bounding area, number of rows and cols, and cell size\n"; color=198)
    printstyled("seem self consistent. It may help to resample the precipitation\n"; color=198)
    printstyled("grid with similar extent to the DEM to avoid any challenges that\n"; color=198)
    printstyled("may result from the poor approximation of a large grid as a\n"; color=198)
    printstyled("simple rectangle.\n"; color=198)
    rethrow()
  end

  # along with an ordered sequence of random numbers
  rs = sort(rand(Float32, num_sources))

  # we now determine at what index in the normalized cumulative
  # rainfall distribution each of our random numbers falls. these
  # indices will be distributed in proportion to the rainfall rate.
  indices = DataStructures.MutableLinkedList{Int}()
  rs_idx = 1
  dist_idx = 1
  while rs_idx < num_sources
    while dist_idx < length(dist) && dist[dist_idx] < rs[rs_idx]
      dist_idx += 1
    end
    append!(indices, dist_idx)
    rs_idx += 1
  end

  # next we need to turn these linear indices into 2D grid locations.
  # additionally, we want to randomly shuffle them so that we can simply
  # step through them in order, but have new source locations generated
  # pseudo-randomly with likelihood proportional to the local rainfall.
  # this prevents us from having to generate a random number many times
  # within the simulation, which is relatively costly.
  indices_2d = Random.shuffle!(collect(Index, map(to_table_indices(width), indices)))

  # finally, let's compute the total precipitation flux across the field
  cropped_area = cropped_registration.nrows * cropped_registration.ncols * cropped_registration.cell_size_meters^2
  to_area = to_registration.nrows * to_registration.ncols * to_registration.cell_size_meters^2
  relative_area = to_area/cropped_area
  # debug: print(" relative area=$(relative_area)")
  scale_factor = is_horizontal ? cropped_registration.cell_size_meters : cropped_registration.cell_size_meters^2 * relative_area
  total_flux = sum(cropped_flux_field) * scale_factor

  # convert indices to locations in target georegistration
  target_indices = convert_source_locations(indices_2d, cropped_registration, to_registration, true)

  # we return the set of spatial source locations and the total flux
  # for this period.
  (target_indices, total_flux)
end



"""
    draw_indices_from_weighted_array(array_of_weights::Array{T}, num::Int) where {T <: Real}

    TBW
"""
function draw_indices_from_weighted_array(array_of_weights::Array{T}, num::Int) where {T <: Real}

  # Create a normalized, flattened, cumulative distribution of
  # relative rainfall rates
  dist = normalized_cumulative(array_of_weights)

  # along with an ordered sequence of random numbers
  rs = sort(rand(Float32, num))

  # we now determine at what index in the normalized cumulative
  # rainfall distribution each of our random numbers falls. these
  # indices will be distributed in proportion to the rainfall rate.
  indices = DataStructures.MutableLinkedList{Int}()
  rs_idx = 1
  dist_idx = 1
  while rs_idx < num
    while dist_idx < length(dist) && dist[dist_idx] < rs[rs_idx]
      dist_idx += 1
    end
    append!(indices, dist_idx)
    rs_idx += 1
  end

  Random.shuffle!(collect(indices))
end


"""
    to_table_indices(width: Int)::Function

    Poor man's (i.e. Julia developer's) curried function to compute
    2D indices from linear location in an array, the calculation of
    which depends on the width of the implicit raster. Returns a
    function which takes a linear index as argument and returns
    a 2D spatial location on a grid.
"""
function to_table_indices(width::Int)::Function
  (idx) -> Index( div(idx-1,width)+1, mod(idx-1, width)+1)  # plus ones for Julia's one-offset arrays
end


"""
    convert_source_locations(indices::Array{Index}, from_reg::GeoRegistration, to_reg::GeoRegistration)

    Converts the set of indices from one projection to another. Only the points that lie within the
    boundary of the second georegistered area are returned.
"""
function convert_source_locations(
  indices::Array{Index},
  from_reg::GeoRegistration,
  to_reg::GeoRegistration,
  do_jiggle::Bool = false;
  filter_in_bounds = true )::Array{Index}

  # convert the indices to real-world lat/longs using from_reg
  latlongs = map(idx -> latlongof(idx, from_reg), indices)

  # may not want, for example, precipitation sources to all end up at cell
  # centroids if precipitation grid is more coarse than DEM. let's jiggle
  # them uniformly within the cell size.
  dx = from_reg.cell_size / 2.0
  latlongs = 
    if do_jiggle
      map(
        ll -> [jiggle(ll[1], dx), jiggle(ll[2], dx)], 
        latlongs
      )
    else
      latlongs
    end

  # then convert back to indices using the to_reg
  to_indices = map(ll -> indexof(ll, to_reg), latlongs)

  # keep only those points that lie within the boundary of the to_reg if requested
  if filter_in_bounds
    filter(idx -> !is_outside_boundary(idx, to_reg), to_indices)
  else
    to_indices
  end

end


"""
    jiggle(x::Real, dx::Real)

    Adds a uniformly distributed value within the range of +/- dx to the value of x
"""
function jiggle(x::Real, dx::Real)
  x + 2.0 * dx * (rand() - 0.5)
end


"""
    open_precipitation_series(
      filename_pattern::String,
      start_index::Int32,
      end_index::Int32,
      num_sources::Int32,
      time_step::Float64,
      interval::Float64,
      write_distributions::Bool,
      output_directory::String,
      to_registration::GeoRegistration
    )

    Return a precipitation time series from a set of precipitation
    rate rasters.
"""
function open_precipitation_series(
  filename_pattern::String,
  start_index::Int,
  end_index::Int,
  num_sources::Int,
  time_step::Float64,
  interval::Float64,
  to_registration::GeoRegistration
) :: PrecipitationTimeSeries

  # initialize height, width, and GeoRegistration variables that we
  # can fill in later
  height = 0
  width = 0
  from_registration = GeoRegistration()

  # create a queue to hold new precipitation periods
  events = DataStructures.Queue{PrecipitationPeriod}()

  # create a progress bar to chart precipitation load and generation
  prog = Progress(end_index-start_index; desc="Precipitation       ")

  # precipitation rate raster files are expected to be numbered
  # sequentially so we just step through them in order here.
  for i in start_index:end_index

    # generate a filename based on the provided pattern and index
    filename = Printf.format(Printf.Format(filename_pattern), i)

    # and open the grid file
    precip = open_raster(filename)

    # fill in variables initialized above from the registration information
    from_registration = precip.registration
    height = from_registration.nrows
    width = from_registration.ncols

    # Create a new spatial distribution of source locations corresponding to the
    # spatial rainfall rate.
    (_, evt) = time_computation("") do 

      # create a precipitation period corresponding to this source distribution
      # the third argument is the function that defines when (in terms of
      # simulation time steps) this source distribution is valid. note that
      # we're only really specifying an upper bound here since it's assumed
      # that the expiration of the prior distribution will set the lower bound.
      p = create_period(num_sources,
        precip.data .* mmPerHourToMetersPerSecond,
        from_registration,
        to_registration,
        false,
        true,
        (j) -> (j < ((i-start_index+1) * trunc(Int, interval/time_step)))
      )

      # debug: print("sources=$(length(p.sources)), flux=$(p.total_flux) ")
      update!(prog, i-start_index; showvalues = [("file", short_filename(filename)), ("sources", length(p.sources)), ("flux", "$(p.total_flux) m^3/s")])

      p

    end

    # add the new precipitation period event to the queue
    DataStructures.enqueue!(events, evt)

  end

  # finish progress bar
  finish!(prog)

  # we'll also want to create a default event that can be used if the simulation
  # time extends beyond the precipitation periods we have defined. currently the
  # default is to have no precipitation.
  default_event = create_period(num_sources, zeros(Float32, height, width), from_registration, to_registration, false, true, (j)->true)
  DataStructures.enqueue!(events, default_event)
  
  # once the default event is added, the set of precipitation events is wrapped
  # and returned as a PrecipitationTimeSeries.
  PrecipitationTimeSeries(events)

end


"""
    open_multiband_precipitation(
      filename::String,
      min_index::Int,
      max_index::Int,
      num_sources::Int,
      time_step::Float64,
      interval::Float64,
      scale_factor::Float64,
      write_distributions::Bool,
      output_directory::String,
      to_registration::GeoRegistration
    ) :: PrecipitationTimeSeries

TBW
"""
function open_multiband_precipitation(
  filename::String,
  min_index::Int,
  max_index::Int,
  num_sources::Int,
  time_step::Float64,
  interval::Float64,
  scale_factor::Float64,
  to_registration::GeoRegistration
) :: PrecipitationTimeSeries

  # initialize height, width, and GeoRegistration variables that we
  # can fill in later
  height = 0
  width = 0
  from_registration = GeoRegistration()

  # create a queue to hold new precipitation periods
  events = DataStructures.Queue{PrecipitationPeriod}()

  # read precipitation grids from multiband geotiff
  iter = open_multiband_stream(filename)

  # precipitation rate raster files are expected to be numbered
  # sequentially so we just step through them in order here.
  lower = if min_index>0 min_index else 1 end
  upper = if (max_index>0 && max_index<=length(iter)) max_index else length(iter) end

  # create a progress bar to chart precipitation load and generation
  prog = Progress(upper-lower; desc="Multi-band Precip   ")

  for i in lower:upper

    # grab the grid of interest (we're not treating our poor iterator
    # well here, but should be fine the way it's defined)
    (precip, _) = iterate(iter, i)

    # fill in variables initialized above from the registration information
    from_registration = precip.registration
    height = from_registration.nrows
    width = from_registration.ncols

    # Create a new spatial distribution of source locations corresponding to the
    # spatial rainfall rate.
    (_, evt) = time_computation("Creating source distribution for grid $(i)...") do 

      # let's handle the potential for the precipitation data to include a no data value
      precip_field = broadcast( x -> begin
        isapprox(x, precip.registration.no_data_value) || isnan(x) ? 0f0 : x
      end, precip.data)

      # and let's apply the scale factor (this defaults to 1.0 if none is supplied)
      precip_field = broadcast(x -> x * Float32(scale_factor), precip_field)

      # create a precipitation period corresponding to this source distribution
      # the third argument is the function that defines when (in terms of
      # simulation time steps) this source distribution is valid. note that
      # we're only really specifying an upper bound here since it's assumed
      # that the expiration of the prior distribution will set the lower bound.
      p = create_period(num_sources,
        precip_field .* mmPerHourToMetersPerSecond,
        from_registration,
        to_registration,
        false,
        true,
        (j) -> (j < ((i-lower+1) * trunc(Int, interval/time_step)))
      )

      p
    end

    # update progress bar status.
    update!(prog, i-lower; showvalues = [("band", i), ("sources", length(evt.sources)), ("flux", evt.total_flux)])

    # add the new precipitation period event to the queue
    DataStructures.enqueue!(events, evt)

  end

  # finish progress bar
  finish!(prog)

  # we'll also want to create a default event that can be used if the simulation
  # time extends beyond the precipitation periods we have defined. currently the
  # default is to have no precipitation.
  default_event = create_period(num_sources, zeros(Float32, height, width), from_registration, to_registration, false, true, (j)->true)
  DataStructures.enqueue!(events, default_event)
  
  # once the default event is added, the set of precipitation events is wrapped
  # and returned as a PrecipitationTimeSeries.
  PrecipitationTimeSeries(events)

end



"""
    normalized_cumulative(xs::Array{T}) where {T <: Real}

    Shortcut to create a normalized cumulative distribution from
    an array of xs. The xs can be in the form of either a 1D or 2D array.
"""
function normalized_cumulative(xs::Array{T}) where {T <: Real}
  cumulative(normalize(xs))
end


"""
    normalize(xs::Array{T, 1}) where {T <: Real}

    Normalize xs so that they sum to 1.0
"""
function normalize(xs::Array{T, 1}) where {T <: Real}
  n = length(xs)
  if n == 0
    Array{T}(undef, 0)
  elseif n == 1
    [one(T)]
  else
    total = sum(xs)
    xs ./ total
  end
end


"""
    normalize(data::Array{T,2}) where {T <: Real}

    Returns a flattened, 1D, normalized version of the supplied 2D array, data.
"""
function normalize(data::Array{T,2})::Array{T,1} where {T <: Real}
  data_size = size(data)
  width = data_size[2]
  n = data_size[1] * data_size[2]
  flat = Array{T}(undef, n)
  if n == 0
    Array{T}(undef, 0)
  elseif n == 1
    [one(T)]
  else
    # min_val = data[1,1]
    for i in 1:data_size[2]
      for j in 1:data_size[1]
        x = data[j,i]
        flat[(j-1)*width + i] = x
        # min_val = x < min_val ? x : min_val
      end
    end
    # offset = flat .- min_val
    total = sum(flat)
    flat ./ total
  end

end


"""
    cumulative(xs::Array{T}) where {T <: Real}

    Creates a cumulative tally of the elements of xs.
"""
function cumulative(xs::Array{T}) where {T <: Real}
  cumul = Array{T}(undef, length(xs))
  cumul[1] = xs[1]
  for i in 2:length(xs)
    cumul[i] = cumul[i-1] + xs[i]
  end
  cumul
end


"""
    open_nwm_data(filename::String, time_step::Float64)::PrecipitationTimeSeries

    Opens data from the national water model and returns a PrecipitationTimeSeries
    instance. NWM data is expected to be held in a CSV file, the first column of
    which is the latitude and the second the longitude. In each row then, the
    remaining columns will represent a time series of the lateral influx of
    precipitation at that location in m^3/s. For example, then, if you have 100 time
    steps you would have 2 + 100 columns.
"""
function open_nwm_data(
  nwm_filename::String,
  dem_registration::GeoRegistration,
  time_step::Float64,
  interval::Float64,
  num_sources::Int) :: PrecipitationTimeSeries

  # open data from the national water model that has been recast
  # to the format described above 
  table = open_csv(Float64, nwm_filename)

  # filter to only those points that lie within the DEM grid
  latlongs = slicematrix(table[:,1:2])
  all_grid_locations = indexof.(latlongs, Ref(dem_registration))
  keepers = .!(is_outside_boundary.(all_grid_locations, Ref(dem_registration)))
  bounded_table = table[keepers, :]
  bounded_grid_locations = all_grid_locations[keepers]

  # step through each time period specified in the NWM data
  events = DataStructures.Queue{PrecipitationPeriod}()
  for t in 3:size(bounded_table)[2]

    # compute the total flux during that period of time
    total_flux = sum(bounded_table[:,t]) 

    # choose a random set of source locations from those in the NWM data,
    # with each source locations likelihood of being chosen proportional
    # to the lateral flux there
    indices = draw_indices_from_weighted_array(bounded_table[:,t], num_sources)
    source_locations = bounded_grid_locations[indices]

    # create a new precipitation period from this data
    evt = PrecipitationPeriod(
      source_locations,
      total_flux,
      (j) -> (j < ((t-2) * trunc(Int, interval/time_step)))
    )

    # and add that period to the precipitation event queue
    DataStructures.enqueue!(events, evt)
  end

  # we'll also want to create a default event that can be used if the simulation
  # time extends beyond the precipitation periods we have defined. currently the
  # default is to have no precipitation.
  default_event = PrecipitationPeriod([], 0.0, (j)->true)
  DataStructures.enqueue!(events, default_event)
  
  # once the default event is added, the set of precipitation events is wrapped
  # and returned as a PrecipitationTimeSeries.
  PrecipitationTimeSeries(events)

end



"""
    generate_dam_failure(
      config::Dict{String,Any},
      registration::GeoRegistration,
      iteration::Int) :: Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

TBW
"""
function generate_dam_failure(
  config::Dict{String,Any},
  registration::GeoRegistration,
  iteration::Int) :: Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

  if haskey(config, "dam-failure")
    time_computation("Generating dam failure hydrograph...\n ") do 

      # parse configuration for dam failure parameters

      latitude = config["dam-failure"]["latitude"]
      longitude = config["dam-failure"]["longitude"]
      index = indexof(latitude, longitude, registration)
      if is_outside_boundary(index, registration)
        error("Dam failure location specified lies outside of DEM boundary.")
      end

      breach_width_top = rand_value(parse_distribution(config["dam-failure"]["breach-width-top"]))
      breach_width_bottom = rand_value(parse_distribution(config["dam-failure"]["breach-width-bottom"]))
      reservoir_volume_initial = rand_value(parse_distribution(config["dam-failure"]["reservoir-volume-initial"]))
      reservoir_depth_curve_info = config["dam-failure"]["reservoir-depth-curve"]
      dam_height_initial = rand_value(parse_distribution(config["dam-failure"]["dam-height-initial"]))
      dam_height_final = rand_value(parse_distribution(config["dam-failure"]["dam-height-final"]))
      failure_period = rand_value(parse_distribution(config["dam-failure"]["failure-period"]))
      time_step = config["time-step-seconds"]
      num_steps = config["max-steps"]

      # parse the reservoir depth curve information

      reservoir_curve_function = parse_reservoir_depth_curve(reservoir_depth_curve_info)

      # generate the hydrograph, which will necessarily have one entry
      # for each simulation time step

      (fluxes, _) = hydrograph(
        breach_width_top = breach_width_top,
        breach_width_bottom = breach_width_bottom,
        initial_reservoir_volume = reservoir_volume_initial,
        reservoir_depth_curve = reservoir_curve_function,
        initial_dam_height = dam_height_initial,
        final_dam_height = dam_height_final,
        failure_period = failure_period,
        time_step = time_step,
        num_steps = num_steps
      )

      # save fluxes to a CSV file in the output directory for reference
      flux_table = [[i*time_step, fluxes[i]] for i in 1:length(fluxes)]
      padded_iteration = Printf.@sprintf("%03d", iteration)
      save_csv(config["output-directory"] * "dam-breach-flux-$padded_iteration.csv", flux_table)

      # create precipitation periods for each step
      events = DataStructures.Queue{PrecipitationPeriod}()
      for i in 1:length(fluxes)
        evt = PrecipitationPeriod(
          [index],
          fluxes[i],
          (j) -> (j<i+1)
        )
        DataStructures.enqueue!(events, evt)
      end
      PrecipitationTimeSeries(events)
    end
  else
    (0.0, nothing)
  end

end


"""
    parse_reservoir_depth_curve(curve_info)

Parses the depth `curve_info` based on type and returns
a function that maps reservoir volume to depth.  
"""
function parse_reservoir_depth_curve(curve_info)

  if typeof(curve_info) == String
    table = open_csv(Float64, curve_info)
    curve = [table[i,:] for i in 1:size(table,1)]
    reservoir_depth_curve(curve)
  elseif typeof(curve_info) == Float64
    reservoir_depth_curve(curve_info)
  elseif typeof(curve_info) <: Vector
    # the json parser treats the array of arrays as a Vector{Any}
    # the following step is needed to explicitly convert to A
    # Vector{Vector{Float64}}
    vec = [Float64.(row) for row in curve_info]
    reservoir_depth_curve(vec)
  else
    error("Unexpected type parsing 'reservoir-depth-curve', found $(typeof(curve_info))")
  end

end


"""
    open_bounding_flood_series(
      filename_pattern::String,
      start_index::Int,
      end_index::Int,
      index_step::Int,
      index_offset::Int,
      seconds_per_step::Float64,
      time_step_seconds::Float64,
      surrounding_dem_filename::String,
      surrounding_dem_units::String,
      num_sources::Int,
      inset_from_border::Int,
      write_distributions::Bool,
      output_directory::String,
      inner_registration::GeoRegistration,
      manning_coef::Float64
    ) :: PrecipitationTimeSeries

TBW
"""
function open_bounding_flood_series(
  filename_pattern::String,
  start_index::Int,
  end_index::Int,
  index_step::Int,
  index_offset::Int,
  seconds_per_step::Float64,
  time_step_seconds::Float64,
  surrounding_dem_filename::String,
  surrounding_dem_units::String,
  num_sources::Int,
  inset_from_border::Int,
  inner_registration::GeoRegistration,
  manning_coef::Float64
) :: PrecipitationTimeSeries

  # create a queue to hold new precipitation periods
  events = DataStructures.Queue{PrecipitationPeriod}()

  # open surrounding DEM
  surrounding_dem = open_raster(surrounding_dem_filename, false; units=surrounding_dem_units); println()

  # in most cases we won't have written an initial grid of zero flood depths
  # so we need to add an initial default zero depth event to be used until
  # the first flood boundary condition become applicable
  default_event = create_period(
    num_sources,
    zeros(Float32, inner_registration.nrows, inner_registration.ncols),
    inner_registration,
    inner_registration,
    true,
    false,
    (j)->(j < (start_index-index_offset) * Int(seconds_per_step/time_step_seconds))
  )
  DataStructures.enqueue!(events, default_event)

  # create a progress bar to chart precipitation load and generation
  prog = Progress(end_index-start_index; desc="Boundary Conditions ")

  # step through the flood depth files in order
  for i in start_index:index_step:end_index

    # generate a filename based on the provided pattern and index
    filename = Printf.format(Printf.Format(filename_pattern), i)

    # update progress bar status.
    update!(prog, i-start_index; showvalues = [("file", short_filename(filename))])

    # and open the grid file
    surrounding_flood = open_raster(filename)

    # Create a new spatial distribution of source locations on the boundary
    # corresponding to the estimated incoming flux.
    (_, evt) = time_computation("Creating boundary distribution...") do

      flux_grid = create_boundary_flux_grid(surrounding_flood, surrounding_dem, manning_coef, inner_registration, inset_from_border)

      p = create_period(num_sources,
        flux_grid,
        inner_registration,  # the flux should already be in the registration of the inner, hi-res DEM
        inner_registration,
        true,
        false,
        (j) -> (j < (i+index_step-index_offset) * Int(seconds_per_step/time_step_seconds))
      )

      # print(" sources=$(length(p.sources)), flux=$(p.total_flux), valid till step=$((i+index_step-index_offset) * Int(seconds_per_step/time_step_seconds))")
      p
      
    end

    # add the new precipitation period event to the queue
    DataStructures.enqueue!(events, evt)

  end

  # finish progress bar
  finish!(prog)

  # we'll also want to create a default event that can be used if the simulation
  # time extends beyond the precipitation periods we have defined. currently the
  # default is to have no precipitation.
  default_event = create_period(
    num_sources,
    zeros(Float32, inner_registration.nrows, inner_registration.ncols),
    inner_registration,
    inner_registration,
    true,
    false,
    (j)->true
  )

  DataStructures.enqueue!(events, default_event)
  
  # once the default event is added, the set of precipitation events is wrapped
  # and returned as a PrecipitationTimeSeries.
  PrecipitationTimeSeries(events)

end



"""
    create_boundary_flux_grid(
      surrounding_flood :: Grid,
      surrounding_dem :: Grid,
      manning_coef :: Float64,
      inner_registration :: GeoRegistration,
      inset_from_border :: Int = 3
    ) :: Array{Float32, 2}

TBW
"""
function create_boundary_flux_grid(
  surrounding_flood :: Grid,
  surrounding_dem :: Grid,
  manning_coef :: Float64,
  inner_registration :: GeoRegistration,
  inset_from_border :: Int = 3
) :: Array{Float32, 2}

  offset = inset_from_border + 1

  # initialize a grid of zeros to hold estimated boundary fluxes
  flux_field = zeros(Float32, inner_registration.nrows, inner_registration.ncols)

  # top boundary
  inner_indices = [Index(offset,j) for j in offset:inner_registration.ncols-offset]
  fill_fluxes!((1,0), inner_indices, inner_registration, surrounding_flood, surrounding_dem, manning_coef, flux_field)

  # bottom boundary
  inner_indices = [Index(inner_registration.nrows-offset,j) for j in offset:inner_registration.ncols-offset]
  fill_fluxes!((-1,0), inner_indices, inner_registration, surrounding_flood, surrounding_dem, manning_coef, flux_field)

  # left boundary
  inner_indices = [Index(j,offset) for j in offset:inner_registration.nrows-offset]
  fill_fluxes!((0,1), inner_indices, inner_registration, surrounding_flood, surrounding_dem, manning_coef, flux_field)

  # right boundary
  inner_indices = [Index(j,inner_registration.ncols-offset) for j in offset:inner_registration.nrows-offset]
  fill_fluxes!((0,-1), inner_indices, inner_registration, surrounding_flood, surrounding_dem, manning_coef, flux_field)

  flux_field

end


"""
    fill_fluxes!(
      direction::Tuple{Int,Int},
      inner_indices::Array{Index,1},
      inner_registration::GeoRegistration,
      surrounding_flood::Grid,
      surrounding_dem::Grid,
      manning_coef::Float64,
      flux_field::Array{Float32, 2}
    )

    TBW
"""
function fill_fluxes!(
  direction::Tuple{Int,Int},
  inner_indices::Array{Index,1},
  inner_registration::GeoRegistration,
  surrounding_flood::Grid,
  surrounding_dem::Grid,
  manning_coef::Float64,
  flux_field::Array{Float32, 2}
)

  # figure out how many cells of the inner grid, linearly, per cell of the outer grid
  # we'll spread the influx across a boundary this many cells wide
  scale_factor = round(Int, surrounding_flood.registration.cell_size_meters / inner_registration.cell_size_meters)

  # convert indices on the boundary of the inset sim into the corresponding indices on the
  # outer flood grid in which the nested simulation is embedded.
  target_indices = convert_source_locations(inner_indices, inner_registration, surrounding_flood.registration, false)

  # compute the flux at each of the locations in the outer flood grid where it bounds
  # the inset simulation
  flux_vector = [
    flux_per_unit_length(surrounding_dem, surrounding_flood, idx, direction, manning_coef)
    for idx in target_indices
  ]

  # note that the '+=' ensures that any flux filled near corners from other boundaries
  # remains included
  for j in eachindex(inner_indices)
    for i in 1:scale_factor
      flux_field[inner_indices[j].row + direction[1]*(i-1), inner_indices[j].col + direction[2]*(i-1)] += flux_vector[j]/scale_factor
    end
  end

end


"""
    flux_per_unit_length(
      dem::Grid,
      flood::Grid,
      idx::Index,
      direction::Tuple{Int,Int},
      manning_coef::Float64
    ) :: Float32

Determine the flux per unit length through the specified grid cell.
Just need to multiply the returned value, then, by the grid-cell width
of the target grid to determine flux.
"""
function flux_per_unit_length(
  dem::Grid,
  flood::Grid,
  idx::Index,
  direction::Tuple{Int,Int},
  manning_coef::Float64
) :: Float32

  h = flood[idx]
  flux = 0.0

  # compute the surface water elevation slope
  try
    ea = flood[idx] + dem[idx]
    eb = flood[idx.row+direction[1], idx.col+direction[2]] + dem[idx.row+direction[1], idx.col+direction[2]]
    s = (ea-eb) / flood.registration.cell_size_meters
    s = s < 0.0 ? 0.0 : s

    # estimate local component of velocity in the specified direction
    # from manning's formula
    v = h^0.67 * sqrt(s) / manning_coef

    flux = v * h

  catch
    println("Not enough info to compute flux on boundary. Using 0.0.")
  end

  flux

end

