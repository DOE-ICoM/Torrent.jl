
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


const DS = DataStructures

"""
Data structure to represent a contaminated area within a DEM.
"""
struct ContaminatedArea

  "Latitude of lower-left corner of contaminated area."
  ll_latitude :: Float32

  "Longitude of lower-left corner of contaminated area."
  ll_longitude :: Float32

  "Latitude of upper-right corner of contaminated area."
  ur_latitude :: Float32

  "Longitude of upper-right corner of contaminated area."
  ur_longitude :: Float32

  "Time step at which contamination becomes active."
  start_step :: Int

  "Time step at which contamination becomes inactive."
  end_step :: Int

  "Density of contaminant that is expected to be picked up by each rivulet traversing this area."
  contamination_rate :: Float32

end

"""
Holds info on a rivulet that is being tracked.
"""
mutable struct RivuletTrack

  "ID of the tracked rivulet."
  id :: Int

  "Track of the rivulet."
  track :: DS.MutableLinkedList{Tuple{Index, Int}}

  "Time step at which tracking of the rivulet began."
  initial_step :: Int

  "Latest time step at which rivulet has been tracked."
  latest_step :: Int

end


"""
    RivuletTrack()

Construct and initial RivuletTrack with an empty instance of a MutableLinkedList and
default values of -1 for the other fields.
"""
function RivuletTrack()
  RivuletTrack(-1, DS.MutableLinkedList{Tuple{Index,Int}}(), -1, -1)
end


"""
The Simulation structure holds information about the context of a simulation,
general characteristics of the rivulets that will be used in the simulation,
and references to the Digital Elevation Model, precipitation time series, and
resulting depth raster.
"""
mutable struct Simulation

  "The time series of geospatial precipitation rates that will drive the flooding."
  source_series::Vector{PrecipitationTimeSeries}

  "Where to store simulation results."
  output_directory::String

  "How frequently (time steps) to write a snapshot of the current depths. [Ex: 240]"
  write_depth_every::Int

  "How frequently (time steps) to update the peak depth for each grid cell. [Ex: 240]"
  compute_max_every::Int

  "Number of threads to use (Scala only). Must be set from the command line in Julia."
  number_of_threads::Int
  
  "Length of each rivulet in grid cells. [Ex: 200]"
  rivulet_length :: Int

  "Thickness of each rivulet in meters. Thinner rivulets meann finer depth resolution. [Ex: 0.05]"
  rivulet_thickness :: Float64

  "Total volume of a rivulet (using the available constructor this will be computed for you)."
  rivulet_volume::Float64

  "Duration of each time step (seconds). [Ex: 60.0]"
  time_step :: Float64

  "Reference to the Digital Elevation Model to use for the simulation."
  dem :: Grid

  "Reference to depth grid (using the available constructor this will be initialized for you)."
  depth :: Array{Float32, 2}

  "Height of the DEM and depth grids."
  height :: Int

  "Width of the DEM and depth grids."
  width :: Int

  "Whether to stop tracking rivulets when they hit a no_data_value cell in the DEM."
  exclude_no_data_cells :: Bool

  "Manning's coefficient may be specified in the configuration file, or if not, will default to 0.04."
  manning_coef :: Union{Float64,Grid}

  "A reentrant lock that may be used to avoid race conditions."
  lk :: ReentrantLock

  "Relative location of neighbors that will be check for next rivulet move (initialized by the constructor)."
  neighbors :: Vector{Tuple{Index,Float64}}

  "An optional instance of a contaminated area."
  contaminated_area :: Union{Nothing,ContaminatedArea}

  "Grid to hold contamination level within a cell."
  contamination :: Union{Nothing,Array{Float32, 2}}

  "Whether to track only contaminated rivulets."
  rivulet_tracking_only_contaminated :: Bool

  "Array of rivulets IDs to track."
  rivulet_tracking_ids :: Vector{Int}

  "Paths of rivulets being tracked."
  rivulet_tracks :: Vector{RivuletTrack}

  "Number of rivulets being tracked so far."
  rivulets_num_tracking :: Int

  "How frequently to track position of rivulets."
  rivulets_tracking_every_time_steps :: Int

end


"""
    Simulation(precipitation::PrecipitationTimeSeries,
      output_directory::String, write_depth_every::Int,
      compute_max_every::Int, number_of_threads::Int,
      rivulet_length::Int, rivulet_thickness::Float64,
      time_step::Float64, manning_representation::Union{Float64,String}, dem::Grid
    )

    A convenience constructor that fills in many of the fields for you, such as a
    zero-initialized depth grid and list of relative neighbor locations.
"""
function Simulation(
  source_series::Vector{PrecipitationTimeSeries},
  output_directory::String, write_depth_every::Int,
  compute_max_every::Int, number_of_threads::Int,
  rivulet_length::Int, rivulet_thickness::Float64,
  time_step::Float64, manning_representation::Union{Float64,String}, dem::Grid,
  exclude_no_data_cells::Bool,
  contaminated_area::Union{Nothing,ContaminatedArea},
  rivulet_tracking_total_num::Int,
  rivulet_tracking_only_contaminated::Bool,
  rivulet_track_every_time_steps::Int
)

  # Depth grid initialized to zero everywhere.
  depth = zeros(Float32, dem.registration.nrows, dem.registration.ncols)

  # Contamination grid initialized to zero everywhere, if a contaminated area has been defined.
  contamination = isnothing(contaminated_area) ? nothing : zeros(Float32, dem.registration.nrows, dem.registration.ncols)

  # Volume of a rivulet
  vol = rivulet_thickness * rivulet_length * dem.registration.cell_size_meters^2.0

  # Set of nearest neighbor locations, including diagonals and their associated
  # center to center distance from the current location in grid cells.
  nbrs = [
    (Index(-1, -1), 1.41),
    (Index(-1, 0), 1.0),
    (Index(-1, 1), 1.41),
    (Index(0, -1), 1.0),
    (Index(0, 1), 1.0),
    (Index(1, -1), 1.41),
    (Index(1, 0), 1.0),
    (Index(1, 1), 1.41)
  ]

  # parses the representation of the manning coefficient. this may be either
  # a constant float value or a string filename representing a raster grid.
  manning_coef = parse_manning(manning_representation, dem.registration.units)

  # setup arrays to hold rivulet tracking information
  rivulet_tracking_ids = fill(-1, rivulet_tracking_total_num)
  rivulet_tracking_paths = [RivuletTrack() for i in 1:rivulet_tracking_total_num]

  # Call the full constructor with these initialized values.
  Simulation(source_series, output_directory, write_depth_every,
    compute_max_every, number_of_threads, rivulet_length, rivulet_thickness,
    vol, time_step, dem, depth, dem.registration.nrows, dem.registration.ncols,
    exclude_no_data_cells, manning_coef, ReentrantLock(), nbrs,
    contaminated_area,
    contamination,
    rivulet_tracking_only_contaminated,
    rivulet_tracking_ids,
    rivulet_tracking_paths,
    0,
    rivulet_track_every_time_steps
  )

end


"""
The rivulet structure holds information about each individual rivulet as it
moves through the depth/dem grids.
"""
mutable struct Rivulet

  "A unique identifier for each rivulet."
  id :: Int

  "The current index of the head position within the path array."
  head_idx :: Int

  "The number of steps the rivulet has traveled so far."
  steps_taken :: Int

  """An array holding the patial locations of each element of the rivulet, ie, the path it has followed.
      going back length(path) steps.
  """
  path :: Array{Union{Index, Nothing}}

  "Whether or not the rivulet head has run into a boundary yet."
  head_has_hit_boundary :: Bool

  "Density of any contaminant within the rivulet."
  contaminant_density :: Float32

end


"""
    Rivulet(sim::Simulation, id::Int, location::Index, rivulet_length::Int)

    Convenience constructor for Rivulet information.
"""
function Rivulet(sim::Simulation, id::Int, location::Index, contaminant_density::Float32 = 0f0)
  head_idx = 1
  path = Array{Union{Index,Nothing}}(undef, sim.rivulet_length)
  fill!(path, nothing)
  path[head_idx] = location  # add initial location to the path at the head idx
  Rivulet(id, head_idx, 0, path, is_on_boundary(sim, location), contaminant_density)
end


"""
    move_head!(rivulet :: Rivulet) :: Union{T,Nothing}

    Move the rivulet forward one step. This has the side effect of updating
    the depth grid. The new location of the rivulet head is returned, unless
    the rivulet has previously hit the boundary, in which case nothing is returned.
"""
function move_head!(sim::Simulation, rivulet::Rivulet, velocity::Float64, sim_step::Int)::Union{Index,Nothing}

  if(rivulet.head_has_hit_boundary)
    return nothing
  else

    new_location = lowest_surface_nbr_abs(sim, rivulet)

    # if this neighbor is on the boundary we need to remember that so
    # we don't try to move the head further in the future, instead just
    # letting the tail complete its evolution.
    if is_on_boundary(sim, new_location)
      rivulet.head_has_hit_boundary = true
    end

    # update the depth grid with additional rivulet_thickness depth
    # at the new head location.
    array_update(sim.depth, new_location, sim.lk) do h
      h + sim.rivulet_thickness
    end

    # note whether or not we move through a contaminated area on this step
    hit_contaminated_area = false

    # handle contamination if a contamination area has been specified
    if !isnothing(sim.contamination)

      # if we're in a contaminated area, we need to add additional contamination to the rivulet
      if is_contaminated_area(new_location, sim_step, sim.contaminated_area, sim.dem.registration)
        rivulet.contaminant_density += sim.contaminated_area.contamination_rate * sim.time_step * sim.rivulet_thickness
        hit_contaminated_area = true
      end

      # then need to add aditional contamination to the contaminant grid at this location
      array_update(sim.contamination, new_location, sim.lk) do c
        c + rivulet.contaminant_density 
      end

    end

    # initiate tracking of this rivulet if rivulet and sim state meet criteria
    maybe_initiate_rivulet_tracking(sim, rivulet, hit_contaminated_area, sim_step)  # note, this function has many side effects on the sim structure

    # update rivulet tracking if necessary (this could become an onerous test if the number
    # of rivulets being tracked is large). for most debugging purposes n should be small,
    # though, making the computational burden negligible.
    if mod(sim_step, sim.rivulets_tracking_every_time_steps) == 0
      i = findfirst(x -> x==rivulet.id, sim.rivulet_tracking_ids)
      if !isnothing(i)
        DS.append!(sim.rivulet_tracks[i].track, (new_location, sim_step) )
        sim.rivulet_tracks[i].latest_step = sim_step
      end
    end

    # return the location of the new head
    new_location
  end
end


"""
    maybe_initiate_rivulet_tracking(sim::Simulation, rivulet::Rivulet, in_contaminated_area::Bool)

Test whether the conditions for tracking this rivulet are met, i.e., it is or isn't in a contaminated
area and we haven't yet met the overall quota for the number of rivulets we're supposed to track.
"""
function maybe_initiate_rivulet_tracking(sim::Simulation, rivulet::Rivulet, in_contaminated_area::Bool, sim_step::Int)

  is_time = if in_contaminated_area && !isnothing(sim.contaminated_area)
    sim_step >= sim.contaminated_area.start_step && sim_step <= sim.contaminated_area.end_step
  else
    true
  end

  if ( sim.rivulets_num_tracking < size(sim.rivulet_tracking_ids,1) &&
    ( (sim.rivulet_tracking_only_contaminated && in_contaminated_area) || 
      !sim.rivulet_tracking_only_contaminated ) &&
      is_time
  )

    if !(rivulet.id in sim.rivulet_tracking_ids)
      # println("\n===========================")
      # println("Tracking new rivulet: $(rivulet.id)")
      # println("===========================\n")
      sim.rivulets_num_tracking += 1
      sim.rivulet_tracking_ids[sim.rivulets_num_tracking] = rivulet.id
      sim.rivulet_tracks[sim.rivulets_num_tracking].id = rivulet.id
      sim.rivulet_tracks[sim.rivulets_num_tracking].initial_step = sim_step
    end
  end

end


"""
    dry_tail!(sim::Simulation, rivulet::Rivulet, tail_location::Index, tail_idx::Int) :: Bool

    Moves the tail of the rivulet forward one step. This has the side effect
    updating the depth grid.
"""
function dry_tail!(sim::Simulation, rivulet::Rivulet, tail_location::Index, tail_idx::Int) :: Bool

  # update the depth grid, subtracting rivulet_thickness
  # from the location of the rivulet tail.
  if array_read(sim.depth, tail_location, sim.lk) >= sim.rivulet_thickness
    array_update(sim.depth, tail_location, sim.lk) do h
      h - sim.rivulet_thickness
    end
  end

  # and, if we're managing contamination, need to subtract that off from the contaminant grid
  # could end up in a situation where a rivulet has accumulated more contamination at the head
  # than it originally contributed to a contaminant grid cell. We'll constrain the grid to not
  # be less than zero for the moment.
  if(!isnothing(sim.contamination) && rivulet.contaminant_density > 0.0)  # used to be 0.0f which introduced a bug that was then inadvertently swallowed by the try/catch block in SimulationRunner
    array_update(sim.contamination, tail_location, sim.lk) do c
      c - (rivulet.contaminant_density <= c ? rivulet.contaminant_density : c)
    end
  end

  # in the case where we're drying the tail of a rivulet that has hit
  # the boundary, it's important that we set this location to nothing
  # so we know to skip over it and dry the next element in the path in
  # the next iteration
  rivulet.path[tail_idx] = nothing

  # return whether the tail is now on the boundary
  is_on_boundary(sim, tail_location)  # || tail_location == head_location(rivulet)

end


"""
    lowest_surface_nbr_relative(rivulet::Rivulet) :: Index

    Determine the relative indices of the neighbor with the lowest slope from
    the current rivulet head location.
"""
function lowest_surface_nbr_relative(sim::Simulation, rivulet::Rivulet) :: Index
  
  # figure out which neighbor has the minimum slope from the current head location
  row, col = tuplefrom(rivulet.path[rivulet.head_idx])
  min_nbr = (min_by(sim.neighbors) do nbr  # typeof(nbr) == (Index, Float64)

    # the neighbor location is the first element of
    # the tuple and the relative distance the second
    j, i = tuplefrom(nbr[1])
    distance = nbr[2]

    # if the neighbor location is outside the boundary make the apparent
    # slope super high so we don't choose to step there
    if is_outside_boundary(sim, Index(row+j, col+i))
      1000000.0

    # otherwise compute the relative slope
    else
      ((sim.dem[row+j, col+i] + sim.depth[row+j, col+i]) - 
        (sim.dem[row, col] + sim.depth[row, col])) / distance
    end

  end)[1]  # and retain the Index of the neighbor with the minimum slope

  # grab this neighbor location
  dy, dx = tuplefrom(min_nbr)

  # and compute its surface elevation
  nbr_surface_elevation = sim.dem[row+dy, col+dx] + sim.depth[row+dy, col+dx]

  # as well as the local surface elevation. technically, we could probably
  # just pass this out of the min_by calculation above someway to avoid
  # recomputing, but i'm not sure that's necessary
  local_surface_elevation = sim.dem[row, col] + sim.depth[row, col]

  # return the neighbor location unless we're in a dip, in which case
  # fill in the current location of the head first
  nbr_surface_elevation < local_surface_elevation ? min_nbr : Index(0,0)

end


"""
    lowest_surface_nbr_abs(rivulet::Rivulet) :: Index

    Determine the absolute indices of the neighbor with the lowest slope from
    the current rivulet head location.
"""
function lowest_surface_nbr_abs(sim::Simulation, rivulet::Rivulet) :: Index
  nbr = lowest_surface_nbr_relative(sim, rivulet)
  Index(head_location(rivulet).row + nbr.row, head_location(rivulet).col + nbr.col)
end


"""
    min_by(f::Function, v::Vector{T}) where T

    Returns the element of v which gives the minimum value when passed to f.
    Neither Julia's min or minimum functions act in quit this way.
"""
function min_by(f::Function, v::Vector{T}) where T

  # initialize minimum element and values to first element of v
  min_elem = v[1]
  min_val = f(v[1])

  # step through other elements, apply f, and determine which one,
  # if any, results in the lowest value
  for i in 2:length(v)
    x = f(v[i])
    if x < min_val
      min_elem = v[i]
      min_val = x
    end
  end

  # return the element of v that had the smallest value when
  # f was applied to it
  min_elem
end


"""
    head_location(rivulet :: Rivulet) :: Index

    Get grid indices of current rivulet head location.
"""
function head_location(rivulet :: Rivulet) :: Index
  rivulet.path[rivulet.head_idx]
end


"""
    tail_location(sim::Simulation, rivulet::Rivulet) :: Index

    Get grid indices of current rivulet tail location.
"""
function tail_location(sim::Simulation, rivulet::Rivulet) :: Index
  rivulet.path[path_idx_after(sim, rivulet.head_idx)]
end


"""
    path_idx_after(sim::Simulation, i::Int) :: Int

    Determine the index of the path element that follows i. This
    function takes care of wrapping when you approach the end of
    the array.
"""
function path_idx_after(sim::Simulation, i::Int) :: Int
  i == sim.rivulet_length ? 1 : i+1
end


"""
    path_idx_before(sim::Simulation, i::Int) :: Int

    Determine the index of the path element that precedes i. This
    function takes care of wrapping when you approach the end of
    the array.
"""
function path_idx_before(sim::Simulation, i::Int) :: Int
  i == 1 ? sim.rivulet_length : i-1
end



"""
    is_on_boundary(sim::Simulation, location::Index) :: Bool

    Determines whether location is on a simulation boundary row/col.
"""
function is_on_boundary(sim::Simulation, location::Index) :: Bool
  location.row == 1 ||
  location.col == 1 ||
  location.row == size(sim.dem.data, 1) ||
  location.col == size(sim.dem.data, 2) ||
  (sim.exclude_no_data_cells && typeof(sim.dem.registration.no_data_value) == Float64 ? sim.dem[location] == sim.dem.registration.no_data_value : false)
end


"""
    is_outside_boundary(sim::Simulation, location::Index)::Bool

    Determines whether locations is outside the simulation boundary.
"""
function is_outside_boundary(sim::Simulation, location::Index)::Bool
  location.row < 1 || location.col < 1 || location.row > sim.height || location.col > sim.width
end



"""
    grid_cells_to_move(sim::Simulation, a::Index, ob::Union{Index,Nothing}) :: Int

    Determines the number of grid cells to move the head based on the current
    location, a, and previous location, ob. If nothing is passed for ob then the
    rivulet is moved a default of 1 cell. A tuple is returned containing the number
    of grid cells to move and the actual, calculated velocity.
"""
function grid_cells_to_move(sim::Simulation, a::Index, ob::Union{Index,Nothing}) :: Tuple{Int,Float64}
  
  # if there's no second index supplied, just return 1
  isnothing(ob) && return (1, sim.dem.registration.cell_size_meters / sim.time_step)

  # compute the surface elevations and slope
  sea = sim.dem[a] + array_read(sim.depth, a, sim.lk)
  seb = sim.dem[ob] + array_read(sim.depth, ob, sim.lk)
  s = -(sea - seb) / sim.dem.registration.cell_size_meters
  s = s < 0.0 ? 0.0 : s

  # velocity from manning's formula
  v = array_read(sim.depth, a, sim.lk)^0.67 * sqrt(s) / manning(sim.manning_coef,a)

  # this next line deserves a lot of discussion and additional research. it is
  # to account for the fact that ...
  # v = v < 1.0 ? 1.0 : v  # where 1.0 = 1.0 meters per second

  # compute the ideal number of grid cells to move based on the floating point velocity
  ideal_number = v * sim.time_step / sim.dem.registration.cell_size_meters

  # if the ideal number is less than one cell, move at least one cell
  # (may want to REMOVE THIS, but I think it was added to keep us from getting
  # weird structures forming at some places in the east of the Boulder scenario)
  actual_number = if ideal_number < 1.0
    1

  # otherwise, we'll deterministically move the full integer number of steps
  # and stochastically move one additional step with probability determined by
  # the fractional part of the floating point velocity. this ensures that on average
  # we move the right number of grid cells, without having to deal with tracking
  # partial, intra-cell positions.
  else
    i = trunc(Int, ideal_number)
    frac = ideal_number - i
    rand() < frac ? i+1 : i
  end

  # return this number of grid cells to move
  (actual_number, v)

end


"""
    step(sim::Simulation, rivulet::Rivulet) :: Bool

    Perform a one time step update of the rivulet's location and path.
"""
function step(sim::Simulation, rivulet::Rivulet, sim_step::Int) :: Tuple{Bool, Float64}

  # get the number of grid cells to move based on the local slope
  # and manning velocity formula
  cells, vel = grid_cells_to_move(
    sim,
    head_location(rivulet),  # current head location
    rivulet.steps_taken >= 1 ? rivulet.path[path_idx_before(sim, rivulet.head_idx)] : nothing  # previous head location
  )

  # move that many cells forward this time step
  cells_moved = 0
  tail_on_boundary = false
  while cells_moved < cells && !tail_on_boundary

    # if we're at least a rivulet length into this rivulet's travels
    # or the rivulet's head has already hit the boundary, we need to dry
    # the tail
    if rivulet.steps_taken > sim.rivulet_length-1 || rivulet.head_has_hit_boundary

      # find presumptive location of the tail
      tail_idx = path_idx_after(sim, rivulet.head_idx)
      tail = rivulet.path[tail_idx]

      # manage case where head hits boundary before full rivulet
      # length has unfolded by scanning through to find actual
      # current end. note that, as written, this essentially
      # removes water from the field sooner than it should be
      # removed. might improve upon this later by letting full
      # length of rivulet emerge before we start erasing it.
      while isnothing(tail) && tail_idx != rivulet.head_idx
        tail_idx = path_idx_after(sim, tail_idx)
        tail = rivulet.path[tail_idx]
      end

      # then dry tail
      tail_on_boundary = dry_tail!(sim, rivulet, tail, tail_idx)

    end  # if

    # now that the tail has been dried, we can move the head of
    # the rivulet to its next location
    new_head = move_head!(sim, rivulet, vel, sim_step)

    # so long as the head hasn't hit the boundary, in which case
    # move_head! would return a value of nothing
    if !isnothing(new_head)

      # update the index of the head in the path array
      rivulet.head_idx = path_idx_after(sim, rivulet.head_idx)

      # and store the new head location in this element
      rivulet.path[rivulet.head_idx] = new_head
    end

    # update the number of cells moved and the total number of
    # steps the rivulet has taken so far
    cells_moved += 1
    rivulet.steps_taken += 1

  end  # while

  # return a boolean reflecting whether or not the rivulet tail
  # has reached the boundary. if it has, the rivulet can be killed
  # off in the next simulation time step.
  (tail_on_boundary, vel)

end


"""
    parse_manning(x::Float64)::Float64

    Simply returns the input value.
"""
function parse_manning(x::Float64, units::String)::Float64
  x
end


"""
    parse_manning(filename::String)::Grid

    Opens the file, `filename`, which is assumed to represent a raster,
    and returns a Grid instance.
"""
function parse_manning(filename::String, units::String)::Grid
  open_raster(filename; units=units)
end


"""
    manning(value_holder::Float64, row::Int, col::Int)::Float64

    If the value_holder is a Float64, the value_holder is returned.
"""
function manning(value_holder::Float64, row::Int, col::Int)::Float64
  value_holder
end


"""
    manning(value_holder::Float64, index::Index)::Float64

    If the value_holder is a Float64, the value_holder is returned.
"""
function manning(value_holder::Float64, index::Index)::Float64
  value_holder
end


"""
    manning(value_holder::Grid, row::Int, col::Int)::Float64

    If the value_holder is a grid, the row/col indices are used to extract
    the manning coefficient at the location specified.
"""
function manning(value_holder::Grid, row::Int, col::Int)::Float64
  value_holder[row, col]
end


"""
    manning(value_holder::Grid, index::Index)::Float64

    If the value_holder is a grid, the indices, index, are used to extract
    the manning coefficient at the location specified.
"""
function manning(value_holder::Grid, index::Index)::Float64
  value_holder[index]
end


"""
    is_in_contaminated_area(location::Index, contaminated_area::ContaminatedArea)

Whether `location` lies within a contaminated area and time.
"""
function is_contaminated_area(
  location::Index,
  sim_step::Int,
  ca::ContaminatedArea,
  registration::GeoRegistration
) :: Bool

  if sim_step < ca.start_step || sim_step > ca.end_step
    false
  else
    lat, lon = latlongof(location, registration)
    lat >= ca.ll_latitude && lat <= ca.ur_latitude && lon >= ca.ll_longitude && lon <= ca.ur_longitude
  end
end