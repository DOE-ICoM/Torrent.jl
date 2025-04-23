
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
    torrent(config_filename::String)

Run a Torrent simulation based on the configuration information in
the referenced JSON-format configuration file.
"""
function torrent(config_filename::String)

  # open and parse the configuration file
  config_io = JSON.open(config_filename)
  config = JSON.parse(config_io)

  # open the base Digital Elevation Model
  println();
  dem = open_raster(config["dem"]["filename"], false; units=config["dem"]["units"]); println()
  grid_size_message = Printf.format(Printf.Format("Grid Resolution: %.1f meters ($(size(dem.data,1)) x $(size(dem.data,2)) cells)\n"), dem.registration.cell_size_meters)
  printstyled(grid_size_message; color=39)

  # FastFlood is stochastic. Multiple runs can be requested using the same
  # set of configuration parameters to generate an ensemble that characterizes
  # model forecast uncertainty.
  num = if haskey(config, "num-realizations") config["num-realizations"] else 1 end
  for iter in 1:num

    println()
    printstyled("========================================\n"; color=214)
    printstyled("    ITERATION: $iter\n"; color=214)
    printstyled("========================================\n"; color=214)
    println()
    realization(config, dem, iter)
  end
end


"""
    realization(config::Dict{String,Any}, dem::Grid, iteration::Int64)

    Performs a single realization of a Simulation.
"""
function realization(config::Dict{String,Any}, dem::Grid, iteration::Int)

  # default manning coefficient to use if none is specified in the configuration file
  default_manning_coef = 0.04

  # load the precipitation data for the full simulation time.
  (timing_of_precip, precipitation) = read_precipitation(config, dem.registration)

  # create sources based on a dam failure hydrograph
  (timing_of_dam_failure, dam_failure) = generate_dam_failure(config, dem.registration)

  # load any surrounding flood data that should be used as a boundary condition
  (timing_of_bounding_flood, bounding_flood) = read_boundary_conditions(config, dem.registration, default_manning_coef)

  # create a list of sources where any flood boundary conditions are optional
  source_series :: Vector{PrecipitationTimeSeries} = filter( x -> !isnothing(x), [precipitation, bounding_flood, dam_failure])

  # create an instance of a contaminated area if defined in the configuration
  contaminated_area :: Union{Nothing,ContaminatedArea} = contamination(config)

  # create a Simulation instance to store all of the relevant data
  sim = Simulation(
    source_series,
    config["output-directory"],
    config["save-every-steps"],
    config["compute-max-every-steps"],
    1,  # num threads is no longer used in the Julia version of the code.
    config["rivulet-length"],
    config["rivulet-thickness"],
    config["time-step-seconds"],
    haskey(config, "manning-coef") ? config["manning-coef"] : default_manning_coef,
    dem,
    haskey(config, "exclude-no-data-cells") ? config["exclude-no-data-cells"] : true,
    contaminated_area,
    haskey(config, "rivulet-tracking") ? config["rivulet-tracking"]["num-rivulets"] : 0,
    haskey(config, "rivulet-tracking") ? config["rivulet-tracking"]["only-contaminated"] : false,
    haskey(config, "rivulet-tracking") ? config["rivulet-tracking"]["track-every-time-steps"] : -1,
    haskey(config, "interpolate-output") ? config["interpolate-output"] : false
  )

  # run the simulation
  (timing_of_run, time_spent_saving) = time_computation("Running simulation...") do
    run(sim, config["max-steps"], iteration)
  end

  # report timing info to stdout as well as saving to a file in the output directory
  write_timing_info(config, timing_of_precip, timing_of_bounding_flood, timing_of_dam_failure, time_spent_saving, timing_of_run, iteration)

end  # function realization


"""
    run(sim::Simulation, num_time_steps::Int)

Run the simulation with the parameters specified by the Simulation instance.
"""
function run(sim::Simulation, num_time_steps::Int, iteration::Int)

  # variable to track the current simulation step
  sim_step = 0

  # the index of the next flux source to use, one for each flux type.
  # these have already been randomly spatially shuffled so we can just step
  # through them to get a pseudo-random distribution of rain fall at each
  # time step
  current_source_index_array = ones(Int, length(sim.source_series))

  # create an array to hold the set of currently operational rivulet instances
  rivulets = Array{Rivulet}(undef, 0)

  # track number of rivulets that we're finished simulating because
  # they've left the simulation region
  num_rivulet_escapes = 0

  # a unique identifier for the rivulets and tracker of the total
  # number of rivulets spawned so far
  rivulet_counter = 1

  # initialize a table to hold a running tabulation of current and
  # cumulative rivulet numbers, the associated volume, and related flux
  flux_trace = Matrix{Float64}(undef, num_time_steps, 6)

  # initialize a grid to track peak flood depth over the course
  # of the simulation
  max_depth = zeros(Float32, size(sim.depth))

  # initialize grid to track peak contamination over the course
  # of the simulation
  max_contamination = if !isnothing(sim.contaminated_area)
    zeros(Float32, size(sim.depth))
  else
    nothing
  end

  # track total time spent saving so we can separate that from algorithm
  # simulation time
  total_save_time = 0.0

  # support reporting a rough processing rate
  prev_time = time()

  # track the number of new rivulets created
  num_rivulet_starts = 0

  # debug
  num_escapes_by_error = 0

  # create a progress bar to chart precipitation load and generation
  prog = Progress(num_time_steps; desc="Torrent             ")

  # this loop is the guts of the simulation where we run the simulation
  # through the requested number of time steps
  while sim_step < num_time_steps

    # every so often, we'll update the run status for the user

    # let's compute the rivulet processing rate; typically, the more
    # rivulets the more copying that needs to be done of rivulet
    # arrays, and the slower the processing rate
    proc_rate = length(rivulets)/(time()-prev_time)
    proc_rate_str = Printf.format(Printf.Format("%.1e"), proc_rate)

    # now we can update the status bar
    update!(prog, sim_step; showvalues = [
      ("step", sim_step),
      ("rivulets", length(rivulets)),
      ("cumulative rivulets", rivulet_counter),
      ("processing rate (riv/sec)", proc_rate_str)
    ])
    prev_time = time()

    # and then save some flux information that we can report out at the
    # end for the user
    flux_trace[sim_step+1, 1] = length(rivulets)
    flux_trace[sim_step+1, 2] = length(rivulets) * sim.rivulet_volume
    flux_trace[sim_step+1, 3] = rivulet_counter
    flux_trace[sim_step+1, 4] = rivulet_counter * sim.rivulet_volume
    flux_trace[sim_step+1, 5] = num_rivulet_starts
    flux_trace[sim_step+1, 6] = num_rivulet_starts * sim.rivulet_volume / sim.time_step

    # if you only wanted to run a single thread, the single line below
    # would suffice to do what we needed on a Set of rivulets
    #   filter!((riv) -> !step(sim.structure, riv), rivulets)

    # when multithreading in Julia, though, we need to track rivulets
    # with an array and use a for-loop to iterate through the invocation
    # of each rivulet's step function

    # while doing this we'll need to track whether or not the rivulet
    # has left the simulation region
    rivulet_escaped = Array{Bool}(undef, length(rivulets))

    # as well as the total number of escapees (by tracking this on the 
    # fly we don't have to tally it again later, thereby avoiding a
    # second loop through the array)
    num_rivulet_escapes = 0

    # use multithreading to invoke the step function on each rivulet
    # in the rivulets array
    @threads for j in 1:length(rivulets)

      # When run without the @threads macro, this loop works fine.
      # with the @threads macro it randomly throws an UndefRefError
      # on the last element of the array. don't know why. at this point
      # we're simply catching the error and marking that rivulet as
      # as escaped so that it can be discarded the next time step.
      # this really just seems to effect one rivulet every few tens
      # of time steps. given that we often have tens or hundreds of
      # thousands of rivulets, the impact should be minimal.
      try

        # step a rivulet forward in time
        did_escape, velocity = step(sim, rivulets[j], sim_step)

        # make note of whether or not the rivulet escaped and add to
        # the total number of escapees if appropriate
        rivulet_escaped[j] = did_escape
        if did_escape
          num_rivulet_escapes += 1
        end

      # just swallow the UndefRefError at this point
      catch err
        if isa(err, UndefRefError)
          # println("DEBUG: undefined found at $j of $(length(rivulets)) rivulets")
          rivulet_escaped[j] = true
          num_rivulet_escapes += 1
          num_escapes_by_error += 1
        end
      end

    end  # multithreaded for-loop

    # now that we've updated the state of existing rivulets, we'll start to
    # look toward setting things up for the next iteration. we get all source
    # events corresponding to the current simulation time step
    source_events = map(x -> get_current_period(x, sim_step), sim.source_series)

    # and compute the number of new rivulets that should be created to
    # account for each of the flux events (e.g., rain, boundary flux).
    num_rivulet_starts_by_event = map(source_events) do evt

      # calculate ideal number of new rivulets based on flux, time step,
      # and rivulet volume
      ideal_rivulet_starts = sim.time_step * evt.total_flux / sim.rivulet_volume

      # truncate to the integer number of whole rivulets
      int_part_of_ideal = trunc(Int, ideal_rivulet_starts)

      # and then use the fraction part to randomly add one more so that
      # statistically we should get the right total flux over time and
      # not lose volume flux to truncation (which builds over time)
      frac = ideal_rivulet_starts - int_part_of_ideal
      rand() < frac ? int_part_of_ideal+1 : int_part_of_ideal
    end

    # calculate the total number of rivulet starts across all ongoing events
    num_rivulet_starts = sum(num_rivulet_starts_by_event)

    # now that we know how many new rivulets will be initiated we can
    # create the new rivulet array. this will directly replace the prior 
    # rivulet array. its length is the number of rivulets holding over from
    # the prior time step (those that didn't leave the simulation area) plus
    # the number of new rivulet starts 
    new_rivulet_array = Array{Rivulet}(undef, length(rivulets) - num_rivulet_escapes + num_rivulet_starts)

    # first copy the existing rivulets that didn't escape to this new array
    count = 1
    for j in 1:length(rivulets)
      if !rivulet_escaped[j]
        new_rivulet_array[count] = rivulets[j]
        count += 1
      end
    end

    # now, let's step through each current, ongoing flux event
    flux_sources_so_far = 0
    for j in eachindex(source_events)

      # create new rivulet instances for each one and add them to the new rivulet array
      for i in (count+flux_sources_so_far) : (count+flux_sources_so_far) + num_rivulet_starts_by_event[j] - 1

        # remember the current_source_index steps us through our randomly shuffled
        # source array. just want to make sure the counter stays in bounds.
        if current_source_index_array[j] >= length(source_events[j].sources)
          current_source_index_array[j] = 1
        end

        # get the location of the next rain source
        location = source_events[j].sources[current_source_index_array[j]]
        current_source_index_array[j] = current_source_index_array[j] < length(source_events[j].sources) ? current_source_index_array[j]+1 : 1

        # and create a new rivulet instance at that location, assigning
        # it to the new_rivulet_array
        new_rivulet_array[i] = Rivulet(sim, rivulet_counter, location)

        # update our total rivulet counter
        rivulet_counter += 1

      end

      # note number of new rivulets for this event in the total so far
      flux_sources_so_far += num_rivulet_starts_by_event[j]

    end

    # finally we can update the rivulets reference to point to the new array
    # note that the single threaded version in Julia and multithreaded version
    # in scala use Set containers to avoid copying arrays, but multithreading
    # in Julia doesn't make that approach so easy so we end up with an additional
    # array copy.
    rivulets = new_rivulet_array

    # write current snapshot of depth and contaminant grid(s) if requested
    if sim.write_depth_every > 0 && sim_step > 0

      if mod(sim_step, sim.write_depth_every) == 0
        (save_time, _) = time_computation("") do 

          # compute contaminant concentration rather than absolute amount
          concentration = if !isnothing(sim.contamination)
            broadcast((c,h) -> h > 0.0 ? Float32(c/h) : 0.0f0, sim.contamination, sim.depth)
          else
            nothing
          end

          filename_index = Printf.@sprintf("%03d-%05d", iteration, sim_step)

          if strip(sim.dem.registration.proj_string) != "" 
            save_geotiff(sim.output_directory*"depth-$filename_index.tif", sim.depth, sim.dem.registration, sim.interpolate_output)
            if !isnothing(sim.contamination)
              save_geotiff(sim.output_directory*"concentration-$filename_index.tif", concentration, sim.dem.registration, sim.interpolate_output)
              save_geotiff(sim.output_directory*"contamination-$filename_index.tif", sim.contamination, sim.dem.registration, sim.interpolate_output)
            end
          else
            save_esri_asc_file(sim.output_directory*"depth-$filename_index.asc", sim.depth, sim.dem.registration, sim.interpolate_output)
            if !isnothing(sim.contamination)
              save_esri_asc_file(sim.output_directory*"concentration-$filename_index.asc", concentration, sim.dem.registration, sim.interpolate_output)
              save_esri_asc_file(sim.output_directory*"contamination-$filename_index.asc", sim.contamination, sim.dem.registration, sim.interpolate_output)
            end
          end
        end        

        # and track the total time spent saving
        total_save_time += save_time
      end

    end

    # we may also need to recompute the peak depth reached so far
    # don't want to do this too often (say, every time step) as it's 
    # fairly computationally intensive
    if mod(sim_step, sim.compute_max_every) == 0
      max_in_place!(max_depth, sim.depth)
      if !isnothing(max_contamination); max_in_place!(max_contamination, sim.contamination) end
    end

    # update the simulation time step we're on
    sim_step += 1

  end  # while

  finish!(prog)

  padded_iteration = Printf.@sprintf("%03d", iteration)

  # save the final peak depth and contamination grids
  (save_time, _) = time_computation("Saving peak depth and contamination ", false) do 

    # compute contaminant concentration rather than absolute amount
    max_concentration = if !isnothing(sim.contamination)
      broadcast((c,h) -> h > 0.0 ? Float32(c/h) : 0.0f0, max_contamination, max_depth)
    else
      nothing
    end

    if strip(sim.dem.registration.proj_string) != "" 
      save_geotiff(sim.output_directory * "peak-depth-$padded_iteration.tif", max_depth, sim.dem.registration, sim.interpolate_output)
      if !isnothing(sim.contamination)
        # save_geotiff(sim.output_directory*"peak-concentration-$iteration.tif", max_concentration, sim.dem.registration, true)
        save_geotiff(sim.output_directory*"peak-contamination-$padded_iteration.tif", max_contamination, sim.dem.registration, sim.interpolate_output)
      end
    else
      save_esri_asc_file(sim.output_directory * "peak-depth-$padded_iteration.asc", max_depth, sim.dem.registration, sim.interpolate_output)
      if !isnothing(sim.contamination)
        # save_esri_asc_file(sim.output_directory*"peak-concentration-$iteration.asc", max_concentration, sim.dem.registration, true)
        save_esri_asc_file(sim.output_directory*"peak-contamination-$padded_iteration.asc", max_contamination, sim.dem.registration, sim.interpolate_output)
      end
    end
  end

  # save the final paths of any rivulets being tracked
  (rivulet_tracks_save_time, _) = time_computation("Saving any rivulet tracks ", false) do 
    save_rivulet_tracks(sim, iteration)
  end

  # save the rivulet and volume related time series
  (flux_trace_save_time, _) = time_computation("Saving flux trace ", false) do 
    save_csv(sim.output_directory * "flux-$padded_iteration.csv", flux_trace,
      ["Current Rivulets", "Current Volume", "Cumulative Rivulets", "Cumulative Volume", "New Rivulets", "Flux"]
    )
  end

  # update the total save time
  total_save_time += save_time + rivulet_tracks_save_time + flux_trace_save_time

  # return the total time spent saving results as the output
  # of the run function
  total_save_time

end


"""
    save_rivulet_tracks(sim::Simulation)

TBD
"""
function save_rivulet_tracks(sim::Simulation, iteration::Int)
  if length(sim.rivulet_tracking_ids) > 0
    num_rows = maximum(rt -> DS.length(rt.track), sim.rivulet_tracks)
    num_tracks = length(sim.rivulet_tracking_ids)
    locations = fill(0.0, num_rows, 3*num_tracks)
    for j in 1:num_tracks
      println("Tracked rivulet id: $(sim.rivulet_tracks[j].id), initial_step: $(sim.rivulet_tracks[j].initial_step), latest_step: $(sim.rivulet_tracks[j].latest_step)")
      i = 1
      for (elem, t) in sim.rivulet_tracks[j].track
        lat, lon = latlongof(elem, sim.dem.registration)
        locations[i, 3*j] = t
        locations[i, 3*j-1] = lat
        locations[i, 3*j-2] = lon
        i += 1
      end
    end
    save_csv(sim.output_directory*"rivulet-tracks-$iteration.csv", locations)
  end
end


"""
    read_precipitation(
      config::Dict{String,Any}, 
      registration::GeoRegistration
    )::Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

TBW
"""
function read_precipitation(
  config::Dict{String,Any}, 
  registration::GeoRegistration
)::Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

  time_computation("Opening precipitation time series...\n ") do 

    if haskey(config, "rain-time-series")
      open_precipitation_series(
        config["rain-time-series"]["filename-pattern"],
        config["rain-time-series"]["min-index"],
        config["rain-time-series"]["max-index"],
        config["num-sources"],
        config["time-step-seconds"],
        config["rain-time-series"]["file-interval-seconds"],
        config["write-precipitation-distributions"],
        config["output-directory"],
        registration
      )
    elseif haskey(config, "rain-nwm")
      open_nwm_data(
        config["rain-nwm"]["filename"],
        registration,
        Float64(config["time-step-seconds"]),
        config["rain-nwm"]["flux-interval-seconds"],
        config["num-sources"]
      )
    elseif haskey(config, "rain-multiband-geotiff")
      lower = if haskey(config["rain-multiband-geotiff"], "min-index") config["rain-multiband-geotiff"]["min-index"] else -1 end
      upper = if haskey(config["rain-multiband-geotiff"], "max-index") config["rain-multiband-geotiff"]["max-index"] else -1 end
      open_multiband_precipitation(
        config["rain-multiband-geotiff"]["filename"],
        lower,
        upper,
        config["num-sources"],
        config["time-step-seconds"],
        config["rain-multiband-geotiff"]["band-interval-seconds"],
        config["write-precipitation-distributions"],
        config["output-directory"],
        registration
      )
    else
      nothing
    end
  end
end


"""
    read_boundary_conditions(
      config::Dict{String,Any},
      registration::GeoRegistration,
      default_manning_coef::Float64
    )::Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

TBW
"""
function read_boundary_conditions(
  config::Dict{String,Any},
  registration::GeoRegistration,
  default_manning_coef::Float64
)::Tuple{Float64,Union{PrecipitationTimeSeries,Nothing}}

  time_computation("Opening any bounding flood series...\n ") do 

    if haskey(config, "boundary-conditions-time-series")
      open_bounding_flood_series(
        config["boundary-conditions-time-series"]["filename-pattern"],
        config["boundary-conditions-time-series"]["min-index"],
        config["boundary-conditions-time-series"]["max-index"],
        config["boundary-conditions-time-series"]["index-step"],
        config["boundary-conditions-time-series"]["index-offset"],
        config["boundary-conditions-time-series"]["seconds-per-step"],
        config["time-step-seconds"],
        config["boundary-conditions-time-series"]["dem"]["filename"],
        config["boundary-conditions-time-series"]["dem"]["units"],
        config["num-sources"],
        config["boundary-conditions-time-series"]["inset-from-border"],
        config["write-precipitation-distributions"],
        config["output-directory"],
        registration,
        default_manning_coef
      )
    else
      nothing
    end
  end
end



"""
    write_timing_info(
      config::Dict{String,Any},
      timing_of_precip::Float64,
      timing_of_bounding_flood::Float64,
      time_spent_saving::Float64,
      timing_of_run::Float64
    )

Save information about the timing of different aspects of the simulation
realization, including a copy of the detailed configuration parameters used
in the simulation.
"""
function write_timing_info(
  config::Dict{String,Any},
  timing_of_precip::Float64,
  timing_of_bounding_flood::Float64,
  timing_of_dam_failure::Float64,
  time_spent_saving::Float64,
  timing_of_run::Float64,
  iteration::Int
)

  # and report some timing information
  println("Time to analyze precipitation: $timing_of_precip")
  println("Time to analyze flood boundary conditions: $timing_of_bounding_flood")
  println("Time to generate dam breach hydrograph: $timing_of_dam_failure")
  println("Time spent saving: $time_spent_saving")
  println("Algorithm run time: $(timing_of_run-time_spent_saving)")

  timing_flux_sources = timing_of_precip + timing_of_bounding_flood + timing_of_dam_failure

  # we'll also save this to a JSON file so that it can easily be parsed
  # for analytical purposes
  date = Dates.now()
  format = "yyyy-mm-dd_HH-MM-SS"
  timing_filename = config["output-directory"] * "timing-for-$iteration-" * Dates.format(date, format) * ".json"
  timing_text =
    """
    {
      "configuration": $(JSON.json(config, 2)),
      "timing-flux-sources": $timing_flux_sources,
      "timing-output": $time_spent_saving,
      "timing-algorithm": $(timing_of_run - time_spent_saving),
      "timing-runtime": $(timing_flux_sources + timing_of_run)
    }
    """
  open(timing_filename, "w") do io
    write(io, timing_text)
  end

end


"""
    contamination(config::Dict{String,Any}) :: Union{Nothing,ContaminatedArea}

Parse "contaminated-area" portion of the configuration, if it exists.
"""
function contamination(config::Dict{String,Any}) :: Union{Nothing,ContaminatedArea}
  if haskey(config, "contaminated-area")
    ContaminatedArea(
      config["contaminated-area"]["ll-latitude"],
      config["contaminated-area"]["ll-longitude"],
      config["contaminated-area"]["ur-latitude"],
      config["contaminated-area"]["ur-longitude"],
      config["contaminated-area"]["start-step"],
      config["contaminated-area"]["end-step"],
      config["contaminated-area"]["contamination-rate"]
    )
  else
    nothing
  end
end