module Torrent

import Base.Threads.@threads
import Dates
import JSON
import Statistics
import DataStructures
import Printf
import Random
import ArchGDAL as AG
import LinearAlgebra
import Statistics
import DSP
import CSV
import Tables

using ProgressMeter

# main entry point
export torrent

# Grid related capabilities
export GeoRegistration
export Grid
export indexof
export latlongof
export is_outside_boundary
export is_on_boundary
export open_esri_asc_file
export save_esri_asc_file
export open_geotiff
export open_multiband_geotiff
export open_raster
export save_geotiff
export save_raster

include("helpers.jl")
include("grids.jl")
include("precipitation.jl")
include("rivulets.jl")
include("simulation_runner.jl")

end # module Torrent
