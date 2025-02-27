
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


module Torrent

import Base.Threads.@threads
import Base.iterate
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
import Combinatorics

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
export crop_geotiff

include("helpers.jl")
include("grids.jl")
include("precipitation.jl")
include("rivulets.jl")
include("simulation_runner.jl")

end # module Torrent
