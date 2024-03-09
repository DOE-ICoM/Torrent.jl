
"The georegistration information for a grid."
struct GeoRegistration

  "Number of columns in the raster."
  ncols::Int64

  "Number of rows in the raster."
  nrows::Int64

  "X location of the lower left of the grid."
  xll::Float64

  "Y location of the lower left of the grid."
  yll::Float64

  "Grid cell size."
  cell_size::Float64

  "Value that should be interpreted as no data."
  no_data_value::Union{Float64,Nothing}

  "The cell size of the grid in meters (versus say, degrees)."
  cell_size_meters::Float64

  "A projection string defining the CRS if the grid was created from a GeoTIFF."
  proj_string::String

  "Whether the units of xll, yll, and cell_size are in degrees or meters."
  units::String

end

"Default instance of GeoRegistration"
function GeoRegistration()
  GeoRegistration(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "")
end


"""
    GeoRegistration(ncols, nrows, xll, yll, cell_size, no_data_value, cell_size_meters)

GeoRegistration constructor where the projection string, proj_string, is assumed to be empty.
"""
function GeoRegistration(ncols::Int, nrows::Int, xll::Float64, yll::Float64, cell_size::Float64, no_data_value::Float64, cell_size_meters::Float64, units::String)
  GeoRegistration(ncols, nrows, xll, yll, cell_size, no_data_value, cell_size_meters, "", units)
end


"""
    GeoRegistration(reg::GeoRegistration; no_data_value::Union{Float64,Nothing})

Copy constructor where the no data value may be changed.
"""
function GeoRegistration(reg::GeoRegistration; no_data_value::Union{Float64,Nothing})
  GeoRegistration(reg.ncols, reg.nrows, reg.xll, reg.yll, reg.cell_size, no_data_value, reg.cell_size_meters, reg.proj_string, reg.units)
end


"A georegistered grid."
struct Grid

  "The georegistration information."
  registration::GeoRegistration

  "The actual raster data array."
  data::Array{Float32, 2}
end


"""
    cell_size_conversion_factor(units::String, yll::Float64)

Scale factor to use to convert from georegistration units to meters.
"""
function cell_size_conversion_factor(units::String, yll::Float64) :: Float64
  if lowercase(units) == "meters"
    1.0
  elseif lowercase(units) == "degrees"
    deg2rad = π / 180.0
    a = 6378.1370 # radius at equator in kilometers
    e2 = 0.0066943799901413165
    oneminus_e2 = 1.0 - 0.0066943799901413165
    deg2radX1000 = deg2rad * 1000.0
    lat_rad = yll * deg2rad
    u = 1.0 - e2 * sin(lat_rad) * sin(lat_rad)
    # latmetersperdegree = (oneminus_e2 / u) * (a / sqrt(u)) * deg2radX1000
    lonmetersperdegree = cos(lat_rad) * (a / sqrt(u)) * deg2radX1000
    lonmetersperdegree
  else
    throw(ArgumentError("Units must either be 'meters' or 'degrees'."))
  end

end


"""
  Base.getindex(g::Matrix, idx::Index)

  Overriding the getindex function allows essentially overrides
  the bracket opearator so that a grid element can be accessed as `grid[idx]`.
"""
function Base.getindex(arr::Matrix, idx::Index)
  arr[idx.row, idx.col]
end


"""
  Base.getindex(g::Grid, idx::Index)

  Overriding the getindex function overrides the bracket opearator
  so that a grid element can be accessed as `grid[Index]`.
"""
function Base.getindex(g::Grid, idx::Index)
  g.data[idx.row, idx.col]
end


"""
  Base.getindex(g::Grid, row::Int, col::Int)

  Overriding the getindex function overrides the bracket opearator
  so that a grid element can be accessed as `grid[row, col]`.
"""
function Base.getindex(g::Grid, row::Int, col::Int)
  g.data[row, col]
end



"""
    indexof(lat::Float64, lon::Float64, reg::GeoRegistration)

    Get the grid indices of a supplied latitude and longitude.
"""
function indexof(lat::Float64, lon::Float64, reg::GeoRegistration) :: Index
  col = round(Int, (lon-reg.xll) / reg.cell_size)
  row = reg.nrows - round(Int, (lat-reg.yll) / reg.cell_size)
  Index(row, col)
end


"""
    indexof(lat::Float64, lon::Float64, g::Grid)

    Get the grid indices of a supplied latitude and longitude.
"""
function indexof(lat::Float64, lon::Float64, g::Grid) :: Index
  indexof(lat, lon, g.registration)
end


"""
    indexof(latlongs::Vector{Float64}, reg::GeoRegistration)

    Get the grid indices of a latitude and longitude supplied as a Vector.
"""
function indexof(latlongs::Vector{Float64}, reg::GeoRegistration) :: Index
  lat, lon = latlongs
  indexof(lat, lon, reg)
end


"""
    indexof(latlong::Vector{Float64}, g::Grid)

    Get the grid indices of a latitude and longitude supplied as a Vector.
"""
function indexof(latlong::Vector{Float64}, g::Grid) :: Index
  lat, lon = latlongs
  indexof(lat, lon, g)
end


"""
    latlongof(row::Int, col::Int, reg::GeoRegistration)

    Get the latitude and longitude of cell.
"""
function latlongof(row::Int, col::Int, reg::GeoRegistration) :: Vector{Float64}
  lon = col * reg.cell_size + reg.xll
  lat = (reg.nrows - row) * reg.cell_size + reg.yll
  [lat, lon]
end


"""
    latlongof(idx::Index, reg::GeoRegistration)

    Get the latitude and longitude of cell.
"""
function latlongof(idx::Index, reg::GeoRegistration) :: Vector{Float64}
  latlongof(idx.row, idx.col, reg)
end


"""
    is_outside_boundary(idx::Index, reg::GeoRegistration)

    Whether the index lies outside the boundary specified by the
    nrows/ncols of the georegistration.
"""
function is_outside_boundary(idx::Index, reg::GeoRegistration)::Bool
  idx.row < 1 || idx.col < 1 || idx.row > reg.nrows || idx.col > reg.ncols
end


"""
    is_outside_boundary(idx::Index, g::Grid)

    Whether the index lies outside the boundary specified by the
    nrows/ncols of the georegistration.
"""
function is_outside_boundary(idx::Index, g::Grid)::Bool
  is_outside_boundary(idx, g.registration)
end


"""
    is_on_boundary(idx::Index, reg::GeoRegistration)::Bool

    Whether the index lies on the boundary specified by the
    nrows/ncols of the georegistration.
"""
function is_on_boundary(idx::Index, reg::GeoRegistration)::Bool
  idx.row == 1 || idx.col == 1 || idx.row == reg.nrows || idx.col == reg.ncols
end


"""
    is_on_boundary(idx::Index, g::Grid)::Bool

    Whether the index lies on the boundary specified by the
    nrows/ncols of the georegistration.
"""
function is_on_boundary(idx::Index, g::Grid)::Bool
  is_on_boundary(idx, g.registration)
end


"""
    open_esri_asc_file(filename::String, field_separator::String = "\\s+")::Grid

    Read a grid from an ESRI ASC raster file and return it as a Grid instance.
"""
function open_esri_asc_file_old(filename::String, field_separator::String = "\\s+")::Grid

  _, (registration, data) = time_computation("Loading grid from $filename ... ") do

    # open the file and grab the array of lines
    lines = open(filename) do file
      readlines(file)
    end

    # extract the georegistration from the first six lines
    reg = extract_registration(lines[1:6])

    # then extract the data from the rest
    nonheaderlines = view(lines, 7:length(lines))

    # first splitting all of the lines into string fields based on the provided separator
    substrings = ((line)->split(strip(line), Regex(field_separator))).(nonheaderlines)

    # then parsing the strings as instances of floats
    data = ((vec)->(((ss)->parse(Float32, String(ss))).(vec))).(substrings)

    # the result of the above computation is a vector of vectors which where
    # want to turn into an actual two-dimensional array instance. this explicit
    # copy approach is much faster than methods suggested that utilize reduce
    # hcat (it also gives us the correct orientation, which would otherwise need 
    # to be followed by a transposition)
    matrix = Array{Float32, 2}(undef, length(data), length(data[1]))
    for i in 1:length(data[1])
      for j in 1:length(data)
        matrix[j,i] = data[j][i]
      end
    end

    reg, matrix

  end

  Grid(registration, data)

end


"""
    open_esri_asc_file(filename::String, field_separator::String = "\\s+")::Grid

    Read a grid from an ESRI ASC raster file and return it as a Grid instance.
    This version should be a little faster and, likely, a lot less memory
    intensive than the version above.
"""
function open_esri_asc_file(filename::String, silent::Bool=true, field_separator::String = "\\s+"; units::String = "degrees")::Grid

  _, (registration, data) = time_computation("Loading grid from $(short_filename(filename)) ", silent) do

    open(filename) do file

      # first we just want to grab the first size lines and extract
      # the georegistration information from them
      reg_lines = Array{String}(undef, 6)
      count = 1
      while count <= 6 && !eof(file)
        reg_lines[count] = readline(file)
        count += 1
      end

      # extract the georegistration
      reg = extract_registration(reg_lines, units)

      # now we can extract the data from the rest of the file
      # let's create an array to hold the data
      arr = Array{Float32, 2}(undef, reg.nrows, reg.ncols)

      # then we'll stream through line by line and populate that table
      count = 1
      while !eof(file) && count <= reg.nrows
        line = readline(file)
        stripped_line = strip(line)
        if length(stripped_line) > 0
          fields = split(stripped_line, Regex(field_separator))
          for j in 1:length(fields)
            arr[count,j] = parse(Float32, fields[j])
          end
          count += 1
        end
      end

      # pass the registration and data back out of the open do-block
      # and essentially through the computation timer
      reg, arr

    end  # open file
  end  # time_computation

  Grid(registration, data)

end



"""
    extract_registration(lines::Array{String, 1}, field_separator::String = "\\s+")::GeoRegistration

    Extract the GeoRegistration information from the ESRI ASC header lines.
"""
function extract_registration(lines::Array{String, 1}, units::String, field_separator::String = "\\s+")::GeoRegistration

  # extract dictionary of registration information from header lines
  info = Dict(map(function (line)
                    pair = split(strip(line), Regex(field_separator))
                    (uppercase(pair[1]), pair[2])
                  end, lines))

  yll = parse(Float64, get(info, "YLLCENTER", info["YLLCORNER"]))

  GeoRegistration(
    parse(Int64, info["NCOLS"]),
    parse(Int64, info["NROWS"]),
    parse(Float64, get(info, "XLLCENTER", info["XLLCORNER"])),
    yll,
    parse(Float64, info["CELLSIZE"]),
    parse(Float64, info["NODATA_VALUE"]),
    parse(Float64, info["CELLSIZE"]) * cell_size_conversion_factor(units, yll),
    units
  )                

end


"""
    save_esri_asc_file(
      filename::String, 
      data,
      registration::GeoRegistration,
      interpolate_at_cell_centers::Bool
    )

    Save a raster grid to an ESRI ASC grid file.
"""
function save_esri_asc_file(
  filename::String, 
  data,
  registration::GeoRegistration,
  interpolate_at_cell_centers::Bool
)
  # we may want to interpolate the output to cell centers for comparison with
  # shallow water equation based results, such as from RIFT. if so the size
  # of the grid will change slightly as will the position of the lower left
  # corner of the raster
  decr = interpolate_at_cell_centers ? 1 : 0
  cell_offset = interpolate_at_cell_centers ? registration.cell_size/2.0 : 0.0

  # create the header lines from the registration information and potential
  # cell offsets
  header =
    """
    NCOLS $(registration.ncols - decr)
    NROWS $(registration.nrows - decr)
    XLLCORNER $(registration.xll + cell_offset)
    YLLCORNER $(registration.yll + cell_offset)
    CELLSIZE $(registration.cell_size)
    NODATA_VALUE $(registration.no_data_value)
    """

  time_computation("Saving grid to $(short_filename(filename))...") do 
    
    # perform the interpolation to cell centers if requested
    data_table = if interpolate_at_cell_centers
      arr = interpolate_at_centers(data)
    else
      data
    end

    # then save the file
    open(filename, "w") do io
      
      # writing the header lines first
      print(io, header)
      
      # then the raster table. note that we're retaining only three digits of
      # precision here which should help reduce the file size, while still being
      # sufficient for our current needs (this basically limits our depth
      # resolution to 1mm).
      for j in 1:size(data_table)[1]
        for i in 1:size(data_table)[2]-1
          print(io, Printf.@sprintf("%.3f", data_table[j,i]))
          print(io, " ")
        end
        print(io, Printf.@sprintf("%.3f", data_table[j,size(data_table)[2]]))
        print(io, "\n")
      end
    end
  end
end


"""
    interpolate_at_centers(data::Array{T, 2}) where {T <: Number}

    Interpolate the data array at cell centers resulting in a reduction of one
    cell in size along each axis. 
"""
function interpolate_at_centers(data::Array{T, 2}) where {T <: Number}
  interp = Array{T}(undef, size(data) .- 1)
  for i in 1:size(data)[2]-1
    for j in 1:size(data)[1]-1
      interp[j,i] = 0.25 * (data[j,i] + data[j,i+1] + data[j+1,i] + data[j+1,i+1])
    end
  end
  interp
end


"""
    open_geotiff(filename::String)::Grid

    Uses the ArchGDAL library to load a GeoTIFF. Extracts the first band as well
    as the georegistration information from the file and returns the pair as a
    Grid object.
"""
function open_geotiff(filename::String, silent::Bool=true; units::String="degrees")::Grid

  _, grid = time_computation("Loading GeoTIFF from $(short_filename(filename)) ", silent) do
    
    # open the geotiff and get dataset and band references
    dataset = AG.read(filename)
    band = AG.getband(dataset, 1)

    # grab registration information and translate into the form
    # expected by the GeoRegistration struct
    nrows = AG.height(dataset)
    ncols = AG.width(dataset)
    nodata = AG.getnodatavalue(band)
    xll, dx, _, yul, _, dy = AG.getgeotransform(dataset)
    dy = abs(dy)
    yll = yul - nrows * dy
    cell_size_meters = dy * cell_size_conversion_factor(units, yll)

    # read the data from the band and transpose
    arr = AG.read(band)
    trans = Array{Float32}(undef, size(arr,2), size(arr,1))
    LinearAlgebra.transpose!(trans, arr)

    # create registration and grid instance and return
    reg = GeoRegistration(ncols, nrows, xll, yll, dy, nodata, cell_size_meters, AG.getproj(dataset), units)
    Grid(reg, trans)
  end

  grid

end


"""
    open_multiband_geotiff(filename::String)::Vector{Grid}

    TBW
"""
function open_multiband_geotiff(filename::String, silent::Bool=true; units::String="degrees")::Vector{Grid}

  _, grids = time_computation("Loading multiband GeoTIFF from $(short_filename(filename)) ", silent) do

    dataset = AG.read(filename)
    num_bands = AG.nraster(dataset)
    nrows = AG.height(dataset)
    ncols = AG.width(dataset)

    # get georegistration information from the dataset
    xll, dx, _, yul, _, dy = AG.getgeotransform(dataset)
    dy = abs(dy)
    yll = yul - nrows * dy
    cell_size_meters = dy * cell_size_conversion_factor(units, yll)
    
    gridset = [begin
      band = AG.getband(dataset, i)
      nodata = AG.getnodatavalue(band)
      arr = AG.read(band)
      trans = Array{Float32}(undef, size(arr,2), size(arr,1))
      LinearAlgebra.transpose!(trans, arr)

      reg = GeoRegistration(ncols, nrows, xll, yll, dy, nodata, cell_size_meters, AG.getproj(dataset), units)
      Grid(reg, trans)
    end for i in 1:num_bands]
    gridset
  end
  grids
end


struct BandIterator
  dataset::AG.IDataset
  reg::GeoRegistration
end


"""
    length(iter::BandIterator)

Number of bands in the iterator.
"""
function Base.length(iter::BandIterator)
  AG.nraster(iter.dataset)
end


"""
    interate(iter::BandIterator) :: Union{Tuple{Grid, Int}, Nothing}

Iterate over the bands in a multiband geotiff.
"""
function interate(iter::BandIterator) :: Union{Tuple{Grid, Int}, Nothing}
  iterate(iter, 1)
end


"""
    iterate(iter::BandIterator, band_num::Int) :: Union{Tuple{Grid, Int}, Nothing}

Get the requested band from the iterator and return updated state, i.e., next band number.
"""
function iterate(iter::BandIterator, band_num::Int) :: Union{Tuple{Grid, Int}, Nothing}
  if band_num > AG.nraster(iter.dataset)
    nothing
  else
    band = AG.getband(iter.dataset, band_num)
    nodata = AG.getnodatavalue(band)
    arr = AG.read(band)
    trans = Array{Float32}(undef, size(arr,2), size(arr,1))
    LinearAlgebra.transpose!(trans, arr)
    reg = GeoRegistration(iter.reg; no_data_value = nodata)
    (Grid(reg, trans), band_num+1)
  end
end


"""
    open_multiband_stream(filename::String)::BandIterator

Open a multiband geotiff for streaming of the bands sequentially.
"""
function open_multiband_stream(filename::String; units::String = "degrees")::BandIterator

  dataset = AG.read(filename)
  nrows = AG.height(dataset)
  ncols = AG.width(dataset)

  # get georegistration information from the dataset
  xll, dx, _, yul, _, dy = AG.getgeotransform(dataset)
  dy = abs(dy)
  yll = yul - nrows * dy
  cell_size_meters = dy * cell_size_conversion_factor(units, yll)

  reg = GeoRegistration(ncols, nrows, xll, yll, dy, nothing, cell_size_meters, AG.getproj(dataset), units)
  BandIterator(dataset, reg)
end



"""
    open_raster(filename::String)::Grid

    Opens a raster file and returns a Grid instance. The type of file open_esri_asc_file_old
    depends on the filename extension. ".tif" and ".asc" are valid extensions.
"""
function open_raster(filename::String, silent::Bool = true; units::String="degrees")::Grid
  ext_idx = findlast('.', filename)
  ext = lowercase(filename[ext_idx+1:end])
  if ext == "tif" || ext == "tiff" || ext == "geotiff"
    open_geotiff(filename, silent; units=units)
  elseif ext == "asc" || ext == "txt"
    open_esri_asc_file(filename, silent; units=units)
  else
    throw(ArgumentError("Unrecognized file extension opening a raster. Expected either .tif or .asc."))
  end
end


"""
    save_geotiff(filename::String, grid::Grid)

    Save a data table as a GeoTIFF.
"""
function save_geotiff(filename::String, grid::Grid, interpolate_at_cell_centers::Bool, silent::Bool=true)
  save_geotiff(filename, grid.data, grid.registration, interpolate_at_cell_centers, silent)
end


"""
    save_geotiff(filename::String, data::Matrix{Float32}, registration::GeoRegistration)

    Save a data table as a GeoTIFF.
"""
function save_geotiff(filename::String, data::Matrix{Float32}, registration::GeoRegistration, interpolate_at_cell_centers::Bool, silent::Bool=true)

  # make sure there's a valid projection associated with the registration if
  # we're saving as a geotiff
  @assert strip(registration.proj_string) != ""

  time_computation("Saving grid to $(short_filename(filename)) ", silent) do

    # transpose the input data set back to what the julia libraries are expecting
    # with column major
    trans = Array{Float32}(undef, size(data,2), size(data,1))
    LinearAlgebra.transpose!(trans, data)

    # we may want to interpolate the output to cell centers for comparison with
    # shallow water equation based results, such as from RIFT. if so the size
    # of the grid will change slightly as will the position of the lower left
    # corner of the raster

    trans = if interpolate_at_cell_centers
      interpolate_at_centers(trans)
    else
      trans
    end

    decr = interpolate_at_cell_centers ? 1 : 0
    cell_offset = interpolate_at_cell_centers ? registration.cell_size/2.0 : 0.0

    ncols = registration.ncols - decr
    nrows = registration.nrows - decr
    xll = registration.xll + cell_offset
    yll = registration.yll + cell_offset

    # let's translate back to the upper left reference point that geotiffs like
    # and create the geotransform array
    yul = yll + nrows * registration.cell_size
    gt = [xll, registration.cell_size, 0.0, yul, 0.0, -registration.cell_size]

    AG.create(filename,
      driver = AG.getdriver("GTiff"),
      width = ncols,
      height = nrows,
      nbands = 1,
      dtype = Float32,
      options = ["COMPRESS=LZW"]
    ) do dataset
      AG.write!(dataset, trans, 1)
      AG.setgeotransform!(dataset, gt)
      AG.setproj!(dataset, registration.proj_string)
      band = AG.getband(dataset, 1)
      if !isnothing(registration.no_data_value); AG.setnodatavalue!(band, registration.no_data_value) end
    end
  end  # time_computation
end


"""
    save_geotiff(filename::String, grids::Vector{Grid}, interpolate_at_cell_centers::Bool)

Save a multiband geotiff. 
"""
function save_geotiff(
  filename::String, 
  grids::Vector{Grid};
  interpolate_at_cell_centers::Bool,
  silent::Bool=true,
  alt_proj_string::String="",
  sub_sample_by::Int = 1
)

  # make sure we didn't enter a negative number to subsample by
  @assert sub_sample_by >= 1

  # make sure we have at least one band to save
  @assert length(grids) > 0

  # make sure all the grids are the same size
  @assert length(unique(map(x -> x.registration.nrows, grids))) == 1 && length(unique(map(x -> x.registration.ncols, grids))) == 1

  # if the first grid instance has a non-empty projection string, we'll use it,
  # otherwise we'll use the alt_proj_string (which defaults to an empty string).
  proj_string = if(strip(grids[1].registration.proj_string) !="") 
    grids[1].registration.proj_string
  else
    alt_proj_string
  end

  # make sure one of the two contained a non-empty string
  @assert strip(proj_string) != ""

  # figure out the registration information
  registration = grids[1].registration
  decr = interpolate_at_cell_centers ? 1 : 0
  cell_offset = interpolate_at_cell_centers ? registration.cell_size/2.0 * sub_sample_by : 0.0
  ncols = div(registration.ncols, sub_sample_by) - decr
  nrows = div(registration.nrows, sub_sample_by) - decr
  xll = registration.xll + cell_offset
  yll = registration.yll + cell_offset

  # let's translate back to the upper left reference point that geotiffs like
  # and create the geotransform array
  yul = yll + nrows * registration.cell_size * sub_sample_by
  gt = [xll, registration.cell_size * sub_sample_by, 0.0, yul, 0.0, -registration.cell_size * sub_sample_by]

  # now we're ready to actually create the file and save the data
  time_computation("Saving multi-band grid to $(short_filename(filename)) ", silent) do

      AG.create(filename,
        driver = AG.getdriver("GTiff"),
        width = ncols,
        height = nrows,
        nbands = length(grids),
        dtype = Float32,
        options = ["COMPRESS=LZW"]
      ) do dataset

        # transpose, possibly interpolate, and then write each grid as a band to the dataset
        for i in 1:length(grids)

          # transpose the input data set back to what the julia libraries are expecting
          # with column major
          trans = Array{Float32}(undef, size(grids[i].data,2), size(grids[i].data,1))
          LinearAlgebra.transpose!(trans, grids[i].data)

          # we may want to interpolate the output to cell centers for comparison with
          # shallow water equation based results, such as from RIFT. if so, the size
          # of the grid will change slightly as will the position of the lower left
          # corner of the raster (which we already took into account when doing the
          # registration calculations above).
          trans = if interpolate_at_cell_centers
            interpolate_at_centers(trans)
          else
            trans
          end

          # sub-sample the grid if requested
          trans = if sub_sample_by > 1
            trans[1:sub_sample_by:end, 1:sub_sample_by:end]
          else
            trans
          end

          # write the transposed data to the the appropriate band
          AG.write!(dataset, trans, i)

          # set a no data value for each band
          band = AG.getband(dataset, i)
          AG.setnodatavalue!(band, registration.no_data_value)
        end

        # add projection information
        AG.setgeotransform!(dataset, gt)
        AG.setproj!(dataset, proj_string)
        
      end  # create multiband geotiff

  end  # time_computation

end


"""
    save_raster(filename::String, grid::Grid, interpolate_at_centers::Bool = false)

    Save a raster in the format specified by the filename extension.
"""
function save_raster(filename::String, grid::Grid, interpolate_at_centers::Bool = false)
  ext_idx = findlast('.', filename)
  ext = lowercase(filename[ext_idx+1:end])
  if ext == "tif" || ext == "tiff" || ext == "geotiff"
    save_geotiff(filename, grid, interpolate_at_centers)
  elseif ext == "asc" || ext == "txt"
    save_esri_asc_file(filename, grid.data, grid.registration, interpolate_at_centers)
  else
    throw(ArgumentError("Unrecognized file extension saving a raster. Expected either .tif or .asc."))
  end
end


"""
  convert(idx::Index, from_registration::GeoRegistration, to_registration::GeoRegistration) :: Index

  TBW
"""
function convert(idx::Index, from_registration::GeoRegistration, to_registration::GeoRegistration) :: Index
  indexof(latlongof(idx.row, idx.col, from_registration), to_registration)
end


"""
    rerasterize{T<:Real}(
      data::Array{T,2},
      from_registration::GeoRegistration,
      to_registration::GeoRegistration) :: Array{T,2}

Currently assumes that the 'from' raster is bigger than the 'to' raster.
"""
function rerasterize(
  data::Array{Float32,2},
  from_registration::GeoRegistration,
  to_registration::GeoRegistration) :: Array{Float32,2}

  new_array = zeros(Float32, to_registration.nrows, to_registration.ncols)
  for j in 1:to_registration.nrows
    for i in 1:to_registration.ncols
      idx = convert(Index(j,i), to_registration, from_registration)
      new_array[j,i] = data[idx.row, idx.col]
    end
  end
  new_array
end


"""
    smooth(grid::Grid, radius::Int) :: Grid

    TBW
"""
function smooth(grid::Grid, radius::Int) :: Grid
  kern = fill(Float32(1.0/((2*radius+1.0)^2)), 2*radius+1, 2*radius+1)
  smoothed = DSP.conv(grid.data, kern)[2*radius+1:grid.registration.nrows-2*radius-1, 2*radius+1:grid.registration.ncols-2*radius-1]
  Grid(grid.registration, smoothed)
end


"""
    smooth(raster_filename::String, radius::Int, output_filename::String)

    TBW
"""
function smooth(raster_filename::String, radius::Int, output_filename::String)
  grid = open_raster(raster_filename)
  smoove = smooth(grid, radius)
  save_raster(output_filename, smoove, false)
end



# function fractal_dimension(grid::Grid, scales::Vector{Int}) :: Grid

#   sqpi = sqrt(2.0 / pi)
#   height = grid.registration.nrows
#   weight = grid.registration.ncols
#   s0 = scales[1]
#   s1 = scales[2]
#   ls0 = log(s0 * grid.registration.cell_size_meters)
#   ls1 = log(s1 * grid.registration.cell_size_meters)
#   dls = ls1 - ls0

#   fractal_dimension_grid = zeros(Float32, height, width)
#   vertical_scale_grid = zeros(Float32, height, width)

#   for j in 1:height
#     for i in 1:width

#       mean_diffs = map( scale -> begin
#           other_points = filter(x -> x[3] == true, fractal_measurement_points(grid, Index(j,i), scale))
#           elevation_differences = map( x -> abs(grid[x[1],x[2]] - grid[j,i]), other_points)
#           mean_elevation_difference = sum(elevation_differences) / length(other_points)
#           log(mean_elevation_difference)
#         end,
#         scales
#       )

      

#     end
#   end


#   fractal_dimension_grid
# end


# """
#     fractal_measurement_points(
#       grid::Grid, 
#       central_point::Index, 
#       scale::Int, 
#       num::Int = 8)::Vector{Tuple{Int, Int, Bool}}

#     TBW
# """
# function fractal_measurement_points(
#   grid::Grid, 
#   central_point::Index, 
#   scale::Int, 
#   num::Int = 8)::Vector{Tuple{Int, Int, Bool}}

#   height = grid.registration.nrows
#   width = grid.registration.ncols

#   list = [begin
#     theta = 2.0 * pi * Float64(i) / num
#     [round(Int, scale * sin(theta)) + central_point.row,
#       round(Int, scale * cose(theta)) + central_point.col]
#   end for i in 0:num-1]

#   map( x -> (x[1], x[2], x[1]>=1 && x[2]>=0 && x[1]<=height && x[2]<=width), list)

# end