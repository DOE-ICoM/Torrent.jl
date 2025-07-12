
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



"The 2D index of a grid location."
struct Index
  row::Int64
  col::Int64
end


"A file pattern and range that can be used to read in a series of files."
struct FileRange
  filename_pattern :: String
  range :: AbstractRange{Int64}
end


"""
    array_update(f::Function, arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)

    Wrapper function for a basic array update that can be replaced by a thread-safe
    version that avoids race conditions if necessary.
"""
function array_update(f::Function, arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)
  x = arr[idx.row, idx.col]
  arr[idx.row, idx.col] = f(x)
end


"""
    safe_array_update(f::Function, arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)

    The thread-safe version of the array update function. In my initial experience, this
    seemed to introduce a LOT of overhead (like, tripling the algorithm run time). It's
    not clear that race conditions are frequent or catastrophic for the current application
    either. As such, we'll have to do further investigation about the trade off between
    computational complexity and results.
"""
function safe_array_update(f::Function, arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)
  lock(lk)
  try
    x = arr[idx.row, idx.col]
    arr[idx.row, idx.col] = f(x)
  finally
    unlock(lk)
  end
end


"""
    array_read(arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)

    Wrapper function for a basic array read that can be replaced by a thread-safe
    version that avoids race conditions if necessary.
"""
function array_read(arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)
  arr[idx.row, idx.col]
end


"""
    safe_array_read(arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)

    The thread-safe version of the array read function. In my initial experience, this
    seemed to introduce a LOT of overhead (like, tripling the algorithm run time). It's
    not clear that race conditions are frequent or catastrophic for the current application
    either. As such, we'll have to do further investigation about the trade off between
    computational complexity and results.
"""
function safe_array_read(arr::Array{Float32, 2}, idx::Index, lk::ReentrantLock)
  lock(lk)
  try
    arr[idx.row, idx.col]
  finally
    unlock(lk)
  end
end


"""
    array_read(arr::Array{Float32, 2}, row::Int, col::Int, lk::ReentrantLock)

    Wrapper function for a basic array read that can be replaced by a thread-safe
    version that avoids race conditions if necessary.
"""
function array_read(arr::Array{Float32, 2}, row::Int, col::Int, lk::ReentrantLock)
  arr[row, col]
end


"""
    safe_array_read(arr::Array{Float32, 2}, row::Int, col::Int, lk::ReentrantLock)

    The thread-safe version of the array read function. In my initial experience, this
    seemed to introduce a LOT of overhead (like, tripling the algorithm run time). It's
    not clear that race conditions are frequent or catastrophic for the current application
    either. As such, we'll have to do further investigation about the trade off between
    computational complexity and results.
"""
function safe_array_read(arr::Array{Float32, 2}, row::Int, col::Int, lk::ReentrantLock)
  lock(lk)
  try
    arr[row, col]
  finally
    unlock(lk)
  end
end


"""
    tuplefrom(idx::Index)

    Create a tuple from an Index. This allows one to assign both variables at once
    in a single line of Julia code, e.g.: j, i = tuplefrom(idx)
"""
function tuplefrom(idx::Index)
  (idx.row, idx.col)
end


"""
    time_computation(f::Function, message::String)

Times the supplied computation and returns a Tuple of the computation time and
result. Normally, you'll want to use this with Julia's do-block syntax. 
"""
function time_computation(f::Function, message::String, silent::Bool=true)
  if !silent
    printstyled(message; color=39)
  end
  start = time()
  result = f()
  dt = 1000.0 * (time() - start)
  if !silent
    str = Printf.format(Printf.Format(" [%.0f ms]\r"), dt)
    printstyled(str; color=39)
  end
  (dt, result)
end


"""
    save_csv(filename::String, data, header=nothing)

Save a CSV file with a header row.
"""
function save_csv(filename::String, data, header::Union{Nothing,Array{String}}=nothing)
  
  (timing, _) = time_computation("Saving CSV to $(short_filename(filename))... ") do

    # if input data is a vector of vectors, let's convert to a matrix first
    # may be a better way to do this, including just having two separate implementations,
    # but for the moment we'll keep things simple, though it involves so inefficient copying
    matrix = if isa(data, Vector)  # unfortunately a more comprehensive check such as Vector{Vector{<:Number}} doesn't seem to work
      reduce(hcat, data)'
    elseif isa(data, Array{<:Number, 2})
      data
    else
      throw(ArgumentError("While saving CSV, data must be either Vector{Vector{T}} or Array{T,2}."))
    end

    # open the file for writing
    open(filename, "w") do io

      # if an array of column names has been specified in the header
      # we'll write those out first
      if !isnothing(header)
        for i in 1:length(header)-1
          print(io, header[i] * ",")
        end
        print(io, header[length(header)])
        print(io, "\n")
      end

      # then write out the actual table. note that we could probably do
      # this whole thing more "functionally", but it would almost certainly
      # involve array copies rather than simply and efficiently streaming
      # out the table to disk.
      for j in 1:size(matrix)[1]
        for i in 1:size(matrix)[2]-1
          print(io, matrix[j,i])
          print(io, ",")
        end
        print(io, matrix[j,size(matrix)[2]])
        print(io, "\n")
      end
    end  # open
  end  # time_computation
  timing
end


"""
    open_csv(
      filename::String,
      field_separator::String = ",",
      skip_first_line::Boolean = false,
      skip_blank_lines::Boolean = true
    )

  Open a CSV file.
"""
function open_csv(
  field_type::Type{T},
  filename::String,
  skip_first_line::Bool = false,
  default_value::T = zero(T),
  field_separator::String = ","
)::Array{T, 2} where {T}

  (timing, result) = time_computation("Opening CSV file, $(short_filename(filename)), ... ") do 

    # use a linked list to collect parsed rows initially since
    # we won't know the file length from the outset
    table = DataStructures.MutableLinkedList{Array{T}}()
    line = ""
    count = 0

    open(filename, "r") do io

      # if skipping first line just read it in and throw away
      if skip_first_line && !eof(io)
        readline(io)
      end

      # then read the rest of the file
      while !eof(io)
        line = strip(readline(io))
        if length(line) != 0
          fields = split(line, field_separator)
          typed_fields = if length(fields) == 1 && fields[1] == ""
              T[]
            else
              map(fields) do s
                try
                  parse(T, s)
                catch e
                  if isa(e, ArgumentError)
                    default_value
                  end
                end
              end
            end  # if
          DataStructures.append!(table, typed_fields)
        end
      end  # while
    end  # open
    
    # want to return an Array{T,2} so have to copy the
    # linked list over to an array
    vector = collect(table)
    arr = Array{T, 2}(undef, length(vector), length(first(vector)))
    for j in 1:length(vector)
      for i in 1:length(vector[1])
        arr[j, i] = vector[j][i]
      end
    end
    arr
  end  # time_computation

  result
end




"""
    max_in_place(a::Array{T,2}, b::Array{T,2}) where {T <: Number}

    Element by element maximum where the matrix a is updated with the result.
    This is used to update the peak flood depth grid every so often.
"""
function max_in_place!(a::Array{T,2}, b::Array{T,2}) where {T <: Number}
  @assert size(a) == size(b)
  @threads for i in 1:size(a)[2]
    for j in 1:size(a)[1]
      a[j,i] = max(a[j,i], b[j,i])
    end
  end
end


"""
    short_filename(filename::String)::String

    Return a simple version of the filename, minus all the path information
"""
function short_filename(filename::String)::String
  idx = findlast('/', filename)
  filename[idx+1:end]
end


"""
    slicematrix(A::AbstractMatrix)

    Convert a matrix into an array of arrays.
"""
function slicematrix(A::AbstractMatrix)
  return [A[i, :] for i in 1:size(A,1)]
end


"""
    masked!(data::Matrix{Float32},
  mask::Matrix{Float32};
  threshold::Real=0.0,
  replacement_value::Real=0.0
)

TBW
"""
function masked!(
  data::Matrix{Float32},
  mask::Matrix{Float32};
  threshold::Real=0.0,
  replacement_value::Real=0.0
)
  for i in 1:size(data,2)
    for j in 1:size(data,1)
      data[j,i] = mask[j,i]>threshold ? data[j,i] : replacement_value
    end
  end
  data
end



"""
    raster_math(
      data_a::Matrix{Float32},
      data_b::Matrix{Float32},
      f::Function
    )

TBW
"""
function raster_math(
  data_a::Matrix{Float32},
  data_b::Matrix{Float32},
  f::Function
)
  result = similar(data_a)
  for i in 1:size(data_a,2)
    for j in 1:size(data_a,1)
      result[j,i] = f(data_a[j,i], data_b[j,i])
    end
  end
  result 
end


"""
    parse_distribution(x)

Parses a parameter value that may represent a distribution. If a
`Float64` value is passed the value is assumed to represent a delta-
function distribution centered at that value (i.e. to be deterministic).
If a `Dict` is passed it is assumed to have to elements, either `mean`
and `std` (in which case a normal distribution is returned); or `lower`
and `upper` (in which case a uniform distribution is returned).
"""
function parse_distribution(x)
  if typeof(x) == Float64
    x
  elseif typeof(x) == Dict{String,Any}
    if haskey(x, "mean")
      Distributions.Normal(x["mean"], x["std"])
    elseif haskey(x, "lower")
      Distributions.Uniform(x["lower"], x["upper"])
    else
      error("Distribution specification must include either mean/std or lower/upper bounds.")
    end
  else
    error("Unexpected type parsing distribution-related parameter value.")
  end
end


"""
    rand_value(x)

If `dist` is a distribution, returns a random value selected from that
distribution. If `dist` is a `Float64`, the number is simply returned.
"""
function rand_value(dist)
  if typeof(dist) == Float64
    dist
  else
    Distributions.rand(dist,1)[1]
  end
end