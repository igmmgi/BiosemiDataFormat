# channels.jl - Channel manipulation functions for BiosemiDataFormat

"""
    select_channels_bdf!(bdf, channels)

Select specific channels from BDF data (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `channels::Vector{<:Union{Int,String}}`: Channels to keep
  - Vector of integers: channel indices (1-based, excluding status channel)
  - Vector of strings: channel labels (e.g., ["Fp1", "Cz"])
  - Mix of types allowed (e.g., [1, "Fp1"])

# Returns
- `Nothing`: Modifies `bdf` in-place

# Examples
```julia
# Select channels by index
select_channels_bdf!(dat, [1, 3, 5])

# Select channels by label
select_channels_bdf!(dat, ["Fp1", "Cz", "A1"])

# Mix of channel types
select_channels_bdf!(dat, [1, "Fp1"])
```

# Notes
- Modifies the original data structure
- The status/trigger channel is always included automatically
- Updates header information (num_channels, channel_labels, etc.)
- Channel indices are 1-based
- Use `select_channels_bdf` for non-mutating version

# See also
- `select_channels_bdf`: Non-mutating version
- `delete_channels_bdf!`: Remove specific channels
- `read_bdf`: Read with channel selection
"""
function select_channels_bdf!(bdf::BiosemiData, channels::Union{Vector{Int}, Vector{String}})
  @info "Selecting channels: $channels"
  channels = channel_index(bdf.header.channel_labels, channels)
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end

# Convenience methods for single inputs
select_channels_bdf!(bdf_in::BiosemiData, channels::Union{Int,String}) = select_channels_bdf!(bdf_in, [channels])


"""
    select_channels_bdf(bdf_in, channels)

Select specific channels from BDF data (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to keep

# Returns
- `BiosemiData`: New data structure with selected channels

# Examples
```julia
# Select and get new structure
dat_selected = select_channels_bdf(dat, [1, 3, 5])
dat_selected = select_channels_bdf(dat, ["Fp1", "Cz"])
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `select_channels_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `select_channels_bdf!`: In-place version
- `delete_channels_bdf`: Remove specific channels
"""
function select_channels_bdf(bdf_in::BiosemiData, channels::Union{Vector{Int}, Vector{String}})
  bdf_out = deepcopy(bdf_in)
  select_channels_bdf!(bdf_out, channels)
  return bdf_out
end

# Convenience methods for single inputs
select_channels_bdf(bdf_in::BiosemiData, channels::Union{Int,String}) = select_channels_bdf(bdf_in, [channels])


"""
    delete_channels_bdf!(bdf, channels)

Remove specific channels from BDF data (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to remove

# Returns
- `Nothing`: Modifies `bdf` in-place

# Examples
```julia
# Remove channels by index
delete_channels_bdf!(dat, [2, 4])

# Remove channels by label
delete_channels_bdf!(dat, ["Fp2", "F8"])

# Mix of channel types
delete_channels_bdf!(dat, [2, "Fp2"])
```

# Notes
- Modifies the original data structure
- The status/trigger channel cannot be removed
- Updates header information (num_channels, channel_labels, etc.)
- Channel indices are 1-based
- Use `delete_channels_bdf` for non-mutating version

# See also
- `delete_channels_bdf`: Non-mutating version
- `select_channels_bdf!`: Keep specific channels
- `read_bdf`: Read with channel selection
"""
function delete_channels_bdf!(bdf::BiosemiData, channels::Union{Vector{Int}, Vector{String}})
  @info "Deleting channels: $channels"
  channels = channel_index(bdf.header.channel_labels, channels)
  channels = filter(x -> !(x in channels[1:end-1]), 1:length(bdf.header.channel_labels))
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end

# Convenience methods for single inputs
delete_channels_bdf!(bdf_in::BiosemiData, channels::Union{Int,String}) = delete_channels_bdf!(bdf_in, [channels])


"""
    delete_channels_bdf(bdf_in, channels)

Remove specific channels from BDF data (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to remove

# Returns
- `BiosemiData`: New data structure with channels removed

# Examples
```julia
# Remove channels and get new structure
dat_reduced = delete_channels_bdf(dat, [2, 4])
dat_reduced = delete_channels_bdf(dat, ["Fp2", "F8"])
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `delete_channels_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `delete_channels_bdf!`: In-place version
- `select_channels_bdf`: Keep specific channels
"""
function delete_channels_bdf(bdf_in::BiosemiData, channels::Union{Vector{Int}, Vector{String}})
  bdf_out = deepcopy(bdf_in)
  delete_channels_bdf!(bdf_out, channels)
  return bdf_out
end

# Convenience methods for single inputs
delete_channels_bdf(bdf_in::BiosemiData, channels::Union{Int,String}) = delete_channels_bdf(bdf_in, [channels])


"""
    channel_index(labels, channels)

Convert channel specifications to channel indices.

# Examples
```julia
# Single channel by label
idx = channel_index(["A1", "A2", "A3"], "A2")  # Returns [2, 3]

# Single channel by index
idx = channel_index(["A1", "A2", "A3"], 2)     # Returns [2, 3]

# Multiple channels by index
idx = channel_index(["A1", "A2", "A3"], [1, 3])  # Returns [1, 3, 3]

# Multiple channels by label
idx = channel_index(["A1", "A2", "A3"], ["A1", "A3"])  # Returns [1, 3, 3]

# Mixed types
idx = channel_index(["A1", "A2", "A3"], [1, "A2"])  # Returns [1, 2, 3]

# Trigger channel only
idx = channel_index(["A1", "A2", "A3"], [-1])  # Returns [3, 3]
```

# Arguments
- `labels::Vector{<:AbstractString}`: Available channel labels
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channel specifications

# Returns
- `Vector{Int}`: Channel indices including the status channel

# Notes
- This is an internal function used by other functions
- Channel indices are 1-based
- The status channel is always included automatically
- Returns indices in the order specified
- Useful for channel selection and deletion operations
"""
function channel_index(labels, channels::Array{String})
  channels = [findfirst(x .== labels) for x in channels]
  any(channels .=== nothing) && error("Requested channel label is not in the bdf file!")
  return sort(unique(append!(channels, length(labels))))
end

function channel_index(labels, channels::Array{Int})
  trigSelected = findall(x -> x == -1, channels)
  if length(trigSelected) > 0
    channels[trigSelected] = repeat([length(labels)], length(trigSelected))
  end
  any(channels .> length(labels)) && error("Requested channel number greater than number of channels in file!")
  any(channels .< 1) && error("Requested channel number less than 1!")
  return sort(unique(append!(channels, length(labels))))
end

# Handle SubString types and other AbstractString types
function channel_index(labels::Vector{<:AbstractString}, channels::Array{String})
  channel_idx = Int[]
  for chan in channels
    idx = findfirst(x -> x == chan, labels)
    idx === nothing && error("Channel label $chan not found!")
    push!(channel_idx, idx)
  end
  # add status channel if not already included
  if length(labels) ∉ channel_idx
    push!(channel_idx, length(labels))
  end
  return sort(unique(channel_idx))
end

function channel_index(labels::Vector{<:AbstractString}, channels::Array{Int})
  trigSelected = findall(x -> x == -1, channels)
  if length(trigSelected) > 0
    channels[trigSelected] = repeat([length(labels)], length(trigSelected))
  end
  any(channels .> length(labels)) && error("Requested channel number greater than number of channels in file!")
  any(channels .< 1) && error("Requested channel number less than 1!")
  # add status channel if not already included
  if length(labels) ∉ channels
    push!(channels, length(labels))
  end
  return sort(unique(channels))
end

# Convenience methods for single inputs
channel_index(labels, channels::String) = channel_index(labels, [channels])
channel_index(labels, channels::Int) = channel_index(labels, [channels])
channel_index(labels::Vector{<:AbstractString}, channels::String) = channel_index(labels, [channels])
channel_index(labels::Vector{<:AbstractString}, channels::Int) = channel_index(labels, [channels])


