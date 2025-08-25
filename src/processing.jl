# processing.jl - Data processing functions for BiosemiDataFormat

"""
    crop_bdf!(bdf, crop_type, val)

Reduce the length of BDF data by cropping (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `crop_type::String`: Cropping method
  - `"records"`: Crop by data record numbers
  - `"triggers"`: Crop between trigger events
- `val::Vector{Int}`: Cropping parameters
  - For `"records"`: start_record to end_record (1-based)
  - For `"triggers"`: start_trigger to end_trigger

# Returns
- `Nothing`: Modifies `bdf` in-place

# Examples
```julia
# Crop between records 10-20
crop_bdf!(dat, "records", [10 20])

# Crop between first occurrence of trigger 1 and last of trigger 2
crop_bdf!(dat, "triggers", [1 2])
```

# Notes
- Modifies the original data structure
- Updates header information (num_data_records)
- Recalculates time vector and trigger information
- Records are 1-based indexing
- Use `crop_bdf` for non-mutating version

# See also
- `crop_bdf`: Non-mutating version
- `downsample_bdf!`: Reduce sampling rate
- `merge_bdf`: Combine multiple files
"""
function crop_bdf!(bdf::BiosemiData, crop_type::String, val::Array{Int})

  @info "Cropping data: $crop_type, $val"
  crop_check(crop_type, val, bdf.triggers.count.keys, bdf.header.num_data_records)
  sample_rate = bdf.header.sample_rate[1]
  nsamples = size(bdf.data, 1)
  idxStart, idxEnd = find_crop_index(bdf.triggers, crop_type, val, sample_rate, nsamples)

  # crop
  bdf.header.num_data_records = Int(((idxEnd - idxStart) + 1) / sample_rate)
  bdf.data = bdf.data[idxStart:idxEnd, :]
  bdf.time = time_range(bdf.header.sample_rate[1], bdf.header.num_data_records)
  bdf.status = bdf.status[idxStart:idxEnd]

  # recaculate trigger information
  trig = bdf.triggers.raw[idxStart:idxEnd]
  bdf.triggers = trigger_info(trig, bdf.header.sample_rate[1])

end


"""
    crop_bdf(bdf_in, crop_type, val)

Reduce the length of BDF data by cropping (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `crop_type::String`: Cropping method ("records" or "triggers")
- `val::Vector{Int}`: Cropping parameters

# Returns
- `BiosemiData`: New cropped data structure

# Examples
```julia
# Crop and get new structure
dat_cropped = crop_bdf(dat, "records", [10 20])
dat_cropped = crop_bdf(dat, "triggers", [1 2])
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `crop_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `crop_bdf!`: In-place version
- `downsample_bdf`: Reduce sampling rate
"""
function crop_bdf(bdf_in::BiosemiData, crop_type::String, val::Array{Int})
  bdf_out = deepcopy(bdf_in)
  crop_bdf!(bdf_out, crop_type, val)
  return bdf_out
end


"""
    downsample_bdf!(bdf, dec)

Reduce the sampling rate of BDF data by an integer factor (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `dec::Int`: Downsampling factor (must be power of 2)

# Returns
- `Nothing`: Modifies `bdf` in-place

# Downsampling Process
- Applies anti-aliasing filter using DSP.resample
- Reduces data length by factor `dec`
- Updates header information (sample_rate, num_samples)
- Recalculates time vector
- Adjusts trigger indices to new sampling rate

# Examples
```julia
# Downsample by factor of 2 (256 Hz → 128 Hz)
downsample_bdf!(dat, 2)

# Downsample by factor of 4 (2048 Hz → 512 Hz)
downsample_bdf!(dat, 4)
```

# Notes
- `dec` must be a power of 2 (2, 4, 8, 16, ...)
- Uses padding to handle edge effects
- Trigger timing is automatically adjusted
- Modifies the original data structure
- Use `downsample_bdf` for non-mutating version

# See also
- `downsample_bdf`: Non-mutating version
- `crop_bdf!`: Reduce data length
- `merge_bdf`: Combine multiple files
"""
function downsample_bdf!(bdf::BiosemiData, dec::Int)

  !ispow2(dec) && error("dec should be power of 2!")
  @info "Downsampling data by factor: $dec"

  # padding at start/end
  nsamp = dec * 20  # enough samples?
  ndec = div(nsamp, dec)

  data = Matrix{Float32}(undef, div(size(bdf.data, 1), dec), size(bdf.data, 2))
  for i in 1:size(bdf.data, 2)  
    start_pad = reverse(bdf.data[1:nsamp, i])
    end_pad = reverse(bdf.data[end-(nsamp-1):end, i])
    padded_data = vcat(start_pad, bdf.data[:, i], end_pad)
    tmp_dat = resample(padded_data, 1 / dec)
    data[:, i] = convert(Vector{Float32}, tmp_dat[ndec+1:end-ndec])
  end
  bdf.data = data

  bdf.header.sample_rate = div.(bdf.header.sample_rate, dec)
  bdf.header.num_samples = div.(bdf.header.num_samples, dec)
  bdf.time = (0:size(bdf.data, 1)-1) / bdf.header.sample_rate[1]

  # update triggers
  bdf.triggers.raw = zeros(Int16, size(bdf.data, 1))  
  bdf.triggers.idx = convert(Vector{Int64}, round.(bdf.triggers.idx / dec))
  valid_indices = filter(x -> 1 <= x <= size(bdf.data, 1), bdf.triggers.idx)
  if length(valid_indices) > 0
    bdf.triggers.raw[valid_indices] = bdf.triggers.val[1:length(valid_indices)]
  end

end


"""
    downsample_bdf(bdf_in, dec)

Reduce the sampling rate of BDF data by an integer factor (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `dec::Int`: Downsampling factor (must be power of 2)

# Returns
- `BiosemiData`: New downsampled data structure

# Examples
```julia
# Downsample and get new structure
dat_ds = downsample_bdf(dat, 2)
dat_ds = downsample_bdf(dat, 4)
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `downsample_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `downsample_bdf!`: In-place version
- `crop_bdf`: Reduce data length
"""
function downsample_bdf(bdf_in::BiosemiData, dec::Int)
  bdf_out = deepcopy(bdf_in)
  downsample_bdf!(bdf_out, dec)
  return bdf_out
end


"""
    merge_bdf(bdfs)

Merge multiple BioSemi BDF data structures into a single structure.

# Arguments
- `bdfs::Vector{BiosemiData}`: Vector of BDF data structures to merge

# Returns
- `BiosemiData`: Merged data structure

# Requirements
All input structures must have:
- Same number of channels
- Same channel labels
- Same sampling rate

# Merging Process
- Data records are concatenated in order
- Trigger information is recomputed from merged data
- Time vector is recalculated for the full duration
- Header information is updated with new record count

# Examples
```julia
# Merge two files
dat1 = read_bdf("file1.bdf")
dat2 = read_bdf("file2.bdf")
merged = merge_bdf([dat1, dat2])

# Merge multiple files
files = ["file1.bdf", "file2.bdf", "file3.bdf"]
data = [read_bdf(f) for f in files]
merged = merge_bdf(data)
```

# Notes
- Files are merged in the order they appear in the vector
- Total duration = sum of individual file durations
- Trigger counts are recalculated from the merged data
- Original data structures are not modified
- Useful for combining recordings from the same session

# See also
- `read_bdf`: Read individual BDF files
- `write_bdf`: Write merged data to file
- `crop_bdf`: Reduce length of merged data if needed
"""
function merge_bdf(bdfs::Array{BiosemiData})

  file_names = join([x.filename for x in bdfs], ", ")
  @info "Merging files: $file_names"

  # check data structs to merge have same number of channels, channel labels + sample rate
  num_chans = (x -> x.header.num_channels).(bdfs)
  !all(x -> x == num_chans[1], num_chans) && error("1+ files have different number of channels!")
  chan_labels = (x -> x.header.channel_labels).(bdfs)
  !all(y -> y == chan_labels[1], chan_labels) && error("1+ files have different channel labels!")
  sample_rate = (x -> x.header.sample_rate).(bdfs)
  !all(y -> y == sample_rate[1], sample_rate) && error("1+ files have different sample rates!")

  # make copy so that bdf_in is not altered
  bdf_out = deepcopy(bdfs[1])

  # merge data
  bdf_out.header.num_data_records = sum((x -> x.header.num_data_records).(bdfs))
  bdf_out.data = vcat((x -> x.data).(bdfs)...)
  bdf_out.status = vcat((x -> x.status).(bdfs)...)

  # recaculate trigger information
  trig = vcat((x -> x.triggers.raw).(bdfs)...)
  bdf_out.triggers = trigger_info(trig, bdf_out.header.sample_rate[1])

  # merged time 
  bdf_out.time = time_range(bdf_out.header.sample_rate[1], bdf_out.header.num_data_records)

  return bdf_out

end


"""
    time_range(sample_rate, num_data_records)

Generate time vector for BDF data.

# Arguments
- `sample_rate::Int`: Sampling rate in Hz
- `num_data_records::Int`: Number of data records

# Returns
- `StepRangeLen{Float64}`: Time vector starting at 0

# Time Calculation
- Time starts at 0 seconds
- Increments by 1/sample_rate for each sample
- Ends at (num_data_records - 1/sample_rate) seconds
- Total duration = num_data_records seconds

# Examples
```julia
# 60 seconds of data at 256 Hz
time = time_range(256, 60)
@assert length(time) == 256 * 60
@assert first(time) == 0.0
@assert last(time) == 59.99609375  # (60 - 1/256)
```

# Notes
- This is an internal function used by other functions
- Time vector length matches the number of data samples
- Useful for plotting and time-based operations
"""
function time_range(sample_rate::Int, num_data_records::Int)
  return 0:1/sample_rate:(num_data_records-(1/sample_rate))
end


"""
    trigger_info(trig_raw, sample_rate)

Extract trigger event information from raw trigger channel data.

# Arguments
- `trig_raw::Vector{Int16}`: Raw trigger values for each sample
- `sample_rate::Int`: Sampling rate in Hz

# Returns
- `BiosemiTriggers`: Structured trigger information containing:
  - `raw`: Original trigger values
  - `idx`: Sample indices where triggers occur
  - `val`: Trigger values at trigger events
  - `count`: Count of each trigger value
  - `time`: Trigger timing matrix [value, time_since_previous]

# Algorithm
1. Finds sample indices where trigger values change (diff ≥ 1)
2. Extracts trigger values at those indices
3. Calculates time intervals between consecutive triggers
4. Counts occurrences of each trigger value
5. Creates timing matrix with trigger values and intervals

# Notes
- This is an internal function used by `read_bdf`
- Trigger events are detected when the trigger value increases
- Time intervals are calculated in seconds
- The first trigger has time interval 0
"""
function trigger_info(trig_raw, sample_rate)

  # trigger events + time
  trig_idx = findall(diff(trig_raw) .>= 1) .+ 1
  trig_val = trig_raw[trig_idx]
  trig_time = hcat(trig_val, pushfirst!(diff(trig_idx), 0) / sample_rate)

  # unique triggers + count
  trig_count = OrderedDict(i => 0 for i in sort!(collect(Set(trig_val))))
  for val in trig_val
    trig_count[val] += 1
  end

  triggers = BiosemiTriggers(trig_raw, trig_idx, trig_val, trig_count, trig_time)

  return triggers

end


# Internal helper functions for cropping
function crop_check(crop_type::String, val::Array{Int}, triggers, num_data_records)
  crop_type ∉ ["triggers", "records"] && error("crop_type not recognized!")
  length(val) != 2 && error("val should be of length 2!")
  if crop_type == "triggers"
    val[1] ∉ triggers && error("val[1] not available trigger!")
    val[2] ∉ triggers && error("val[2] not available trigger!")
  elseif crop_type == "records"
    val[1] < 1 && error("val[1] less than 1")
    val[2] > num_data_records && error("val[2] > number of data records!")
  end
end

function find_crop_index(triggers::BiosemiTriggers, crop_type::String, val::Array{Int}, sample_rate::Int, data_length::Int)
  # find idxStart/idxEnd 
  if crop_type == "triggers"
    borders = 1:sample_rate:data_length
    trigStart = findfirst(triggers.val .== val[1])
    trigEnd = findlast(triggers.val .== val[2])
    idxStart = triggers.idx[trigStart]
    idxStart = findfirst(borders .>= idxStart) * sample_rate
    idxEnd = triggers.idx[trigEnd]
    idxEnd = (findlast(borders .<= idxEnd) * sample_rate) - 1
  elseif crop_type == "records"
    idxStart = ((val[1] - 1) * sample_rate) + 1
    idxEnd = (val[2] * sample_rate)
  end
  return idxStart, idxEnd
end
