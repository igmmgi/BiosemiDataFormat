"""
    BiosemiDataFormat

Julia package for reading, writing, and processing BioSemi 24-bit EEG data files (BDF format).

This package provides functionality to:
- Read BDF files into Julia data structures
- Write Julia data structures back to BDF format
- Crop data by time or trigger events
- Downsample data by integer factors
- Select or reduce the number of channels
- Merge multiple BDF files
- Process trigger and status channel information

# File Format
BioSemi BDF files store 24-bit EEG data with metadata including channel information,
sampling rates, and trigger events. See the [BioSemi BDF specification](https://www.biosemi.com/faq_file_format.htm)
for detailed format information.

# Quick Start
```julia
using BiosemiDataFormat

# Read a BDF file
dat = read_bdf("data.bdf")

# Read only header information
hdr = read_bdf("data.bdf", header_only=true)

# Select specific channels
dat = read_bdf("data.bdf", channels=[1, 3, 5])

# Crop data between records 10-20
cropped = crop_bdf(dat, "records", [10 20])

# Downsample by factor of 2
downsampled = downsample_bdf(dat, 2)

# Write processed data
write_bdf(downsampled, "processed.bdf")
```

# Main Functions
- `read_bdf`: Read BDF files with optional channel selection
- `write_bdf`: Write data structures to BDF format
- `crop_bdf`: Reduce data length by time or triggers
- `downsample_bdf`: Reduce sampling rate by integer factors
- `select_channels_bdf`: Keep only specified channels
- `delete_channels_bdf`: Remove specified channels
- `merge_bdf`: Combine multiple BDF files
- `channel_index`: Convert channel labels to indices

# Data Structures
- `BiosemiHeader`: File header information
- `BiosemiData`: Complete data structure with header, data, and metadata
- `BiosemiTriggers`: Trigger event information

# Examples
See the test files for comprehensive usage examples.
"""
module BiosemiDataFormat

using DSP
using Logging
using OrderedCollections

# BDF file format constants
const BDF_HEADER_SIZE = 256
const BDF_ID_BYTES = 8
const BDF_TEXT_BYTES = 80
const BDF_CHANNEL_LABEL_BYTES = 16
const BDF_TRANSDUCER_BYTES = 80
const BDF_UNIT_BYTES = 8
const BDF_VALUE_BYTES = 8
const BDF_FILTER_BYTES = 80
const BDF_SAMPLES_BYTES = 8
const BDF_RESERVED_BYTES = 32
const BDF_DATA_FORMAT_BYTES = 44
const BDF_DATE_BYTES = 8
const BDF_TIME_BYTES = 8
const BDF_CHANNEL_COUNT_BYTES = 4
const BDF_RECORD_COUNT_BYTES = 8
const BDF_DURATION_BYTES = 8
const BDF_HEADER_SIZE_BYTES = 8

# BDF data format constants
const BDF_SAMPLES_PER_BYTE = 3  # 24-bit samples stored as 3 bytes
const BDF_STATUS_CHANNEL_OFFSET = 1  # Status channel is always last

export
  crop_bdf!,
  crop_bdf,
  delete_channels_bdf!,
  delete_channels_bdf,
  downsample_bdf!,
  downsample_bdf,
  merge_bdf,
  read_bdf,
  select_channels_bdf!,
  select_channels_bdf,
  write_bdf

"""
    BiosemiHeader

Data structure containing BioSemi BDF file header information.

# Fields
- `id1::Vector{UInt8}`: File identifier (first byte)
- `id2::Vector{UInt8}`: File identifier (remaining bytes)
- `text1::String`: Subject information
- `text2::String`: Recording information
- `start_date::String`: Recording start date (DD.MM.YYYY)
- `start_time::String`: Recording start time (HH.MM.SS)
- `num_bytes_header::Int`: Header size in bytes
- `data_format::String`: Data format (typically "24BIT")
- `num_data_records::Int`: Number of data records
- `duration_data_records::Int`: Duration of each record in seconds
- `num_channels::Int`: Number of channels (including status channel)
- `channel_labels::Vector{String}`: Channel names/labels
- `transducer_type::Vector{String}`: Transducer type for each channel
- `channel_unit::Vector{String}`: Physical units for each channel
- `physical_min::Vector{Int}`: Physical minimum values
- `physical_max::Vector{Int}`: Physical maximum values
- `digital_min::Vector{Int}`: Digital minimum values
- `digital_max::Vector{Int}`: Digital maximum values
- `pre_filter::Vector{String}`: Pre-filtering information
- `num_samples::Vector{Int}`: Number of samples per record per channel
- `reserved::Vector{String}`: Reserved header space
- `scale_factor::Vector{Float32}`: Scale factors for data conversion
- `sample_rate::Vector{Int}`: Sampling rate for each channel

# Notes
- The last channel is always the status/trigger channel
- `scale_factor` is automatically calculated from physical and digital ranges
- `sample_rate` is calculated from `num_samples` and `duration_data_records`
"""
mutable struct BiosemiHeader
  id1::Vector{UInt8}
  id2::Vector{UInt8} 
  text1::String
  text2::String
  start_date::String
  start_time::String
  num_bytes_header::Int
  data_format::String
  num_data_records::Int
  duration_data_records::Int
  num_channels::Int
  channel_labels::Vector{String}
  transducer_type::Vector{String}
  channel_unit::Vector{String}
  physical_min::Vector{Int}
  physical_max::Vector{Int}
  digital_min::Vector{Int}
  digital_max::Vector{Int}
  pre_filter::Vector{String}
  num_samples::Vector{Int}
  reserved::Vector{String}
  scale_factor::Vector{Float32}
  sample_rate::Vector{Int}
end

"""
    BiosemiTriggers

Data structure containing trigger and status channel information.

# Fields
- `raw::Vector{Int16}`: Raw trigger values for each sample
- `idx::Vector{Int}`: Sample indices where triggers occur
- `val::Vector{Int}`: Trigger values at trigger events
- `count::OrderedDict{Int,Int}`: Count of each trigger value
- `time::Matrix{Float64}`: Trigger timing information (value × time)

# Notes
- `raw` contains trigger values for every sample (0 for no trigger)
- `idx` contains sample indices where trigger values change
- `time` matrix has columns: [trigger_value, time_since_previous_trigger]
"""
mutable struct BiosemiTriggers
  raw
  idx
  val
  count
  time
end

"""
    BiosemiData

Complete data structure containing BioSemi BDF file data and metadata.

# Fields
- `filename::String`: Source filename
- `header::BiosemiHeader`: File header information
- `data::Matrix{Float32}`: EEG data matrix (samples × channels)
- `time::StepRangeLen{Float64}`: Time vector for each sample
- `triggers::BiosemiTriggers`: Trigger event information
- `status::Vector{Int16}`: Status channel values for each sample

# Notes
- `data` matrix excludes the status channel (stored separately in `status`)
- `time` vector starts at 0 and increments by 1/sample_rate
- `status` contains the status channel values (typically 0)
"""
mutable struct BiosemiData
  filename
  header
  data
  time
  triggers
  status
end

function Base.show(io::IO, hdr::BiosemiHeader)
  println(io, "Number of Channels: ", length(hdr.channel_labels) - 1)
  println(io, "Channel Labels: ", join(hdr.channel_labels[1:end-1], ", "))
  println(io, "Number of Data Records: ", hdr.num_data_records)
  println(io, "Sample Rate: ", hdr.sample_rate[1])
end

function Base.show(io::IO, trig::BiosemiTriggers)
  println(io, "Triggers (Value => Count): ", join(trig.count, ", "))
end

function Base.show(io::IO, dat::BiosemiData)
  println(io, "Filename: $(dat.filename)")
  print(io, dat.header)
  println(io, "Data Size: $(size(dat.data))")
  println(io, dat.triggers)
end


"""
    read_bdf(filename; header_only=false, channels=[])

Read BioSemi Data Format (BDF) files into Julia data structures.

# Arguments
- `filename::String`: Path to the BDF file
- `header_only::Bool=false`: If `true`, only read header information and return `BiosemiHeader`
- `channels::Vector{<:Union{Int,String}}=[]`: Specific channels to read
  - Empty vector (default): read all channels
  - Vector of integers: channel indices (1-based, excluding status channel)
  - Vector of strings: channel labels (e.g., ["Fp1", "Cz"])
  - Use `-1` to include only the trigger/status channel
  - Mix of types allowed (e.g., [1, "Fp1", -1])

# Returns
- `BiosemiData`: Complete data structure with header, data, time, triggers, and status
- `BiosemiHeader`: If `header_only=true`

# Data Structure
The returned `BiosemiData` contains:
- `header`: File metadata and channel information
- `data`: EEG data matrix (samples × selected_channels)
- `time`: Time vector starting at 0
- `triggers`: Trigger event information
- `status`: Status channel values

# Examples
```julia
# Read entire file
dat = read_bdf("data.bdf")

# Read only header
hdr = read_bdf("data.bdf", header_only=true)

# Read specific channels by index
dat = read_bdf("data.bdf", channels=[1, 3, 5])

# Read specific channels by label
dat = read_bdf("data.bdf", channels=["Fp1", "Cz", "A1"])

# Read only trigger channel
dat = read_bdf("data.bdf", channels=[-1])

# Mix of channel types
dat = read_bdf("data.bdf", channels=[1, "Fp1", -1])
```

# Notes
- Channel indices are 1-based
- The status/trigger channel is always included automatically
- Data is automatically scaled using header calibration information
- Trigger information is extracted from the status channel
- File format follows [BioSemi BDF specification](https://www.biosemi.com/faq_file_format.htm)

# See also
- `write_bdf`: Write data back to BDF format
- `crop_bdf`: Reduce data length
- `select_channels_bdf`: Select channels after reading
"""
function read_bdf(filename::String; header_only::Bool=false, channels::Union{Array{Any},Array{Int},Array{String}}=[])

  @info "Reading file: $filename"
  if !isfile(filename)
    error("File $filename does not exist!")
  end
  fid = open(filename, "r")

  # create header dictionary
  id1 = read!(fid, Array{UInt8}(undef, 1))
  id2 = read!(fid, Array{UInt8}(undef, BDF_ID_BYTES - 1))
  text1 = String(read!(fid, Array{UInt8}(undef, BDF_TEXT_BYTES)))
  text2 = String(read!(fid, Array{UInt8}(undef, BDF_TEXT_BYTES)))
  start_date = String(read!(fid, Array{UInt8}(undef, BDF_DATE_BYTES)))
  start_time = String(read!(fid, Array{UInt8}(undef, BDF_TIME_BYTES)))
  num_bytes_header = parse(Int, String(read!(fid, Array{UInt8}(undef, BDF_HEADER_SIZE_BYTES))))
  data_format = strip(String(read!(fid, Array{UInt8}(undef, BDF_DATA_FORMAT_BYTES))))
  num_data_records = parse(Int, String(read!(fid, Array{UInt8}(undef, BDF_RECORD_COUNT_BYTES))))
  duration_data_records = parse(Int, String(read!(fid, Array{UInt8}(undef, BDF_DURATION_BYTES))))
  num_channels = parse(Int, String(read!(fid, Array{UInt8}(undef, BDF_CHANNEL_COUNT_BYTES))))
  channel_labels = split(String(read!(fid, Array{UInt8}(undef, BDF_CHANNEL_LABEL_BYTES * num_channels))))
  transducer_type = split(String(read!(fid, Array{UInt8}(undef, BDF_TRANSDUCER_BYTES * num_channels))), r"[ ]{2,}", keepempty=false)
  channel_unit = split(String(read!(fid, Array{UInt8}(undef, BDF_UNIT_BYTES * num_channels))))
  physical_min = parse.(Int, [String(read!(fid, Array{UInt8}(undef, BDF_VALUE_BYTES))) for _ in 1:num_channels])
  physical_max = parse.(Int, [String(read!(fid, Array{UInt8}(undef, BDF_VALUE_BYTES))) for _ in 1:num_channels])
  digital_min = parse.(Int, [String(read!(fid, Array{UInt8}(undef, BDF_VALUE_BYTES))) for _ in 1:num_channels])
  digital_max = parse.(Int, [String(read!(fid, Array{UInt8}(undef, BDF_VALUE_BYTES))) for _ in 1:num_channels])
  pre_filter = split(String(read!(fid, Array{UInt8}(undef, BDF_FILTER_BYTES * num_channels))), r"[ ]{2,}", keepempty=false)
  num_samples = parse.(Int, [String(read!(fid, Array{UInt8}(undef, BDF_SAMPLES_BYTES))) for _ in 1:num_channels])
  reserved = split(String(read!(fid, Array{UInt8}(undef, BDF_RESERVED_BYTES * num_channels))))
  scale_factor = convert(Array{Float32}, ((physical_max .- physical_min) ./ (digital_max .- digital_min)))
  sample_rate = convert(Array{Int}, num_samples ./ duration_data_records)

  hd = BiosemiHeader(id1, id2, text1, text2, start_date, start_time, num_bytes_header, data_format,
    num_data_records, duration_data_records, num_channels, channel_labels, transducer_type,
    channel_unit, physical_min, physical_max, digital_min, digital_max, pre_filter,
    num_samples, reserved, scale_factor, sample_rate)

  if header_only
    close(fid)
    return hd
  end

  # read data
  bdf = read!(fid, Array{UInt8}(undef, BDF_SAMPLES_PER_BYTE * (num_data_records * num_channels * num_samples[1])))
  close(fid)

  channels = !isempty(channels) ? channel_index(channel_labels, channels) : 1:num_channels

  dat, time, trig, status = bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples[1], sample_rate[1])
  channels != 1:num_channels && update_header_bdf!(hd, channels)

  triggers = trigger_info(trig, sample_rate[1])

  return BiosemiData(filename, hd, dat, time, triggers, status)

end


"""
    bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples, sample_rate)

Convert raw BDF binary data to Julia data matrix and extract trigger information.

# Arguments
- `bdf::Vector{UInt8}`: Raw binary data from BDF file
- `num_channels::Int`: Total number of channels in file
- `channels::Vector{Int}`: Selected channel indices (including status channel)
- `scale_factor::Vector{Float32}`: Scale factors for each channel
- `num_data_records::Int`: Number of data records
- `num_samples::Int`: Number of samples per record per channel
- `sample_rate::Int`: Sampling rate in Hz

# Returns
- `dat_chans::Matrix{Float32}`: Scaled EEG data matrix (samples × channels)
- `time::StepRangeLen{Float64}`: Time vector for each sample
- `trig_chan::Vector{Int16}`: Raw trigger values for each sample
- `status_chan::Vector{Int16}`: Status channel values for each sample

# Notes
- Converts 24-bit BDF format to 32-bit float data
- Automatically applies scale factors for physical units
- Extracts trigger and status information from the last channel
- Time vector starts at 0 and increments by 1/sample_rate
- This is an internal function used by `read_bdf`
"""
function bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples, sample_rate)

  dat_chans = Matrix{Float32}(undef, (num_data_records * num_samples), length(channels) - 1)
  time = time_range(sample_rate, num_data_records)
  trig_chan = Vector{Int16}(undef, num_data_records * num_samples)
  status_chan = Vector{Int16}(undef, num_data_records * num_samples)

  # Use Set for faster channel lookup
  channels_set = Set(channels)

  pos = 1
  for rec = 0:(num_data_records-1)
    offset = rec * num_samples
    chan_idx = 1
    for chan = 1:num_channels
      if chan in channels_set  # selected channel
        if chan < num_channels
          for samp = 1:num_samples
            @inbounds dat_chans[offset+samp, chan_idx] = Float32(((Int32(bdf[pos]) << 8) | (Int32(bdf[pos+1]) << 16) | (Int32(bdf[pos+2]) << 24)) >> 8) * scale_factor[chan]
            pos += BDF_SAMPLES_PER_BYTE
          end
        else  # last channel is always Status channel
          for samp = 1:num_samples
            @inbounds trig_chan[offset+samp] = Int16(bdf[pos]) | (Int16(bdf[pos+1]) << 8)
            @inbounds status_chan[offset+samp] = Int16(bdf[pos+2])
            pos += BDF_SAMPLES_PER_BYTE
          end
        end
        chan_idx += 1
      else # channel not selected
        pos += num_samples * BDF_SAMPLES_PER_BYTE
      end
    end
  end

  return dat_chans, time, trig_chan, status_chan

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


"""
    write_bdf(bdf_in, filename="")

Write BioSemi BDF data structure to a BDF file.

# Arguments
- `bdf_in::BiosemiData`: Data structure to write
- `filename::String=""`: Output filename. If empty, uses `bdf_in.filename`

# File Format
Writes data in BioSemi BDF 24-bit format with:
- 256-byte header containing metadata
- 24-bit data samples for each channel
- Status channel as the last channel
- Proper scaling and calibration information

# Header Information
The header includes:
- File identification and metadata
- Channel information (labels, units, ranges)
- Sampling rate and data record information
- Pre-filtering and transducer information

# Examples
```julia
# Write to specified filename
write_bdf(dat, "output.bdf")

# Write using filename from data structure
dat.filename = "processed.bdf"
write_bdf(dat)
```

# Notes
- Automatically converts data back to 24-bit BDF format
- Applies inverse scaling to restore original digital values
- Preserves all header metadata and channel information
- File format follows [BioSemi BDF specification](https://www.biosemi.com/faq_file_format.htm)
- Overwrites existing files without warning

# See also
- `read_bdf`: Read BDF files
- `crop_bdf`: Reduce data length before writing
- `downsample_bdf`: Reduce sampling rate before writing
"""
function write_bdf(bdf_in::BiosemiData, filename::String="")

  if isempty(filename) 
    filename = bdf_in.filename
  end
  fid = open(filename, "w") 

  write(fid, 0xff)
  [write(fid, UInt8(i)) for i in bdf_in.header.id2]
  [write(fid, UInt8(i)) for i in bdf_in.header.text1]
  [write(fid, UInt8(i)) for i in bdf_in.header.text2]
  [write(fid, UInt8(i)) for i in bdf_in.header.start_date]
  [write(fid, UInt8(i)) for i in bdf_in.header.start_time]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_bytes_header), BDF_HEADER_SIZE_BYTES)]
  [write(fid, UInt8(i)) for i in rpad(bdf_in.header.data_format, BDF_DATA_FORMAT_BYTES)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_data_records), BDF_RECORD_COUNT_BYTES)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.duration_data_records), BDF_DURATION_BYTES)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_channels), BDF_CHANNEL_COUNT_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.channel_labels for j in rpad(i, BDF_CHANNEL_LABEL_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.transducer_type for j in rpad(i, BDF_TRANSDUCER_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.channel_unit for j in rpad(i, BDF_UNIT_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.physical_min for j in rpad(i, BDF_VALUE_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.physical_max for j in rpad(i, BDF_VALUE_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.digital_min for j in rpad(i, BDF_VALUE_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.digital_max for j in rpad(i, BDF_VALUE_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.pre_filter for j in rpad(i, BDF_FILTER_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.num_samples for j in rpad(i, BDF_SAMPLES_BYTES)]
  [write(fid, UInt8(j)) for i in bdf_in.header.reserved for j in rpad(i, BDF_RESERVED_BYTES)]

  data = round.(Int32, (bdf_in.data ./ transpose(bdf_in.header.scale_factor[1:end-1])))
  trigs = bdf_in.triggers.raw
  status = bdf_in.status
  num_data_records = bdf_in.header.num_data_records
  num_samples = bdf_in.header.num_samples[1]
  num_channels = bdf_in.header.num_channels

  bdf = matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)

  @info "Writing file: $filename"
  write(fid, Array{UInt8}(bdf))
  close(fid)

end


"""
    matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)

Convert Julia data matrix back to BDF 24-bit binary format.

# Arguments
- `data::Matrix{<:Number}`: EEG data matrix (samples × channels)
- `trigs::Vector{<:Integer}`: Trigger values for each sample
- `status::Vector{<:Integer}`: Status channel values for each sample
- `num_data_records::Int`: Number of data records
- `num_samples::Int`: Number of samples per record per channel
- `num_channels::Int`: Total number of channels (including status)

# Returns
- `Vector{UInt8}`: Binary data in BDF 24-bit format

# Format
- Each sample is stored as 3 bytes (24-bit)
- Data is written in record-major order (record × channel × sample)
- Status channel is written last for each record
- This is an internal function used by `write_bdf`

# Notes
- Converts float data back to integer format
- Applies inverse scaling to restore original digital values
- Maintains BDF file structure and byte ordering
"""
function matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)
  bdf = Array{UInt8}(undef, BDF_SAMPLES_PER_BYTE * (num_data_records * num_channels * num_samples))
  pos = 1
  for rec = 0:(num_data_records-1)
    for chan = 1:num_channels
      if chan < num_channels
        for samp = 1:num_samples
          data_val = data[rec*num_samples+samp, chan]
          bdf[pos] = (data_val % UInt8)
          bdf[pos+1] = ((data_val >> 8) % UInt8)
          bdf[pos+2] = ((data_val >> 16) % UInt8)
          pos += BDF_SAMPLES_PER_BYTE
        end
      else  # last channel is Status channel
        for samp = 1:num_samples
          trig_val = trigs[rec*num_samples+samp]
          status_val = status[rec*num_samples+samp]
          bdf[pos] = trig_val % UInt8
          bdf[pos+1] = (trig_val >> 8) % UInt8
          bdf[pos+2] = (status_val) % UInt8
          pos += BDF_SAMPLES_PER_BYTE
        end
      end
    end
  end
  return bdf
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
    delete_channels_bdf!(bdf, channels)

Remove specified channels from a BioSemi BDF data structure (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to remove
  - Integer: channel index (1-based)
  - String: channel label
  - Vector: multiple channels to remove

# Returns
- `Nothing`: Modifies `bdf` in-place

# Examples
```julia
# Remove single channel by index
delete_channels_bdf!(dat, 1)

# Remove single channel by label
delete_channels_bdf!(dat, "Fp1")

# Remove multiple channels
delete_channels_bdf!(dat, [1, 3, 5])
delete_channels_bdf!(dat, ["Fp1", "Fp2"])
```

# Notes
- Modifies the original data structure
- Updates header information (num_channels, labels, etc.)
- Recalculates header size and metadata
- Status channel is always preserved
- Use `delete_channels_bdf` for non-mutating version

# See also
- `delete_channels_bdf`: Non-mutating version
- `select_channels_bdf!`: Keep only specified channels
- `read_bdf`: Read data with channel selection
"""
function delete_channels_bdf!(bdf::BiosemiData, channels::Union{Array{Int},Array{String}})
  @info "Deleting channels: $channels"
  channels = channel_index(bdf.header.channel_labels, channels)
  channels = filter(x -> !(x in channels[1:end-1]), 1:length(bdf.header.channel_labels))
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end
delete_channels_bdf!(bdf_in::BiosemiData, channels::Union{Int,String}) = delete_channels_bdf!(bdf_in, [channels])

"""
    delete_channels_bdf(bdf_in, channels)

Remove specified channels from a BioSemi BDF data structure (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to remove

# Returns
- `BiosemiData`: New data structure with channels removed

# Examples
```julia
# Remove channels and get new structure
dat_clean = delete_channels_bdf(dat, [1, 2])
dat_clean = delete_channels_bdf(dat, "Fp1")
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `delete_channels_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `delete_channels_bdf!`: In-place version
- `select_channels_bdf`: Keep only specified channels
"""
function delete_channels_bdf(bdf_in::BiosemiData, channels::Union{Array{Int},Array{String}})
  bdf_out = deepcopy(bdf_in)
  delete_channels_bdf!(bdf_out, channels)
  return bdf_out
end
delete_channels_bdf(bdf_in::BiosemiData, channels::Union{Int,String}) = delete_channels_bdf(bdf_in, [channels])


"""
    select_channels_bdf!(bdf, channels)

Keep only specified channels in a BioSemi BDF data structure (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to keep
  - Integer: channel index (1-based)
  - String: channel label
  - Vector: multiple channels to keep

# Returns
- `Nothing`: Modifies `bdf` in-place

# Examples
```julia
# Keep single channel by index
select_channels_bdf!(dat, 1)

# Keep single channel by label
select_channels_bdf!(dat, "Cz")

# Keep multiple channels
select_channels_bdf!(dat, [1, 3, 5])
select_channels_bdf!(dat, ["Fp1", "Cz", "O1"])
```

# Notes
- Modifies the original data structure
- Updates header information (num_channels, labels, etc.)
- Recalculates header size and metadata
- Status channel is always preserved
- Use `select_channels_bdf` for non-mutating version

# See also
- `select_channels_bdf`: Non-mutating version
- `delete_channels_bdf!`: Remove specified channels
- `read_bdf`: Read data with channel selection
"""
function select_channels_bdf!(bdf::BiosemiData, channels::Union{Array{Int},Array{String}})
  @info "Selecting channels: $channels"
  channels = channel_index(bdf.header.channel_labels, channels)
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end
select_channels_bdf!(bdf_in::BiosemiData, channels::Union{Int,String}) = select_channels_bdf!(bdf_in, [channels])

"""
    select_channels_bdf(bdf_in, channels)

Keep only specified channels in a BioSemi BDF data structure (non-mutating).

# Arguments
- `bdf_in::BiosemiData`: Input data structure
- `channels::Union{Vector{<:Union{Int,String}}, Int, String}`: Channels to keep

# Returns
- `BiosemiData`: New data structure with only specified channels

# Examples
```julia
# Select channels and get new structure
dat_selected = select_channels_bdf(dat, [1, 3, 5])
dat_selected = select_channels_bdf(dat, ["Fp1", "Cz"])
```

# Notes
- Returns a new data structure (original unchanged)
- Calls `select_channels_bdf!` internally
- Useful when you want to preserve the original data

# See also
- `select_channels_bdf!`: In-place version
- `delete_channels_bdf`: Remove specified channels
"""
function select_channels_bdf(bdf_in::BiosemiData, channels::Union{Array{Int},Array{String}})
  bdf_out = deepcopy(bdf_in)
  select_channels_bdf!(bdf_out, channels)
  return bdf_out
end
select_channels_bdf(bdf_in::BiosemiData, channels::Union{Int,String}) = select_channels_bdf(bdf_in, [channels])


"""
    crop_bdf!(bdf, crop_type, val)

Reduce the length of BDF data by cropping (in-place).

# Arguments
- `bdf::BiosemiData`: Data structure to modify
- `crop_type::String`: Cropping method
  - `"records"`: Crop by data record numbers
  - `"triggers"`: Crop between trigger events
- `val::Vector{Int}`: Cropping parameters
  - For `"records"`: [start_record, end_record] (1-based)
  - For `"triggers"`: [start_trigger, end_trigger]

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


"""
downsample_bdf(bdf_in::BiosemiData, dec::Int)
  Reduce the sampling rate within a BiosemiData struct by an integer factor (dec).
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat2 = downsample_bdf(dat1, 2)
  ```
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
    update_header_bdf!(hd, channels)

Update header information after channel selection/deletion.

# Arguments
- `hd::BiosemiHeader`: Header to modify
- `channels::Vector{Int}`: Selected channel indices

# Returns
- `Nothing`: Modifies `hd` in-place

# Updated Fields
- `num_channels`: Number of selected channels
- `channel_labels`: Labels for selected channels
- `transducer_type`: Transducer types for selected channels
- `channel_unit`: Units for selected channels
- `physical_min/max`: Physical ranges for selected channels
- `digital_min/max`: Digital ranges for selected channels
- `pre_filter`: Pre-filtering info for selected channels
- `num_samples`: Sample counts for selected channels
- `reserved`: Reserved space for selected channels
- `sample_rate`: Sample rates for selected channels
- `scale_factor`: Scale factors for selected channels
- `num_bytes_header`: Recalculated header size

# Notes
- This is an internal function used by channel selection functions
- Automatically recalculates header size based on channel count
- Preserves all metadata for selected channels
- Status channel is always included automatically
"""
function update_header_bdf!(hd::BiosemiHeader, channels::Array{Int})
  hd.num_channels = length(channels)
  hd.physical_min = hd.physical_min[channels]
  hd.physical_max = hd.physical_max[channels]
  hd.digital_min = hd.digital_min[channels]
  hd.digital_max = hd.digital_max[channels]
  hd.scale_factor = hd.scale_factor[channels]
  hd.transducer_type = hd.transducer_type[channels]
  hd.num_samples = hd.num_samples[channels]
  hd.channel_unit = hd.channel_unit[channels]
  hd.reserved = hd.reserved[channels]
  hd.sample_rate = hd.sample_rate[channels]
  hd.channel_labels = hd.channel_labels[channels]
  hd.pre_filter = hd.pre_filter[channels]
  hd.num_bytes_header = (length(channels) + BDF_STATUS_CHANNEL_OFFSET) * BDF_HEADER_SIZE
end


"""
channel_index(labels::Array{AbstractString}, channels::Array{AbstractString})
  Return channel index given labels and desired selection.
"""
function channel_index(labels, channels::Array{String})
  channels = [findfirst(x .== labels) for x in channels]
  any(channels .=== nothing) && error("Requested channel label is not in the bdf file!")
  return sort(unique(append!(channels, length(labels))))
end
channel_index(labels, channels::String) = channel_index(labels, [channels])


"""
channel_index(labels::Array{AbstractString}, channels::Array{Int})
  Return channel index given labels and desired selection.
"""
function channel_index(labels, channels::Array{Int})
  trigSelected = findall(x -> x == -1, channels)
  if length(trigSelected) > 0
    channels[trigSelected] = repeat([length(labels)], length(trigSelected))
  end
  any(channels .> length(labels)) && error("Requested channel number greater than number of channels in file!")
  any(channels .< 1) && error("Requested channel number less than 1!")
  return sort(unique(append!(channels, length(labels))))
end
channel_index(labels, channels::Int) = channel_index(labels, [channels])

end # module
