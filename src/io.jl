# io.jl - Input/Output functions for BiosemiDataFormat

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
function read_bdf(filename::String; header_only::Bool=false, channels=[])

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



# Arguments
- `bdf_in::BiosemiData`: Data structure to write
- `filename::String=""`: Output filename (uses original filename if empty)

# Returns
- `Nothing`: Writes file to disk

# Examples
```julia
# Write to a new file
write_bdf(dat, "processed_data.bdf")

# Write using original filename
dat.filename = "my_data.bdf"
write_bdf(dat)

# Write processed data
dat_cropped = crop_bdf(dat, "triggers", [100, 200])
write_bdf(dat_cropped, "cropped_data.bdf")
```

# Notes
- Creates a new BDF file in standard BioSemi format
- Automatically applies scale factors and converts to 24-bit format
- Preserves all header information and metadata
- If no filename is provided, uses the filename stored in the data structure
- Overwrites existing files without warning

# See also
- `read_bdf`: Read BDF files
- `crop_bdf`: Reduce data before writing
- `select_channels_bdf`: Select channels before writing
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
    update_header_bdf!(hd, channels)

Update header information after channel selection/deletion.

# Arguments
- `hd::BiosemiHeader`: Header to modify
- `channels::Vector{Int}`: Selected channel indices (including status channel)

# Returns
- `Nothing`: Modifies `hd` in-place

# Updated Fields
- `num_channels`: Number of selected channels
- `channel_labels`: Labels of selected channels
- `transducer_type`: Transducer types for selected channels
- `channel_unit`: Units for selected channels
- `physical_min`: Physical minimum values for selected channels
- `physical_max`: Physical maximum values for selected channels
- `digital_min`: Digital minimum values for selected channels
- `digital_max`: Digital maximum values for selected channels
- `pre_filter`: Pre-filtering information for selected channels
- `num_samples`: Number of samples for selected channels
- `reserved`: Reserved header space for selected channels
- `scale_factor`: Scale factors for selected channels
- `sample_rate`: Sampling rates for selected channels
- `num_bytes_header`: Updated header size

# Notes
- This is an internal function used by channel manipulation functions
- Automatically recalculates `num_bytes_header` based on channel count
- Preserves the status channel as the last channel
- Updates all channel-specific header fields
"""
function update_header_bdf!(hd::BiosemiHeader, channels::Array{Int})

  hd.num_channels = length(channels)
  hd.channel_labels = hd.channel_labels[channels]
  hd.transducer_type = hd.transducer_type[channels]
  hd.channel_unit = hd.channel_unit[channels]
  hd.physical_min = hd.physical_min[channels]
  hd.physical_max = hd.physical_max[channels]
  hd.digital_min = hd.digital_min[channels]
  hd.digital_max = hd.digital_max[channels]
  hd.pre_filter = hd.pre_filter[channels]
  hd.num_samples = hd.num_samples[channels]
  hd.reserved = hd.reserved[channels]
  hd.scale_factor = hd.scale_factor[channels]
  hd.sample_rate = hd.sample_rate[channels]
  hd.num_bytes_header = (length(channels) + BDF_STATUS_CHANNEL_OFFSET) * BDF_HEADER_SIZE

end
