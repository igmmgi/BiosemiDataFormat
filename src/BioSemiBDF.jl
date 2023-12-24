module BioSemiBDF

using DSP, OrderedCollections

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

mutable struct BioSemiHeader
  id1
  id2
  text1
  text2
  start_date
  start_time
  num_bytes_header
  data_format
  num_data_records
  duration_data_records
  num_channels
  channel_labels
  transducer_type
  channel_unit
  physical_min
  physical_max
  digital_min
  digital_max
  pre_filter
  num_samples
  reserved
  scale_factor
  sample_rate
end

mutable struct BioSemiTriggers
  raw
  idx
  val
  count
  time
end

mutable struct BioSemiData
  filename
  header
  data
  time
  triggers
  status
end

function Base.show(io::IO, hdr::BioSemiHeader)
  println(io, "Number of Channels: ", length(hdr.channel_labels) - 1)
  println(io, "Channel Labels: ", join(hdr.channel_labels[1:end-1], ", "))
  println(io, "Number of Data Records: ", hdr.num_data_records)
  println(io, "Sample Rate: ", hdr.sample_rate[1])
end

function Base.show(io::IO, trig::BioSemiTriggers)
  println(io, "Triggers (Value => Count): ", join(trig.count, ", "))
end

function Base.show(io::IO, dat::BioSemiData)
  println(io, "Filename: $(dat.filename)")
  print(io, dat.header)
  println(io, "Data Size: $(size(dat.data))")
  println(io, dat.triggers)
end


"""
read_bdf(filename::String; header_only::Bool=false, channels::Union{Array{Any}, Array{Int}, Array{String}}}=[])

  Reads BioSemi Data Format (bdf) files.
  See https://www.biosemi.com/faq_file_format.htm for file format details.
  ### Inputs:
  * filename String
  * header_only Bool=false
  * channels Array Int/String
  ### Outputs:
  * BioSemiData struct with the following fields
  * header BioSemiHeader
  * data (samples * channels) Matrix
  * time
  * triggers BioSemiTriggers
  * status Array
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat1 = read_bdf("filename1.bdf", header_only = true)
  dat1 = read_bdf("filename1.bdf", channels = [1,3,5])
  dat1 = read_bdf("filename1.bdf", channels = ["Fp1", "Cz"])
  dat1 = read_bdf("filename1.bdf", channels = [-1])  # trigger channel only
  ```
"""
function read_bdf(filename::String; header_only::Bool=false, channels::Union{Array{Any},Array{Int},Array{String}}=[])

  fid = open(filename, "r")

  # create header dictionary
  id1 = read!(fid, Array{UInt8}(undef, 1))
  id2 = read!(fid, Array{UInt8}(undef, 7))
  text1 = String(read!(fid, Array{UInt8}(undef, 80)))
  text2 = String(read!(fid, Array{UInt8}(undef, 80)))
  start_date = String(read!(fid, Array{UInt8}(undef, 8)))
  start_time = String(read!(fid, Array{UInt8}(undef, 8)))
  num_bytes_header = parse(Int, String(read!(fid, Array{UInt8}(undef, 8))))
  data_format = strip(String(read!(fid, Array{UInt8}(undef, 44))))
  num_data_records = parse(Int, String(read!(fid, Array{UInt8}(undef, 8))))
  duration_data_records = parse(Int, String(read!(fid, Array{UInt8}(undef, 8))))
  num_channels = parse(Int, String(read!(fid, Array{UInt8}(undef, 4))))
  channel_labels = split(String(read!(fid, Array{UInt8}(undef, 16 * num_channels))))
  transducer_type = split(String(read!(fid, Array{UInt8}(undef, 80 * num_channels))), r"[ ]{2,}", keepempty=false)
  channel_unit = split(String(read!(fid, Array{UInt8}(undef, 8 * num_channels))))
  physical_min = parse.(Int, [String(read!(fid, Array{UInt8}(undef, 8))) for _ in 1:num_channels])
  physical_max = parse.(Int, [String(read!(fid, Array{UInt8}(undef, 8))) for _ in 1:num_channels])
  digital_min = parse.(Int, [String(read!(fid, Array{UInt8}(undef, 8))) for _ in 1:num_channels])
  digital_max = parse.(Int, [String(read!(fid, Array{UInt8}(undef, 8))) for _ in 1:num_channels])
  pre_filter = split(String(read!(fid, Array{UInt8}(undef, 80 * num_channels))), r"[ ]{2,}", keepempty=false)
  num_samples = parse.(Int, [String(read!(fid, Array{UInt8}(undef, 8))) for _ in 1:num_channels])
  reserved = split(String(read!(fid, Array{UInt8}(undef, 32 * num_channels))))
  scale_factor = convert(Array{Float32}, ((physical_max .- physical_min) ./ (digital_max .- digital_min)))
  sample_rate = convert(Array{Int}, num_samples ./ duration_data_records)

  hd = BioSemiHeader(id1, id2, text1, text2, start_date, start_time, num_bytes_header, data_format,
    num_data_records, duration_data_records, num_channels, channel_labels, transducer_type,
    channel_unit, physical_min, physical_max, digital_min, digital_max, pre_filter,
    num_samples, reserved, scale_factor, sample_rate)

  if header_only
    close(fid)
    return hd
  end

  # read data
  bdf = read!(fid, Array{UInt8}(undef, 3 * (num_data_records * num_channels * num_samples[1])))
  close(fid)

  channels = !isempty(channels) ? channel_idx(channel_labels, channels) : 1:num_channels

  dat, time, trig, status = bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples[1], sample_rate[1])
  channels != 1:num_channels && update_header_bdf!(hd, channels)

  triggers = triggerInfo(trig, sample_rate[1])

  return BioSemiData(filename, hd, dat, time, triggers, status)

end



"""
bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples)

  Internal functon used within read_bdf to read BioSemi 24bit data representation
  into julia data array/matrix
"""
function bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples, sample_rate)

  dat_chans = Matrix{Float32}(undef, (num_data_records * num_samples), length(channels) - 1)
  time = time_range(sample_rate, num_data_records)
  trig_chan = Array{Int16}(undef, num_data_records * num_samples)
  status_chan = Array{Int16}(undef, num_data_records * num_samples)

  pos = 1
  for rec = 0:(num_data_records-1)
    offset = rec * num_samples
    chan_idx = 1
    for chan = 1:num_channels
      if chan in channels  # selected channel
        if chan < num_channels
          for samp = 1:num_samples
            @inbounds dat_chans[offset+samp, chan_idx] = Float32(((Int32(bdf[pos]) << 8) | (Int32(bdf[pos+1]) << 16) | (Int32(bdf[pos+2]) << 24)) >> 8) * scale_factor[chan]
            pos += 3
          end
        else  # last channel is always Status channel
          for samp = 1:num_samples
            @inbounds trig_chan[offset+samp] = ((Int16(bdf[pos])) | (Int16(bdf[pos+1]) << 8))
            @inbounds status_chan[offset+samp] = Int16(bdf[pos+2])
            pos += 3
          end
        end
        chan_idx += 1
      else # channel not selected
        pos += num_samples * 3
      end
    end
  end

  return dat_chans, time, trig_chan, status_chan

end




"""
triggerInfo(trig)
  Internal functon used within read_bdf to convert raw trigger line values to a summary of trigger information 
"""
function triggerInfo(trig_raw, sample_rate)

  # trigger events + time
  trig_idx = findall(diff(trig_raw) .>= 1) .+ 1
  trig_val = trig_raw[trig_idx]
  trig_time = hcat(trig_val, pushfirst!(diff(trig_idx), 0) / sample_rate)

  # unique triggers + count
  trig_count = OrderedDict(i => 0 for i in sort!(collect(Set(trig_val))))
  for val in trig_val
    trig_count[val] += 1
  end

  triggers = BioSemiTriggers(trig_raw, trig_idx, trig_val, trig_count, trig_time)

  return triggers

end


"""
write_bdf(bdf_in::BioSemiData)

  Write BioSemiRaw structs to *.bdf file.
  See https://www.biosemi.com/faq_file_format.htm for file format details.
  ### Inputs:
  * BioSemiData struct
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  write_bdf(dat1)
"""
function write_bdf(bdf_in::BioSemiData, filename::String="")

  isempty(filename) ? fid = open(bdf_in.filename, "w") : fid = open(filename, "w")

  write(fid, 0xff)
  [write(fid, UInt8(i)) for i in bdf_in.header.id2]
  [write(fid, UInt8(i)) for i in bdf_in.header.text1]
  [write(fid, UInt8(i)) for i in bdf_in.header.text2]
  [write(fid, UInt8(i)) for i in bdf_in.header.start_date]
  [write(fid, UInt8(i)) for i in bdf_in.header.start_time]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_bytes_header), 8)]
  [write(fid, UInt8(i)) for i in rpad(bdf_in.header.data_format, 44)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_data_records), 8)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.duration_data_records), 8)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header.num_channels), 4)]
  [write(fid, UInt8(j)) for i in bdf_in.header.channel_labels for j in rpad(i, 16)]
  [write(fid, UInt8(j)) for i in bdf_in.header.transducer_type for j in rpad(i, 80)]
  [write(fid, UInt8(j)) for i in bdf_in.header.channel_unit for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.physical_min for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.physical_max for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.digital_min for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.digital_max for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.pre_filter for j in rpad(i, 80)]
  [write(fid, UInt8(j)) for i in bdf_in.header.num_samples for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header.reserved for j in rpad(i, 32)]

  data = round.(Int32, (bdf_in.data ./ transpose(bdf_in.header.scale_factor[1:end-1])))
  trigs = bdf_in.triggers.raw
  status = bdf_in.status
  num_data_records = bdf_in.header.num_data_records
  num_samples = bdf_in.header.num_samples[1]
  num_channels = bdf_in.header.num_channels

  bdf = matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)

  write(fid, Array{UInt8}(bdf))
  close(fid)

end


"""
matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)

  Internal functon used within write_bdf to write Julia BioSemiBDF data matrix
  to bdf 24bit file format.
"""
function matrix2bdf(data, trigs, status, num_data_records, num_samples, num_channels)

  bdf = Array{UInt8}(undef, 3 * (num_data_records * num_channels * num_samples))
  pos = 1
  for rec = 0:(num_data_records-1)
    for chan = 1:num_channels
      if chan < num_channels
        for samp = 1:num_samples
          data_val = data[rec*num_samples+samp, chan]
          bdf[pos] = (data_val % UInt8)
          bdf[pos+1] = ((data_val >> 8) % UInt8)
          bdf[pos+2] = ((data_val >> 16) % UInt8)
          pos += 3
        end
      else  # last channel is Status channel
        for samp = 1:num_samples
          trig_val = trigs[rec*num_samples+samp]
          status_val = status[rec*num_samples+samp]
          bdf[pos] = trig_val % UInt8
          bdf[pos+1] = (trig_val >> 8) % UInt8
          bdf[pos+2] = (status_val) % UInt8
          pos += 3
        end
      end
    end
  end
  return bdf
end


"""
merge_bdf(bdfs::Array{BioSemiData}, filename::String="merged.bdf")

  Merge BioSemiRaw structs to single BioSemiRaw struct. Checks that the
  input BioSemiRaw structs have the same number of channels, same channel
  labels and that each channel has the same sample rate.
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat2 = read_bdf("filename2.bdf")
  dat3 = merge_bdf([dat1, dat2])
  ```
"""
function merge_bdf(bdfs::Array{BioSemiData})

  # check data structs to merge have same number of channels, channel labels + sample rate
  num_chans = (x -> x.header.num_channels).(bdfs)
  !all(x -> x == num_chans[1], num_chans) && error("Different number of channels in bdf_in")
  chan_labels = (x -> x.header.channel_labels).(bdfs)
  !all(y -> y == chan_labels[1], chan_labels) && error("Different channel labels bdf_in")
  sample_rate = (x -> x.header.sample_rate).(bdfs)
  !all(y -> y == sample_rate[1], sample_rate) && error("Different sample rate in bdf_in")

  # make copy so that bdf_in is not altered
  bdf_out = deepcopy(bdfs[1])

  # merge data
  bdf_out.header.num_data_records = sum((x -> x.header.num_data_records).(bdfs))
  bdf_out.data = vcat((x -> x.data).(bdfs)...)
  bdf_out.status = vcat((x -> x.status).(bdfs)...)

  # recaculate trigger information
  trig = vcat((x -> x.triggers.raw).(bdfs)...)
  bdf_out.triggers = triggerInfo(trig, bdf_out.header.sample_rate[1])

  # merged time 
  println(sample_rate)
  bdf_out.time = time_range(bdf_out.header.sample_rate[1], bdf_out.header.num_data_records)

  return bdf_out

end

function time_range(sample_rate::Int, num_data_records::Int)
  return 0:1/sample_rate:(num_data_records-(1/sample_rate))
end

"""
delete_channels_bdf(bdf_in::BioSemiData, channels::Union{Array{Int}, Array{String}})

  Delete specific channels from BioSemiData struct. Channels can be specified
  using channel numbers or channel labels.
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat1 = delete_channels_bdf(dat, [1, 2])
  dat1 = delete_channels_bdf(dat, ["Fp1", "F1"])
  ```
"""

function delete_channels_bdf(bdf_in::BioSemiData, channels::Union{Array{Int},Array{String}})
  bdf_out = deepcopy(bdf_in)
  channels = channel_idx(bdf_out.header.channel_labels, channels)
  channels = filter(x -> !(x in channels[1:end-1]), 1:(length(bdf_in.header.channel_labels)))
  update_header_bdf!(bdf_out.header, channels)
  bdf_out.data = bdf_out.data[:, channels[1:end-1]]
  return bdf_out
end

function delete_channels_bdf!(bdf::BioSemiData, channels::Union{Array{Int},Array{String}})
  channels = channel_idx(bdf.header.channel_labels, channels)
  channels = filter(x -> !(x in channels[1:end-1]), 1:length(bdf.header.channel_labels))
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end


"""
select_channels_bdf(bdf_in::BioSemiData, channels::Union{Array{Int}, Array{String}})

  Select specific channels from BioSemiData struct. Channels can be specified
  using channel numbers or channel labels.
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat1 = select_channels_bdf(dat, [1, 2])
  dat1 = select_channels_bdf(dat, ["Fp1", "F1"])
  ```
"""
function select_channels_bdf(bdf_in::BioSemiData, channels::Union{Array{Int},Array{String}})
  bdf_out = deepcopy(bdf_in)
  channels = channel_idx(bdf_out.header.channel_labels, channels)
  update_header_bdf!(bdf_out.header, channels)
  bdf_out.data = bdf_out.data[:, channels[1:end-1]]
  return bdf_out
end

function select_channels_bdf!(bdf::BioSemiData, channels::Union{Array{Int},Array{String}})
  channels = channel_idx(bdf.header.channel_labels, channels)
  update_header_bdf!(bdf.header, channels)
  bdf.data = bdf.data[:, channels[1:end-1]]
end


function crop_check(crop_type::String, val::Array{Int}, triggers, num_data_records)
  crop_type ∉ ["triggers", "records"] && error("crop_type not recognized!")
  length(val) != 2 && error("val should be of length 2")
  if crop_type == "triggers"
    val[1] ∉ triggers && error("val[1] not available trigger")
    val[2] ∉ triggers && error("val[2] not available trigger")
  elseif crop_type == "records"
    val[1] < 1 && error("val[1] less than 1")
    val[2] > num_data_records && error("val[2] > number of data records")
  end
end


"""
crop_bdf(bdf_in::BioSemiData, crop_type::tString, val::Array{Int}, filename::String)

  Recuce the length of the recorded data. The border upon which to crop the bdf file can be defined using either
  a start and end trigger ("triggers") or a start and end record ("records").
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat2 = crop_bdf(dat1, "triggers", [1 2])   # between first trigger 1 and last trigger 2
  dat3 = crop_bdf(dat1, "records",  [1 100]) # data records 1 to 100 inclusive
  ```
"""
function crop_bdf(bdf_in::BioSemiData, crop_type::String, val::Array{Int})

  crop_check(crop_type, val, bdf_in.triggers.count.keys, bdf_in.header.num_data_records)

  sample_rate = bdf_in.header.sample_rate[1]
  nsamples = size(bdf_in.data, 1)

  idxStart, idxEnd = find_crop_index(bdf_in.triggers, crop_type, val, sample_rate, nsamples)

  # copy data and crop
  bdf_out = deepcopy(bdf_in)
  bdf_out.header.num_data_records = Int(((idxEnd - idxStart) + 1) / sample_rate)
  bdf_out.data = bdf_out.data[idxStart:idxEnd, :]
  bdf_out.time = time_range(bdf_out.header.sample_rate[1], bdf_out.header.num_data_records)
  bdf_out.status = bdf_out.status[idxStart:idxEnd]

  # recaculate trigger information
  trig = bdf_out.triggers.raw[idxStart:idxEnd]
  bdf_out.triggers = triggerInfo(trig, bdf_out.header.sample_rate[1])

  return bdf_out

end

function crop_bdf!(bdf::BioSemiData, crop_type::String, val::Array{Int})

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
  bdf.triggers = triggerInfo(trig, bdf.header.sample_rate[1])

end


function find_crop_index(triggers::BioSemiTriggers, crop_type::String, val::Array{Int}, sample_rate::Int, data_length::Int)

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
downsample_bdf(bdf_in::BioSemiData, dec::Int)

  Reduce the sampling rate within a BioSemiData struct by an integer factor (dec).
  ### Examples:
  ```julia
  dat1 = read_bdf("filename1.bdf")
  dat2 = downsample_bdf(dat1, 2)
  ```
"""
function downsample_bdf(bdf_in::BioSemiData, dec::Int)

  !ispow2(dec) && error("dec should be power of 2!")

  # padding at start/end
  nsamp = dec * 20  # enough samples?
  ndec = div(nsamp, dec)

  bdf_out = deepcopy(bdf_in)
  data = Matrix{Float32}(undef, div(size(bdf_out.data, 1), dec), size(bdf_out.data, 2))
  for i in 1:size(bdf_out.data, 2)
    tmp_dat = resample([reverse(bdf_out.data[1:nsamp, i]); bdf_out.data[:, i]; reverse(bdf_out.data[end-(nsamp-1):end, i])], 1 / dec)
    data[:, i] = convert(Array{Float32}, tmp_dat[ndec+1:end-ndec])
  end
  bdf_out.data = data

  bdf_out.header.sample_rate = div.(bdf_out.header.sample_rate, dec)
  bdf_out.header.num_samples = div.(bdf_out.header.num_samples, dec)
  bdf_out.time = (0:size(bdf_out.data, 1)-1) / bdf_out.header.sample_rate[1]

  # update triggers
  bdf_out.triggers.raw = zeros(Int16, 1, size(bdf_out.data, 1))
  bdf_out.triggers.idx = convert(Array{Int64}, round.(bdf_out.triggers.idx / dec))
  bdf_out.triggers.raw[bdf_out.triggers.idx] = bdf_out.triggers.val

  return bdf_out

end

function downsample_bdf!(bdf::BioSemiData, dec::Int)

  !ispow2(dec) && error("dec should be power of 2!")

  # padding at start/end
  nsamp = dec * 20  # enough samples?
  ndec = div(nsamp, dec)

  data = Matrix{Float32}(undef, div(size(bdf.data, 1), dec), size(bdf.data, 2))
  for i in 1:size(bdf.data, 1)
    tmp_dat = resample([reverse(bdf.data[1:nsamp, :]); bdf.data[:, i]; reverse(bdf.data[end-(nsamp-1):end, i])], 1 / dec)
    data[:, i] = convert(Array{Float32}, tmp_dat[ndec+1:end-ndec])
  end
  bdf.data = data

  bdf.header.sample_rate = div.(bdf.header.sample_rate, dec)
  bdf.header.num_samples = div.(bdf.header.num_samples, dec)
  bdf.time = (0:size(bdf.data, 1)-1) / bdf.header.sample_rate[1]

  # update triggers
  bdf.triggers.raw = zeros(Int16, 1, size(bdf.data, 1))
  bdf.triggers.idx = convert(Array{Int64}, round.(bdf.triggers.idx / dec))
  bdf.triggers.raw[bdf.triggers.idx] = bdf.triggers.val

end


"""
update_header_bdf(hd::Dict, channels::Array{Int})

  Updates header Dict within BioSemiData struct following the selection
  of specific channels in read_bdf or select_channels_bdf.
"""
function update_header_bdf!(hd::BioSemiHeader, channels::Array{Int})
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
  hd.num_bytes_header = (length(channels) + 1) * 256
end


"""
channel_idx(labels::Array{AbstractString}, channels::Array{AbstractString})

  Return channel index given labels and desired selection.
"""
function channel_idx(labels, channels::Array{String})
  channels = [findfirst(x .== labels) for x in channels]
  any(channels .=== nothing) && error("A requested channel label is not in the bdf file!")
  return sort(unique(append!(channels, length(labels))))
end


"""
channel_idx(labels::Array{AbstractString}, channels::Array{Int})

  Return channel index given labels and desired selection.
"""
function channel_idx(labels, channels::Array{Int})
  trigSelected = findall(x -> x == -1, channels)
  if length(trigSelected) > 0
    channels[trigSelected] = repeat([length(labels)], length(trigSelected))
  end
  any(channels .> length(labels)) && error("A requested channel number greater than number of channels in file!")
  any(channels .< 1) && error("A requested channel number less than 1!")
  return sort(unique(append!(channels, length(labels))))
end

channel_idx(labels, channels::Int) = channel_idx(labels, [channels])

end


