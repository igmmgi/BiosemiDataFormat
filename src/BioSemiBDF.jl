module BioSemiBDF

using
StatsBase,
DSP

export
crop_bdf,
downsample_bdf,
merge_bdf,
select_channels_bdf,
read_bdf,
write_bdf

mutable struct BioSemiData
  header::Dict
  data::Matrix
  time::Array
  triggers::Dict
  status::Array
end


"""
read_bdf(filename::String; header_only::Bool=false, channels::Union{Array{Any}, Array{Int}, Array{String}}}=[])
Reads BioSemi Data Format (bdf) files.
See https://www.biosemi.com/faq_file_format.htm for file format details.
### Inputs:
* filename String
* header_only Bool=false
* header_only Bool=false
* channels Array Int/String
### Outputs:
* BioSemiData struct with the following fields
* header Dict
* data (channels * samples) Matrix
* time
* triggers Dict
* status Array
### Examples:
```julia
dat1 = read_bdf("filename1.bdf")
dat1 = read_bdf("filename1.bdf", header_only = true)
dat1 = read_bdf("filename1.bdf", channels = [1,3,5])
dat1 = read_bdf("filename1.bdf", channels = ["Fp1", "Cz"])
```
"""
function read_bdf(filename::String; header_only::Bool=false, channels::Union{Array{Any}, Array{Int}, Array{String}}=[])

  fid = open(filename, "r")

  # create header dictionary
  hd1 = Dict{String, Any}(
                          "filename"              => filename,
                          "id1"                   => read!(fid, Array{UInt8}(undef, 1)),
                          "id2"                   => ascii(String(read!(fid, Array{UInt8}(undef, 7)))),
                          "text1"                 => ascii(String(read!(fid, Array{UInt8}(undef, 80)))),
                          "text2"                 => ascii(String(read!(fid, Array{UInt8}(undef, 80)))),
                          "start_date"            => ascii(String(read!(fid, Array{UInt8}(undef, 8)))),
                          "start_time"            => ascii(String(read!(fid, Array{UInt8}(undef, 8)))),
                          "num_bytes_header"      => parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef, 8)))))),
                          "data_format"           => strip(ascii(String(read!(fid, Array{UInt8}(undef, 44))))),
                          "num_data_records"      => parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef, 8)))))),
                          "duration_data_records" => parse(Int, ascii(String(read!(fid, Array{UInt8}(undef, 8))))),
                          "num_channels"          => parse(Int, ascii(String(read!(fid, Array{UInt8}(undef, 4)))))
                         )
  hd2 = Dict{String, Any}(
                          "channel_labels"  => [String(strip(ascii(String(read!(fid, Array{UInt8}(undef,  16)))))) for _ in 1:hd1["num_channels"]],
                          "transducer_type" => [String(strip(ascii(String(read!(fid, Array{UInt8}(undef, 80)))))) for _ in 1:hd1["num_channels"]],
                          "channel_unit"    => [String(strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "physical_min"    => [parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "physical_max"    => [parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "digital_min"     => [parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "digital_max"     => [parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "pre_filter"      => [String(strip(ascii(String(read!(fid, Array{UInt8}(undef,  80)))))) for _ in 1:hd1["num_channels"]],
                          "num_samples"     => [parse(Int, strip(ascii(String(read!(fid, Array{UInt8}(undef,  8)))))) for _ in 1:hd1["num_channels"]],
                          "reserved"        => [String(strip(ascii(String(read!(fid, Array{UInt8}(undef,  32)))))) for _ in 1:hd1["num_channels"]]
                         )
  hd3 = Dict{String, Any}(
                          "scale_factor" => convert(Array{Float32}, ((hd2["physical_max"] .- hd2["physical_min"]) ./ (hd2["digital_max"] .- hd2["digital_min"]))),
                          "sample_rate"  => convert(Array{Int}, hd2["num_samples"] ./ hd1["duration_data_records"])
                         )
  hd = merge(hd1, hd2, hd3)
  header_only && return hd

  # read data
  bdf = read!(fid, Array{UInt8}(undef, 3*(hd["num_data_records"]*hd["num_channels"]*hd["num_samples"][1])))
  close(fid)

  channels = !isempty(channels) ? channel_idx(hd["channel_labels"], channels) : collect(1:hd["num_channels"])

  dat, time, trig, status = bdf2matrix(bdf, hd["num_channels"], channels, hd["scale_factor"], hd["num_data_records"], hd["num_samples"], hd["sample_rate"])
  channels != collect(1:hd["num_channels"]) && update_header_bdf!(hd, channels)

  # events
  trig_idx = findall(diff(trig) .>= 1) .+ 1
  trig_val = trig[trig_idx]

  # create triggers dictionary
  triggers = Dict{String, Any}(
                               "raw"   => trig,
                               "idx"   => trig_idx,
                               "val"   => trig_val,
                               "count" => sort(countmap(trig_val)),
  "time"  => hcat(trig_val, pushfirst!(diff(trig_idx), 0) / hd["sample_rate"][1])
 )

  return BioSemiData(hd, dat, time, triggers, status)

end


"""
bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples)
Internal functon used within read_bdf to read BioSemi 24bit data representation
into julia data array/matrix
"""
function bdf2matrix(bdf, num_channels, channels, scale_factor, num_data_records, num_samples, sample_rate)

  dat_chans   = Matrix{Float32}(undef, length(channels)-1, (num_data_records*num_samples[1]))
  trig_chan   = Array{Int16}(undef, num_data_records*num_samples[1])
  status_chan = Array{Int16}(undef, num_data_records*num_samples[1])
  pos = 1;
  for rec = 0:(num_data_records-1)
    idx = 1
    for chan = 1:num_channels
      if chan in channels  # selected channel
        if chan < num_channels
          for samp = 1:num_samples[1]
            dat_chans[idx, rec*num_samples[1]+samp] = Float32(((Int32(bdf[pos]) << 8) | (Int32(bdf[pos+1]) << 16) | (Int32(bdf[pos+2]) << 24)) >> 8) * scale_factor[chan]
            pos += 3
          end
        else  # last channel is always Status channel
          for samp = 1:num_samples[1]
            trig_chan[rec*num_samples[1]+samp] = ((Int16(bdf[pos])) | (Int16(bdf[pos+1]) << 8))
            status_chan[rec*num_samples[1]+samp] = Int16(bdf[pos+2])
            pos += 3
          end
        end
        idx += 1
      else # channel not selected
        pos += num_samples[1]*3
      end
    end
  end
  time = collect(0:size(dat_chans, 2) - 1) / sample_rate[1]

  return dat_chans, time, trig_chan, status_chan

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

  isempty(filename) ? fid = open(bdf_in.header["filename"], "w") : fid = open(filename, "w")

  write(fid, 0xff)
  [write(fid, UInt8(i)) for i in bdf_in.header["id2"]]
  [write(fid, UInt8(i)) for i in bdf_in.header["text1"]]
  [write(fid, UInt8(i)) for i in bdf_in.header["text2"]]
  [write(fid, UInt8(i)) for i in bdf_in.header["start_date"]]
  [write(fid, UInt8(i)) for i in bdf_in.header["start_time"]]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header["num_bytes_header"]), 8)]
  [write(fid, UInt8(i)) for i in rpad(bdf_in.header["num_bytes_header"], 44)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header["num_data_records"]), 8)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header["duration_data_records"]), 8)]
  [write(fid, UInt8(i)) for i in rpad(string(bdf_in.header["num_channels"]), 4)]
  [write(fid, UInt8(j)) for i in bdf_in.header["channel_labels"] for j in rpad(i, 16)]
  [write(fid, UInt8(j)) for i in bdf_in.header["transducer_type"] for j in rpad(i, 80)]
  [write(fid, UInt8(j)) for i in bdf_in.header["channel_unit"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["physical_min"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["physical_max"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["digital_min"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["digital_max"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["pre_filter"] for j in rpad(i, 80)]
  [write(fid, UInt8(j)) for i in bdf_in.header["num_samples"] for j in rpad(i, 8)]
  [write(fid, UInt8(j)) for i in bdf_in.header["reserved"] for j in rpad(i, 32)]

  data             = round.(Int32, (bdf_in.data ./ bdf_in.header["scale_factor"][1:end-1]))
  trigs            = bdf_in.triggers["raw"]
  status           = bdf_in.status
  num_data_records = bdf_in.header["num_data_records"]
  num_samples      = bdf_in.header["num_samples"][1]
  num_channels     = bdf_in.header["num_channels"]

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

  bdf = Array{UInt8}(undef, 3*(num_data_records*num_channels*num_samples))
  pos = 1
  for rec = 0:(num_data_records-1)
    for chan = 1:num_channels
      if chan < num_channels
        for samp = 1:num_samples
          data_val   = data[chan, rec*num_samples + samp]
          bdf[pos  ] = (data_val % UInt8)
          bdf[pos+1] = ((data_val >> 8) % UInt8)
          bdf[pos+2] = ((data_val >> 16) % UInt8)
          pos += 3
        end
      else  # last channel is Status channel
        for samp = 1:num_samples
          trig_val   = trigs[rec*num_samples + samp]
          status_val = status[rec*num_samples + samp]
          bdf[pos  ] = trig_val % UInt8
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
merge_bdf(bdf_in::Array{BioSemiData}, filename::String="merged.bdf")
Merge BioSemiRaw structs to single BioSemiRaw struct. Checks that the
input BioSemiRaw structs have the same number of channels, same channel
labels and that each channel has the same sample rate.
### Examples:
```julia
dat1 = read_bdf("filename1.bdf")
dat2 = read_bdf("filename2.bdf")
dat3 = merge_bdf([dat1, dat2], "filename_merged.bdf")
```
"""
function merge_bdf(bdf_in::Array{BioSemiData}, filename::String)

  # check data structs to merge have same number of channels, channel labels + sample rate
  num_chans = (x -> x.header["num_channels"]).(bdf_in)
  !all(x -> x == num_chans[1], num_chans) && error("Different number of channels in bdf_in")
  chan_labels = (x -> x.header["channel_labels"]).(bdf_in)
  !all(y -> y == chan_labels[1], chan_labels) && error("Different channel labels bdf_in")
  sample_rate = (x -> x.header["sample_rate"]).(bdf_in)
  !all(y -> y == sample_rate[1], sample_rate) && error("Different sample rate in bdf_in")

  # make copy so that bdf_in is not altered
  bdf_out = deepcopy(bdf_in[1])
  bdf_out.header["filename"] = filename
  bdf_out.header["num_data_records"] = sum((x -> x.header["num_data_records"]).(bdf_in))

  # merged dat_chan Matrix (channels x samples)
  bdf_out.data = hcat((x -> x.data).(bdf_in)...)

  # merged time ans status array
  bdf_out.time   = collect(0:size(bdf_out.data, 2) -1) / bdf_in[1].header["sample_rate"][1]
  bdf_out.status = vcat((x -> x.status).(bdf_in)...)

  # merged triggers dict with offset idx
  for bdf in bdf_in[2:end]
    idx_offset = size(bdf.data, 2)
    bdf_out.triggers["idx"] = vcat(bdf_out.triggers["idx"], bdf.triggers["idx"] .+ idx_offset)
  end
  bdf_out.triggers["raw"]   = vcat((x -> x.triggers["raw"]).(bdf_in)...)
  bdf_out.triggers["val"]   = vcat((x -> x.triggers["val"]).(bdf_in)...)
  bdf_out.triggers["count"] = sort(countmap(bdf_out.triggers["val"]))

  return bdf_out

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
function select_channels_bdf(bdf_in::BioSemiData, channels::Union{Array{Int}, Array{String}})
  bdf_out = deepcopy(bdf_in)
  channels = channel_idx(bdf_out.header["channel_labels"], channels)
  update_header_bdf!(bdf_out.header, channels)
  bdf_out.data = bdf_out.data[channels[1:end-1], :]
  return bdf_out
end


"""
crop_bdf(bdf_in::BioSemiData, crop_type::tString, val::Array{Int}, filename::String)
Recuce the length of the recorded data. The border upon which to crop the bdf file can be defined using either
a start and end trigger ("triggers") or a start and end record ("records").
### Examples:
```julia
dat1 = read_bdf("filename1.bdf")
dat2 = crop_bdf(dat1, "triggers", [1 2],   "filename1_cropped.bdf") # between first trigger 1 and last trigger 2
dat3 = crop_bdf(dat1, "records",  [1 100], "filename1_cropped.bdf") # data records 1 to 100 inclusive
```
"""
function crop_bdf(bdf_in::BioSemiData, crop_type::String, val::Array{Int}, filename::String)

  length(val) != 2 && error("val should be of length 2")

  sample_rate = bdf_in.header["sample_rate"][1]
  if crop_type == "triggers"

    # find trigger value index
    trigStart = findfirst(bdf_in.triggers["val"] .== val[1])
    trigEnd   = findlast(bdf_in.triggers["val"]  .== val[2])
    idxStart  = bdf_in.triggers["idx"][trigStart]
    idxEnd    = bdf_in.triggers["idx"][trigEnd]

    # need to find boundardy equal to record breaks
    borders  = collect(1:sample_rate:size(bdf_in.data, 2))
    idxStart =  findfirst(borders .>= idxStart) * sample_rate
    idxEnd   = (findlast(borders  .<= idxEnd)   * sample_rate) - 1

  elseif crop_type == "records"

    # find trigger value index
    idxStart  = ((val[1]-1) * sample_rate) + 1
    idxEnd    =  (val[2]    * sample_rate)
    trigStart = findfirst(bdf_in.triggers["idx"] .>= idxStart)
    trigEnd   = findlast(bdf_in.triggers["idx"]  .<= idxEnd)

  else
    error("crop_type not recognized!")
  end

  # copy data and crop
  bdf_out = deepcopy(bdf_in)
  bdf_out.header["filename"] = filename
  bdf_out.header["num_data_records"] = Int(((idxEnd - idxStart)+1) / sample_rate)
  bdf_out.data = bdf_out.data[:, idxStart:idxEnd]
  bdf_out.time = collect(0:size(bdf_out.data, 2) -1) / bdf_out.header["sample_rate"][1]

  # update triggers
  bdf_out.triggers["raw"]   = bdf_out.triggers["raw"][idxStart:idxEnd]
  bdf_out.triggers["idx"]   = bdf_out.triggers["idx"][trigStart:trigEnd]
  bdf_out.triggers["val"]   = bdf_out.triggers["val"][trigStart:trigEnd]
  bdf_out.triggers["count"] = sort(countmap(bdf_out.triggers["val"]))

  return bdf_out

end


"""
downsample_bdf(bdf_in::BioSemiData, dec::Int, filename::String)
Reduce the sampling rate within a BioSemiData struct by an integer factor (dec).
### Examples:
```julia
dat1 = read_bdf("filename1.bdf")
dat2 = downsample_bdf(dat1, 2, "filename1_downsampled.bdf")
```
"""
# TO DO: make separate mirror function
function downsample_bdf(bdf_in::BioSemiData, dec::Int, filename::String)

  !ispow2(dec) && error("dec should be power of 2!")

  nsamp = dec * 20  # enough samples?
  ndec = div(nsamp, dec)

  bdf_out = deepcopy(bdf_in)
  data = Matrix{Float32}(undef, size(bdf_out.data, 1), div(size(bdf_out.data, 2), dec))
  for i in 1:size(bdf_out.data, 1)
    tmp_dat    = resample([reverse(bdf_out.data[i, 1:nsamp]); bdf_out.data[i, :]; reverse(bdf_out.data[i, end-(nsamp-1):end])], 1/dec)
    data[i, :] = convert(Array{Float32}, tmp_dat[ndec+1:end-ndec])
  end
  bdf_out.data = data

  bdf_out.header["sample_rate"] = div.(bdf_out.header["sample_rate"], dec)
  bdf_out.header["num_samples"] = div.(bdf_out.header["num_samples"], dec)
  bdf_out.time = collect(0:size(bdf_out.data, 2) -1) / bdf_out.header["sample_rate"][1]

  # update triggers
  bdf_out.triggers["raw"] = zeros(Int16, 1, size(bdf_out.data, 2))
  bdf_out.triggers["idx"] = convert(Array{Int64}, round.(bdf_out.triggers["idx"] / dec))
  bdf_out.triggers["raw"][bdf_out.triggers["idx"]] = bdf_out.triggers["val"]

  return bdf_out

end

"""
update_header_bdf(hd::Dict, channels::Array{Int})
Updates header Dict within BioSemiData struct following the selection
of specific channels in read_bdf or select_channels_bdf.
"""
function update_header_bdf!(hd::Dict, channels::Array{Int})
  hd["num_channels"]     = length(channels)
  hd["physical_min"]     = hd["physical_min"][channels]
  hd["physical_max"]     = hd["physical_max"][channels]
  hd["digital_min"]      = hd["digital_min"][channels]
  hd["digital_max"]      = hd["digital_max"][channels]
  hd["scale_factor"]     = hd["scale_factor"][channels]
  hd["transducer_type"]  = hd["transducer_type"][channels]
  hd["num_samples"]      = hd["num_samples"][channels]
  hd["channel_unit"]     = hd["channel_unit"][channels]
  hd["reserved"]         = hd["reserved"][channels]
  hd["sample_rate"]      = hd["sample_rate"][channels]
  hd["channel_labels"]   = hd["channel_labels"][channels]
  hd["pre_filter"]       = hd["pre_filter"][channels]
  hd["num_bytes_header"] = (length(channels)+1) * 256
end

"""
channel_idx(labels::Array{String}, channels::Array{String})
Return channel index given labels and desired selection.
"""
function channel_idx(labels::Array{String}, channels::Array{String})
  channels = [findfirst(x .== labels) for x in channels]
  any(channels .== nothing) && error("A requested channel label is not in the bdf file!")
  println("Selecting channels:", labels[channels])
  return unique(append!(channels, length(labels)))
end

"""
channel_idx(labels::Array{String}, channels::Array{Int})
Return channel index given labels and desired selection.
"""
function channel_idx(labels::Array{String}, channels::Array{Int})
  any(channels .> length(labels)) && error("Requested channel number greater than number of channels in file!")
  any(channels .< 1)              && error("Requested channel number less than 1!")
  println("Selecting channels:", labels[channels])
  return unique(append!(channels, length(labels)))
end

end
