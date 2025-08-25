# Examples

## Basic Usage

### Reading and Writing BDF Files

```julia
using BiosemiDataFormat

# Read a BDF file
dat = read_bdf("eeg_data.bdf")

# Check basic information
println("Number of channels: ", dat.header.num_channels)
println("Sampling rate: ", dat.header.sample_rate[1], " Hz")
println("Duration: ", dat.header.duration_data_records, " seconds")

# Write the data to a new file
write_bdf(dat, "output.bdf")
```

### Channel Selection

```julia
# Select specific channels by name
channels = ["Fp1", "Fp2", "Cz", "O1", "O2"]
dat_selected = select_channels_bdf(dat, channels)

# Select channels by index
dat_selected = select_channels_bdf(dat, [1, 3, 5, 7])

# Mix of names and indices
dat_selected = select_channels_bdf(dat, [1, "Cz", 5])
```

### Data Cropping

```julia
# Crop by time range (in seconds)
dat_cropped = crop_bdf(dat, "time", [10.0, 20.0])

# Crop by data records
dat_cropped = crop_bdf(dat, "records", [100, 200])

# Crop by trigger values
dat_cropped = crop_bdf(dat, "triggers", [100, 200])
```

### Downsampling

```julia
# Downsample by factor of 2
dat_downsampled = downsample_bdf(dat, 2)

# Downsample by factor of 4
dat_downsampled = downsample_bdf(dat, 4)

# Check the new sampling rate
println("Original sampling rate: ", dat.header.sample_rate[1], " Hz")
println("New sampling rate: ", dat_downsampled.header.sample_rate[1], " Hz")
```

### Merging Multiple Files

```julia
# Merge two BDF files
file1 = "session1.bdf"
file2 = "session2.bdf"

dat1 = read_bdf(file1)
dat2 = read_bdf(file2)

dat_merged = merge_bdf([dat1, dat2])
```

### Channel Manipulation

```julia
# Delete specific channels
dat_reduced = delete_channels_bdf(dat, ["Fp1", "Fp2"])

# Check remaining channels
println("Remaining channels: ", dat_reduced.header.channel_labels)

```


