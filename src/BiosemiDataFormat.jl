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
dat = read_bdf("eeg_data.bdf")

# Select specific channels
dat_selected = select_channels_bdf(dat, ["Fp1", "Cz", "O1"])

# Crop data to specific time range
dat_cropped = crop_bdf(dat, "triggers", [100, 200])

# Write modified data
write_bdf(dat_cropped, "processed_data.bdf")
```

# License
This package is licensed under the MIT License.
"""
module BiosemiDataFormat

using DSP
using Logging
using OrderedCollections

# Include organized module files
include("types.jl")
include("processing.jl")
include("channels.jl")
include("io.jl")

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

end # module
