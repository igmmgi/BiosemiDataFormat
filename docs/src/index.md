# BiosemiDataFormat

A Julia package for reading, writing, and processing BioSemi 24-bit EEG data files (BDF format).

## Overview

BioSemi BDF (BioSemi Data Format) files store 24-bit EEG data with metadata including channel information, sampling rates, and trigger events. This package provides comprehensive functionality to work with these files in Julia.

## Features

- **File I/O**: Read and write Biosemi (\*.bdf) files
- **Basic Data Processing**: Crop, downsample, and merge data
- **Channel Management**: Select, delete, and manipulate channels
- **Trigger Analysis**: Extract and analyze trigger events
- **Status Channel**: Handle BioSemi status channel information

## Quick Start

```julia
using BiosemiDataFormat

# Read a BDF file
dat = read_bdf("eeg_data.bdf")

# Select specific channels
dat_selected = select_channels_bdf(dat, ["Fp1", "Cz", "O1"])

# Crop data to specific time range
dat_cropped = crop_bdf(dat, "triggers", [100, 500])

# Downsample data
dat_downsampled = downsample_bdf(dat, 2)

# Write modified data
write_bdf(dat_downsampled, "processed_data.bdf")
```

## Installation

```julia
using Pkg
Pkg.add("BiosemiDataFormat")
```

## Documentation

- [API Reference](@ref)
