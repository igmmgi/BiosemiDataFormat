# BiosemiDataFormat

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igmmgi.github.io/BiosemiDataFormat.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igmmgi.github.io/BiosemiDataFormat.jl/dev/)
[![Build Status](https://github.com/igmmgi/BiosemiDataFormat.jl/workflows/Documentation/badge.svg)](https://github.com/igmmgi/BiosemiDataFormat.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Julia code for BioSemi 24 bit EEG files. The code can be used for:

- reading files into Julia data struct
- cropping file length
- reducing the sample rate
- selecting/reducing the number of channels
- writing a Julia data struct to a bdf fileformat

## Installation

```julia
] # julia pkg manager
add https://github.com/igmmgi/BiosemiDataFormat.git # install from  GitHub
test BiosemiDataFormat # optional
```

## Functions

- crop_bdf
- downsample_bdf
- merge_bdf
- select_channels_bdf
- read_bdf
- write_bdf

## Basic Example

```julia
using BiosemiDataFormat

dat1 = read_bdf("filename1.bdf")
dat2 = read_bdf("filename2.bdf")
dat3 = merge_bdf([dat1, dat2])
write_bdf(dat3, "filename3.bdf")

```
