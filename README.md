# BioSemiBDF
Julia code for BioSemi 24 bit EEG files. The code can be used for:

* reading bdf files into Julia data struct
* cropping the length of recorded bdf files
* reducing the sample rate of bdf files
* selecting/reducing the number of channels in a bdf file
* writing a Julia data struct to a bdf fileformat

## Installation
``` julia
] # julia pkg manager
add https://github.com/igmmgi/BioSemiBDF.git # install from  GitHub
test BioSemiBDF # optional
```

## Functions
* crop_bdf 
* downsample_bdf 
* merge_bdf
* select_channels_bdf
* read_bdf
* write_bdf

## Basic Example
``` julia
using BioSemiBDF

dat1 = read_bdf("filename1.bdf")
dat2 = read_bdf("filename2.bdf")
dat3 = merge_bdf([dat1, dat2], "filename3.bdf")
write_bdf(dat3)

```

