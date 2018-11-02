# BioSemiBDF
Julia code to read BioSemi 24 bit EEG files.

## Installation
``` julia
] # julia pkg manager
add https://github.com/igmmgi/BioSemiBDF.git # install from  GitHub
```

## Basic Example
``` julia
using BioSemiBDF

dat1 = read_bdf("filename1.bdf")
dat2 = read_bdf("filename2.bdf")
dat3 = merge_bdf([dat1, dat2], "filename3.bdf")
write_bdf(dat3)

```

