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

bdf = read_bdf("filename.bdf")

```

