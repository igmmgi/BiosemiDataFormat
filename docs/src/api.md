# API Reference

## Module

```@docs
BiosemiDataFormat
```

## Data Structures

```@docs
BiosemiDataFormat.BiosemiHeader
BiosemiDataFormat.BiosemiTriggers
BiosemiDataFormat.BiosemiData
```

## File I/O Functions

### Reading BDF Files

```@docs
BiosemiDataFormat.read_bdf
```

### Writing BDF Files

```@docs
BiosemiDataFormat.write_bdf
```

## Data Processing Functions

### Cropping and Slicing

```@docs
BiosemiDataFormat.crop_bdf
BiosemiDataFormat.crop_bdf!
BiosemiDataFormat.time_range
BiosemiDataFormat.trigger_info
```

### Downsampling

```@docs
BiosemiDataFormat.downsample_bdf
BiosemiDataFormat.downsample_bdf!
```

### Merging

```@docs
BiosemiDataFormat.merge_bdf
```

## Channel Management Functions

### Channel Selection

```@docs
BiosemiDataFormat.select_channels_bdf
BiosemiDataFormat.select_channels_bdf!
```

### Channel Deletion

```@docs
BiosemiDataFormat.delete_channels_bdf
BiosemiDataFormat.delete_channels_bdf!
```

### Channel Utilities

```@docs
BiosemiDataFormat.channel_index
```
