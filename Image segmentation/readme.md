# Cilia Image Processing Pipeline

This Python script provides an automated image processing pipeline for analyzing cilia structures in microscopy images. The pipeline performs segmentation of nuclei, cilia, and basal bodies, then combines the results into a multi-channel output image.

## Overview

The script processes multi-channel TIFF images containing:
- Nucleic acid staining (nuclei)
- Protein staining 
- Cilia structures
- Basal body structures

It applies segmentation algorithms to identify and mask these cellular components, then outputs a processed multi-channel image with both original and segmented data.

## Features

- **Multi-channel image processing**: Handles 4-channel input images
- **Automated segmentation**: Segments nuclei, cilia, and basal bodies
- **Cilia refinement**: Improves cilia segmentation using basal body information
- **ImageJ compatibility**: Outputs images in ImageJ-compatible format

## Requirements

The script requires the following Python packages:
- numpy
- tifffile
- pathlib
- os
- time (optional, currently commented out)

## Usage

### Basic Usage

```python
python cilia_processing.py
```

### Configuration

Modify the following variables in the script:

```python
local_folder = "/cilia_images/"  # Input directory containing TIFF files
save_path = Path('/save_dir/')   # Output directory for processed images
```

### Input Requirements

- Input images must be multi-channel TIFF files (`.tif` extension)
- Images should contain 4 channels in the expected order:
  1. Nucleic acid staining
  2. Protein staining
  3. Cilia structures
  4. Basal body structures

### Output

For each input image `image.tif`, the script generates:
- `image_processed.tif` - Multi-channel output with 7 channels:
  1. Cilia segmentation mask
  2. Basal body segmentation mask
  3. Nuclei segmentation mask
  4. Original protein channel
  5. Original cilia channel
  6. Original basal body channel
  7. Original nucleic acid channel

## File Structure

```
project/
├── cilia_processing.py
├── errors.txt (generated during processing)
├── /cilia_images/ (input directory)
│   ├── image1.tif
│   ├── image2.tif
│   └── ...
└── /save_dir/ (output directory)
    ├── image1_processed.tif
    ├── image2_processed.tif
    └── ...
```

## Notes

- The script automatically creates the output directory if it doesn't exist
- Processing time monitoring is available but currently commented out
- Output images are saved in ImageJ-compatible format with proper metadata
- All segmentation masks are converted to 8-bit format (0-255 range)
