# Cilia Intensity Profile Analysis

A Python script for processing and analyzing cilia intensity profiles from microscopy data. This tool extracts intensity measurements along cilia structures, normalizes spatial coordinates, and generates processed datasets for downstream analysis.

## Overview

This script processes microscopy data files containing cilia measurements, performing spatial interpolation and data aggregation to create standardized intensity profiles. It's designed to work with tab-delimited measurement files from image analysis software and generates pickle files containing processed data ready for visualization and statistical analysis.

## Features

- **Automated File Processing**: Batch processes multiple microscopy measurement files
- **Spatial Normalization**: Standardizes cilia measurements to fixed coordinate grids
- **Dual-Region Analysis**: Separates basal body region (0-1.16 μm) from ciliary shaft
- **Threshold Integration**: Applies background correction using intensity thresholds
- **Cell Line Annotation**: Links measurements to cell line information and experimental conditions
- **Data Export**: Saves processed data as pickle files and generates summary CSV

## Requirements

```python
pandas>=1.3.0
numpy>=1.20.0
webdav3>=3.14.0
pathlib
pickle
glob
csv
```

## Input Files

### Required Data Files
- **`.CQ/` directory**: Contains individual measurement files (`*_CQl.txt`)
- **`Filtered_Staining_List.xlsx`**: Excel file with plate layout and cell line information
- **`AllCQs.txt`**: Tab-delimited file containing intensity thresholds for background correction

### File Format Requirements
- Measurement files must contain columns: `ID`, `Arc length [micron]`, `Intensity A`
- Excel file must contain columns: `Plate id`, `Position`, `Cell line`, plus condition columns (positions 40+)
- Threshold file must contain: `File name`, `ID`, `Intensity threshold A`

## Processing Workflow

### 1. Data Loading and Validation
- Reads measurement files from `.CQ/` directory
- Matches files with plate layout information
- Validates data integrity and handles missing values

### 2. Cilia Segmentation
- **Basal Body Region**: 0 to 1.16 μm (20 interpolated points)
- **Ciliary Shaft**: 1.16 μm to maximum length (80 interpolated points)
- Creates 100-point normalized intensity profiles per cilium

### 3. Background Correction
- Retrieves appropriate intensity thresholds from reference file
- Matches thresholds by filename and cilia ID
- Handles missing threshold data gracefully

### 4. Data Aggregation
- Groups measurements by plate ID and well position
- Tracks cell line information and experimental conditions
- Counts basal bodies with and without cilia structures

## Output Files

### Pickle Files
- **`intensity_profiles_A_all.pkl`**: Dictionary of normalized intensity profiles
- **`cilia_lengths_all.pkl`**: Dictionary of maximum cilia lengths per measurement
- **`annotations_all.pkl`**: Dictionary of cell line and condition annotations
- **`thresholds_all.pkl`**: Dictionary of background intensity thresholds

### CSV Summary
- **`path_list.csv`**: Summary statistics including basal body counts per well

## Data Structure

### Dictionary Keys
All pickle files use consistent keys: `"PlateID_WellPosition"` (e.g., `"2167_B1"`)

### Profile Format
- Each intensity profile contains 100 normalized spatial points
- Values represent background-corrected fluorescence intensities
- Spatial coordinates span from basal body to ciliary tip

### Annotation Format
- Cell line information stored as strings
- Experimental conditions stored as lists of condition names
- Derived from columns 40+ in the Excel metadata file

## Usage Example

```python
import pickle
import pandas as pd

# Load processed data
intensity_profiles = pickle.load(open('intensity_profiles_A_all.pkl', 'rb'))
cilia_lengths = pickle.load(open('cilia_lengths_all.pkl', 'rb'))
annotations = pickle.load(open('annotations_all.pkl', 'rb'))

# Access data for specific well
well_key = "2167_B1"
profiles = intensity_profiles[well_key]
lengths = cilia_lengths[well_key]
cell_line_info = annotations[well_key]

print(f"Number of cilia in well {well_key}: {len(profiles)}")
print(f"Cell line: {cell_line_info[0] if cell_line_info else 'Unknown'}")
```

## Error Handling

- **Missing Files**: Continues processing when measurement files are missing
- **Data Validation**: Handles "?" placeholder values in measurement data
- **Threshold Matching**: Uses zero threshold when matching fails
- **Empty Regions**: Fills ciliary shaft with zeros when no measurements exist

## Troubleshooting

### Common Issues
- **File encoding errors**: Ensure measurement files use ISO-8859-1 encoding
- **Missing thresholds**: Check filename matching between measurement and threshold files
- **Excel formatting**: Verify plate ID and position columns contain expected data types

## License

This script is designed for research use. Please cite appropriately when used in publications.