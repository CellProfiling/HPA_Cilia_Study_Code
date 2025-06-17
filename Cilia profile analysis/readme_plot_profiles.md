# Cilia Intensity Profile Analysis

A Python script for analyzing and visualizing cilia intensity profiles from microscopy data.

## Description

This script processes cilia intensity measurements from text files to generate normalized intensity profiles along the length of cilia. The analysis includes:

- Reading intensity measurements from tab-separated files
- Normalizing intensity profiles by subtracting background thresholds
- Interpolating data points along standardized length grids
- Generating heatmaps to visualize intensity patterns across multiple cilia

## Features

- **Data Processing**: Handles missing values (marked as "?") in intensity and length measurements
- **Length Partitioning**: Splits cilia into two segments at 1.16 microns for detailed analysis
- **Intensity Normalization**: Subtracts background thresholds and normalizes to maximum values
- **Gaussian Smoothing**: Applies smoothing filter to reduce noise in profiles
- **Visualization**: Creates publication-ready heatmaps using seaborn

## Requirements

```
numpy
pandas
matplotlib
seaborn
scipy
```

## Installation

Install the required dependencies:

```bash
pip install numpy pandas matplotlib seaborn scipy
```

## File Structure

```
project/
├── script.py                    # Main analysis script
├── GPR75_002/                   # Directory containing cilia measurement files
│   ├── *_CQl.txt               # Individual cilia measurement files
├── GPR75_002_AllCQs.txt        # Threshold data file
└── GPR75_002.svg               # Generated heatmap output
```

## Input Data Format

### Cilia Measurement Files (`*_CQl.txt`)
Tab-separated files with columns:
- `ID`: Unique identifier for each cilium
- `Arc length [micron]`: Distance along cilium
- `Intensity A`: Fluorescence intensity values

### Threshold File (`GPR75_002_AllCQs.txt`)
Tab-separated file with columns:
- `File name`: Original image filename
- `ID`: Cilium identifier
- `Intensity threshold A`: Background threshold value

## Usage

1. **Prepare your data**: Ensure measurement files are in the `GPR75_002/` directory
2. **Update file paths**: Modify the script to point to your data directory and threshold file
3. **Run the analysis**:

```bash
python script.py
```

## Output

The script generates:
- **Console output**: Statistics on cilia with and without profiles
- **SVG heatmap**: `GPR75_002.svg` showing normalized intensity profiles

## Key Parameters

- **Length partition**: 1.16 microns (separates proximal and distal regions)
- **Grid resolution**: 20 points (proximal) + 80 points (distal) = 100 total
- **Smoothing sigma**: 5 (Gaussian filter parameter)
- **Colormap**: "magma" for heatmap visualization

## Functions

### `intensity_profile(file_path, df2)`
Processes a single cilia measurement file and returns intensity profiles.

**Parameters:**
- `file_path`: Path to the measurement file
- `df2`: DataFrame containing threshold data

**Returns:**
- `cilia_intensity`: Dictionary of interpolated intensity profiles
- `thresholds_local`: Dictionary of threshold values for each cilium

## Data Processing Steps

1. Load measurement and threshold data
2. Group measurements by cilium ID
3. Handle missing values and data type conversion
4. Interpolate intensity values along standardized length grids
5. Apply background subtraction using thresholds
6. Normalize profiles to common scale
7. Apply Gaussian smoothing
8. Generate and save heatmap visualization

## Troubleshooting

- **Missing threshold data**: Script will use threshold = 0 and print warning
- **Empty length segments**: Distal region filled with zeros if no data beyond 1.16μm
- **File encoding**: Uses ISO-8859-1 encoding for compatibility with international characters


## License

This script is provided as-is for research purposes.