# Deep Learning Pipeline for Cell Cycle Analysis

This repository contains a collection of Python scripts for training and running deep learning models on cell cycle analysis data, with specialized utilities for cellular image segmentation and analysis.

## Overview

The pipeline is designed for training neural networks (particularly DPN-UNet architectures) on cellular imaging data, with specific focus on:
- Cell cycle phase prediction
- Cellular component segmentation (nuclei, cilia, basal bodies)
## Scripts Description

### Core Scripts

#### `run_cell_cycle_general.py`
**Main execution script** for running the cell cycle analysis pipeline.
- Configures and launches training or inference workflows
- Supports multiple model architectures (DPNUNet, VQVAE, AE)
- Handles GPU/CPU device allocation
- Configurable for different gene datasets (default: CKAP2)

#### `callbacks.py`
**Training callbacks and utilities** for model training workflow.
- `ModelSaver`: Saves best models during training
- `CheckpointSaver`: Regular checkpoint saving functionality  
- `TensorBoard`: Logging training metrics to TensorBoard
- `PredictionSaver`: Saves prediction outputs during training
- `HistorySaver`: Tracks and visualizes training history
- `ModelFreezer`: Implements gradual unfreezing of model layers

#### `configs.py` & `configs_parser.py`
**Configuration management** for training parameters.
- Default configuration settings for model training
- Command-line argument parsing
- Hyperparameter definitions (learning rate, batch size, etc.)
- Dataset path and channel configuration

#### `utils.py`
**Utility functions** for image processing and model management.
- Image preprocessing functions (normalization, scaling)
- Model loading and GPU allocation
- Watershed segmentation utilities
- File I/O operations for training data


## Installation

### Requirements
```bash
pip install torch torchvision
pip install opencv-python
pip install scikit-image
pip install tifffile
pip install tensorboardX
pip install imageio
pip install matplotlib
pip install scipy
pip install namedlist
pip install GPUtil
```

### Additional Dependencies
- `imgseg` (for GeoJSON mask generation)
- CUDA-compatible PyTorch for GPU acceleration

## Usage

### Training a Model
```bash
python run_cell_cycle_general.py
```

### Configuration
Modify the configuration in `run_cell_cycle_general.py`:

```python
# Example configuration
config.phase = "train"  # or "test"
config.architecture = "DPNUNet"
config.dataset_path = "path/to/your/dataset"
config.input_channel_files = ["cilia.png", "basal.png", "nuclei.png"]
config.target_channel_files = [None, "protein.png", None]
config.batch_size = 16
config.lr = 5e-4
```

### Input Data Format
The pipeline expects datasets organized as:
```
dataset_path/
├── train/
│   ├── sample1/
│   │   ├── cilia.png
│   │   ├── basal.png
│   │   ├── nuclei.png
│   │   └── protein.png
│   └── sample2/
│       └── ...
└── valid/
    └── ...
```


## Model Architectures

### Supported Models
- **DPNUNet**: Dual Path Network with U-Net architecture
- **VQVAE**: Vector Quantized Variational Autoencoder
- **AE**: Standard Autoencoder

### Loss Functions
- Binary Cross Entropy (BCE)
- Dice Loss
- Mean Squared Error (MSE)
- Combined loss functions

## Output Structure
```
output_directory/
├── weights/           # Model checkpoints
├── logs/             # Training logs and TensorBoard files
├── outputs/          # Prediction outputs and visualizations
└── history.pkl       # Training history data
```

## Key Features

### Training Features
- **Multi-GPU support** with automatic GPU allocation
- **Mixed precision training** support
- **Gradual unfreezing** for transfer learning
- **Comprehensive logging** with TensorBoard integration
- **Automatic checkpointing** and best model saving

### Image Processing Features
- **Multi-channel support** (up to 4 channels)
- **3D image stack processing** for Z-stack microscopy
- **Advanced segmentation** algorithms for cellular components
- **Flexible preprocessing** pipelines

### Analysis Capabilities
- **Cell cycle phase prediction**
- **Protein localization analysis**
- **Cellular component segmentation**
- **Temporal tracking** across image sequences

## Configuration Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `batch_size` | Training batch size | 16 |
| `lr` | Learning rate | 5e-4 |
| `nb_epoch` | Number of training epochs | 20 |
| `target_rows/cols` | Input image dimensions | 512 |
| `input_channel_num` | Number of input channels | 3 |
| `scale_factor` | Image scaling factor | 1.0 |
| `warmup` | Epochs for gradual unfreezing | 0 |

## Custom Preprocessing
```python
from utils import get_preprocessing

# Available preprocessing options:
# - "identity": No preprocessing
# - "min-max": Min-max normalization
# - "bit-depth": Bit depth normalization
# - "rescale99": 99th percentile rescaling
# - "equalize-adapthist": Adaptive histogram equalization

preprocessor = get_preprocessing("min-max")
```

### Custom Callbacks
```python
from callbacks import Callback

class CustomCallback(Callback):
    def on_epoch_end(self, epoch):
        # Custom logic here
        pass
```

## License

This project is designed for research purposes.
