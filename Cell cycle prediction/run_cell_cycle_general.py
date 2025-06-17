import os
from utils import allocate_gpu

try:
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("Using device: " + device)
except:
    print("No GPU available, will use CPU")
    device = "cpu"


from configs import get_default_config
from helpers import (
    start_training,
    start_regression_inference,
    start_regression_inference_cellcycle,
)

GENE_NAME = "CKAP2"

config = get_default_config()
config.phase = "test"
config.architecture = "DPNUNet" #"VQVAE" #"DPNUNet" #"AE" #"DPNUNet"
#config.folder = "/data/cellcycle/DPN-training/Fucci-dataset-v3-minaug"
config.folder = "{}_output".format(GENE_NAME)
config.dataset_path = GENE_NAME
#config.dataset_path = "/data/cellcycle/datasets/Fucci-dataset-v3_filtered"
config.outputs_dir = os.path.join(config.dataset_path, config.folder, "outputs")
config.weights_dir = os.path.join(config.dataset_path, config.folder, "weights")
config.logs_dir = os.path.join(config.dataset_path, config.folder, "logs")
config.load_model_from = r"weights/dpn_unet_best.pth"
# config.load_model_from = None
config.nb_epoch = 10001
config.batch_size = 16
config.epoch_size = 1
config.target_cols = 512
config.target_rows = 512
config.lr = 5e-4
config.ignore_target_size = True
config.loss = {"type":"mse"}# {"type": "mse+bce", "bce":0.4, "dice":0.6}
config.activation = "linear"
config.input_channel_files = ["cilia.png", "basal.png", "nuclei.png"]
config.target_channel_files = [None, "protein.png", None]
config.input_channel_num = 3
config.target_channel_num = 3
config.scale_factor = 0.5
config.num_workers = 0
config.save_interval = 10
config.training_augmentations = "victor"#"minimal" #"victor"
config.inputs_preprocessing = "bit-depth" #"min-max" #  # "min-max"
config.targets_preprocessing = "rescale99"#"equalize-adapthist"#"bit-depth"
config.test_pad = 0
config.device = device

if config.phase == "train":
    start_training(config)
else:
    start_regression_inference_cellcycle(config)
