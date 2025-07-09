"""
Processes multi-channel microscopy images to extract nuclei, cilia, and protein channels for cell cycle prediction.
"""

import numpy as np
import os
from skimage import filters
from scipy import ndimage as ndi
from skimage.morphology import closing, ball
from skimage.measure import label, regionprops
from PIL import Image
from tifffile import imread, imwrite


def load_image_stack(file_path):
    """Load multi-frame TIFF image into numpy array."""
    img = imread(file_path)
    return np.array(img)


def normalize_image(img):
    """Normalize image to 0-255 range."""
    img = img / 255
    img = img / np.max(img) * 255 if np.max(img) > 0 else img
    return img.astype(np.uint8)


def process_nuclei_segmentation(c1_stack):
    """Segment nuclei and find optimal z-plane."""
    z_length = c1_stack.shape[0]
    
    # Smooth and threshold
    smooth = filters.gaussian(c1_stack, sigma=3)
    thresh = filters.threshold_otsu(smooth)
    binary = smooth > thresh
    binary = closing(binary, ball(3))
    
    # Fill holes slice by slice
    nucleus_seg = []
    for i in range(z_length):
        nucleus_seg.append(ndi.binary_fill_holes(binary[i]))
    nucleus_seg = np.array(nucleus_seg).astype(np.uint8)
    
    # Label the binary image
    nucleus_seg = label(nucleus_seg)
    
    # Find optimal z-plane based on weighted centroids
    nucleus_props = regionprops(nucleus_seg, intensity_image=c1_stack)
    z_planes = []
    for prop in nucleus_props:
        z_planes.append(prop.weighted_centroid[0])
    
    # Remove nan values and calculate average
    z_planes = [z for z in z_planes if not np.isnan(z)]
    z_avg = round(np.mean(z_planes)) if z_planes else z_length // 2
    
    return z_avg


def main():
    """Main processing function."""
    GENE_NAME = "CKAP2"
    image_dir = f"./data/{GENE_NAME}"
    n_slices = 21  # number of slices per channel
    
    for file in os.listdir(image_dir):
        print(f"Processing {file}")
        
        # Load image stack
        img_stack = load_image_stack(os.path.join(image_dir, file))
        
        # Split channels
        c1_stack = img_stack[:n_slices, :, :]        # cell nuclei
        c2_stack = img_stack[n_slices:n_slices*2, :, :]    # cilia
        c3_stack = img_stack[n_slices*2:n_slices*3, :, :]  # not used
        c4_stack = img_stack[n_slices*3:n_slices*4, :, :]  # protein
        c5_stack = img_stack[n_slices*4:n_slices*5, :, :]  # basal bodies
        
        # Process nuclei and find optimal z-plane
        z_avg = process_nuclei_segmentation(c1_stack)
        
        # Extract images at optimal z-plane
        z_avg_img = normalize_image(c1_stack[z_avg, :, :])
        c4_avg_img = normalize_image(c4_stack[z_avg, :, :])
        
        # Create maximum projections
        c2_max = normalize_image(np.max(c2_stack, axis=0))
        c5_max = normalize_image(np.max(c5_stack, axis=0))
        
        # Print intensity statistics
        print(f"Max intensities - cilia: {c2_max.max()}, basal: {c5_max.max()}, "
              f"nuclei: {z_avg_img.max()}, protein: {c4_avg_img.max()}")
        
        # Create output directory
        output_dir = os.path.join(f"{GENE_NAME}/test", file.split(".")[0])
        os.makedirs(output_dir, exist_ok=True)
        
        # Save processed images
        Image.fromarray(z_avg_img).save(os.path.join(output_dir, "nuclei.png"))
        Image.fromarray(c2_max).save(os.path.join(output_dir, "cilia.png"))
        Image.fromarray(c5_max).save(os.path.join(output_dir, "basal.png"))
        Image.fromarray(c4_avg_img).save(os.path.join(output_dir, "protein.png"))
        
    print("Processing complete")


if __name__ == "__main__":
    main()