from tifffile import imread, imwrite
import numpy as np
from cilia_utils import *
# from HPA_segmentor_test import hpa_segmentor
from pathlib import Path
import os


def stack_(nucleic_img_stack, protein_img_stack, cilia_img_stack, bb_img_stack, nucleic_seg_stack, cilia_seg_stack,
           bb_seg_stack):
    # stack the images
    new_stack = np.stack(
        [cilia_seg_stack, bb_seg_stack, nucleic_seg_stack, protein_img_stack, cilia_img_stack, bb_img_stack,
         nucleic_img_stack], axis=-1)

    return new_stack

def main(local_folder, img_name):
    # start = time.time()
    
    # img_id = img_id[:-1]
    save_path = Path('/save_dir/')
    if not save_path.exists():
        save_path.mkdir()


    nucleic_img_stack, protein_img_stack, cilia_img_stack, bb_img_stack = read_img_stack(os.path.join(local_folder, img_name))
    #-----------------------------------------------------------------------------------------
    # Run the pipeline
    nucleic_mask = nucleic_segmentation(nucleic_img_stack, save_path)

    nucleic_mask = nucleic_mask.astype(np.uint8)
    cilia_mask = cilia_segmentation(cilia_img_stack, save_path)
    bb_mask = basal_body_segmentation(bb_img_stack, save_path)
    bb_mask = 255 * bb_mask.astype(np.uint8)

    cilia_mask = cilia_refine(cilia_mask, bb_mask)
    cilia_mask = 255 * cilia_mask.astype(np.uint8)

    # check returns are not none
    try:
        new_stack = stack_(nucleic_img_stack, protein_img_stack, cilia_img_stack, bb_img_stack, nucleic_mask, cilia_mask, bb_mask)
        # move the channel axis to the second dimension
        new_stack = np.moveaxis(new_stack, 3, 1)
        imwrite(os.path.join(save_path, img_name.replace(".tif", "_processed.tif")), new_stack, metadata={'axes': 'ZCYX'},
                imagej=True)
        print("Image processed: ", img_name)
        return True
    except Exception as e:
        with open('errors.txt', 'a') as f:
            f.write(img_name + ":  images don't have the same dimensions")
            f.write('\n')
        return False

if __name__ == "__main__":
    local_folder = "/cilia_images/"

    for imfile in os.listdir(local_folder):
        if imfile.endswith(".tif"):
            main(local_folder, imfile)
