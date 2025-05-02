from tifffile import imread, imsave
import numpy as np
from skimage import filters
from pathlib import Path
from scipy import ndimage as ndi

from skimage.morphology import closing, ball, remove_small_objects, remove_small_holes, dilation, disk, ball, opening, closing
from skimage import exposure
from skimage.measure import label, regionprops
from skimage.filters import sato

import os

def imread_(path):
    img = imread(path)
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i, j].sum() > 0:
                return img[i, j]
            

def adjust_image(img, min_size, connectivity, sigma, z_length):
    binary_cleaned = remove_small_objects(img > 0, min_size=min_size, connectivity=connectivity)
    img_new = binary_cleaned.astype("int") * img
    img_new = filters.gaussian(img_new, sigma=sigma)

    p2, p90 = np.percentile(img_new, (2, 98))
    img_new = exposure.rescale_intensity(img_new, in_range=(p2, p90))

    # Equalization
    img_new = exposure.equalize_hist(img_new)

    # Adaptive Equalization
    img_new = exposure.equalize_adapthist(img_new, clip_limit=0.03)

    img_new = ndi.maximum_filter(img_new, size=1, mode='constant')
    return img_new


def nucleic_segmentation(nucleic_img_stack, save_dir):
    z_length = nucleic_img_stack.shape[0]

    smooth = filters.gaussian(nucleic_img_stack, sigma=3)

    thresh = filters.threshold_otsu(smooth)
    binary = smooth > 0.75 * thresh
    binary = closing(binary, ball(3))

    nucleus_seg = []
    for i in range(z_length):
        nucleus_seg.append(ndi.binary_fill_holes(binary[i]))
    nucleus_seg = np.array(nucleus_seg).astype(np.uint8)

    nucleus_seg = label(nucleus_seg)
    
    return nucleus_seg


def generate_er_mt(nucleic_img_stack, protein_img_stack, save_dir):
    z_length = nucleic_img_stack.shape[0]

    smooth = filters.gaussian(nucleic_img_stack, sigma=3)

    thresh = filters.threshold_otsu(smooth)
    binary = smooth > 0.75 * thresh
    binary = closing(binary, ball(3))

    nucleus_seg = []
    for i in range(z_length):
        nucleus_seg.append(ndi.binary_fill_holes(binary[i]))
    nucleus_seg = np.array(nucleus_seg).astype(np.uint8)

    ### cell segmentation
    img_stack = protein_img_stack.astype("f")

    # Adjust the image
    img_stack = adjust_image(img_stack, min_size=2, connectivity=3, sigma=2, z_length=z_length)

    # img_stack_a = img_stack_a[5]
    # img_stack_b = img_stack_b[5]

    thresh = filters.threshold_otsu(img_stack)
    binary = img_stack > thresh

    binary = remove_small_holes(binary, 1000)
    # find connected components
    label_objects, nb_labels = ndi.label(binary)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > 25000
    mask_sizes[0] = 0
    binary_cleaned = mask_sizes[label_objects]

    imsave(os.path.join(save_dir, "binary_cleaned.tif"), binary_cleaned.astype(np.uint8))

    nucleus = smooth
    nucleus = nucleus - np.min(nucleus)
    nucleus = nucleus / np.max(nucleus) * 233



    temp = img_stack * filters.gaussian(np.where(nucleus_seg == 1, 0.1, 1), sigma = 5, preserve_range=True)
    pseudo_mask = filters.gaussian(binary_cleaned, sigma = 5, preserve_range=True)
    temp *= pseudo_mask
    # smoothing in z
    temp = filters.gaussian(temp, sigma = [3, 0, 0], preserve_range=True)
    temp /= np.max(temp)
    temp *= 250


    mt = np.zeros_like(temp)
    for i in range(z_length):
        mt[i] = sato(temp[i], black_ridges = False, sigmas = [3])
    mt = filters.gaussian(mt, sigma = [3, 0, 0], preserve_range=True)
    mt /= np.max(mt)
    mt *= 250

    image_list = [[], [], []]
    # save per slice
    for i in range(z_length):
        imsave(os.path.join(save_dir, f"mt_" + str(i) + ".tif"), mt[i].astype(np.uint8))
        imsave(os.path.join(save_dir, f"er_" + str(i) + ".tif"), temp[i].astype(np.uint8))
        imsave(os.path.join(save_dir, f"nucleus_" + str(i) + ".tif"), nucleus[i].astype(np.uint8))
        image_list[0].append(os.path.join(save_dir, f"mt_" + str(i) + ".tif"))
        image_list[1].append(os.path.join(save_dir, f"er_" + str(i) + ".tif"))
        image_list[2].append(os.path.join(save_dir, f"nucleus_" + str(i) + ".tif"))

    return image_list


def refine_v2(nucleic_label, cell_label):
    z_length = nucleic_label.shape[0]

    cell_id_current = np.max(cell_label[0]) + 1

    for i in range(z_length - 1):
        cell_label_current = cell_label[i]
        cell_label_next = cell_label[i + 1]
       
        cell_label_next_binary = np.zeros_like(cell_label_next)

        for j in range(1, np.max(cell_label_next) + 1):
            overlaps = []
            for k in range(1, np.max(cell_label_current) + 1):
                region = cell_label_next == j
                region = region.astype("int")
                region_prev = cell_label_current == k
                overlap = np.sum(np.logical_and(region, region_prev))
                overlaps.append(overlap)
            if len(overlaps) == 0:
                continue
            min_index = np.argmax(overlaps)
            if overlaps[min_index] < 1000:
                continue

            region = cell_label_next == j
            region = region.astype("int")
            # opening
            region = dilation(region, disk(5))
            region_prev = cell_label_current == min_index + 1
            region_prev = region_prev.astype("int")
            region_prev = dilation(region_prev, disk(3))

            cell_label_next_binary[np.logical_and(region, region_prev) > 0] = min_index + 1

        # for j in range(1, np.max(cell_label_next) + 1):
        #     if j not in cell_label_next_binary:
        #         cell_label_next_binary[cell_label_next == j] = cell_id_current 
        #         cell_id_current += 1
        
        cell_label[i + 1] = cell_label_next_binary

    return nucleic_label, cell_label


def refine(nucleic_label, cell_label):
    z_length = nucleic_label.shape[0]

    for i in range(z_length - 1):
        cell_label_current = cell_label[i]
        cell_label_next = cell_label[i + 1]
        # get centers
        centers = []
        for j in range(1, np.max(cell_label_current) + 1):
            centers.append(np.mean(np.where(cell_label_current == j), axis=1))
        if len(centers) == 0:
            continue
        cell_label_next_binary = cell_label_next > 0
        cell_label_next_binary = cell_label_next_binary.astype("int")
        # relabel cell_label_next_binary by centers
        for j in range(1, np.max(cell_label_next) + 1):
            center = np.mean(np.where(cell_label_next == j), axis=1)
            distances = np.linalg.norm(np.array(centers) - center, axis=1)
            min_index = np.argmin(distances)

            region = cell_label_next == j
            region = region.astype("int")
            # opening
            region = dilation(region, disk(5))
            cell_label_next_binary[region > 0] = min_index + 1

        
        for j in range(1, np.max(cell_label_current) + 1):
            if j not in cell_label_next_binary:
                cell_label_next_binary[cell_label_current == j] = j

        # cell_label_next_binary = np.where(cell_label_next > 0, cell_label_next_binary, 0)

        k = 1
        # check for overlap
        for j in range(1, np.max(cell_label_next_binary) + 1):
            region = cell_label_next_binary == j
            nucleic_region = nucleic_label[i + 1] > 0
            if np.sum(region) < 1000 or np.sum(np.logical_and(region, nucleic_region)) < 10:
                cell_label_next_binary[region > 0] = 0
                continue
            region_2 = cell_label_current == j
            overlap = np.sum(np.logical_and(region, region_2))
            if overlap == 0:
                cell_label_next_binary[region > 0] = np.max(cell_label_next_binary) + k
                k += 1

        cell_label[i + 1] = cell_label_next_binary


    # check individual slices
    for i in range(1, z_length - 1):
        # identify components in the current slice
        cell_label_current = cell_label[i]
        # check overlap with the previous and next slices
        cell_label_previous = cell_label[i - 1]
        cell_label_next = cell_label[i + 1]
        for j in range(1, np.max(cell_label_current) + 1):
            region = cell_label_current == j
            region = region.astype("int")
            overlap = np.sum(np.logical_and(region, cell_label_previous > 0))

            overlap_next = np.sum(np.logical_and(region, cell_label_next > 0))
            if overlap < 10 or overlap_next < 10:
                cell_label_current[region > 0] = 0

    # check the last slice
    cell_label_current = cell_label[-1]
    cell_label_previous = cell_label[-2]
    for j in range(1, np.max(cell_label_current) + 1):
        region = cell_label_current == j
        region = region.astype("int")
        overlap = np.sum(np.logical_and(region, cell_label_previous > 0))
        if overlap < 10:
            cell_label_current[region > 0] = 0


        
    nucleic_label_binary = nucleic_label > 0
    nucleic_label = np.multiply(nucleic_label_binary > 0, cell_label)

    return nucleic_label, cell_label


def read_img_stack(data_path):
    data_path = Path(data_path)
    # save_path = Path(r'output')
    # if not save_path.exists():
    #     save_path.mkdir()

    # Load the images
    nucleic_img_stack = []
    for i in range(30):
        suffix = "*c0_z" + str(i) + ".tif"
        directory = data_path

        for i in directory.rglob(suffix):
            nucleic_img_stack.append(imread_(i))
    nucleic_img_stack = np.array(nucleic_img_stack)


    protein_img_stack = []
    for i in range(30):
        suffix = "*c1_z" + str(i) + ".tif"
        directory = data_path

        for i in directory.rglob(suffix):
            protein_img_stack.append(imread_(i))
    protein_img_stack = np.array(protein_img_stack)

    cilia_img_stack = []
    for i in range(30):
        suffix = "*c2_z" + str(i) + ".tif"
        directory = data_path

        for i in directory.rglob(suffix):
            cilia_img_stack.append(imread_(i))
    cilia_img_stack = np.array(cilia_img_stack)

    bb_img_stack = []
    for i in range(30):
        suffix = "*c3_z" + str(i) + ".tif"
        directory = data_path

        for i in directory.rglob(suffix):
            bb_img_stack.append(imread_(i))
    bb_img_stack = np.array(bb_img_stack)

    return nucleic_img_stack, protein_img_stack, cilia_img_stack, bb_img_stack


def delete_temp_files(save_dir):
    for f in os.listdir(save_dir):
        if "mt_" in f or "er_" in f or "nucleus_" in f or "binary_cleaned" in f:
            os.remove(os.path.join(save_dir, f))

def basal_body_segmentation(bb_img_stack, save_path):
    blurred = filters.gaussian(bb_img_stack, sigma=1) - filters.gaussian(bb_img_stack, sigma=10)
    # set negative values to 0
    blurred[blurred < 0] = 0
    # enhanced = sato(blurred, sigmas = [5])
    thresh = filters.threshold_triangle(blurred)
    # if thresh < 0.015:
    #     thresh = 0.015
    binary = blurred > 10 * thresh

    # image open
    binary = opening(binary, ball(1))

    # remove small objects
    binary = remove_small_objects(binary, min_size=75)

    # # save the mask
    # imsave(os.path.join(save_path, "bb_mask.tif"), 255 * binary.astype(np.uint8))
    return binary


def cilia_segmentation(cilia_img_stack, save_path):
    blurred = filters.gaussian(cilia_img_stack, sigma=1) - filters.gaussian(cilia_img_stack, sigma=5)
    # set negative values to 0
    blurred[blurred < 0] = 0
    # enhanced = sato(blurred, sigmas = [5])
    thresh = filters.threshold_otsu(blurred)
    if thresh < 0.025:
        thresh = 0.025
    binary = blurred > thresh

    # image open
    binary = opening(binary, ball(1))

    # remove small objects
    binary = remove_small_objects(binary, min_size=200)

    # # save the mask
    # imsave(os.path.join(save_path, "cilia_mask.tif"), 255 * binary.astype(np.uint8))

    return binary

def cilia_refine(cilia_mask, bb_mask):
    # label cilia
    label_cilia = label(cilia_mask)
    # get bounding box of cilia
    props = regionprops(label_cilia)
    margin = 15
    for prop in props:
        minz, minr, minc, maxz, maxr, maxc = prop.bbox
        minr = max(0, minr - margin)
        minc = max(0, minc - margin)
        maxr = min(label_cilia.shape[1], maxr + margin)
        maxc = min(label_cilia.shape[2], maxc + margin)



        sum_bb = np.sum(bb_mask[:, minr:maxr, minc:maxc])
        if sum_bb == 0:
            # this label is not a cilia
            cilia_mask[label_cilia == prop.label] = 0

    return cilia_mask


# if __name__ == "__main__":
#     nuclei_mask_stack = imread(r"E:\cellX_dev\output\nuclei_mask.tif")
#     cell_mask_stack = imread(r"E:\cellX_dev\output\cell_mask.tif")
#     nuclei_mask_stack, cell_mask_stack = refine(nuclei_mask_stack, cell_mask_stack)
#     imsave(r"E:\cellX_dev\output\nuclei_mask_refined.tif", nuclei_mask_stack)
#     imsave(r"E:\cellX_dev\output\cell_mask_refined.tif", cell_mask_stack)