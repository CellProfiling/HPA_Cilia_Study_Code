import glob

from webdav3.client import Client
import os
import numpy as np
from pathlib import Path
import pandas as pd
import csv
import pickle




# Change this with your credentials
def intensity_profile(file_path, df2):
    # file_path = r"C:\Users\sunh\Downloads\2131_G11_88_new_stack_CQl.txt"

    # Read the file into a DataFrame, skipping unnecessary metadata lines
    df = pd.read_csv(file_path, sep='\t', header=0, encoding='ISO-8859-1')  # Assuming the metadata is on the first line
    n_BB_wo_cilia = 0
    n_BB_w_cilia = 0
    cilia_intensity = {}
    cilia_length = {}
    thresholds_local = {}

    for i in np.unique(df["ID"]):
        # get rows with the same ID
        temp = df[df["ID"] == i]
        if len(temp) == 1:
            n_BB_wo_cilia += 1
            continue
        n_BB_w_cilia += 1
        length = temp["Arc length  [micron]"].values
        # replace ? with 0
        length = np.where(length == "?", 0, length)
        length = length.astype(float)
        intensity = temp["Intensity A"].values
        # replace ? with 0
        intensity = np.where(intensity == "?", 0, intensity)
        intensity = intensity.astype(float)

        # get the threshold by matching Filename and ID
        row2 = df2[(df2["File name"] == filename) & (df2["ID"] == i)] 
        if len(row2) == 0:
            print("No matching row found in df2 for ID:", i)
            threshold = 0
        else:
            threshold = row2["Intensity threshold A"].values[0]
        if threshold == "Intensity threshold A":
            threshold = 0

        # partition the length into 2 parts by 1.16 micron
        length1 = length[length < 1.16]
        intensity1 = intensity[length < 1.16]
        grid1 = np.linspace(np.min(length1), 0.2, 20)
        intensity_interp1 = np.interp(grid1, length1, intensity1)
        length2 = length[length >= 1.16]
        if len(length2) == 0:
            intensity_interp2 = np.zeros(80)
        else:
            intensity2 = intensity[length >= 1.16]
            grid2 = np.linspace(1.16, np.max(length2), 80)
            intensity_interp2 = np.interp(grid2, length2, intensity2)
        # concatenate the two parts
        intensity_interp = np.concatenate((intensity_interp1, intensity_interp2))



        # smmoth the interpolated intensity
        # intensity_interp = scipy.ndimage.gaussian_filter1d(intensity_interp, 5)

        # store the interpolated intensity
        cilia_intensity[i] = intensity_interp
        cilia_length[i] = np.max(length)
        thresholds_local[i] = threshold

    
    return cilia_intensity, cilia_length, thresholds_local, n_BB_wo_cilia, n_BB_w_cilia




df = pd.read_excel("Filtered_Staining_List.xlsx", header=0) 

cq_dir = ".CQ/"

f_csv = open("path_list.csv", "w")
csv_writer = csv.writer(f_csv, lineterminator="\n")
csv_writer.writerow(["Plate id", "Position", "n_BB_wo_cilia", "n_BB_w_cilia"])


df2 = pd.read_csv(r'AllCQs.txt', sep='\t', header=0, encoding='ISO-8859-1')


intensity_profiles = {}
cilia_lengths = {}
annotations = {}
thresholds = {}

for init_file in os.listdir(cq_dir):
    if "CQl" not in init_file:
        continue
    plate_id = int(init_file.split('_')[0])
    well_id = init_file.split('_')[1]

    filename = init_file.replace("_CQl.txt", ".tif")

    # get the row with the same plate_id and well_id in df
    row = df[(df["Plate id"] == plate_id) & (df["Position"] == well_id)]

    if len(row) == 0:
        cell_line = "Unknown"
        temp = ["Unknown"]
    else:
        cell_line = row["Cell line"].values[0]
        temp = []
        for i in range(40, len(row.columns)):
            if row.iloc[0, i] == 1:
                # get the name of the column
                temp.append(row.columns[i])

    cilia_intensity, cilia_lengths_, thresholds_local, n_BB_wo_cilia, n_BB_w_cilia = intensity_profile(os.path.join(cq_dir, filename), df2)

    assert len(cilia_intensity) == len(cilia_lengths_), "cilia_intensity and cilia_length have different lengths"

    if str(plate_id) + "_" + well_id not in intensity_profiles:
        intensity_profiles[str(plate_id) + "_" + well_id] = []
        cilia_lengths[str(plate_id) + "_" + well_id] = []
        annotations[str(plate_id) + "_" + well_id] = []
        thresholds[str(plate_id) + "_" + well_id] = []

    for key in cilia_intensity:
        intensity_profiles[str(plate_id) + "_" + well_id].append(cilia_intensity[key])
        cilia_lengths[str(plate_id) + "_" + well_id].append(cilia_lengths_[key])
        annotations[str(plate_id) + "_" + well_id].append(temp)
        thresholds[str(plate_id) + "_" + well_id].append(thresholds_local[key])


    csv_writer.writerow([plate_id, well_id, n_BB_wo_cilia, n_BB_w_cilia])

f_csv.close()

pickle.dump(intensity_profiles, open("intensity_profiles_A_all.pkl", "wb"))
pickle.dump(cilia_lengths, open("cilia_lengths_all.pkl", "wb"))
pickle.dump(annotations, open("annotations_all.pkl", "wb"))
pickle.dump(thresholds, open("thresholds_all.pkl", "wb"))