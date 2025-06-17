import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy


# Change this with your credentials
def intensity_profile(file_path, df2):
    # file_path = r"C:\Users\sunh\Downloads\2131_G11_88_new_stack_CQl.txt"

    # Read the file into a DataFrame, skipping unnecessary metadata lines
    df = pd.read_csv(file_path, sep='\t', header=0, encoding='ISO-8859-1')  # Assuming the metadata is on the first line
    n_BB_wo_cilia = 0
    n_BB_w_cilia = 0
    cilia_intensity = {}
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

        cilia_intensity[i] = intensity_interp
        thresholds_local[i] = threshold
    return cilia_intensity, thresholds_local

all_profiles = []
all_thresholds = []
df2 = pd.read_csv(r'GPR75_002_AllCQs.txt', sep='\t', header=0, encoding='ISO-8859-1')


for init_file in os.listdir("GPR75_002"):
    filename = init_file.replace("_CQl.txt", ".tif")
    cilia_intensity, thresholds = threshold = intensity_profile(os.path.join("GPR75_002", init_file), df2)
    for key in cilia_intensity:
        all_profiles.append(cilia_intensity[key])
        all_thresholds.append(thresholds[key])

# normalize the profiles
all_profiles = np.array(all_profiles)
for i in range(len(all_profiles)):
    all_profiles[i] = all_profiles[i] - all_thresholds[i]
all_profiles[all_profiles < 0] = 0

max_ = max(10000, np.max(all_profiles))
profiles = all_profiles / max_

profiles = scipy.ndimage.gaussian_filter1d(profiles, axis=1, sigma=5)
# plot it using seaborn heatmap
plt.figure(figsize=(4, 10))
sns.heatmap(profiles, cmap="magma", cbar=True)

# no ticks on both axes
plt.xticks([])
plt.yticks([])

# save the figure as svg
plt.savefig("GPR75_002.svg", dpi=300, bbox_inches='tight', pad_inches=0.1)