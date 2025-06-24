# HPA Cilia Study - Code repository
This repository contains instructions and programming code for analyses of the data presented in the preprint ["Intrinsic Heterogeneity In Primary Cilia Revealed Through Spatial Proteomics (DOI: 10.1101/2024.10.20.619273)"](https://www.biorxiv.org/content/10.1101/2024.10.20.619273).

## (1) Downloading data from the protein atlas
1. Query the XML data for the ENSEMBL gene ID of your interest through a link composed of ```https://www.proteinatlas.org/``` + ENSEMBL ID (e.g., ```ENSG00000137691```) + ```.xml```, e.g., https://www.proteinatlas.org/ENSG00000137691.xml for CFAP300.
2. In the XML search for the subassay with subtype ("ciliated cell lines"): ```<subAssay type="human" subtype="ciliated cell lines"> ... </subAssay>```

![image](https://github.com/user-attachments/assets/d36c0976-48eb-4c22-b280-77ef9b9c0a4c)


3. Nested in the subAssay element you will find different "data" elements, which each represent the images for a specific cell line.

![image](https://github.com/user-attachments/assets/675443fd-7e28-4ee8-b95e-bf596999fde8)

4. This allows you to see the links for all z planes. E.g., for the image in the screenshot, this reveals 21 different z planes with their corresponding links:
   1. For z_index=1: ```https://images.proteinatlas.org/38585/2146_D7_42_blue_red_green.jpg```
   2. For z_index=2: ```https://images.proteinatlas.org/38585/2146_D7_41_blue_red_green.jpg```
   3. ...
   4. For z_index=20: ```https://images.proteinatlas.org/38585/2146_D7_24_blue_red_green.jpg```
   5. For z_index=21: ```https://images.proteinatlas.org/38585/2146_D7_23_blue_red_green.jpg```

5. Extract a download id for each z plane from the link:
   1. Remove the beginning (```https://images.proteinatlas.org/```) from the link. E.g., ```https://images.proteinatlas.org/38585/2146_D7_42_blue_red_green.jpg``` becomes ```/38585/2146_D7_42_blue_red_green.jpg```
   2. Remove the ending ( ```blue_red_green.jpg```) from the remaining link. E.g.,  ```/38585/2146_D7_42_blue_red_green.jpg``` becomes ```/38585/2146_D7_42_```
   3. The final image id (```/38585/2146_D7_42_```) remains

6. Recombine image id for each z plane to create links for downloading the individual tif images for each channel.
   1. For the blue DAPI / Nuclei channel, the download link will be: ```https://www.proteinatlas.org/download_file.php?filename=``` + the image id (e.g., ```/38585/2146_D7_42_```) + ```_blue&format=tif.gz```, so e.g., ```https://www.proteinatlas.org/download_file.php?filename=/38585/2146_D7_42_blue&format=tif.gz```
   2. For the red Cilia marker channel, the download link will be: ```https://www.proteinatlas.org/download_file.php?filename=``` + the image id (e.g., ```/38585/2146_D7_42_```) + ```_red&format=tif.gz```, so e.g., ```https://www.proteinatlas.org/download_file.php?filename=/38585/2146_D7_42_red&format=tif.gz```
   3. For the yellow Basal body marker channel, the download link will be: ```https://www.proteinatlas.org/download_file.php?filename=``` + the image id (e.g., ```/38585/2146_D7_42_```) + ```_yellow&format=tif.gz```, so e.g., ```https://www.proteinatlas.org/download_file.php?filename=/38585/2146_D7_42_yellow&format=tif.gz```
   4. For the green Protein of interest marker channel, the download link will be: ```https://www.proteinatlas.org/download_file.php?filename=``` + the image id (e.g., ```/38585/2146_D7_42_```) + ```_green&format=tif.gz```, so e.g., ```https://www.proteinatlas.org/download_file.php?filename=/38585/2146_D7_42_green&format=tif.gz```

7. Download for all z planes all four channel .tif.gz files and place them all together into one folder

8. Extract all .tif.gz files so they become .tif files and rename each file based on the z plane and the channel to match the following scheme: <image id after removing the front part between / and /> + ```_c``` + <channel number> + ```_z``` + <z_index> + ```.tif```, such as, e.g., ```2146_D7_42_c0_z0.tif``` for the DAPI channel image of the first plane for the example images above.
   1. Use channel number 0 for the DAPI channel
   2. Use channel number 1 for the Cilia channel
   3. Use channel number 3 for the Protein of interest channel
   4. Use channel number 4 for the Basal body channel

9. Assemble all images into a multi-channel multi-plane tif stack using the following script: XXX

## (2) Segmenting the images
Images downloaded and assembled as described above can be subjected to the [segmentation script](https://github.com/CellProfiling/HPA_Cilia_Study_Code/tree/main/Image%20segmentation) to create cilia, basal body, and nucleus segmentations.

## (3)  Run CiliaQ analysis on the segmented images
1. Run CiliaQ analysis (CiliaQ version V0.2.1) on all 7-channel images with the following CiliaQ settings

*TODO - add screenshots*

2. Collect all CiliaQ output files ending with ```_CQs.txt``` in a folder.
3. Add the table legend file provided [here](https://github.com/CellProfiling/HPA_Cilia_Study_Code/blob/main/Merging_CiliaQ_output_files/0000_Header_CQs.txt) to the folder with all CiliaQ's ```_CQs.txt``` output files. Make sure that this legend file is the first item when alphabetically sorting all files in the folder (if not rename to have it become the first item while keeping the file ending ```_CQs.txt``` intact)
4. Finally, assemble all ```_CQs.txt```-ending files produced by CiliaQ through concatenating all files ending with ```_CQs.txt``` into a single text file, e.g. by adding [this](https://github.com/CellProfiling/HPA_Cilia_Study_Code/blob/main/Merging_CiliaQ_output_files/Combiner.bat) Windows batch file to the folder and executing it; it will produce a new file called ```AllCQsFilesMerged.txt``` with all concatenated ```CQs.txt``` files, which can then be used for analyzing statistics on the cilia like the length or orientation angle (as shown in Figure 2 in the preprint).
 
## (4) Extract and normalize intensity profiles
1. Use the ```AllCQsFilesMerged.txt``` file from the previous step for this analysis.
2. Collect all ```_CQl.txt``` files created by CiliaQ in the previous steps in a folder.
3. Run the analysis as described in the [readme file](https://github.com/CellProfiling/HPA_Cilia_Study_Code/blob/main/Cilia%20profile%20analysis/readme_create_cilia_profiles.md) and use the scripts in [this repository](https://github.com/CellProfiling/HPA_Cilia_Study_Code/tree/main/Cilia%20profile%20analysis)
   1. Note that you need to have a specific excel file for this that lists all the images that you want to include and has specific columns available as explained in the readme file.
   2. The output files from this analysis can be further used to cluster and analyze intensity profiles.
4. To plot intensity profiles follow this [readme file](https://github.com/CellProfiling/HPA_Cilia_Study_Code/blob/main/Cilia%20profile%20analysis/readme_plot_profiles.md) and use the scripts in [this repository](https://github.com/CellProfiling/HPA_Cilia_Study_Code/tree/main/Cilia%20profile%20analysis)

## (5) Run cell cycle predictions on images

## (6) Analyze cell cycle prediction data

## (7) Other bioinformatic analysis of protein lists
- [Functional enrichment analysis of protein lists](https://github.com/CellProfiling/HPA_Cilia_Study_Code/tree/main/functional-enrichment-analysis)
- [Test if the distribution of the number of subcellular locations of ciliary proteins is significantly different from all proteins in the whole cell](https://github.com/CellProfiling/HPA_Cilia_Study_Code/blob/main/supplementary-figure-s2/statistical_analysis_protein_multilocalization.ipynb)
