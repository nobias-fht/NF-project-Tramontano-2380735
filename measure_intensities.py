import numpy as np
import matplotlib.pyplot as plt
import tifffile
import czifile
import os
import pandas as pd


# ------------------------- CHECK THIS PARAMETERS AND CHANGE THEM IF NEEDED -------------------------
# Channels order
cells_ch = 1
VP35_ch = 2
IRF3_ch = 3
dapi_ch = 4

# Folder paths
raw_images_folder_path = './input_Images'
cell_masks_folder_path = './cellpose/cells_cellpose_masks'
nuclei_masks_folder_path = './cellpose/nuclei_cellpose_masks'
output_folder_path = './output'
# ---------------------------------------------------------------------------------------------------


# functions

# Function to select a specific mask ID from a given mask
def select_maskID(mask, maskID):
    """
    Selects a specific mask ID from a given mask.

    Parameters:
    mask (numpy.ndarray): The input mask array.
    maskID (int): The ID of the mask to select.

    Returns:
    numpy.ndarray: A binary mask where the selected mask ID is 1 and all other values are 0.
    """
    mask_ids = np.atleast_1d(maskID)
    selected_mask = np.isin(mask, mask_ids)
    return selected_mask

# Function to get the unique values and their counts from an array, 
# keeping just the most connon and those that are within a certain 
# threshold of the most common count. 
# In this case, we are interested in the unique nuclei IDs within a cell mask 
# and we want to keep just one nucleus or multinuclear if they are of 
# comparable size (within a certain threshold set at 0.7 here).
def unique_counts(x, threshold_ratio=0.70):
    counts_dict = {}

    for v in x:
        if v == 0:
            continue
        counts_dict[v] = counts_dict.get(v, 0) + 1

    if not counts_dict:
        return [], []

    most_common_count = max(counts_dict.values())
    threshold = threshold_ratio * most_common_count

    filtered_items = [
        (v, c) for v, c in counts_dict.items()
        if c >= threshold
    ]

    filtered_items.sort(key=lambda item: item[1], reverse=True)

    values = [v for v, _ in filtered_items]
    counts = [c for _, c in filtered_items]

    return values, counts

# Create output folder if it doesn't exist
if not os.path.exists(output_folder_path):
    os.makedirs(output_folder_path) 

output_folder_masks_path = os.path.join(output_folder_path, "masks")
if not os.path.exists(output_folder_masks_path):    
    os.makedirs(output_folder_masks_path)



# Main script to measure intensities

# Initialize lists to store results
filenames = []
cell_IDs = []
x_positions = []
y_positions = []

cytoplasm_intensity_mean_IRF3 = []
cytoplasm_intensity_median_IRF3 = []
cytoplasm_intensity_25_IRF3 = []
cytoplasm_intensity_75_IRF3 = []
nuclei_intensity_mean_IRF3 = []
nuclei_intensity_median_IRF3 = []
nuclei_intensity_25_IRF3 = []
nuclei_intensity_75_IRF3 = []

cytoplasm_intensity_mean_VP35 = []
cytoplasm_intensity_median_VP35 = []
cytoplasm_intensity_25_VP35 = []
cytoplasm_intensity_75_VP35 = []
nuclei_intensity_mean_VP35 = []
nuclei_intensity_median_VP35 = []
nuclei_intensity_25_VP35 = []
nuclei_intensity_75_VP35 = []

# Load filenames
filenames_list = [ filename for filename in os.listdir(raw_images_folder_path) if filename.endswith('.czi') ]

for filename in filenames_list:
    print(f"Processing file: {filename}")
    # liad the image
    img = czifile.imread(os.path.join(raw_images_folder_path, filename))
    img = img.squeeze()

    # Load cell and nuclei masks
    cell_mask = tifffile.imread(os.path.join(cell_masks_folder_path, filename.replace('.czi', '.tif')))
    nuclei_mask = tifffile.imread(os.path.join(nuclei_masks_folder_path, filename.replace('.czi', '.tif')))

    # initialize a mask for all cytoplasm
    all_cyto_masks = np.zeros_like(cell_mask, dtype=bool)
    all_nuceli_masks = np.zeros_like(cell_mask, dtype=bool)

    # loop over all the cells IDs in the image
    for cell_id in np.unique(cell_mask):
        #print(f"Processing cell ID {cell_id} in file {filename}...")
        if cell_id == 0:
            continue

        # select current cell mask
        cell_mask_ID = select_maskID(cell_mask, cell_id)
        y_center, x_center = np.argwhere(cell_mask_ID > 0).mean(axis=0)

        # muclei masks within the cell
        nuclei_in_cell, _ = unique_counts(nuclei_mask[cell_mask_ID>0])
        nuclei_mask_ID = select_maskID(nuclei_mask, nuclei_in_cell)

        # cytoplasm mask is the cell mask minus the nuclei mask
        cytoplasm_mask = np.logical_and(cell_mask_ID, np.logical_not(nuclei_mask_ID))

        all_cyto_masks = all_cyto_masks+ cytoplasm_mask*cell_id
        all_nuceli_masks = all_nuceli_masks + nuclei_mask_ID*cell_id
        #all_cyto_masks = np.logical_or(all_cyto_masks, cytoplasm_mask)

        if np.sum(cytoplasm_mask) == 0 or np.sum(nuclei_mask_ID) == 0:
            continue

        # measure intensities in the cytoplasm and nuclei
        filenames.append(filename)
        cell_IDs.append(cell_id)
        x_positions.append(x_center)
        y_positions.append(y_center)

        cytoplasm_intensity_mean_IRF3.append(np.mean(img[IRF3_ch,cytoplasm_mask]))
        cytoplasm_intensity_median_IRF3.append(np.median(img[IRF3_ch,cytoplasm_mask]))
        cytoplasm_intensity_25_IRF3.append(np.quantile(img[IRF3_ch,cytoplasm_mask], 0.25))
        cytoplasm_intensity_75_IRF3.append(np.quantile(img[IRF3_ch,cytoplasm_mask], 0.75))
        nuclei_intensity_mean_IRF3.append(np.mean(img[IRF3_ch,nuclei_mask_ID]))
        nuclei_intensity_median_IRF3.append(np.median(img[IRF3_ch,nuclei_mask_ID]))
        nuclei_intensity_25_IRF3.append(np.quantile(img[IRF3_ch,nuclei_mask_ID], 0.25))
        nuclei_intensity_75_IRF3.append(np.quantile(img[IRF3_ch,nuclei_mask_ID], 0.75))

        cytoplasm_intensity_mean_VP35.append(np.mean(img[VP35_ch,cytoplasm_mask]))
        cytoplasm_intensity_median_VP35.append(np.median(img[VP35_ch,cytoplasm_mask]))
        cytoplasm_intensity_25_VP35.append(np.quantile(img[VP35_ch,cytoplasm_mask], 0.25))
        cytoplasm_intensity_75_VP35.append(np.quantile(img[VP35_ch,cytoplasm_mask], 0.75))
        nuclei_intensity_mean_VP35.append(np.mean(img[VP35_ch,nuclei_mask_ID]))
        nuclei_intensity_median_VP35.append(np.median(img[VP35_ch,nuclei_mask_ID]))
        nuclei_intensity_25_VP35.append(np.quantile(img[VP35_ch,nuclei_mask_ID], 0.25))
        nuclei_intensity_75_VP35.append(np.quantile(img[VP35_ch,nuclei_mask_ID], 0.75))


    data_table = pd.DataFrame({
        'filename': filenames,
        'cell_ID': cell_IDs,
        'x_position': x_positions,
        'y_position': y_positions,

        'cytoplasm_intensity_mean_IRF3': cytoplasm_intensity_mean_IRF3,
        'cytoplasm_intensity_median_IRF3': cytoplasm_intensity_median_IRF3,
        'cytoplasm_intensity_25_IRF3': cytoplasm_intensity_25_IRF3,
        'cytoplasm_intensity_75_IRF3': cytoplasm_intensity_75_IRF3,
        'nuclei_intensity_mean_IRF3': nuclei_intensity_mean_IRF3,
        'nuclei_intensity_median_IRF3': nuclei_intensity_median_IRF3,
        'nuclei_intensity_25_IRF3': nuclei_intensity_25_IRF3,
        'nuclei_intensity_75_IRF3': nuclei_intensity_75_IRF3,

        'cytoplasm_intensity_mean_VP35': cytoplasm_intensity_mean_VP35,
        'cytoplasm_intensity_median_VP35': cytoplasm_intensity_median_VP35,
        'cytoplasm_intensity_25_VP35': cytoplasm_intensity_25_VP35,
        'cytoplasm_intensity_75_VP35': cytoplasm_intensity_75_VP35,
        'nuclei_intensity_mean_VP35': nuclei_intensity_mean_VP35,
        'nuclei_intensity_median_VP35': nuclei_intensity_median_VP35,
        'nuclei_intensity_25_VP35': nuclei_intensity_25_VP35,
        'nuclei_intensity_75_VP35': nuclei_intensity_75_VP35
    })

    # save the data table as a CSV file
    data_table.to_csv(os.path.join(output_folder_path, f"data_table.csv"), index=False)
    
    # save the masks as TIFF files
    tifffile.imwrite(os.path.join(output_folder_masks_path, f"all_cyto_masks_{filename.replace('.czi', '.tif')}"), all_cyto_masks.astype(np.uint16))
    tifffile.imwrite(os.path.join(output_folder_masks_path, f"all_nuceli_masks_{filename.replace('.czi', '.tif')}"), all_nuceli_masks.astype(np.uint16))

print("Processing completed for all files.")