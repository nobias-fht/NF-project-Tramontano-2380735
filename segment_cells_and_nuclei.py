import numpy as np
import matplotlib.pyplot as plt
from cellpose import models, io, utils
import tifffile
import czifile
import os

# ------------------------- CHECK THIS PARAMETERS AND CHANGE THEM IF NEEDED -------------------------
# Channels
cells_ch = 1
dapi_ch = 4

# Folder paths
raw_images_folder_path = './input_Images'
cell_model_path = './cellpose/CP_20260710_154647'
cell_output_masks_folder_path = './cellpose/cells_cellpose_masks'
nuclei_output_masks_folder_path = './cellpose/nuclei_cellpose_masks'
# ---------------------------------------------------------------------------------------------------

# create output folders if they don't exist 
if not os.path.exists(cell_output_masks_folder_path):
    os.makedirs(cell_output_masks_folder_path)
if not os.path.exists(nuclei_output_masks_folder_path):
    os.makedirs(nuclei_output_masks_folder_path)

model_nuclei = models.CellposeModel(model_type='cyto3', gpu=False)
model_cells = models.CellposeModel(model_type=cell_model_path, gpu=False)

for filename in os.listdir(raw_images_folder_path):
    if filename.endswith('.czi'):
        # Load the image
        img = czifile.imread(os.path.join(raw_images_folder_path, filename))
        img = np.squeeze(img, axis=None)
        print(f"Processing {filename}")
  
        # Run Cellpose to segment nuceli with pretrained model
        nuclei_im = img[dapi_ch-1]
        nuclei_masks, flows, styles = model_nuclei.eval(nuclei_im, diameter=200, flow_threshold=0.4, cellprob_threshold=0.0)
        print(f"- Segmented {np.max(nuclei_masks)} nuclei")

        # Run Cellpose to segment cells with fine-tuned model
        cells_im = img[[cells_ch-1,dapi_ch-1],:,:]
        cells_masks, flows, styles = model_cells.eval(cells_im, channels=[0,1])
        print(f"- Segmented {np.max(cells_masks)} cells")
        
        # Save the masks as a TIFF file
        output_filename = os.path.splitext(filename)[0] + '.tif'
        tifffile.imwrite(os.path.join(nuclei_output_masks_folder_path, output_filename), nuclei_masks.astype(np.uint16))
        tifffile.imwrite(os.path.join(cell_output_masks_folder_path, output_filename), cells_masks.astype(np.uint16))

print("Segmentation completed for all images.")