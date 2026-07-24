# Installation instructions

If you don't have python please follow installation instructions here: https://www.python.org/downloads/

If you don't have uv please follow installation instructions here: https://docs.astral.sh/uv/getting-started/installation/

1. Clone the repository:

```bash
git clone https://github.com/nobias-fht/Data-Reducer-2009821-DHA-PILOT.git
cd NF-project-Tramontano-2380735
```

2. Create a virtual environment with `uv`:

```bash
uv sync
```


# User Guide

## Segmentation of cells and nuclei with cellpose
Open the file **"segment_cells_and_nuclei.py"** and edit the parameters section if needed. 
Then on the terminal window type:

```bash
uv run segment_cells_and_nuclei.py
```
Two folders called 'cells_cellpose_masks' and 'nuclei_cellpose_masks' will be created containing the segmentation masks for cells and nuclei of each raw image. Depending on your machine each image may take approximatelly 5-10 minutes to be processed.

---

## Quantify fluorescence intensities

Open the file **"measure_intensities.py"** and edit the parameters section if needed. 
Then on the terminal window type:

```bash
uv run measure_intensities.py
```
An output folder will be creted containing: 
- A **".csv file"** with single cell fluorescence intensity maeasurements on the VP35 and IRF-3 channels within the nucleus and cytoplasm. 
- A new **"masks folder"** containing the masks for nuclei and cytoplasm with the same cellID/label found in the .csv file.