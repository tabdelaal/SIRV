# SIRV
## Spatial inference of RNA velocity at the single-cell resolution

### Implementation description

Python implementation can be found in the 'SIRV' folder. The ```SIRV``` function takes as input 1) spatial transcriptomics data as ```scanpy/scvelo Anndata``` object, 2) scRNA-seq data as ```scanpy/scvelo Anndata``` object, having spliced and unspliced expressions and (optionally) cell label/metadata annotations, 3) number of principal vectors *(PVs)* needed for integration, 4) column names (identifiers) of the label/metadata annotations to be transferred from scRNA-seq to spatial data (optional).

For full description, please check the ```SIRV``` function description in ```main.py```.

### Results reproducibility

The ```DevelopingMouseBrain.py``` script shows how to reproduce the results obtained in the preprint when testing SRV on the Developing Mouse Brain Atlas data.

### Datasets

The processed spatial (HybISS) and scRNA-seq datasets of the Developing Mouse Brain can be downloaded as ```scanpy Anndata .h5``` files from [Surfdrive](https://surfdrive.surf.nl/files/index.php/s/11LsBOAmRhDEZOJ)

For citation and further information please refer to: "SIRV: Spatial inference of RNA velocity at the single-cell resolution", [bioRxiv2021](https://www.biorxiv.org/content/10.1101/2021.07.26.453774v)
