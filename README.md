# Spatial Converter
[Seurat by satijalab](https://satijalab.org/seurat/), [AnnData by scverse](https://anndata.readthedocs.io/en/stable/) and [BioConductor's SpatialExperiment](https://www.bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) are some of the most commonly used data objects for spatial omics data. These functions provide the ability to convert between Seurat and the other data objects, using csv files as an intermediary. 

# Dependencies
R:
- [Seurat](https://satijalab.org/seurat/)
- [readr](https://readr.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [stringr](https://stringr.tidyverse.org/)
- [SpatialExperiment](https://www.bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) (for conversion between Seurat and SpatialExperiment)
- [SPIAT](https://trigosteam.github.io/SPIAT/) (for conversion between Seurat and SPIAT-compatible SpatialExperiment)

Python: (if converting to and from AnnData)
- [pandas](https://pandas.pydata.org/)
- [anndata](https://anndata.readthedocs.io/en/stable/)
- [numpy](https://numpy.org/)

# Syntax
## AnnData to Seurat
### Python:
```
from anndata_functions import * 
Anndata_to_csv(adata, export_dir)
```
**Arguments:**
- *adata:* AnnData object in workspace
- *export_dir:* String. Directory for csv files to be written to

**Returns:**

None

### R:
```
source("seurat_functions.R")
csv_to_Seurat(import_dir, assay_name)
```
**Arguments:**
- *import_dir:* String. Directory where csv files are located
- *assay_name:* String. Desired name for Seurat assay.

**Returns:**

Seurat object

## Seurat to AnnData
### R:
```
source("seurat_functions.R")
Seurat_to_csv(seurat_obj, export_dir)
```
**Arguments:**
- *seurat_obj:* Seurat object, with spatial coordinates in the ``@image$centroids`` slot.
- *export_dir:* String. Directory for csv files to be written to

**Returns:**

None

### Python:
```
from anndata_functions import *
csv_to_Anndata(import_dir)
```
**Arguments:**
- *import_dir:* String. Directory where csv files are located

**Returns:**

AnnData Object

## Seurat to SpatialExperiment
### R:
```
source("seurat_functions.R")
Seurat_to_SPE(seurat_obj)
```
**Arguments:**
- *seurat_obj:* Seurat object, with spatial coordinates in the ``@image$centroids`` slot.

**Returns:**

SpatialExperiment Object

## Seurat to SPIAT-compatible SpatialExperiment
### R:
```
source("seurat_functions.R")
Seurat_to_SPIATSPE(seurat_obj)
```
**Arguments:**
- *seurat_obj:* Seurat object, with spatial coordinates in the ``@image$centroids`` slot.

**Returns:**

SpatialExperiment Object in format suitable for use with SPIAT functions

# Note:
- When converting from Seurat to AnnData:
  - The function assumes that the original Seurat object includes spatial data, in the image@centroids slot.
- When converting from AnnData to Seurat:
	- Feature/Gene metadata will be lost, due to Seurat having no consistent method for storing this information. If required, please add the information manually.
  - If there are HTML color codes in the metadata/obs dataframe, please remove before conversion, as R struggles to read in hexadecimal numerals.
- When converting from Seurat to SpatialExperiment, some metadata, such as dimension reduction will be lost.
- When converting from Seurat to SPIAT-compatible SpatialExperiment, all cell metadata, except Cell Identity will be lost.
