# BIDCell to Seurat

This repository contains an enabler script to convert the BIDcell segmentation output to a Seurat object for the processing of 10X Genomics _Xenium_ spatial transcriptomics data.

You can find detailed information on BIDcell at the following citation
> Fu, X., Lin, Y., Lin, D., Mechtersheimer, D., Wang, C., Ameen, F., Ghazanfar, S., Patrick, E., Kim, J., & Yang, J. Y. H. BIDCell: Biologically-informed self-supervised learning for segmentation of subcellular spatial transcriptomics data. Nat Commun 15, 509 (2024). https://doi.org/10.1038/s41467-023-44560-w

And at the following github repository: https://github.com/SydneyBioX/BIDCell

For Seurat, multiple citations and help resources are available via their dedicated website: https://satijalab.org/seurat/

# Getting started

## Environment

A conda environment with both BIDcell and Seurat installed is recommended.

## Running the script

> [!IMPORTANT]  
> run BIDcell beforehand!

1. Source the BIDcell_to_Seurat function in the BIDcell_to_Seurat.R script to your own script/session like so:

```R
source("path/to/BIDcell_to_Seurat.R")
```

2. Run the conversion function with the following parameters
    * name = output name of file and field of view (FOV)
    * BIDCell_outs_folder = location of the BIDCell output
    * output_folder = location where Seurat rds is supposed to be saved

# Contact

Made by [Markus Boesch](markus.boesch@kuleuven.be) and [Francisco Pestana](francisco.pestana@kuleuven.be)

# Citation
10.5281/zenodo.13619779

