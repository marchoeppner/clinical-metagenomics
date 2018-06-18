![](images/ikmb_bfx_logo.png)
# IKMB Metagenomic profiling pipeline

## Overview

This pipelines analyses short reads and identifies the most likely species in the respective sample. 

## Executing the pipeline

All dependencies for this pipeline are provided through Bioconda. To enable to pipeline environment, please do:
`conda env create -f environment.yml
source activate metagenomics
`

Please note that Pathoscope, Metaphlan and Kaiju require reference databases which are *not* included with the Bioconda packages. On RZCluster, these are
available automatically through the included config file. 

To run the pipeline, do:

`nextflow -c nextflow.config run main.nf --folder /path/to/reads`

where nextflow.config and main.nf need to be passed with their full path to the local git clone. 

If you wish to run the pipeline directly from the git repository on our server:

`nextflow run bfx-core/NF-metagenomic-profiling -hub ikmb --folder /path/to/folder`




