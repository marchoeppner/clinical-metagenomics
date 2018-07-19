![](images/ikmb_bfx_logo.png)
# IKMB Metagenomic profiling pipeline

## Overview

This pipelines analyses short reads and identifies the most likely species in the respective sample. 

## Executing the pipeline

This pipeline requires Nextflow >= 0.30.2. All dependencies are provided through Bioconda. Please make sure that conda/miniconda2 are available before starting the pipeline. 

Please note that Pathoscope, Metaphlan and Kaiju require reference databases which are *not* included with the Bioconda packages. On RZCluster, these are
available automatically through the included config file. 

This pipeline uses a CSV formatted file as input, using the following columns (separated by ';'):

* patientID
* sampleID
* sample type (e.g. Blood)
* read type (e.g. 2x150bp)
* platform (e.g. HiSeq 4000)

A simple script to create such a basic input format from a folder of FastQ files is included in the bin/ subfolder.

To run the pipeline, do:

`nextflow -c nextflow.config run main.nf --samples Samples.csv`

where nextflow.config and main.nf need to be passed with their full path to the local git clone. 




