# Single-cell mapping of transcriptional and lineage landscape in developing mouse brain

![image](https://github.com/Y0NEKO/ChenLab-CREST/blob/main/figure1.jpg)

* CREST is a Cre-driven molecular recorder for single cell lineage tracing in vivo.  
* CREST reveals a spatiotemporal lineage landscape of developing ventral midbrain.  
* Single clonal analysis revealed 7 archetypes of ventral midbrain progenitors with
distinct fates.  
* Single clonal analysis reveals divergent differentiation routes and corresponding
* transcriptional dynamics of floor plate cells.  
* Semi-prospective clonal tracing links transcriptome states of progenitors to their
terminal fates across stages.  

# Preporcess raw single cell RNA data:
ex: sh 0.preprocess.sh /R1.fastq /R2.fastq outpath  
Produce CB/UMI and array sequence  

# CREST Barcode recovering
## Install and using our R package CRESTBuild:
```
devtools::install_github("Y0NEKO/ChenLab-CREST/CRESTBuild")
```
## Usage:
```
mycrest = LoadFiles("data/samplesheet.txt")
thread = 10
mycrest = FindScar(mycrest, thread, match = 1, mismatch = -3, gapOpening = 6, gapExtension = 1)
mycrest = INDELChangeForm(mycrest, thread)
mycrest = INDELCons(mycrest, thread)
write.csv(mycrest@finalIndel, file = "final_scarform.csv", quote = F, row.names = F)
```

Note: samplesheet.txt is a sample information file contains:reference file(fasta), cutsite, CB_UMI
(Please see our example data)


# License
The CRESTBuild package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See the GNU General Public License for more details.
A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3

## other analysis
other analysis method is in BarcodeRecover.R and BarcodeAnalysis.R
