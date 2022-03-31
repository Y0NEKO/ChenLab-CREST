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

# Aanlysis overview
## preporcess raw data:
ex: sh 0.preprocess.sh /R1.fastq /R2.fastq outpath  
Produce CB/UMI and array sequence  

## array recovering
ex: Rscript array_recovery_main.R PreData2 samplesheet.txt outpath 8  
where,samplesheet.txt is a sample information file contains:reference file(fasta), cutsite, CB_UMI  
then call consensus
ex: Rscript array_recovery_main.R CallScar samplesheet.txt outpath 8  

or use our package LinTInd to recover array, example of usage can check array_recovery_with_LinTInd.R.

## Lineage tree reconstruction
Analyze the steps and methods in script lineage_tree_reconstruction_method.R

## other analysis
other analysis method is in main_figure_analysis.R
