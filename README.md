# Single-cell mapping of transcriptional and lineage landscape in developing mouse brain

![image]https://github.com/Y0NEKO/ChenLab-CREST/blob/main/figure1.jpg

CREST is a Cre-driven molecular recorder for single cell lineage tracing in vivo.
CREST reveals a spatiotemporal lineage landscape of developing ventral midbrain.
Single clonal analysis revealed 7 archetypes of ventral midbrain progenitors with
distinct fates.
Single clonal analysis reveals divergent differentiation routes and corresponding
transcriptional dynamics of floor plate cells.
Semi-prospective clonal tracing links transcriptome states of progenitors to their
terminal fates across stages.

# Aanlysis overview
step1--preporcess data produced from cellranger to get UMI and R2 information:
ex: sh 0.preprocess.sh /R1.fastq /R2.fastq outpath

outfiles:

"CB_UMI" reads information for next steps "cb.umi.tsv" UMI information

step2--call scar from "CB_UMI"

if you have only one scar,just use "PreData1" command

ex:

Rscript main.R PreData1 samplesheet.txt outpath 8

if you have two scars to split,use "PreData2" command

ex:

Rscript main.R PreData2 samplesheet.txt outpath 8

where,samplesheet.txt is a sample information file contains:reference file(fasta), cutsite, CB_UMI

step3--Scar formation transform

ex:

Rscript main.R CallScar samplesheet.txt outpath 8

step4--some basic analysis

ex:

Rscript main.R ScarAnalysis samplesheet.txt outpath 8

Bulk:

ex:

sh step_bulk.sh /home/xyz2020/scar/bulk/dose/HetV4M/1.rawdata/HetV4M2_FKDL202572705-1a/HetV4M2_FKDL202572705-1a_1.fq.gz /home/xyz2020/scar/bulk/dose/HetV4M/1.rawdata/HetV4M2_FKDL202572705-1a/HetV4M2_FKDL202572705-1a_2.fq.gz samplesheet3.txt bulktest2
