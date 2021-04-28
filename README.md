

**Aim preprocess scar single cell data and basic analysis**

Author W&L

create data 2021.04.28**

USAGE

**step1--preporcess data produced from cellranger to get UMI and R2 information:**

ex:
sh 0.preprocess.sh /R1.fastq /R2.fastq outpath

outfiles:

"CB_UMI" reads information for next steps
"cb.umi.tsv" UMI information

**step2--call scar from "CB_UMI"**

if you have only one scar,just use "PreData1" command

ex:

Rscript main.R PreData1 samplesheet.txt outpath 8


if you have two scars to split,use "PreData2" command

ex:

Rscript main.R PreData2 samplesheet.txt outpath 8

where,samplesheet.txt is a sample information file contains:reference file(fasta), cutsite, CB_UMI


**step3--Scar formation transform**

ex:

Rscript main.R CallScar samplesheet.txt outpath 8


**step4--some basic analysis**

ex:

Rscript main.R ScarAnalysis samplesheet.txt outpath 8



**Bulk:**

ex:

sh step_bulk.sh /home/xyz2020/scar/bulk/dose/HetV4M/1.rawdata/HetV4M2_FKDL202572705-1a/HetV4M2_FKDL202572705-1a_1.fq.gz /home/xyz2020/scar/bulk/dose/HetV4M/1.rawdata/HetV4M2_FKDL202572705-1a/HetV4M2_FKDL202572705-1a_2.fq.gz samplesheet3.txt bulktest2

