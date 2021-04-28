#source /home/csdn/miniconda3/bin/activate
#conda activate r4-base

fq1=$1
fq2=$2
sheet1=$3
sheet2=$4
outpath=$5
thread=$6
#sh 0.preprocess_bulk_single.sh ${fq1} ${fq2} ${outpath}
sh 0.preprocess_bulk.sh ${fq1} ${fq2} ${outpath} ${thread}
Rscript 1.main.R PreDataBulk ${sheet1} ${outpath}/v1 ${thread}
Rscript 1.main.R CallScar ${sheet1} ${outpath}/v1 ${thread}
Rscript 1.main.R ScarAnalysis ${sheet1} ${outpath}/v1 ${thread}

Rscript 1.main.R PreDataBulk ${sheet2} ${outpath}/v2 ${thread}
Rscript 1.main.R CallScar ${sheet2} ${outpath}/v2 ${thread}
Rscript 1.main.R ScarAnalysis ${sheet2} ${outpath}/v2 ${thread}
