source /home/csdn/miniconda3/bin/activate
conda activate r4-base

input=$1
sheet=$2
outpath=$3
thread=$4
sh 0.preprocess_sc.sh ${input} ${outpath}

#1.1 for two arrays 
#find array and separate v1 / v2
Rscript 1.main.R PreData2 ${sheet} ${outpath} ${thread}
Rscript 1.main.R CallScar ${sheet} ${outpath}/v1 ${thread}
Rscript 1.main.R CallScar ${sheet} ${outpath}/v2 ${thread}
Rscript 1.main.R ScarAnalysis ${sheet} ${outpath}/v1 ${thread}
Rscript 1.main.R ScarAnalysis ${sheet} ${outpath}/v2 ${thread}

#1.2 for one array
Rscript 1.main.R PreData1 ${sheet} ${outpath} ${thread}
Rscript 1.main.R CallScar ${sheet} ${outpath} ${thread}
Rscript 1.main.R ScarAnalysis ${sheet} ${outpath} ${thread}


