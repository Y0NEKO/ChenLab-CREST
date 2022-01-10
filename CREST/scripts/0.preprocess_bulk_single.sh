
#QC cut adapter
fq1=$1
fq2=$2
outpath=$3
cutadapt -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -e 0.1 -O 5 -m 50  -n  2 --pair-filter=both -o ${outpath}/R1.fq.gz -p ${outpath}/R2.fq.gz ${fq1} ${fq2}
fastqc ${outpath}/R1.fq.gz
fastqc ${outpath}/R2.fq.gz


#assemble
flash -M 250 -r 250 -f 250 -s 25 -O -d ${outpath} 2>&1 ${fq1} ${fq2}

#select reads and UMIs (might be array)
sed -n '2~4p' ${outpath}/out.extendedFrags.fastq > ${outpath}/reads.txt
#awk -F "GAACGC|ATGGCC" '{print$2}' ${dirpath}/reads.txt > ${dirpath}/umi.txt
awk -F "ATGGCCATCAT" '{print substr($1,length($1)-15)}' ${outpath}/reads.txt > ${outpath}/umi.txt
paste ${outpath}/reads.txt ${outpath}/umi.txt > ${outpath}/reads_UMI_pre.txt
awk 'length($2)<=18 && length($2)>=6' ${outpath}/reads_UMI_pre.txt > ${outpath}/CB_UMI




