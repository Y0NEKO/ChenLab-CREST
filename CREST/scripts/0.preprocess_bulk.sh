
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
mkdir ${outpath}/v1
mkdir ${outpath}/v2
sed -n '2~4p' ${outpath}/out.extendedFrags.fastq | grep "AGTTCT" | grep "ATTTAA"  > ${outpath}/v1/reads_v1.txt
#awk -F "GAACGC|ATGGCC" '{print$2}' ${dirpath}/reads.txt > ${dirpath}/umi.txt
sed -n '2~4p' ${outpath}/out.extendedFrags.fastq | grep "AGTTCT" | grep "TTATAA"  > ${outpath}/v2/reads_v2.txt
awk -F "ATGGCCATCAT" '{print substr($1,length($1)-15)}' ${outpath}/v1/reads_v1.txt > ${outpath}/v1/umi_v1.txt
awk -F "ATGGCCATCAT" '{print substr($1,length($1)-15)}' ${outpath}/v2/reads_v2.txt > ${outpath}/v2/umi_v2.txt
paste ${outpath}/v1/reads_v1.txt ${outpath}/v1/umi_v1.txt > ${outpath}/v1/reads_UMI_pre_v1.txt
paste ${outpath}/v2/reads_v2.txt ${outpath}/v2/umi_v2.txt > ${outpath}/v2/reads_UMI_pre_v2.txt
awk 'length($2)<=18 && length($2)>=6' ${outpath}/v1/reads_UMI_pre_v1.txt > ${outpath}/v1/CB_UMI
awk 'length($2)<=18 && length($2)>=6' ${outpath}/v2/reads_UMI_pre_v2.txt > ${outpath}/v2/CB_UMI




