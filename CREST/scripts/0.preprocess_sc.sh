#!/bin/bash
R1=$1
R2=$2
outpath=$3
#sh /picb/sysgenomics2/projects/wangluyue/Xie/Dual_2nd/cellranger_count.sh ${outpath} fq

umi_tools whitelist --stdin ${R1} --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --set-cell-number=10000 --log2stderr > ${outpath}/whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ${R1} \
                  --stdout ${outpath}/R1_extracted.fastq.gz \
                  --read2-in ${R2} \
                  --read2-out=${outpath}/R2_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist=${outpath}/whitelist.txt


mkdir ${outpath}/UMI_CB
zcat ${outpath}/R2_extracted.fastq.gz | sed -n '1~4p' | awk -F "_" '{print $1,$2,$3}' | awk '{print$1,$2,$3}' > ${outpath}/UMI_CB/ID_CB_UMI
zcat ${outpath}/R2_extracted.fastq.gz | sed -n '2~4p' > ${outpath}/UMI_CB/Reads
paste ${outpath}/UMI_CB/ID_CB_UMI ${outpath}/UMI_CB/Reads | awk -v OFS="\t" '{print$1,$4,$2,$3}' > ${outpath}/UMI_CB/cb.umi.tsv
awk '$3!="" && $4!=""' ${outpath}/UMI_CB/cb.umi.tsv > ${outpath}/UMI_CB/CB_UMI
sed -i '1i\Read.ID\tRead.Seq\tCell.BC\tUMI' ${outpath}/UMI_CB/CB_UMI
echo "DONE"

#samtools view ${input} > ${outpath}/possorted_genome_bam.sam
#mkdir ${outpath}/UMI_CB
#gawk -f 0.reads_umi.sh ${outpath}/possorted_genome_bam.sam > ${outpath}/UMI_CB/cb.umi.tsv
#awk '$3!="" && $4!=""' ${outpath}/UMI_CB/cb.umi.tsv > ${outpath}/UMI_CB/CB_UMI
#rm ${outpath}/possorted_genome_bam.sam


