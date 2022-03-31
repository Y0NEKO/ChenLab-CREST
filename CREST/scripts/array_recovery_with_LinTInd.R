library(LinTInd)

data <- "../data/CB_UMI/CB_UMI"
fafile<-"../data/CB_UMI/V1.fasta"
cutsite<-"../data/CB_UMI/V1.cutSites"
data<-read.table(data,sep="\t",header=TRUE)
ref<-ReadFasta(fafile)
cutsite<-read.table(cutsite,col.names = c("indx","start","end"))
celltype<-read.table(celltype,header=TRUE,stringsAsFactors=FALSE)

scarinfo<-FindIndel(data=data,scarfull=ref,scar=cutsite,indel.coverage="All",type="test",cln=1)
scarinfo<-IndelForm(scarinfo,cln=1)
scarinfodf<-do.call(rbind.data.frame, your_list)
write.csv(scarinfodf,file = "final_scarform_v1.csv")
