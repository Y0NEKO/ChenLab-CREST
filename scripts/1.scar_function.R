#Create on 8 April, 2021
#stat and plot
library(ggplot2)
library(Biostrings)
library(stringdist)
library(rlist)
library(dplyr) 
library(circlize)
library(RColorBrewer)
library(parallel)
library(Seurat)
library(dplyr)



#input data
ReadSheet = function(filename){
  filels = read.table(filename,sep = "\t",header = T)
  files = filels$Path
  names(files) = filels$File
  return(files)
}
ReadCutsite = function(filename){
  segref = read.table(filename)
  colnames(segref) = c("segid","start","end")
  
  scar = NULL
  type = NULL
  for (i in 1:nrow(segref)) {
    scar = c(scar,segref[i,]$start:segref[i,]$end)
    type = c(type,rep(segref[i,]$segid,(segref[i,]$end-segref[i,]$start)+1))
  }
  scarseg = data.frame("scar" = scar,"type" = as.character(type))
  
  return(scarseg)
}
ReadFasta = function(filename){
  sv = read.table(filename)
  scarfull = DNAString(sv[2,1])
  scar = DNAString(sv[4,1])
  return(c("scarfull" = scarfull,"scar" = scar))
}

#find_barcode
FindScar = function(data,scarfull,scar,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  point = length(scarfull)-length(scar)-length(scar)-6
  find_barcode<-function(data){
    #s3<-reverseComplement(DNAString(as.character(data[3])))
    s3<-DNAString(as.character(data))
    #alig3<-pairwiseAlignment(sv3,s3,substitutionMatrix = mat,type="global",gapOpening = 6, gapExtension= 1)
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)

    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos<-matchPattern(scar,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read<-as.character(subject(alig))
      if(length(scar_pos)!=0){
        r_scar<-testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
        stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
        p <- as.character(alignedPattern(alig)[[1L]])
        s <- as.character(alignedSubject(alig)[[1L]])
        p <- strsplit(p,"")[[1]]
        s <- strsplit(s,"")[[1]]
        lenp <- length(p)
        index <- 1
        i <- 0
        del_flag <- F
        in_flag <- F
        del_start<-c()
        del_end<-c()
        ins_start<-c()
        ins_width<-c()
        ins_editseq<-c()
        while(i < lenp){
          i <- i + 1
          if(s[[i]] == '-'){
            if(!del_flag){
              del_flag <- T
              del_start<-c(del_start,index)
              width <- 1
            }
            else{
              width <- width + 1
            }
          }
          else{
            if(del_flag){
              del_flag <- F
              del_end<-c(del_end,index)
              #print(paste("del stop width", width))
            }
          }
          if(p[[i]] == '-'){
            if(!in_flag){
              in_flag <- T
              ins_start<-c(ins_start,index)
              width <- 1
              
              editseq<-s[[i]]
            }
            else{
              width <- width + 1
              
              editseq<-paste0(editseq,s[[i]])
            }
          }
          else{
            if(in_flag){
              in_flag <- F
              ins_width<-c(ins_width,width)
              ins_editseq<-c(ins_editseq,editseq)
              #print(paste("in stop width", width))
            }
          }
          if(p[[i]] != '-')
            index <- index + 1
        }
        if(del_flag){
          del_flag <- F
          del_end<-c(del_end,index)
        }
        if(in_flag){
          in_flag <- F
          ins_width<-c(ins_width,width)
          ins_editseq<-c(ins_editseq,editseq)
          #print(paste("in stop width", width))
        }
        ins<-IRanges(start = ins_start,width = ins_width)
        mcols(ins)$seq <- ins_editseq
        del<-IRanges(start = del_start,end = del_end)
        
      }else{
        r_scar<-"unknown"
        del<-NA
        ins<-NA
      }
    }
    type<-"-"
    fin_dat<-data.frame(new.reads=r_read,scar.BC=r_scar,type=type)
    return(list(del,ins,fin_dat))
  }
  
  
  cl = makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  environment(point) <- .GlobalEnv
  environment(scarfull) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  environment(data) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(find_barcode) <- .GlobalEnv
  environment(testFunction) <- .GlobalEnv
  clusterExport(cl,c('point','scarfull','scar','data','mat','find_barcode','testFunction'),envir = environment())
  scar_BC = parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)
  
  #output
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data[,-5]
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1))
}

FindScarBulk = function(data,scarfull,scar,cln){
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction <- function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  point<-length(scarfull)-length(scar)-length(scar)-6
  find_barcode<-function(data){
    #s3<-reverseComplement(DNAString(as.character(data[3])))
    s3<-DNAString(as.character(data))
    #alig3<-pairwiseAlignment(sv3,s3,substitutionMatrix = mat,type="global",gapOpening = 6, gapExtension= 1)
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos<-matchPattern(scar,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read<-as.character(subject(alig))
      if(length(scar_pos)!=0){
        r_scar<-testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
        stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
        p <- as.character(alignedPattern(alig)[[1L]])
        s <- as.character(alignedSubject(alig)[[1L]])
        p <- strsplit(p,"")[[1]]
        s <- strsplit(s,"")[[1]]
        lenp <- length(p)
        index <- 1
        i <- 0
        del_flag <- F
        in_flag <- F
        del_start<-c()
        del_end<-c()
        ins_start<-c()
        ins_width<-c()
        ins_editseq<-c()
        while(i < lenp){
          i <- i + 1
          if(s[[i]] == '-'){
            if(!del_flag){
              del_flag <- T
              del_start<-c(del_start,index)
              width <- 1
            }
            else{
              width <- width + 1
            }
          }
          else{
            if(del_flag){
              del_flag <- F
              del_end<-c(del_end,index)
              #print(paste("del stop width", width))
            }
          }
          if(p[[i]] == '-'){
            if(!in_flag){
              in_flag <- T
              ins_start<-c(ins_start,index)
              width <- 1
              
              editseq<-s[[i]]
            }
            else{
              width <- width + 1
              
              editseq<-paste0(editseq,s[[i]])
            }
          }
          else{
            if(in_flag){
              in_flag <- F
              ins_width<-c(ins_width,width)
              ins_editseq<-c(ins_editseq,editseq)
              #print(paste("in stop width", width))
            }
          }
          if(p[[i]] != '-')
            index <- index + 1
        }
        if(del_flag){
          del_flag <- F
          del_end<-c(del_end,index)
        }
        if(in_flag){
          in_flag <- F
          ins_width<-c(ins_width,width)
          ins_editseq<-c(ins_editseq,editseq)
          #print(paste("in stop width", width))
        }
        ins<-IRanges(start = ins_start,width = ins_width)
        mcols(ins)$seq <- ins_editseq
        del<-IRanges(start = del_start,end = del_end)
        
      }else{
        r_scar<-"unknown"
        del<-NA
        ins<-NA
      }
    }
    type<-"-"
    fin_dat<-data.frame(new.reads=r_read,scar.BC=r_scar,type=type)
    return(list(del,ins,fin_dat))
  }
  library(parallel)
  cl<-makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  environment(point) <- .GlobalEnv
  environment(scarfull) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  environment(data) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(find_barcode) <- .GlobalEnv
  environment(testFunction) <- .GlobalEnv
  clusterExport(cl,c('point','scarfull','scar','data','mat','find_barcode','testFunction'),envir = environment())
  scar_BC<-parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)
  
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data[,-5]
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1))
}

SplitScar = function(data,scarfull1,scar1,scarfull2,scar2,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  point1 = length(scarfull1)-length(scar1)-length(scar1)-6
  point2 = length(scarfull2)-length(scar2)-length(scar2)-6
  find_barcode<-function(data){
    #s3<-reverseComplement(DNAString(as.character(data[3])))
    s3<-DNAString(as.character(data))
    #alig3<-pairwiseAlignment(sv3,s3,substitutionMatrix = mat,type="global",gapOpening = 6, gapExtension= 1)
    alig1<-pairwiseAlignment(scarfull1,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    alig2<-pairwiseAlignment(scarfull2,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    if(score(alig1)>score(alig2)){
      alig<-alig1
      scar<-scar1
      type<-"V1"
      point<-point1
    }else if(score(alig1)<score(alig2)){
      alig<-alig2
      scar<-scar2
      type<-"V2"
      point<-point2
    }else{
      alig<-alig2
      scar<-scar2
      r_read<-"unknown"
      r_scar<-"unknown"
      type<-"unknown"
      point<-point2
      del=NA
      ins=NA
    }
    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos<-matchPattern(scar,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read<-as.character(subject(alig))
      if(length(scar_pos)!=0){
        r_scar<-testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
        stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
        p <- as.character(alignedPattern(alig)[[1L]])
        s <- as.character(alignedSubject(alig)[[1L]])
        p <- strsplit(p,"")[[1]]
        s <- strsplit(s,"")[[1]]
        lenp <- length(p)
        index <- 1
        i <- 0
        del_flag <- F
        in_flag <- F
        del_start<-c()
        del_end<-c()
        ins_start<-c()
        ins_width<-c()
        ins_editseq<-c()
        while(i < lenp){
          i <- i + 1
          if(s[[i]] == '-'){
            if(!del_flag){
              del_flag <- T
              del_start<-c(del_start,index)
              width <- 1
            }
            else{
              width <- width + 1
            }
          }
          else{
            if(del_flag){
              del_flag <- F
              del_end<-c(del_end,index)
              #print(paste("del stop width", width))
            }
          }
          if(p[[i]] == '-'){
            if(!in_flag){
              in_flag <- T
              ins_start<-c(ins_start,index)
              width <- 1
              
              editseq<-s[[i]]
            }
            else{
              width <- width + 1
              
              editseq<-paste0(editseq,s[[i]])
            }
          }
          else{
            if(in_flag){
              in_flag <- F
              ins_width<-c(ins_width,width)
              ins_editseq<-c(ins_editseq,editseq)
              #print(paste("in stop width", width))
            }
          }
          if(p[[i]] != '-')
            index <- index + 1
        }
        if(del_flag){
          del_flag <- F
          del_end<-c(del_end,index)
        }
        if(in_flag){
          in_flag <- F
          ins_width<-c(ins_width,width)
          ins_editseq<-c(ins_editseq,editseq)
          #print(paste("in stop width", width))
        }
        ins<-IRanges(start = ins_start,width = ins_width)
        mcols(ins)$seq <- ins_editseq
        del<-IRanges(start = del_start,end = del_end)
        
      }else{
        r_scar<-"unknown"
        del<-NA
        ins<-NA
      }
    }
    fin_dat<-data.frame(new.reads=r_read,scar.BC=r_scar,type=type)
    return(list(del,ins,fin_dat))
  }
  
  
  cl = makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  environment(point1) <- .GlobalEnv
  environment(point2) <- .GlobalEnv
  environment(scarfull1) <- .GlobalEnv
  environment(scarfull2) <- .GlobalEnv
  environment(scar1) <- .GlobalEnv
  environment(scar2) <- .GlobalEnv
  environment(data) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(find_barcode) <- .GlobalEnv
  environment(testFunction) <- .GlobalEnv
  clusterExport(cl,c('point1','point2','scarfull1','scar1','scarfull2','scar2','data','mat','find_barcode','testFunction'),envir = environment())
  scar_BC = parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)

  #output
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data[data$type=="V1",-5]
  data_v2<-data[data$type=="V2",-5]
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"v1/UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data_v2<-data_v2[data_v2$scar.BC!="unknown",]
  write.table(data_v2,"v2/UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  scar_BC_v1<-scar_BC[data$type=="V1"]
  scar_BC_v2<-scar_BC[data$type=="V2"]
  saveRDS(scar_BC_v1,"v1/indel.rds")
  saveRDS(scar_BC_v2,"v2/indel.rds")
  return(list("V1_INDEL" = scar_BC_v1,"V2_INDEL" = scar_BC_v2,"V1_Scar" = data_v1,"V2_Scar" = data_v2))
}

#data preprocess
INDELChangeForm = function(INDEL_ranges,data,scarref,outpath,cln){
  change_form_stat<-function(indel){
    indel<-indel[c(1,2)]
    indel<-unlist(indel)
    if(length(indel)==0){
      return("unknown")
    }else{
      ins<-data.frame(indel[2])
      ins$seq=as.character(indel[[2]]@elementMetadata$seq)
      site_ins<-apply(ins,1,function(x){c(x[1]:x[2])})
      if(dim(ins)[1]==1){
        site_ins<-list(as.numeric(site_ins))
      }
      cutsite_ins<-lapply(site_ins,function(x){unique(scarref$type[scarref$scar %in% x])})
      tag_ins<-apply(ins,1,function(x){paste0(x[3],"I+",x[1],x[4])})
      tag_ins<-lapply(tag_ins,function(x){rep(x,length(cutsite_ins[[which(tag_ins==x)]]))})
      tag_ins<-unlist(tag_ins)
      cutsite_ins<-unlist(cutsite_ins)
      del<-data.frame(indel[1])
      site_del<-apply(del,1,function(x){c(x[1]:x[2])})
      if(dim(del)[1]==1){
        site_del<-list(as.numeric( site_del))
      }
      cutsite_del<-lapply(site_del,function(x){unique(scarref$type[scarref$scar %in% x])})
      tag_del<-apply(del,1,function(x){paste0((x[3]-1),"D+",x[1])})
      tag_del<-lapply(tag_del,function(x){rep(x,length(cutsite_del[[which(tag_del==x)]]))})
      tag_del<-unlist(tag_del)
      cutsite_del<-unlist(cutsite_del)
      tag<-c(tag_del,tag_ins)
      cutsite<-c(cutsite_del,cutsite_ins)
      tag_all<-rep("NONE",11)
      if(length(tag)==0){
        return(paste(tag_all,collapse = "_"))
      }else{
        for(x in c(1:length(tag))){
          if(tag_all[as.numeric(cutsite[x])]=="NONE"){
            tag_all[as.numeric(cutsite[x])]<-tag[x]
          }else{
            tag_all[as.numeric(cutsite[x])]<-paste(tag_all[as.numeric(cutsite[x])],tag[x],sep="&")
          }
        }
      }
      return(paste(tag_all,collapse = "_"))
    }
  }
  cl<-makeCluster(cln)
  environment(change_form_stat) <- .GlobalEnv
  environment(INDEL_ranges) <- .GlobalEnv
  environment(scarref) <- .GlobalEnv
  clusterExport(cl,c('INDEL_ranges','change_form_stat',"scarref"), envir = environment())
  scar_form_p<-parLapply(cl,INDEL_ranges,change_form_stat)
  stopCluster(cl)
  scar_form<-unlist(scar_form_p)
  scar_form<-gsub(" ","",scar_form)
  data$scar_form_all<-scar_form
  data$scar_form<-unlist(lapply(scar_form,function(x){paste(unlist(strsplit(x,"_"))[2:10],collapse = "_")}))
  write.csv(data,paste0(outpath,"/indel_pattern.csv"),quote=F,row.names = F)
  return(data)
}
INDELCons = function(INDEL_ranges,data,scarref,outpath,cln){
  Cell.BC<-data.frame(table(data$Cell.BC))
  Cell.BC<-Cell.BC[Cell.BC$Freq>1,]
  INDEL_ranges<-INDEL_ranges[data$Cell.BC %in% Cell.BC$Var1]
  data<-data[data$Cell.BC %in% Cell.BC$Var1,]
  data_1<-data[,c(3,4,7,9)]
  max_reads_stat<-function(x,dat){
    temreads<-dat[dat$Cell.BC==x,]
    read_data<-unique(temreads[,c(2,4)])
    read_data<-data.frame(table(read_data$scar_form))
    read_data<-read_data[order(-read_data$Freq),]
    if(read_data$Freq[1]==1){
      fin_read="unknown"
      UMI_pro=0
      UMI_num=0
    }else{
      fin_read=read_data$Var1[1]
      UMI_num=read_data$Freq[1]
      UMI_pro=UMI_num/sum(read_data$Freq)
    }
    fin_dat<-data.frame(Cell.BC=x,UMI_num=UMI_num,scar_form=fin_read,UMI_pro=UMI_pro)
    return(fin_dat)
  }
  cl<-makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(stringdist))
  environment(data_1) <- .GlobalEnv
  environment(max_reads_stat) <- .GlobalEnv
  environment(Cell.BC) <- .GlobalEnv
  clusterExport(cl,c('data_1','max_reads_stat','Cell.BC'), envir = environment())
  data_con<-parLapply(cl,as.character(Cell.BC$Var1),function(x)tryCatch(max_reads_stat(x,dat=data_1),error=function(e) NULL))
  stopCluster(cl)
  data_con<-data_con[!sapply(data_con,is.null)]
  data_con<-do.call("rbind",data_con)
  data_con<-data_con[data_con$scar_form!="unknown",]
  
  write.csv(data_con,paste0(outpath,"/consensus.csv"),quote=F,row.names = F)
  
  INDEL_ranges_man<-list()
  for(scarform in data_con$scar_form){
    index=which(data$scar_form==scarform)[1]
    INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
  }
  saveRDS(INDEL_ranges_man,paste0(outpath,"/image.rds"))
  return(INDEL_ranges_man)
}

ScarToINDELRange = function(data,scar,out = F){
  if(nrow(data) == 0){
    print("No available Scar")
    return(0)
  }
  Cell.BC = data.frame(table(data$Cell.BC))
  Cell.BC = Cell.BC[Cell.BC$Freq>1,]
  data = data[data$Cell.BC %in% Cell.BC$Var1,]
  data_1 = data[,c("Cell.BC","UMI","scar_f")]
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  
  max_reads_stat<-function(x,dat){
    temreads<-dat[dat$Cell.BC==x,]
    read_data<-data.frame(table(temreads$scar_f))
    read_data<-read_data[order(-read_data$Freq),]
    if(read_data$Freq[1]==1){
      fin_read="unknown"
      reads_pro=0
      reads_num=0
      #UMI=0
      del=NA
      ins=NA
    }else{
      #reads_pro=read_data$Freq[1]/sum(read_data$Freq)
      #fin_read=as.character(read_data$Var1[1])
      #UMI=length(unique(temreads$UMI[temreads$scar_f==fin_read]))
      #reads_num=read_data$Freq[1]
      #if(reads_pro<0.6){
      scarstrdist<-stringdistmatrix(as.character(read_data$Var1),as.character(read_data$Var1))
      scarindex<-which(apply(scarstrdist,1,function(x){sum(read_data$Freq[which(x<9)])>(sum(read_data$Freq)/2)}))
      fin_read<-consensusString(temreads$scar_f[temreads$scar_f %in% as.character(read_data$Var1)[scarindex]])
      reads_pro=round(length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))/sum(read_data$Freq),4)
      reads_num=length(temreads$scar_f %in%as.character(read_data$Var1)[scarindex])
      #UMI=length(unique(temreads$UMI[temreads$scar_f %in% as.character(read_data$Var1)[scarindex]]))
      fin_read<-gsub("\\?","",fin_read)
      #}
      s1<-DNAString(fin_read)
      alig<-pairwiseAlignment(scar,s1,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension = 1)
      stopifnot(is(alig, "PairwiseAlignments"), length(x) == 1L)
      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      p <- strsplit(p,"")[[1]]
      s <- strsplit(s,"")[[1]]
      lenp <- length(p)
      index <- 1
      i <- 0
      del_flag <- F
      in_flag <- F
      del_start<-c()
      del_end<-c()
      ins_start<-c()
      ins_width<-c()
      ins_editseq<-c()
      while(i < lenp){
        i <- i + 1
        if(s[[i]] == '-'){
          if(!del_flag){
            del_flag <- T
            del_start<-c(del_start,index)
            width <- 1
          }
          else{
            width <- width + 1
          }
        }
        else{
          if(del_flag){
            del_flag <- F
            del_end<-c(del_end,index)
            #print(paste("del stop width", width))
          }
        }
        if(p[[i]] == '-'){
          if(!in_flag){
            in_flag <- T
            ins_start<-c(ins_start,index)
            #print(paste("in start", index))
            width <- 1
            
            editseq<-s[[i]]
          }
          else{
            width <- width + 1
            
            editseq<-paste0(editseq,s[[i]])
          }
        }
        else{
          if(in_flag){
            in_flag <- F
            ins_width<-c(ins_width,width)
            ins_editseq<-c(ins_editseq,editseq)
            #print(paste("in stop width", width))
          }
        }
        if(p[[i]] != '-')
          index <- index + 1
      }
      if(del_flag){
        del_flag <- F
        del_end<-c(del_end,index)
      }
      if(in_flag){
        in_flag <- F
        ins_width<-c(ins_width,width)
        ins_editseq<-c(ins_editseq,editseq)
        #print(paste("in stop width", width))
      }
      ins<-IRanges(start = ins_start,width = ins_width)
      mcols(ins)$seq <- ins_editseq
      del<-IRanges(start = del_start,end = del_end)
    }
    fin_dat<-data.frame(UMI=x,scar=fin_read,reads_pro=reads_pro,reads_num=reads_num)
    return(list(fin_dat,del,ins))
  }
  
  cl = makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(stringdist))
  environment(data_1) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  environment(max_reads_stat) <- .GlobalEnv
  environment(Cell.BC) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  clusterExport(cl,c('data_1','scar','mat','max_reads_stat','Cell.BC'),envir=environment())
  INDEL_ranges = parLapply(cl,as.character(Cell.BC$Var1),function(x)tryCatch(max_reads_stat(x,dat=data_1), error=function(e) NULL))
  stopCluster(cl)
  INDEL_ranges = INDEL_ranges[!sapply(INDEL_ranges,is.null)]
  saveRDS(INDEL_ranges,paste0(out,"/image.rds"))

  return(INDEL_ranges)
}

#scar data stat
ScarFilt = function(INDEL_ranges){
  read_counts = do.call("rbind",sapply(INDEL_ranges,function(x){return(x[1])}))
  INDEL_ranges = INDEL_ranges[read_counts$scar!="unknown" & read_counts$reads_num>1]
  return(INDEL_ranges)
}

ScarDf = function(INDEL_ranges){
  read_counts = do.call("rbind",sapply(INDEL_ranges,function(x){return(x[1])}))
  read_counts = read_counts[read_counts$scar!="unknown" & read_counts$reads_num>1,]
  write.table(read_counts,"new_max_res.txt",quote=F,row.names=F)
  write.table(read_counts,"new_max_res_full.txt",quote=F,row.names=F)
  return(read_counts)
}

#edit number stat
EditStat = function(INDEL_ranges, read_counts, out = F){
  edit_num_stat = function(x){
    dat = c(x[[2]],x[[3]])
    return(length(dat))
  }
  edit_num = unlist(lapply(INDEL_ranges,edit_num_stat))
  edit_num = data.frame(read_counts$UMI,edit_num)
  write.csv(edit_num,"edit_num.csv",quote=F,row.names=F)

  return(edit_num)
}

#indel split
INDELSplitDF = function(INDEL_ranges){
  del_ranges = unlist(lapply(INDEL_ranges,function(x){x[2]}))
  ins_ranges = unlist(lapply(INDEL_ranges,function(x){x[3]}))
  del_ranges_df = del_ranges[[1]]
  for(i in 2:length(del_ranges)){
    del_ranges_df = c(del_ranges_df,del_ranges[[i]])
  }
  ins_ranges_df = ins_ranges[[1]]
  for(i in 2:length(ins_ranges)){
    ins_ranges_df = c(ins_ranges_df,ins_ranges[[i]])
  }
  return(list("del_ranges" = del_ranges,"ins_ranges" = ins_ranges,
              "del_ranges_df" = del_ranges_df,"ins_ranges_df" = ins_ranges_df))
}

#plot indel size
INDELSizePlot = function(indel_r,out = F){
  if(length(indel_r$del_ranges_df) > 0){
    del_tab = data.frame(table(width(indel_r$del_ranges_df)))
  }else{
    del_tab = 0
  }
  if(length(indel_r$ins_ranges_df) > 0){
    ins_tab = data.frame(table(width(indel_r$ins_ranges_df)))
    ins_tab$Freq =  -(ins_tab$Freq)
  }else{
    ins_tab = 0
  }
  
  if(is.null(nrow(del_tab)) != T){
    del_tab$type = "Deletion"
  }
  if(is.null(nrow(ins_tab)) != T){
    ins_tab$type = "Insertion"
  }
  
  del_tab$Var1 = factor((as.numeric(as.character(del_tab$Var1))-1))
  tab = rbind(del_tab,ins_tab)
  
  
  tab$Var1 = as.numeric(as.character(tab$Var1))
  p = ggplot(data=tab, aes(x=Var1, y=Freq, fill=type)) + 
    geom_bar(stat="identity")+theme_bw() +
    scale_fill_manual(values=c('lightpink1','lightblue2')) +
    xlab("Indel size (bp)")+ ylab("Occurrences") +
    scale_x_continuous(limits=c(0,max(tab$Var1)+1),breaks=seq(0,max(tab$Var1)+1,20))
  print(p)
  ggsave(p,filename = out,width = 10,height = 6)
  return(p)
}

#delete frequency ribbon plot
INDELRibbonPlot = function(del_ranges_df,scarseg,out){
  all_site_stat = function(i,dat){
    return(c(start(dat[i]):end(dat[i])))
  }
  all_site_per_stat = function(i,dat){
    return(length(which(dat==i))/length(dat))
  }
  
  del_allsite = unlist(lapply(1:length(del_ranges_df),all_site_stat,dat=del_ranges_df))
  del_allsite_per = data.frame("Site" = unique(del_allsite),
                               "Freq" = unlist(lapply(unique(del_allsite),all_site_per_stat,dat=del_allsite)))
  siteplus = c(1:239)[!c(1:239 %in% del_allsite_per$Site)]
  siteplus = data.frame(Site=siteplus,Freq=rep(0,length(siteplus)))
  del_allsite_per = rbind(del_allsite_per,siteplus)
  p = ggplot()+geom_ribbon(data = scarseg,aes(x=scar,ymin=0,ymax=max(del_allsite_per$Freq),fill=type),alpha=0.1) +
    geom_line(data = del_allsite_per,aes(x = Site,y = Freq),size=2,colour="lightskyblue") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks=c()) +
    xlab("")+ylab("Frequency")+theme_bw()
  print(p)
  ggsave(p, filename = out,width = 8,height = 2.5)
  return(p)
}

#change scar form and plot density
INDELToAllReadCounts = function(read_counts,INDEL_ranges,scarseg){
  change_form_stat = function(indel){
    indel = indel[c(2,3)]
    indel = unlist(indel)
    if(length(indel)==0){
      return("unknown")
    }else{
      ins = data.frame(indel[2])
      ins$seq=as.character(indel[[2]]@elementMetadata$seq)
      site_ins = apply(ins,1,function(x){c(x[1]:x[2])})
      if(dim(ins)[1]==1){
        site_ins = list(as.numeric( site_ins))
      }
      cutsite_ins = lapply(site_ins,function(x){unique(scarseg$type[scarseg$scar %in% x])})
      tag_ins = apply(ins,1,function(x){paste0(x[3],"I+",x[1],x[4])})
      tag_ins = lapply(tag_ins,function(x){rep(x,length(cutsite_ins[[which(tag_ins==x)]]))})
      tag_ins = unlist(tag_ins)
      cutsite_ins = unlist(cutsite_ins)
      del = data.frame(indel[1])
      site_del = apply(del,1,function(x){c(x[1]:x[2])})
      if(dim(del)[1]==1){
        site_del = list(as.numeric( site_del))
      }
      cutsite_del = lapply(site_del,function(x){unique(scarseg$type[scarseg$scar %in% x])})
      tag_del = apply(del,1,function(x){paste0((x[3]-1),"D+",x[1])})
      tag_del = lapply(tag_del,function(x){rep(x,length(cutsite_del[[which(tag_del==x)]]))})
      tag_del = unlist(tag_del)
      cutsite_del = unlist(cutsite_del)
      tag = c(tag_del,tag_ins)
      cutsite = c(cutsite_del,cutsite_ins)
      tag_all = rep("NONE",9)
      if(length(tag)==0){
        return(paste(tag_all,collapse = "_"))
      }else{
        for(x in c(1:length(tag))){
          if(tag_all[as.numeric(cutsite[x])]=="NONE"){
            tag_all[as.numeric(cutsite[x])] = tag[x]
          }else{
            tag_all[as.numeric(cutsite[x])] = paste(tag_all[as.numeric(cutsite[x])],tag[x],sep="&")
          }
        }
      }
      return(paste(tag_all,collapse = "_"))
    }
  }
  scar_form = unlist(lapply(INDEL_ranges,change_form_stat))
  read_counts$scar_form = scar_form
  scar_tab = data.frame(table(scar_form))
  scar_tab = scar_tab[order(-scar_tab$Freq),]
  read_counts = read_counts[!scar_form %in% c("unknown"),]
  read_counts = read_counts[order(-as.numeric(as.character(read_counts$reads_num))),]
  read_counts$scar_form = gsub(" ","",read_counts$scar_form)
  write.table(read_counts,"allReadCounts",row.names=F,sep="\t",quote=F)
  return(read_counts)
}

#select consensus scar and write to allReadCounts file
INDELScarRanges = function(INDEL_ranges,read_counts){
  Scar_ranges = INDEL_ranges[!scar_form %in% c("unknown")]
  Scar_ranges = Scar_ranges[order(-as.numeric(as.character(read_counts$reads_num)))]
  names(Scar_ranges) = as.character(read_counts$UMI)
  if(outpath != F){
    save(read_counts,Scar_ranges,file="Scar_ranges.RData")
  }
  return(Scar_ranges)
}

#chord plot
DelChordPlot = function(del_ranges_df,model_names,pos_set,output){
  
  #pos_set = c(1,20,28,47,55,74,82,101,109,128,136,155,163,182)
  
  #set color
  rgb = data.frame("r" = c(0,155,1,243,51,244,211,150,50,230),
                   "g" = c(172,89,152,156,73,207,84,30,130,130),
                   "b" = c(172,182,218,19,94,64,0,100,200,55))
  
  #draw chord
  {
    #Initialize
    data_pos = as.data.frame(del_ranges_df)
    #you can change data size here
    data_pos$from = 0
    data_pos$to = 0
    data_pos$value = 1
    
    #build data for chord drawing
    Classify = function(pos){
      i = 1
      axis = 0
      while (i < length(pos_set)) {
        if(pos>=pos_set[i] && pos<=pos_set[i+1])
        {
          axis = (i+1)/2
        }
        i = i + 1
      }
      return(axis)
    }
    
    for (i in 1:nrow(data_pos)) {
      start = data_pos[i,1]
      end = data_pos[i,2]
      x = Classify(start)
      y = Classify(end)
      if(x != 0 & y != 0){
        data_pos[i,"from"] = model_names[x]
        data_pos[i,"to"] = model_names[y]
        data_pos[i,1] = data_pos[i,1] - pos_set[x*2-1]
        data_pos[i,2] = data_pos[i,2] - pos_set[y*2-1]
      }
    }
    data_pos = data_pos[which(data_pos$from != 0),]
    
    #set color
    grid.col = rgb(rgb$r,rgb$g,rgb$b,max = 255)
    line.col = rgb(rgb$r,rgb$g,rgb$b,alpha=100,max = 255)
    grid.col = grid.col[1:length(model_names)]
    line.col = line.col[1:length(model_names)]
    
    #rgb(231,104,21, maxColorValue = 255)
    df = data.frame("names" = factor(model_names, levels = model_names), 
                    "color" = grid.col,"lcol" = line.col)
    df$color = as.character(df$color)
    df$lcol = as.character(df$lcol)
    
    #plot and save
    circos.clear()
    pdf(file=output, width=8, height=5, pointsize=8)
    circos.initialize(df$names, xlim = c(pos_set[1]-1,pos_set[2]))
    circos.track(ylim = c(0, 1),  track.height=0.1,
                 panel.fun = function(x, y){
                   sector.index = get.cell.meta.data("sector.index")
                   xlim = get.cell.meta.data("xlim")
                   ylim = get.cell.meta.data("ylim")
                   i = get.cell.meta.data("sector.numeric.index")
                   circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                               col = df$color[i], border=df$color[i])
                   circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, 
                               niceFacing = TRUE, col="black", facing="bending.inside", font=2, 
                               adj=c(0.5,0.5))
                   
                 }, bg.border = NA)
    
    for (j in 1:nrow(data_pos)) {
      circos.link(data_pos$from[j], data_pos$start[j], data_pos$to[j], data_pos$end[j], border = FALSE,
                  col = df[which(df$names == data_pos$from[j]), "lcol"])
    }
    dev.off()
  }
  
  
}

