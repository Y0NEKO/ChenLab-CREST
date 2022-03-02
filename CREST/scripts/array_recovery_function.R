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
library(pheatmap)
library(reshape2)
library(ggtree)

library(ggvenn)
library(data.tree)
library(networkD3)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(ggtree)
library(stringr)
library(ggstance)
library(ggnewscale)
library(grid)
library(gtable)
library(rlist)
library(purrr)

#input data
ReadSheet = function(filename){
  filels = read.table(filename,sep = "\t",header = T)
  files = filels$Path
  names(files) = filels$File
  return(files)
}

ReadCutsite = function(segref,reftype=NULL){
    colnames(segref) = c("indx","start","end")   
    scar = NULL
    type = NULL
    if(is.null(reftype)){
      for (i in 1:nrow(segref)) {
        scar = c(scar,segref[i,]$start:segref[i,]$end)
        type = c(type,rep(segref[i,]$indx,(segref[i,]$end-segref[i,]$start)+1))
      }
      scarseg = data.frame("scar" = scar,"type" = as.character(type))   
    }else{
      endsite<-NA
      for (i in 2:nrow(segref)) {
        endsite<-c(endsite,segref[["end"]][i]+segref[["start"]][1])
      }
      endsite[nrow(segref)]<-segref[["end"]][1]
      scarseg = data.frame("scar" = c(1:endsite[nrow(segref)]),"type" = NA)
      #endsite<-is.na(endsite)
      for (i in rev(endsite)) {
        if(!is.na(i)){
          scarseg$type[1:i]<-which(endsite==i)-1
         }else{
          break
         }
      }      
    }
     return(scarseg)
}

ReadFasta = function(filename){
  sv = read.table(filename)
  scarfull = DNAString(sv[2,1])
  return(c("scarfull" = scarfull))
}

#find_barcode
align_to_range = function(p,s,cut){
  pn <- strsplit(p,"")[[1]]
  sn <- strsplit(s,"")[[1]]
  lenp <- length(pn)
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
    if(sn[[i]] == '-'){
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
    if(pn[[i]] == '-'){
      if(!in_flag){
        in_flag <- T
        ins_start<-c(ins_start,index)
        width <- 1
        
        editseq<-sn[[i]]
      }
      else{
        width <- width + 1
        
        editseq<-paste0(editseq,sn[[i]])
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
    if(pn[[i]] != '-')
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
  
  ins_start = ins_start-cut
  ins_end = ins_start + ins_width
  del_start = del_start-cut
  del_end = del_end-cut
  
  ins<-IRanges(start = ins_start,end = ins_end)
  mcols(ins)$seq <- ins_editseq
  del<-IRanges(start = del_start,end = del_end)
  
  
  return(list("del" = del,"ins" = ins))
  
}
FindScar = function(data,scarfull,scar,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  point = length(scarfull)- 2*(scar["end"] - scar["start"])-6
  find_barcode<-function(data){
    type = "none"
    s3<-DNAString(as.character(data))
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    scarshort = subseq(as.character(scarfull),scar["start"],scar["end"])
    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos = matchPattern(scarshort,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read = as.character(subject(alig))
      if(length(scar_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
      
      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      delins = align_to_range(p,s,scar["start"])
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat<-data.frame(new.reads=r_read, scar.BC=r_scar,type=type)
    return(list("del" = del,"ins" = ins,fin_dat))
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
  clusterExport(cl,c('point','scarfull','scar','data','mat','find_barcode','testFunction',"align_to_range"),envir = environment())
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
  data_v1<-data[,-4]
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1))
  
}

FindScarBulk = function(data,scarfull,scar,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  point = length(scarfull)- 2*(scar["end"] - scar["start"])-6
  find_barcode<-function(data){
    s3<-DNAString(as.character(data))
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    scarshort = subseq(as.character(scarfull),scar["start"],scar["end"])
    type = "none"
    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos = matchPattern(scarshort,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read = as.character(subject(alig))
      if(length(scar_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
      
      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      delins = align_to_range(p,s,scar["start"])
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat<-data.frame(new.reads=r_read,scar.BC=r_scar,type=type)
    return(list("del" = del,"ins" = ins,fin_dat))
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
  environment(align_to_range) <- .GlobalEnv
  clusterExport(cl,c('point','scarfull','scar','data','mat','find_barcode','testFunction',"align_to_range"),
                envir = environment())
  scar_BC<-parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)
  
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data[,-4]
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1))
}



SplitScar = function(data,scarfull1,scarfull2,scar1,scar2,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  
  point1 = length(scarfull1)- 2*(scar1["end"] - scar1["start"])-6-50
  point2 = length(scarfull2)- 2*(scar2["end"] - scar2["start"])-6-50
  find_barcode<-function(data){
    s3<-DNAString(as.character(data))
    alig1<-pairwiseAlignment(scarfull1,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    alig2<-pairwiseAlignment(scarfull2,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    if(score(alig1)>score(alig2)){
      alig<-alig1
      type<-"V1"
      point<-point1
      scar = scar1
      scarshort = subseq(as.character(scarfull1),scar1["start"],scar1["end"])
    }else if(score(alig1)<score(alig2)){
      alig<-alig2
      type<-"V2"
      point<-point2
      scar = scar2
      scarshort = subseq(as.character(scarfull2),scar2["start"],scar2["end"])
    }else{
      alig<-alig1
      r_read<-"unknown"
      r_scar<-"unknown"
      type<-"unknown"
      point<-point1
      scar = scar1
      scarshort = subseq(as.character(scarfull1),scar1["start"],scar1["end"])
      del=NA
      ins=NA
    }
    if(score(alig)<=point){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos = matchPattern(scarshort,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read = as.character(subject(alig))
      if(length(scar_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
      
      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      delins = align_to_range(p,s,scar["start"])
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat<-data.frame(new.reads=r_read,scar.BC=r_scar,type=type)
    return(list("del" = del,"ins" = ins,fin_dat))
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
  clusterExport(cl,c('point1','point2','scarfull1','scarfull2','scar1','scar2',
                     'data','mat','find_barcode','testFunction',"align_to_range"),envir = environment())
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
    tag_all<-rep("NONE",length(unique(scarref$type)))
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
INDELChangeForm = function(INDEL_ranges,data,scarref,outpath,cln){
  cl<-makeCluster(cln)
  environment(change_form_stat) <- .GlobalEnv
  environment(INDEL_ranges) <- .GlobalEnv
  environment(scarref) <- .GlobalEnv
  clusterExport(cl,c('INDEL_ranges','change_form_stat',"scarref"), envir = environment())
  scar_form_p<-parLapply(cl,INDEL_ranges,change_form_stat)
  stopCluster(cl)
  scar_form<-unlist(scar_form_p)
  scar_form<-gsub(" ","",scar_form)
  data$scar_form<-scar_form
  write.csv(data,paste0(outpath,"/indel_pattern.csv"),quote=F,row.names = F)
  return(data)
  
}

INDELCons = function(INDEL_ranges,data,scarref,outpath,scarfull,scar,cln){
  
  Cell.BC<-data.frame(table(data$Cell.BC))
  Cell.BC<-Cell.BC[Cell.BC$Freq>1,]
  data<-data[data$Cell.BC %in% Cell.BC$Var1,]
  data_1<-data[,c("Cell.BC","UMI","scar_f","scar_form")]
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  # x = as.character(Cell.BC$Var1)[2]
  # dat =data_1
  max_reads_stat = function(x,dat){
    temreads = dat[dat$Cell.BC==x,]
    
    #consensus
    read_data = data.frame(table(as.character(temreads$scar_f)))
    read_data = read_data[order(-read_data$Freq),]
    scar_data = data.frame(table(as.character(temreads$scar_form)))
    scar_data = scar_data[order(-scar_data$Freq),]
    
    if(scar_data$Freq[1]==1){
      reads_num=0
      #UMI=0
      del=NA
      ins=NA
    }else{
      #consensus:
      scarstrdist = stringdistmatrix(as.character(read_data$Var1),as.character(read_data$Var1))
      scarindex = which(apply(scarstrdist,1,function(x){sum(read_data$Freq[which(x<9)])>(sum(read_data$Freq)/3)}))
      if(length(scarindex) > 0){
        fin_read = consensusString(temreads$scar_f[temreads$scar_f %in% as.character(read_data$Var1)[scarindex]])
        reads_pro_cons = round(length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))/sum(read_data$Freq),4)
        fin_read_cons = gsub("\\?","",fin_read)
        s1 = DNAString(fin_read_cons)
        aligc = pairwiseAlignment(scarfull,s1,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension = 1)
        stopifnot(is(aligc, "PairwiseAlignments"), length(x) == 1L)
        p <- as.character(alignedPattern(aligc)[[1L]])
        s.cons <- as.character(alignedSubject(aligc)[[1L]])
        indel.cons = align_to_range(p,s.cons,scar["start"])
        pat.cons =  change_form_stat(indel.cons)
      }else{
        reads_pro_cons = 0
        pat.cons = "unknown"
      }
      
      
      #reads main
      pat.main = as.character(scar_data$Var1[1])
      reads_pro_main = round(scar_data$Freq[1]/sum(scar_data$Freq),4)
      
      #umi main
      read_data_umi = unique(temreads)
      read_data_umi = data.frame(table(read_data_umi$scar_form))
      read_data_umi = read_data_umi[order(-read_data_umi$Freq),]
      pat.umi = as.character(read_data_umi$Var1[1])
      umi_pro = round(read_data_umi$Freq[1]/length(unique(temreads)$UMI),4)
      reads_pro_umi = round(scar_data[which(scar_data$Var1 == pat.umi),"Freq"]/sum(scar_data$Freq),4)
       
      #
      reads_num = nrow(temreads)
      umi_num = length(unique(temreads)$UMI)
      
      
    }
    fin_line = data.frame("Cell.BC" = x,"cons" = pat.cons,"main" = pat.main,"umim" = pat.umi,
                 "reads_pro.cons" = reads_pro_cons,
                 "reads_pro.main" = reads_pro_main,
                 "reads_pro.umim" = reads_pro_umi,
                 "umi_pro" = umi_pro,
                 "reads_num" = reads_num,
                 "umi_num" = umi_num, stringsAsFactors = F)
    return(fin_line)
  }
  cl<-makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(stringdist))
  environment(data_1) <- .GlobalEnv
  environment(max_reads_stat) <- .GlobalEnv
  environment(Cell.BC) <- .GlobalEnv
  environment(scarref) <- .GlobalEnv
  environment(align_to_range) <- .GlobalEnv
  environment(change_form_stat) <- .GlobalEnv
  environment(scarfull) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  clusterExport(cl,c('data_1','max_reads_stat','Cell.BC',"align_to_range","scarref","change_form_stat","scarfull","mat","scar"), envir = environment())
  data_con<-parLapply(cl,as.character(Cell.BC$Var1),function(x)tryCatch(max_reads_stat(x,dat=data_1),error=function(e) NULL))
  stopCluster(cl)
  
  
  data_con<-data_con[!sapply(data_con,is.null)]
  data_con<-do.call("rbind",data_con)
  
  write.csv(data_con,paste0(outpath,"/final_scarform.csv"),quote=F,row.names = F)
  
  
  INDEL_ranges<-INDEL_ranges[data$Cell.BC %in% Cell.BC$Var1]
  INDEL_ranges_man<-list()
  for(scarform in data_con$cons){
    index=which(data$scar_form==scarform)[1]
    INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
  }
  saveRDS(INDEL_ranges_man,paste0(outpath,"/image_cons.rds"))
  
  INDEL_ranges_man<-list()
  for(scarform in data_con$main){
    index=which(data$scar_form==scarform)[1]
    INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
  }
  saveRDS(INDEL_ranges_man,paste0(outpath,"/image_main.rds"))
  
  INDEL_ranges_man<-list()
  for(scarform in data_con$umim){
    index=which(data$scar_form==scarform)[1]
    INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
  }
  saveRDS(INDEL_ranges_man,paste0(outpath,"/image_umim.rds"))
  
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
  del_ranges = unlist(lapply(INDEL_ranges,function(x){x[1]}))
  ins_ranges = unlist(lapply(INDEL_ranges,function(x){x[2]}))
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
  del_allsite = del_allsite[which(del_allsite >= segref[1,2] & del_allsite <= segref[nrow(segref),3])]
  del_allsite_per = data.frame("Site" = unique(del_allsite),
                               "Freq" = unlist(lapply(unique(del_allsite),all_site_per_stat,dat=del_allsite)))
  siteplus = c(segref[1,2]:segref[7,3])[!c(segref[1,2]:segref[7,3] %in% del_allsite_per$Site)]
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

TagDataProcess = function(data,prefix,Cells){
  
  TagStat = function(x) {
    x = as.character(x)
    umi = x[10]
    CB = x[1]
    x = x[11]
    x = unlist(strsplit(x,"_|&"))
    x = x[!x%in%c("NONE")]
    x = unique(x)
    if(length(x) == 0){
      return(NA)
    }else{
      return(data.frame(Cell.BC = CB,Reads_num = umi,Tag = x))
    }
  }
  
  tag = NULL
  common.CB = NULL
  if((length(data)>1)){
    for (i in 1:(length(data)-1)) {
      common.CB = c(common.CB,intersect(data[[i]]$Cell.BC, data[[i+1]]$Cell.BC))
    }
  }else{
    common.CB = data[[1]]$Cell.BC
  }
  
  for (i in 1:length(data)) {
    data[[i]]=data[[i]][data[[i]]$Cell.BC %in% Cells$Cell.BC,]
    data[[i]]=data[[i]][data[[i]]$Cell.BC %in% common.CB,]
    tagi = apply(data[[i]],1,TagStat)
    tagi = do.call("rbind",tagi)
    tagi = na.omit(tagi)
    tabi = data.frame(table(tagi$Tag)/length(as.character(unique(tagi$Cell.BC))))
    #black list filter
    
    tagi$Tag = paste(prefix[i], tagi$Tag, sep = "")
    
    tag = rbind(tag,tagi)
  }
  
  return(tag)
}

ArrayDataProcess = function(data,prefix,Cells){
  
  tag = NULL
  common.CB = NULL
  if((length(data)>1)){
    for (i in 1:(length(data)-1)) {
      common.CB = c(common.CB,intersect(data[[i]]$Cell.BC, data[[i+1]]$Cell.BC))
    }
  }else{
    common.CB = data[[1]]$Cell.BC
  }
  
  if(length(data)>1){
    tagi = Reduce(function(x,y) merge(x = x[,c(1,10,11)], y = y[,c(1,10,11)], by = "Cell.BC"),data)
  }else{
    tagi = data[[1]][,c(1,10,11)]
  }
  tagi = tagi[tagi$Cell.BC %in% Cells$Cell.BC,]
  
  Tag = NULL
  num = 0
  for (i in 1:length(data)) {
    tagi[,1+(2*i)] = paste0(prefix[i],tagi[,1+(2*i)])
    Tag = paste0(Tag,tagi[,1+(2*i)])
    num = num + tagi[,2*i]
    if(i < length(data)){
      Tag = paste0(Tag,"_")
    }
  }
  tag = data.frame(Cell.BC = tagi$Cell.BC,Reads_num = num,Tag = Tag)
  
  return(tag)
}

CellMap1 = function(tag,Cells){
  #1&2 allreadcount 3 celltype
  
  # data_2<-data_2[paste0(data_2$Cell.BC,"_1_2")%in%Cells$X,]
  # data_1<-data_1[data_1$Cell.BC%in% intersect(data_1$Cell.BC,data_2$Cell.BC),]
  # data_2<-data_2[data_2$Cell.BC%in% intersect(data_1$Cell.BC,data_2$Cell.BC),]
  

  tag$celltype = Cells$Cell.type[match(tag$Cell.BC,Cells$Cell.BC)]

  one_jac_stat2<-function(l,VBC1,clone_stat_data){
    VBC2<-as.numeric(clone_stat_data[,l])
    num <- sum(sapply(1:length(VBC1), function(x)(min(VBC1[x],VBC2[x]))))
    den <- sum(sapply(1:length(VBC1), function(x)(max(VBC1[x],VBC2[x]))))
    return(num/den)
  }
  one_jac_stat1<-function(i,clone_stat_data){
    VBC1<-as.numeric(clone_stat_data[,i])
    return(unlist(lapply(clu,one_jac_stat2,VBC1=VBC1,clone_stat_data=clone_stat_data)))
  }
  one_op_stat2<-function(y,x,jac,jac_sample){
    jac_real<-jac[x,y]
    jac_pred<-as.numeric(unlist(lapply(jac_sample,function(z){z[x,y]})))
    ob<-jac_real/mean(jac_pred)
    zscore<-(jac_real-mean(jac_pred))/sd(jac_pred)
    p<-pnorm(zscore,lower.tail = F)
    return(c(ob,p))
  }
  one_op_stat1<-function(x,jac,jac_sample,clu){
    return(lapply(clu,one_op_stat2,x=x,jac=jac,jac_sample=jac_sample))
  }
  
  
  clone_tab<-acast(tag,Tag~celltype)
  clu<-colnames(clone_tab)
  clone_tab<-clone_tab[apply(clone_tab,1,sum)>=2,]
  
  #plot tag size
  pdf("tag_size.pdf",width = 4,height = 3)
  plot(density(log2(table(apply(clone_tab,1,sum)))),main = "Tag size (>=2cells)",xlab = "log2 (cells)")
  dev.off()
  
  
  annotation_col = data.frame(Group=factor(substring(rownames(clone_tab),first=1,last = 2)))
  row.names(annotation_col)<-rownames(clone_tab)
  pheatmap(t(clone_tab),cluster_cols = F,cluster_rows = F,show_colnames = F,border=FALSE,annotation_col = annotation_col,scale = "row",filename = "tag_heatmap_scale_all.pdf",width = 5,height = 3)
  
  all_jac<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab)),ncol=length(clu)))
  names(all_jac)<-clu
  row.names(all_jac)<-clu
  
  jac_sample<-list()
  for(t in c(1:500)){
    tag$tags_sample<-sample(tag$Tag,length(tag$Tag))
    clone_tab_sample_sub<-acast(tag,tags_sample~celltype)
    jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab_sample_sub)),ncol=length(clu)))
    names(jac_sample_sub)<-clu
    row.names(jac_sample_sub)<-clu
    jac_sample<-c(jac_sample,list(jac_sample_sub))
  }
  
  saveRDS(jac_sample,paste0("jac_sample_all.rds"))
  
  #plot pvalue
  all_ob_p<-lapply(clu,one_op_stat1,jac=all_jac,jac_sample=jac_sample,clu=clu)
  all_ob<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==1]]}))
  all_ob<-data.frame(all_ob)
  all_p<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==0]]}))
  all_p<-data.frame(all_p)
  names(all_ob)<-clu
  row.names(all_ob)<-clu
  names(all_p)<-clu
  row.names(all_p)<-clu
  
  all_ob[is.na(all_ob)] = 0
  pheatmap(all_ob/max(all_ob),border=FALSE,cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("dodgerblue4","white","darkred"))(100),filename = "obp_all.pdf",width = 3.3,height = 3)
  pheatmap(all_p,display_numbers=T,cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(50)),file=paste0("pvalue_all.pdf"),width = 3.3,height = 3)
  pheatmap(all_jac,display_numbers=T,cluster_rows = F,cluster_cols = F,file=paste0("jaccard_all.pdf"),width = 3.3,height = 3)
  
  write.csv(all_ob,paste0("all_ob.csv"),quote=F)
  write.csv(all_p,paste0("all_pvalue.csv"),quote=F)
  write.csv(all_jac,paste0("all_jac.csv"),quote=F)
  
  #tag filt
  sample_list<-list()
  for(i in c(1:500)){
    tag$Tag.sample<-sample(tag$Tag,length(tag$Tag))
    clone_stat_sample<-acast(tag,Tag.sample~celltype)
    sample_list<-c(sample_list,list(clone_stat_sample))
  }
  
  virus_list<-rownames(clone_tab)
  check_stat2<-function(c){
    return(unlist(lapply(virus_list,function(x){check_stat1(x=x,c=c)})))
  }
  
  check_stat1<-function(x,c){
    s_clone<-unlist(lapply(sample_list,function(i){return(i[x,c])}))
    real_clone<-clone_tab[x,c]
    zscore<-(real_clone-mean(s_clone))/sd(s_clone)
    return(pnorm(zscore,lower.tail = F))
  }
  
  check_res<-lapply(clu,check_stat2)
  clone_check<-do.call("cbind",check_res)
  colnames(clone_check)<-clu
  rownames(clone_check)<-virus_list
  clone_tab[!clone_check<=0.05]<-0
  clone_tab<-clone_tab[apply(clone_tab,1,sum)>1,]
  pheatmap(t(clone_tab),cluster_cols = F,cluster_rows = F,show_colnames = F,border=FALSE,annotation_col = annotation_col,scale = "row",filename = "tag_heatmap_scale_filter.pdf",width = 5,height = 3)
  write.csv(clone_tab,paste0("filter_clone_tab.csv"),quote=F)
  
  #plot filtered tag
  clone_num<-apply(clone_tab,1,sum)
  pdf(paste0("tag_size_filt.pdf"),width = 4,height = 3)
  plot(density(log2(clone_num)),main = "Tag size (>=2cells)",xlab = "log2 (cells)")
  dev.off()
  
  #jac filtered
  all_jac<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab)),ncol=length(clu)))
  names(all_jac)<-clu
  row.names(all_jac)<-clu
  jac_sample<-list()
  for(t in c(1:500)){
    tag$tags_sample<-sample(tag$Tag,length(tag$Tag))
    clone_tab_sample_sub<-acast(tag,tags_sample~celltype)
    jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab_sample_sub)),ncol=length(clu)))
    names(jac_sample_sub)<-clu
    row.names(jac_sample_sub)<-clu
    jac_sample<-c(jac_sample,list(jac_sample_sub))
  }
  saveRDS(jac_sample,paste0("jac_sample_filte.rds"))
  
  all_ob_p<-lapply(clu,one_op_stat1,jac=all_jac,jac_sample=jac_sample,clu=clu)
  all_ob<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==1]]}))
  all_ob<-data.frame(all_ob)
  all_p<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==0]]}))
  all_p<-data.frame(all_p)
  names(all_ob)<-clu
  row.names(all_ob)<-clu
  names(all_p)<-clu
  row.names(all_p)<-clu
  
  all_ob[is.na(all_ob)] = 0
  pheatmap(all_ob/max(all_ob),border=FALSE,cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("dodgerblue4","white","darkred"))(100),filename = "obp_filte.pdf",width = 3.3,height = 3)
  pheatmap(all_p,display_numbers=T,cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(50)),file=paste0("pvalue_filte.pdf"),width = 3.3,height = 3)
  pheatmap(all_jac,display_numbers=T,cluster_rows = F,cluster_cols = F,file=paste0("jaccard_filte.pdf"),width = 3.3,height = 3)
  write.csv(all_ob,paste0("filte_ob.csv"),quote=F)
  write.csv(all_p,paste0("filte_pvalue.csv"),quote=F)
  write.csv(all_jac,paste0("filte_jac.csv"),quote=F)

}



#tree build
#build tree
BuildTagTree = function(tag,Cells,outname){
  
  #tag stat
  Tag = data.frame(table(tag$Tag))
  tag$Cell.num = Tag$Freq[match(tag$Tag,Tag$Var1)]
  tag_tab = acast(tag,Cell.BC~Tag)
  tag_tab[!is.na(tag_tab)] = 1
  tag_tab[is.na(tag_tab)] = 0
  
  #tag integrate 
  cell_tab = data.frame(table(tag$Cell.BC))
  cell_tab = cell_tab[order(-cell_tab$Freq),]
  tags_all = lapply(as.character(cell_tab$Var1),function(x){sort(as.character(tag$Tag[tag$Cell.BC == x]))})
  tags_paste = sapply(tags_all,function(x){paste(x,collapse = "_")})
  tags_tab = data.frame(table(tags_paste))
  tags_tab$num = unlist(lapply(as.character(tags_tab$tags_paste),function(x){length(strsplit(x,split = "_")[[1]])}))
  tags_tab = tags_tab[order(-tags_tab$num,-tags_tab$Freq),]
  tags_uni = lapply(as.character(tags_tab$tags_paste),function(x){strsplit(x,split = "_")[[1]]})
  Tag_1 = Tag[Tag$Var1 %in% unlist(tags_uni[tags_tab$num == 1]),]
  Tag_1 = Tag_1[order(-Tag_1$Freq),]
  tags_uni[tags_tab$num == 1] = as.list(as.character(Tag_1$Var1))
  tags_tab[tags_tab$num == 1,] = tags_tab[tags_tab$num == 1,][match(as.character(Tag_1$Var1),
                                                                    as.character(tags_tab$tags_paste[tags_tab$num==1])),]
  
  #node build
  cluster_stat = function(i){
    x = tags_uni[[i]]
    n = tags_tab$num[i]
    tags_belone = NA
    for(y_ind in which(tags_tab$num<n)){
      y=tags_uni[[y_ind]]
      if(length(intersect(x,y))>0 & length(setdiff(y,x))==0){
        tags_belone = y_ind
        break
      }else{
        next
      }
    }
    return(tags_belone)
  }
  
  belons = sapply(which(tags_tab$num > 1),cluster_stat)
  belons[tags_tab$num == 1]= NA
  
  Tags = as.list(which(tags_tab$num == 1))
  names(Tags) = which(tags_tab$num == 1)
  nodes = as.list(c(1:length(tags_tab$num)))
  
  for(i in c(1:length(tags_tab$num))){
    n = belons[[i]]
    while(!is.na(n)){
      nodes[[i]] = c(n,nodes[[i]])
      n = belons[[n]]
    }
  }
  nodes_len = sapply(nodes,length)
  
  for(i in c(1:length(tags_tab$num))){
    nodes[[i]] = as.character(tags_tab$tags_paste)[nodes[[i]]]
    if(length(nodes[[i]])<max(nodes_len)){
      nodes[[i]][(length(nodes[[i]])+1):max(nodes_len)] = NA
    }else{
      next
    }
  }
  nodes=data.frame(do.call("rbind",nodes))
  names(nodes) = paste("N",as.character(1:ncol(nodes)),sep = "")
  nodes$pathString = do.call(paste, c("N0",nodes, sep="/"))
  nodes$pathString = gsub("/NA","",nodes$pathString)
  
  #save tree figure and rds
  population = as.Node(nodes)
  saveNetwork(diagonalNetwork(ToListExplicit(population, unname = TRUE),
                              margin = 10,fontSize = 8,height = 7000),
              file = paste0(outname,".html"))
  
  saveRDS(population,paste0(outname,".rds"))
  
  #save celltype tab
  cell_tab$tags = tags_paste
  cell_tab$celltype = Cells$Cell.type[match(as.character(cell_tab$Var1),Cells$Cell.BC)]
  write.csv(cell_tab,paste0(outname,"_cell_tab.csv"),row.names = F,quote = F)
  return(list(population,cell_tab))
  
}
#plot tag tree
PlotTagTree = function(tree,cm,prefix,outname){
  #0.require packages
  # library(data.tree)
  # library(dplyr)
  # library(networkD3)
  # library(ggplot2)
  # library(ggtree)
  # library(reshape2)
  # library(stringr)
  # library(ggstance)
  # library(grid)
  # library(gtable)
  # library(rlist)
  # library(purrr)
  
  #1.trans data.tree to newick with internal node
  td = ToDataFrameNetwork(tree)
  from_node = unique(td$from)
  td[which(td$to %in% from_node),"to"] = paste("node.",td[which(td$to %in% from_node),"to"],sep = "")
  td = rbind(data.frame("from" = from_node, "to" = from_node),td)
  td$from = paste("node.",td$from,sep = "")
  td = td[which(td$to != "N0"),]
  td_node = FromDataFrameNetwork(td)
  tree_nwk = ToNewick(td_node)
  write(tree_nwk,paste(outname,"nwk",sep = "."))
  
  
  #2.data build
  #(1)tree_data for plotting tree
  tree_data = read.tree(paste(outname,"nwk",sep = "."))
  tl = tree_data$tip.label
  
  #(2)tree_cell_plot for plotting cell composition 
  tree_cell = cm %>% group_by(tags, celltype) %>% summarise(count = sum(Freq))
  tree_cell = dcast(tree_cell,tags~celltype,value.var = "count")
  tree_cell[is.na(tree_cell)] = 0
  #normlize
  tree_cell[,-1] = log(tree_cell[,-1]+1)
  tree_cell[,-1] = t(apply(tree_cell[,-1], 1, function(x) {x/sum(x)}))
  tree_cell_l = melt(tree_cell)
  tree_cell_plot = tree_cell_l
  
  #option for cell composition bar plot 
  # tree_cell_l = melt(tree_cell)
  # tree_cell_plot = left_join(data.frame(id=tl),tree_cell_l,by=c("id" = "tags"))
  
  #option for 
  # tree_cell_plot$variable = as.character(tree_cell_plot$variable)
  # tree_cell_plot[which(tree_cell_plot$value>0),3] = tree_cell_plot[which(tree_cell_plot$value>0),2]
  # tree_cell_plot[which(tree_cell_plot$value==0),3] = NA
  
  #(3)indel for plotting pattern
  #library(stringr)
  
  indel = data.frame(id=NA,start=NA,end=NA,width=NA,type=NA)
  indel = indel[-1,]
  #build indel data frame
  for(i in c(1:length(tl[-1]))){
    x = tl[i]
    x_str = strsplit(x,split = "_")
    del = NULL
    ins = NULL
    for (k in 1:length(prefix)) {
      x_pre = x_str[[1]][grepl(prefix[k], x_str[[1]])]
      if(length(x_pre)>0){
        x_pre_del = x_pre[grep("D",x_pre)]
        if(length(x_pre_del) > 0){
          x_pre_del = str_extract_all(x_pre_del, "[0-9]+")
          for (j in 1:length(x_pre_del)) {
            line = as.numeric(x_pre_del[[j]])
            del = rbind(del, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                        end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                        type="deletion"))
          }
          
        }
        x_pre_ins = x_pre[grep("I",x_pre)]
        if(length(x_pre_ins) > 0){
          x_pre_ins = str_extract_all(x_pre_ins, "[0-9]+")
          for (j in 1:length(x_pre_ins)) {
            line = as.numeric(x_pre_ins[[j]])
            ins = rbind(ins, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                        end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                        type="insertion"))
          }
        }
        
      }
    }
    indel = rbind(indel,del,ins)
  }
  
  #(4)reorder tree
  SortTree = function(label){
    
    #construct data structure
    taglst = list()
    for (i in 1:length(label)) {
      tchr = label[i]
      tline = strsplit(tchr,"_")
      element = list("labels" = tchr,"tags" = tline[[1]], "tag_number" = length(tline[[1]]))
      taglst[[i]] = element
    }
    taglst.sort = taglst[list.order(taglst,-tag_number)]
    for (i in 2:(length(taglst.sort) - 1)) {
      point = taglst.sort[[i]]
      score = 0
      for (j in (i+1):length(taglst.sort)) {
        queue = taglst.sort[[j]]
        qscore = length(intersect(point$tags,queue$tags))
        
        if(qscore > score){
          score = qscore
          taglst.sort[[j]] = NULL
          taglst.sort = list.insert(taglst.sort,i+1,queue)
        }
      }
    }
    
    
    #result
    label.sort = map_chr(taglst.sort, 1)
  }
  label.sort =  SortTree(tree_data$tip.label)
  tree_data_sort = ape::rotateConstr(tree_data, label.sort)
  
  #3.plot integrated tree
  #set expand can change width of the tree
  p = ggtree(tree_data_sort,size = 0.1, ladderize=F) + geom_tiplab(size = 0.3) + 
    xlim_expand(c(0, 200),panel = "Tree")
  
  p1 = p + geom_facet(panel = "cells",data = tree_cell_plot,
                      geom = geom_tile,
                      mapping = aes(x = as.numeric(as.factor(variable)),fill = value,color = variable)) +
    #    scale_fill_viridis_d(option="D", name="discrete\nvalue") +
    scale_fill_gradient(low = "white",high = "#440130")  +
    new_scale_color() +
    geom_facet(panel = "indel pattern", data = indel,
               geom = geom_segment,
               mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
               size = 0.3)
  p1
  
  # p1 = facet_plot(p,panel = "cell composition",data = tree_cell_plot,
  #                 geom = geom_barh,
  #                 mapping = aes(x = value,fill = variable),
  #                 stat="identity") %>%
  #   facet_plot(panel = "indel pattern", data = indel,
  #              geom = geom_segment,
  #              mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
  #              size = 1)
  
  
  #change facet grid
  gt = ggplot_gtable(ggplot_build(p1))
  gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
  gt$widths[9] = 0.5*gt$widths[9] # in this case it was colmun 7 - reduce the width by a half
  grid.draw(gt) # plot with grid draw
  
  ggsave(gt,filename = paste(outname,"pdf",sep = "."),height = 100,width = 15,units = "cm")
}
#Main
TagTreeMain = function(data,prefix,Cells,outname){
  BuildTagTree(data,Cells,outname)
  tree = readRDS(paste0(outname,".rds")) #tag tree rds file
  cm = read.csv(paste0(outname,"_cell_tab.csv")) #cell matrix
  PlotTagTree(tree,cm,prefix,outname)
}


