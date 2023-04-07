#stat and plot
library(Biostrings)
library(stringdist)
library(rlist)
library(dplyr) 
library(parallel)
library(reshape2)
library(stringr)

setClass(
  "CRESTData",
  slots = list(samplesheet = "character",
               refN = "numeric",
               # refFasta = "DNAString",
               refFasta = "list",
               cutsite = "list",
               cutsitelong = "list",
               cutsiterange = "list",
               rawCBUMI = "data.frame",
               indels = "list",
               indelPattern = "data.frame",
               finalIndel = "data.frame",
               celltype = "data.frame")
)

#Build Input data
LoadFiles = function(filename){
  #load sample sheet
  filels = read.table(filename,sep = "\t",header = T)
  files = filels$Path
  names(files) = filels$File
  
  #load cutsite
  cutsiteraw = read.table(files["CutSite"])
  cutsiteranges = list()
  cutsites = list()
  cutsitels = list()
  for (i in 1:length(unique(cutsiteraw$V4))) {
    cutsiterawi = cutsiteraw[cutsiteraw$V4 == unique(cutsiteraw$V4)[i],]
    cutsiterange = c("start" = cutsiterawi[1,2], "end" = cutsiterawi[1,3])
    cutsite = cutsiterawi[-1,]
    colnames(cutsite) = c("segid","start","end","group")
    
    #transfer cutsite to long segment
    cutsitel = NULL
    for (j in 1:nrow(cutsite)) {
      scar = cutsite[j,]$start:cutsite[j,]$end
      type = rep(cutsite[j,]$segid,(cutsite[j,]$end-cutsite[j,]$start)+1)
      cutsitel = rbind(cutsitel, data.frame("scar" = scar,"type" = as.character(type)))
    }
    cutsiteranges[[i]] = cutsiterange
    cutsites[[i]] = cutsite
    cutsitels[[i]] = cutsitel
  }
  names(cutsiteranges) = names(cutsites) = names(cutsitels) = unique(cutsiteraw$V4)
  #load refseq
  refseqsraw = read.table(files["Fasta"])
  refseqs = list()
  for (i in 1:as.integer(nrow(refseqsraw)/2)) {
    refseqs[[i]] = DNAString(refseqsraw[i*2,1])
    names(refseqs)[i] = substr(refseqsraw[i*2 - 1,1],2,nchar(refseqsraw[i*2 - 1, 1]))
  }
  
  #load raw CBUMI
  rawCBUMI = read.table(files["CB_UMI"], stringsAsFactors=F, header=T)
  mycrest = new("CRESTData",
                samplesheet = files,
                refN = length(refseqs),
                # refFasta = "DNAString",
                refFasta = refseqs,
                cutsite = cutsites,
                cutsitelong = cutsitels,
                cutsiterange = cutsiteranges,
                rawCBUMI = rawCBUMI)
  return(mycrest)
}
# mycrest = LoadFiles("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/myscript/CRESTLineage/test_data/samplesheet.txt")

FindScar = function(mycrest, cln = 10, match = 1, 
                    mismatch = -3, gapOpening = 6,
                    gapExtension = 1, minscore = NULL){
  #build alignment matrix
  mat  =  nucleotideSubstitutionMatrix(match = match, mismatch = mismatch)
  #get threshold of each reference
  # calThred = function(refFasta, cutrange){
  #   return(length(refFasta) - (cutrange["end"] - cutrange["start"]) + 10)
  # }
  refFasta = mycrest@refFasta
  cutrange = mycrest@cutsiterange
  # point1 = length(scarfull1)- 2*(scar1["end"] - scar1["start"])-6-50
  
  thredls = NULL
  for (i in 1:length(refFasta)) {
    thredls = c(thredls, length(refFasta[[i]])- 2*(cutrange[[i]]["end"] - cutrange[[i]]["start"]) - 6 - 50)
  }
  if(!is.null(minscore)){
    thred = minscore
  }else{
    thred = min(thredls)
  }
  
  align_to_range = function(p,s,cut){
    pn = strsplit(p,"")[[1]]
    sn = strsplit(s,"")[[1]]
    lenp = length(pn)
    index = 1
    i = 0
    del_flag = F
    in_flag = F
    del_start=c()
    del_end=c()
    ins_start=c()
    ins_width=c()
    ins_editseq=c()
    while(i < lenp){
      i = i + 1
      if(sn[[i]] == '-'){
        if(!del_flag){
          del_flag = T
          del_start = c(del_start,index)
          width = 1
        }
        else{
          width = width + 1
        }
      }
      else{
        if(del_flag){
          del_flag = F
          del_end = c(del_end,index)
          #print(paste("del stop width", width))
        }
      }
      if(pn[[i]] == '-'){
        if(!in_flag){
          in_flag = T
          ins_start = c(ins_start,index)
          width = 1
          
          editseq=sn[[i]]
        }
        else{
          width = width + 1
          
          editseq = paste0(editseq,sn[[i]])
        }
      }
      else{
        if(in_flag){
          in_flag = F
          ins_width = c(ins_width,width)
          ins_editseq = c(ins_editseq,editseq)
          #print(paste("in stop width", width))
        }
      }
      if(pn[[i]] != '-')
        index = index + 1
    }
    if(del_flag){
      del_flag = F
      del_end = c(del_end,index)
    }
    if(in_flag){
      in_flag = F
      ins_width = c(ins_width,width)
      ins_editseq = c(ins_editseq,editseq)
      #print(paste("in stop width", width))
    }
    
    ins_start = ins_start-cut
    ins_end = ins_start + ins_width
    del_start = del_start-cut
    del_end = del_end-cut
    
    ins = IRanges(start = ins_start,end = ins_end)
    mcols(ins)$seq = ins_editseq
    del = IRanges(start = del_start,end = del_end)
    
    
    return(list("del" = del,"ins" = ins))
    
  }
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  findBarcode = function(seq){
    scores = NULL
    aligls = list()
    for (i in 1:mycrest@refN) {
      s3 = DNAString(as.character(seq))
      alig = pairwiseAlignment(refFasta[[i]],
                               s3,
                               substitutionMatrix = mat,
                               type="global-local", 
                               gapOpening = gapOpening,
                               gapExtension= gapExtension)
      aligls[[i]] = alig
      scores = c(scores, score(alig))
    }
    
    maxscore = max(scores)
    maxi = which.max(scores)
    if(maxscore <= thred){
      r_read = "unknown"
      r_scar = "unknown"
      del = NA
      ins = NA
      group = NA
    }else{
      group = names(refFasta)[maxi]
      refshort = subseq(as.character(refFasta[[maxi]]),cutrange[[maxi]]["start"],cutrange[[maxi]]["end"])
      align_pos = matchPattern(refshort, as.character(pattern(aligls[[maxi]])), with.indels = T, max.mismatch = 10)
      r_read = as.character(subject(aligls[[maxi]]))
      if(length(align_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(aligls[[maxi]])),start=start(align_pos),end=end(align_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(aligls[[maxi]], "PairwiseAlignments"), length(seq) == 1L)
      p = as.character(alignedPattern(aligls[[maxi]])[[1L]])
      s = as.character(alignedSubject(aligls[[maxi]])[[1L]])
      delins = align_to_range(p,s,cutrange[[maxi]]["start"])
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat = data.frame(new.reads = r_read, scar.BC = r_scar, type = group)
    return(list("del" = del,"ins" = ins, fin_dat))
  }
  
  #parallel process
  cl = makeCluster(cln,type = "FORK")
  clusterEvalQ(cl,library(Biostrings))
  environment(thred) = .GlobalEnv
  environment(refFasta) = .GlobalEnv
  environment(cutrange) = .GlobalEnv
  environment(mat) = .GlobalEnv
  environment(mycrest) = .GlobalEnv
  environment(gapExtension) = .GlobalEnv
  environment(gapOpening) = .GlobalEnv
  environment(findBarcode) = .GlobalEnv
  environment(testFunction) = .GlobalEnv
  environment(align_to_range) = .GlobalEnv
  clusterExport(cl,c('thred','refFasta','cutrange',
                     'mycrest','gapExtension','gapOpening','mat','findBarcode',
                     'testFunction','align_to_range'),
                envir = environment())
  rawindel = parLapply(cl, mycrest@rawCBUMI$Read.Seq, findBarcode)
  stopCluster(cl)
  
  #output
  scarcol = do.call("rbind",sapply(rawindel,function(x){return(x[3])}))
  resdata = cbind(mycrest@rawCBUMI, scarcol)
  resdataf = resdata[resdata$scar.BC != "unknown",]
  rawindel = rawindel[resdata$scar.BC != "unknown"]
  resdataf$scar_f = gsub("[-]", "", as.character(resdataf$scar.BC))
  mycrest@indels = list("INDEL" = rawindel,"Scar" = resdataf)
  # write.table(resdataf,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  # saveRDS(rawindel,"indel.rds")
  
  return(mycrest)
}
# mycrest = FindScar(mycrest)

#data preprocess
INDELChangeForm = function(mycrest, cln = 10){
  INDEL_ranges = mycrest@indels$INDEL
  Scar =  mycrest@indels$Scar
  cutsitelong = mycrest@cutsitelong
  
  #change single indel form
  change_form_stat = function(indel){
    indel = indel[c(1,2)]
    indel = unlist(indel)
    if(length(indel)==0){
      return("unknown")
    }else{
      ins = data.frame(indel[2])
      ins$seq=as.character(indel[[2]]@elementMetadata$seq)
      site_ins = apply(ins,1,function(x){c(x[1]:x[2])})
      if(dim(ins)[1]==1){
        site_ins = list(as.numeric(site_ins))
      }
      cutsite_ins = lapply(site_ins,function(x){unique(scarref$type[scarref$scar %in% x])})
      tag_ins = apply(ins,1,function(x){paste0(x[3],"I+",x[1],x[4])})
      tag_ins = lapply(tag_ins,function(x){rep(x,length(cutsite_ins[[which(tag_ins==x)]]))})
      tag_ins = unlist(tag_ins)
      cutsite_ins = unlist(cutsite_ins)
      del = data.frame(indel[1])
      site_del = apply(del,1,function(x){c(x[1]:x[2])})
      if(dim(del)[1]==1){
        site_del = list(as.numeric( site_del))
      }
      cutsite_del = lapply(site_del,function(x){unique(scarref$type[scarref$scar %in% x])})
      tag_del = apply(del,1,function(x){paste0((x[3]-1),"D+",x[1])})
      tag_del = lapply(tag_del,function(x){rep(x,length(cutsite_del[[which(tag_del==x)]]))})
      tag_del = unlist(tag_del)
      cutsite_del = unlist(cutsite_del)
      tag = c(tag_del,tag_ins)
      cutsite = c(cutsite_del,cutsite_ins)
      tag_all = rep("NONE",length(unique(scarref$type)))
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
  cl = makeCluster(cln)
  environment(change_form_stat)  =  .GlobalEnv
  environment(INDEL_ranges)  =  .GlobalEnv
  environment(cutsitelong)  =  .GlobalEnv
  
  types = names(cutsitelong)
  Scar$scar_form = ""
  for (i in 1:length(types)) {
    scarref = cutsitelong[names(cutsitelong) == types[i]][[1]]
    environment(scarref)  =  .GlobalEnv
    clusterExport(cl,c('INDEL_ranges','change_form_stat',"cutsitelong","scarref"), 
                  envir = environment())
    scar_form_p = parLapply(cl,INDEL_ranges[Scar$type == types[i]],
                            change_form_stat)
    scar_form = unlist(scar_form_p)
    scar_form = gsub(" ","",scar_form)
    Scar[Scar$type == types[i],]$scar_form = scar_form
  }
  stopCluster(cl)
  
  mycrest@indelPattern = Scar
  return(mycrest)
}
# mycrest = INDELChangeForm(mycrest)
# 
# for (oi in unique(mycrest@indelPattern$type)) {
#   scari = mycrest@indelPattern[mycrest@indelPattern$type == oi,]
#   dir.create(paste0(outpath,"/",oi))
#   write.csv(scari,paste0(outpath,"/",oi,"/indel_pattern.csv"),quote=F,row.names = F)
# }
INDELCons = function(mycrest){
  data = mycrest@indelPattern
  Cell.BC = data.frame(table(data$Cell.BC))
  Cell.BC = Cell.BC[Cell.BC$Freq>1,]
  data = data[data$Cell.BC %in% Cell.BC$Var1,]
  data_1 = data[,c("Cell.BC","UMI","scar_f","scar_form")]
  max_reads_stat = function(x,dat){
    temreads = dat[dat$Cell.BC==x,]
    #consensus
    read_data = data.frame(table(as.character(temreads$scar_f)))
    read_data = read_data[order(-read_data$Freq),]
    scar_data = data.frame(table(as.character(temreads$scar_form)))
    scar_data = scar_data[order(-scar_data$Freq),]
    
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
    fin_line = data.frame("Cell.BC" = x,
                          # "cons" = pat.cons,
                          "main" = pat.main,"umim" = pat.umi,
                          # "reads_pro.cons" = reads_pro_cons,
                          "reads_pro.main" = reads_pro_main,
                          "reads_pro.umim" = reads_pro_umi,
                          "umi_pro" = umi_pro,
                          "reads_num" = reads_num,
                          "umi_num" = umi_num, stringsAsFactors = F)
    return(fin_line)
  }
  
  # 
  data_con = NULL
  for (i in 1:length(as.character(Cell.BC$Var1))) {
    line = max_reads_stat(as.character(Cell.BC$Var1)[i],dat=data_1)
    data_con = rbind(data_con,line)
  }
  mycrest@finalIndel = data_con
  return(mycrest)
  
}
# mycrest = INDELCons(mycrest)

