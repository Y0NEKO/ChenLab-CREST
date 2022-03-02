#Recover array edit pattern
#2022/03/02 by lhx

source("array_recovery_function.R")

#arg1 command arg2 sheetfile arg3 outpath
args=commandArgs(T)
command = args[1]
sheetfile = args[2]
outpath = args[3]

if(length(args)>3){
  cln = as.numeric(args[4])
}else{
  cln = 2
}

PreData1 = function(filename,scarfull,scar,cln){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  print("Find scar...")
  data = FindScar(data0,scarfull,scar,cln)
  print("Done")
}
PreDataBulk = function(filename,scarfull,scar,cln){
  data0 = read.table(filename,stringsAsFactors=F)
  # data0 = read.table(files["CB_UMI"],stringsAsFactors=F)
  names(data0) = c("Read.Seq","UMI")
  data0$Cell.BC = data0$UMI
  print("Find scar...")
  data = FindScar(data0,scarfull,scar,cln)
  print("Done")
}
PreData2 = function(filename,scarfull1,scarfull2,scar1,scar2,cln){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  if(!dir.exists("v1")){
    dir.create("v1")
  }
  if(!dir.exists("v2")){
    dir.create("v2")
  }
  print("Split and Find scar...")
  data = SplitScar(data0,scarfull1,scarfull2,scar1,scar2,cln)
  print("Done")
}
CallScar = function(INDEL,Scar,scarref_all,scarfull,scar,cln){
  print("Scar to INDEL Range...")
  data = INDELChangeForm(INDEL,Scar,scarref_all,".",cln)
  print("Done")
  print("Select Consensus reads...")
  INDELCons(INDEL,data,scarref,".",scarfull,scar,cln)
  print("Done")
}

ScarAnalysis = function(imagels,scarref,segref,consensus){
  for (i in 1:length(imagels)) {
    print(paste0("plot for ",names(imagels)[i],":"))
    
    print("INDEL trans...")
#    edit_num = EditStat(INDEL_ranges,read_counts)
    indel_r = INDELSplitDF(imagels[[i]])
    #plot
    print("INDEL Size plot...")
    tryCatch(INDELSizePlot(indel_r,out = paste0(names(imagels)[i],"_edit_size.pdf")), error=function(e) "plot error")
    
    print("INDEL Ribbon plot...")
    tryCatch(INDELRibbonPlot(indel_r$del_ranges_df,scarref,out = paste0(names(imagels)[i],"_edit_ribbon.pdf")), error=function(e) "plot error")
    pdf(paste0(names(imagels)[i],"_reads_proportion.pdf"),width=4,height=3)
    tryCatch(plot(density(consensus[,4+i]),main = paste0(names(imagels)[i]," scar density plot")), error=function(e) "plot error")
    dev.off()
    
    model_names = segref$V1
    pos_set = as.vector(t(segref[,2:3]))
    print("INDEL Chord plot...")
    tryCatch(DelChordPlot(indel_r$del_ranges_df,model_names,pos_set,
                 output = paste0(names(imagels)[i],"_chord.pdf")), error=function(e) "plot error")
  }
  
  pdf("umi_proportion.pdf",width=4,height=3)
  plot(density(consensus[,8]),main = "umi scar density plot")
  dev.off()
  vdata = list("consensus" = paste0(as.character(1:nrow(consensus)),consensus$cons),
               "reads main" = paste0(as.character(1:nrow(consensus)),consensus$main),
               "umi main" = paste0(as.character(1:nrow(consensus)),consensus$umim))
  pdf("venn.pdf",width=7,height=7)
  ggvenn(vdata)
  dev.off()
}

LineageAnalysis1 = function(Cells,scarref,segref,consensus){
  
  if(!dir.exists("LineageAnalysis")){
    dir.create("LineageAnalysis")
  }
  setwd("LineageAnalysis")
  
  prefix = c("cons","readm","umim")
  for (i in 1:3) {
    consensus$pattern = consensus[,i+1]
    hetdata = consensus[which(consensus$pattern != "unknown"),]
    data <- list(hetdata)
    tag = TagDataProcess(data, "", Cells)
    
    print("HeatMap plot...")
    if(!dir.exists(prefix[i])){
      dir.create(prefix[i])
    }
    setwd(prefix[i])
    
    if(!dir.exists("jaccard")){
      dir.create("jaccard")
    }
    setwd("jaccard/")
    tryCatch(CellMap1(tag,Cells), error=function(e) "plot error")
    setwd("../")
    
    print("tree plot...")
    if(!dir.exists("tagtree")){
      dir.create("tagtree")
    }
    setwd("tagtree/")
    tryCatch(TagTreeMain(data = tag,
                         prefix = pl,Cells = Cells,outname = "tag_tree"), error=function(e) "plot error")
    setwd("../")
    
    setwd("../")
  }
  
  print("DONE")
}

LineageAnalysis2 = function(Cells,scarref,segref,consensus1,consensus2){
  
  if(!dir.exists("LineageAnalysis")){
    dir.create("LineageAnalysis")
  }
  setwd("LineageAnalysis")
  
  prefix = c("cons","readm","umim")
  for (i in 1:3) {
    consensus1$pattern = consensus1[,i+1]
    hetdata1 = consensus1[which(consensus1$pattern != "unknown"),]
    consensus2$pattern = consensus2[,i+1]
    hetdata2 = consensus2[which(consensus2$pattern != "unknown"),]
    data <- list(hetdata1,hetdata2)
    pl = c("v1.","v2.")
    tag = TagDataProcess(data, pl, Cells)
    
    print("HeatMap plot...")
    if(!dir.exists(prefix[i])){
      dir.create(prefix[i])
    }
    setwd(prefix[i])
    
    if(!dir.exists("jaccard")){
      dir.create("jaccard")
    }
    setwd("jaccard/")
    tryCatch(CellMap1(tag,Cells), error=function(e) "plot error")
    setwd("../")
    
    print("tree plot...")
    if(!dir.exists("tagtree")){
      dir.create("tagtree")
    }
    setwd("tagtree/")
    tryCatch(TagTreeMain(data = tag,
                prefix = pl,Cells = Cells,outname = "tag_tree"), error=function(e) "plot error")
    setwd("../")
    
    setwd("../")
  }
  
  print("DONE")
}

LineageAnalysis1x = function(Cells,scarref,segref,consensus){
  
  if(!dir.exists("LineageAnalysis")){
    dir.create("LineageAnalysis")
  }
  setwd("LineageAnalysis")
  
  prefix = c("cons","readm","umim")
  for (i in 1:3) {
    consensus$pattern = consensus[,i+1]
    hetdata = consensus[which(consensus$pattern != "unknown"),]
    data <- list(hetdata)
    tag1 = ArrayDataProcess(data, "", Cells)
    tag2 = TagDataProcess(data, "", Cells)
    
    print("HeatMap plot...")
    if(!dir.exists(prefix[i])){
      dir.create(prefix[i])
    }
    setwd(prefix[i])
    
    if(!dir.exists("jaccard")){
      dir.create("jaccard")
    }
    setwd("jaccard/")
    tryCatch(CellMap1(tag1,Cells), error=function(e) "plot error")
    setwd("../")
    
    print("tree plot...")
    if(!dir.exists("tagtree")){
      dir.create("tagtree")
    }
    setwd("tagtree/")
    tryCatch(TagTreeMain(data = tag2,
                prefix = "",Cells = Cells,outname = "tag_tree"), error=function(e) "plot error")
    setwd("../")
    
    setwd("../")
  }
  
  print("DONE")
}

LineageAnalysis2x = function(Cells,scarref,segref,consensus1,consensus2){
  
  if(!dir.exists("LineageAnalysis")){
    dir.create("LineageAnalysis")
  }
  setwd("LineageAnalysis")
  
  prefix = c("cons","readm","umim")
  for (i in 1:3) {
    consensus1$pattern = consensus1[,i+1]
    hetdata1 = consensus1[which(consensus1$pattern != "unknown"),]
    consensus2$pattern = consensus2[,i+1]
    hetdata2 = consensus2[which(consensus2$pattern != "unknown"),]
    data <- list(hetdata1,hetdata2)
    pl = c("v1.","v2.")
    tag1 = ArrayDataProcess(data, pl, Cells)
    tag2 = TagDataProcess(data, pl, Cells)
    
    print("HeatMap plot...")
    if(!dir.exists(prefix[i])){
      dir.create(prefix[i])
    }
    setwd(prefix[i])
    
    if(!dir.exists("jaccard")){
      dir.create("jaccard")
    }
    setwd("jaccard/")
    tryCatch(CellMap1(tag1,Cells), error=function(e) "plot error")
    setwd("../")
    
    print("tree plot...")
    if(!dir.exists("tagtree")){
      dir.create("tagtree")
    }
    setwd("tagtree/")
    tryCatch(TagTreeMain(data = tag2,
                prefix = pl,Cells = Cells,outname = "tag_tree"), error=function(e) "plot error")
    setwd("../")
    
    setwd("../")
  }
  
  print("DONE")
}



#initialize
files = ReadSheet(sheetfile)
setwd(outpath)
cutsite = read.table(files["CutSite"])
segref = cutsite[-1,]
scarref = ReadCutsite(segref)
scarref_all<-ReadCutsite(cutsite,reftype="ALL")

scar = c("start" = cutsite[1,2], "end" = cutsite[1,3])


if(command == "PreData1"){
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  PreData1(files["CB_UMI"],scarfull,scar,cln)
  print("call scar:")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  CallScar(INDEL,Scar,scarref_all,scarfull,scar,cln)
  
}else if(command == "PreDataBulk"){
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  PreDataBulk(files["CB_UMI"],scarfull,scar,cln)
  print("call scar:")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  CallScar(INDEL,Scar,scarref_all,scarfull,scar,cln)
  
}else if(command == "PreData2"){
  cutsite1 = read.table(files["CutSite1"])
  segref1 = cutsite1[-1,]
  scarref1 = ReadCutsite(segref1)
  scarref1_all<-ReadCutsite(cutsite1,reftype="ALL")
  scar1 = c("start" = cutsite1[1,2], "end" = cutsite1[1,3])
  
  cutsite2 = read.table(files["CutSite2"])
  segref2 = cutsite2[-1,]
  scarref2 = ReadCutsite(segref2)
  scarref2_all<-ReadCutsite(cutsite2,reftype="ALL")
  scar2 = c("start" = cutsite2[1,2], "end" = cutsite2[1,3])
  
  sv1 = ReadFasta(files["Fasta1"])
  scarfull1 = sv1$scarfull
  sv2 = ReadFasta(files["Fasta2"])
  scarfull2 = sv2$scarfull
  PreData2(files["CB_UMI"],scarfull1,scarfull2,scar1,scar2,cln)
  
  print("call scar from v1:")
  setwd("v1")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta1"])
  scarfull = sv$scarfull1
  CallScar(INDEL,Scar,scarref1_all,scarfull1,scar1,cln)
  
  print("call scar from v2:")
  setwd("../v2")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta2"])
  scarfull = sv$scarfull1
  CallScar(INDEL,Scar,scarref2_all,scarfull2,scar2,cln)
  
}else if(command == "CallScar2"){
  cutsite1 = read.table(files["CutSite1"])
  segref1 = cutsite1[-1,]
  scarref1 = ReadCutsite(segref1)
  scarref1_all<-ReadCutsite(cutsite1,reftype="ALL")
  scar1 = c("start" = cutsite1[1,2], "end" = cutsite1[1,3])
  
  cutsite2 = read.table(files["CutSite2"])
  segref2 = cutsite2[-1,]
  scarref2 = ReadCutsite(segref2)
  scarref2_all<-ReadCutsite(cutsite2,reftype="ALL")
  scar2 = c("start" = cutsite2[1,2], "end" = cutsite2[1,3])
  
  sv1 = ReadFasta(files["Fasta1"])
  scarfull1 = sv1$scarfull
  sv2 = ReadFasta(files["Fasta2"])
  scarfull2 = sv2$scarfull
  
  print("call scar from v1:")
  setwd("v1")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta1"])
  scarfull = sv$scarfull1
  CallScar(INDEL,Scar,scarref1_all,scarfull1,scar1,cln)
  
  print("call scar from v2:")
  setwd("../v2")
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta2"])
  scarfull = sv$scarfull1
  CallScar(INDEL,Scar,scarref2_all,scarfull2,scar2,cln)
  
}else if(command == "CallScar"){
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T,stringsAsFactors = F,fill = T)
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  CallScar(INDEL,Scar,scarref_all,scarfull,scar,cln)
}else if(command == "ScarAnalysis"){
  image_cons = readRDS("image_cons.rds")
  image_main = readRDS("image_main.rds")
  image_umim = readRDS("image_umim.rds")
  imagels = list("consensus" = image_cons, "main reads" = image_main,"main umi" = image_umim)
  consensus = read.csv("final_scarform.csv",header = T,stringsAsFactors = F)
  ScarAnalysis(imagels,scarref,segref,consensus)
}else if(command == "LineageAnalysis1"){
  
  Cells<-read.csv(files["Cells"],header = T,stringsAsFactors = F)
  colnames(Cells) = c("Cell.BC","Cell.type")
  Cells$Cell.BC = substring(Cells$Cell.BC,1,16)
  consensus = read.csv("final_scarform.csv",header = T,stringsAsFactors = F)
  LineageAnalysis1(Cells,scarref,segref,consensus)
  
}else if(command == "LineageAnalysis2"){
  
  Cells<-read.csv(files["Cells"],header = T,stringsAsFactors = F)
  colnames(Cells) = c("Cell.BC","Cell.type")
  Cells$Cell.BC = substring(Cells$Cell.BC,1,16)
  consensus1 = read.csv("v1/final_scarform.csv",header = T,stringsAsFactors = F)
  consensus2 = read.csv("v2/final_scarform.csv",header = T,stringsAsFactors = F)
  LineageAnalysis2(Cells,scarref,segref,consensus1,consensus2)
  
}else if(command == "LineageAnalysis1x"){
  
  Cells<-read.csv(files["Cells"],header = T,stringsAsFactors = F)
  colnames(Cells) = c("Cell.BC","Cell.type")
  Cells$Cell.BC = substring(Cells$Cell.BC,1,16)
  consensus = read.csv("final_scarform.csv",header = T,stringsAsFactors = F)
  LineageAnalysis1x(Cells,scarref,segref,consensus)
  
}else if(command == "LineageAnalysis2x"){
  
  Cells<-read.csv(files["Cells"],header = T,stringsAsFactors = F)
  colnames(Cells) = c("Cell.BC","Cell.type")
  Cells$Cell.BC = substring(Cells$Cell.BC,1,16)
  consensus1 = read.csv("v1/final_scarform.csv",header = T,stringsAsFactors = F)
  consensus2 = read.csv("v2/final_scarform.csv",header = T,stringsAsFactors = F)
  LineageAnalysis2x(Cells,scarref,segref,consensus1,consensus2)
  
}


