source("1.scar_function.R")

#arg1 command arg2 sheetfile arg3 outpath
args=commandArgs(T)
command = args[1]
sheetfile = args[2]
outpath = args[3]
cln = as.numeric(args[4])

PreData1 = function(filename,scarfull,scar,cln){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  print("Find scar...")
  data = FindScar(data0,scarfull,scar,cln)
  print("Done")
}
PreDataBulk = function(filename,scarfull,scar,cln){
  data0 = read.table(filename,stringsAsFactors=F)
  names(data0) = c("Read.Seq","UMI")
  data0$Cell.BC = data0$UMI
  print("Find scar...")
  data = FindScar(data0,scarfull,scar,cln)
  print("Done")
}
PreData2 = function(filename,scarfull1,scar1,scarfull2,scar2,cln){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  if(!dir.exists("v1")){
    dir.create("v1")
  }
  if(!dir.exists("v2")){
    dir.create("v2")
  }
  print("Split and Find scar...")
  data = SplitScar(data0,scarfull1,scar1,scarfull2,scar2,cln)
  print("Done")
}
CallScar = function(INDEL,Scar,segref,cln){
  print("Scar to INDEL Range...")
  data = INDELChangeForm(INDEL,Scar,segref,".",cln)
  print("Done")
  print("Select Consensus reads...")
  INDELCons(INDEL,data,segref,".",cln)
  print("Done")
}

ScarAnalysis = function(image,scarseg,segref,read_counts){
  INDEL_ranges = readRDS(image)
  print("INDEL trans...")
  #edit_num = EditStat(INDEL_ranges,read_counts)
  indel_r = INDELSplitDF(INDEL_ranges)
  #plot
  print("INDEL Size plot...")
  INDELSizePlot(indel_r,out = "edit_size.pdf")
  print("INDEL Ribbon plot...")
  INDELRibbonPlot(indel_r$del_ranges_df,scarseg,out = "edit_ribbon.pdf")
  pdf("reads_proportion.pdf",width=4,height=3)
  plot(density(read_counts$UMI_pro),main = "scar density plot")
  dev.off()
  
  model_names = segref$V1
  pos_set = as.vector(t(segref[,2:3]))
  print("INDEL Chord plot...")
  DelChordPlot(indel_r$del_ranges_df,model_names,pos_set,
               output = "chord.pdf")
  print("DONE")
  
}


#initialize
files = ReadSheet(sheetfile)
setwd(outpath)


if(command == "PreData1"){
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  scar = sv$scar
  PreData1(files["CB_UMI"],scarfull,scar,cln)
}else if(command == "PreDataBulk"){
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  scar = sv$scar
  PreDataBulk(files["CB_UMI"],scarfull,scar,cln)
}else if(command == "PreData2"){
  sv1 = ReadFasta(files["Fasta1"])
  scarfull1 = sv1$scarfull
  scar1 = sv1$scar
  sv2 = ReadFasta(files["Fasta2"])
  scarfull2 = sv2$scarfull
  scar2 = sv2$scar
  PreData2(files["CB_UMI"],scarfull1,scar1,scarfull2,scar2,cln)
}else if(command == "CallScar"){
  segrf = ReadCutsite(files["CutSite"])
  INDEL = readRDS("indel.rds")
  Scar = read.table("UMI_reads_scar_full.txt",header = T)
  CallScar(INDEL,Scar,segref,cln)
}else if(command == "ScarAnalysis"){
  image = readRDS("image.rds")
  scarseg = ReadCutsite(files["CutSite"])
  read_counts = read.csv("consensus.csv",header = T)
  segref = read.table(files["CutSite"])
  ScarAnalysis(image,scarseg,segref,read_counts)
}


