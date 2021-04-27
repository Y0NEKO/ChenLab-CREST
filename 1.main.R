source("D:/LAB/p6_chen_lab/myscript/newest_4_23/1.scar_function.R")

#arg1 command arg2 sheetfile arg3 outpath
args=commandArgs(T)
command = args[1]
sheetfile = args[2]
outpath = args[3]


PreData1 = function(filename,scarfull,scar,segref){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  print("Find scar...")
  data = FindScar(data0,scarfull,scar)
  print("Done")
  print("Scar to INDEL Range...")
  INDEL_ranges = INDELChange(data$INDEL,data$Scar,segref,".")
  print("Done")
}
PreDataBulk = function(filename,scarfull,scar,segref){
  data0 = read.table(filename,stringsAsFactors=F)
  names(data0) = c("Read.Seq","UMI")
  data0$Cell.BC = data0$UMI
  print("Find scar...")
  data = FindScar(data0,scarfull,scar)
  print("Done")
  print("Scar to INDEL Range...")
  INDEL_ranges = INDELChange(data$INDEL,data$Scar,segref,".")
  print("Done")
}

PreData2 = function(filename,scarfull1,scar1,scarfull2,scar2,segref1,segref2){
  data0 = read.table(filename, stringsAsFactors=F, header=T)
  if(!dir.exists("v1")){
    dir.create("v1")
  }
  if(!dir.exists("v2")){
    dir.create("v2")
  }
  print("Split and Find scar...")
  data = SplitScar(data0,scarfull1,scar1,scarfull2,scar2)
  print("Done")
  print("Change Scar form for V1...")
  INDEL_ranges1 = INDELChange(data$V1_INDEL,data$V1_Scar,segref1,"v1")
  print("Done")
  print("Change Scar form for V2...")
  INDEL_ranges2 = INDELChange(data$V2_INDEL,data$V2_Scar,segref2,"v2")
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
  segref = ReadCutsite(files["CutSite"])
  PreData1(files["CB_UMI"],scarfull,scar,segref)
}else if(command == "PreDataBulk"){
  sv = ReadFasta(files["Fasta"])
  scarfull = sv$scarfull
  scar = sv$scar
  segref = ReadCutsite(files["CutSite"])
  PreDataBulk(files["CB_UMI"],scarfull,scar,segref)
}else if(command == "PreData2"){
  sv1 = ReadFasta(files["Fasta1"])
  scarfull1 = sv1$scarfull
  scar1 = sv1$scar
  sv2 = ReadFasta(files["Fasta2"])
  scarfull2 = sv2$scarfull
  scar2 = sv2$scar
  segref1 = ReadCutsite(files["CutSite1"])
  segref2 = ReadCutsite(files["CutSite2"])
  PreData2(files["CB_UMI"],scarfull1,scar1,scarfull2,scar2,segref1,segref2)
}else if(command == "ScarAnalysis"){
  image = files["INDEL"]
  scarseg = ReadCutsite(files["CutSite"])
  read_counts = read.csv(files["Cons"],header = T)
  segref = read.table(files["CutSite"])
  ScarAnalysis(image,scarseg,segref,read_counts)
}



