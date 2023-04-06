#CREST analysis script 2022.12.01 by Y0NEKO
#qc with input 'final_scarform.csv' file and build data for next analysis
setwd("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/")
.libPaths(c("/home/liuhengxin/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library",
            "/picb/sysgenomics2/people/liuhengxin/software/anaconda3/envs/r4-base/lib/R/library/",
            "/picb/sysgenomics2/people/liuhengxin/software/miniconda3/envs/r4.1/lib/R/library/"))
source("myscript/anascipt/tree_method.R")


#load all sample
{
  #bulk RNA
  samplels = list.files(c("bulk_result/X.bulk_in_vivo/nkx_p1-1/",
                          "bulk_result/X.bulk_in_vivo/nkx_p1-2/",
                          "bulk_result/X.bulk_in_vivo/nkx_p1-3/"),
                        full.names = T)
  v1lsbnk = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsbnk = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsbnk) = names(v2lsbnk) = list.files(c("bulk_result/X.bulk_in_vivo/nkx_p1-1/",
                                                 "bulk_result/X.bulk_in_vivo/nkx_p1-2/",
                                                 "bulk_result/X.bulk_in_vivo/nkx_p1-3/"),
                                               full.names = T)
  nkx1v1 = v1lsbnk[1:8]
  nkx2v1 = v1lsbnk[9:16]
  nkx3v1 = v1lsbnk[17:24]
  nkx1v2 = v2lsbnk[1:8]
  nkx2v2 = v2lsbnk[9:16]
  nkx3v2 = v2lsbnk[17:24]
  MyCombine = function(nkx1v1){
    nkx1v1t = NULL
    for (i in 1:length(nkx1v1)) {
      nkx1v1t = rbind(nkx1v1t,nkx1v1[[i]])
    }
    return(nkx1v1t)
  }
  nkx1v1t = MyCombine(nkx1v1)
  nkx2v1t = MyCombine(nkx2v1)
  nkx3v1t = MyCombine(nkx3v1)
  nkx1v2t = MyCombine(nkx1v2)
  nkx2v2t = MyCombine(nkx2v2)
  nkx3v2t = MyCombine(nkx3v2)
  
  samplels = list.files("bulk_result/X.bulk_in_vivo/",full.names = T)
  samplels = samplels[grep("Dual",samplels)]
  v1lsb = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsb = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsb) = names(v2lsb) = unlist(lapply(strsplit(samplels,"/"), "[[",4)) 
  v1lsb$Dual_E15.5_1 = rbind(v1lsb$Dual_E15.5_1a,v1lsb$Dual_E15.5_1b)
  v1lsb$Dual_E15.5_2 = rbind(v1lsb$Dual_E15.5_2a,v1lsb$Dual_E15.5_2b)
  v1lsb$Dual_E15.5_3 = rbind(v1lsb$Dual_E15.5_3a,v1lsb$Dual_E15.5_3a_less)
  v1lsb = v1lsb[which(!names(v1lsb) %in% c("Dual_E15.5_1a","Dual_E15.5_1b",
                                           "Dual_E15.5_2a","Dual_E15.5_2b",
                                           "Dual_E15.5_3a","Dual_E15.5_3a_less"))]
  v2lsb$Dual_E15.5_1 = rbind(v2lsb$Dual_E15.5_1a,v2lsb$Dual_E15.5_1b)
  v2lsb$Dual_E15.5_2 = rbind(v2lsb$Dual_E15.5_2a,v2lsb$Dual_E15.5_2b)
  v2lsb$Dual_E15.5_3 = rbind(v2lsb$Dual_E15.5_3a,v2lsb$Dual_E15.5_3a_less)
  v2lsb = v2lsb[which(!names(v2lsb) %in% c("Dual_E15.5_1a","Dual_E15.5_1b",
                                           "Dual_E15.5_2a","Dual_E15.5_2b",
                                           "Dual_E15.5_3a","Dual_E15.5_3a_less"))]
  
  v1lsb$nkx1 = nkx1v1t;v2lsb$nkx1 = nkx1v2t
  v1lsb$nkx2 = nkx2v1t;v2lsb$nkx2 = nkx2v2t
  v1lsb$nkx3 = nkx3v1t;v2lsb$nkx3 = nkx3v2t
  for (i in 1:length(v1lsb)) {
    v1lsb[[i]]$pattern = v1lsb[[i]]$main
    v2lsb[[i]]$pattern = v2lsb[[i]]$main
  }
  qsave(list("v1" = v1lsb,"v2" = v2lsb),"bulk_result/X.bulk_in_vivo/bulk_total.qs")
  bulkraw = qread("bulk_result/X.bulk_in_vivo/bulk_total.qs")
  v1lsb = bulkraw$v1
  v2lsb = bulkraw$v2
  
  #scRNA
  samplels =list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE15/pandaE15_rep1/",full.names = T)[-1]
  ScLoad = function(filepath,cellpath){
    samplels = filepath
    v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
    v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
    names(v1lsx) = names(v2lsx) = unlist(lapply(strsplit(filepath,"/"),"[[",length(strsplit(filepath,"/")[[1]])))
    cell = read.csv(cellpath)
    return(list("v1lsx" = v1lsx,"v2lsx" = v2lsx, "cell" = cell))
  }
  E15raw1 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE15/pandaE15_rep1/",full.names = T)[-1],"sc_result/CREST_final_20221204/E15_identityV1_20221204.csv")
  E15raw2 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE15/pandaE15_rep2/",full.names = T)[-1],"sc_result/CREST_final_20221204/E15_identityV1_20221204.csv")
  E15raw3 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE15/pandaE15_rep3/",full.names = T)[-1],"sc_result/CREST_final_20221204/E15_identityV1_20221204.csv")
  
  E11raw1 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE11/rep1/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  E11raw2 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE11/rep2/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  E11raw3 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/snapE11/rep3/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")

  E11praw1 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaE11/rep1/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  E11praw2 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaE11/rep2/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  E11praw3 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaE11/rep3/",full.names = T),"sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  
  OGNraw1 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaOGN/rep1/",full.names = T),"sc_result/CREST_final_20221204/Organoid_identity_20221204.csv")
  OGNraw2 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaOGN/rep2/",full.names = T),"sc_result/CREST_final_20221204/Organoid_identity_20221204.csv")
  OGNraw3 = ScLoad(list.files("sc_result/CREST_final_20221204/CREST RECOVERY/pandaOGN/rep3/",full.names = T),"sc_result/CREST_final_20221204/Organoid_identity_20221204.csv")
  
  e15cell = read.csv("sc_result/CREST_final_20221204/E15_all_cellident_group_final_20221207.csv")
  e11cell = read.csv("sc_result/CREST_final_20221204/E11_identity_20221204.csv")
  ogncell = read.csv("sc_result/CREST_final_20221204/pandaOGN_identity_group_20230108.csv")
  
  # cellls = list("E15_1" = e15cell1,"E15_2" = e15cell2, "E11" = e11cell,"OGN" = ogncell)
  cellls = list("E15" = e15cell,"E11" = e11cell,"OGN" = ogncell)
  for (i in 1:length(cellls)) {
    cellls[[i]]$Cell.BC = substr(cellls[[i]]$X,1,16)
    colnames(cellls[[i]])[2] = "Cell.type"
  }
  cellls$E15$Cell.BC = cellls$E15$X
  cellls$E11_rep1 = cellls$E11[which(cellls$E11$group == "snapE11.5-rep1"),]
  cellls$E11_rep2 = cellls$E11[which(cellls$E11$group == "snapE11.5-rep2"),]
  cellls$E11_rep3 = cellls$E11[which(cellls$E11$group == "snapE11.5-rep3"),]
  cellls$pandaE11_rep1 = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep1"),]
  cellls$pandaE11_rep2 = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep2"),]
  cellls$pandaE11_rep3 = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep3"),]
  cellls$E15_rep1 = cellls$E15[which(cellls$E15$group == "snapE15.5-rep1"),]
  cellls$E15_rep2 = cellls$E15[which(cellls$E15$group == "snapE15.5-rep2"),]
  cellls$E15_rep3 = cellls$E15[which(cellls$E15$group == "snapE15.5-rep3"),]
  cellls$OGN_rep1 = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep1"),]
  cellls$OGN_rep2 = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep2"),]
  cellls$OGN_rep3 = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep3"),]
  # cellls$E15_1$Cell.BC = cellls$E15_1$X
  # cellls$E15_2$Cell.BC = cellls$E15_2$X
  
  E15raw1$cell = cellls$E15[which(cellls$E15$group == "snapE15.5-rep1"),]
  E15raw2$cell = cellls$E15[which(cellls$E15$group == "snapE15.5-rep2"),]
  E15raw3$cell = cellls$E15[which(cellls$E15$group == "snapE15.5-rep3"),]

  # E15raw1$cell2 = cellls$E15_2[which(cellls$E15_2$group == "snapE15.5-rep1"),]
  # E15raw2$cell2 = cellls$E15_2[which(cellls$E15_2$group == "snapE15.5-rep2"),]
  # E15raw3$cell2 = cellls$E15_2[which(cellls$E15_2$group == "snapE15.5-rep3"),]
  
  E11raw1$cell = cellls$E11[which(cellls$E11$group == "snapE11.5-rep1"),]
  E11raw2$cell = cellls$E11[which(cellls$E11$group == "snapE11.5-rep2"),]
  E11raw3$cell = cellls$E11[which(cellls$E11$group == "snapE11.5-rep3"),]
  
  E11praw1$cell = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep1"),]
  E11praw2$cell = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep2"),]
  E11praw3$cell = cellls$E11[which(cellls$E11$group == "pandaE11.5-rep3"),]
  
  OGNraw1$cell = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep1"),]
  OGNraw2$cell = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep2"),]
  OGNraw3$cell = cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep3"),]
  
  # intersect(OGNraw2$v1lsx$`Scarlib-1`$Cell.BC,cellls$E11[which(cellls$E11$group == "pandaE11.5-rep3"),"Cell.BC"])
  # intersect(E11praw3$v1lsx$`Scarlib-1`$Cell.BC,cellls$OGN[which(cellls$OGN$group == "pandaOGN-rep2"),"Cell.BC"])
  
  scrawls = list("E11_rep1" = E11raw1,"E11_rep2" = E11raw2,
                 "E11_rep3" = E11raw3,
                 "pandaE11_rep1" = E11praw1,"pandaE11_rep2" = E11praw2,"pandaE11_rep3" = E11praw3,
                 "E15_rep1" = E15raw1,"E15_rep2" = E15raw2,"E15_rep3" = E15raw3,
                 "OGN_rep1" = OGNraw1,"OGN_rep2" = OGNraw2,"OGN_rep3" = OGNraw3)
  library(qs)
  qsave(scrawls,"sc_result/CREST_final_20221204/scrawls_01_08.qs")
  scrawls =qread("sc_result/CREST_final_20221204/scrawls_12_14.qs")
  # v1lsc = 
  
  n = 12
  myrbind = function(x,y){
    return(rbind(x[c("Cell.BC","umim")],y[c("Cell.BC","umim")]))
  }
  for (i in 1:length(scrawls)) {
    v1lsl = Reduce(myrbind,scrawls[[i]]$v1lsx)
    v2lsl = Reduce(myrbind,scrawls[[i]]$v2lsx)
    v1lsb[[n+1]] = v1lsl
    v2lsb[[n+1]] = v2lsl
    names(v1lsb)[n+1] = names(v2lsb)[n+1] = names(scrawls)[i]
    v1lsb[[n+1]]$pattern = v1lsb[[n+1]]$umim
    v2lsb[[n+1]]$pattern = v2lsb[[n+1]]$umim
    n = n +1
  }
  
  qsave(list("v1" = v1lsb,"v2" = v2lsb),"sc_result/CREST_final_20221204/allrawls_12_05.qs")
  
}

#blacklist select
{
  MergeAllArray = function(v1lsb){
    arraydft = NULL
    for (i in 2:length(v1lsb)) {
      arraydf = as.data.frame(table(v1lsb[[i]]$pattern))
      arraydf = arraydf[which(!arraydf$Var1 %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      
      arraydf = arraydf[order(-arraydf$Freq),]
      arraydf$Freq = arraydf$Freq/sum(arraydf$Freq)
      
      colnames(arraydf)[2] = names(v1lsb)[i]
      if(i == 2){
        arraydft = arraydf
      }else{
        arraydft = merge(arraydft, arraydf, by = "Var1", all = T)
      }
    }
    
    arraydft[is.na(arraydft)] = 0
    # blackls = arraydft[which(rowSums(arraydft[-1]>0.001)>1),]
    blackls = arraydft
    return(blackls)
    
  }
  MergeAllIndel = function(v1lsb){
    tagdft = NULL
    for (i in 2:length(v1lsb)) {
      print(i)
      df = v1lsb[[i]]$pattern
      tag = unlist(lapply(strsplit(df,"_|&"),unique))
      tag = tag[which(!tag == "NONE")]
      
      tagdf = as.data.frame(table(tag))
      tagdf[,2] = tagdf[,2]/sum(tagdf[,2])
      
      colnames(tagdf)[2] = names(v1lsb)[i]
      if(i == 2){
        tagdft = tagdf
      }else{
        tagdft = merge(tagdft, tagdf, by = "tag", all = T)
      }
    }
    
    tagdft[is.na(tagdft)] = 0
    # blackls = tagdft[which(rowSums(tagdft[-1]>0.005)>1),]
    blackls = tagdft
    return(blackls)
  }
  
  v1blayb = MergeAllArray(v1lsb)
  v2blayb = MergeAllArray(v2lsb)
  # v1bltgb = MergeAllIndel(v1lsb)
  # v2bltgb = MergeAllIndel(v2lsb)
  
  
  blt = list("v1_array" = v1blayb,"v2_array" = v2blayb)
  qsave(blt,file = "sc_result/CREST_final_20221204/total_sample_stat_12_05.qs")
  blt = qread(file = "sc_result/CREST_final_20221204/total_sample_stat_12_05.qs")
  # blt = readRDS("final_result/x.blacklist_total_2_24.rds")
  tmp2 = blt$v2_array
  bl = list()
  for (i in 1:2) {
    line = blt[[i]]
    # bl[[i]] = line[which(rowSums(line[13:15]>0.001)>1 | rowSums(line[10:12]>0.005)>1 | rowSums(line[2:9]>0.001)>1),]
    # bl[[i]] = line[which(rowSums(line[2:17]>0.001)>1 | rowSums(line[2:ncol(line)]>0)==(ncol(line)-1)),]
    bl[[i]] = line[which(rowSums(line[c(2:14,16:21)]>0.01)>1 | rowSums(line[c(2:14,18:24)]>0.01)>1 | 
                           rowSums(line[2:ncol(line)]>0)==(ncol(line)-1) |
                           rowSums(line[2:ncol(line)]>0.3)>0),]
    # bl[[i]] = line[which(rowSums(line[2:ncol(line)]>0.1)>0),]
    names(bl)[i] = names(blt)[i]
  }
  
  nrow(bl$v1_array)
  nrow(bl$v2_array)
  nrow(bl$v2_tag)
  nrow(bl$v1_tag)
  
  dir.create("final_result/blacklist")
  qsave(bl,"final_result/blacklist/x.blacklist_0.01_12_14_final_all_sample.qs")
  # qsave(bl,"final_result/blacklist/x.blacklist_0.01_12_09_bulk_all_sample.qs")
  # qsave(bl,"final_result/blacklist/x.blacklist_0.001_12_06_bulk_e11_sample.qs")
  bl = qread("final_result/blacklist/x.blacklist_0.01_12_14_final_all_sample.qs")
  # bl = readRDS("final_result/x.blacklist_0.001_2_25.rds")
  
  
}


#Filter data
{
  
  FilterRecover = function(v1lsx,cell,thres,title){
    library(ggvenn)
    library(ggplot2)
    library(ggpubr)
    rawcbls = list()
    rawfcbls = list()
    cmfcbls = list()
    
    v1lsxf = list()
    #hard filter: umi counts > 3
    for (i in 1:length(v1lsx)) {
      rawcbls[[i]] = v1lsx[[i]]$Cell.BC
      v1lsxf[[i]] = v1lsx[[i]][which(v1lsx[[i]]$umi_num > thres),]
      rawfcbls[[i]] = v1lsxf[[i]]$Cell.BC
    }
    comBCraw = Reduce(intersect,rawcbls)
    comBC = Reduce(intersect,rawfcbls)
    comBCwithC = comBC[which(comBC %in% cell$Cell.BC)]
    
    
    for (i in 1:length(v1lsxf)) {
      v1lsxf[[i]] = v1lsxf[[i]][which(v1lsxf[[i]]$Cell.BC %in% comBC),]
      cmfcbls[[i]] = paste0(v1lsxf[[i]]$Cell.BC,",",v1lsxf[[i]]$umim)
    }
    names(rawcbls) = names(rawfcbls) =  names(cmfcbls) = names(v1lsxf) = names(v1lsx)
    
    comarray1 = Reduce(intersect,cmfcbls)
    comarray1withC = comarray1[which(unlist(lapply(strsplit(comarray1,","), "[[",1)) %in% comBCwithC)]
    result = v1lsx[[1]]
    result$bcarray = paste0(result$Cell.BC,",",result$umim)
    result = result[which(result$bcarray %in% comarray1),]
    # resultst = result %>% group_by(umim) %>% summarise(prop = n()/nrow(result))
    # resultf = result[which(result$umim %in% resultst[which(resultst$prop < 0.3),]$umim),]
    result_withcellan = merge(result[,c("Cell.BC","umim")],cell[c("Cell.BC","Cell.type")],by="Cell.BC");colnames(result)[2] = "pattern";colnames(result_withcellan)[2] = "pattern"
    stat = data.frame("process" = c("common_BCraw","common_BC","common_BCwithCell",
                                    "common_array","common_arraywithCell","common_with_cell"), 
                      "count" = c(length(comBCraw),length(comBC),length(comBCwithC),
                                  length(comarray1),length(comarray1withC),
                                  nrow(result_withcellan)))
    
    p1.0 = ggvenn(rawcbls,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[1])
    p1.1 = ggvenn(rawfcbls,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[2])
    p1.2 = ggvenn(cmfcbls,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[3])
    p = ggarrange(p1.0,p1.1,p1.2,nrow = 3)
    
    return(list("resultraw" = result,"result_withcellan" = result_withcellan,
                "rawcblist" = rawcbls,"rawfcblist" = rawfcbls,
                "cmfcblist" = cmfcbls,"stat" = stat,"venn_plot" = p))
  }
  # scrawls$E15_rep1$cell = scrawls$E15_rep1$cell1
  # scrawls$E15_rep2$cell = scrawls$E15_rep2$cell1
  # scrawls$E15_rep3$cell = scrawls$E15_rep3$cell1
  
  datalsv1 = list()
  datalsv2 = list()
  {
    for (i in 1:length(scrawls)) {
      v1title = c("V1 raw BC recover rate","V1 UMIs num(3) filtered BC recover rate",
                  "V1 UMIs num(3) filtered BC-pattern recover rate")
      v2title = c("V2 raw BC recover rate","V2 UMIs num(3) filtered BC recover rate",
                  "V2 UMIs num(3) filtered BC-pattern recover rate")
      datalsv1[[i]] = FilterRecover(scrawls[[i]]$v1lsx,scrawls[[i]]$cell,3,v1title)
      datalsv2[[i]] = FilterRecover(scrawls[[i]]$v2lsx,scrawls[[i]]$cell,3,v2title)
      names(datalsv1)[i] = names(datalsv2)[i] = names(scrawls)[i]
    }
    # datalsv1$E15_rep1_v2 = FilterRecover(scrawls$E15_rep1$v1lsx,scrawls$E15_rep1$cell2,3,v1title)
    # datalsv2$E15_rep1_v2 = FilterRecover(scrawls$E15_rep1$v2lsx,scrawls$E15_rep1$cell2,3,v2title)
    # datalsv1$E15_rep2_v2 = FilterRecover(scrawls$E15_rep2$v1lsx,scrawls$E15_rep2$cell2,3,v1title)
    # datalsv2$E15_rep2_v2 = FilterRecover(scrawls$E15_rep2$v2lsx,scrawls$E15_rep2$cell2,3,v2title)
    # datalsv1$E15_rep3_v2 = FilterRecover(scrawls$E15_rep3$v1lsx,scrawls$E15_rep3$cell2,3,v1title)
    # datalsv2$E15_rep3_v2 = FilterRecover(scrawls$E15_rep3$v2lsx,scrawls$E15_rep3$cell2,3,v2title)
    
  }
  
  for (i in 1:length(datalsv1)) {
    print(names(datalsv1)[i])
    print(datalsv1[[i]]$stat)
    print(datalsv2[[i]]$stat)
  }
  datals$v1 = datalsv1
  datals$v2 = datalsv2
  # datals = list("v1" = datalsv1,"v2" = datalsv2)
  datals$cellls = cellls
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_01_08_raw.qs")
  datals = qread(file = "final_result/CREST_final_analysis/array_recover_data_total_12_14.qs")
  
  #filter with blacklist
  bl = datals$blacklist
  # bl = qread("final_result/blacklist/x.blacklist_0.01_12")
  {
    datals$arraydata = list()
    datals$arraydata_bl = list()
    for (i in 1:length(datals$v1)) {
      colnames(datals$v1[[i]]$result_withcellan)[2] = "pattern"
      v1f = datals$v1[[i]]
      # result_withcellan_blf_v1 = v1f$result_withcellan[which(!v1f$result_withcellan$pattern %in% bl$v1_array$Var1),]
      
      colnames(datals$v2[[i]]$result_withcellan)[2] = "pattern"
      v2f = datals$v2[[i]]
      # result_withcellan_blf_v2 = v2f$result_withcellan[which(!v2f$result_withcellan$pattern %in% bl$v2_array$Var1),]
      
      mergedf = merge(v1f$result_withcellan, v2f$result_withcellan, by = "Cell.BC",all = T, suffixes = c(".v1",".v2"))
      mergedf$Cell.type = mergedf$Cell.type.v1
      mergedf[which(is.na(mergedf$Cell.type)),"Cell.type"] = mergedf[which(is.na(mergedf$Cell.type)),"Cell.type.v2"]
      mergedf = mergedf[-c(3,5)]
      mergedf$pattern = paste0(mergedf$pattern.v1,"-",mergedf$pattern.v2)
      
      mergedff = mergedf[which((!mergedf$pattern.v1 %in% bl$v1_array$Var1 &
                                  !is.na(mergedf$pattern.v1)) | 
                                 (!mergedf$pattern.v2 %in% bl$v2_array$Var1) &
                                 !is.na(mergedf$pattern.v2)),]
      mergedff$pattern = paste0(mergedff$pattern.v1,"-",mergedff$pattern.v2)
      
      tmp = mergedff %>% group_by(pattern) %>% summarise(prop = n()/nrow(mergedff))
      # mergedff = mergedff[which(!mergedff$pattern %in% tmp[which(tmp$prop > 0.1),]$pattern),]
      
      print(names(datals$v1)[i])
      print(max(tmp$prop))
      # datals$v2[[i]]$stat = rbind(datals$v2[[i]]$stat,data.frame("process" = "blacklist_filter","count" = nrow(result_withcellan_blf_v2)))
      # datals$v1[[i]]$stat = rbind(datals$v1[[i]]$stat,data.frame("process" = "blacklist_filter","count" = nrow(result_withcellan_blf_v1)))
      datals$arraydata[[i]] = mergedf
      datals$arraydata_bl[[i]] = mergedff
      names(datals$arraydata)[i] = names(datals$arraydata_bl)[i] = names(datals$v1)[i]
    }
    for (i in 1:length(datals$arraydata_bl)) {
      datals$arraydata_bl[[i]]$sample = names(datals$arraydata_bl)[i]
      datals$arraydata_bl[[i]]$pattern = paste0(names(datals$arraydata_bl)[i],"-",
                                                 datals$arraydata_bl[[i]]$pattern)
      
      datals$arraydata[[i]]$sample = names(datals$arraydata)[i]
      datals$arraydata[[i]]$pattern = paste0(names(datals$arraydata)[i],"-",
                                                 datals$arraydata[[i]]$pattern)
    }
    datals$arraydata_bl$E11 = rbind(datals$arraydata_bl$E11_rep1,datals$arraydata_bl$E11_rep2,datals$arraydata_bl$E11_rep3)
    datals$arraydata_bl$E15 = rbind(datals$arraydata_bl$E15_rep1,datals$arraydata_bl$E15_rep2,datals$arraydata_bl$E15_rep3)
    # datals$arraydata_bl$E15v2 = rbind(datals$arraydata_bl$E15_rep1_v2,datals$arraydata_bl$E15_rep2_v2,datals$arraydata_bl$E15_rep3_v2)
    datals$arraydata_bl$pandaE11 = rbind(datals$arraydata_bl$pandaE11_rep1,datals$arraydata_bl$pandaE11_rep2,datals$arraydata_bl$pandaE11_rep3)
    datals$arraydata_bl$OGN = rbind(datals$arraydata_bl$OGN_rep1,datals$arraydata_bl$OGN_rep2,datals$arraydata_bl$OGN_rep3)
    
    
    datals$arraydata$E11 = rbind(datals$arraydata$E11_rep1,datals$arraydata$E11_rep2,datals$arraydata$E11_rep3)
    datals$arraydata$E15 = rbind(datals$arraydata$E15_rep1,datals$arraydata$E15_rep2,datals$arraydata$E15_rep3)
    # datals$arraydata$E15v2 = rbind(datals$arraydata$E15_rep1_v2,datals$arraydata$E15_rep2_v2,datals$arraydata$E15_rep3_v2)
    datals$arraydata$pandaE11 = rbind(datals$arraydata$pandaE11_rep1,datals$arraydata$pandaE11_rep2,datals$arraydata$pandaE11_rep3)
    datals$arraydata$OGN = rbind(datals$arraydata$OGN_rep1,datals$arraydata$OGN_rep2,datals$arraydata$OGN_rep3)
    
    datals$v1$E15_rep2$resultraw = datals$v1$E15_rep2$resultraw[-c(3,5)]
    datals$v2$E15_rep2$resultraw = datals$v2$E15_rep2$resultraw[-c(3,5)]
    datals$v1$E11$resultraw = rbind(datals$v1$E11_rep1$resultraw,datals$v1$E11_rep2$resultraw,datals$v1$E11_rep3$resultraw)
    datals$v1$pandaE11$resultraw = rbind(datals$v1$pandaE11_rep1$resultraw,datals$v1$pandaE11_rep2$resultraw,datals$v1$pandaE11_rep3$resultraw)
    datals$v1$E15$resultraw = rbind(datals$v1$E15_rep1$resultraw,datals$v1$E15_rep2$resultraw,datals$v1$E15_rep3$resultraw)
    datals$v1$OGN$resultraw = rbind(datals$v1$OGN_rep1$resultraw,datals$v1$OGN_rep2$resultraw,datals$v1$OGN_rep3$resultraw)
    
    datals$v2$E11$resultraw = rbind(datals$v2$E11_rep1$resultraw,datals$v2$E11_rep2$resultraw,datals$v2$E11_rep3$resultraw)
    datals$v2$pandaE11$resultraw = rbind(datals$v2$pandaE11_rep1$resultraw,datals$v2$pandaE11_rep2$resultraw,datals$v2$pandaE11_rep3$resultraw)
    datals$v2$E15$resultraw = rbind(datals$v2$E15_rep1$resultraw,datals$v2$E15_rep2$resultraw,datals$v2$E15_rep3$resultraw)
    datals$v2$OGN$resultraw = rbind(datals$v2$OGN_rep1$resultraw,datals$v2$OGN_rep2$resultraw,datals$v2$OGN_rep3$resultraw)
    
    datals$v1$E11$result_withcellan = rbind(datals$v1$E11_rep1$result_withcellan,datals$v1$E11_rep2$result_withcellan,datals$v1$E11_rep3$result_withcellan)
    datals$v1$pandaE11$result_withcellan = rbind(datals$v1$pandaE11_rep1$result_withcellan,datals$v1$pandaE11_rep2$result_withcellan,datals$v1$pandaE11_rep3$result_withcellan)
    datals$v1$E15$result_withcellan = rbind(datals$v1$E15_rep1$result_withcellan,datals$v1$E15_rep2$result_withcellan,datals$v1$E15_rep3$result_withcellan)
    datals$v1$OGN$result_withcellan = rbind(datals$v1$OGN_rep1$result_withcellan,datals$v1$OGN_rep2$result_withcellan,datals$v1$OGN_rep3$result_withcellan)
    
    datals$v2$E11$result_withcellan = rbind(datals$v2$E11_rep1$result_withcellan,datals$v2$E11_rep2$result_withcellan,datals$v2$E11_rep3$result_withcellan)
    datals$v2$pandaE11$result_withcellan = rbind(datals$v2$pandaE11_rep1$result_withcellan,datals$v2$pandaE11_rep2$result_withcellan,datals$v2$pandaE11_rep3$result_withcellan)
    datals$v2$E15$result_withcellan = rbind(datals$v2$E15_rep1$result_withcellan,datals$v2$E15_rep2$result_withcellan,datals$v2$E15_rep3$result_withcellan)
    datals$v2$OGN$result_withcellan = rbind(datals$v2$OGN_rep1$result_withcellan,datals$v2$OGN_rep2$result_withcellan,datals$v2$OGN_rep3$result_withcellan)
    
  }
  
  lapply(datals$arraydata, nrow)
  lapply(datals$arraydata_bl, nrow)
  lapply(datals$v1, function(x) print(x$stat))
  lapply(datals$v2, function(x) print(x$stat))
  
  # dir.create("final_result/CREST_final_analysis")
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_12_14_blallsample0.01_0.3single.qs")
  
  
}

#data build
{
  datals$figure = list()
  #base stat, heatmap and lint tree construct
  #qc plot
  statdf = NULL
  for (i in 1:12) {
    datals$v1[[i]]$stat$sample = names(datals$v1)[i]
    datals$v1[[i]]$stat$group = "v1"
    statdf = rbind(statdf,datals$v1[[i]]$stat)
    
    datals$v2[[i]]$stat$sample = names(datals$v2)[i]
    
    datals$v2[[i]]$stat$group = "v2"
    statdf = rbind(statdf,datals$v2[[i]]$stat)
    statdf = rbind(statdf, data.frame("process" = c("v1_v2_combine","v1_v2_combine_blfil"),
                                      "count" = c(nrow(datals$arraydata[[i]]),nrow(datals$arraydata_bl[[i]])),
                                      "sample" = names(datals$v1)[i],"group" = "v1_v2"))
  }
  statdf$process = factor(statdf$process,levels = unique(statdf$process))
  library(ggsci)
  p1.0 = ggplot(statdf,aes(x = process, y = count, fill = group)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    geom_text(aes(x = process, y = count, label = count,angle = 90),size= 2,
              position = position_dodge(width=1)) +
    facet_wrap(~sample,nrow = 3) + theme_bw() + scale_fill_jama()+
    theme(axis.text.x =  element_text(angle = 90))
  p1.0
  # dir.create("final_result/CREST_final_analysis/0.qc")
  ggsave(p1.0,filename = "final_result/CREST_final_analysis/0.qc/qc_stat_barplot.pdf",
         width = 12,height = 8)
  
  datals$statdf = statdf
  
  # datals$figure$heatmap$E11$abspearson@column_order
  # datals$figure$heatmap$E15$abspearson
  # datals$figure$heatmap$OGN$abspearson
  
  
  
  arraydata  = datals$arraydata_bl$E11
  #CalHeatmap
  MyHeatmapCal = function(arraydata,title,order = NULL){
    library(amap)
    library(reshape2)
    library(ComplexHeatmap)
    library(RColorBrewer)
    cell_count = dcast(arraydata,pattern~Cell.type);
    rownames(cell_count) = cell_count$pattern;cell_count = cell_count[-1]
    cell_count = as.matrix(cell_count)
    
    cell_norm = cell_count[which(rowSums(cell_count>0)>1),]
    cell_norm = log(cell_norm+1)
    cell_norm = t(apply(cell_norm, 1, function(x) {x/sum(x)}))
    
    # pdf(paste0(outname,"_abspearson.pdf"),width = 10,height = 8)
    library(psych)
    mxpv = corr.test(cell_norm)
    mxpv = mxpv$p
    dd = Dist(t(cell_norm), method = "abspearson")
    mx = as.matrix(1 - dd)
    
    mxl = melt(mx)
    mxpvl = melt(mxpv)
    colnames(mxl)[3] = "corr"
    colnames(mxpvl)[3] = "pvalue"
    mxl = merge(mxl,mxpvl)
    
    for (i in 1:nrow(mx)) {
      mx[i,i] = 1
    }
    dd[is.na(dd)] = 0

    if(length(order)>0){
      p1 = Heatmap(mx,heatmap_legend_param = list(title = "abspearson"),row_title = title,
                   row_order = order[which(order %in% rownames(mx))],
                   column_order = order[which(order %in% colnames(mx))],
                   # cell_fun = function(j, i, x, y, w, h, fill){
                   #   if(mxpv[i, j] < 0.001) {
                   #     grid.text("***", x, y)
                   #   } else if(mxpv[i, j] < 0.01) {
                   #     grid.text("**", x, y)
                   #   } else if(mxpv[i, j] < 0.05) {
                   #     grid.text("*", x, y)
                   #   }
                   # },
                   # cluster_rows = F,cluster_columns = F,
                   col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)))
      mxl$Var1 = factor(mxl$Var1,levels = order)
      mxl$Var2 = factor(mxl$Var2,levels = rev(order))
    }else{
      p1 = Heatmap(mx,heatmap_legend_param = list(title = "abspearson"),row_title = title,
                   # cell_fun = function(j, i, x, y, w, h, fill){
                   #   if(mxpv[i, j] < 0.01) {
                   #     grid.text("*", x, y)
                   #   }
                   # },
                   # clustering_distance_rows = "spearman",
                   col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
                   clustering_distance_rows = function(x){x = dd; return(x)},
                   # clustering_distance_columns = "spearman",-
                   clustering_distance_columns = function(x){x = dd; return(x)},
                   clustering_method_columns = "complete",clustering_method_rows = "complete")
      
    }
    
    p2 = ggplot(mxl, aes(x=Var1,y=Var2, fill=corr,  height=sqrt(1-pvalue),  width=sqrt(1-pvalue))) +
      geom_tile()+
      scale_fill_distiller(palette = "RdYlGn", direction=-1, limits=c(0,1),values = c(0,0.05,0.15,1),name="abpearson") +
      scale_color_distiller(palette = "RdYlGn", direction=-1, limits=c(0,1),values = c(0,0.05,0.15,1),name="abpearson") +
      coord_fixed() + xlab("") + ylab("")+
      theme(axis.text.x=element_text(angle=90, vjust=0.5),
            panel.background=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
      )
    
    hc = hclust(dd, method = "complete")
    
    
    distfigure = list("abspearson" = p1,"abspearsonhc" = hc,"abspearson_withpv" = p2,"cormx" = mxl)
    return(distfigure)
  }

  heatmapls = list()
  for (i in 1:length(datals$arraydata_bl)) {
    heatmapls[[i]] = MyHeatmapCal(datals$arraydata_bl[[i]],names(datals$arraydata_bl)[i])
    names(heatmapls)[i] = names(datals$arraydata_bl)[i]
  }
  # dir.create("final_result/CREST_final_analysis/1.heatmap_tree")
  ggexport(heatmapls,filename = "final_result/CREST_final_analysis/1.heatmap_tree/heatmap_tree/total_heatmap.pdf",
           width = 10,height = 8)
  
  heatmapls$E11$abspearson_withpv
  heatmapls$E15$abspearson_withpv
  heatmapls$OGN$abspearson
  heatmapls$pandaE11$abspearson
  # datals = datals
  heatmapls_order = list()
  heatmaplsv1_order = list()
  heatmaplsv2_order = list()
  heatmapls_ordercm = list()
  heatmaplsv1_ordercm = list()
  heatmaplsv2_ordercm = list()
  
  
  orderE11 = c('NPBL','NbGABA0','NbGABABL','NbGLUAL1',
               'NbGABABI','NbGLUAL2','NPBM','NbBM0','NbBM1',
               'OMTN','Rgl1','NbFP','GLUAL','Peric')
  orderOGN = c('GABABL1','GABABL4','GABABL2','GLUAL1','GABABI','GLUMHB2','Rgl2BP',
               'NbGABABI','Rgl3','NbDA','RglUnk','NbGLU','DA',
               'GABA_unk',
               'GLUBM2','GLUBM1','GLUFP','NbBM1','NbFP')
  orderE15 = c('NbGABABL','GABABL2','GABABL3','GABABL4','GABABL1','Rgl3',
               'GLUAL4','GLUAL1','GLUAL2','GLUAL3','Rgl2AL','GABAAL','NbGABAAL',
               'Rgl2BP','GLUMHB1','GLUMHB2','GABAMHB','GLUBM2','GLUBM1',
               'NbBM1','NbGABABI','GABABI',
               'NbDA','DA','NbFP','Rgl1','NbGLU','GLUFP',
               'GLUMHB3','GLUMHB4')
  orderls = list(orderE11,orderE11,orderE11,orderE11,orderE11,orderE11,orderE15,orderE15,orderE15,
                 orderOGN,orderOGN,orderOGN,orderE11,orderE15,orderE11,orderOGN)
  names(orderls) = names(datals$arraydata_bl)
  
  commonctE11 = Reduce(intersect,list(datals$arraydata_bl$E11_rep1$Cell.type,datals$arraydata_bl$E11_rep2$Cell.type,datals$arraydata_bl$E11_rep3$Cell.type))
  commonctE15 = Reduce(intersect,list(datals$arraydata_bl$E15_rep1$Cell.type,datals$arraydata_bl$E15_rep2$Cell.type,datals$arraydata_bl$E15_rep3$Cell.type))
  commonctOGN = Reduce(intersect,list(datals$arraydata_bl$OGN_rep1$Cell.type,datals$arraydata_bl$OGN_rep2$Cell.type,datals$arraydata_bl$OGN_rep3$Cell.type))
  commonctpE11 = Reduce(intersect,list(datals$arraydata_bl$pandaE11_rep1$Cell.type,datals$arraydata_bl$pandaE11_rep2$Cell.type,datals$arraydata_bl$pandaE11_rep3$Cell.type))
  commonctls = list(commonctE11,commonctE11,commonctE11,commonctpE11,commonctpE11,commonctpE11,commonctE15,commonctE15,commonctE15,
                    commonctOGN,commonctOGN,commonctOGN,commonctE11,commonctE15,commonctpE11,commonctOGN)
  
  for (i in 1:length(datals$arraydata_bl)) {
    v1i = datals$v1[[names(datals$arraydata_bl)[i]]]$result_withcellan
    v1i = v1i[which(!v1i$pattern %in% bl$v1_array$Var1),]
    v2i = datals$v2[[names(datals$arraydata_bl)[i]]]$result_withcellan
    v2i = v2i[which(!v2i$pattern %in% bl$v2_array$Var1),]
    heatmaplsv1_order[[i]] = MyHeatmapCal(v1i,names(datals$arraydata_bl)[i],
                                        orderls[[i]])
    heatmaplsv2_order[[i]] = MyHeatmapCal(v2i,names(datals$arraydata_bl)[i],
                                          orderls[[i]])
    heatmapls_order[[i]] = MyHeatmapCal(datals$arraydata_bl[[i]],names(datals$arraydata_bl)[i],
                                        orderls[[i]])
    #common cell
    
    heatmaplsv1_ordercm[[i]] = MyHeatmapCal(v1i[which(v1i$Cell.type %in% commonctls[[i]]),],names(datals$arraydata_bl)[i],
                                          orderls[[i]])
    heatmaplsv2_ordercm[[i]] = MyHeatmapCal(v2i[which(v2i$Cell.type %in% commonctls[[i]]),],names(datals$arraydata_bl)[i],
                                          orderls[[i]])
    heatmapls_ordercm[[i]] = MyHeatmapCal(datals$arraydata_bl[[i]][which(datals$arraydata_bl[[i]]$Cell.type %in% commonctls[[i]]),],names(datals$arraydata_bl)[i],
                                        orderls[[i]])
    names(heatmapls_order)[i] = names(heatmaplsv1_order)[i] = names(heatmaplsv2_order)[i] =
      names(heatmapls_ordercm)[i] = names(heatmaplsv1_ordercm)[i] = names(heatmaplsv2_ordercm)[i] =
      names(datals$arraydata_bl)[i]
  }
  ggexport(heatmapls_order,filename = "final_result/CREST_final_analysis/1.heatmap_tree/total_heatmap_ordered.pdf",
           width = 10,height = 8)
  ggexport(heatmaplsv1_order,filename = "final_result/CREST_final_analysis/1.heatmap_tree/total_heatmap_ordered_v1.pdf",
           width = 10,height = 8)
  ggexport(heatmaplsv2_order,filename = "final_result/CREST_final_analysis/1.heatmap_tree/total_heatmap_ordered_v2.pdf",
           width = 10,height = 8)
  # datals$figure$heatmap_order = heatmapls_order
  datals$figure$heatmap = heatmapls
  datals$figure$heatmap_order = heatmapls_order
  datals$figure$heatmap_order_v1 = heatmaplsv1_order
  datals$figure$heatmap_order_v2 = heatmaplsv2_order
  datals$figure$heatmap_order_cm = list("v1" = heatmaplsv1_ordercm,"v2" = heatmaplsv2_ordercm,
                                        "v1_v2" = heatmapls_ordercm)
  datals$figure$qc_process_stat = p1.0
  datals$cellorder = orderls
  datals$blacklist = bl
  datals$cellls = cellls
  datals$arraydata_bl$E15
  
  heatmap_order_cm = list("v1" = heatmaplsv1_ordercm,"v2" = heatmaplsv2_ordercm,
                          "v1_v2" = heatmapls_ordercm)
  qsave(heatmap_order_cm,file = "final_result/CREST_final_analysis/heatmapls_order_cm.qs")
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_withfigure_01_10_edit.qs")
  
  datals = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_19_edit.qs")

  #indel tree
  oldtree =  readRDS("final_result/xE15_noDoublet/E15_bltg_tree.rds")
  tags = unique(oldtree$cm$tags)
  tagsls = strsplit(tags,"_|&")
  tags = data.frame("tags" = rep(tags,unlist(lapply(tagsls, length))),
                    "indels" = unlist(tagsls))
  indeldf = merge(oldtree$cm,tags,by = "tags")
  colnames(indeldf)[c(4,5)] = c("Cell.type","pattern")
  
  ordertmp = c("GluAL3","GluAL2","GluAL1", "GluAL4","Rgl3","GabaAL","Rgl2AL",
               "NProgAL","GabaProgBL",
               "GabaBL2","GabaBL3","GabaBL1","GabaBL4","GabaBI","Rgl2BL" ,"GabaProgBI",
               "GluBM1","NbBM1" ,"GluBM2","GluFP","NbGlu","NbFP","DA","NbDA","Rgl1")
  tmp = MyHeatmapCal(indeldf,"E15.5 InDel",order =ordertmp)
  pdf("final_result/CREST_final_analysis/rebbutal_analysis/E15_old_InDel_heatmap_ordered.pdf",
      width = 10,height = 8)
  print(tmp$abspearson)
  dev.off()
  
}


