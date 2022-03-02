#lineage reconstruction analysis
#2022/03/02 by lhx

source("lineage_tree_reconstruction_method.R")
#read data
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
  
  samplels = list.files("bulk_result/X.bulk_in_vivo/",full.names = T)[c(1:8,9,10,12:13)]
  v1lsb = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsb = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsb) = names(v2lsb) = list.files("bulk_result/X.bulk_in_vivo/")[c(1:8,9,10,12:13)]
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
  
  
  #scRNA
  filepath = list.files("sc_result/x.TAQ_E15/",full.names = T)[c(1,5,6)]
  cellpath = "sc_result/x.TAQ_E15/E15clustername20211126 v2.csv"
  ScLoad = function(filepath,cellpath){
    samplels = filepath
    v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
    v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
    names(v1lsx) = names(v2lsx) = unlist(lapply(strsplit(filepath,"/"),"[[",length(strsplit(filepath,"/")[[1]])))
    cell = ReadCell(cellpath)
    return(list("v1lsx" = v1lsx,"v2lsx" = v2lsx, "cell" = cell))
  }
  E15data = ScLoad(list.files("sc_result/x.TAQ_E15/",full.names = T)[c(1,5,6)],"sc_result/x.TAQ_E15/E15clustername20211126 v2.csv")
  E11data = ScLoad(list.files("sc_result/X.E11",full.names = T)[c(7,8,9)],"sc_result/X.E11/E11_singlet_renamed_celltype.csv")
  E11sdata = ScLoad(list.files("sc_result/STF-E11",full.names = T)[1:3],"sc_result/STF-E11/STF_E11_singlet_renamed_celltype.csv")
  Divdata = ScLoad(list.files("sc_result/STF-DIV7",full.names = T)[1:3],"sc_result/STF-DIV7/STF-DIV7_celltype.csv")
  
  v1lsb$scE15 = MyCombine(E15data$v1lsx)
  v1lsb$scE11 = MyCombine(E11data$v1lsx)
  v1lsb$scSTF = rbind(MyCombine(E11sdata$v1lsx),MyCombine(Divdata$v1lsx))
  
  v2lsb$scE15 = MyCombine(E15data$v2lsx)
  v2lsb$scE11 = MyCombine(E11data$v2lsx)
  v2lsb$scSTF = rbind(MyCombine(E11sdata$v2lsx),MyCombine(Divdata$v2lsx))
  
  for (i in (length(v1lsb)-2):length(v1lsb)) {
    v1lsb[[i]]$pattern = v1lsb[[i]]$umim
    v2lsb[[i]]$pattern = v2lsb[[i]]$umim
  }
  
  
}




#blacklist select
{
  SelectBlacklist1 = function(v1lsx){
    arraydft = NULL
    for (i in 2:length(v1lsx)) {
      arraydf = as.data.frame(table(v1lsx[[i]]$pattern))
      arraydf = arraydf[which(!arraydf$Var1 %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      
      arraydf = arraydf[order(-arraydf$Freq),]
      arraydf$Freq = arraydf$Freq/sum(arraydf$Freq)
      
      colnames(arraydf)[2] = names(v1lsx)[i]
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
  SelectBlacklist2 = function(v1lsx){
    tagdft = NULL
    for (i in 2:length(v1lsx)) {
      print(i)
      tag = NULL
      df = v1lsx[[i]]$pattern
      for (j in 1:length(df)) {
        x = df[j]
        x_str = strsplit(x,split = "_")
        x_pre = x_str[[1]]
        x_pre = unique(x_pre[!(x_pre %in% "NONE")])
        tag = c(tag,x_pre)
      }
      
      tagdf = as.data.frame(table(tag))
      tagdf[,2] = tagdf[,2]/sum(tagdf[,2])
      
      colnames(tagdf)[2] = names(v1lsx)[i]
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
  
  v1blayb = SelectBlacklist1(v1lsb)
  v2blayb = SelectBlacklist1(v2lsb)
  v1bltgb = SelectBlacklist2(v1lsb)
  v2bltgb = SelectBlacklist2(v2lsb)
  
  
  blt = list("v1_array" = v1blayb,"v2_array" = v2blayb,
             "v1_tag" = v1bltgb,"v2_tag" = v2bltgb)
  bl = list()
  for (i in 1:4) {
    line = blt[[i]]
    bl[[i]] = line[which(rowSums(line[13:15]>0.001)>1 | rowSums(line[10:12]>0.001)>1 | rowSums(line[2:9]>0.005)>1),]
    names(bl)[i] = names(blt)[i]
  }
  
  saveRDS(bl,"final_result/x.blacklist_0.001_2_25.rds")
  
  
}


#tree reconstruction
{

  FilterRecover = function(v1lsx,cell,thres,title){
  library(ggvenn)
  rawlist = list()
  rflist = list()
  baflist = list()
  for (i in 1:length(v1lsx)) {
    rawlist[[i]] = v1lsx[[i]]$Cell.BC
    v1lsx[[i]] = v1lsx[[i]][which(v1lsx[[i]]$reads_num >thres),]
    rflist[[i]] = v1lsx[[i]]$Cell.BC
  }
  
  
  comBC = Reduce(intersect,rflist)
  for (i in 1:length(v1lsx)) {
    v1lsx[[i]] = v1lsx[[i]][which(v1lsx[[i]]$Cell.BC %in% comBC),]
    baflist[[i]] = paste0(v1lsx[[i]]$Cell.BC,"-",v1lsx[[i]]$umim)
  }
  names(rawlist) = names(rflist) =  names(baflist) = names(v1lsx)
  
  # p1.0 = ggvenn(rawlist,stroke_alpha = 0.5,stroke_size = 0.5)
  # p1.1 = ggvenn(rflist,stroke_alpha = 0.5,stroke_size = 0.5)
  # p1.2 = ggvenn(baflist,stroke_alpha = 0.5,stroke_size = 0.5)
  
  comarray1 = Reduce(intersect,baflist)
  
  result = v1lsx[[1]]
  result$bcarray = paste0(result$Cell.BC,"-",result$umim)
  result = result[which(result$bcarray %in% comarray1),]
  result = merge(result[,c("Cell.BC","umim")],cell,by="Cell.BC");colnames(result)[2] = "pattern"
  stat = data.frame("process" = c("common_BC","common_array","common_with_cell"), 
                    "count" = c(length(comBC),length(comarray1),nrow(result)))
  
  p1.0 = ggvenn(rawlist,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[1])
  p1.1 = ggvenn(rflist,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[2])
  p1.2 = ggvenn(baflist,stroke_alpha = 0,stroke_size = 0.5) + ggtitle(title[3])
  p = ggarrange(p1.0,p1.1,p1.2,nrow = 3)
  
  return(list("result" = result,"rawlist" = rawlist,"rflist" = rflist,"baflist" = baflist,"stat" = stat,"venn_plot" = p))
  }
  
  v1tmp = E15data$v1lsx$AT
  v2tmp = E15data$v2lsx$AT
  v1tmp$pattern = v1tmp$umim
  v2tmp$pattern = v2tmp$umim
  bl1 = bl$v1_array
  bl2 = bl$v2_array
  test = BlacklistFilter2(v1tmp,v2tmp,bl1,bl2)
  BlacklistFilter2 = function(v1tmp,v2tmp,bl1,bl2){
    newdf = NULL
    df = merge(v1tmp[,c("Cell.BC","pattern")],v2tmp[,c("Cell.BC","pattern")],by = "Cell.BC",all = T)
    # df = df[which(df$reads_pro.umi >= 0.5 & df$umi_num >= 3),]
    for (j in 1:nrow(df)) {
      flag = 0
      x1 = df$pattern.x[j]
      flag = length(bl1[,1][which(bl1[,1] == x1)])
      if(flag == 0 & !is.na(x1)){
        x1_str = strsplit(x1,split = "_")
        x1_pre = x1_str[[1]]
        x1_pre = unique(x1_pre[!(x1_pre %in% "NONE")])
        flag = length(intersect(x1_pre,bl1[,1]))
      }
      if(flag==0 & !is.na(x1)){
        newdf = rbind(newdf,df[j,])
      }else{
        x2 = df$pattern.y[j]
        flag = length(bl2[,1][which(bl2[,1] == x2)])
        if(flag == 0 & !is.na(x2)){
          x2_str = strsplit(x2,split = "_")
          x2_pre = x2_str[[1]]
          x2_pre = unique(x2_pre[!(x2_pre %in% "NONE")])
          flag = length(intersect(x2_pre,bl2[,1]))
        }
        if(flag==0 & !is.na(x2)){
          newdf = rbind(newdf,df[j,])
        }
      }
      
    }
    newdata = list("v1" = v1tmp[which(v1tmp$Cell.BC %in% newdf$Cell.BC),],"v2" = v2tmp[which(v2tmp$Cell.BC %in% newdf$Cell.BC),])
    return(newdata)
    
  }
  XBulidTree = function(v1lsx,v2lsx,cell11,outpath,bl){
    print("qc")
    prefix = c("v1.","v2.")
    
    v1title = c("V1 raw BC recover rate","V1 reads num(10) filtered BC recover rate",
                "V1 reads num(10) filtered BC-pattern recover rate")
    v1test11 = FilterRecover(v1lsx,cell11,10,v1title)
    v1test11$stat
    # v1test11$venn_plot
    
    v2title = c("V2 raw BC recover rate","V2 reads num(10) filtered BC recover rate",
                "V2 reads num(10) filtered BC-pattern recover rate")
    v2test11 = FilterRecover(v2lsx,cell11,10,v2title)
    v2test11$stat
    # v2test11$venn_plot
    saveRDS(v1test11,file = paste0(outpath,"_v1qc.rds"))
    saveRDS(v2test11,file = paste0(outpath,"_v2qc.rds"))
    # v1test11 = readRDS("final_result/STF_E11/e11_v1qc.rds")
    # v2test11 = readRDS("final_result/STF_E11/e11_v2qc.rds")
    
    ggsave(ggarrange(v1test11$venn_plot,v2test11$venn_plot,nrow = 1),
           filename = paste0(outpath,"_recover_rate_stat.pdf"),
           width = 14,height = 10)
    
    v1testf11 = BlacklistFilter(v1test11$result,bl$v1_array)
    v2testf11 = BlacklistFilter(v2test11$result,bl$v2_array)
    
    print("tree construct")
    dataE11 = BlacklistFilter2(v1test11$result,v2test11$result,bl$v1_array,bl$v2_array)
    outnameE11 = paste0(outpath,"_blay")
    E11tree = TreeConstruct(dataE11,prefix,outnameE11,cell11)
    E11tree = CirclePlot1(E11tree,outnameE11)
    E11tree = MyHeatmap(E11tree,outnameE11)
    
    saveRDS(E11tree,file = paste0(outpath,"_blay_tree.rds"))
    
    v1testf11tg = BlacklistFilter(v1test11$result,bl$v1_tag)
    v2testf11tg = BlacklistFilter(v2test11$result,bl$v2_tag)
    dataE11tg = BlacklistFilter2(v1test11$result,v2test11$result,bl$v1_tag,bl$v2_tag)
    outnameE11 = paste0(outpath,"_bltg")
    E11treetg = TreeConstruct(dataE11tg,prefix,outnameE11,cell11)
    E11treetg = CirclePlot1(E11treetg,outnameE11)
    E11treetg = MyHeatmap(E11treetg,outnameE11)
    saveRDS(E11treetg,file = paste0(outpath,"_bltg_tree.rds"))
    
    
    prefix = c("v1.")
    dataE11v1 = list(v1testf11)
    outnameE11v1 = paste0(outpath,"_blay_v1")
    E11treev1 = TreeConstruct(dataE11v1,prefix,outnameE11v1,cell11)
    E11treev1 = CirclePlot1(E11treev1,outnameE11v1)
    E11treev1 = MyHeatmap(E11treev1,outnameE11v1)
    saveRDS(E11treev1,file = paste0(outpath,"_blay_v1_tree.rds"))
    
    prefix = c("v2.")
    dataE11v2 = list(v2testf11)
    outnameE11v2 = paste0(outpath,"_blay_v2")
    E11treev2 = TreeConstruct(dataE11v2,prefix,outnameE11v2,cell11)
    E11treev2 = CirclePlot1(E11treev2,outnameE11v2)
    # E11treev2 = readRDS("final_result/STF_E11/E11_blay_v2_tree.rds")
    E11treev2 = MyHeatmap(E11treev2,outnameE11v2)
    saveRDS(E11treev2,file = paste0(outpath,"_blay_v2_tree.rds"))
    return(E11tree)
  }
  
  test = BlacklistFilter(v1test11$result,bl$v1_array)
  
  outpath = "final_result/xdraft/tree_fil_clonenum1_3_2/E15_bl0.001"
  E15tree = XBulidTree(E15data$v1lsx,E15data$v2lsx,E15data$cell,outpath,bl)
  E15tree2 = readRDS("final_result/xE15.5/E15.5_blay_taq_tree.rds")
  E11tree2 = readRDS("final_result/xE11/E11_12_21_blay_tree.rds")
  outpath = "final_result/xdraft/tree_fil_clonenum1_3_2/E11_bl0.001"
  E11tree = XBulidTree(E11data$v1lsx,E11data$v2lsx,E11data$cell,outpath,bl)
  outpath = "final_result/xdraft/tree_fil_clonenum1_3_2/E11_STF_bl0.001"
  E11stree = XBulidTree(E11sdata$v1lsx,E11sdata$v2lsx,E11sdata$cell,outpath,bl)
  outpath = "final_result/xdraft/tree_fil_clonenum1_3_2/Div_STF_bl0.001"
  Divtree = XBulidTree(Divdata$v1lsx,Divdata$v2lsx,Divdata$cell,outpath,bl)
  
  
  
}

