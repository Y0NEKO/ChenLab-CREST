#2022.11.26
#revise article
setwd("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/")
.libPaths(c("/home/liuhengxin/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library",
            "/picb/sysgenomics2/people/liuhengxin/software/anaconda3/envs/r4-base/lib/R/library/",
            "/picb/sysgenomics2/people/liuhengxin/software/miniconda3/envs/r4.1/lib/R/library/"))
#0. Basic Function
{
  VocanoFigure = function(DA.markers,title,highlight,thres){
    DA.markers$reg = "NoSig"
    DA.markers[which(DA.markers$p_val<0.05 & DA.markers$avg_log2FC > thres),"reg"] = "Up"
    DA.markers[which(DA.markers$p_val<0.05 & DA.markers$avg_log2FC < -thres),"reg"] = "Down"
    DA.markers$gene =rownames(DA.markers)
    DA.markers$reg = factor(DA.markers$reg,levels = c("NoSig","Down","Up"))
    # library(ggrepel)
    p = ggplot(DA.markers,aes(x = avg_log2FC,y = -log10(p_val),color = reg)) + geom_point() + 
      scale_color_manual(values=c("#808080","#DC143C","#00008B"))+#确定点的颜色
      ggrepel::geom_text_repel(
        data = DA.markers[which((DA.markers$p_val<0.05 & abs(DA.markers$avg_log2FC) > thres)),],
        # data = DA.markers[which(rownames(DA.markers) %in% highlight),],
        aes(label = gene),
        size = 3,
        segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
      theme_bw()+#修改图片背景
      theme(
        legend.title = element_blank()#不显示图例标题
      )+
      ylab('-log10 (pvalue)')+#修改y轴名称
      xlab('log2 (FoldChange)')+#修改x轴名称
      ggtitle(title) +
      geom_vline(xintercept=c(-thres,thres),lty=3,col="black",lwd=thres) +#添加横线|FoldChange|>2
      geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=thres)#添加竖线padj<0.05
    return(p)
  }
  
  #upset pvalue
  CaldmUpalue = function(e11dmup){
    e11pv = colSums(e11dmup[-1])/nrow(e11dmup)
    rownames(e11dmup) = e11dmup$X;e11dmup = e11dmup[-1]
    library(tidyr)
    e11pvmx = crossing(AL = 0:1, BL = 0:1, BI = 0:1, BM = 0:1, FP = 0:1)
    e11pvmx$pe = 0
    e11pvmx$count = 0
    e11pvmx$expect = 0
    e11pvmx$pvalue = 0
    sum(e11pvmx$pe)
    for (i in 2:nrow(e11pvmx)) {
      line = e11pvmx[i,1:5]
      linei = which(line>0)
      pe = prod(e11pv[linei],1-e11pv[-linei])
      vect = e11dmup[which(rowSums(e11dmup[linei] > 0) == length(linei) & rowSums(e11dmup[-linei] > 0) == 0),]
      e11pvmx$pe[i] = pe
      e11pvmx$count[i] = nrow(vect)
      
      
      
    }
    
    e11pvmx$pe = e11pvmx$pe/sum(e11pvmx$pe)
    e11pvmx$expect = e11pvmx$pe*nrow(e11dmup)
    
    for (i in 2:nrow(e11pvmx)) {
      tmp = binom.test(e11pvmx$count[i], nrow(e11dmup), p = e11pvmx$pe[i], 
                       alternative = "two.sided",conf.level= 0.95)
      e11pvmx$pvalue[i] = tmp$p.value
    }
    e11pvmx$FDP = p.adjust(e11pvmx$pvalue)
    return(e11pvmx)
  }
}

#1. rebuttal before 11.23------------
{
  divv1 = readRDS("final_result/STF_DIV7/DIV7_v1qc.rds")
  divv2 = readRDS("final_result/STF_DIV7/DIV7_v2qc.rds")
  E11Sv1 = readRDS("final_result/STF_E11/singlet_E11_v1qc.rds")
  E11Sv2 = readRDS("final_result/STF_E11/singlet_E11_v2qc.rds")
  e11v1 = readRDS("final_result/xE11_noDoublet/E11_12_21_v1qc.rds")
  e11v2 = readRDS("final_result/xE11_noDoublet/E11_12_21_v2qc.rds")
  e15v1 = readRDS("final_result/xE15_noDoublet/E15_v1qc.rds")
  e15v2 = readRDS("final_result/xE15_noDoublet/E15_v2qc.rds")
  celle11 = ReadCell("sc_result/X.E11/E11_singlet_renamed_celltype.csv")
  celldiv = ReadCell("sc_result/STF-DIV7/STF-DIV7_celltype.csv")
  cells11 = ReadCell("sc_result/STF-E11/STF_E11_singlet_renamed_celltype.csv")
  celle15 = ReadCell("sc_result/x.TAQ_E15/E15clustername20211126 v2.csv")
  doublet = read.csv("sc_result/x.TAQ_E15/E15.5_doublet_BC.csv")$x
  celle15 = celle15[which(!celle15$Cell.BC %in% substr(doublet,1,16)),]
  
  v1qc1 = readRDS("raw_data/E15_vitro_scar/lib1_v1qc.rds")
  v2qc1 = readRDS("raw_data/E15_vitro_scar/lib1_v2qc.rds")
  v1qc2 = readRDS("raw_data/E15_vitro_scar/lib2_v1qc.rds")
  v2qc2 = readRDS("raw_data/E15_vitro_scar/lib2_v2qc.rds")
  cellvitro = read.csv("raw_data/E15_vitro_scar/220923_metadata.csv")
  colnames(cellvitro)[c(1,13)] = c("Cell.BC","Cell.type")
  cellvitro1 = cellvitro[which(substr(cellvitro$Cell.BC,18,18)=="1"),]
  cellvitro2 = cellvitro[which(substr(cellvitro$Cell.BC,18,18)=="2"),]
  cellvitro1$Cell.BC = substr(cellvitro1$Cell.BC,1,16)
  cellvitro2$Cell.BC = substr(cellvitro2$Cell.BC,1,16)
  
  #(1)stat barcode edit
  tedit = c(e11v2$result$pattern,e15v2$result$pattern,divv2$result$pattern,E11Sv2$result$pattern)
  editv1 = list(e11v1$result$pattern,e15v1$result$pattern,divv1$result$pattern,E11Sv1$result$pattern)
  editv2 = list(e11v2$result$pattern,e15v2$result$pattern,divv2$result$pattern,E11Sv2$result$pattern)
  bl = readRDS("final_result/x.blacklist_0.001_2_25.rds")
  Editstat = function(tedit){
    edit_stat = NULL
    edit = unique(tedit)
    for(i in 1:length(edit)){
      line = edit[i]
      line_size = length(tedit[which(tedit == line)])
      line_str = unique(unlist(str_split(line, "_|&")))
      line_str = line_str[which(line_str != "NONE")]
      if(length(line_str)>0){
        line_del = line_str[grep("D",line_str)]
        line_ins = line_str[grep("I",line_str)]
        line_deln = unlist(str_extract_all(line_del, "[0-9]+D"))
        line_insn = unlist(str_extract_all(line_ins, "[0-9]+I"))
        line_start = c(unlist(str_extract_all(line_del, "\\+[0-9]+")),unlist(str_extract_all(line_ins, "\\+[0-9]+")))
        line_start = as.numeric(substr(line_start,2,nchar(line_start)))
        line_deln = as.numeric(substr(line_deln,1,nchar(line_deln)-1))
        line_insn = as.numeric(substr(line_insn,1,nchar(line_insn)-1))
        edit_stat =rbind(edit_stat,data.frame("array" = line,"clone_size" = line_size,
                                              "edit_num" = length(line_str),'width' = c(line_deln,line_insn),"start" = line_start,
                                              "type" = rep(c("deletion","insertion"),c(length(line_del),length(line_ins)))))
      }
      
    }
    return(edit_stat)
  }
  
  # summary(unique(edit_stat[1:2])$clone_size)
  # summary(unique(edit_stat[which(edit_stat$type == "deletion"),])$width)
  # summary(unique(edit_stat[which(edit_stat$type == "insertion"),])$width)
  # table(edit_stat$type)/nrow(edit_stat)
  
  #(2) relation ship of edit num and clone size
  library(ggExtra)
  names(editv1) = names(editv2) = c("E11","E15","DIV","pandaE11")
  for (i in 1:4) {
    edit_stat1 = Editstat(editv1[[i]])
    edit_stat2 = Editstat(editv2[[i]])
    edit_stat1 = edit_stat1[which(edit_stat1$type == "deletion"),] %>% group_by(array) %>% 
      summarise(clone_size = unique(clone_size),width = sum(width),edit_num = unique(edit_num))
    edit_stat1 = edit_stat1[which(!edit_stat1$array %in% bl$Var1),]
    edit_stat2 = edit_stat2[which(edit_stat2$type == "deletion"),] %>% group_by(array) %>% 
      summarise(clone_size = unique(clone_size),width = sum(width),edit_num = unique(edit_num))
    edit_stat2 = edit_stat2[which(!edit_stat2$array %in% bl$Var1),]
    
    p2.4.1 = ggstatsplot::ggbetweenstats(data = edit_stat1,x = "edit_num", y = "clone_size") + 
      ggtitle(paste0(names(editv1)[i],"V1"))
    p2.4.2 = ggstatsplot::ggbetweenstats(data = edit_stat2,x = "edit_num", y = "clone_size") + 
      ggtitle(paste0(names(editv1)[i],"V2"))
    p2.4 = ggarrange(p2.4.1,p2.4.2,nrow = 2) 
    ggsave(p2.4,filename = paste0("final_result/xtotal/",names(editv1)[i],"_editnum_clonesize_barplot.pdf"),
           width = 10,height = 10)
  }
  
  
  edit_stat$clone_size = log10(edit_stat$clone_size)
  
  
  p2.0 = ggplot(edit_stat,aes(x = edit_num, y =  log10(clone_size),size = width)) + 
    scale_size(range = c(1,3)) +
    geom_jitter(alpha = 0.5) + theme(legend.position = "") +
    theme_bw()
  ggMarginal(p2.0, type="density")
  
  p2.1 = ggplot(edit_stat,aes(x = width, y =  log10(clone_size),size = edit_num,color = type)) + 
    scale_size(range = c(1,3)) +
    geom_jitter(alpha = 0.5) + theme(legend.position = "left") +
    theme_bw()
  ggMarginal(p2.1, type="density")
  
  p2.2 = ggplot(edit_stat,aes(x = width, y =  edit_num,size = log10(clone_size),color = type)) + 
    scale_size(range = c(1,3)) +
    geom_jitter(alpha = 0.5) + theme(legend.position = "left") +
    theme_bw()
  ggMarginal(p2.2, type="density")
  ggarrange(ggMarginal(p2.0, type="density"),ggMarginal(p2.2, type="density"),nrow = 2) 
  
  
  #E15.5 vitro recover--
  #load lib
  {
    cellvitro = read.csv("raw_data/E15_vitro_scar/220923_metadata.csv")
    samplels1 = list.files("raw_data/E15_vitro_scar/",full.names = T)[1:3]
    colnames(cellvitro)[c(1,13)] = c("Cell.BC","Cell.type")
    cellvitro1 = cellvitro[which(substr(cellvitro$Cell.BC,18,18)=="1"),]
    cellvitro2 = cellvitro[which(substr(cellvitro$Cell.BC,18,18)=="2"),]
    cellvitro1$Cell.BC = substr(cellvitro1$Cell.BC,1,16)
    cellvitro2$Cell.BC = substr(cellvitro2$Cell.BC,1,16)
    v1lsx1 = lapply(paste0(samplels1,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
    v2lsx1 = lapply(paste0(samplels1,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
    names(v1lsx1) = names(v2lsx1) = list.files("raw_data/E15_vitro_scar/")[1:3]
    v1title = c("E15 vitro lib1 V1 raw BC recover rate","E15 vitro lib1 V1 reads num(10) filtered BC recover rate",
                "E15 vitro lib1 V1 reads num(10) filtered BC-pattern recover rate")
    v1qc1 = FilterRecover(v1lsx1,cellvitro1,10,v1title)
    v1qc1$stat
    # v1test11$venn_plot
    v2title = c("E15 vitro lib2 V1 raw BC recover rate","E15 vitro lib2 V1 reads num(10) filtered BC recover rate",
                "E15 vitro lib2 V1 reads num(10) filtered BC-pattern recover rate")
    v2qc1 = FilterRecover(v2lsx1,cellvitro1,10,v2title)
    v2qc1$stat
    # v2test11$venn_plot
    saveRDS(v1qc1,file = paste0("raw_data/E15_vitro_scar/lib1","_v1qc.rds"))
    saveRDS(v2qc1,file = paste0("raw_data/E15_vitro_scar/lib1","_v2qc.rds"))
    
    #lib2
    samplels2 = list.files("raw_data/E15_vitro_scar/",full.names = T)[6:8]
    v1lsx2 = lapply(paste0(samplels2,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
    v2lsx2 = lapply(paste0(samplels2,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
    names(v1lsx2) = names(v2lsx2) = list.files("raw_data/E15_vitro_scar/")[6:8]
    v1title = c("E15 vitro lib2 V1 raw BC recover rate","E15 vitro lib2 V1 reads num(10) filtered BC recover rate",
                "E15 vitro lib2 V1 reads num(10) filtered BC-pattern recover rate")
    v1qc2 = FilterRecover(v1lsx2,cellvitro2,10,v1title)
    v1qc2$stat
    # v1test11$venn_plot
    v2title = c("E15 vitro lib2 V2 raw BC recover rate","E15 vitro lib2 V2 reads num(10) filtered BC recover rate",
                "E15 vitro lib2 V2 reads num(10) filtered BC-pattern recover rate")
    v2qc2 = FilterRecover(v2lsx2,cellvitro2,10,v2title)
    v2qc2$stat
    # v2test11$venn_plot
    saveRDS(v1qc2,file = paste0("raw_data/E15_vitro_scar/lib2","_v1qc.rds"))
    saveRDS(v2qc2,file = paste0("raw_data/E15_vitro_scar/lib2","_v2qc.rds"))
    
  }
  Recoverstat = function(v1qc1,v2qc1,cellvitro1,sample){
    v1bc1 = Reduce(intersect,v1qc1$rawlist)
    v2bc1 = Reduce(intersect,v2qc1$rawlist)
    v1v2bc1 = intersect(v1bc1,v2bc1)
    reBC1 = data.frame("BC" = cellvitro1$Cell.BC,
                       "recover" = "unrecover") 
    reBC1[which(reBC1$BC %in% v1bc1),"recover"] = "V1"
    reBC1[which(reBC1$BC %in% v2bc1),"recover"] = "V2"
    reBC1[which(reBC1$BC %in% v1v2bc1),"recover"] = "V1_V2"
    reBC1 = merge(reBC1,cellvitro1[c("Cell.BC","Cell.type")],by.x = "BC",by.y = "Cell.BC")
    reBC1.st = reBC1 %>% group_by(Cell.type,recover) %>% summarise(cell.counts = n())
    reBC1.st$Cell.type = as.character(reBC1.st$Cell.type)
    reBC1.st$sample = sample
    return(reBC1.st)
  }
  reBC1.st = Recoverstat(v1qc1,v2qc1,cellvitro1,"lib1")
  reBC2.st = Recoverstat(v1qc2,v2qc2,cellvitro2,"lib2")
  reBCe11.st = Recoverstat(e11v1,e11v2,celle11,"E11")
  reBCe15.st = Recoverstat(e15v1,e15v2,celle15,"E15")
  reBCes11.st = Recoverstat(E11Sv1,E11Sv2,cells11,"STF_E11")
  reBCediv.st = Recoverstat(divv1,divv2,celldiv,"DIV7")
  reBCrt.st = rbind(reBC1.st,reBC2.st)
  reBCrt.st$sample = "Retro_E15"
  
  reBC = rbind(reBCe11.st,reBCe15.st,reBCes11.st,reBCediv.st,reBCrt.st)
  # reBC$Cell.type = factor(reBC$Cell.type,levels = as.character(0:26))
  p4.0 = ggplot(reBC,aes(x = Cell.type,y = cell.counts,fill = recover)) + 
    geom_bar(position="fill", stat="identity",width = 0.5) +
    facet_grid(~sample,scales = "free") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90))
  p4.0
  
  reBCt = reBC %>% group_by(sample,recover) %>% 
    summarise(cell.counts = sum(cell.counts))
  reBCt = reBCt %>% group_by(sample) %>% 
    summarise(cell.counts = cell.counts, recover = recover,
              cell.pos = cell.counts/sum(cell.counts))
  
  # reBCt$cell.pos = reBCt$cell.counts/sum(reBCt$cell.counts)
  reBCt$cell.pos.label = paste0(round(reBCt$cell.pos,3)*100,"%")
  reBCt$sample = factor(reBCt$sample,levels = c("E11","E15","STF_E11","DIV7","Vitro_E15"))
  p4.1 = ggplot(reBCt,aes(x = sample, y = cell.counts, fill = recover)) + 
    geom_bar(position="fill", stat="identity",width = 0.5) +
    geom_text(aes(label=cell.pos.label),position = position_fill(vjust = 0.5)) +
    theme_bw()
  p4.1
  
  library(Seurat)
  p4.t = ggarrange(p4.0+NoLegend(),p4.1,widths = c(2,1))
  ggsave(p4.0,filename = "final_result/xtotal/rebuttal_total_recover_rate_unfiltered_per_cell.pdf",width = 12,
         height = 5)
  ggsave(p4.1,filename = "final_result/xtotal/rebuttal_total_recover_rate_unfiltered.pdf",width = 8,
         height = 5)
  
  
  #heatmap
  # celldiv = ReadCell("sc_result/STF-DIV7/STF-DIV7_celltype.csv")
  outpath1 = "final_result/xtotal/rebuttallib1/rebuttal_E15_vitrolib1"
  outpath2 = "final_result/xtotal/rebuttallib2/rebuttal_E15_vitrolib2"
  lib1tree = XBulidTree(v1lsx1,v2lsx1,cellvitro1,outpath1,bl)
  lib2tree = XBulidTree(v1lsx2,v2lsx2,cellvitro2,outpath2,bl)
  
  #edit pos and indel ribbon
  #画每个segment的编辑diversity
  StatIndelLen = function(v1ls,v2ls){
    indelen = NULL
    indelseg = list()
    for (i in 1:length(v1ls)) {
      v1 = v1ls[[i]]
      v2 = v2ls[[i]]
      v1indelseg = ScartoIndel(v1$pattern)
      v2indelseg = ScartoIndel(v2$pattern)
      v1indelseg$sample = v2indelseg$sample = names(v1ls)[i]
      v1indelseg$group = "V1";v2indelseg$group = "V2"
      line2 = rbind(v1indelseg,v2indelseg)
      indelseg[[i]] = line2
      names(indelseg)[i] = names(v1ls)[i]
      # v1 = v1[v1$umim != "NONE_NONE_NONE_NONE_NONE_NONE_NONE",]
      # v2 = v2[v2$umim != "NONE_NONE_NONE_NONE_NONE_NONE_NONE",]
      v1spl = strsplit(v1$pattern,"_|&")
      v1taglen = unlist(lapply(v1spl, function(x){length(unique(x[which(x!="NONE")]))}))
      v2spl = strsplit(v2$pattern,"_|&")
      v2taglen = unlist(lapply(v2spl, function(x){length(unique(x[which(x!="NONE")]))}))
      line = data.frame(sample = names(v1ls)[i],indelenth = c(v1taglen,v2taglen),
                        group = rep(c("V1","V2"),c(length(v1taglen),length(v2taglen))))
      
      indelen = rbind(indelen,line)
      
      
    }
    return(list(indelen,indelseg))
    
  }
  ScartoIndel = function(tl){
    indel = NULL
    #build indel data frame
    for(i in c(1:length(tl))){
      print(i)
      x = tl[i]
      x_str = strsplit(x,split = "_")[[1]]
      for (j in 1:length(x_str)) {
        x_str_seg = strsplit(x_str[j],split = "&")[[1]]
        segid = j
        # if(x_str_seg[1] != "NONE"){
        if(x_str_seg[1] != "NONE"){
          for (k in 1:length(x_str_seg)) {
            # x_pre_num = str_extract_all(x_str_seg[k], "[0-9]+")
            # x_pre_num = as.numeric(x_pre_num[[1]])
            # start = x_pre_num[length(x_pre_num)]
            # end = start + x_pre_num[1]
            # segid.end = min(cutsite[which(cutsite$V3>end),"V1"])
            if(length(x_str_seg[k][grep("I",x_str_seg[k])])>0){
              indeltype = "Insertion"
            }else if(j == 7){
              if(x_str_seg[k] %in% strsplit(x_str[j-1],split = "&")[[1]]){
                indeltype = "MultipleDeletion"
              }else{
                indeltype = "SingleDeletion"
              }
              
            }else if(j == 1){
              if(x_str_seg[k] %in% strsplit(x_str[j+1],split = "&")[[1]]){
                indeltype = "MultipleDeletion"
              }else{
                indeltype = "SingleDeletion"
              }
              
            }else if(x_str_seg[k] %in% strsplit(x_str[j+1],split = "&")[[1]] |
                     x_str_seg[k] %in% strsplit(x_str[j-1],split = "&")[[1]]){
              indeltype = "MultipleDeletion"
            }else{
              indeltype = "SingleDeletion"
            }
            
            indelline = data.frame(id = j,array = x,indel = x_str_seg[k],indeltype = indeltype)
            indel = rbind(indel,indelline)
            
          }
        }else{
          indelline = data.frame(id = segid,array = x,indel = "NONE",indeltype = "Unedit")
          indel = rbind(indel,indelline)
        }
        
        
      }
    }
    return(indel)
  }
  
  samplels = list.files("bulk_result/X.bulk_in_vivo/",full.names = T)[c(1:8,9,10,12:13)]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = names(v2lsx) = list.files("bulk_result/X.bulk_in_vivo/")[c(1:8,9,10,12:13)]
  names(v1lsx) = names(v2lsx) = c("ctrl","E10-1","E10-2","E12","E15-1a","E15-1b",
                                  "E15-2a","E15-2b","E15-3a","E15-3c","E9","P0")
  v1lsx = v1lsx[c(5:10)]
  v2lsx = v2lsx[c(5:10)]
  v1lsx$E15_1 = rbind(v1lsx$`E15-1a`,v1lsx$`E15-1b`)
  v1lsx$E15_2 = rbind(v1lsx$`E15-2a`,v1lsx$`E15-2b`)
  v1lsx$E15_3 = rbind(v1lsx$`E15-3a`,v1lsx$`E15-3b`)
  v1lsx$E15_1$pattern = v1lsx$E15_1$umim;v1lsx$E15_2$pattern = v1lsx$E15_2$umim;v1lsx$E15_3$pattern = v1lsx$E15_3$umim;
  
  v2lsx$E15_1 = rbind(v2lsx$`E15-1a`,v2lsx$`E15-1b`)
  v2lsx$E15_2 = rbind(v2lsx$`E15-2a`,v2lsx$`E15-2b`)
  v2lsx$E15_3 = rbind(v2lsx$`E15-3a`,v2lsx$`E15-3b`)
  v2lsx$E15_1$pattern = v2lsx$E15_1$umim;v2lsx$E15_2$pattern = v2lsx$E15_2$umim;v2lsx$E15_3$pattern = v2lsx$E15_3$umim;
  
  v1ls = list("E11" = e11v1$result,"E15" = e15v1$result,"STF_E11" = E11Sv1$result,"DIV" = divv1$result,
              "Vitro_E15" = rbind(v1qc1$result,v1qc2$result),"Bulk_E15_1" = v1lsx$E15_1,"Bulk_E15_2" = v1lsx$E15_2,"Bulk_E15_3" = v1lsx$E15_3)
  v2ls = list("E11" = e11v2$result,"E15" = e15v2$result,"STF_E11" = E11Sv2$result,"DIV" = divv2$result,
              "Vitro_E15" = rbind(v2qc1$result,v2qc2$result),"Bulk_E15_1" = v2lsx$E15_1,"Bulk_E15_2" = v2lsx$E15_2,"Bulk_E15_3" = v2lsx$E15_3)
  indelseg = StatIndelLen(v1ls,v2ls)
  indelsegbulk = StatIndelLen(v1ls[6:8],v2ls[6:8])
  indelseg = readRDS("final_result/xtotal/indelseg_raw.rds")
  indelseg[[2]] = append(indelseg[[2]], indelsegbulk[[2]])
  indelsegf = indelseg[[2]]
  for (i in 1:length(indelsegf)) {
    indelsegf[[i]] = indelsegf[[i]][which(!indelsegf[[i]]$array %in% c(as.character(bl$v1_array$Var1), as.character(bl$v2_array$Var1))),]
  }
  # indelseg = StatIndelLen(v1ls,v2ls)
  
  saveRDS(indelseg,"final_result/xtotal/indelseg_raw.rds")
  indelseg = readRDS("final_result/xtotal/indelseg_raw.rds")
  
  
  EditTypeStat = function(indelseg,barpos = "stack"){
    p2.3 = list()
    for (i in 1:length(indelseg)) {
      line = indelseg[[i]]
      line$indeltype = factor(line$indeltype,levels = c("Unedit","Insertion","MultipleDeletion",
                                                        "SingleDeletion"))
      # line = line[which(line$indel != 'NONE'),]
      segidden = line %>% group_by(group,id,indeltype) %>%summarise(diversity = length(unique(indel)))
      
      p2.3.0 = ggplot(segidden,aes(x = id, y = diversity,fill = indeltype)) + 
        geom_bar(color = "black",stat = "identity",position = barpos) + facet_grid(~group,scales = "free") + 
        theme_pubr() + scale_x_continuous(breaks = c(unique(indelseg[[i]]$id))) +
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127"))+ xlab("") +
        ggtitle(names(indelseg)[i])
      p2.3.0
      p2.3[[i]] = p2.3.0
    }
    return(p2.3)
  }
  p2.3.1 = EditTypeStat(indelseg[[2]])
  p2.3.2 = EditTypeStat(indelseg[[2]],barpos = "fill")
  p2.3.3 = EditTypeStat(indelsegf)
  p2.3.4 = EditTypeStat(indelsegf,barpos = "fill")
  ggexport(p2.3.1,filename = "final_result/xtotal/segment_indel_diversity_stat.pdf",width = 6,height = 4)
  ggexport(p2.3.2,filename = "final_result/xtotal/segment_indel_diversity_stat_fill.pdf",width = 6,height = 4)
  ggexport(p2.3.3,filename = "final_result/xtotal/segment_indel_diversity_stat_blfiltered.pdf",width = 6,height = 4)
  ggexport(p2.3.4,filename = "final_result/xtotal/segment_indel_diversity_stat_blfiltered_fill.pdf",width = 6,height = 4)
  
  #plot ribbon
  tmp = Editstat(v1ls$Bulk_E15_1$pattern)
  e11vf = BlacklistFilter2(e11v1$result,e11v2$result,bl$v1_array,bl$v2_array)
  INDELRibbon = function(pattern,bl,scarrefv1){
    edit_stat_e11 = Editstat(pattern)
    ind_allsite_per = NULL
    for (i in 1:nrow(edit_stat_e11)) {
      print(i)
      line = edit_stat_e11[i,]
      line = line[rep(1,line$width),]
      line$pos = unique(line$start):(unique(line$start)+unique(line$width)-1)
      # linerep = line[rep(1:nrow(line),unique(line$clone_size)),]
      ind_allsite_per = rbind(ind_allsite_per,line)
    }
    
    ind_allsite_sum_t = ind_allsite_per %>% 
      group_by(pos,type) %>% summarise(num = sum(clone_size))
    max_editnum = length(pattern)
    ind_allsite_sum_t$pop = ind_allsite_sum_t$num/max_editnum
    
    pt = ggplot()+geom_ribbon(data = scarrefv1,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
      geom_line(data = ind_allsite_sum_t,aes(x = pos,y = pop,color = type),size=1) +
      scale_x_continuous(breaks=c()) +
      xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
    pt
    
    ind_allsite_sum = ind_allsite_per[which(!ind_allsite_per$array %in% bl$Var1),] %>% 
      group_by(pos,type) %>% summarise(num = sum(clone_size))
    max_editnum = length(pattern[which(!pattern %in% bl$Var1)])
    ind_allsite_sum$pop = ind_allsite_sum$num/max_editnum
    
    pf = ggplot()+geom_ribbon(data = scarrefv1,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
      geom_line(data = ind_allsite_sum,aes(x = pos,y = pop,color = type),size=1) +
      scale_x_continuous(breaks=c()) +
      xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
    pf
    
    return(list("ribbon" = pt, "ribbon_fil" = pf,"ind_allsite_per" = ind_allsite_per))
  }
  
  ribpl = list()
  ribplf = list()
  for (i in 1:length(v1ls)) {
    v1rib = INDELRibbon(v1ls[[i]]$pattern,bl$v1_array,scarrefv1)
    v2rib = INDELRibbon(v2ls[[i]]$pattern,bl$v2_array,scarrefv2)
    # ind_allsite_sum_t = v1rib$ind_allsite_per %>% 
    #   group_by(pos,type) %>% summarise(num = sum(clone_size))
    # max_editnum = length(v1ls[[i]]$pattern)
    # ind_allsite_sum_t$pop = ind_allsite_sum_t$num/max_editnum
    # pt1 = ggplot()+geom_ribbon(data = scarrefv1,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
    #   geom_line(data = ind_allsite_sum_t,aes(x = pos,y = pop,color = type),size=1) +
    #   scale_x_continuous(breaks=c()) +
    #   xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
    # 
    # ind_allsite_sum_t = v2rib$ind_allsite_per %>% 
    #   group_by(pos,type) %>% summarise(num = sum(clone_size))
    # max_editnum = length(v2ls[[i]]$pattern)
    # pt2 = ggplot()+geom_ribbon(data = scarrefv2,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
    #   geom_line(data = ind_allsite_sum_t,aes(x = pos,y = num/max_editnum,color = type),size=1) +
    #   scale_x_continuous(breaks=c()) +
    #   xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
    ribpl[[i]] = v1rib$ribbon + ggtitle(paste0(names(v1ls)[i]," V1")) + v2rib$ribbon + ggtitle(paste0(names(v1ls)[i]," V2"))
    ribplf[[i]] = v1rib$ribbon_fil + ggtitle(paste0(names(v1ls)[i]," V1"))  + v2rib$ribbon_fil + ggtitle(paste0(names(v1ls)[i]," V2"))
  }
  
  ggexport(ribpl,filename = "final_result/xtotal/segment_indel_cutsite_stat.pdf",width = 8,height = 4)
  ggexport(ribplf,filename = "final_result/xtotal/segment_indel_cutsite_stat_blfiltered.pdf",width = 8,height = 4)
  
  
  
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
  
  segrefv1 = read.table("raw_data/ref/Cutsite CREST/V1.cutSites")
  segrefv1 = segrefv1[-1,]
  scarrefv1 = ReadCutsite(segrefv1)
  
  segrefv2 = read.table("raw_data/ref/Cutsite CREST/V2.cutSites")
  segrefv2 = segrefv2[-1,]
  scarrefv2 = ReadCutsite(segrefv2)
  
  
  
}

#2.benchmark----------
{
  
  dir.create("final_result/xtotal/benchmark")
  library(Seurat)
  e11tree = readRDS("final_result/xE11_noDoublet/E11_12_21_blay_tree.rds")
  e11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  
  e15tree = readRDS("final_result/xE15_noDoublet/E15_blay_tree.rds")
  e15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  
  #1. DCLEAR
  .libPaths(c("/home/liuhengxin/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library",
              "/picb/sysgenomics2/people/liuhengxin/software/anaconda3/envs/r4-base/lib/R/library/",
              "/picb/sysgenomics2/people/liuhengxin/software/miniconda3/envs/r4.1/lib/R/library/",
              "/picb/sysgenomics2/people/liuhengxin/software/miniconda3/envs/r4.1.2/lib/R/library/"))
  
  .libPaths(c("/picb/sysgenomics2/people/liuhengxin/software/miniconda3/envs/r4.1.2/lib/R/library/"))
  #Analysis of increasing number of tree clusters with glu nb nbglu cor
  CalSpatialCor2 = function(cmmxst,hubcell){
    library(amap)
    
    cmmxst = as.matrix(cmmxst)
    cmmxst = log(cmmxst+1)
    cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
    
    dd = Dist(t(cmmxst), method = "abspearson")
    mymx = as.matrix(1 - dd)
    for (i in 1:nrow(mymx)) {
      mymx[i,i] = 1
    }
    dd[is.na(dd)] = 0
    
    mymx[upper.tri(mymx)] = NA
    mymxl = melt(mymx)
    mymxl = mymxl[!is.na(mymxl$value),];mymxl = mymxl[which(mymxl$Var1 != mymxl$Var2),]
    
    mymxi = melt(mymx[hubcell,hubcell])
    mymxi = mymxi[!is.na(mymxi$value),];mymxi = mymxi[which(mymxi$Var1 != mymxi$Var2),]
    
    mymxi$Var1 = as.character(mymxi$Var1);mymxi$Var2 = as.character(mymxi$Var2)
    mymxl$Var1 = as.character(mymxl$Var1);mymxl$Var2 = as.character(mymxl$Var2)
    mymxi$other = 0
    for (i in 1:nrow(mymxi)) {
      line = mymxi[i,]
      mymxr = mymxl[which((mymxl$Var1 == line$Var1 & mymxl$Var2 != line$Var2)|
                            (mymxl$Var1 != line$Var1 & mymxl$Var2 == line$Var2)),]
      mymxi$other[i] = mean(mymxr$value)
    }
    
    spatcor = mymxi
    spatcor$cor = spatcor$value/mean(spatcor$other)
    spatcor$group = paste0(spatcor$Var1,"-",spatcor$Var2)
    return(spatcor)
  }
  
  #Analysis of increasing number of tree clusters
  tree = dctree$hclust
  hubcell = c("GluFP","NbFP","NbGlu","DA","NbDA")
  SplitHclustTreeCor = function(tree,cm,hubcell){
    spatcor.dlt = NULL
    maxclone = length(tree$labels)
    for (n in c(seq(50,as.integer(maxclone/50)*50,50),maxclone)) {
      node.cluster <- as.data.frame(cutree(tree, k=n))
      colnames(node.cluster) = "group"
      
      clonegroup = node.cluster
      colnames(clonegroup) = "group"
      cmmx = merge(cm,clonegroup, by.x = "tags",by.y = "row.names")
      
      cmmxst = dcast(cmmx,group~celltype);
      rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
      
      spatcor.dl = CalSpatialCor2(cmmxst,hubcell)
      spatcor.dl$sample = n
      spatcor.dlt = rbind(spatcor.dlt,spatcor.dl)
    }
    return(spatcor.dlt)
  }
  
  #Analysis of rep correlation
  treels = castreels2
  tree =treels[[1]]
  cm = datacmbls[[1]]
  SplitHclustTreeRepCor = function(treels,datacmbls,nsample){
    repcor = NULL
    maxclone = min(length(treels[[1]]$labels),length(treels[[2]]$labels),length(treels[[3]]$labels))
    CalCormx = function(tree,cm,n){
      node.cluster <- as.data.frame(cutree(tree, k=n))
      colnames(node.cluster) = "group"
      cmmx = merge(cm,node.cluster, by.x = "tags",by.y = "row.names")
      
      cmmxst = dcast(cmmx,group~celltype);
      rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
      cmmxst = as.matrix(cmmxst)
      cmmxst = log(cmmxst+1)
      cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
      
      dd = Dist(t(cmmxst), method = "abspearson")
      mymx = as.matrix(1 - dd)
      for (i in 1:nrow(mymx)) {
        mymx[i,i] = 1
      }
      return(mymx)
    }
    
    #c(seq(50,as.integer(maxclone/50)*50,50))
    for (n in nsample) {
      mx1 = CalCormx(treels[[1]],datacmbls[[1]],n)
      mx2 = CalCormx(treels[[2]],datacmbls[[2]],n)
      mx3 = CalCormx(treels[[3]],datacmbls[[3]],n)
      
      mx1 = melt(mx1)
      mx2 = melt(mx2)
      mx3 = melt(mx3)
      
      mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),list(mx1,mx2,mx3))
      colnames(mxt)[3:5] = c("rep1","rep2","rep3")
      mxt[is.na(mxt)] = 0
      mxt = mxt[which(mxt$Var1 != mxt$Var2),]
      # mpg.model = lm(rep3~rep1+rep2,data = mxt)
      
      repcorl = data.frame("group" = c("r1_r2","r1_r3","r2_r3"),
                           "cor" = c(cor(mxt$rep1,mxt$rep2),cor(mxt$rep1,mxt$rep3),cor(mxt$rep2,mxt$rep3)),
                           "sample" = n)
      repcorl$cocor = sqrt( ( (repcorl$cor[1])^2 + (repcorl$cor[2])^2 + (repcorl$cor[3])^2 ) - 
                              ( 2 * repcorl$cor[1] * repcorl$cor[2] * repcorl$cor[3]) )
      repcor = rbind(repcor,repcorl)
      
    }
    mx1 = CalCormx(treels[[1]],datacmbls[[1]],length(treels[[1]]$labels))
    mx2 = CalCormx(treels[[2]],datacmbls[[2]],length(treels[[2]]$labels))
    mx3 = CalCormx(treels[[3]],datacmbls[[3]],length(treels[[3]]$labels))
    mx1 = melt(mx1)
    mx2 = melt(mx2)
    mx3 = melt(mx3)
    mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),list(mx1,mx2,mx3))
    colnames(mxt)[3:5] = c("rep1","rep2","rep3")
    mxt[is.na(mxt)] = 0
    mxt = mxt[which(mxt$Var1 != mxt$Var2),]
    
    repcorl = data.frame("group" = c("r1_r2","r1_r3","r2_r3"),
                         "cor" = c(cor(mxt$rep1,mxt$rep2),cor(mxt$rep1,mxt$rep3),cor(mxt$rep2,mxt$rep3)),
                         "sample" = "single")
    repcorl$cocor = sqrt( ( (repcorl$cor[1])^2 + (repcorl$cor[2])^2 + (repcorl$cor[3])^2 ) - 
                            ( 2 * repcorl$cor[1] * repcorl$cor[2] * repcorl$cor[3]) )
    repcor = rbind(repcor,repcorl)
    return(repcor)
    
  }
  
  nsample = c(1:10,seq(50,800,50))
  #DCLEAR analysis
  {
    dctree1 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_rep1_tree.rds")
    dctree2 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_rep2_tree.rds")
    dctree3 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_rep3_tree.rds")
    dctreels = list(dctree1$hclust,dctree2$hclust,dctree3$hclust)
    dctreels2 = list(dctree1$hclust2,dctree2$hclust2,dctree3$hclust2)
    datacmbls = list(e15treels$r1$datacmb,e15treels$r2$datacmb,e15treels$r3$datacmb)
    library(ggtree)
    plot(dctree1$NJ)
    
    pdc = list()
    for (i in 1:3) {
      dctree = dctreels2[[i]]
      # p1.0 = ggtree(dctree$NJ,size= 0.2) + geom_tree(size = 0.1) + geom_tiplab(size= 0.5,align = T,linesize = 0.1) + 
      #   ggtitle(paste0("DCLEAR NJ tree of E15 rep ",i))
      # n = 50
      dctree$node.cluster <- as.data.frame(cutree(dctree$hclust, k=n))
      colnames(dctree$node.cluster) = "group"
      dctree$color <- dctree$node.cluster[match(dctree$NJ$tip.label, rownames(dctree$node.cluster)), 'group']
      p1.1 = ggtree(dctree$NJ,size =  0.2) + geom_tree(size = 0.1) + geom_tiplab(size = 0.5,color = dctree$color,align = T,
                                                                                 linesize = 0.1) + theme_tree() +
        ggtitle(paste0("DCLEAR NJ tree with 50 clones of E15 rep ",i))
      pdc[[i]] = p1.1
    }
    ggexport(pdc,filename = "final_result/xtotal/benchmark/DCLEAR_tree.pdf",width = 8,height = 12)
    #example
    {
      n = 50
      dctree$node.cluster <- as.data.frame(cutree(dctree$hclust, k=n))
      colnames(dctree$node.cluster) = "group"
      dctree$color <- dctree$node.cluster[match(dctree$NJ$tip.label, rownames(dctree$node.cluster)), 
                                          'group']
      
      p1.1 = ggtree(dctree$NJ,size =  0.2) + geom_tree(size = 0.1) + geom_tiplab(size = 0.5,color = dctree$color,align = T,
                                                                                 linesize = 0.1) + theme_tree() + ggtitle("NJ tree of DCLEAR")
      p1.1
      p1.2 = ggtree(dctree$hclust,size =  0.2) + geom_tree(size = 0.1) + geom_tiplab(size = 0.5,color = dctree$color,align = T,
                                                                                     linesize = 0.1) + 
        theme_tree2() + ggtitle("hclust tree of DCLEAR")
      p1.2
      ggexport(list(p1.1,p1.2),filename = "final_result/xtotal/benchmark/DCLEAR/DCLEAR_tree_e15_plot.pdf",
               width = 8,height = 10)
      
      clonegroup = dctree$node.cluster
      colnames(clonegroup) = "group"
      cmmx = e15tree$cm
      cmmx = merge(cmmx,clonegroup, by.x = "tags",by.y = "row.names")
      
      cmmxst = dcast(cmmx,group~celltype);
      rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
      cmmxst = as.matrix(cmmxst)
      # scarst[which(scarst>0)] = 1
      
      # cmmxst = cmmxst[which(rowSums(cmmxst>0)>1 & rowSums(cmmxst>0)< 0.5*ncol(cmmxst)),]
      cmmxst = log(cmmxst+1)
      cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
      
      dd = Dist(t(cmmxst), method = "abspearson")
      mx = as.matrix(1 - dd)
      for (i in 1:nrow(mx)) {
        mx[i,i] = 1
      }
      dd[is.na(dd)] = 0
      p1 = Heatmap(mx,heatmap_legend_param = list(title = "abspearson"),
                   # clustering_distance_rows = "spearman",
                   col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
                   clustering_distance_rows = function(x){x = dd; return(x)},
                   # clustering_distance_columns = "spearman",-
                   clustering_distance_columns = function(x){x = dd; return(x)},
                   clustering_method_columns = "complete",clustering_method_rows = "complete")
      print(p1)
      
    }
    
    repcor.hm = SplitHclustTreeRepCor(dctreels,datacmbls,nsample)
    repcor.dc = SplitHclustTreeRepCor(dctreels2,datacmbls,nsample)
    
    p2 = ggplot(spatcor.dlt2,aes(x = sample, y = cor,group = group)) + geom_line(aes(color = group),size = 0.3) +
      scale_x_continuous(breaks = c(seq(50,550,50),length(unique(cmmx$tags)))) +
      xlab("tree cluster numbers") + ylab("correlation compared to other pairs") +
      theme_pubr()
    p2
    ggsave(p2,filename = "final_result/xtotal/benchmark/DCLEAR/DCLEAR_cluster_cor.pdf",width = 8,height = 7)
    
    ggplot(repcor.dl,aes(x = sample,y = cor,color = group)) + geom_point() + theme_bw() +
      ggplot(repcor.dl,aes(x = sample,y = meancor)) + geom_point() + theme_bw()
  }
  
  #Cassiopeia analysis
  {
    #Build Cassiopeia data
    e15treels = readRDS("final_result/xE15_replicate/e15rep_total_tree.rds")
    
    INDELFormTransTest = function(datacmb,outpath){
      
      idstat = NULL
      for (n in 1:nrow(datacmb)) {
        #v1
        indel = unlist(strsplit(datacmb[n,"tags"],"-"))[2]
        if(indel!="NA"){
          indels = unlist(strsplit(indel,"_"))
          indels = paste0("v1.",indels)
          idstat = rbind(idstat,data.frame(Cell.BC = datacmb[n,"Var1"],indel = unique(indels)))
        }
        
        #v2
        indel = unlist(strsplit(datacmb[n,"tags"],"-"))[3]
        if(indel!="NA"){
          indels = unlist(strsplit(indel,"_"))
          indels = paste0("v2.",indels)
          idstat = rbind(idstat,data.frame(Cell.BC = datacmb[n,"Var1"],indel = unique(indels)))
        }
      }
      idstatf = idstat %>% group_by(indel) %>% summarise(count = n())
      idstatf$freq = idstatf$count/length(unique(idstat$Cell.BC))
      # idstatf = idstatf[which(idstatf$freq < 0.9),]
      idstatf = idstatf[order(-idstatf$freq),]
      idstatf$id = as.numeric(as.factor(idstatf$indel))
      idstatf[which(idstatf$freq > 0.9 | idstatf$indel == "v1.NONE" | idstatf$indel == "v2.NONE"),"id"] = 0
      # idstatf[which(idstatf$freq > 0.9 | idstatf$indel == "v1.NONE" | idstatf$indel == "v2.NONE"),"id"] = 0
      
      idnew = NULL
      for (n in 1:nrow(datacmb)) {
        #v1
        indel = unlist(strsplit(datacmb[n,"tags"],"-"))[2]
        if(indel=="NA"){
          v1idnew = as.data.frame(rbind(rep(0,9)))
        }else{
          indels = unlist(strsplit(indel,"_"))
          indels = paste0("v1.",indels)
          indels = idstatf$id[match(indels,idstatf$indel)]
          indelsrm = indels
          for (j in 2:length(indels)) {
            if(indels[j] == indels[j-1] & indels[j] != 0) indelsrm[j] = -1
          }
          v1idnew = as.data.frame(rbind(indelsrm))
        }
        #v2
        indel = unlist(strsplit(datacmb[n,"tags"],"-"))[3]
        if(indel=="NA"){
          v2idnew = as.data.frame(rbind(rep(0,8)))
        }else{
          indels = unlist(strsplit(indel,"_"))
          indels = paste0("v2.",indels)
          indels = idstatf$id[match(indels,idstatf$indel)]
          indelsrm = indels
          for (j in 2:length(indels)) {
            if(indels[j] == indels[j-1] & indels[j] != 0) indelsrm[j] = -1
          }
          v2idnew = as.data.frame(rbind(indelsrm))
        }
        idnewline = cbind(v1idnew,v2idnew)
        idnew = rbind(idnew,idnewline)
        
      }
      
      colnames(idnew) = paste0("r",1:ncol(idnew))
      datacmbn = cbind(datacmb[1],idnew)
      colnames(datacmbn)[1] = "cellBC"
      
      priors = NULL
      for (m in 2:ncol(datacmbn)) {
        vline = datacmbn[,m]
        vlinerm = vline
        vlinerm[which(!vlinerm %in% c(0,-1))] = as.numeric(as.factor(vline[which(!vline %in% c(0,-1))]))
        datacmbn[,m] = vlinerm
        iddf = unique(data.frame("oldid" = vline[which(!vline %in% c(0,-1))],
                                 "newid" = vlinerm[which(!vlinerm %in% c(0,-1))] ))
        propor = merge(idstatf,iddf,by.x = "id",by.y = "oldid")
        propor$col = m
        priors = rbind(priors,propor)
        
      }
      
      write.csv(datacmbn,file = paste0(outpath,"allel_table.csv"),quote = F,row.names = F)
      write.csv(priors,
                file = paste0(outpath,"prior.csv"),quote = F,row.names = F)
      
      #stat indel freq
      
      return(list("allel_tab" = datacmbn,"priors" = priors))
    }
    cas_r1df = INDELFormTransTest(e15treels$r1$datacmb,"final_result/CREST_final_analysis/benchmark/Cassiopeia/e15rep1")
    cas_r2df = INDELFormTransTest(e15treels$r2$datacmb,"final_result/CREST_final_analysis/benchmark/Cassiopeia/e15rep2")
    cas_r3df = INDELFormTransTest(e15treels$r3$datacmb,"final_result/CREST_final_analysis/benchmark/Cassiopeia/e15rep3")
    
    #analysis
    castree1 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/Cassiopeia-master/mytest/e15r1_greedy_tree.nex")
    castree2 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/Cassiopeia-master/mytest/e15r2_greedy_tree.nex")
    castree3 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/Cassiopeia-master/mytest/e15r3_greedy_tree.nex")
    
    
    CasToHclust = function(e15r1tree_cas,datacmb){
      edge = as.data.frame(e15r1tree_cas$edge)
      colnames(edge) = c("from","to")
      edgeterm = datacmb[match(e15r1tree_cas$tip.label,datacmb$Var1),"tags"]
      edgeterm = data.frame("cb" = e15r1tree_cas$tip.label,"clone" = edgeterm,"id" = 1:length(edgeterm))
      edgeclone = merge(edge,edgeterm,by.x = "to",by.y = "id")
      edgeclone = unique(edgeclone[c("from","clone")])
      length(unique(edgeclone$clone))
      edgeclone = edgeclone %>% group_by(clone) %>% summarise(from = mean(from))
      edgeclonemx = NULL
      for (i in 1:nrow(edgeclone)) {
        line =  edgeclone
        line$clone2 = edgeclone[i,]$clone
        line$from2 = edgeclone[i,]$from
        line$dist = abs(line$from2 - line$from)
        edgeclonemx = rbind(edgeclonemx, line)
      }
      edgeclonemx = dcast(edgeclonemx,clone~clone2,value.var = "dist")
      rownames(edgeclonemx) = edgeclonemx$clone;edgeclonemx = edgeclonemx[-1]
      hc = hclust(as.dist(edgeclonemx))
      return(hc)
    }
    castreels = list(CasToHclust(castree1,e15treels$r1$datacmb),
                     CasToHclust(castree2,e15treels$r2$datacmb),
                     CasToHclust(castree3,e15treels$r3$datacmb))
    castreels2 = list(castree1,castree2,castree3)
    pcas = list()
    for (i in 1:3) {
      castree = castreels2[[i]]
      # p1.0 = ggtree(dctree$NJ,size= 0.2) + geom_tree(size = 0.1) + geom_tiplab(size= 0.5,align = T,linesize = 0.1) + 
      #   ggtitle(paste0("DCLEAR NJ tree of E15 rep ",i))
      # n = 50
      node.cluster <- as.data.frame(cutree(castreels[[i]], k=n))
      colnames(node.cluster) = "group"
      node.cluster = merge(node.cluster,datacmbls[[i]],by.x = "row.names",by.y = "tags")
      castree$color <- node.cluster[match(castree$tip.label, node.cluster$Var1), 'group']
      p2.1 = ggtree(castree,size =  0.05) + geom_tree(size = 0.05) + geom_tiplab(size = 0.1,color = castree$color,align = T,
                                                                                 linesize = 0.05) + theme_tree() +
        ggtitle(paste0("Cassiopeia greedy tree with 50 clones of E15 rep ",i))
      p2.1
      pcas[[i]] = p2.1
    }
    ggexport(pcas,filename = "final_result/xtotal/benchmark/Cas_tree.pdf",width = 30,height =100)
    
    datacmbls = list(e15treels$r1$datacmb,e15treels$r2$datacmb,e15treels$r3$datacmb)
    repcor.cas = SplitHclustTreeRepCor(castreels,datacmbls,nsample)
    write.csv(repcor.cas,file = "final_result/CREST_final_analysis/benchmark/E15_repcor_cas.csv",quote = F,row.names = F)
    
    ggplot(repcor.cas,aes(x = sample,y = cor,color = group)) + geom_point() + theme_bw()
    
    mycolor = e15r1dfhc$hcgroup[match(e15r1tree_cas$tip.label,e15r1dfhc$Cell.BC)]
    p1.1 = ggtree(e15r1tree_cas,size =  0.2) + geom_tree(size = 0.1) + 
      geom_tiplab(size = 0.1,color = mycolor,align = T,linesize = 0.1) + theme_tree() + ggtitle("tree of Cas")
    p1.1
    
  }
  
  #LinTind tree
  {
    tree_lt = e15treels$r1$vertice
    LinTToHclust = function(tree_lt,cm){
      edgeclone = data.frame("clone" = tree_lt$tags,"id" = 1:nrow(tree_lt))
      edgeclone = edgeclone[edgeclone$clone %in% cm$tags,]
      edgeclonemx = NULL
      for (i in 1:nrow(edgeclone)) {
        line =  edgeclone
        line$clone2 = edgeclone[i,]$clone
        line$id2 = edgeclone[i,]$id
        line$dist = abs(line$id2 - line$id)
        edgeclonemx = rbind(edgeclonemx, line)
      }
      edgeclonemx = dcast(edgeclonemx,clone~clone2,value.var = "dist")
      rownames(edgeclonemx) = edgeclonemx$clone;edgeclonemx = edgeclonemx[-1]
      hc = hclust(as.dist(edgeclonemx))
      return(hc)
    }
    lttreels = list(LinTToHclust(e15treels$r1$vertice,e15treels$r1$cm),
                    LinTToHclust(e15treels$r2$vertice,e15treels$r2$cm),
                    LinTToHclust(e15treels$r3$vertice,e15treels$r3$cm))
    cmls = list(e15treels$r1$cm,e15treels$r2$cm,e15treels$r3$cm)
    repcor.lt = SplitHclustTreeRepCor(lttreels,cmls)
    
    ggplot(repcor.lt,aes(x = sample,y = cor,color = group)) + geom_point() + theme_bw() +
      ggplot(repcor.lt,aes(x = sample,y = cocor)) + geom_point() + theme_bw()
  }
  
  
  treels = list("DCLEAR" = dctreels,"Cassiopeia" = castreels,"LinTind" = lttreels)
  repcor.hm$method = "hamming"
  repcor.dc$method = "DCLEAR"
  repcor.cas$method = "Cassiopeia"
  repcor = rbind(repcor.hm,repcor.dc,repcor.cas)
  write.csv(repcor,file = "final_result/CREST_final_analysis/benchmark/E15_repcor_total.csv")
  
  
  
  #simulate control data
  {
    
    GenerateEdit = function(time,prograte,difrate,editratee,editratel,isdiratee,isdiratel){
      ctsim = list(c(1:5),c(4:8),c(1,9:12))
      progsim = data.frame("celltype" = sample(1:5,200,replace = T),"r0" = 0)
      progsim$celltype = paste0("prog",progsim$celltype)
      # progsim$r1 = 1:nrow(progsim)
      allcell = progsim 
      #proliferation
      for (t in 1:time) {
        print(t)
        allcell = cbind(allcell,0)
        colnames(allcell)[t+2] = paste0("r",t)
        progsim = allcell[which(substr(allcell$celltype,1,4) == "prog"),]
        nerosim = allcell[which(substr(allcell$celltype,1,4) != "prog"),]
        newcell = NULL
        for (i in 1:nrow(progsim)) {
          if(substr(progsim[i,]$celltype,1,4) == "prog"){
            isdiri = runif(1)
            if(t < time/2){isdiratei = isdiratee}else{isdiratei = isdiratel}
            if(isdiri < isdiratei){
              newcell = rbind(newcell,progsim[rep(i,sample(1:prograte,1)),])
            }else{
              celli = progsim[i,]
              progid = as.numeric(substr(celli$celltype,nchar(celli$celltype),
                                         nchar(celli$celltype)))
              if(progid < 4){
                newcellt = sample(ctsim[[progid]], sample(1:difrate,1), replace = T)
                celli = celli[rep(1,length(newcellt)),]
                celli$celltype = paste0("neuro",newcellt)
                newcell = rbind(newcell,celli)
              }
              
            }
          }else{
            newcell = rbind(newcell,progsim[i,])
          }
          
        }
        ifedit = runif(nrow(newcell))
        
        if(t < time/2){editrate =editratee}else{editrate = editratel}
        newcell[ifedit < editrate,t+2] = 1:nrow(newcell[ifedit < editrate,])
        allcell = rbind(nerosim,newcell)
      }
      allcell = allcell[-2]
      return(allcell)
    }
    simls1 = list()
    simls2 = list()
    for (i in 1:100) {
      print(paste0("i:",i))
      simls1[[i]] = GenerateEdit(8,3,1,0.5,0.8,0.8,0.3)
      simls2[[i]] = GenerateEdit(8,3,3,0.8,0.05,0.8,0.3)
    }
    library(qs)
    qsave(list("sim1" = simls1,"sim2" = simls2),file = "final_result/CREST_final_analysis/benchmark/simdatals_2_4.qs")
    sim1 = GenerateEdit(8,3,1,0.5,0.8,0.8,0.3)
    sim2 = GenerateEdit(8,3,1,0.5,0.8,0.8,0.3)
    sim3 = GenerateEdit(8,3,1,0.5,0.8,0.8,0.3)
    
    sim21 = GenerateEdit(8,3,3,0.8,0.05,0.8,0.3)
    sim22 = GenerateEdit(8,3,3,0.8,0.05,0.8,0.3)
    sim23 = GenerateEdit(8,3,3,0.8,0.05,0.8,0.3)
    
    #sim analysis
    dir.create("final_result/CREST_final_analysis/benchmark/simulate")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim1")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim2")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim1/cas")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim2/cas")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim1/hamming")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim2/hamming")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim1/dclear")
    dir.create("final_result/CREST_final_analysis/benchmark/simulate/sim2/dclear")
    
    outfile = "final_result/CREST_final_analysis/benchmark/simulate/sim1"
    SimTreeBuild = function(simls1,outfile)
    {
      library(reshape2)
      library(dplyr)
      simdclls = list()
      hmtreels = list()
      for (i in 1:length(simls1)) {
        print(i)
        #tree build
        simi = simls1[[i]]
        
        #hamming tree
        hmtreehc = BuildHammingdata(simi)
        hmtreels[[i]] = hmtreehc
        
        #cas
        simi$celltype = paste0(simi$celltype,"_",1:nrow(simi))
        simil = melt(simi)
        priori = simil %>% group_by(variable,value) %>% summarise(count = length(unique(celltype)),
                                                                  freq = count/nrow(simi))
        priori = priori[which(priori$value != 0),]
        colnames(priori)[1:2] = c("col","id")
        priori$col = as.numeric(substr(priori$col,2,nchar(as.character(priori$col))))
        write.csv(simi,file = paste0(outfile,"/simu2_allel_table_",i,".csv"),quote = F,row.names = F)
        write.csv(priori,file = paste0(outfile,"/simu2_priors_",i,".csv"),quote = F,row.names = F)
        
        #dclear
        simscari = simls1[[i]]
        simscari$tags = apply(simls1[[i]][-1], 1, function(x){return(paste0(sprintf("%05d",x),collapse = ""))})
        simscari = simscari[c(1,ncol(simscari))]
        simdclls[[i]] = list("scardata" = simscari,"states" = unique(unlist(as.matrix(simi[-1]))))
      }
      qsave(hmtreels,file = paste0(outfile,"/hamming/hmtreels.qs"))
      qsave(simdclls,file = paste0(outfile,"/dclear/simdclls.qs"))
    }
    SimTreeBuild(simls1,"final_result/CREST_final_analysis/benchmark/simulate/sim1")
    SimTreeBuild(simls2,"final_result/CREST_final_analysis/benchmark/simulate/sim2")
    
    #simu tree analysis
    {
      #load data
      hamtreels1 = qread("final_result/CREST_final_analysis/benchmark/simulate/sim1/hamming/hmtreels.qs")
      hamtreels2 = qread("final_result/CREST_final_analysis/benchmark/simulate/sim2/hamming/hmtreels.qs")
      dcltreels1 = qread("final_result/CREST_final_analysis/benchmark/simulate/sim1/dclear/sim1dctree.qs")
      dcltreels2 = qread("final_result/CREST_final_analysis/benchmark/simulate/sim2/dclear/sim2dctree.qs")
      
      castreels1 = lapply(list.files("final_result/CREST_final_analysis/benchmark/simulate/sim1/cas/",pattern = "*.nex",full.names = T),
                         read.tree)
      cellclls1 = lapply(list.files("final_result/CREST_final_analysis/benchmark/simulate/sim1/",pattern = "simu2_allel_table_*",full.names = T),
                         read.csv)
      for (i in 1:length(cellclls1)) {
        colnames(cellclls1[[i]])[1] = "Var1"
        cellclls1[[i]]$tags = apply(cellclls1[[i]][-1], 1, function(x){return(paste0(x,collapse = "_"))})
        cellclls1[[i]]$celltype = unlist(lapply(strsplit(cellclls1[[i]]$Var1,"_"),"[[",1))
      }
      
      castreels2 = lapply(list.files("final_result/CREST_final_analysis/benchmark/simulate/sim2/cas/",pattern = "*.nex",
                                     full.names = T),read.tree)
      cellclls2 = lapply(list.files("final_result/CREST_final_analysis/benchmark/simulate/sim2/",pattern = "simu2_allel_table_*",full.names = T),
                         read.csv)
      for (i in 1:length(cellclls2)) {
        colnames(cellclls2[[i]])[1] = "Var1"
        cellclls2[[i]]$tags = apply(cellclls2[[i]][-1], 1, function(x){return(paste0(x,collapse = "_"))})
        cellclls2[[i]]$celltype = unlist(lapply(strsplit(cellclls2[[i]]$Var1,"_"),"[[",1))
      }
      
      #tree analysis
      nsample = c(1,2,3,4,5,10,50,100,300,500,600,700,800)
      # nsample = c(1,2,3,4,5,10,50,100,300,500,600,700,800,"single")
      
      #hm
      SplitHclustTreeRepCorHm = function(treels,nsample){
        repcort = NULL
        # maxclone = min(unlist(lapply(treels, function(x){return(length(x$labels))})))
        CalCormx2 = function(treei,n){
          # ggtree(treei) + geom_tiplab(size= 0.5)
          if(length(treei$labels) >= n){
            node.cluster = data.frame(cutree(treei, k=n),treei$labels)
          }else{
            node.cluster = data.frame(cutree(treei, k=length(treei$labels)),treei$labels)
          }
          
          colnames(node.cluster) = c("group","celltype")
          
          cmmxst = dcast(node.cluster,group~celltype,fun.aggregate = length)
          rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
          cmmxst = as.matrix(cmmxst)
          cmmxst = log(cmmxst+1)
          cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
          
          dd = Dist(t(cmmxst), method = "abspearson")
          mymx = as.matrix(1 - dd)
          for (i in 1:nrow(mymx)) {
            mymx[i,i] = 1
          }
          return(mymx)
        }
        #c(seq(50,as.integer(maxclone/50)*50,50))
        for (n in nsample) {
          mxls = list()
          for (i in 1:length(treels)) {
            mxi = CalCormx2(treels[[i]],n)
            mxls[[i]] = melt(mxi)
          }
          mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),mxls)
          colnames(mxt)[3:ncol(mxt)] = paste0("rep",1:length(treels))
          mxt[is.na(mxt)] = 0
          mxt = mxt[which(mxt$Var1 != mxt$Var2),]
          # mpg.model = lm(rep3~rep1+rep2,data = mxt)
          
          repcor = cor(as.matrix(mxt[-c(1,2)]))
          repcorl = as.data.frame(melt(repcor))
          repcorl = repcorl[which(repcorl$Var1 != repcorl$Var2),]
          colnames(repcorl)[3] = "cor"
          repcorl$group = paste0(repcorl$Var1, "_", repcorl$Var2)
          repcorl$sample = n
          # repcorl = data.frame("group" = c("r1_r2","r1_r3","r2_r3"),
          #                      "cor" = c(cor(mxt$rep1,mxt$rep2),cor(mxt$rep1,mxt$rep3),cor(mxt$rep2,mxt$rep3)),
          #                      "sample" = n)
          # repcorl$cocor = sqrt( ( (repcorl$cor[1])^2 + (repcorl$cor[2])^2 + (repcorl$cor[3])^2 ) - 
          #                         ( 2 * repcorl$cor[1] * repcorl$cor[2] * repcorl$cor[3]) )
          repcort = rbind(repcort,repcorl)
          
        }
        return(repcort)
      }
      repcorhm1 = SplitHclustTreeRepCorHm(hamtreels1, nsample)
      repcorhm2 = SplitHclustTreeRepCorHm(hamtreels2, nsample)
      
      #cas
      CasToHclust2 = function(castree1,cellcl1){
        edge = as.data.frame(castree1$edge)
        colnames(edge) = c("from","to")
        edgeterm = cellcl1[match(castree1$tip.label,cellcl1$Var1),"tags"]
        edgeterm = data.frame("cb" = castree1$tip.label,"clone" = edgeterm,"id" = 1:length(edgeterm))
        edgeclone = merge(edge,edgeterm,by.x = "to",by.y = "id")
        edgeclone = unique(edgeclone[c("from","clone")])
        length(unique(edgeclone$clone))
        edgeclone = edgeclone %>% group_by(clone) %>% summarise(from = mean(from))
        edgeclonetmp = edgeclone$from
        names(edgeclonetmp) = edgeclone$clone
        edgeclonemxs = dist(edgeclonetmp)
        hc = hclust(edgeclonemxs)
        return(hc)
      }
      simcashcls1 = list()
      simcashcls2 = list()
      for (i in 1:length(castreels1)) {
        print(i)
        simcashcls1[[i]] = CasToHclust2(castreels1[[i]],cellclls1[[i]])
        simcashcls2[[i]] = CasToHclust2(castreels2[[i]],cellclls2[[i]])
      }
      SplitHclustTreeRepCor = function(treels,datacmbls,nsample){
        repcort = NULL
        # maxclone = min(length(treels[[1]]$labels),length(treels[[2]]$labels),length(treels[[3]]$labels))
        CalCormx = function(tree,cm,n){
          if(length(treels[[i]]$labels) >= n){
            node.cluster = as.data.frame(cutree(tree, k=n))
          }else{
            node.cluster = as.data.frame(cutree(tree, k=length(treels[[i]]$labels)))
          }
          colnames(node.cluster) = "group"
          cmmx = merge(cm,node.cluster, by.x = "tags",by.y = "row.names")
          
          cmmxst = dcast(cmmx,group~celltype);
          rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
          cmmxst = as.matrix(cmmxst)
          cmmxst = log(cmmxst+1)
          cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
          
          dd = Dist(t(cmmxst), method = "abspearson")
          mymx = as.matrix(1 - dd)
          for (i in 1:nrow(mymx)) {
            mymx[i,i] = 1
          }
          return(mymx)
        }
        
        #c(seq(50,as.integer(maxclone/50)*50,50))
        for (n in nsample) {
          mxls = list()
          for (i in 1:length(treels)) {
            mxi = CalCormx(treels[[i]],datacmbls[[i]],n)
            mxls[[i]] = melt(mxi)
          }
          mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),mxls)
          colnames(mxt)[3:ncol(mxt)] = paste0("rep",1:length(mxls))
          mxt[is.na(mxt)] = 0
          mxt = mxt[which(mxt$Var1 != mxt$Var2),]
          
          repcor = cor(as.matrix(mxt[-c(1,2)]))
          repcorl = as.data.frame(melt(repcor))
          repcorl = repcorl[which(repcorl$Var1 != repcorl$Var2),]
          colnames(repcorl)[3] = "cor"
          repcorl$group = paste0(repcorl$Var1, "_", repcorl$Var2)
          repcorl$sample = n
          
          repcort = rbind(repcort,repcorl)
          
        }
        
        mxls = list()
        for (i in 1:length(treels)) {
          mxi = CalCormx(treels[[i]],datacmbls[[i]],length(treels[[i]]$labels))
          mxls[[i]] = melt(mxi)
        }
        mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),mxls)
        colnames(mxt)[3:ncol(mxt)] = paste0("rep",1:length(mxls))
        mxt[is.na(mxt)] = 0
        mxt = mxt[which(mxt$Var1 != mxt$Var2),]
        
        repcor = cor(as.matrix(mxt[-c(1,2)]))
        repcorl = as.data.frame(melt(repcor))
        repcorl = repcorl[which(repcorl$Var1 != repcorl$Var2),]
        colnames(repcorl)[3] = "cor"
        repcorl$group = paste0(repcorl$Var1, "_", repcorl$Var2)
        repcorl$sample = "single"
        repcort = rbind(repcort,repcorl)
        
        return(repcort)
        
      }
      repcorcas1 = SplitHclustTreeRepCor(simcashcls1, cellclls1, nsample)
      repcorcas2 = SplitHclustTreeRepCor(simcashcls2, cellclls2, nsample)
      
      #dclear
      sim1dc = qread("final_result/CREST_final_analysis/benchmark/simulate/sim1/dclear/simdclls.qs")
      sim2dc = qread("final_result/CREST_final_analysis/benchmark/simulate/sim2/dclear/simdclls2.qs")
      dctreesimls1 = list()
      dctreesimls2 = list()
      dclcellls1 = list()
      dclcellls2 = list()
      hmtreesimls1 = list()
      hmtreesimls2 = list()
      
      for (i in 1:length(dcltreels1)) {
        dctreesimls1[[i]] = dcltreels1[[i]]$hclust2
        dctreesimls2[[i]] = dcltreels2[[i]]$hclust2
        hmtreesimls1[[i]] = dcltreels1[[i]]$hclust
        hmtreesimls2[[i]] = dcltreels2[[i]]$hclust
        dclcellls1[[i]] = sim1dc[[i]]$scardata
        dclcellls2[[i]] = sim2dc[[i]]$scardata
      }
      repcorhm1 = SplitHclustTreeRepCor(hmtreesimls1,dclcellls1,nsample)
      repcorhm2 = SplitHclustTreeRepCor(hmtreesimls2,dclcellls2,nsample)
      repcordc1 = SplitHclustTreeRepCor(dctreesimls1,dclcellls1,nsample)
      repcordc2 = SplitHclustTreeRepCor(dctreesimls2,dclcellls2,nsample)
      
      #plot
      repcorcas1$method = "Cassiopeia"
      repcorhm1$method = "hamming"
      repcordc1$method = "DCLEAR"
      repcorsim1 = rbind(repcorcas1,repcorhm1,repcordc1)
      
      repcorcas2$method = "Cassiopeia"
      repcorhm2$method = "hamming"
      repcordc2$method = "DCLEAR"
      repcorsim2 = rbind(repcorcas2,repcorhm2,repcordc2)
      
      repcore15 = read.csv(file = "final_result/CREST_final_analysis/benchmark/E15_repcor_total.csv")
      repcore15 = repcore15[c("sample","group","cor","method")]
      repcorsim1 = repcorsim1[c("sample","group","cor","method")]
      repcorsim2 = repcorsim2[c("sample","group","cor","method")]
      colnames(repcore15)[1] = colnames(repcorsim1)[1] = colnames(repcorsim2)[1] = "clonenum"
      repcore15$sample = "E15.5(n = 3)"
      repcorsim1$sample = "simulate2(dynamic barcoding, n = 100)"
      repcorsim2$sample = "simulate1(static barcoding, n = 100)"
      
      repcort = rbind(repcore15,repcorsim1,repcorsim2)
      
      library(ggsci)
      library(ggpubr)
      repcort[which(repcort$clonenum == "single"),"method"] = "single clone"
      repcort$clonenum = factor(repcort$clonenum,levels = unique(repcort$clonenum))
      psample = c(1,2,3,4,5,10,50,100,300,500,600,700,800,"single")
      
      repcort.plot = unique(repcort[which(repcort$clonenum %in% psample),])
      repcort.plot = repcort.plot[!duplicated(repcort.plot[c("clonenum","sample","method","group")]),]
      
      p4.0 = ggbarplot(repcort.plot,
                x = "clonenum", y = "cor", 
                fill = "method",size = 0.1,
                width = 0.7,
                position = position_dodge(width = 0.9),
                add = c("mean_sd"),add.params = list(size = 0.2),
                facet.by = "sample",ncol = 1) + 
        # scale_y_continuous(limits = c(0,1)) +
        scale_fill_npg() +
        theme_pubr() + 
        xlab("Clone number") + 
        ylab("Correlation coefficient among replicates")
      p4.0
      ggsave(p4.0,filename = "final_result/CREST_final_analysis/benchmark/benchmark_tree_repcor_diffclonenumber_withsim_2_6.pdf",width = 7,height = 5)
      qsave(repcort,file = "final_result/CREST_final_analysis/benchmark/benchmark_tree_repcor_diffclonenumber_withsim_2_6.qs")
      
    }
    
    
    #hamming tree test
    {
      library(phangorn)
      BuildHammingdata = function(simi){
        scar = apply(simi[-1], 1, function(x){return(paste0(sprintf("%05d",x),collapse = ""))})
        df = data.frame("cellcb" = simi$celltype,"edit" = scar)
        tip_names <- df[, 1]
        states <- unique(unlist(as.matrix(simi[-1])))
        sequence_length <- 50
        # states2num = 1:length(states)
        # names(states2num) = states
        df <- do.call('rbind', strsplit(as.character(df[, 2]), ''))
        rownames(df) <- tip_names
        
        df <- df %>% phyDat(type = 'USER', levels = states)
        D_h = df %>% dist.hamming()
        # D_h = df %>% dist_replacement(reps = 20L, k = 2L)
        # hmtree = NJ(D_h)
        hmtreehc = hclust(D_h)
        remove(df)
        # qsave(hmdatals,"final_result/CREST_final_analysis/benchmark/hmtreehc.qs")
        return(hmtreehc)
      }
      
      hmdatals
      hmtreehc1 = BuildHammingdata(sim21)
      hmtreehc2 = BuildHammingdata(sim22)
      hmtreehc3 = BuildHammingdata(sim23)
      treels = list(hmtreehc1,hmtreehc2,hmtreehc3)
      
      # nsample = c(1,50,100,500,750,1000,1200)
      # nsample = c(5,10,20,30,40,50,100,200,300,500,600,800)
      SplitHclustTreeRepCorHm = function(treels,nsample){
        repcor = NULL
        maxclone = min(length(treels[[1]]$labels),length(treels[[2]]$labels),length(treels[[3]]$labels))
        CalCormx2 = function(treei,n){
          ggtree(treei) + geom_tiplab(size= 0.5)
          node.cluster <- data.frame(cutree(treei, k=n),treei$labels)
          colnames(node.cluster) = c("group","celltype")
          
          cmmxst = dcast(node.cluster,group~celltype);
          rownames(cmmxst) = cmmxst$group;cmmxst = cmmxst[-1]
          cmmxst = as.matrix(cmmxst)
          cmmxst = log(cmmxst+1)
          cmmxst = t(apply(cmmxst, 1, function(x) {x/sum(x)}))
          
          dd = Dist(t(cmmxst), method = "abspearson")
          mymx = as.matrix(1 - dd)
          for (i in 1:nrow(mymx)) {
            mymx[i,i] = 1
          }
          return(mymx)
        }
        
        #c(seq(50,as.integer(maxclone/50)*50,50))
        for (n in nsample) {
          
          mx1 = CalCormx2(treels[[1]],n)
          mx2 = CalCormx2(treels[[2]],n)
          mx3 = CalCormx2(treels[[3]],n)
          
          mx1 = melt(mx1)
          mx2 = melt(mx2)
          mx3 = melt(mx3)
          
          mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),list(mx1,mx2,mx3))
          colnames(mxt)[3:5] = c("rep1","rep2","rep3")
          mxt[is.na(mxt)] = 0
          mxt = mxt[which(mxt$Var1 != mxt$Var2),]
          # mpg.model = lm(rep3~rep1+rep2,data = mxt)
          
          repcorl = data.frame("group" = c("r1_r2","r1_r3","r2_r3"),
                               "cor" = c(cor(mxt$rep1,mxt$rep2),cor(mxt$rep1,mxt$rep3),cor(mxt$rep2,mxt$rep3)),
                               "sample" = n)
          repcorl$cocor = sqrt( ( (repcorl$cor[1])^2 + (repcorl$cor[2])^2 + (repcorl$cor[3])^2 ) - 
                                  ( 2 * repcorl$cor[1] * repcorl$cor[2] * repcorl$cor[3]) )
          repcor = rbind(repcor,repcorl)
          
        }
        return(repcor)
        
      }
      repcor = SplitHclustTreeRepCorHm(treels,c(1:10,seq(50,1600,50)))
      repcor$sample = factor(repcor$sample,levels = unique(repcor$sample))
      ggplot(repcor,aes(x = sample,y = cor, color = group)) + 
        geom_bar(stat = "identity",position = "dodge",fill = "white") + 
        theme_pubr() + 
        scale_color_lancet()+
        xlab("log10(clone number)")
      
      qsave(list("rep1" = sim1,"rep2" = sim2,"rep3" = sim3),
            file = "final_result/CREST_final_analysis/benchmark/simdatals.qs")
      qsave(list("rep1" = hmtreehc1,"rep2" = hmtreehc2,"rep3" = hmtreehc3),
            file = "final_result/CREST_final_analysis/benchmark/hammingtree.qs")
    }
    
    #cas analysis
    {
      simls = list(sim21,sim22,sim23)
      for (i in 1:3) {
        simls[[i]]$celltype = paste0(simls[[i]]$celltype,"_",1:nrow(simls[[i]]))
      }
      priorls = list()
      for (i in 1:length(simls)) {
        priori = NULL
        simi = simls[[i]]
        simil = melt(simi)
        priori = simil %>% group_by(variable,value) %>% summarise(count = length(unique(celltype)),
                                                                  freq = count/nrow(simi))
        priori = priori[which(priori$value != 0),]
        colnames(priori)[1:2] = c("col","id")
        priori$col = as.numeric(substr(priori$col,2,nchar(as.character(priori$col))))
        priorls[[i]] = priori
      }
      
      for (i in 1:3) {
        simls[[i]] = simls[[i]]
        write.csv(simls[[i]],file = paste0("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu2_allel_table_",i,".csv"),
                  quote = F,row.names = F)
        write.csv(priorls[[i]],file = paste0("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu2_priors_",i,".csv"),
                  quote = F,row.names = F)
      }
      
      #simulate cas data analysis
      {
        library(ggtree)
        library(data.tree)
        castree1 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_1.nex")
        cellcl1 = read.csv("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_allel_table_1.csv")
        colnames(cellcl1)[1] = "Var1"
        cellcl1$tags = apply(cellcl1[-1], 1, function(x){return(paste0(x,collapse = "_"))})
        cellcl1$celltype = unlist(lapply(strsplit(cellcl1$Var1,"_"),"[[",1))
        
        castree2 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_2.nex")
        cellcl2 = read.csv("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_allel_table_2.csv")
        colnames(cellcl2)[1] = "Var1"
        cellcl2$tags = apply(cellcl2[-1], 1, function(x){return(paste0(x,collapse = "_"))})
        cellcl2$celltype = unlist(lapply(strsplit(cellcl2$Var1,"_"),"[[",1))
        
        castree3 = read.tree("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_3.nex")
        cellcl3 = read.csv("final_result/CREST_final_analysis/benchmark/Cassiopeia/simu_allel_table_3.csv")
        colnames(cellcl3)[1] = "Var1"
        cellcl3$tags = apply(cellcl3[-1], 1, function(x){return(paste0(x,collapse = "_"))})
        cellcl3$celltype = unlist(lapply(strsplit(cellcl3$Var1,"_"),"[[",1))
        
        #cas tree
        CasToHclust2 = function(castree1,cellcl1){
          # colnames(cellcl1)[1] = "Cell.BC"
          # cellcl1$clone = apply(cellcl1[-1], 1, function(x){return(paste0(x,collapse = "_"))})
          # cellcl1$celltype = unlist(lapply(strsplit(cellcl1$Cell.BC,"_"),"[[",1))
          # 
          edge = as.data.frame(castree1$edge)
          colnames(edge) = c("from","to")
          edgeterm = cellcl1[match(castree1$tip.label,cellcl1$Var1),"tags"]
          edgeterm = data.frame("cb" = castree1$tip.label,"clone" = edgeterm,"id" = 1:length(edgeterm))
          edgeclone = merge(edge,edgeterm,by.x = "to",by.y = "id")
          edgeclone = unique(edgeclone[c("from","clone")])
          length(unique(edgeclone$clone))
          edgeclone = edgeclone %>% group_by(clone) %>% summarise(from = mean(from))
          edgeclonetmp = edgeclone$from
          names(edgeclonetmp) = edgeclone$clone
          edgeclonemxs = dist(edgeclonetmp)
          hc = hclust(edgeclonemxs)
          # edgeclonemx = NULL
          # for (i in 1:nrow(edgeclone)) {
          #   print(i)
          #   line =  edgeclone
          #   line$clone2 = edgeclone[i,]$clone
          #   line$from2 = edgeclone[i,]$from
          #   line$dist = abs(line$from2 - line$from)
          #   edgeclonemx = rbind(edgeclonemx, line)
          # }
          # edgeclonemxs = dcast(edgeclonemx,clone~clone2,value.var = "dist")
          # rownames(edgeclonemxs) = edgeclonemxs$clone;edgeclonemxs = edgeclonemxs[-1]
          # hc = hclust(as.dist(edgeclonemxs))
          return(hc)
        }
        simcashc1 = CasToHclust2(castree1,cellcl1)
        simcashc2 = CasToHclust2(castree2,cellcl2)
        simcashc3 = CasToHclust2(castree3,cellcl3)
        qsave(list(simcashc1,simcashc2,simcashc3),file = "final_result/CREST_final_analysis/benchmark/Cassiopeia/simcashcls.qs")
        
        castreels1 = list(castree1,castree2,castree3)
        castreels2 = list(simcashc1,simcashc2,simcashc3)
        cellclls = list(cellcl1,cellcl2,cellcl3)
        pcas = list()
        ptmp = ggtree(simcashc1) + geom_tiplab(size= 0.5)
        ggsave(ptmp,filename = "final_result/CREST_final_analysis/benchmark/Cassiopeia/test_tree.pdf",width = 10,height = 20)
        
        for (i in 1:3) {
          castree = castreels1[[i]]
          node.cluster <- as.data.frame(cutree(castreels2[[i]], k=n))
          colnames(node.cluster) = "group"
          node.cluster = merge(node.cluster,cellclls[[i]],by.x = "row.names",by.y = "tags")
          castree$color <- node.cluster[match(castree$tip.label, node.cluster$Var1), 'group']
          p2.1 = ggtree(castree,size =  0.05) + geom_tree(size = 0.05) + geom_tiplab(size = 0.1,color = castree$color,align = T,
                                                                                     linesize = 0.05) + theme_tree() +
            ggtitle(paste0("Cassiopeia greedy tree with 50 clones of E15 rep ",i))
          p2.1
          pcas[[i]] = p2.1
        }
        ggexport(pcas,filename = "final_result/xtotal/benchmark/Cas_tree.pdf",width = 30,height =100)
        
        repcorsim.cas = SplitHclustTreeRepCor(castreels2, cellclls, nsample)
        # repcor.simcas2$sample = factor(repcor.simcas2$sample,levels = unique(repcor.simcas2$sample))
        
        ggplot(repcor.simcas2,aes(x = sample,y = cor, color = group)) + 
          geom_bar(stat = "identity",position = "dodge",fill = "white") + 
          theme_pubr() + 
          scale_color_lancet()+
          xlab("log10(clone number)")
        
      }
    }
    
    #haming&dclear anslysis
    {
      #save sim data
      simls = list(sim1,sim2,sim3)
      simscarls =list() 
      for (i in 1:3) {
        simscari = simls[[i]]
        simscari$tags = apply(simls[[i]][-1], 1, function(x){return(paste0(sprintf("%05d",x),collapse = ""))})
        print(length(unique(simscari$tags)))
        simscari = simscari[c(1,ncol(simscari))]
        simscarls[[i]] = list("scardata" = simscari,"states" = unique(unlist(as.matrix(simi[-1]))))
      }
      qsave(simscarls,"final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_sim2_data.qs")
      
      dctreesim1 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_sim_rep1_tree.rds")
      dctreesim2 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_sim_rep2_tree.rds")
      dctreesim3 = readRDS("final_result/CREST_final_analysis/benchmark/DCLEAR/DCLEAR_sim_rep3_tree.rds")
      dctreesimls = list(dctreesim1$hclust,dctreesim2$hclust,dctreesim3$hclust)
      dctreesimls2 = list(dctreesim1$hclust2,dctreesim2$hclust2,dctreesim3$hclust2)
    
      simcellls = list(simscarls[[1]]$scardata,simscarls[[2]]$scardata,simscarls[[3]]$scardata)
      repcorsim.hm = SplitHclustTreeRepCor(dctreesimls,simcellls,nsample)
      repcorsim.dc = SplitHclustTreeRepCor(dctreesimls2,simcellls,nsample)
      
      repcorsim.cas$method = "Cassiopeia"
      repcorsim.hm$method = "hamming"
      repcorsim.dc$method = "DCLEAR"
      
      repcorsim = rbind(repcorsim.hm,repcorsim.dc,repcorsim.cas)
      write.csv(repcorsim,file = "final_result/CREST_final_analysis/benchmark/E15_repcor_total_sim.csv",quote = F,
                row.names = F)
      
      
    }
    
    
  }
  
  #three method benchmark
  {
    
    colnames(repcor)[3] = colnames(repcorsim)[3] = colnames(repcorsim2)[3] = "clonenum"
    repcor$sample = "E15.5"
    repcorsim$sample = "simulate1"
    repcorsim2$sample = "simulate2"
    repcort = rbind(repcor,repcorsim,repcorsim2)
    library(ggsci)
    repcort[which(repcort$clonenum == "single"),"method"] = "single clone"
    repcort$clonenum = factor(repcort$clonenum,levels = unique(repcort$clonenum))
    psample = c(1,2,3,4,5,10,50,100,300,500,600,700,800,"single")
    
    
    repcort.plot = unique(repcort[which(repcort$clonenum %in% psample),][-c(1,2)])
    repcort.plot = repcort.plot[!duplicated(repcort.plot[c("clonenum","sample","method")]),]
    
    p4.0 = ggplot(repcort.plot,
                  aes(x = clonenum,y = cocor,fill = method)) + 
      # scale_color_lancet() +
      scale_fill_npg() +
      geom_bar(stat = "identity",position = position_dodge2(width = 0.9, preserve = "single")) + 
      theme_pubr() + 
      facet_wrap(~sample,ncol = 1) +
      xlab("Clone number") + 
      ylab("Correlation coefficient among replicates")
    p4.0
    
    benchmarkls = list("repcor" = repcort.plot,"figure" = p4.0)
    qsave(benchmarkls,"final_result/CREST_final_analysis/benchmark/benchmarkls.qs")
    datals$figure$benchmarkls = benchmarkls
    qsave(datals,"final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_17_edit.qs")
    ggsave(p4.0,filename = "final_result/CREST_final_analysis/benchmark/benchmark_tree_repcor_diffclonenumber_withsim.pdf",width = 7,height = 5)
    
  }
}

#all figure reanalysis 12.06-----------------------------
#load data
{
  # datalstmp = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_07.qs")
  datals = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_14.qs")
}

#fig5. Glu-DA fate analysis----------
{
  E15trans = qread("final_result/CREST_final_analysis/annotation_20221207/snapE15_FP_pseudotime.qs")
  clonecl = readRDS("final_result/CREST_final_analysis/annotation_20221207/e15.5 clone cluster/snapE15_clonecluster_20221209_cellcounts_cluster.rds")
  # clonecl1 = readRDS("final_result/CREST_final_analysis/annotation_20221207/e15.5 clone cluster/snapE15_V1ident_cellcounts_cluster.rds")
  # clonecl2 = readRDS("final_result/CREST_final_analysis/annotation_20221207/e15.5 clone cluster/snapE15_V2ident_cellcounts_cluster.rds")
  library(Seurat)
  hubc = c("Rgl1", "DA", "GLUFP", "NbDA", "NbGLU", "NbFP")
  clonecl1$cluster.fig
  clonecl2$cluster.fig
  #check data
  # DimPlot(E15trans)
  arraydf = datals$arraydata_bl$E15
  hubclone = names(clonecl$cluster[which(clonecl$cluster == "c6")])
  hubcb = names(E15trans$identsV1[which(E15trans$identsV1 %in% hubc)])
  Idents(E15trans) = E15trans$identsV1
  
  # arraydf = datals$arraydata_bl$E15v2
  # hubclone = names(clonecl2$cluster[which(clonecl2$cluster == "c6")])
  # hubcb = names(E15trans$identsV2[which(E15trans$identsV2 %in% hubc)])
  # Idents(E15trans) = E15trans$identsV2
  
  arraydfhub = arraydf[which(arraydf$pattern %in% hubclone),]
  hubst = dcast(arraydfhub,pattern~Cell.type)
  rownames(hubst) = hubst$pattern;hubst = hubst[-1]
  
  #clone plot
  {
    #2. 把两群细胞合并
    merhub = hubst
    merhub$DAt = merhub$DA + merhub$NbDA
    merhub$Glut = merhub$NbGLU + merhub$GLUFP
    merhub = merhub[c("NbFP","Rgl1","DAt","Glut")]
    merhub.up = as.matrix(merhub)
    merhub.up[which(merhub.up>0)] = T
    merhub.up[which(merhub.up==0)] = F
    
    #3. upset图
    library(UpSetR)
    merhubf = merhub[which(rowSums(merhub[3:4])>0),1:4]
    merhubf.up = as.matrix(merhubf)
    merhubf.up[which(merhubf.up>0)] = T
    merhubf.up[which(merhubf.up==0)] = F
    p3.1 = upset(as.data.frame(merhubf.up),nsets = 4, keep.order = T)
    p3.1
    # dir.create("final_result/CREST_final_analysis/DA_GLU_Clone_analysis")
    # ggexport(p3.1,filename = "final_result/CREST_final_analysis/DA_GLU_Clone_analysis/E15.5_DA_Glu_direction_stat_upset.pdf",width = 8,height = 6)
    #4. 对分类结果作图
    group.da = merhubf[which(merhubf$DAt > 0  & merhubf$Glut == 0),]
    group.glu = merhubf[which(merhubf$DAt == 0 & merhubf$Glut > 0),]
    group.all = merhubf[which(merhubf$DAt > 0 & merhubf$Glut > 0),]
    merhubfl = as.data.frame(t(t(as.matrix(merhubf))/colSums(merhubf)))
    merhubfl = cbind(rownames(merhubfl),merhubfl)
    merhubfl = melt(merhubfl)
    colnames(merhubfl) = c("tags","cell","proportion")
    merhubfl = merhubfl %>% group_by(tags) %>% summarise(cell = cell,proportion = proportion,
                                                         abspro = proportion/sum(proportion))
    merhubfl$group = "unclass"
    merhubfl[which(merhubfl$tags %in% rownames(group.da)),"group"] = "DA-bias"
    merhubfl[which(merhubfl$tags %in% rownames(group.glu)),"group"] = "Glu-bias"
    merhubfl[which(merhubfl$tags %in% rownames(group.all)),"group"] = "Dual-fate"
    merhubfl$prog = FALSE
    merhubfl[which(merhubfl$abspro > 0 & merhubfl$cell == "NbFP"),"prog"] = TRUE
    merhubfl = merhubfl[order(merhubfl$group,-merhubfl$prog),]
    merhubfl$tags = factor(merhubfl$tags,levels = unique(merhubfl$tags))
    merhubfl$cell = factor(merhubfl$cell,levels = c("DAt","Glut","NbFP","Rgl1"))
    p4.1 = ggplot(merhubfl,aes(x = tags,y = cell,fill = proportion)) + 
      geom_point(aes(size = abspro,color = group),shape = 21,alpha = 0.8) + 
      scale_size_continuous(range = c(0, 4)) +
      scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
      xlab("") + ylab("") +
      theme_bw() +
      theme(aspect.ratio = 1/6,
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "bottom")
    p4.1
    # ggsave(p4.1,filename = "final_result/CREST_final_analysis/DA_GLU_Clone_analysis/E15.5_DA_Glu_direction_stat_point_v1.pdf",width = 10,height = 3)
    
  }
  
  #DA portion
  {
    nbfpp = hubst[which(hubst$NbFP > 0 & (hubst$DA + hubst$NbDA) > 0 &
                              (hubst$GLUFP + hubst$NbGLU) == 0),]
    nbfpm = hubst[which(hubst$NbFP == 0 & (hubst$DA + hubst$NbDA) > 0 &
                              (hubst$GLUFP + hubst$NbGLU) == 0),]
    nbfppcb = arraydf[which(arraydf$pattern %in% rownames(nbfpp) & arraydf$Cell.type != "NbFP"),"Cell.BC"]
    nbfpmcb = arraydf[which(arraydf$pattern %in% rownames(nbfpm) & arraydf$Cell.type != "NbFP"),"Cell.BC"]
    nbfppcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfppcb)]
    nbfpmcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfpmcb)]
    nbfpst = data.frame("DA_pseudotime" = c(nbfppcb,nbfpmcb),
                        "type" = rep(c("NBFP+","NBFP-"),c(length(nbfppcb),length(nbfpmcb))))
    # nbfpp = hubst[which(hubst$NbFP > 0 & (hubst$DA + hubst$NbDA) > 0 ),]
    # nbfpm = hubst[which(hubst$NbFP == 0 & (hubst$DA + hubst$NbDA) > 0 ),]
    nbfpp$DA_percentage = (nbfpp$DA)/(nbfpp$DA + nbfpp$NbDA)
    nbfpm$DA_percentage = (nbfpm$DA)/(nbfpm$DA + nbfpm$NbDA)
    nbfpst2 = data.frame("DA_percentage" = c(nbfpp$DA_percentage,nbfpm$DA_percentage),
                        "type" = rep(c("NBFP+","NBFP-"),c(nrow(nbfpp),nrow(nbfpm))))
    
    
    my_comparisons <- list( c("NBFP+", "NBFP-"))
    # nbfpst = read.csv("final_result/xdraft/main_figure_2_25/E15.5_DA_percentage_withornot_NbFP.csv")
    
    p5.1 = ggboxplot(nbfpst, x = "type", y = "DA_pseudotime",
                   fill = "type", palette = "jco",
                   width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(method = "t.test",label.y = max(nbfpst$DA_pseudotime),comparisons = my_comparisons)    
    p5.1
    
    p5.2 = ggboxplot(nbfpst2, x = "type", y = "DA_percentage",
                     fill = "type", palette = "jco",
                     width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(method = "t.test",label.y = 1,comparisons = my_comparisons)    
    p5.2
    # ggsave(p4,filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_percentage_withornot_NbFP_3_10.pdf",width = 4,height = 4)
    # write.csv(nbfpst,"final_result/xdraft/main_figure_2_25/E15.5_DA_percentage_withornot_NbFP.csv")
    
    nbfpgp = hubst[which(hubst$NbFP > 0 & (hubst$DA + hubst$NbDA) == 0 &
                               (hubst$GLUFP + hubst$NbGLU) > 0),]
    nbfpgm = hubst[which(hubst$NbFP == 0 & (hubst$DA + hubst$NbDA) == 0 &
                               (hubst$GLUFP + hubst$NbGLU) > 0),]
    nbfpgpcb = arraydf[which(arraydf$pattern %in% rownames(nbfpgp)),"Cell.BC"]
    nbfpgmcb = arraydf[which(arraydf$pattern %in% rownames(nbfpgm)),"Cell.BC"]
    nbfpgpcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfpgpcb)]
    nbfpgmcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfpgmcb)]
    nbfpgst = data.frame("GLU_pseudotime" = c(nbfpgpcb,nbfpgmcb),
                        "type" = rep(c("NBFP+","NBFP-"),c(length(nbfpgpcb),length(nbfpgmcb))))
    
    nbfpgp$Glu_percentage = (nbfpgp$GLUFP)/(nbfpgp$GLUFP + nbfpgp$NbGLU)
    nbfpgm$Glu_percentage = (nbfpgm$GLUFP)/(nbfpgm$GLUFP + nbfpgm$NbGLU)
    nbfpgst2 = data.frame("Glu_percentage" = c(nbfpgp$Glu_percentage,nbfpgm$Glu_percentage),
                         "type" = rep(c("NBFP+","NBFP-"),c(nrow(nbfpgp),nrow(nbfpgm))))
    
    my_comparisons <- list( c("NBFP+", "NBFP-"))
    # nbfpgst = read.csv("final_result/xdraft/main_figure_2_25/E15.5_Glu_percentage_withornot_NbFP.csv")
    p5.3 = ggboxplot(nbfpgst, x = "type", y = "GLU_pseudotime",
                     fill = "type", palette = "jco",
                     width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(method = "t.test",label.y = max(nbfpgst$GLU_pseudotime),comparisons = my_comparisons)    
    p5.4 = ggboxplot(nbfpgst2, x = "type", y = "Glu_percentage",
                     fill = "type", palette = "jco",
                     width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(method = "t.test",label.y = 1,comparisons = my_comparisons)    
    p5.3+p5.4    
    # ggsave(p4.1,filename = "final_result/xdraft/main_figure_2_25/E15.5_Glu_percentage_withornot_NbFP_3_10.pdf",width = 4,height = 4)
    # write.csv(nbfpgst,"final_result/xdraft/main_figure_2_25/E15.5_Glu_percentage_withornot_NbFP.csv",row.names = F)
    
    # p5.2 = ggbetweenstats(
    #   data = nbfpst,
    #   x = type,
    #   y = DA_percentage,
    #   pairwise.display = "significant"
    # ) + xlab("") + ylab("DA_percentage")
    # 
    # p5
    # p4.1 = ggbetweenstats(
    #   data = nbfpgst,
    #   x = type,
    #   y = Glu_percentage,
    #   pairwise.display = "significant"
    # ) + xlab("") + ylab("Glu_percentage")
    nbfpap = hubst[which(hubst$NbFP > 0),]
    nbfpam = hubst[which(hubst$NbFP == 0),]
    nbfpgpcb = arraydf[which(arraydf$pattern %in% rownames(nbfpgp)),"Cell.BC"]
    nbfpgmcb = arraydf[which(arraydf$pattern %in% rownames(nbfpgm)),"Cell.BC"]
    nbfpgpcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfpgpcb)]
    nbfpgmcb = E15trans$pseudotime[which(names(E15trans$pseudotime) %in% nbfpgmcb)]
    nbfpgst = data.frame("GLU_pseudotime" = c(nbfpgpcb,nbfpgmcb),
                         "type" = rep(c("NBFP+","NBFP-"),c(length(nbfpgpcb),length(nbfpgmcb))))
  }
  
  #DEG find
  {
    trans.da = arraydf[which(arraydf$pattern %in% rownames(group.da) & arraydf$Cell.type == "NbFP"),"Cell.BC"]
    trans.glu = arraydf[which(arraydf$pattern %in% rownames(group.glu) & arraydf$Cell.type == "NbFP"),"Cell.BC"]
    trans.all = arraydf[which(arraydf$pattern %in% rownames(group.all) & arraydf$Cell.type == "NbFP"),"Cell.BC"]
    
    cell_names.da = trans.da
    cell_names.glu = trans.glu
    cell_names.all = trans.all
    
    p6.4.0 = DimPlot(E15trans,label = TRUE,label.size = 3) + NoLegend()
    p6.4.1 = DimPlot(E15trans,label = TRUE, label.size = 3,cells.highlight = cell_names.da, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.2 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.glu, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.3 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.all, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4 = p6.4.0 + p6.4.1 + p6.4.2 + p6.4.3
    p6.4
    
    
    #DEG calculate
    E15trans.hub = subset(E15trans, idents = c("NbFP"))
    tmp = Idents(E15trans.hub)
    tmp = as.character(tmp);names(tmp) = names(Idents(E15trans.hub))
    tmp[which(names(tmp) %in% trans.da)] = "DA_d"
    tmp[which(names(tmp) %in% trans.glu)] = "Glu_d"
    tmp[which(names(tmp) %in% trans.all)] = "All_d"
    tmp = as.factor(tmp)
    Idents(E15trans.hub) = tmp
    E15trans.hub = subset(E15trans.hub, idents = c("DA_d","Glu_d","All_d"))
    DA.markers <- FindMarkers(E15trans.hub, ident.1 = "DA_d",min.pct = 0.25)
    Glu.markers <- FindMarkers(E15trans.hub, ident.1 = "Glu_d", min.pct = 0.25)
    All.markers <- FindMarkers(E15trans.hub, ident.1 = "All_d", min.pct = 0.25)
    DA_Glu.markers <- FindMarkers(E15trans.hub, ident.1 = "DA_d",ident.2 = "Glu_d", min.pct = 0.25)
    levels(E15trans.hub) = c("DA_d","All_d","Glu_d")
    pvc = VocanoFigure(DA_Glu.markers,title = "DA biased vs Glu biased NbFP cell",thres = 1)
    pvc
    #PCA Plot
    E15trans.hub <- FindVariableFeatures(E15trans.hub, selection.method = "vst", nfeatures = 2000)
    E15trans.hub <- RunPCA(E15trans.hub, features = VariableFeatures(object = E15trans.hub))
    p6.6.1 = DimPlot(E15trans.hub, reduction = "pca")
    p6.6.2 = DimHeatmap(E15trans.hub, dims = 1:2,fast = F, balanced = TRUE)
    pcmx = as.data.frame(E15trans.hub@reductions$pca@cell.embeddings)
    pcmx$class = Idents(E15trans.hub)
    pcmx = pcmx[,c("PC_1","PC_2","class")]
    pcmxl = melt(pcmx)
    p6.6.3.0 = ggplot(pcmx,aes(x = PC_1, y = PC_2, col = class)) + geom_point() + 
      theme_pubr()+
      theme(legend.position = c(0.9, 0.9))
    dens1 <- ggplot(pcmx, aes(x = PC_1, fill = class)) + 
      geom_density(alpha = 0.5) + 
      theme_void() + 
      theme(legend.position = "none")
    dens2 <- ggplot(pcmx, aes(x = PC_2, fill = class)) + 
      geom_density(alpha = 0.5) + 
      theme_void() + 
      theme(legend.position = "none") + 
      coord_flip()
    library(patchwork)
    p6.6.3.0 = dens1 + plot_spacer() + p6.6.3.0 + dens2 + 
      plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    p6.6.3.0
    
    library(ggpubr)
    p6.6.3.1 = ggplot(pcmx,aes(y = PC_1, x = class)) + 
      geom_boxplot(width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = list(c("DA_d","All_d"), c("DA_d","Glu_d"), c("All_d","Glu_d")),
                         label = "p.signif",
                         method = "t.test") +
      theme(legend.position = c(5, 20))
    
    p6.6.3.2 = ggplot(pcmx,aes(y = PC_2, x = class)) + 
      geom_boxplot(width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = list(c("DA_d","All_d"), c("DA_d","Glu_d"), c("All_d","Glu_d")),
                         label = "p.signif", method = "t.test")
    p6.6.3.3 = p6.6.3.1+p6.6.3.2
    
    dagluid = as.data.frame(Idents(E15trans.hub))
    
  }
  figurels = list("f.DA_Glu_upset" = p3.1,"e.DA_Glu_clone_point_stat" = p4.1,
                  "i.DA_Glu_NbFP_PCA_point" = p6.6.3.0,"j.DA_Glu_NbFP_PCA1" = p6.6.3.1,
                  "k.DA_Glu_NbFP_PCA2" = p6.6.3.2,"l.DA_Glu_NbFP_DEG_VocanoPlot" = pvc,
                  "g.DA_Glu_NbFP_DA_psudotime" = p5.1,"g.DA_Glu_NbFP_DA_percentage" = p5.2,
                  "g.DA_Glu_NbFP_Glu_psudotime" = p5.3,"g.DA_Glu_NbFP_Glu_percentage" = p5.4,
                  "DA_Glu_NbFP_SeuratDEGPlot" = list("Variable_gene" = p6.6.2,
                                                     "clone_dimplot" = p6.4))
  degdatals = list("Cell_group" = dagluid,
                   "DEG_marker" = list("DA.bias" = DA.markers,"Glu.bias" = Glu.markers,
                                       "Dual.bias" = All.markers,"DA_vs_Glu" = DA_Glu.markers))
  #NW model
  #模型构建
  {
    #输入数据集
    E15trans.hubdg = subset(E15trans.hub, idents = c("DA_d","Glu_d"))
    E15trans.hubdg = FindVariableFeatures(E15trans.hubdg, selection.method = "vst", nfeatures = 32285)
    DGdegs = FindMarkers(E15trans.hubdg,ident.1 = "DA_d",min.pct = 0.25)
    # write.csv(DGdegs,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_DAvsGlu.csv",row.names = T,quote = F)
    
    library(pROC)
    
    trans.da = subset(E15trans.hubdg, idents = "DA_d")
    trans.glu = subset(E15trans.hubdg, idents = "Glu_d")
    damx = as.data.frame(trans.da[["RNA"]]@counts)
    glumx = as.data.frame(trans.glu[["RNA"]]@counts)
    
    NeurBuild = function(markers){
      damxm = damx[rownames(damx) %in% markers,]
      glumxm = glumx[rownames(glumx) %in% markers,]
      damxm = as.data.frame(t(damxm))
      glumxm = as.data.frame(t(glumxm))
      
      Normmx = function(x){
        return((x-min(x))/(max(x)-min(x)))
      }
      damxmn =  as.data.frame(t(apply(damxm, 1, Normmx)))
      glumxmn =  as.data.frame(t(apply(glumxm, 1, Normmx)))
      damxmn$type =  0
      glumxmn$type = 1
      
      dadata = rbind(damxmn,glumxmn)
      colnames(dadata)[1:(ncol(dadata)-1)] = paste0("x",1:(ncol(dadata)-1))
      
      testid = sample(1:nrow(dadata),as.integer(nrow(dadata)/2))
      testdata = dadata[testid,]
      traindata = dadata[setdiff(1:nrow(dadata),testid),]
      
      f <- as.formula(paste0("type ~ ",
                             paste(names(traindata)[1:(ncol(traindata) - 1)],
                                   collapse = " + ")))
      library(neuralnet)
      neur <- neuralnet(f, data = traindata, hidden = 10,
                        act.fct = "logistic", linear.output = T,
                        err.fct = "sse", rep = 1)
      test.hat = compute(neur,testdata[-ncol(testdata)])$net.result
      
      test.predict = test.hat[,1]
      modelroc = roc(testdata$type,test.predict)
      
      return(modelroc)
    }
    marker1 = rownames(DGdegs)
    marker2 = head(VariableFeatures(E15trans.hubdg), length(marker1))
    marker3 = rownames(E15trans.hubdg)[sample(1:length(rownames(E15trans.hubdg)),length(marker1))]
    marker4 = read.delim("raw_data/tf.name.txt",header = F)
    marker4 = marker4$V2
    marker4 = head(VariableFeatures(E15trans.hubdg)[VariableFeatures(E15trans.hubdg) %in% marker4],length(marker1))
    
    aucdf = NULL
    for (i in 1:100) {
      roc1 = NeurBuild(marker1)
      roc2 = NeurBuild(marker2)
      roc3 = NeurBuild(marker3)
      roc4 = NeurBuild(marker4)
      aucdf = rbind(aucdf,data.frame(sample = i, markers = c(paste0("Differentially expressed(n=",length(marker1),")"),
                                                             paste0("Highly variable genes(n=",length(marker2),")"),
                                                             paste0("Random genes(n=",length(marker3),")"),
                                                             paste0("Transcription factors(n=",length(marker4),")")),
                                     AUC = c(roc1$auc,roc2$auc,roc3$auc,roc4$auc)))
    }
    
    aucdf$markers = factor(aucdf$markers,levels = c(paste0("Differentially expressed(n=",length(marker1),")"),
                                                    paste0("Highly variable genes(n=",length(marker2),")"),
                                                    paste0("Random genes(n=",length(marker3),")"),
                                                    paste0("Transcription factors(n=",length(marker4),")")))
    p7 = ggboxplot(aucdf, x = "markers", y = "AUC",
                   fill = "markers", palette = "jco",
                   width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      stat_compare_means(method = "anova", label.y = 1) +
      stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = paste0("Differentially expressed(n=",length(marker1),")"))    
    p7
    aucst = aucdf %>% group_by(markers) %>% summarise(mean(AUC))

  }
  figurels$o.DA_Glu_neuralnet_auc_boxplot = p7
  degdatals$DA_Glu_neuralnet_aucdf = aucdf
  degdatals$DA_Glu_neuralnet_aucdf_stat = aucst
  
  DAvsGluls = list("figurels" = figurels, "datals" = degdatals)
  qsave(DAvsGluls,file = "final_result/CREST_final_analysis/DA_GLU_Clone_analysis/DAvsGlulsV1.qs")
  DAvsGluls$figurels$e.DA_Glu_clone_point_stat
  DAvsGluls$figurels$o.DA_Glu_neuralnet_auc_boxplot
  
  datals$figure$DAvsGluls_f5 = DAvsGluls
  # DAvsGlulsv1 = qread(file = "final_result/CREST_final_analysis/DA_GLU_Clone_analysis/DAvsGlulsV1.qs")
  # 
  # DAvsGlulsv1$figurels$o.DA_Glu_neuralnet_auc_boxplot + DAvsGluls$figurels$o.DA_Glu_neuralnet_auc_boxplot
  # DAvsGlulsv1$figurels$l.DA_Glu_NbFP_DEG_VocanoPlot + DAvsGluls$figurels$l.DA_Glu_NbFP_DEG_VocanoPlot
  # DAvsGlulsv1$figurels$i.DA_Glu_NbFP_PCA_point | DAvsGluls$figurels$i.DA_Glu_NbFP_PCA_point
  # DAvsGlulsv1$figurels$j.DA_Glu_NbFP_PCA1 | DAvsGluls$figurels$j.DA_Glu_NbFP_PCA1
  datals$figure$DAvsGluls_f5$figurels$f.DA_Glu_upset
  datals$figure$DAvsGluls_f5$figurels$l.DA_Glu_NbFP_DEG_VocanoPlot
  datals$figure$DAvsGluls_f5$figurels$g.DA_Glu_NbFP_DA_psudotime
  
}

#fig3. E11 E15桑基图------------
{
  library(qs)
  datals = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_28_edit.qs")
  # arraydf = datals$arraydata_bl$E15

  earlyid11 = c('NPBL','NPBM','Rgl1')
  earlyid15 = c('Rgl2AL','NbGABAAL','Rgl3','NbGABABL','Rgl2BP','NbGABABI','NbBM1','NbFP',"Rgl1")
  
  
  cellorder15 = c('NbGABAAL','Rgl2AL','NbGABABL','NbGABABI','NbBM1','NbFP',"Rgl1",'Rgl2BP','Rgl3',
                  'GLUAL1','GLUAL2','GLUAL3','GLUAL4','GABAAL','GABABL1',
                  'GABABL2','GABABL3','GABABL4','GABABI','GLUBM1','GLUBM2','DA','NbDA','GLUFP','NbGLU',
                  'GLUMHB1','GLUMHB2','GABAMHB','GLUMHB3','GLUMHB4')
  domain15 = read.csv("final_result/CREST_final_analysis/annotation_20221207/E15_spatial_value.csv")
  domain11 = read.csv("final_result/CREST_final_analysis/annotation_20221207/E11_spatial_value.csv")
  colnames(domain15)[c(2,3,4)] = c("Cell.type","domain","domainid")
  colnames(domain11)[c(2,3,4)] = c("Cell.type","domain","domainid")
  # paste0(domain15[order(domain15$domainid),"Cell.type"],collapse = "','")
  BCEnrichment2gen = function(arraydf,domain,earlyid){
    cmtn = merge(arraydf,domain,by = "Cell.type")
    lateid = setdiff(unique(cmtn$Cell.type),earlyid)
    cmtn$phase = "late"
    cmtn[which(cmtn$Cell.type %in% earlyid),"phase"] = "early"
    cmtn = cmtn[order(cmtn$phase,cmtn$domain),]
    {
      cmtnp = cmtn %>% group_by(pattern,phase,Cell.type) %>% 
        summarise(counts = sum(n()),domain = unique(domain))
      #
      cellmx.e = dcast(cmtnp[which(cmtnp$phase == "early"),c("pattern","Cell.type","counts")],
                       pattern~Cell.type,value.var = "counts",fun.aggregate = sum)
      cellmx.l = dcast(cmtnp[which(cmtnp$phase == "late"),c("pattern","Cell.type","counts")],
                       pattern~Cell.type,value.var = "counts",fun.aggregate = sum)
      cellft = merge(cellmx.e,cellmx.l,by = "pattern")
      
    }
    cellft = cellft[which(rowSums(cellft[-1]>0) < ncol(cellft[-1])/3),]
    library(ggalluvial)
    phasest = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                         "phase2" = rep(lateid,length(earlyid)), "Freq.from" = 0, "Freq.to" = 0)
    for (i in 1:nrow(cellft)) {
      line = cellft[i,]
      line.e = unlist(line[,which(colnames(line) %in% earlyid)])
      line.e = line.e[which(line.e>0)]
      line.l = unlist(line[,which(colnames(line) %in% lateid)])
      line.l = line.l[which(line.l>0)]
      for (j in 1:length(line.e)) {
        for (k in 1:length(line.l)) {
          seg.from = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.from"]
          seg.to = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.to"]
          phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.from"] = seg.from + line.e[j]
          phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.to"] = seg.to + line.l[k]
        }
      }
    }
    #method2
    rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
      vec <- rnorm(N, M/N, sd)
      if (abs(sum(vec)) < 0.01) vec <- vec + 1
      vec <- round(vec / sum(vec) * M)
      deviation <- M - sum(vec)
      for (. in seq_len(abs(deviation))) {
        vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
      }
      if (pos.only) while (any(vec < 0)) {
        negs <- vec < 0
        pos  <- vec > 0
        vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
        vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
      }
      vec
    }
    
    phasestml = list()
    for (t in 1:100) {
      cellftm = cellft
      for (i in 2:ncol(cellftm)) {
        cellftm[i] = rand_vect(nrow(cellft),sum(cellft[i]))
      }
      phasestm = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                            "phase2" = rep(lateid,length(earlyid)), "Freq.from" = 0, "Freq.to" = 0)
      for (i in 1:nrow(cellftm)) {
        line = cellftm[i,]
        line.e = unlist(line[,which(colnames(line) %in% earlyid)])
        line.e = line.e[which(line.e>0)]
        line.l = unlist(line[,which(colnames(line) %in% lateid)])
        line.l = line.l[which(line.l>0)]
        for (j in 1:length(line.e)) {
          for (k in 1:length(line.l)) {
            seg.from = phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.from"]
            seg.to = phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.to"]
            phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.from"] = seg.from + line.e[j]
            phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.to"] = seg.to + line.l[k]
          }
        }
      }
      
      phasestml[[t]] = phasestm
    }
    phasestmle = phasestmll = list()
    for (j in 1:nrow(phasestml[[1]])) {
      phasestmle[[j]] = phasestml[[1]][j,3]
      phasestmll[[j]] = phasestml[[1]][j,4]
    }
    
    for (i in 2:length(phasestml)) {
      line = phasestml[[i]]
      for (j in 1:nrow(line)) {
        phasestmle[[j]] = c(phasestmle[[j]],line[j,3])
        phasestmll[[j]] = c(phasestmll[[j]],line[j,4])
      }
    }
    return(list("cmtn" = cmtn,"cellft" = cellft,"phasest" = phasest,"phasestmle" = phasestmle,
                "phasestmll" = phasestmll))
  }
  BCEnrichment2 = function(allu,cellorder = NULL,domaincol){
    #绘制桑基图
    {
      phasestf = NULL
      phasest = allu$phasest
      phasestmle = allu$phasestmle
      phasestmll = allu$phasestmll
      cmtn = allu$cmtn
      for (i in 1:nrow(phasest)) {
        line = phasest[i,]
        if(line$Freq.from > quantile(phasestmle[[i]],probs = 0.1) | 
           line$Freq.to > quantile(phasestmll[[i]],probs = 0.1)){
          phasestf = rbind(phasestf,line)
        }
        
      }
      
      domain_df = unique(cmtn[c("Cell.type","domainid","domain")])
      domain_df = domain_df[order(domain_df$domain),]
      # phasestf = phasest[which(phasest$Freq.from > 0.001 | phasest$Freq.to > 0.001),]
      # phasestf = phasestf[which(phasestf$Freq>20),]
      phasest_lod = to_lodes_form(phasestf,
                                  key = "phase", value = "Cell.type", id = "alluvium",
                                  axes = 1:2)
      phasest_lod$Freq = 0
      freqfrom = phasest_lod[which(phasest_lod$phase == "phase1"),"Freq.from"]
      freqto = phasest_lod[which(phasest_lod$phase == "phase2"),"Freq.to"]
      phasest_lod[which(phasest_lod$phase == "phase1"),"Freq"] = freqfrom/sum(freqfrom)
      phasest_lod[which(phasest_lod$phase == "phase2"),"Freq"] = freqto/sum(freqto)
      
      phasest_lod = merge(phasest_lod,domain_df,
                          by = "Cell.type")
      phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domainid),]
      if(length(cellorder)== 0){
        phasest_lod$Cell.type = factor(phasest_lod$Cell.type,
                                       levels = unique(as.character(phasest_lod$Cell.type)))
      }else{
        phasest_lod$Cell.type = factor(phasest_lod$Cell.type,
                                       levels = cellorder[which(cellorder %in% unique(as.character(phasest_lod$Cell.type)))])
      }
      
      # pv.phasestt  = allu$pv.phasestt
      # pv.phasestt$sig = "N"
      # pv.phasestt[which(),]
      
      pf1 = ggplot(phasest_lod,
                   aes(x = phase, stratum = Cell.type, alluvium = alluvium,
                       y = Freq,
                       fill = domain,
                       # color = Cell.type,
                       label = Cell.type)) +
        scale_x_discrete(expand = c(.1, .1)) +
        scale_fill_manual(values = domaincol) +
        # geom_flow(aes(color = Cell.type)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      pf1
      
      
    }
    allu$figure_alluv = pf1
    return(allu)
  }
  CalAlluPvalue = function(cellallu11){
    cellallu11$phasest$pve = 1
    cellallu11$phasest$pvl = 1
    
    for (i in 1:length(cellallu11$phasestmle)) {
      ldiste = cellallu11$phasestmle[[i]]
      ldistl = cellallu11$phasestmll[[i]]
      ldiste = ldiste[order(ldiste)]
      ldistl = ldiste[order(ldistl)]
      pve = length(ldiste[which(ldiste > cellallu11$phasest$Freq.from[i])])/length(ldiste)
      pvl = length(ldiste[which(ldistl > cellallu11$phasest$Freq.to[i])])/length(ldistl)
      cellallu11$phasest$pve[i] = pve
      cellallu11$phasest$pvl[i] = pvl
    }
    
    # install.packages("ggpol")
    library(gghalves)
    library(ggpol)
    
    nphasest = cellallu11$phasest[which((cellallu11$phasest$Freq.from + cellallu11$phasest$Freq.to)>0),]
    phasestt = rbind(data.frame("from" = nphasest$phase1,
                                "to" = nphasest$phase2,
                                "size" = nphasest$Freq.from,
                                "pvalue" = nphasest$pve,
                                "phase" = "from"),
                     data.frame("from" = nphasest$phase1,
                                "to" = nphasest$phase2,
                                "size" = nphasest$Freq.to,
                                "pvalue" = nphasest$pvl,
                                "phase" = "to"))
    
    library(ggforce)
    df <- data.frame(start = rep(c(-pi/2, pi/2), 3),
                     type = rep(c("Investors", "Assignees"), 3),
                     country = rep(c("Japan", "Germany", "Korea"), each = 2),
                     x = rep(c(1, 2, 3), each = 2),
                     y = rep(c(3, 1, 2), each = 2),
                     count = c(19419, 1132, 8138, 947, 8349, 436))
    
    
    phasestt$start = -1.570796
    phasestt[which(phasestt$phase == "to"),"start"] = 1.570796
    phasestt$x = as.numeric(as.factor(phasestt$from))
    phasestt$y = as.numeric(as.factor(phasestt$to))
    
    r <- 0.5
    scale <- r/max(sqrt(phasestt$size))
    phasestt.sig = phasestt[which(phasestt$pvalue < 0.05),]
    p1.2 = ggplot(phasestt) + 
      geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = scale*sqrt(size + .05), 
                       alpha = 1-pvalue,
                       start = start, end = start + pi, fill = phase),
                   color = "white") +
      geom_text(data = phasestt.sig,
                aes(label = "*", x = x, y = y - start/1.570796*scale*sqrt(size)/2),
                size =10/.pt, vjust = 0.7)+ 
      scale_fill_manual(values = c("#D9534F","#FFAD60")) +
      # scale_fill_viridis(discrete = T) +
      # guides(fill = guide_legend(title = "Type", reverse = T)) +
      xlab("") + ylab("") +
      scale_x_continuous(breaks = 1:max(phasestt$x),
                         labels = levels(as.factor(phasestt$from))) +
      scale_y_continuous(breaks = 1:max(phasestt$y),
                         labels = levels(as.factor(phasestt$to))) +
      coord_fixed() +
      theme_bw() + theme(axis.text.x = element_text(angle = 45,vjust = 0.5)) 
    
    p1.2
    cellallu11$pv.phasestt = phasestt
    cellallu11$figure_allupvalue = p1.2
    return(cellallu11)
  }
  
  names(datals$arraydata_bl)
  sample_todo = c("E11_rep1","E11_rep2","E11_rep3","E11","E15_rep1","E15_rep2","E15_rep3","E15")
  earlyidls = list(earlyid11,earlyid11,earlyid11,earlyid11,earlyid15,earlyid15,earlyid15,earlyid15)
  domainls = list(domain11,domain11,domain11,domain11,domain15,domain15,domain15,domain15)
  
  # cellalluls = list()
  for (i in 1:4) {
    print(i)
    arraydf = datals$arraydata_bl[[sample_todo[i]]]
    # cellallui = cellalluls[[i]]
    cellallui = BCEnrichment2gen(arraydf,domainls[[i]],earlyidls[[i]])
    cellallui = BCEnrichment2(cellallui,domaincol)
    # cellallui = BCEnrichment2(cellallui,cellorder15)
    cellalluls[[i]] = CalAlluPvalue(cellallui)
    names(cellalluls)[i] = sample_todo[i]
    # cellalluls[[i]]$figure_alluv = cellallui$figure_alluv
  }
  
  # dir.create("final_result/CREST_final_analysis/E11_E15_alluv_lineage")
  qsave(cellalluls,"final_result/CREST_final_analysis/E11_E15_alluv_lineage/alluls_2.qs")
  datals$figure$alluls_f3 = cellalluls
  
  
  #predict E11 to E15
  {
    earlyid11 = c('NPBL','NPBM','Rgl1')
    earlyid15 = c('NbGABABL','NbGABABI','NbBM1','NbFP')
    c('Rgl2AL','NbGABAAL','Rgl3','NbGABABL','Rgl2BP','NbGABABI','NbBM1','NbFP',"Rgl1")
    
    oldid15 = c('GLUAL1','GLUAL2','GLUAL3','GLUAL4','GABAAL','GABABL1',
                'GABABL2','GABABL3','GABABL4','GABABI','GLUBM1','GLUBM2',
                'DA','NbDA','GLUFP','NbGLU',
                'GLUMHB1','GLUMHB2','GABAMHB','GLUMHB3','GLUMHB4')
    arraydf = datals$arraydata_bl$E11
    arraydf = arraydf[which(arraydf$Cell.type %in% c(earlyid11,earlyid15)),]
    cellallui = BCEnrichment2gen(arraydf,domain11,earlyid11)
    cellallui = BCEnrichment2(cellallui,NULL,domaincol)
    cellallui = CalAlluPvalue(cellallui)
    
    arraydf15 = datals$arraydata_bl$E15
    arraydf15 = arraydf15[which(arraydf15$Cell.type %in% c(earlyid15,unique(datals$arraydata_bl$OGN$Cell.type))),]
    cellallui15 = BCEnrichment2gen(arraydf15,domain15,earlyid15)
    cellallui15 = BCEnrichment2(cellallui15,NULL,domaincol)
    cellallui15 = CalAlluPvalue(cellallui15)
    pfpre = cellallui$figure_alluv + NoLegend() + cellallui15$figure_alluv
    pfpre
    ggsave(pfpre,filename = "final_result/CREST_final_analysis/E11_E15_alluv_lineage/E11_E15_alluv_predicti.pdf",
           width = 10,height = 8)
    
  }
  
  
  
}

#fig4. E11 E15 domain & clone stat-----------
{
  domaincol = c("AL" = "#D0D520","BL" = "#EEAF31","BI" = "#B71724","BM" = "#571F7C","FP" = "#192C78","unknown" = "grey")
  arraydf = datals$arraydata_bl$E15
  domain15 = read.csv("final_result/CREST_final_analysis/annotation_20221207/E15_spatial_value.csv")
  domain11 = read.csv("final_result/CREST_final_analysis/annotation_20221207/E11_spatial_value.csv")
  colnames(domain15)[c(2,3,4)] = c("Cell.type","domain","domainid")
  colnames(domain11)[c(2,3,4)] = c("Cell.type","domain","domainid")
  domain15[is.na(domain15$domain),"domain"] = "unknown"
  domain11[is.na(domain11$domain),"domain"] = "unknown"
  
  PseudoDisper = function(arraydf,domain){
    library(scatterpie)
    library(viridis)
    library(hrbrthemes)
    # cmtn = merge(cm,pseudo,by.x = "Var1",by.y = "BC")
    cmtn = merge(arraydf,domain,by = "Cell.type")
    cmtn = cmtn[order(cmtn$domainid),]
    cmtn$Cell.type = factor(as.character(cmtn$Cell.type),levels = unique(cmtn$Cell.type))
    unique(domain$domain)
    
    #
    # cmtnf = cmtn[which(cmtn$domain != "unknown"),]
    cmtnf = cmtn[which(cmtn$domain != "unknown"),]
    cmtnf$Cell.type = factor(cmtnf$Cell.type,levels = unique(cmtnf$Cell.type))
    cmtnf = cmtnf[order(cmtnf$domainid),]
    cmtnf$domainid = factor(cmtnf$domainid)
    # cmtnf$domainid = factor(as.numeric(factor(cmtnf$domain,levels = unique(cmtnf$domain))))
    
    root_spr  =NULL
    root_spr_mod  =NULL
    for (i in 1:length(unique(cmtnf$pattern))) {
      root = unique(cmtnf$pattern)[i]
      slice = cmtnf[which(cmtnf$pattern == root),]
      # pseudom = mean(slice$pseudo)
      
      {
        cc_m = table(slice$domainid)
        sdm_m = cc_m[which(cc_m>0)]
        sdm_m = sdm_m[order(-sdm_m)]
        sdm_m = sdm_m/sum(sdm_m)
        H_m = -sum(sdm_m*log(sdm_m))
        if(H_m != 0){
          Hwi_m = H_m/log2(length(sdm_m))
        }else{
          Hwi_m = 0
        }
        bias = 1
        
        sdid_m = names(sdm_m)
        if("unknown" %in% slice$domain){
          sdid_m = sdid_m[which(sdid_m!="unkown")]
          if(length(sdid_m)>1){
            bias = 1 + 0.1*(sd(as.numeric(sdid_m))+1)*(length(sdid_m)+1) 
          }else if(length(sdid_m)==1){
            bias = 1 + 0.1*1*2 
          }
        }else if(length(sdid_m)>1){
          bias = 1 + 0.1*sd(as.numeric(sdid_m))*length(sdid_m)
        }
        
        
        Hb = H_m*bias
        
        mainModel_m = names(table(slice$domain))[which.max(table(slice$domain))]
        meanModel = mean(as.numeric(as.character(slice$domainid)))
        line_m = data.frame(group = root,H = H_m,Hwi = Hwi_m,Hb = Hb,
                            main.domain = mainModel_m,
                            mean.model = meanModel,
                            # pseudo.mean = pseudom,
                            count = nrow(slice),
                            stringsAsFactors = F)
        cc_m = as.data.frame.list(cc_m)
        line_m = cbind(line_m,cc_m)
        root_spr_mod = rbind(root_spr_mod,line_m)
        
        mainModel = names(table(slice$Cell.type))[which.max(table(slice$Cell.type))]
        
        
        
        #cell H
        cc = table(slice$Cell.type)
        sdm = cc[which(cc>0)]
        sdm = sdm[order(-sdm)]
        sdm = sdm/sum(sdm)
        H_c = -sum(sdm*log(sdm))
        if(H_c != 0){
          Hwi_c = H_c/log2(length(sdm))
        }else{
          Hwi_c = 0
        }
        
        line = data.frame(group = root,H = H_m,Hwi = Hwi_m,Hb = Hb,
                          Hc = H_c,Hwic = Hwi_c,
                          main.cell = mainModel,
                          mean.model = meanModel,
                          # pseudo.mean = pseudom,
                          count = nrow(slice),
                          stringsAsFactors = F)
        cc = table(slice$Cell.type)
        cc = as.data.frame.list(cc)
        line = cbind(line,cc)
        root_spr = rbind(root_spr,line)
      }
      
    }
    colnames(root_spr_mod)[-(1:7)] = unique(cmtnf$domain)
    
    root_sprn = root_spr
    root_sprn[-(1:9)] = t(t(as.matrix(root_sprn[-c(1:9)]))/colSums(root_sprn[-(1:9)]))
    root_sprn = root_sprn[which(root_sprn$count >1),]
    root_sprn_mod = root_spr_mod
    root_sprn_mod[-(1:7)] = t(t(as.matrix(root_sprn_mod[-c(1:7)]))/colSums(root_sprn_mod[-(1:7)]))
    root_sprn_mod = root_sprn_mod[which(root_sprn_mod$count >1),]
    
    # tem2 <- as.data.frame(spline(root_sprn$mean.model, root_sprn$Hb,n = nrow(root_sprn),ties = mean))
    p5.1 = ggplot() + 
      geom_scatterpie(data = root_sprn,aes(x=mean.model, y=Hb, group=group,
                                           r = log10(count)/15,alpha = 0.5),color=NA,
                      cols=colnames(root_sprn)[-(1:9)]) + 
      # geom_histogram(data = root_sprn,aes(mean.model),bins = 100) +
      coord_fixed() +
      # geom_smooth(data = tem2, aes(x=x, y=y), se = F, method = 'loess',color = "black") +
      scale_fill_viridis(discrete=TRUE) +
      # scale_fill_manual(values = domaincol) +
      scale_y_continuous(breaks = c(0,1,2,3)) + 
      scale_x_continuous(breaks = c(1,2,3,4,5)) + 
      ylab("Dispersion score(DS)") + xlab("")+
      theme_classic() + theme(legend.position="none")
    p5.1
    p5.1.1 = ggplot(root_sprn,aes(x=mean.model)) + 
      geom_histogram(bins = 200,fill = "black") + 
      theme_bw() +
      theme(axis.title.x = element_blank(),panel.grid = element_blank())
    
    p5.1.1
    p5.1.2 = ggplot(root_sprn,aes(x=Hb)) + 
      geom_histogram(bins = 200,fill = "black") + 
      xlab("Dispersion score") + ylab("Clone Count") + 
      theme_classic()+
      coord_flip()
    p5.1.2
    
    p5.1.t = ggarrange(p5.1,p5.1.2,p5.1.1,heights = c(1,0.2),widths = c(1,0.4),
                       ncol = 2, nrow = 2)
    
    p5.2 = ggplot() + 
      geom_scatterpie(data = root_sprn_mod,aes(x=mean.model, y=Hb, group=group,
                                           r = log10(count)/15,alpha = 0.5),color=NA,
                      cols=colnames(root_sprn_mod)[-(1:7)]) + 
      # geom_histogram(data = root_sprn,aes(mean.model),bins = 100) +
      coord_fixed() +
      # geom_smooth(data = tem2, aes(x=x, y=y), se = F, method = 'loess',color = "black") +
      # scale_fill_viridis(discrete=TRUE) +
      scale_fill_manual(values = domaincol) +
      scale_y_continuous(breaks = c(0,1,2,3)) + 
      scale_x_continuous(breaks = c(1,2,3,4,5)) + 
      xlab("") + ylab("Dispersion score(DS)") + 
      theme_classic() + theme(legend.position="none")
    p5.2
    p5.2.1 = ggplot(root_sprn_mod,aes(x=mean.model)) + 
      # scale_fill_continuous() +
      geom_histogram(bins = 200,fill = "black") + 
      theme_bw() + ylab("Clone count") +
      theme(axis.title.x = element_blank(),panel.grid = element_blank())
    # axis.text.x = element_blank(),)
    p5.2.1
    
    p5.2.2 = ggplot(root_sprn_mod,aes(x=Hb)) + 
      # scale_fill_continuous() +
      geom_histogram(bins = 200,fill = "black") + 
      xlab("Dispersion score") + ylab("Clone Count") + 
      theme_classic()+theme(panel.grid = element_blank())+
      coord_flip()
    # axis.text.x = element_blank(),)
    p5.2.2
    
    
    p5.2.t = ggarrange(p5.2,p5.2.2,p5.2.1,heights = c(1,0.2),widths = c(1,0.4),
                       ncol = 2, nrow = 2)
    
    root_spr_modn = root_spr_mod[which(root_spr_mod$count>1),]
    root_spr_modn[-(1:7)] = t(t(as.matrix(root_spr_modn[-c(1:7)]))/colSums(root_spr_modn[-(1:7)]))
    
    #5.2
    root_spr_modn$max = 0
    for (i in 1:nrow(root_spr_modn)) {
      line = root_spr_modn[i,]
      root_spr_modn[i,"main.domain"] = names(which.max(line[-(1:7)]))
      root_spr_modn$max[i] = max(line[-(1:6)]/sum(line[-(1:7)]))
    }
    p5.3 = ggplot(root_spr_modn,aes(x = main.domain, y = Hb, fill = main.domain)) +
      geom_boxplot(width = 0.6) + 
      scale_fill_manual(values = domaincol) +
      theme_ipsum(base_family = "sans")
    p5.3
    
    
    return(list("pvs_cell" = root_sprn,"pvs_domain" = root_spr_modn,
                "pvs_cell_count" = root_spr,"pvs_domain_count" = root_spr_mod,
                "figure_pvs_cell_pie" = p5.1,
                "figure_pvs_cell_clonesize" = p5.1.1,
                "figure_pvs_cell_disper" = p5.1.2,
                # "figure_pvs_cell_with" = p5.1.t,
                "figure_pvs_domain_pie" = p5.2,
                "figure_pvs_domain_clonesize" = p5.2.1,
                "figure_pvs_domain_disper" = p5.2.2,
                # "figure_pvs_domain" = p5.2.t,
                "figure_pvs_domain_box" = p5.3))
    
  }
  DisperHeatmap = function(pvs){
    ardm = pvs$pvs_domain
    ardm = ardm[-c(1:7,ncol(ardm))]
    ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
    mx = cor(ardm,method = "pearson")
    # pdf("final_result/xdraft/main_figure_2_25/E15_domain_pearson_heatmap.pdf",width = 6,height = 5)
    
    p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
                 col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
                 column_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
                 row_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
                 })
    pvs$figure_domain_heatmap = p1
    return(pvs)
    
  }
  CloneDomainStatPlot = function(pvs){
    domainst = pvs$pvs_domain_count
    domainst = domainst[which(domainst$count > 1),]
    domainst[(ncol(domainst)-4):ncol(domainst)] = domainst[(ncol(domainst)-4):ncol(domainst)]/domainst$count
    # domainst = domainst[which(rowSums(domainst[(ncol(domainst)-4):ncol(domainst)] >= 0.5)>0),]
    
    library(ggpubr)

    #domain density stat
    p1.0 = ggplot(domainst,aes(x = log2(count),color = main.domain)) + geom_density() +
      scale_color_manual(values = domaincol) +
      labs(color = "dominant\ndomain") +
      scale_x_continuous(limits = c(0,max(log2(domainst$count)))) +
      xlab("log2(clone size)") + 
      theme_classic()
    p1.0
    
    #cell complexity
    cellcmpst = pvs$pvs_cell_count
    cellcmpst =  merge(cellcmpst,domainst[c(1,5)],by = "group")
    r = cor(log2(cellcmpst$count),cellcmpst$Hc)
    p1.1.1 = ggplot(cellcmpst,aes(x = log2(count),y = Hc,color = main.domain)) + geom_point() +
      scale_color_manual(values = domaincol) + 
      xlab("log2(clone size)") + ylab("cell type complexity")+
      labs(color = "dominant ndomain") + 
      annotate(geom="text",label = paste0("r = ",round(r,2)),y = max(cellcmpst$Hc),x = 1.5,color = "black") +
      ggtitle("cell type complexity") +
      theme_classic() + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
      
    p1.1.1
    
    #domain max percentage
    domainprop = pvs$pvs_domain
    domainprop = domainprop[which(domainprop$count > 1),]
    domainprop = domainprop[c(1,5,8:12)]
    domainprop[-c(1,2)] = domainprop[-c(1,2)]/rowSums(domainprop[-c(1,2)])
    domainprop = melt(domainprop)
    domainprop = domainprop %>% group_by(group) %>% summarise(main.domain = unique(main.domain),max = max(value))
    domainprop = domainprop[which(domainprop$max != 1),]
    domainprop = merge(domainprop,pvs$pvs_domain[c(1,4)])
    # domainprop = domainprop[which(domainprop$max > 0.5),]
    r = cor(domainprop$Hb,domainprop$max)
    p1.1.2 = ggplot(domainprop,aes(x = Hb,y = max,color = main.domain)) + geom_point() +
      scale_color_manual(values = domaincol) + xlab("Dispersion score") + ylab("max percentage")+
      labs(color = "dominant ndomain") + theme_classic() + 
      annotate(geom="text",label = paste0("r = ",round(r,2)),y = max(domainprop$max),x = 2,
               color = "black") +
      ggtitle("dominant domain percentage") +
      theme_classic() + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
    p1.1.2
    
    #clone size vs dispersion score
    domainsc = domainst
    domainsc = domainsc[which(domainsc$count>1),]
    r = cor(log2(domainsc$count),domainsc$Hb)
    p1.1.3 = ggplot(domainsc,aes(x = log2(count),y = Hb,color = main.domain)) + geom_point() +
      scale_color_manual(values = domaincol) + xlab("log2(clone size)") + ylab("Dispersion score")+
      labs(color = "dominant ndomain") + theme_classic() + 
      annotate(geom="text",label = paste0("r = ",round(r,2)),y = max(domainsc$Hb),x = 1.5,
               color = "black") +
      ggtitle("dispersion potential") +
      theme_classic() + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
    p1.1.3
    
    
    p1.1 = ggarrange( p1.1.1 ,p1.1.2 ,p1.1.3,nrow = 1, common.legend = TRUE, legend="bottom")
    p1.1
    
    
    #domain clone size boxplot& mann kendall test
    library(trend)
    # domainkend = pvs$pvs_domain
    domainkend = domainst
    domainkendst = domainkend %>% group_by(main.domain) %>% summarise(nc = median(count))
    domainkendst$main.domain = factor(domainkendst$main.domain,levels = colnames(domainkend)[8:12])
    domainkendst = domainkendst[order(domainkendst$main.domain),]
    mkt = mk.test(domainkendst$nc)
    mkt$p.value
    domainkend$main.domain = factor(domainkend$main.domain,levels = colnames(domainkend)[8:12])
    p1.2 = ggplot(domainkend,aes(y = log2(count), x = main.domain, fill = main.domain)) + 
      geom_boxplot(width = 0.5,color = "black") + 
      geom_jitter(width = 0.1,size = 1)+
      stat_compare_means(comparisons = list(c("BM","FP"),c("BI","BM"),c("BL","BI"), c("AL","BL")),
                         label = "p.signif",
                         method = "t.test") +
      scale_fill_manual(values = domaincol) + 
      ylab("log2(clone size)") + xlab("Dominant Domain") + 
      ggtitle(label = paste0("Mann-Kendall trend test:\nz = ",round(unlist(mkt$statistic),4),
                             ", p = ",round(mkt$p.value,4))) +
      theme_classic()+theme(legend.position = "none",plot.title = element_text(hjust = 0.8,size = 10,
                                                                               vjust = -10))
    pvs$figure_domain_basic_stat_list = list("domain_density.f4b" = p1.0,"clone_domain_stat.ef6f" = p1.1,
                                          "domain_mann_kendall_test.ef6e" = p1.2)
    return(pvs)
    
  }
  CaldmUpalue = function(pvs){
    library(UpSetR)
    dmup = pvs$pvs_domain
    dmup = dmup[c(1,8:12)]
    tmp = as.matrix(dmup[-1])
    tmp[which(tmp>0)] = 1
    dmup[-1] = tmp
    rownames(dmup) = dmup$group;dmup = dmup[-1]
    
    dmup.up = as.matrix(dmup)
    dmup.up[which(dmup.up>0)] = T
    dmup.up[which(dmup.up==0)] = F
    
    library(tidyr)
    pv = colSums(dmup)/nrow(dmup)
    pvmx = crossing(AL = 0:1, BL = 0:1, BI = 0:1, BM = 0:1, FP = 0:1)
    pvmx$pe = 0
    pvmx$count = 0
    pvmx$expect = 0
    pvmx$pvalue = 0
    for (i in 2:nrow(pvmx)) {
      line = pvmx[i,1:5]
      linei = which(line>0)
      pe = prod(pv[linei],1-pv[-linei])
      vect = dmup[which(rowSums(dmup[linei] > 0) == length(linei) & rowSums(dmup[-linei] > 0) == 0),]
      pvmx$pe[i] = pe
      pvmx$count[i] = nrow(vect)
    }
    
    pvmx$pe = pvmx$pe/sum(pvmx$pe)
    pvmx$expect = pvmx$pe*nrow(dmup)
    
    for (i in 2:nrow(pvmx)) {
      tmp = binom.test(pvmx$count[i], nrow(dmup), p = pvmx$pe[i], 
                       alternative = "two.sided",conf.level= 0.95)
      pvmx$pvalue[i] = tmp$p.value
    }
    pvmx$FDP = p.adjust(pvmx$pvalue)
    pvmx$elementnum = rowSums(pvmx[1:5])
    pvmx = pvmx[order(pvmx$elementnum,-pvmx$count),]
    pvmx = pvmx[which(pvmx$count>0),]
    upcol = NULL
    for (i in 1:nrow(pvmx)) {
      line = pvmx[i,]
      if(line$elementnum == 1){
        coli = "#E60011"
      }else if(line$elementnum == 2){
        tmpid = which(line[1:5]>0)
        if(tmpid[2]-tmpid[1] == 1){
          coli = "#00A0E9"
        }else{
          coli = "black"
        }
      }else{
        coli = "black"
      }
      
      if(line$FDP < 0.05 & line$count > line$expect){
        colsigi = "red"
      }else{
        colsigi = "black"
      }
      upcol = rbind(upcol,data.frame("colsig" = colsigi,"colbar" = coli))
    }
    
    p1.0 = upset(as.data.frame(dmup.up),nsets = 5, sets = colnames(dmup.up),keep.order = T,
                 # matrix.color = rep(upcol$colsig,2),
                 main.bar.color = upcol$colbar)
    p1.0
    
    pvmx$sig = "N"
    pvmx[which(pvmx$FDP<0 & pvmx$count > pvmx$expect),"sig"] = "Y"
    
    pvs$upset_plot = list("figure" = p1.0,"pvalue" = pvmx)
    
    return(pvs)
  }
  # ggsave(pvs11$figure_pvs,filename = "final_result/xdraft/main_figure_2_25/E11_domain_disper_cmb.pdf",width = 8,height = 6)
  # ggsave(pvs15$figure_pvs,filename = "final_result/xdraft/main_figure_2_25/E15_domain_disper_cmb.pdf",width = 8,height = 6)
  domainls = list(domain11,domain11,domain11,domain11,domain15,domain15,domain15,domain15)
  ddisperls = list()
  # for (i in 1:length(sample_todo)) {
  for (i in 1:8) {
    print(i)
    arraydf = datals$arraydata_bl[[sample_todo[i]]]
    pvs = PseudoDisper(arraydf,domainls[[i]])
    pvs = DisperHeatmap(pvs)
    pvs = CloneDomainStatPlot(pvs)
    pvs = CaldmUpalue(pvs)
    ddisperls[[i]] = pvs
    names(ddisperls)[i] = sample_todo[i]
  }
  
  dir.create("final_result/CREST_final_analysis/Clone_domain_E11_E15_analysis")
  qsave(ddisperls,file = "final_result/CREST_final_analysis/Clone_domain_E11_E15_analysis/ddisperls.qs")
  
  domaincol = c("AL" = "#D0D520","BL" = "#EEAF31","BI" = "#B71724","BM" = "#571F7C","FP" = "#192C78")
  ScaleDisperPie = function(ddisperlsi,width,height)
  {
    pnew = ggplot() + 
      geom_scatterpie(data = ddisperlsi$pvs_domain,aes(x=mean.model*width, y=Hb, group=group,
                                                       r = log10(count)/15,alpha = 0.5),color=NA,
                      cols=colnames(ddisperlsi$pvs_domain)[-c(1:7,13)]) + 
      coord_fixed() +
      scale_fill_manual(values = domaincol) +
      scale_y_continuous(breaks = c(0,1,2,3)*height,labels = c(0,1,2,3)) + 
      scale_x_continuous(breaks = c(1,2,3,4,5)*width,labels = c(1,2,3,4,5)) + 
      xlab("") + ylab("Dispersion score(DS)") + 
      theme_classic() + theme(legend.position="none")
    return(pnew)
  }
  ScaleDisperPie(ddisperlsi,1.5,1)
  
  ddisperls$E11 = CloneDomainStatPlot(ddisperls$E11)
  ddisperls$E15 = CloneDomainStatPlot(ddisperls$E15)
  ddisperls$E11$figure_domain_basic_stat_list$domain_mann_kendall_test.ef6e
  ddisperls$E15_rep1$figure_domain_basic_stat_list$domain_mann_kendall_test.ef6e
  ddisperls$E15$figure_domain_basic_stat_list$domain_mann_kendall_test.ef6e
  ddisperls$E15$figure_domain_basic_stat_list$clone_domain_stat.ef6f
  ddisperls$E15_rep3$figure_domain_basic_stat_list$clone_domain_stat.ef6f
  CloneDomainStatPlot
  ddisperls$E11_rep1$figure_pvs_domain
  ddisperls$E11_rep2$figure_pvs_domain
  ddisperls$E11$figure_pvs_domain
  ddisperls$E15$figure_pvs_domain
  datals$figure$DomainDispersion_f4 = ddisperls
  qsave(ddisperls,"final_result/CREST_final_analysis/Clone_domain_E11_E15_analysis/ddisperls.qs")
  qsave(datals,"final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_17_edit.qs")
  
  
}

#fig6. OGN analysis
{
  datals = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_19_edit.qs")
  #
  myorderpE11 = c("NbGABA0","NPBL","NbGABABL","NbGLUAL1",
                  "Rgl1","NbFP",
                  "NPBM","NbBM0","NbBM1",
                  "NbGABABI","NbGLUAL2","OMTN","GLUAL","Peric")
  # myorderOGN = c("GABABL1","GABABL4","GABABL2","GLUAL1","GABABI","GLUMHB2",
  #                "GLUBM2","GLUBM1","NbBM1",
  #                "NbDA","DA","NbGLU","GLUFP","NbFP",
  #                "NbGABABI","GABA_unk","Rgl1")
  myorderOGN = c('GABABL1','GABABL4','GABABL2','GLUAL1','GABABI','GLUMHB2',
                 'Rgl2BP',
               'NbGABABI','Rgl3','NbDA','NbGLU','DA',
               'GABA_unk',
               'RglUnk',
               'GLUBM2','GLUBM1','GLUFP','NbBM1','NbFP')
  
  progcell = c("E11_NPBL","E11_NPBM","E11_Rgl1")
  
  domain11
  domain15
  
  ht_cluster = data.frame("celltype" = c("E11_NbGABA0","E11_NPBL","E11_NbGABABL","E11_NbGLUAL1",
                                         "E11_NPBM","E11_NbBM0","E11_NbBM1",
                                         "E11_Rgl1","E11_NbFP",
                                         "OGN_GABABL1","OGN_GABABL4","OGN_GABABL2","OGN_GLUAL1","OGN_GABABI","OGN_GLUMHB2",
                                         "OGN_GLUBM2","OGN_GLUBM1","OGN_NbBM1","OGN_Rgl2BP",
                                         "OGN_NbDA","OGN_DA","OGN_NbGLU","OGN_GLUFP",'OGN_RglUnk',"OGN_NbFP"),
                          "domain" = rep(c("c1","c2","c3","c1","c2","c3"),c(4,3,2,6,4,6)))
  
  
  domainp11 = domain11
  domainogn = domain15
  domainp11$Cell.type = paste0("E11_",domainp11$Cell.type)
  domainogn$Cell.type = paste0("OGN_",domainogn$Cell.type)
  ht_cluster = rbind(domainp11,domainogn)
  ht_cluster = ht_cluster[order(ht_cluster$domain),]
  colnames(ht_cluster)[2] = "celltype"
  
  myrow_order = c("E11_NbGABA0","E11_NPBL","E11_NbGABABL","E11_NbGLUAL1",
                  "E11_Rgl1","E11_NbFP",
                  "E11_NPBM","E11_NbBM0","E11_NbBM1",
                  "E11_NbGABABI","E11_NbGLUAL2","E11_OMTN","E11_GLUAL","E11_Peric")
  mycol_order = c("OGN_GABABL1","OGN_GABABL4","OGN_GABABL2","OGN_GLUAL1","OGN_GABABI","OGN_GLUMHB2",
                  "OGN_GLUBM2","OGN_GLUBM1","OGN_NbBM1",
                  "OGN_Rgl2BP",
                  "OGN_NbDA","OGN_DA","OGN_NbGLU","OGN_GLUFP",'OGN_RglUnk',
                  "OGN_NbFP", 
                  "OGN_NbGABABI","OGN_GABA_unk","OGN_Rgl3")
  mycol_order = c(orderE15[which(orderE15 %in% myorderOGN)],"RglUnk","GABA_unk")
  mycol_order = paste0("OGN_",mycol_order)
  
  #Heatmap
  OGNHeatmap = function(e11arydf,ognarydf,title,myorder = T)
  {
    count1 = e11arydf
    count2 = ognarydf
    count1 = count1[,c("pattern",myorderpE11[which(myorderpE11 %in% colnames(count1))])]
    count2 = count2[,c("pattern",myorderOGN[which(myorderOGN %in% colnames(count2))])]
    colnames(count1)[-1] = paste0("E11_",colnames(count1)[-1])
    colnames(count2)[-1] = paste0("OGN_",colnames(count2)[-1])
    e11divcms = merge(count1,count2,by = "pattern")
    
    e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$pattern;e11divcms.pd = e11divcms.pd[-1]
    e11divcms.pd = as.matrix(e11divcms.pd)
    
    split = factor(rep(c("E11","OGN"),c(ncol(count1[-1]),ncol(count2[-1]))),levels = c("E11","OGN"))
    names(split) = colnames(e11divcms.pd)
    
    e11divcms.pd = e11divcms.pd[do.call("order", as.data.frame(-e11divcms.pd)), ]
    
    # pdf("final_result/xE11_DIV7/E11-DIV7_array_cmp_heatmap.pdf",width = 6,height = 10)
    p1.0_ht = Heatmap(log2(e11divcms.pd+1),show_row_names = F,name = "log2 counts",
            # col = c("#23174C","#376BB4","white"),
            col = colorRampPalette(brewer.pal(9, "Blues"))(30),
            column_gap = unit(5, "mm"),
            cluster_row_slices = FALSE,row_order = rownames(e11divcms.pd),
            column_split = split,column_order = colnames(e11divcms.pd))
    p1.0_ht
    #cmpare heatmap
    #方法1
    cell_norm = e11divcms.pd[which(rowSums(e11divcms.pd>0)>1),]
    cell_norm = cell_norm[,which(colSums(cell_norm)>0)]
    cell_norm = log(cell_norm+1)
    cell_norm = t(apply(cell_norm, 1, function(x) {x/sum(x)}))
    # e11divcmf.mx = cor(as.matrix(cell_norm),method = "spearman")
    # e11divcmf.mx = e11divcmf.mx[colnames(count1)[-1][which(colnames(count1)[-1] %in% rownames(mx))],
    #                             colnames(count2)[-1][which(colnames(count2)[-1] %in% colnames(mx))]]
    mx = cor(cell_norm)
    dd = Dist(t(cell_norm), method = "abspearson")
    mx = as.matrix(1 - dd)
    for (i in 1:nrow(mx)) {
      mx[i,i] = 1
    }
    mx = mx[colnames(count1)[-1][which(colnames(count1)[-1] %in% rownames(mx))],
            colnames(count2)[-1][which(colnames(count2)[-1] %in% colnames(mx))]]
    # dd = dd[is.na(dd)] = 0
    
    p1.1_htcor = Heatmap(mx,heatmap_legend_param = list(title = "abspearson"),row_title = title,
            clustering_method_columns = "complete",clustering_method_rows = "complete",
            col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4))
            )
    p1.1_htcor
    
    if(myorder){
      p1.1_htcor_order = Heatmap(mx,heatmap_legend_param = list(title = "abspearson"),row_title = title,
                                 row_order = myrow_order[which(myrow_order%in% rownames(mx))],
                                 column_order = mycol_order[which(mycol_order%in% colnames(mx))],
                                 col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4))
                                 # col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
      )
      p1.1_htcor_order
    }else{
      p1.1_htcor_order = NULL
    }
    
    return(list("cellmx.e" = count1,"cellmx.l" = count2,
                "clonedf" = e11divcms,"e11_ogn_clone_heatmap_figure" = p1.0_ht,
                "e11_ogn_correlation_heatmap_figure" = p1.1_htcor,
                "e11_ogn_correlation_heatmap_order_figure" = p1.1_htcor_order))
  }
  
  #pandaE11-OGN alluv figure
  E11OGNAlluv = function(e11ognlsi,celllevel = NULL)
  {
    cellmx.e = e11ognlsi$cellmx.e
    cellmx.l = e11ognlsi$cellmx.l
    cellft = e11ognlsi$clonedf
    earlyid = colnames(cellmx.e)[-1]
    earlyid = earlyid[which(earlyid %in% ht_cluster$celltype)]
    lateid = colnames(cellmx.l)[-1]
    lateid = lateid[which(lateid %in% ht_cluster$celltype)]
    rownames(cellft) = cellft$pattern;cellft = cellft[-1]
    cellft = cellft[colnames(cellft) %in% ht_cluster$celltype]
    cellft = cellft[which(rowSums(cellft>0) < ncol(cellft)/3),]
    if(length(celllevel) == 0){
      celllevel = ht_cluster$celltype
    }
    
    
    DIVAllu = function(cellft,celllevel,ht_cluster,progcell){
      phasest = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                           "phase2" = rep(lateid,length(earlyid)), "Freq.from" = 0, "Freq.to" = 0)
      for (i in 1:nrow(cellft)) {
        line = cellft[i,]
        line.e = unlist(line[,which(colnames(line) %in% earlyid)])
        line.e = line.e[which(line.e>0)]
        line.l = unlist(line[,which(colnames(line) %in% lateid)])
        line.l = line.l[which(line.l>0)]
        for (j in 1:length(line.e)) {
          for (k in 1:length(line.l)) {
            seg.from = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.from"]
            seg.to = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.to"]
            phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.from"] = seg.from + line.e[j]
            phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq.to"] = seg.to + line.l[k]
          }
        }
      }
      # phasest = phasest[which((phasest$Freq.to + phasest$Freq.from) > 0),]
      #method2
      rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
        vec <- rnorm(N, M/N, sd)
        if (abs(sum(vec)) < 0.01) vec <- vec + 1
        vec <- round(vec / sum(vec) * M)
        deviation <- M - sum(vec)
        for (. in seq_len(abs(deviation))) {
          vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
        }
        if (pos.only) while (any(vec < 0)) {
          negs <- vec < 0
          pos  <- vec > 0
          vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
          vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
        }
        vec
      }
      phasestml = list()
      for (t in 1:100) {
        print(t)
        cellftm = cellft
        for (i in 2:ncol(cellftm)) {
          cellftm[i] = rand_vect(nrow(cellft),sum(cellft[i]),sd = sd(unlist(cellft[i])))
        }
        phasestm = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                              "phase2" = rep(lateid,length(earlyid)), "Freq.from" = 0, "Freq.to" = 0)
        for (i in 1:nrow(cellftm)) {
          line = cellftm[i,]
          line.e = unlist(line[,which(colnames(line) %in% earlyid)])
          line.e = line.e[which(line.e>0)]
          line.l = unlist(line[,which(colnames(line) %in% lateid)])
          line.l = line.l[which(line.l>0)]
          for (j in 1:length(line.e)) {
            for (k in 1:length(line.l)) {
              seg.from = phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.from"]
              seg.to = phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.to"]
              phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.from"] = seg.from + line.e[j]
              phasestm[which(phasestm$phase1 == names(line.e)[j] & phasestm$phase2 == names(line.l)[k]),"Freq.to"] = seg.to + line.l[k]
            }
          }
        }
        
        phasestml[[t]] = phasestm
      }
      phasestmle = phasestmll = list()
      for (j in 1:nrow(phasestml[[1]])) {
        phasestmle[[j]] = phasestml[[1]][j,3]
        phasestmll[[j]] = phasestml[[1]][j,4]
      }
      
      phasestmle = phasestmll = list()
      for (j in 1:nrow(phasestml[[1]])) {
        phasestmle[[j]] = phasestml[[1]][j,3]
        phasestmll[[j]] = phasestml[[1]][j,4]
      }
      
      for (i in 2:length(phasestml)) {
        line = phasestml[[i]]
        for (j in 1:nrow(line)) {
          phasestmle[[j]] = c(phasestmle[[j]],line[j,3])
          phasestmll[[j]] = c(phasestmll[[j]],line[j,4])
        }
      }
      
      
      #filter--------
      phasestf = NULL
      for (i in 1:nrow(phasest)) {
        line = phasest[i,]
        if(line$Freq.from > quantile(phasestmle[[i]],probs = 0.1) | 
           line$Freq.to > quantile(phasestmll[[i]],probs = 0.1)){
          phasestf = rbind(phasestf,line)
        }
      }
      
      # phasest = cellallud7$phasest
      # phasestmle = cellallud7$phasestmle
      # phasestmll = cellallud7$phasestmll
      
      PhasePlot = function(phasestf,domain_df,celllevel){
        # phasestf = phasest
        phasest_lod = to_lodes_form(phasestf,
                                    key = "phase", value = "celltype", id = "alluvium",
                                    axes = 1:2)
        phasest_lod$Freq = 0
        freqfrom = phasest_lod[which(phasest_lod$phase == "phase1"),"Freq.from"]
        freqto = phasest_lod[which(phasest_lod$phase == "phase2"),"Freq.to"]
        phasest_lod[which(phasest_lod$phase == "phase1"),"Freq"] = freqfrom/sum(freqfrom)
        phasest_lod[which(phasest_lod$phase == "phase2"),"Freq"] = freqto/sum(freqto)
        
        phasest_lod = merge(phasest_lod,domain_df,
                            by = "celltype")
        phasest_lod = phasest_lod[order(phasest_lod$domain),]
        phasest_lod$celltype = factor(phasest_lod$celltype,
                                      levels = unique(celllevel))
        
        
        pf1 = ggplot(phasest_lod,
                     aes(x = phase, stratum = celltype, alluvium = alluvium,
                         y = Freq,
                         fill = domain,label = celltype),color = "black") +
          scale_x_discrete(expand = c(.1, .1)) +
          scale_fill_manual(values = domaincol) +
          # geom_flow(aes(color = celltype)) +
          geom_flow(aes()) +
          geom_stratum(alpha = .5) +
          geom_text(stat = "stratum", size = 3) + theme_void()
        pf1
        return(pf1)
        
      }
      # phasest1 = cellallud7$prog$phasest
      # phasest2 = cellallud7$other$phasest
      phasest1 = phasest[which(phasest$phase1 %in% progcell),]
      phasest2 = phasest[which(!phasest$phase1 %in% progcell),]
      phasestf1 = phasestf[which(phasestf$phase1 %in% progcell),]
      phasestf2 = phasestf[which(!phasestf$phase1 %in% progcell),]
      pf1 = PhasePlot(phasestf1,ht_cluster,celllevel)
      pf2 = PhasePlot(phasestf2,ht_cluster,celllevel)
      pft = PhasePlot(phasestf,ht_cluster,celllevel)
      
      cellallud71 = list("cellft" = cellft,"phasest" = phasest1,"phasestf" = phasestf1,"phasestmle" = phasestmle[which(phasest$phase1 %in% progcell)],
                         "phasestmll" = phasestmll[which(phasest$phase1 %in% progcell)],"alluv_figure" = pf1)
      cellallud72 = list("cellft" = cellft,"phasest" = phasest2,"phasestf" = phasestf2,"phasestmle" = phasestmle[which(!phasest$phase1 %in% progcell)],
                         "phasestmll" = phasestmll[which(!phasest$phase1 %in% progcell)],"alluv_figure" = pf2)
      cellallud7t = list("cellft" = cellft,"phasest" = phasest,"phasestf" = phasestf,"phasestmle" = phasestmle,
                         "phasestmll" = phasestmll,"alluv_figure" = pft)
      
      cellallud7 = list("prog" = cellallud71,"other" = cellallud72,"total" = cellallud7t)
      
      return(cellallud7)
      
    }
    
    cellallud7 = DIVAllu(cellft,celllevel,ht_cluster,progcell)
    cellallud7$prog = CalAlluPvalue(cellallud7$prog)
    cellallud7$other = CalAlluPvalue(cellallud7$other)
    cellallud7$total = CalAlluPvalue(cellallud7$total)
    e11ognlsi$Allufigurels = cellallud7
    return(e11ognlsi)
  }
  
  celllevel = unique(ht_cluster$celltype)
  # celllevel = celllevel[c(1:4,8,5:7,9:23)]
  
  pe11id = c("pandaE11_rep1","pandaE11_rep2","pandaE11_rep3","pandaE11")
  ognid = c("OGN_rep1","OGN_rep2","OGN_rep3","OGN")
  e11ognls = list()
  for (i in 1:length(pe11id)) {
    e11arydf = datals$arraydata_bl[[pe11id[i]]]
    ognarydf = datals$arraydata_bl[[ognid[i]]]
    e11arydf = dcast(e11arydf, pattern~Cell.type)
    ognarydf = dcast(ognarydf, pattern~Cell.type)
    e11arydf$pattern = substr(e11arydf$pattern,10,nchar(e11arydf$pattern))
    ognarydf$pattern = substr(ognarydf$pattern,5,nchar(ognarydf$pattern))
    
    e11arydfv1 = datals$v1[[pe11id[i]]]$result_withcellan
    ognarydfv1 = datals$v1[[ognid[i]]]$result_withcellan
    e11arydfv1 = dcast(e11arydfv1, pattern~Cell.type)
    ognarydfv1 = dcast(ognarydfv1, pattern~Cell.type)

    
    e11arydfv2 = datals$v2[[pe11id[i]]]$result_withcellan
    ognarydfv2 = datals$v2[[ognid[i]]]$result_withcellan
    e11arydfv2 = dcast(e11arydfv2, pattern~Cell.type)
    ognarydfv2 = dcast(ognarydfv2, pattern~Cell.type)
    
    e11ognlsiv1 = OGNHeatmap(e11arydfv1,ognarydfv1,paste0(pe11id[i],"-",ognid[i]))
    e11ognlsiv2 = OGNHeatmap(e11arydfv2,ognarydfv2,paste0(pe11id[i],"-",ognid[i]))
    
    e11ognlsi = OGNHeatmap(e11arydf,ognarydf,paste0(pe11id[i],"-",ognid[i]))
    e11ognlsi = E11OGNAlluv(e11ognlsi,celllevel)
    e11ognlsi$v1_heatmap = e11ognlsiv1
    e11ognlsi$v2_heatmap = e11ognlsiv2
    
    e11ognls[[i]] = e11ognlsi
    
    # e11ognls[[i]]$e11_ogn_correlation_heatmap_figure = e11ognlsi$e11_ogn_correlation_heatmap_figure
    # e11ognls[[i]]$e11_ogn_correlation_heatmap_order_figure = e11ognlsi$e11_ogn_correlation_heatmap_order_figure
  }
  # e11ognls = datals$figure$pandaE11_OGNls_f6
  names(e11ognls) = paste0(pe11id,"-",ognid)
  e11ognls$`pandaE11-OGN`$Allufigurels$total$alluv_figure
  e11ognls$`pandaE11_rep1-OGN_rep1`$v1_heatmap$e11_ogn_correlation_heatmap_figure
  e11ognls$`pandaE11_rep2-OGN_rep2`$e11_ogn_correlation_heatmap_order_figure
  e11ognls$`pandaE11_rep3-OGN_rep3`$e11_ogn_correlation_heatmap_order_figure
  e11ognls$`pandaE11-OGN`$e11_ogn_correlation_heatmap_figure
  e11ognls$`pandaE11_rep1-OGN_rep1`$Allufigurels$prog$alluv_figure
  e11ognls$`pandaE11_rep2-OGN_rep2`$Allufigurels$prog$alluv_figure
  e11ognls$`pandaE11_rep3-OGN_rep3`$Allufigurels$prog$alluv_figure
  
  dir.create("final_result/CREST_final_analysis/OGN_analysis")
  qsave(e11ognls,file = "final_result/CREST_final_analysis/OGN_analysis/pandaE11_OGNls_1_9.qs")
  datals$figure$pandaE11_OGNls_f6 = e11ognls
  qsave(datals,"final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_28_edit.qs")
}

#fig4. E11 E15 clone cluster edit
{
  e15cl = readRDS("final_result/CREST_final_analysis/annotation_20221207/e15.5 clone cluster/snapE15_clonecluster_20221209_cellcounts_cluster.rds")
  e11cl = readRDS("final_result/CREST_final_analysis/annotation_20221207/snapE11_20221214_cellcounts_cluster.rds")
  
  cellordern = c('GLUAL1','GLUAL2','GLUAL3','GLUAL4','Rgl2AL','GABAAL','NbGABAAL',
                 'Rgl3','NbGABABL','GABABL1','GABABL2','GABABL3','GABABL4','Rgl2BP',
                 'NbGABABI','GABABI','NbBM1','GLUBM1','GLUBM2',
                 'NbDA','DA','NbFP','Rgl1','NbGLU','GLUFP',
                 'GLUMHB1','GLUMHB2','GABAMHB','GLUMHB3','GLUMHB4')
  arraydf = datals$arraydata_bl$E15
  CloneClusterAppendUni = function(e15cl,arraydf){
    e15cc = dcast(arraydf,pattern~Cell.type)
    rownames(e15cc) = e15cc$pattern;e15cc = e15cc[-1]
    e15cc = e15cc[rowSums(e15cc)>1,]
    e15cc = as.data.frame(t(t(as.matrix(e15cc))/colSums(e15cc)))
    e15cc = e15cc[,names(e15cl$cluster.fig@column_order)]
    e15cc = e15cc[do.call("order", -e15cc),]
    adid = rownames(e15cc)[which(!rownames(e15cc) %in% names(e15cl$cluster))]
    tmp = e15cc[which(rownames(e15cc) %in% adid),]
    # e15_clu$cellcounts = e15_clu$cellcounts[which(rownames(e15_clu$cellcounts) %in% names(nsplit)),]
    tmpord = e15cl$cluster.fig@column_order[order(e15cl$cluster.fig@column_order)]
    e15add.hmap <- Heatmap(as.matrix(tmp[names(tmpord)]), name = "normlize enrichment score",
                           row_order = rownames(tmp),
                           # col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                           show_row_names = F,column_order = names(e15cl$cluster.fig@column_order))
    e15cl$cluster.fig@column_order
    e15clt = e15cl$cluster.fig %v% e15add.hmap
    return(e15clt)
  }
  
  e15clt = CloneClusterAppendUni(e15cl,datals$arraydata_bl$E15)
  e11clt = CloneClusterAppendUni(e11cl,datals$arraydata_bl$E11)
  
  datals$figure$CloneCluster_f4 = list("E11" = list("cc" = e11cl,"cc_withuni" = e11clt),
                                       "E15" = list("cc" = e15cl,"cc_withuni" = e15clt))
  qsave(e11clt,file = "final_result/CREST_final_analysis/E11_E15_basic_stat/E11_clone_cluster_withuni_qs")
  
}

#sfig. replicates analysis----------
{
  #replicates correlation stat
  {
    maxtrixls = list("E11" = list(melt(datals$figure$heatmap$E11_rep1$abspearson@matrix),
                                  melt(datals$figure$heatmap$E11_rep2$abspearson@matrix),
                                  melt(datals$figure$heatmap$E11_rep3$abspearson@matrix)),
                     # "pandaE11" = list(melt(datals$figure$heatmap$pandaE11_rep1$abspearson@matrix),
                     #                   melt(datals$figure$heatmap$pandaE11_rep2$abspearson@matrix),
                     #                   melt(datals$figure$heatmap$pandaE11_rep3$abspearson@matrix)),
                     "E15" = list(melt(datals$figure$heatmap$E15_rep1$abspearson@matrix),
                                  melt(datals$figure$heatmap$E15_rep2$abspearson@matrix),
                                  melt(datals$figure$heatmap$E15_rep3$abspearson@matrix)),
                     # "OGN" = list(melt(datals$figure$heatmap$OGN_rep1$abspearson@matrix),
                     #              melt(datals$figure$heatmap$OGN_rep2$abspearson@matrix),
                     #              melt(datals$figure$heatmap$OGN_rep3$abspearson@matrix)),
                     "pandaE11_OGN" = list(melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep1-OGN_rep1`$e11_ogn_correlation_heatmap_figure@matrix),
                                           melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep2-OGN_rep2`$e11_ogn_correlation_heatmap_figure@matrix),
                                           melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep3-OGN_rep3`$e11_ogn_correlation_heatmap_figure@matrix)))
    
    maxtrixlsv1 = list("E11" = list(melt(datals$figure$heatmap_order_v1$E11_rep1$abspearson@matrix),
                                  melt(datals$figure$heatmap_order_v1$E11_rep2$abspearson@matrix),
                                  melt(datals$figure$heatmap_order_v1$E11_rep3$abspearson@matrix)),
                     "E15" = list(melt(datals$figure$heatmap_order_v1$E15_rep1$abspearson@matrix),
                                  melt(datals$figure$heatmap_order_v1$E15_rep2$abspearson@matrix),
                                  melt(datals$figure$heatmap_order_v1$E15_rep3$abspearson@matrix)),
                     "pandaE11_OGN" = list(melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep1-OGN_rep1`$v1_heatmap$e11_ogn_correlation_heatmap_figure@matrix),
                                           melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep2-OGN_rep2`$v1_heatmap$e11_ogn_correlation_heatmap_figure@matrix),
                                           melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep3-OGN_rep3`$v1_heatmap$e11_ogn_correlation_heatmap_figure@matrix)))
    
    maxtrixlsv2 = list("E11" = list(melt(datals$figure$heatmap_order_v2$E11_rep1$abspearson@matrix),
                                    melt(datals$figure$heatmap_order_v2$E11_rep2$abspearson@matrix),
                                    melt(datals$figure$heatmap_order_v2$E11_rep3$abspearson@matrix)),
                       "E15" = list(melt(datals$figure$heatmap_order_v2$E15_rep1$abspearson@matrix),
                                    melt(datals$figure$heatmap_order_v2$E15_rep2$abspearson@matrix),
                                    melt(datals$figure$heatmap_order_v2$E15_rep3$abspearson@matrix)),
                       "pandaE11_OGN" = list(melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep1-OGN_rep1`$v2_heatmap$e11_ogn_correlation_heatmap_figure@matrix),
                                             melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep2-OGN_rep2`$v2_heatmap$e11_ogn_correlation_heatmap_figure@matrix),
                                             melt(datals$figure$pandaE11_OGNls_f6$`pandaE11_rep3-OGN_rep3`$v2_heatmap$e11_ogn_correlation_heatmap_figure@matrix)))
    
    
  }
  
  
  RepCorCal = function(maxtrixls){
    mxt = Reduce(function(x,y) merge(x,y,by = c("Var1","Var2")),maxtrixls)
    colnames(mxt)[3:ncol(mxt)] = paste0("rep",1:(ncol(mxt)-2))
    mxt[is.na(mxt)] = 0
    mxt$Var1 = as.character(mxt$Var1)
    mxt$Var2 = as.character(mxt$Var2)
    mxt = mxt[which(mxt$Var1 != mxt$Var2),]
    repcorl = NULL
    for (i in 3:ncol(mxt)) {
      for (j in (i+1):ncol(mxt)) {
        if(j <= ncol(mxt) & j!=i){
          repcorl = rbind(repcorl,data.frame("group" = paste0(colnames(mxt)[i],"_",colnames(mxt)[j]),
                                             "cor" = cor(mxt[,i],mxt[,j])))
        }
      }
    }
    
    if(nrow(repcorl) == 3){
      repcorl$cocor = sqrt( ( (repcorl$cor[1])^2 + (repcorl$cor[2])^2 + (repcorl$cor[3])^2 ) - 
                              ( 2 * repcorl$cor[1] * repcorl$cor[2] * repcorl$cor[3]) )
      # mxt$rep1 = log(mxt$rep1+0.000001)
      # mxt$rep2 = log(mxt$rep2+0.000001)
      # mxt$rep3 = log(mxt$rep3+0.000001)
      
      p1.0 = ggscatter(mxt,x = "rep1",y = "rep2",size = 0.5) + 
        ggtitle(paste0("cor = " ,round(cor(mxt$rep1,mxt$rep2),3),
                       ", cor test p = ", round(cor.test(mxt$rep1,mxt$rep2)$p.value,3)))
      p1.1 = ggscatter(mxt,x = "rep1",y = "rep3",size = 0.5) + 
        ggtitle(paste0("cor = " ,round(cor(mxt$rep1,mxt$rep3),3),
                       ", cor test p = ", round(cor.test(mxt$rep1,mxt$rep3)$p.value,3)))
      p1.2 = ggscatter(mxt,x = "rep2",y = "rep3",size = 0.5) + 
        ggtitle(paste0("cor = " ,round(cor(mxt$rep2,mxt$rep3),3),
                       ", cor test p = ", round(cor.test(mxt$rep2,mxt$rep3)$p.value,3)))
      pt = ggarrange(p1.0,p1.1,p1.2,nrow = 1)
      
    }else{
      repcorl$cocor = repcorl$cor
      pt = ggscatter(mxt,x = "rep1",y = "rep2",size = 0.5) + ggtitle(paste0("cor = " ,round(cor(mxt$rep1,mxt$rep2),3)))
    }
    return(list("repcordf" = repcorl,"repcorpp" = pt,"mxt" = mxt))
  }
  repcorls = list()
  repcorlsv1 = list()
  repcorlsv2 = list()
  repcordf = NULL
  for (i in 1:length(maxtrixls)) {
    repcorls[[i]] = RepCorCal(maxtrixls[[i]])
    repcorls[[i]]$repcordf$sample = names(repcorls)[i] = names(maxtrixls)[i]
    repcordfi = repcorls[[i]]$repcordf
    mxt = repcorls[[i]]$mxt
    mxt = melt(mxt[-c(1,2)])
    anvlm = lm(value~variable,mxt)
    repcordfi$p.value = anova(anvlm)$`Pr(>F)`[1]
    repcordf = rbind(repcordf,repcordfi)

  }
  
  p1.0 = ggbarplot(repcordf[which(!repcordf$sample %in% c("pandaE11","OGN")),],
                   x = "sample", y = "cor", fill = "grey",width = 0.5,
            add = c("mean_sd", "jitter")) + 
    geom_text(data = repcordf[which(!repcordf$sample %in% c("pandaE11","OGN")),],
              aes(x = sample, y = 1, label = paste0("p = ", round(p.value,3)))) + 
    scale_y_continuous(limits = c(0,1)) +
    scale_color_grey() + 
    xlab("") + ylab("Replicants lineage correlation of different samples")
  p1.0
  
  arrycordf = NULL
  arrycorppls = list()
  for (i in 1:length(datals$figure$heatmap)) {
    cormxid = names(datals$figure$heatmap)[i]
    cormxt = datals$figure$heatmap[[i]]$cormx
    cormxv1 = datals$figure$heatmap_order_v1[[i]]$cormx
    cormxv2 = datals$figure$heatmap_order_v2[[i]]$cormx
    colnames(cormxt)[3:4] = paste0("total.",colnames(cormxt)[3:4])
    colnames(cormxv1)[3:4] = paste0("V1.",colnames(cormxv1)[3:4])
    colnames(cormxv2)[3:4] = paste0("V2.",colnames(cormxv2)[3:4])
    cormxt = merge(cormxt,cormxv1)
    cormxt = merge(cormxt,cormxv2)
    cormxt[is.na(cormxt)] = 0
    
    pt.0 = ggscatter(cormxt,x = "total.corr",y = "V1.corr",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$total.corr,cormxt$V1.corr),3),
                     ", cor test p = ", round(cor.test(cormxt$total.corr,cormxt$V1.corr)$p.value,3))) +
      xlab("V1&V2") + ylab("V1")
    pt.1 = ggscatter(cormxt,x = "total.corr",y = "V2.corr",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$total.corr,cormxt$V2.corr),3),
                     ", cor test p = ", round(cor.test(cormxt$total.corr,cormxt$V2.corr)$p.value,3))) +
      xlab("V1&V2") + ylab("V2")
    pt.2 = ggscatter(cormxt,x = "V1.corr",y = "V2.corr",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$V1.corr,cormxt$V2.corr),3),
                     ", cor test p = ", round(cor.test(cormxt$V1.corr,cormxt$V2.corr)$p.value,3))) +
      xlab("V1") + ylab("V2")
    pt = ggarrange(pt.0,pt.1,pt.2,nrow = 1)
    arrycorppls[[i]] = pt
    names(arrycorppls)[i] = cormxid
    cordfi = data.frame("cor" = c(cor(cormxt$total.corr,cormxt$V1.corr),
                       cor(cormxt$V1.corr,cormxt$V2.corr),
                       cor(cormxt$total.corr,cormxt$V2.corr)),
                       "group" = c("V1:V1&V2","V1:V2","V2:V1&V2"),
                       "sample" = cormxid)
    arrycordf = rbind(arrycordf,cordfi)
    
  }
  arrycorppls$E11_rep1
  arrycordf$sample_group = unlist(lapply(strsplit(arrycordf$sample,"_"),"[[",1))
  for (i in 1:4) {
    cormxid = names(datals$figure$pandaE11_OGNls_f6)[i]
    cormxt = datals$figure$pandaE11_OGNls_f6[[i]]$e11_ogn_correlation_heatmap_figure@matrix
    cormxv1 = datals$figure$pandaE11_OGNls_f6[[i]]$v1_heatmap$e11_ogn_correlation_heatmap_figure@matrix
    cormxv2 = datals$figure$pandaE11_OGNls_f6[[i]]$v2_heatmap$e11_ogn_correlation_heatmap_figure@matrix
    cormxt = melt(cormxt)
    cormxv1 = melt(cormxv1)
    cormxv2 = melt(cormxv2)
    colnames(cormxt)[3] = paste0("total.",colnames(cormxt)[3])
    colnames(cormxv1)[3] = paste0("V1.",colnames(cormxv1)[3])
    colnames(cormxv2)[3] = paste0("V2.",colnames(cormxv2)[3])
    cormxt = merge(cormxt,cormxv1)
    cormxt = merge(cormxt,cormxv2)
    cormxt[is.na(cormxt)] = 0
    
    
    pt.0 = ggscatter(cormxt,x = "total.value",y = "V1.value",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$total.value,cormxt$V1.value),3),
                     ", cor test p = ", round(cor.test(cormxt$total.value,cormxt$V1.value)$p.value,3))) +
      xlab("V1&V2") + ylab("V1")
    pt.1 = ggscatter(cormxt,x = "total.value",y = "V2.value",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$total.value,cormxt$V2.value),3),
                     ", cor test p = ", round(cor.test(cormxt$total.value,cormxt$V2.value)$p.value,3))) +
      xlab("V1&V2") + ylab("V2")
    pt.2 = ggscatter(cormxt,x = "V1.value",y = "V2.value",size = 0.5) + 
      ggtitle(paste0("cor = " ,round(cor(cormxt$V1.value,cormxt$V2.value),3),
                     ", cor test p = ", round(cor.test(cormxt$V1.value,cormxt$V2.value)$p.value,3))) +
      xlab("V1") + ylab("V2")
    pt = ggarrange(pt.0,pt.1,pt.2,nrow = 1)
    arrycorppls[[i+16]] = pt
    names(arrycorppls)[i] = cormxid
    
    cordfi = data.frame("cor" = c(cor(cormxt$total.value,cormxt$V1.value),
                                  cor(cormxt$V1.value,cormxt$V2.value),
                                  cor(cormxt$total.value,cormxt$V2.value)),
                        "group" = c("V1:V1&V2","V1:V2","V2:V1&V2"),
                        "sample" = cormxid,
                        "sample_group" = "pandaE11_OGN")
    arrycordf = rbind(arrycordf,cordfi)
  }
  
  ggexport(arrycorppls[c(13,14,4)],filename = "final_result/CREST_final_analysis/0.qc/V1_V2_lineage_correlation.pdf",
           width = 12,height = 4)
  qsave(arrycorppls,file = "final_result/CREST_final_analysis/0.qc/V1_V2_lineage_correlation.qs")
  arrycordf = arrycordf[which(arrycordf$sample_group != arrycordf$sample),]
  p1.1 = ggbarplot(arrycordf[which(arrycordf$sample %in% c("E11","E15","pandaE11-OGN")),],
                   x = "sample", y = "cor", color = "group",
                   width = 0.5,
                   # add = c("mean_sd", "jitter"),
                   position = position_dodge(0.8)) +
    geom_text(data = arrycordf[which(arrycordf$sample %in% c("E11","E15","pandaE11-OGN")),],
              aes(x = sample,y = cor, label = round(cor,3)),position = position_dodge2(0.8)) +
    scale_color_jama() + 
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") + ylab("only V1 vs only V2 vs V1&V2 lineage correlation")
  p1.1
  
  ggexport(p1.0,filename = "final_result/CREST_final_analysis/rebbutal_analysis/replicates_lineage_correlation_1_9.pdf",
           width = 5,height = 5)
  ggexport(p1.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/array_type_correlation_1_9.pdf",
           width = 6,height = 5)
  
  replicates_cor = list("repcor_pointls" = repcorls,
                        "array_pointls" = arrycorppls,
                        "repcor_barplot" = p1.0,
                        "repcordf" = repcordf,"arraytypecordf" = arrycordf,
                        "arraytypecor_barplot" = p1.1)
  datals$figure$replicates_cor = list("repcor_pointls" = repcorls,"repcor_barplot" = p1.0,
                                      "repcordf" = repcordf,"arraytypecordf" = arrycordf,
                                      "arraytypecor_barplot" = p1.1)
  ognrevisels = list()
  ognrevisels$replicates_cor = replicates_cor
  ognrevisels$pandaE11_OGNls_f6 = datals$figure$pandaE11_OGNls_f6
  
  
  #ref vs e15 cell
  ggexport(replicates_cor$repcor_pointls,filename = "final_result/CREST_final_analysis/rebbutal_analysis/replicates_cor_1_9.pdf",width = 12,height = 4)
  ggexport(replicates_cor$array_pointls,filename = "final_result/CREST_final_analysis/rebbutal_analysis/replicates_cor_array_1_9.pdf",width = 12,height = 4)
  
  qsave(replicates_cor,"final_result/CREST_final_analysis/benchmark/replicates_cor.qs")
  qsave(datals,"final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_17_edit.qs")
  
}

#sfig. E11 OGN sister compare
{
  E11trans = qread("final_result/CREST_final_analysis/annotation_20221207/E11_subclustered_renamed_UMAPswap_20221209.qs")
  e11_clu = readRDS("final_result/CREST_final_analysis/annotation_20221207/snapE11_clonecluster_20221209_cellcounts_cluster.rds")
  e11cm = datals$arraydata_bl$E11
  e11cm$Cell.BC = paste0(e11cm$Cell.BC,"-1_",substr(e11cm$sample,8,8))
  e11cm = e11cm[which(e11cm$Cell.type != "Peric"),]
  # DimPlot(E11trans,group.by = "orig.ident")
  cmrid = sample(1:nrow(e11cm),nrow(e11cm)/2)
  e11cmar = e11cm[cmrid,]
  e11cmbr = e11cm[setdiff(1:nrow(e11cm),cmrid),]
  cmtagc = intersect(e11cmar$pattern,e11cmbr$pattern)
  e11cma = dcast(e11cmar,Cell.type~pattern)
  e11cmb = dcast(e11cmbr,Cell.type~pattern)
  rownames(e11cma) = e11cma$Cell.type;e11cma = e11cma[-1]
  rownames(e11cmb) = e11cmb$Cell.type;e11cmb = e11cmb[-1]
  progcell = c("Rgl1","NPBL","NPBM")
  e11rna = as.data.frame(GetAssayData(E11trans,assay = "RNA"))
  e11rna = e11rna[,colnames(e11rna) %in% e11cm$Cell.BC]
  
  #analysis
  {
    #function
    SisterUmapAna = function(funseur,funcma,funcmb,ida,idb,cmtag)
    {
      p5.1 = DimPlot(funseur) | DimPlot(funseur,split.by = "group",group.by = "clonecluster")
      funump = as.data.frame(funseur@reductions$umap@cell.embeddings)
      
      distmat <- as.matrix(dist(funump))    #matrix of Euclidean distances between rows
      diag(distmat) <- Inf                      #remove zeros on diagonal
      distvec = NULL
      for (i in 1:length(cmtag)) {
        distvec = c(distvec,distmat[ida[i],ncol(funcma) + idb[i]])
      }
      
      pairid = c(ncol(funcma) + idb,ida)
      names(pairid) = rownames(funump)[c(ida,ncol(funcma) + idb)]
      funump$xend = NA
      funump$yend = NA
      funump[names(pairid),]$xend <- funump$UMAP_1[pairid]  #set end point of segment from each point
      funump[names(pairid),]$yend <- funump$UMAP_2[pairid]
      # funump$highgroup = highgroup
      funump$clonegroup = funseur$clonecluster
      funump$clustergroup = Idents(funseur)
      
      funump$groupdist = NA
      funump$groupdist[c(ida,ncol(funcma) + idb)] = distvec 
      funump$highgroup = "background"
      
      funump = funump[order(funump$groupdist),]
      # for (i in unique(levels(funump$clustergroup))) {
      #   funump[which(funump$clustergroup == i)[1:10],"highgroup"] = as.character(i)
      # }
      
      funump$highgroup = "background"
      funump[which(!is.na(funump$clonegroup)),][1:200,"highgroup"] = funump[which(!is.na(funump$clonegroup)),][1:200,"clustergroup"]
      table(funump$highgroup)
      
      p5.2 = ggplot() + 
        geom_point(data = funump[which(funump$highgroup == "background"),],
                   aes(x = UMAP_1, y = UMAP_2),size = 1,color = "grey") +
        geom_point(data = funump[which(funump$highgroup != "background"),],
                   aes(x = UMAP_1, y = UMAP_2,color = as.character(clustergroup)),
                   size = 1.5,shape = 1) +
        geom_segment(data = funump[which(funump$highgroup != "background"),],
                     aes(x = UMAP_1, y = UMAP_2, xend = xend, yend = yend,color = as.character(clustergroup)),
                     size = 0.5,alpha = 0.5) + theme_bw() + 
        labs(color='UMAPcluster') 
      p5.2
      
      p5.2.1 = ggplot() + 
        geom_point(data = funump[which(funump$highgroup == "background"),],
                   aes(x = UMAP_1, y = UMAP_2),size = 1,color = "grey") +
        geom_point(data = funump[which(funump$highgroup != "background"),],
                   aes(x = UMAP_1, y = UMAP_2,color = as.character(clonegroup)),
                   size = 1.5,shape = 1) +
        geom_segment(data = funump[which(funump$highgroup != "background"),],
                     aes(x = UMAP_1, y = UMAP_2, xend = xend, yend = yend,color = as.character(clonegroup)),
                     size = 0.5,alpha = 0.5) + theme_bw() + 
        labs(color='clonecluster') 
      p5.2.1
      
      
      #dist densityplot
      # distvec
      distrand = c(distmat)[sample(1:length(distmat),length(distvec))]
      distdf = data.frame("dist" = c(distvec,distrand),
                          "group" = rep(c("sister_clone","random_clone"),c(length(distvec),length(distrand))))
      
      p5.3 = ggplot(distdf, aes(x = dist, col = group)) +
        geom_density() +
        scale_color_ipsum() +
        theme_bw(base_family = "sans") +
        # scale_x_continuous(limits = c(0,max(distdf$dist))) +
        ylab("proportion of clone pairs") + xlab("UMAP distance between clone pairs")
      p5.3
      
      #pearson
      nearest <- apply(distmat[1:ncol(funcma),ncol(funcma) + 1:ncol(funcmb)], 1, which.min)   #find index of nearest point to each point
      if(length(VariableFeatures(cellseur)) > 0){
        funcor = cor(funcma[VariableFeatures(funseur),],funcmb[VariableFeatures(funseur),],method = "pearson")
      }else{
        funcor = cor(funcma,funcmb,method = "pearson")
      }
      
      corsist = NULL
      cornear = NULL
      corsistcl = NULL
      for (i in 1:length(ida)) {
        cornear = c(cornear,funcor[i,nearest[i]])
        corsist = c(corsist,funcor[ida[i],idb[i]])
        corsistcl = rbind(corsistcl,data.frame("pearson" = funcor[ida[i],idb[i]],
                                               "group" = "sister",
                                               "subgroup" = funump[rownames(funcor)[ida[i]],"clonegroup"]))
      }
      # table(corsistcl$subgroup)
      corsistcl = corsistcl[which(!is.na(corsistcl$subgroup)),]
      corrand = c(as.matrix(funcor))[sample(1:length(c(as.matrix(funcor))),length(corsist))]
      cordf = data.frame("pearson" = c(corsist,corrand,cornear),
                         "group" = rep(c("sister","random","nearist"),c(length(corsist),length(corsist),length(corsist))),
                         "subgroup" = rep(c("sister total","random","nearist"),c(length(corsist),length(corsist),length(corsist))))
      corsistcl$subgroup = paste0("sister ",corsistcl$subgroup)
      cordf = rbind(cordf,corsistcl)
      # cordf$group = factor(cordf$group,levels = c("random","nearist","sister"))
      cordf$subgroup = factor(cordf$subgroup,levels = c("random","nearist","sister total",paste0("sister ",1:9),"sister NA"))
      
      res.aov <- aov(pearson ~ subgroup, data = cordf[which(cordf$group == "sister"),])
      res.aovp = summary(res.aov)
      resp = res.aovp[[1]]$`Pr(>F)`[1]
      
      p5.4 = ggplot(cordf, aes(x = subgroup,y = pearson, color = group)) +
        geom_boxplot() +
        stat_compare_means(comparisons = list(c("sister total","random"), c("sister total","nearist")),
                           label = "p.signif", method = "t.test") +
        annotate("text", x="sister 5", y=1, label= paste0("sister anova pvalue = ",resp)) + 
        scale_color_ipsum() +
        theme_bw(base_family = "sans") +
        ylab("pearson") + xlab("")
      p5.4
      
      
      #sister clone different cluster
      # cellseur <- FindClusters(cellseur, resolution = 0.1)
      # cellseur <- RunUMAP(cellseur, dims = 1:10)
      # DimPlot(cellseur)
      funump$clustergroup = Idents(funseur)[rownames(funump)]
      head(funump)
      #cluster distance
      funumppair = funump[1:(length(cmtag)*2),]
      funumppair$clustergroup = as.character(funumppair$clustergroup)
      funumprdm = funump[sample(1:nrow(funump),length(cmtag)*2),]
      funumprdm$clustergroup = as.character(funumprdm$clustergroup)
      pairstat = NULL
      pt = 1
      while(pt <= length(cmtag)*2){
        print(pt)
        if(funumppair[pt,"clustergroup"] == funumppair[pt+1,"clustergroup"]){
          pairline = 0
        }
        else{
          line = funumppair[pt,]
          linedist = distmat[rownames(line),]
          linedist = linedist[order(linedist)]
          linecluster = as.character(funump[names(linedist),"clustergroup"])
          clusterdist = unique(linecluster)
          pairline = which(clusterdist == funumppair[pt+1,"clustergroup"]) - 1
        }
        
        if(funumprdm[pt,"clustergroup"] == funumprdm[pt+1,"clustergroup"]){
          rdmline = 0
        }
        else{
          line = funumprdm[pt,]
          linedist = distmat[rownames(line),]
          linedist = linedist[order(linedist)]
          linecluster = as.character(funump[names(linedist),"clustergroup"])
          clusterdist = unique(linecluster)
          rdmline = which(clusterdist == funumprdm[pt+1,"clustergroup"]) - 1
        }
        pairstat = rbind(pairstat,data.frame("pair" = pairline,"random" = rdmline))
        pt = pt + 2
      }
      pairstat
      
      pairstatl = melt(pairstat)
      p5.5 = ggplot(pairstatl,aes(value,col = variable)) + 
        stat_ecdf(geom = "point") +
        # scale_color_manual(values = c("#D9534F","#FFAD60")) + 
        scale_color_manual(values = c("blue","#FFAD60")) + 
        labs(color='group') +
        xlab("Cluster distance") + ylab("ecdf of paired clone") +
        theme_pubr()
      p5.5
      
      
      #joint probability
      clusterid = Idents(funseur)
      clusterid = as.data.frame(clusterid)
      clusterid$group = substr(rownames(clusterid),0,1)
      clusterid$edit = substr(rownames(clusterid),3,nchar(rownames(clusterid)))
      clusterid$clusterid = as.numeric(as.character(clusterid$clusterid))
      
      clsispmx = matrix(0,length(levels(Idents(funseur))),length(levels(Idents(funseur))))
      rownames(clsispmx) = colnames(clsispmx) = as.numeric(levels(Idents(funseur)))
      for (i in as.numeric(levels(Idents(funseur)))) {
        print(i)
        for (j in as.numeric(levels(Idents(funseur)))) {
          if(i == j){
            cluster1 = clusterid[which(clusterid$clusterid == i & clusterid$group == "a"),]
            cluster2 = clusterid[which(clusterid$clusterid == i & clusterid$group == "b"),]
            tedit = nrow(cluster1) + nrow(cluster2)
            pedit = length(intersect(cluster1$edit,cluster2$edit)) * 2
            clsispmx[i+1,j+1] = pedit/tedit
          }else{
            cluster1 = clusterid[which(clusterid$clusterid == i),]
            cluster2 = clusterid[which(clusterid$clusterid == j),]
            tedit = nrow(cluster1) + nrow(cluster2)
            pedit = length(intersect(cluster1$edit,cluster2$edit)) * 2
            clsispmx[i+1,j+1] = pedit/tedit
          }
          
        }
      }
      diag(clsispmx)
      p5.6 = Heatmap(clsispmx,
                     row_order = rownames(clsispmx),
                     column_order = colnames(clsispmx),
                     col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                     cluster_rows = F,cluster_columns = F)
      p5.6
      p5 = list("umap_group" = p5.1,"sisterpair_umap_cluster" = p5.2,
                "sisterpair_clone_cluster" = p5.2.1,
                "umap_distance_density" = p5.3,
                "correlation_cmp" = p5.4,
                "paired_clone_ecdf" = p5.5,
                "paired_clone_colocation_heatmap" = p5.6)
      sisdatals = list("funseur"= funseur,"funump" = funump, "distdf" = distdf,
                    "cordf" = cordf, "pairstatl" = pairstatl,"clsispmx" = clsispmx)
      resultls = list("data" = sisdatals,"fg" = p5)
      
      return(resultls)
      
    }
    SisterUmapAnaPlot = function(cellcmpls){
      funseur = cellcmpls$data$funseur
      funump = cellcmpls$data$funump
      distdf = cellcmpls$data$distdf
      cordf = cellcmpls$data$cordf
      pairstatl = cellcmpls$data$pairstatl
      clsispmx = cellcmpls$data$clsispmx
      
      p5.1 = DimPlot(funseur) | DimPlot(funseur,split.by = "group",group.by = "clonecluster")
      p5.2 = ggplot() + 
        geom_point(data = funump[which(funump$highgroup == "background"),],
                   aes(x = UMAP_1, y = UMAP_2),size = 1,color = "grey") +
        geom_point(data = funump[which(funump$highgroup != "background"),],
                   aes(x = UMAP_1, y = UMAP_2,color = as.character(clustergroup)),
                   size = 1.5,shape = 1) +
        geom_segment(data = funump[which(funump$highgroup != "background"),],
                     aes(x = UMAP_1, y = UMAP_2, xend = xend, yend = yend,color = as.character(clustergroup)),
                     size = 0.5,alpha = 0.5) + theme_bw() + 
        labs(color='UMAPcluster') 
      p5.2
      
      p5.2.1 = ggplot() + 
        geom_point(data = funump[which(funump$highgroup == "background"),],
                   aes(x = UMAP_1, y = UMAP_2),size = 1,color = "grey") +
        geom_point(data = funump[which(funump$highgroup != "background"),],
                   aes(x = UMAP_1, y = UMAP_2,color = as.character(clonegroup)),
                   size = 1.5,shape = 1) +
        geom_segment(data = funump[which(funump$highgroup != "background"),],
                     aes(x = UMAP_1, y = UMAP_2, xend = xend, yend = yend,color = as.character(clonegroup)),
                     size = 0.5,alpha = 0.5) + theme_bw() + 
        labs(color='clonecluster') 
      p5.2.1
      p5.3 = ggplot(distdf, aes(x = dist, col = group)) +
        geom_density() +
        scale_color_ipsum() +
        theme_bw(base_family = "sans") +
        ylab("proportion of clone pairs") + xlab("UMAP distance between clone pairs")
      p5.3
      
      res.aov <- aov(pearson ~ subgroup, data = cordf[which(cordf$group == "sister"),])
      res.aovp = summary(res.aov)
      resp = res.aovp[[1]]$`Pr(>F)`[1]
      p5.4 = ggplot(cordf, aes(x = subgroup,y = pearson, color = group)) +
        geom_boxplot() +
        stat_compare_means(comparisons = list(c("sister total","random"), c("sister total","nearist")),
                           label = "p.signif", method = "t.test") +
        annotate("text", x="sister 5", y=1, label= paste0("sister anova pvalue = ",resp)) + 
        scale_color_ipsum() +
        theme_bw(base_family = "sans") +
        ylab("pearson") + xlab("")
      p5.4
      p5.5 = ggplot(pairstatl,aes(value,col = variable)) + 
        stat_ecdf(geom = "point") +
        scale_color_manual(values = c("#D9534F","#FFAD60")) + 
        labs(color='group') +
        xlab("Cluster distance") + ylab("ecdf of paired clone") +
        theme_bw()
      p5.5
      p5.6 = Heatmap(clsispmx,
                     row_order = rownames(clsispmx),
                     column_order = colnames(clsispmx),
                     col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                     cluster_rows = F,cluster_columns = F)
      p5.6
      fg = list(p5.1,p5.2,p5.2.1,p5.3,p5.4,p5.5,p5.6)
      return(fg)
      
    }
    ClonetoCell = function(e11cmar,e11cmbr,cmtagr){
      smpnamea = colnames(dcast(e11cmar,Cell.type~pattern))[-1]
      smpnameb = colnames(dcast(e11cmbr,Cell.type~pattern))[-1]
      linebca = as.character(e11cmar[which(e11cmar$pattern == smpnamea[1]),"Cell.BC"])
      lclrnaa = e11rna[,colnames(e11rna) %in% linebca]
      if(!is.null(ncol(lclrnaa))){
        lclrnaas = as.data.frame(rowMeans(lclrnaa))
      }else{
        lclrnaas = as.data.frame(lclrnaa)
      }
      e11clrnaa = lclrnaas
      colnames(e11clrnaa)[1] = smpnamea[1]
      
      for (i in 2:length(smpnamea)) {
        print(i)
        linebca = as.character(e11cmar[which(e11cmar$pattern == smpnamea[i]),"Cell.BC"])
        lclrnaa = e11rna[,colnames(e11rna) %in% linebca]
        if(!is.null(ncol(lclrnaa))){
          lclrnaas = as.data.frame(rowMeans(lclrnaa))
        }else{
          lclrnaas = as.data.frame(lclrnaa)
        }
        e11clrnaa = cbind(e11clrnaa,lclrnaas)
        colnames(e11clrnaa)[i] = smpnamea[i]
      }
      
      linebcb = as.character(e11cmbr[which(e11cmbr$pattern == smpnameb[1]),"Cell.BC"])
      lclrnab = e11rna[,colnames(e11rna) %in% linebcb]
      if(!is.null(ncol(lclrnab))){
        lclrnabs = as.data.frame(rowMeans(lclrnab))
      }else{
        lclrnabs = as.data.frame(lclrnab)
      }
      e11clrnab = lclrnabs
      colnames(e11clrnab)[1] = smpnameb[1]
      for (i in 2:length(smpnameb)) {
        print(i)
        linebcb = as.character(e11cmbr[which(e11cmbr$pattern == smpnameb[i]),"Cell.BC"])
        lclrnab = e11rna[,colnames(e11rna) %in% linebcb]
        if(!is.null(ncol(lclrnab))){
          lclrnabs = as.data.frame(rowMeans(lclrnab))
        }else{
          lclrnabs = as.data.frame(lclrnab)
        }
        e11clrnab = cbind(e11clrnab,lclrnabs)
        colnames(e11clrnab)[i] = smpnameb[i]
      }
      
      clrida = match(cmtagr,colnames(e11clrnaa))
      clridb = match(cmtagr,colnames(e11clrnab))
      colnames(e11clrnaa) = paste0("a-",colnames(e11clrnaa))
      colnames(e11clrnab) = paste0("b-",colnames(e11clrnab))
      transeur <- CreateSeuratObject(counts = cbind(e11clrnaa,e11clrnab), min.cells = 0)
      group = rep(c("A","B"),c(ncol(e11clrnaa),ncol(e11clrnab)))
      names(group) = colnames(cbind(e11clrnaa,e11clrnab))
      transeur$group = group
      transeur$clonecluster = clonecluster[names(clonecluster) %in% colnames(transeur)]
      transeur <- ScaleData(transeur)
      transeur <- FindVariableFeatures(transeur, selection.method = "vst", nfeatures = 2000)
      transeur <- RunPCA(transeur, features = VariableFeatures(transeur))
      transeur <- FindNeighbors(transeur, dims = 1:5)
      return(list("seur" = transeur,"cma" = e11clrnaa,
                  "cmb" = e11clrnab, "ida" = clrida, "idb" = clridb))
    }
    tmpfg = SisterUmapAnaPlot(cellcmpls)
    
    e11cmarp = e11cmar[which(e11cmar$Cell.type %in% progcell),]
    e11cmbrp = e11cmbr[which(e11cmbr$Cell.type %in% progcell),]
    cmtagrp = intersect(unique(e11cmarp$pattern), unique(e11cmbrp$pattern))
    e11cmaro = e11cmar[which(!e11cmar$Cell.type %in% progcell),]
    e11cmbro = e11cmbr[which(!e11cmbr$Cell.type %in% progcell),]
    cmtagro = intersect(unique(e11cmaro$pattern), unique(e11cmbro$pattern))
    
    transeurls = ClonetoCell(e11cmar,e11cmbr,cmtagc)
    transeurpls = ClonetoCell(e11cmarp,e11cmbrp,cmtagrp)
    transeurols = ClonetoCell(e11cmaro,e11cmbro,cmtagro)
    
    transeurls$seur <- FindClusters(transeurls$seur, resolution = 0.8)
    transeurls$seur <- RunUMAP(transeurls$seur, dims = 1:10)
    
    transeurpls$seur <- FindClusters(transeurpls$seur, resolution = 0.8)
    transeurpls$seur <- RunUMAP(transeurpls$seur, dims = 1:10)
    
    transeurols$seur <- FindClusters(transeurols$seur, resolution = 0.8)
    transeurols$seur <- RunUMAP(transeurols$seur, dims = 1:10)
    
    # cellcmpls = SisterUmapAna(cellseur, e11cma, e11cmb,
    #                           cmida,cmidb,cmtagc)
    transcmpls = SisterUmapAna(transeurls$seur,transeurls$cma,transeurls$cmb,
                               transeurls$ida,transeurls$idb,cmtagc)
    transcmppls = SisterUmapAna(transeurpls$seur,transeurpls$cma,transeurpls$cmb,
                                transeurpls$ida,transeurpls$idb,cmtagrp)
    transcmpols = SisterUmapAna(transeurols$seur,transeurols$cma,transeurols$cmb,
                                transeurols$ida,transeurols$idb,cmtagro)
    
    transcmpls$data = transcmpls$data[-1]
    transcmppls$data = transcmppls$data[-1]
    transcmpols$data = transcmpols$data[-1]
    
    ggsave(transcmpls$fg[[1]],filename = "final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_umap.pdf",width = 10,height = 5)
    pdf("final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_cmp.pdf",width = 7,height = 5)
    invisible(lapply(transcmpls$fg[-1], print))
    dev.off()
    
    ggsave(transcmppls$fg[[1]],filename = "final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_umap_prog.pdf",width = 10,height = 5)
    pdf("final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_cmp_prog.pdf",width = 7,height = 5)
    invisible(lapply(transcmppls$fg[-1], print))
    dev.off()
    
    ggsave(transcmpols$fg[[1]],filename = "final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_umap_other.pdf",width = 10,height = 5)
    pdf("final_result/CREST_final_analysis/sister_cmp_E11/E11_sister_clone_trans_cmp_other.pdf",width = 7,height = 5)
    invisible(lapply(transcmpols$fg[-1], print))
    dev.off()
    
    # transcmpols$fg$umap_distance_density
    sistercmpls = list("prog" = transcmppls$data,"other" = transcmpols$data,"total" = transcmpls$data)
    # sistercmpls$prog$fg = sistercmpls$prog$fg[-1]
    # sistercmpls$other$fg = sistercmpls$other$fg[-1]
    # sistercmpls$total$fg = sistercmpls$total$fg[-1]
    
  }
  dir.create("final_result/CREST_final_analysis/sister_cmp_E11")
  qsave(sistercmpls,file = "final_result/CREST_final_analysis/sister_cmp_E11/sistercmpls.qs")
  datals$figure$sister_e11_cmp_fs9 = sistercmpls
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_10_edit.qs")
  
}

#sfig2. basic eidt stat bulk---------
{
  #E8.5-P0
  bulkls = qread("bulk_result/X.bulk_in_vivo/bulk_total_12_9.qs")
  bulkstatls = list()
  
  for (i in 1:length(bulkls$v1)) {
    bulkls$v1[[i]] = bulkls$v1[[i]][which(bulkls$v1[[i]]$reads_num > 3 & bulkls$v1[[i]]$reads_pro.main > 0.5),]
    bulkls$v2[[i]] = bulkls$v2[[i]][which(bulkls$v2[[i]]$reads_num > 3 & bulkls$v2[[i]]$reads_pro.main > 0.5),]
    
  }
  ggvenn(list("nkx1" = bulkls$v1$nkx1$main,"nkx2" = bulkls$v1$nkx2$main,"nkx3" = bulkls$v1$nkx3$main))
  
  samplels = c(list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/previous/",full.names = T),
               list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/",full.names = T)[-13])
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
  names(v1lsx) = names(v2lsx) = c(list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/previous/"),
                                  list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/")[-13])
  bulkls = list("v1" = v1lsx,"v2" = v2lsx)
  qsave(bulkls,"bulk_result/X.bulk_in_vivo/bulk_total_1_5.qs")
  
  v1lsx = bulkls$v1[c(12,10,1:9,11)]
  v2lsx = bulkls$v2[c(12,10,1:9,11)]
  colrank = names(v1lsx) = names(v2lsx) = c("E8.5","E9.5","E10-1","E10-2","E12","E15-1a","E15-1b",
                                  "E15-2a","E15-2b","E15-3a","E15-3b","P0")
  
  ScartoIndel = function(tl){
    indel = data.frame(id=NA,start=NA,end=NA,width=NA,type=NA)
    indel = indel[-1,]
    #build indel data frame
    for(i in c(1:length(tl[-1]))){
      x = tl[i]
      x1 = strsplit(x,split = "-")[[1]][1]
      x_str = strsplit(x1,split = "_")
      del = NULL
      ins = NULL
      
      x_pre = unique(x_str[[1]])
      if(length(x_pre)>0){
        x_pre_del = x_pre[grep("D",x_pre)]
        if(length(x_pre_del) > 0){
          x_pre_del = str_extract_all(x_pre_del, "[0-9]+")
          for (j in 1:length(x_pre_del)) {
            line = as.numeric(x_pre_del[[j]])
            del = rbind(del, data.frame(id=x,start=line[length(line)],
                                        end=line[length(line)]+line[length(line)-1],width=line[length(line)-1],
                                        type="deletion"))
          }
          
        }
        x_pre_ins = x_pre[grep("I",x_pre)]
        if(length(x_pre_ins) > 0){
          x_pre_ins = str_extract_all(x_pre_ins, "[0-9]+")
          for (j in 1:length(x_pre_ins)) {
            line = as.numeric(x_pre_ins[[j]])
            ins = rbind(ins, data.frame(id=x,start=line[length(line)],
                                        end=line[length(line)]+line[length(line)-1],width=line[length(line)-1],
                                        type="insertion"))
          }
        }
        
      }
      
      indel = rbind(indel,del,ins)
    }
    return(indel)
  }
  ScarPlot = function(v1lsx,num,colrank){
    
    scardft = NULL
    index = NULL
    
    for (i in 1:length(v1lsx)) {
      
      scardf = as.data.frame(table(v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]$umim))
      scardf = scardf[which(!scardf$Var1 %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      scardf = scardf[order(-scardf$Freq),]
      scardf$Freq = scardf$Freq/sum(scardf$Freq)
      
      scardf = scardf[1:num,]
      colnames(scardf)[2] = names(v1lsx)[i]
      if(i == 1){
        scardft = scardf
      }else{
        scardft = merge(scardft, scardf, by = "Var1", all = T)
        
      }
      
    }
    
    scardft[is.na(scardft)] = 0
    scardft = scardft[do.call("order", c(scardft[-1], list(decreasing=TRUE))),]
    indel = ScartoIndel(as.character(scardft$Var1))
    index = data.frame(id = scardft$Var1,order = nrow(scardft):1)
    indel = merge(indel,index,by = "id")
    scardfti = melt(scardft)
    colnames(scardfti) = c("id","sample","proportion")
    scardfti = merge(scardfti,index,by = "id")
    
    scardfti$sample = factor(scardfti$sample,levels = colrank)
    p3.0 = ggplot() + geom_segment(data = indel,aes(x = start, xend = end, y = order, yend = order,color=type),
                                   size = 0.7) +
      theme_void()  + theme(legend.position = "left") + ylab("") + 
      ggplot() + geom_tile(data = scardfti,aes(x = sample,y = order,fill = proportion,height = 1,width = 1),size = 2) + 
      ylab("") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(),axis.title.x = element_blank()) + 
      # scale_fill_viridis(option = "A")
      scale_fill_gradient2(midpoint = 0.001,low = "blue",high = "red")
    return(list("scardf" = scardfti,"indel" = indel,"figure" = p3.0))
  }
  
  library(stringr)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  scardf1 = ScarPlot(v1lsx,20,colrank)
  scardf2 = ScarPlot(v2lsx,20,colrank)
  bulkstatls$array_heatmap = list("v1" = scardf1,"v2" = scardf2)
  
  #1. 统计每种array在多少个样本中出现
  StatArrayOcc = function(v1lsx){
    arraydft = NULL
    for (i in 2:length(v1lsx)) {
      arraydf = as.data.frame(table(v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]$main))
      arraydf = arraydf[which(!arraydf$Var1 %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                   "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
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
    arraydft$E15_1 = (arraydft$`E15-1a`+ arraydft$`E15-1b`)/2
    arraydft$E15_2 = (arraydft$`E15-2a`+ arraydft$`E15-2b`)/2
    arraydft$E15_3 = (arraydft$`E15-3a`+ arraydft$`E15-3b`)/2
    arraydft = arraydft[!colnames(arraydft) %in% c("E15-1a","E15-1b","E15-2a","E15-2b",
                                                   "E15-3a","E15-3b")]
    
    ayoc = as.data.frame(table(rowSums(arraydft[-1]>0)))
    ayoc$norm = ayoc$Freq/sum(ayoc$Freq)
    return(ayoc)
  }
  
  ayocv1 = StatArrayOcc(v1lsx)
  ayocv2 = StatArrayOcc(v2lsx)
  ayocv1$group = "v1"
  ayocv2$group = "v2"
  ayocv = rbind(ayocv1,ayocv2)
  p1 = ggplot(ayocv,aes(x = Var1, y = norm, fill = group)) + 
    geom_bar(stat = "identity",position = "dodge",width = 0.7) +
    scale_fill_manual(values=c('red','blue')) +
    theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) + theme(panel.grid = element_blank())
  print(p1)
  bulkstatls$array_occurrence = list("arrayocr" = ayocv,"figure" = p1)
  
  #2. bulk indel occurrence histogram
  {
    v1lsx[[6]] = rbind(v1lsx[[6]],v1lsx[[7]])
    v1lsx[[8]] = rbind(v1lsx[[8]],v1lsx[[9]])
    v1lsx[[10]] = rbind(v1lsx[[10]],v1lsx[[11]])
    names(v1lsx)[c(6,8,10)] = c("E15-1","E15-2","E15-3")
    v1lsx = v1lsx[-c(7,9,11)]
    
    v2lsx[[6]] = rbind(v2lsx[[6]],v2lsx[[7]])
    v2lsx[[8]] = rbind(v2lsx[[8]],v2lsx[[9]])
    v2lsx[[10]] = rbind(v2lsx[[10]],v2lsx[[11]])
    names(v2lsx)[c(6,8,10)] = c("E15-1","E15-2","E15-3")
    v2lsx = v2lsx[-c(7,9,11)]
    
    id1 = NULL
    id2 = NULL
    for (i in 2:length(v1lsx)) {
      line1 = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
      line2 = v2lsx[[i]][which(v2lsx[[i]]$reads_pro.main >= 0.5 & v2lsx[[i]]$reads_num >= 3),]
      line.id1 = unique(unlist(strsplit(unique(line1$main),split = "_|&")))
      line.id2 = unique(unlist(strsplit(unique(line2$main),split = "_|&")))
      id1 = c(id1,line.id1)
      id2 = c(id2,line.id2)
    }
    
    idocv1 = as.data.frame(table(id1))
    idocv2 = as.data.frame(table(id2))
    idocv1 = as.data.frame(table(idocv1$Freq))
    idocv2 = as.data.frame(table(idocv2$Freq))
    idocv1$norm = idocv1$Freq / sum(idocv1$Freq)
    idocv2$norm = idocv2$Freq / sum(idocv2$Freq)
    
    idocv1$group = "V1"
    idocv2$group = "V2"
    idocv = rbind(idocv1,idocv2)
    
    p1.1 = ggplot(idocv,aes(x = Var1, y = norm, fill = group)) + 
      geom_bar(stat = "identity",position = "dodge",width = 0.7) +
      scale_fill_manual(values=c('red','blue')) +
      theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
      scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1)) + theme(panel.grid = element_blank())
    print(p1.1)
    bulkstatls$indel_occurrence = list("indelocr" = idocv, "figure" = p1.1)
    
  }
  
  #3. 
  FilterBlacklist = function(v1lsx,bl1){
    v1lsx.f = list()
    for (i in 1:length(v1lsx)){
      df = v1lsx[[i]]
      df = df[which(df$reads_pro.main >= 0.5 & df$reads_num >= 3),]
      newdf = NULL
      for (j in 1:nrow(df)) {
        x = df$umim[j]
        flag = 0
        flag = length(bl1[,1][which(bl1[,1] == x)])
        if(flag == 0){
          x_str = strsplit(x,split = "_")
          x_pre = x_str[[1]]
          x_pre = unique(x_pre[!(x_pre %in% "NONE")])
          flag = length(intersect(x_pre,bl1[,1]))
        }
        
        if(flag==0){
          newdf = rbind(newdf,df[j,])
        }
        
      }
      v1lsx.f[[i]] = newdf
    }
    names(v1lsx.f) = names(v1lsx)
    return(v1lsx.f)
  }
  v1lsx.f = FilterBlacklist(v1lsx,datals$blacklist$v1_array)
  v2lsx.f = FilterBlacklist(v2lsx,datals$blacklist$v2_array)
  
  #编辑size统计
  EditSizePlot = function(v1lsx.f)
  {
    indelt = NULL
    for (i in 2:length(v1lsx.f)) {
      indel = ScartoIndel(as.character(v1lsx.f[[i]]$umim))
      indelt = rbind(indelt,indel)
    }
    indelpd = as.data.frame(table(indelt[,c(5,4)]))
    indelpd$norm = indelpd$Freq/nrow(indelt)
    indelpd[which(indelpd$type == "insertion"),"norm"] = -indelpd[which(indelpd$type == "insertion"),"norm"]
    indelpd$width = as.numeric(as.character(indelpd$width))
    
    p = ggplot(indelpd,aes(x = width, y = norm,fill = type)) + 
      geom_bar(stat="identity")+theme_bw() +
      scale_fill_manual(values=c('lightpink1','lightblue2')) +
      xlab("Indel size (bp)")+ ylab("Occurrences") +
      scale_x_continuous(limits = c(0,max(indelpd$width)+1),breaks=seq(0,max(indelpd$width)+1,20))
    print(p)
    return(p)
  }
  v1edp = EditSizePlot(v1lsx.f)
  v2edp = EditSizePlot(v2lsx.f)
  bulkstatls$Indel_editsize_figure = list("v1" = v1edp,"v2" = v2edp)
  # ggsave(v1edp,filename = "final_result/xBulk/edit_size_v1_2_17_qc.pdf",width = 10,height = 6)
  # ggsave(v2edp,filename = "final_result/xBulk/edit_size_v2_2_17_qc.pdf",width = 10,height = 6)
  
  #edit array proportions
  names(v1lsx)
  EditArrayStat = function(v1lsx){
    propft = NULL
    for (i in 1:length(v1lsx)) {
      arraydf = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
      arraydff = arraydf[which(!arraydf$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                    "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      propl = nrow(arraydff) / nrow(arraydf)
      propft = rbind(propft,data.frame("sample" = names(v1lsx)[i],"edit_prop" = propl))
    }
    
    #intact target sites
    
    tarsite = NULL
    for (i in 1:length(v1lsx)) {
      # arraydf = v1lsx[[i]]
      arraydf = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
      # arraydf = arraydf[which(!arraydf$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
      #                                               "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      linetar = NULL
      for (j in 1:nrow(arraydf)) {
        line = arraydf[j,"main"]
        x1 = strsplit(line,split = "-")[[1]][1]
        x_str = strsplit(x1,split = "_")[[1]]
        linetar = c(linetar,length(x_str[which(x_str == "NONE")]))
      }
      tarsite = c(tarsite,mean(linetar))
    }
    propft$intact_site = tarsite
    return(propft)
  }
  v1propft = EditArrayStat(v1lsx)
  v2propft = EditArrayStat(v2lsx)
  bulkstatls$Edit_array_proportion = list("v1" = v1propft,"v2" = v2propft)
  dir.create("final_result/CREST_final_analysis/bulk_analysis")
  write.table(v1propft,file = "final_result/xBulk/eidt_stat_v1_2_17_qc.txt",row.names = F,
              quote = F,sep = "\t")
  write.table(v2propft,file = "final_result/xBulk/eidt_stat_v2_2_17_qc.txt",row.names = F,
              quote = F,sep = "\t")
  
  #big cut array stat
  CutThroughProp = function(v1lsx){
    ctdf = NULL
    maxn = length(unlist(strsplit(split = "_",v1lsx[[1]]$umim[1])))
    for (i in 1:length(v1lsx)) {
      v1i = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
      v1i = v1i[which(!v1i$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                    "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      indi = strsplit(split = "_|&",v1i$umim)
      indil = unlist(lapply(indi, function(x)return(max(table(x)))))
      ctrate = length(indil[which(indil == maxn)])/length(indil)
      ctdf = rbind(ctdf,data.frame("sample" = names(v1lsx)[i],"cutthrough_rate" = ctrate))
    }
    return(ctdf)
  }
  v1lsx$nkx1 = bulkls$v1$nkx1;v1lsx$nkx2 = bulkls$v1$nkx2;v1lsx$nkx3 = bulkls$v1$nkx3;
  v2lsx$nkx1 = bulkls$v2$nkx1;v2lsx$nkx2 = bulkls$v2$nkx2;v2lsx$nkx3 = bulkls$v2$nkx3;
  ctdfv1 = CutThroughProp(v1lsx)
  ctdfv2 = CutThroughProp(v2lsx)
  ctdfv1$group = "V1"
  ctdfv2$group = "V2"
  ctdf = rbind(ctdfv1[1:9,],ctdfv2[1:9,])
  ctdf$sample = unlist(lapply(strsplit(ctdf$sample,split = "-"),"[[",1))
  ctdf = ctdf %>% group_by(group,sample) %>% summarise(meanct = mean(cutthrough_rate),sd = sd(cutthrough_rate))
  ctdf[is.na(ctdf$sd),"sd"] = 0
  ctdf$sampleid = as.numeric(substr(ctdf$sample,2,nchar(ctdf$sample)))
  ctdf[which(ctdf$sampleid == 0),"sampleid"] = 20.5
  ctdf[which(ctdf$sampleid == 15),"sampleid"] = 15.5
  ctdf[which(ctdf$sampleid == 10),"sampleid"] = 10.5
  ctdf[which(ctdf$sampleid == 12),"sampleid"] = 12.5
  
  pd = position_dodge(0.1)
  yline1 = mean(ctdfv1[10:12,2])
  yline2 = mean(ctdfv2[10:12,2])
  # ctdf$sample = factor(ctdf$sample,levels = unique(ctdf$sample))
  p1.2 = ggplot(ctdf, aes(x = sampleid,y = meanct,color = group,group = group)) + 
    geom_errorbar(aes(ymin=meanct-sd, ymax=meanct+sd,color = group), width=.5, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, aes(shape=group)) + 
    scale_shape_manual(values = c(19,15))+
    geom_hline(yintercept = yline1,color = "red", linetype='dotted') +
    geom_hline(yintercept = yline2,color = "blue", linetype='dotted') +
    scale_x_continuous(breaks = c(8.5,9.5,10.5,12.5,15.5,20.5),labels = c('E8.5','E9.5','E10.5','E12.5','E15.5','P0'))+
    scale_color_manual(values = c("red","blue")) +
    ylab("Cut Through Array rate") + 
    theme_pubr() + xlab("")
    # theme(legend.justification=c(1,0),panel.grid = element_blank(),
    #       legend.position=c(1,0))
  bulkstatls$array_cutthroughRate = list("v1ctdf" = ctdfv1,"v2ctdf" = ctdfv2
                                         ,"ctdftotal" = ctdf,"figure" = p1.2)
  qsave(bulkstatls,file = "final_result/CREST_final_analysis/bulk_analysis/bulkstatls.qs")
  datals$figure$bulkstatls_fs2 = bulkstatls
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_10_edit.qs")
  
  
  # stat added bulk sample
  samplels = c(list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/previous/",full.names = T),
               list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/",full.names = T)[-13])
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
  names(v1lsx) = names(v2lsx) = c(list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/previous/"),
                                  list.files("bulk_result/X.bulk_in_vivo/bulk_20221125/")[-13])
  for (i in 1:length(v1lsx)) {
    v1lsx[[i]]$main = unlist(lapply(strsplit(v1lsx[[i]]$main,"_"),function(x){paste0(x[-4],collapse = "_")}))
  }
  v1propft = EditArrayStat(v1lsx)
  v2propft = EditArrayStat(v2lsx)
  tmp = v1lsx$e7_pos
  write.table(v1propft,file = "final_result/CREST_final_analysis/rebbutal_analysis/eidt_stat_v1_12_19_qc.txt",row.names = F,
              quote = F,sep = "\t")
  write.table(v2propft,file = "final_result/CREST_final_analysis/rebbutal_analysis/eidt_stat_v2_12_19_qc.txt",row.names = F,
              quote = F,sep = "\t")
  
  
  #stat V1 cut site and V2 cut site
  v1lsxi = bulkls$v1[c(1:10)]
  v2lsxi = bulkls$v2[c(1:10)]
  
  for (i in 1:length(v1lsxi)) {
    v1lsxi[[i]]$umim = unlist(lapply(strsplit(v1lsxi[[i]]$umim,"_"),function(x){paste0(x[-4],collapse = "_")}))
  }
  CutNumStat = function(v1lsxi,group){
    ctdf = NULL
    idldf = NULL
    arldf = NULL
    idstdf = NULL
    maxn = length(unlist(strsplit(split = "_",v1lsxi[[1]]$umim[1])))
    for (i in 1:length(v1lsxi)) {
      v1i = v1lsxi[[i]][which(v1lsxi[[i]]$reads_pro.main >= 0.5 & v1lsxi[[i]]$reads_num >= 3),]
      v1i = v1i[which(!v1i$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                       "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      
      arrayi = v1i$umim
      indi = strsplit(split = "_",arrayi)
      indil = unlist(lapply(indi, function(x)return(length(x[which(x != "NONE")]))))
      
      ctrate = indil
      ctdf = rbind(ctdf,data.frame("sample" = names(v1lsxi)[i],"cutrate" = ctrate))
      
      indi2 = strsplit(split = "_|&",arrayi)
      indi2l = unlist(lapply(indi2 , function(x){
        x = unique(x)
        xd = x[grep("D",x)]
        xdn = as.numeric(unlist(lapply(strsplit(xd,"D"),"[[",1)))
        return(xdn)
      }))
      indi2d = unlist(lapply(indi2 , function(x){
        x = unique(x)
        xd = x[grep("D",x)]
        return(xd)
      }))
      idldf = rbind(idldf,data.frame("sample" = names(v1lsxi)[i],"indel" = indi2d,"indelen" = indi2l))
      
      arl = unlist(lapply(indi2 , function(x){
        x = unique(x)
        xd = x[grep("D",x)]
        xdn = as.numeric(unlist(lapply(strsplit(xd,"D"),"[[",1)))
        return(sum(xdn))
      }))
      arldf = rbind(arldf,data.frame("sample" = names(v1lsxi)[i],"array" = arrayi,"arraylen" = arl))
      
      
      indls = unlist(lapply(indi2, function(x)return(table(x))))
      idstdf = rbind(idstdf,data.frame("indel" = names(indls),
                                       "sample" = names(v1lsxi)[i],
                                       "sites" = indls))
      
      
    }
    ctdf$group = group
    arldf$group = group
    idldf$group = group
    idstdf$group = group

    return(list("ctn" = ctdf,"arl" = arldf, "idl" = idldf,"ids" = idstdf))
  }
  
  v1ctst = CutNumStat(v1lsxi,"V1")
  v2ctst = CutNumStat(v2lsxi,"V2")

  cratedf = rbind(v1ctst$ctn,v2ctst$ctn)
  # cratedf$cutrate = cratedf$cutrate/8
  pb.1 = ggboxplot(cratedf, x = "group", y = "cutrate", fill = "group",
                   width = 0.5) + 
    scale_fill_jama() +
    stat_compare_means(aes(group = group),
                       label = "p.signif",
                       method = "t.test") +
    theme(legend.position = "none") +
    xlab("") + ylab("Edited site number in each array") 
  pb.1
  
  v1ctst$idl$lenrate = v1ctst$idl$indelen/239
  v2ctst$idl$lenrate = v2ctst$idl$indelen/212
  idldf = rbind(v1ctst$idl,v2ctst$idl)
  pb.2 = ggboxplot(idldf, x = "group", y = "lenrate", fill = "group",
                   width = 0.5) + 
    scale_fill_jama() +
    stat_compare_means(aes(group = group),
                       label = "p.signif",
                       method = "t.test") +
    theme(legend.position = "none") +
    xlab("") + ylab("Ratio of InDel length to ref array") 
  pb.2
  v1ctst$arl$arraylenr = 233 - v1ctst$arl$arraylen
  v2ctst$arl$arraylenr = 212 - v2ctst$arl$arraylen
  v1ctst$arl$lenrate = 1 - v1ctst$arl$arraylen/233
  v2ctst$arl$lenrate = 1 - v2ctst$arl$arraylen/212
  
  arldf = rbind(v1ctst$arl,v2ctst$arl)
  pb.3 = ggboxplot(arldf, x = "group", y = "arraylenr", fill = "group",
                   width = 0.5) + 
    scale_fill_jama() +
    stat_compare_means(aes(group = group),
                       label = "p.signif",
                       method = "t.test") +
    theme(legend.position = "none") +
    xlab("") + ylab("Remaining barcode length") 
  pb.3
  
  
  idsdf = rbind(v1ctst$ids,v2ctst$ids)
  idsdf = idsdf[grep("D",idsdf$indel),]
  table(idsdf$indel)[order(-table(idsdf$indel))]
  idsdf = idsdf %>% group_by(sample,group,sites) %>% summarise(count = n())
  idsdf = idsdf %>% group_by(sample,group) %>% mutate(norm = count/sum(count))
  tmp = idsdf %>% group_by(group,sites) %>% summarise(mean = mean(norm),sd = sd(norm))

  pb.4 = ggbarplot(idsdf, x = "sites", fill = "group",
                   y = "norm",add = c("mean_sd"),
                   position = position_dodge()) +
    scale_fill_jama() +
    stat_compare_means(aes(group = group),
                       label = "p.signif",
                       method = "t.test") +
    scale_x_continuous(breaks = c(1:8)) +
    theme(legend.position = "bottom")+
    xlab("Indel cut sites number distribution") + ylab("") 
  pb.4
  
  pbt  = ggarrange(ggarrange(pb.1,pb.2,pb.3,nrow = 1),pb.4,nrow = 2)
  pbt
  write.csv(tmp,file = "final_result/CREST_final_analysis/bulk_analysis/bulk_v1_v2_edit_type_cutsite_density.csv",
            quote = F,row.names = F)
  ggsave(pbt,filename = "final_result/CREST_final_analysis/bulk_analysis/bulk_v1_v2_edit_type_cmp.pdf",
         width = 8,height = 7)
  
  
  
}

#sfig5. sc sample v1 v2 recover stat & E11 E15 clone size stat
{
  scbasicstls = qread("final_result/CREST_final_analysis/E11_E15_basic_stat/scbasicstls.qs")
  scbasicstls = list()
  # datals$cellls = cellls
  #统计E11
  cellls = datals$cellls
  #(1) recover stat
  {
    calid = c("E11_rep1","E11_rep2","E11_rep3","E15_rep1","E15_rep2","E15_rep3")
    #array CB in trans/trans CB number
    recBCt = NULL
    for (i in 1:length(names(datals$v1)[1:12])) {
      cali = names(datals$v1)[i]
      v1bc = datals$v1[[cali]]$resultraw$Cell.BC
      v2bc = datals$v2[[cali]]$resultraw$Cell.BC
      v1v2bc = union(v1bc,v2bc)
      v1v2bconly = intersect(v1bc,v2bc)
      transcb = data.frame("BC" = cellls[[cali]]$Cell.BC,"Cell.type" = cellls[[cali]]$Cell.type)
      v1rec = transcb[which(transcb$BC %in% v1bc),]
      v2rec = transcb[which(transcb$BC %in% v2bc),]
      v1v2rec = transcb[which(transcb$BC %in% v1v2bc),]
      unrec = transcb[which(!transcb$BC %in% v1v2bc),]
      v1v2reconly = transcb[which(transcb$BC %in% v1v2bconly),]
      
      v1rec = v1rec %>% group_by(Cell.type) %>% summarise(cell.counts = n())
      v2rec = v2rec %>% group_by(Cell.type) %>% summarise(cell.counts = n())
      v1v2rec = v1v2rec %>% group_by(Cell.type) %>% summarise(cell.counts = n())
      v1v2reconly = v1v2reconly %>% group_by(Cell.type) %>% summarise(cell.counts = n())
      unrec = unrec %>% group_by(Cell.type) %>% summarise(cell.counts = n())
      
      v1rec$recover = "V1_rec";v2rec$recover = "V2_rec";v1v2rec$recover = "V1/V2_rec";unrec$recover = "unrecovered";v1v2reconly$recover = "V1&V2_rec";
      recBC = rbind(v1rec,v2rec,v1v2rec,unrec,v1v2reconly)
      transcb = transcb %>% group_by(Cell.type) %>% summarise(total.counts = n())
      recBC = merge(recBC,transcb)
      recBC$cell.prop = recBC$cell.counts/recBC$total.counts
      recBC$sample = cali
      # reBC11 = data.frame("BC" = cellls[[cali]]$Cell.BC,
      #                     "recover" = "unrecover","Cell.type" = cellls[[cali]]$Cell.type)
      # reBC11[which(reBC11$BC %in% v1bc11),"recover"] = "V1"
      # reBC11[which(reBC11$BC %in% v2bc11),"recover"] = "V2"
      # reBC11[which(reBC11$BC %in% v1v2bc11),"recover"] = "V1_V2"
      # reBC11 = as.data.frame(table(reBC11[2:3]))
      # colnames(reBC11) = c("recover","Cell.type","cell.counts")
      recBCt = rbind(recBCt,recBC)
    }
    # reBCdf$recover = factor(as.character(reBCdf$recover),levels = c("V1","V1_V2","V2","unrecover"))
    # reBC$Cell.type = factor(reBC$Cell.type,levels = as.character(0:26))
    # p4.0 = ggplot(reBCdf,aes(x = Cell.type,y = cell.counts,fill = recover)) + 
    #   geom_bar(position="fill", stat="identity",width = 0.8) +
    #   facet_wrap(~sample,scales = "free",ncol = 3) + 
    #   theme_bw() + theme(axis.text.x = element_text(angle = 90))
    # p4.0
    head(recBCt)
    recBCt$group = unlist(lapply(strsplit(recBCt$sample,split = "_"), "[[", 1))
    
    
    recBCt$recover = factor(recBCt$recover,levels = c("unrecovered","V1_rec","V2_rec","V1&V2_rec",
                                                      "V1/V2_rec"))
    recBCt[which(recBCt$group == "pandaE11"),"group"] = "E11"
    p4.0.1 = ggbarplot(
      recBCt, x = "Cell.type", y = "cell.prop", facet.by = "group",
      scales = "free",
      add = c("mean_sd", "jitter"), 
      color = "recover", 
      position = position_dodge(0.8)
    ) + 
      # facet_wrap(~group,scales = "free",ncol = 2) + 
      scale_color_jama() + 
      scale_y_continuous(limits = c(0,1)) + 
      xlab("") + ylab("Recovered cells proportion") +
      theme(axis.text.x = element_text(angle = 90))
    p4.0.1
    
    
    recBCts = recBCt %>% group_by(sample,recover) %>% 
      summarise(cell.counts = sum(cell.counts),total.counts = sum(total.counts),
                group = unique(group))
    recBCts = recBCts %>% group_by(sample) %>% 
      summarise(cell.counts = cell.counts, recover = recover,
                cell.prop = cell.counts/total.counts,group = unique(group))
    recBCts[which(recBCts$group == "pandaE11"),"group"] = "E11"
    p4.0.2 = ggbarplot(
      recBCts,
      x = "group", y = "cell.prop", 
      add = c("mean_sd", "jitter"), 
      color = "recover", 
      position = position_dodge(0.8)
    ) + 
      scale_color_jama() + 
      scale_y_continuous(limits = c(0,1)) +
      xlab("") + ylab("Recovered cells proportion")
    p4.0.2
 
    #data stat
    head(recBCts)
    recr = recBCts[which(recBCts$recover == "unrecovered"),]
    recr$cell.prop = 1 - recr$cell.prop
    recr$recover = "recovered"
    recBCts2 = rbind(recBCts,recr)
    recBCts_st = recBCts2 %>% group_by(group,recover) %>% summarise(mean.prop = mean(cell.prop),
                                                                   sd = sd(cell.prop))
    
    
    recBCts_st2 = recBCts2 %>% group_by(recover) %>% summarise(mean.recover.prop = mean(cell.prop),
                                                                    sd = sd(cell.prop))
    recoverttls = list("recover_stat_total" = recr,
                       "recover_stat_cells" = recBCt,
                       "recover_stat_sample_replactes" = recBCts,
                       "recover_stat_sample" = recBCts_st,
                       "recover_stat_total" = recBCts_st2,
                       "recover_barplot" = list("cells" = p4.0.1, "total" = p4.0.2))
    write.csv(recoverttls$recover_stat_sample,file = "final_result/CREST_final_analysis/rebbutal_analysis/recover_stat_1_9.csv",
              quote = F,row.names = F)
    ognrevisels$recoverttls = recoverttls
  }
  
  #clone cell type stat
  {
    CloneCnum = function(edrgl1,title){
      cnumst = as.data.frame(table(rowSums(edrgl1>0)))
      p9.1 = ggbarplot(cnumst,x = "Var1",y = "Freq",stat = "identity",fill = "grey",size = 0.2) + 
        # scale_fill_viridis(discrete = T) +
        geom_text(aes(x = Var1,y = Freq,label = Freq),vjust = -1) +
        NoLegend() + xlab("celltype num") + ylab("clone frequence") +
        ggtitle(title)
      p9.1
      return(p9.1)
    }
    ClonePie = function(edrgl1,title){
      cnumst = as.data.frame(table(rowSums(edrgl1>0)))
      cnumst$prop = round(cnumst$Freq/sum(cnumst$Freq),3)
      
      cnumstp = NULL
      propsum = 0
      for (i in 1:nrow(cnumst)) {
        propsum = propsum + cnumst[i,"prop"]
        cnumstp = rbind(cnumstp,cnumst[i,])
        if(propsum > 0.95){
          cnumstp = rbind(cnumstp,data.frame("Var1" = paste0("> ",cnumst[i,"Var1"]),"Freq" = 0,
                                             "prop" = round(1-propsum,3)))
          break
        }
      }
      cnumstp$group = cnumstp$Var1
      labs <- paste0(cnumstp$group, " (", cnumstp$prop*100, "%)")
      
      p9.2 = ggdonutchart(cnumstp, x = "prop", label = rev(labs), color = "white",
                          lab.pos = "in", lab.font = "black",
                          fill = "group") + NoLegend() + 
        # scale_fill_viridis(discrete = T,option = "E") +
        # scale_fill_manual(values = colorRampPalette(brewer.pal(11, "RdBu"))(nrow(cnumstp))) +
        ggtitle(paste0(title," clonenumber = ",nrow(edrgl1)))
      p9.2
      return(p9.2)
    }
    colorRampPalette(brewer.pal(9, "Blues"))(10)
    calid = names(datals$v1)
    celltype_barplot = list()
    celltype_pieplot = list()
    clonesizedf = NULL
    for (i in 1:length(calid)) {
      cts = dcast(datals$arraydata_bl[[calid[i]]],pattern~Cell.type)
      cts = cts[-1]
      cts$count = rowSums(cts)
      cts = cts[which(cts$count>1),]
      p9.1.1 = CloneCnum(cts[-ncol(cts)],calid[i])
      p9.2.1 = ClonePie(cts[-ncol(cts)],calid[i])
      cts$sample = calid[i]
      clonesizedf = rbind(clonesizedf,cts[c("count","sample")])
      celltype_barplot[[i]] = p9.1.1
      celltype_pieplot[[i]] = p9.2.1
      names(celltype_barplot)[i] = names(celltype_pieplot)[i] = calid[i]
    }
    
    clonecelltypestls = list("celltype_barplot" = celltype_barplot,
                             "celltype_pieplot" = celltype_pieplot)
    clonesizedf$group = unlist(lapply(strsplit(clonesizedf$sample,split = "_"),"[[",1))
    clonesizedf$sample = factor(clonesizedf$sample,
                                levels = c(paste0(rep(c("E11_rep","E15_rep","pandaE11_rep",
                                                      "OGN_rep"),3),
                                                rep(1:3,each = 4)),"E11","E15","pandaE11","OGN"))
    
    p1.3.1 = ggplot(clonesizedf[which(clonesizedf$sample != clonesizedf$group),],
                    aes(y = log2(count), x = sample, fill = group)) + 
      geom_boxplot(width = 0.5,color = "black") + 
      geom_jitter(width = 0.1,size = 1)+
      # stat_compare_means(comparisons = list(c("E11_rep1","E15_rep1"),
      #                                       c("E11_rep2","E15_rep2"),
      #                                       c("E11_rep3","E15_rep3")),
      #                    label = "p.signif",
      #                    method = "t.test") +
      scale_fill_jama() +
      ylab("log2(clone size)") + xlab("") + 
      theme_classic()
    p1.3.1
    
    p1.3.2 = ggplot(clonesizedf[which(clonesizedf$sample %in% c("E11","E15")),],
                    aes(y = log2(count), x = sample, fill = group)) + 
      geom_boxplot(width = 0.5,color = "black") + 
      geom_jitter(width = 0.1,size = 1)+
      stat_compare_means(comparisons = list(c("E11","E15")),
                         label = "p.signif",
                         method = "t.test") +
      # scale_fill_jama() +
      scale_fill_manual(values = c("#CCCCFF","#FFE0E0")) +
      ylab("log2(clone size)") + xlab("") + 
      theme_classic()
    p1.3.2
    
    p1.3.3 = ggplot(clonesizedf[which(clonesizedf$sample %in% c("E11","E15")),],
                    aes(x = log2(count), color = sample, fill = sample)) + 
      geom_density(alpha = 0.7) + 
      scale_fill_manual(values = c("#CCCCFF","#FFE0E0")) +
      # scale_color_manual(values = c("#CCCCFF","#FFE0E0")) +
      ylab("log2(clone size)") + xlab("density") + 
      scale_x_continuous(limits = c(0,log2(max(clonesizedf$count))))+
      theme_classic()
    p1.3.3
    
    clonestls = list("clonecelltypestls" = clonecelltypestls,
                     "clonesizestls" = list("clonesizedf" = clonesizedf,
                                            "sampleboxplot1" = p1.3.1,
                                            "sampleboxplot2" = p1.3.2,
                                            "sampledensityplot" = p1.3.3))
  }
  
  
  scbasicstls$recoverstls = recoverttls
  scbasicstls$clonestls = clonestls
  
  # clone cluster edit size stat
  {
    myclcell = datals$arraydata_bl$E11
    Celleditst = function(myclcell,mytitle){
      tmpv1 = myclcell$pattern.v1
      tmpv2 = myclcell$pattern.v2
      tmpv1 = lapply(strsplit(tmpv1,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      tmpv2 = lapply(strsplit(tmpv2,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      editi = unlist(lapply(tmpv1,length)) + unlist(lapply(tmpv2,length))
      editdf = data.frame("celltype" = myclcell$Cell.type,"editnum" = editi)
      editdf = editdf %>% group_by(celltype) %>% mutate(cellnum = n())
      editdf = editdf[which(editdf$cellnum >10),]
      p2.6 = ggstatsplot::ggbetweenstats(data = editdf,x = "celltype", 
                                         y = "editnum",
                                         point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha
                                                           = 0.4, size = 0, stroke = 0),
                                         pairwise.comparisons = T,pairwise.display = "ns") + 
        ggtitle(mytitle) + xlab("celltype")
      p2.6
      return(p2.6)
    }
    
    CloneCleditst = function(mycluster,mytitle){
      cloneid = names(mycluster)
      tmp = strsplit(x = cloneid,split = "-")
      tmpv1 = unlist(lapply(tmp, "[[",2))
      tmpv2 = unlist(lapply(tmp, "[[",3))
      tmpv1 = lapply(strsplit(tmpv1,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      tmpv2 = lapply(strsplit(tmpv2,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      editi = unlist(lapply(tmpv1,length)) + unlist(lapply(tmpv2,length))
      editdf = data.frame("cluster" = mycluster,"editnum" = editi)
      p2.5 = ggstatsplot::ggbetweenstats(data = editdf,x = "cluster", y = "editnum") + 
        ggtitle(mytitle)
      
      return(p2.5)
    }
    p2.5.1 = CloneCleditst(e15cl$cluster,"E15.5 edit number of different clone cluster")
    p2.5.2 = CloneCleditst(e11cl$cluster,"E11.5 edit number of different clone cluster")
    ggsave(p2.5.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15.5_edit_number_clone_cluster.pdf",
           width = 7,height = 5)
    ggsave(p2.5.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E11.5_edit_number_clone_cluster.pdf",
           width = 7,height = 5)
    
    mydomain = datals$figure$DomainDispersion_f4$E11$pvs_domain
    CloneDomaineditst = function(mydomain,mytitle){
      cloneid = mydomain$group
      tmp = strsplit(x = cloneid,split = "-")
      tmpv1 = unlist(lapply(tmp, "[[",2))
      tmpv2 = unlist(lapply(tmp, "[[",3))
      tmpv1 = lapply(strsplit(tmpv1,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      tmpv2 = lapply(strsplit(tmpv2,split = "_|&"), function(x){x = unique(x); x = x[!x %in% c('NA','NONE')]})
      editi = unlist(lapply(tmpv1,length)) + unlist(lapply(tmpv2,length))
      editdf = data.frame("domain" = mydomain$main.domain,"editnum" = editi)
      p2.6 = ggstatsplot::ggbetweenstats(data = editdf,x = "domain", y = "editnum") + 
        scale_color_manual(values = domaincol) +
        ggtitle(mytitle) + xlab("dominant domain")
      p2.6
      return(p2.6)
    }
    p2.6.1 = CloneDomaineditst(datals$figure$DomainDispersion_f4$E15$pvs_domain,"E15.5 edit number of clones with different dominant domain")
    p2.6.2 = CloneDomaineditst(datals$figure$DomainDispersion_f4$E11$pvs_domain,"E11.5 edit number of clones with different dominant domain")
    
    scbasicstls$clone_edit_numst = list("cluster" = list("E11" = p2.5.1,"E15" = p2.5.2),
                                        "domain" = list("E15" = p2.6.1,"E11" = p2.6.2))
    ggsave(p2.6.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15.5_edit_number_domain.pdf",
           width = 7,height = 5)
    ggsave(p2.6.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E11.5_edit_number_domain.pdf",
           width = 7,height = 5)
    datals$cellorder$E11
    datals$arraydata_bl$E11$Cell.type = factor(datals$arraydata_bl$E11$Cell.type,
                                               levels = c("NPBL", "NPBM", "Rgl1","NbBM0","NbBM1",
                                                          "NbGABA0"  ,"NbGABABL" ,"NbGLUAL1" ,"NbGABABI" ,
                                                          "NbGLUAL2",  "OMTN" ,  "NbFP","GLUAL","Peric"))
    myclcell1 = datals$arraydata_bl$E11
    myclcell2 = datals$arraydata_bl$E15
    myclcell1 = myclcell1[which(myclcell1$Cell.type %in% c("NPBL", "NPBM", "Rgl1")),]
    myclcell2 = myclcell2[which(myclcell2$Cell.type %in% c("Rgl1", "Rgl2AL","Rgl2BP", "Rgl3")),]
    # myclcell2$Cell.type = paste0("E15.5-",myclcell2$Cell.type)
    # myclcell = rbind(myclcell1,myclcell2)
    
    p2.7.1 = Celleditst(myclcell1,"E11.5 edit number of different progenitor cell types")
    p2.7.1
    p2.7.2 = Celleditst(myclcell2,"E15.5 edit number of different progenitor cell types")
    ggsave(p2.7.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E11.5_edit_number_celltype_progenitor.pdf",
           width = 7,height = 5)
    ggsave(p2.7.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15.5_edit_number_celltype_progenitor.pdf",
           width = 9,height = 6)
    
  }
  
  # edit type draw
  {
    #(1) stat barcode edit
    bl = datals$blacklist
    Editstat = function(tedit){
      edit_stat = NULL
      edit = unique(tedit)
      for(i in 1:length(edit)){
        line = edit[i]
        line_size = length(tedit[which(tedit == line)])
        line_strt = unique(unlist(str_split(line, "_")))
        line_str = unique(unlist(str_split(line, "_|&")))
        line_str = line_str[which(line_str != "NONE" & line_str != "NA")]
        if(length(line_str)>0){
          line_del = line_str[grep("D",line_str)]
          line_ins = line_str[grep("I",line_str)]
          line_deln = unlist(str_extract_all(line_del, "[0-9]+D"))
          line_insn = unlist(str_extract_all(line_ins, "[0-9]+I"))
          line_start = c(unlist(str_extract_all(line_del, "\\+[0-9]+")),unlist(str_extract_all(line_ins, "\\+[0-9]+")))
          line_start = as.numeric(substr(line_start,2,nchar(line_start)))
          line_deln = as.numeric(substr(line_deln,1,nchar(line_deln)-1))
          line_insn = as.numeric(substr(line_insn,1,nchar(line_insn)-1))
          if(length(line_start) == length(line_str)){
            edit_stat =rbind(edit_stat,data.frame("array" = line,"clone_size" = line_size,
                                                  "edit_num" = length(line_strt),'width' = c(line_deln,line_insn),"start" = line_start,
                                                  "type" = rep(c("deletion","insertion"),c(length(line_del),length(line_ins)))))
          }
        }
        
      }
      return(edit_stat)
    }
    library(ggExtra)
    arrayeditls = list()
    ClonesizeEditNumPlot = function(edit_stat,mytitle){
      edit_stat = edit_stat[which(edit_stat$type == "deletion"),] %>% group_by(array) %>% 
        summarise(clone_size = unique(clone_size),width = sum(width),edit_num = unique(edit_num))
      edit_stat$log2_clone_size = log2(edit_stat$clone_size)
      edit_stat$edit_num_id = as.character(edit_stat$edit_num)
      edit_stat[which(as.numeric(edit_stat$edit_num_id) >= 13),]$edit_num_id = ">=13"
      edit_stat$edit_num_id = factor(edit_stat$edit_num_id,levels = c(as.character(c(1:13)),">=13"))
      edit_stat$log2_clone_size_id = as.character(edit_stat$log2_clone_size)
      
      if(max(edit_stat$log2_clone_size) > 5){
        maxedit = 5
      }else{
        maxedit = as.integer(max(edit_stat$log2_clone_size))+1
      }
      for (i in 0:maxedit) {
        edit_stat[(edit_stat$log2_clone_size >= i & edit_stat$log2_clone_size < i+1),"log2_clone_size_id"] = paste0(i,"-",i+1)
      }
      edit_stat[(edit_stat$log2_clone_size >= 5),"log2_clone_size_id"] = "5+"
      edit_statmd = edit_stat %>% group_by(log2_clone_size_id) %>% 
        summarise(md = median(edit_num))
      library(trend)
      if(length(edit_statmd$md)>2){
        mkt = mk.test(edit_statmd$md)
        p.value = mkt$p.value
        if(is.na(mkt$p.value)){
          p.value = 1
        }
      }else{
        p.value = NA
      }
      
      edit_stat = edit_stat[order(edit_stat$log2_clone_size_id),]
      edit_stat$log2_clone_size_id = factor(edit_stat$log2_clone_size_id,levels = unique(edit_stat$log2_clone_size_id))
      p2.4.2 = ggviolin(data = edit_stat,x = "log2_clone_size_id", y = "edit_num",
                        draw_quantiles = 0.5) +
        geom_jitter(width = 0.2,size = 0.3,alpha = 0.5,height = 0.05)+
        annotate("text",x = maxedit-1, y = max(edit_stat$edit_num)+5,
                 label = paste0("Mann-Kendall trend test: p =", round(p.value,3))) +
        ggtitle(mytitle) + theme_classic() + 
        xlab("log2 clone size range") + ylab("edit num")
      return(p2.4.2)
      
    }
    datals$arraydata_bl$E11_pandaE11 = rbind(datals$arraydata_bl$E11,datals$arraydata_bl$pandaE11)
    for (i in 1:length(datals$arraydata_bl)) {
      print(i)
      arayi = datals$arraydata_bl[[i]]
      arayi = arayi[which(!is.na(arayi$pattern.v1) & !is.na(arayi$pattern.v2)),]
      editv1 = arayi$pattern.v1
      editv2 = arayi$pattern.v2
      
      # library(ggstatsplot)
      # ggbetweenstats(data = edit_stat,x = "log2_clone_size_id", y = "edit_num")
      edit_stat = Editstat(paste0(editv1,"_",editv2))
      p2.4.1 = ClonesizeEditNumPlot(edit_stat,names(datals$arraydata_bl)[i])
      
      
      arayi = arayi[which(!arayi$pattern.v1 %in% as.character(datals$blacklist$v1_array$Var1) &
                            !arayi$pattern.v2 %in% as.character(datals$blacklist$v2_array$Var1)),]
      
      editv1f = arayi$pattern.v1
      editv2f = arayi$pattern.v2
      edit_statf = Editstat(paste0(editv1f,"_",editv2f))
      p2.4.2 = ClonesizeEditNumPlot(edit_statf,paste0(names(datals$arraydata_bl)[i]," filtered"))
      
      edit_statv1 = Editstat(editv1[which(!editv1 %in% as.character(datals$blacklist$v1_array$Var1))])
      p2.4.3 = ClonesizeEditNumPlot(edit_statv1,paste0(names(datals$arraydata_bl)[i]," V1"))
      
      edit_statv2 = Editstat(editv2[which(!editv2 %in% as.character(datals$blacklist$v2_array$Var1))])
      p2.4.4 = ClonesizeEditNumPlot(edit_statv2,paste0(names(datals$arraydata_bl)[i]," V2"))
      
      arrayeditls[[i]] = list(p2.4.1,p2.4.2,p2.4.3,p2.4.4)
      names(arrayeditls)[i] = names(datals$arraydata_bl)[i]
    }
    scbasicstls$editnum_clonesize = arrayeditls
    # dir.create("final_result/CREST_final_analysis/rebbutal_analysis")
    arrayeditls$E11[[1]]
    arrayeditls$E11_pandaE11[[1]]
    arrayeditls$E11[[1]]
    ggsave(arrayeditls$E15,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15_edit_num_clone_size_violin.pdf",
           width = 10,height = 6)
    ggexport(arrayeditls,filename = "final_result/CREST_final_analysis/rebbutal_analysis/edit_num_clone_size_violin_list.pdf",
             width = 8,height = 6)
    #(2) edit pos and indel ribbon
    #segment diversity
    StatIndelLen = function(v1ls,v2ls){
      indelen = NULL
      indelseg = list()
      for (i in 1:length(v1ls)) {
        v1 = v1ls[[i]]
        v2 = v2ls[[i]]
        v1indelseg = ScartoIndel(unique(v1$pattern))
        v2indelseg = ScartoIndel(unique(v2$pattern))
        v1indelseg = merge(v1indelseg,v1[c("Cell.BC","pattern")],by.x = "array",by.y = "pattern")
        v2indelseg = merge(v2indelseg,v2[c("Cell.BC","pattern")],by.x = "array",by.y = "pattern")
        v1indelseg$sample = v2indelseg$sample = names(v1ls)[i]
        v1indelseg$group = "V1";v2indelseg$group = "V2"
        line2 = rbind(v1indelseg,v2indelseg)
        indelseg[[i]] = line2
        names(indelseg)[i] = names(v1ls)[i]
        # v1 = v1[v1$umim != "NONE_NONE_NONE_NONE_NONE_NONE_NONE",]
        # v2 = v2[v2$umim != "NONE_NONE_NONE_NONE_NONE_NONE_NONE",]
        v1spl = strsplit(v1$pattern,"_|&")
        v1taglen = unlist(lapply(v1spl, function(x){length(unique(x[which(x!="NONE")]))}))
        v2spl = strsplit(v2$pattern,"_|&")
        v2taglen = unlist(lapply(v2spl, function(x){length(unique(x[which(x!="NONE")]))}))
        line = data.frame(sample = names(v1ls)[i],indelenth = c(v1taglen,v2taglen),
                          group = rep(c("V1","V2"),c(length(v1taglen),length(v2taglen))))
        
        indelen = rbind(indelen,line)
        
        
      }
      return(list(indelen,indelseg))
      
    }
    ScartoIndel = function(tl){
      indel = NULL
      #build indel data frame
      for(i in c(1:length(tl))){
        print(i)
        x = tl[i]
        x_str = strsplit(x,split = "_")[[1]]
        for (j in 1:length(x_str)) {
          x_str_seg = strsplit(x_str[j],split = "&")[[1]]
          segid = j
          # if(x_str_seg[1] != "NONE"){
          if(x_str_seg[1] != "NONE"){
            for (k in 1:length(x_str_seg)) {
              # x_pre_num = str_extract_all(x_str_seg[k], "[0-9]+")
              # x_pre_num = as.numeric(x_pre_num[[1]])
              # start = x_pre_num[length(x_pre_num)]
              # end = start + x_pre_num[1]
              # segid.end = min(cutsite[which(cutsite$V3>end),"V1"])
              if(length(x_str_seg[k][grep("I",x_str_seg[k])])>0){
                indeltype = "Insertion"
              }else if(j == 7){
                if(x_str_seg[k] %in% strsplit(x_str[j-1],split = "&")[[1]]){
                  indeltype = "MultipleDeletion"
                }else{
                  indeltype = "SingleDeletion"
                }
                
              }else if(j == 1){
                if(x_str_seg[k] %in% strsplit(x_str[j+1],split = "&")[[1]]){
                  indeltype = "MultipleDeletion"
                }else{
                  indeltype = "SingleDeletion"
                }
                
              }else if(x_str_seg[k] %in% strsplit(x_str[j+1],split = "&")[[1]] |
                       x_str_seg[k] %in% strsplit(x_str[j-1],split = "&")[[1]]){
                indeltype = "MultipleDeletion"
              }else{
                indeltype = "SingleDeletion"
              }
              
              indelline = data.frame(id = j,array = x,indel = x_str_seg[k],indeltype = indeltype)
              indel = rbind(indel,indelline)
              
            }
          }else{
            indelline = data.frame(id = segid,array = x,indel = "NONE",indeltype = "Unedit")
            indel = rbind(indel,indelline)
          }
          
          
        }
      }
      return(indel)
    }
    
    #build data
    {
      # bulkls$v1$Dual_E15.5_1a
      # v1lsx = list()
      # v1lsx$E15_1 = rbind(bulkls$v1$Dual_E15.5_1a,bulkls$v1$Dual_E15.5_1b)
      # v1lsx$E15_2 = rbind(bulkls$v1$Dual_E15.5_2a,bulkls$v1$Dual_E15.5_2b)
      # v1lsx$E15_3 = rbind(bulkls$v1$Dual_E15.5_3a,bulkls$v1$Dual_E15.5_3b)
      # v1lsx$E15_1$pattern = v1lsx$E15_1$umim;v1lsx$E15_2$pattern = v1lsx$E15_2$umim;v1lsx$E15_3$pattern = v1lsx$E15_3$umim;
      # 
      # v2lsx = list()
      # v2lsx$E15_1 = rbind(bulkls$v2$Dual_E15.5_1a,bulkls$v2$Dual_E15.5_1b)
      # v2lsx$E15_2 = rbind(bulkls$v2$Dual_E15.5_2a,bulkls$v2$Dual_E15.5_2b)
      # v2lsx$E15_3 = rbind(bulkls$v2$Dual_E15.5_3a,bulkls$v2$Dual_E15.5_3b)
      # v2lsx$E15_1$pattern = v2lsx$E15_1$umim;v2lsx$E15_2$pattern = v2lsx$E15_2$umim;v2lsx$E15_3$pattern = v2lsx$E15_3$umim;
      # 
      # v1ls = list()
      # v2ls = list()
      # for (i in 1:length(datals$v1)) {
      #   v1ls[[i]] = datals$v1[[i]]$resultraw
      #   v2ls[[i]] = datals$v2[[i]]$resultraw
      #   names(v1ls)[i] = names(v2ls)[i] = names(datals$v1)[i]
      # }
      # v1lstmp = list("E11_rep3" = v1ls$E11_rep3)
      # v2lstmp = list("E11_rep3" = v2ls$E11_rep3)
      # indelsegtmp = StatIndelLen(v1lstmp,v2lstmp)
      # indelseg[[2]] = append(indelseg[[2]], indelsegtmp[[2]])
      # indelseg = StatIndelLen(v1ls,v2ls)
      # qsave(indelseg,file = "final_result/CREST_final_analysis/rebbutal_analysis/indelseg.qs")
      # indelsegbulk = StatIndelLen(v1lsx,v2lsx)
      # qsave(indelsegbulk,file = "final_result/CREST_final_analysis/rebbutal_analysis/indelseg.qs")
      # indelseg = readRDS("final_result/xtotal/indelseg_raw.rds")
      # indelseg[[2]] = append(indelseg[[2]], indelsegbulk[[2]])
      qsave(indelseg,file = "final_result/CREST_final_analysis/rebbutal_analysis/indelseg.qs")
    }
    
    indelsegi = indelseg[[2]]
    totalseg = rbind(indelsegf$E11_rep1,indelsegf$E11_rep2,indelsegf$E11_rep3,
                     indelsegf$E15_rep1,indelsegf$E15_rep2,indelsegf$E15_rep3,
                     indelsegf$pandaE11_rep1,indelsegf$pandaE11_rep2,indelsegf$pandaE11_rep3,
                     indelsegf$OGN_rep1,indelsegf$OGN_rep2,indelsegf$OGN_rep3)
    # totalsegf = totalseg[which(!totalseg$array %in% c(as.character(bl$v1_array$Var1), as.character(bl$v2_array$Var1))),]
    totalseg = totalseg[which(!(totalseg$group == "V1" &totalseg$id == 4)),]
    totalseg[which((totalseg$group == "V1" & totalseg$id > 4)),"id"] = totalseg[which((totalseg$group == "V1" & totalseg$id > 4)),"id"] - 1
    
    
    # for (i in 1:length(indelsegf)) {
    #   indelsegf[[i]] = indelsegf[[i]][which(!indelsegf[[i]]$array %in% c(as.character(bl$v1_array$Var1), as.character(bl$v2_array$Var1))),]
    # }
    # 
    indelseg_del4 = indelseg[[2]]
    # indelseg_del4f = list()
    for (i in 1:length(indelseg_del4)) {
      indelseg_del4[[i]] = indelseg_del4[[i]][which(!(indelseg_del4[[i]]$group == "V1" &
                                                        indelseg_del4[[i]]$id == 4)),]
      indelseg_del4[[i]][which((indelseg_del4[[i]]$group == "V1" &
                                  indelseg_del4[[i]]$id > 4)),"id"] = indelseg_del4[[i]][which((indelseg_del4[[i]]$group == "V1" &
                                                                                                  indelseg_del4[[i]]$id > 4)),"id"] - 1
      # indelseg_del4f[[i]] = indelseg_del4[[i]][which(!indelseg_del4[[i]]$array %in% c(as.character(bl$v1_array$Var1), as.character(bl$v2_array$Var1))),]
    }
    # names(indelseg_del4f) = names(indelseg_del4)
    # indelseg = StatIndelLen(v1ls,v2ls)
    
    #bulktype stat
    {
      # oldbulkls = qread("bulk_result/X.bulk_in_vivo/bulk_total_12_9.qs")
      # bulkls$v1$nkx1 = oldbulkls$v1$nkx1
      # bulkls$v2$nkx1 = oldbulkls$v2$nkx1
      # bulkls$v1$nkx2 = oldbulkls$v1$nkx2
      # bulkls$v2$nkx2 = oldbulkls$v2$nkx2
      # bulkls$v1$nkx3 = oldbulkls$v1$nkx3
      # bulkls$v2$nkx3 = oldbulkls$v2$nkx3
      v1lsxb = bulkls$v1
      v2lsxb = bulkls$v2
      for (i in 1:length(v1lsxb)) {
        v1lsxb[[i]]$main = unlist(lapply(strsplit(v1lsxb[[i]]$main,"_"),function(x){paste0(x[-4],collapse = "_")}))
      }
      
      for (i in 1:length(v1lsxb)) {
        v1lsxb[[i]]$pattern = v1lsxb[[i]]$main
        v2lsxb[[i]]$pattern = v2lsxb[[i]]$main
      }
      indelsegbulk = StatIndelLen(v1lsxb,v2lsxb)
      qsave(indelsegbulk,file = "final_result/CREST_final_analysis/rebbutal_analysis/indelsegbulk.qs")
      indelsegbulki = do.call(rbind,indelsegbulk[[2]])
      indelsegbulki = list("En1_CREST_bulk" = indelsegbulki)
      p2.3.3 = EditTypeStat(indelsegbulki,barpos = "fill")
      p2.3.3[[1]]
      ggexport(p2.3.3[[1]],filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_diversity_stat_bulktotal.pdf",width = 6,height = 4)
      
    }
    
    EditTypeStat = function(indelseg,barpos = "stack"){
      p2.3 = list()
      for (i in 1:length(indelseg)) {
        line = indelseg[[i]]
        line$indeltype = factor(line$indeltype,levels = c("Unedit","Insertion","MultipleDeletion",
                                                          "SingleDeletion"))
        # line = line[which(line$indel != 'NONE'),]
        segidden = line %>% group_by(group,id,indeltype) %>% summarise(diversity = length(unique(indel)))
        
        p2.3.0 = ggplot(segidden,aes(x = id, y = diversity,fill = indeltype)) + 
          geom_bar(color = "black",stat = "identity",position = barpos) + facet_grid(~group,scales = "free") + 
          theme_pubr() + scale_x_continuous(breaks = c(unique(indelseg[[i]]$id))) +
          scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127"))+ xlab("") +
          ggtitle(names(indelseg)[i])
        p2.3.0
        p2.3[[i]] = p2.3.0
      }
      return(p2.3)
    }
    p2.3.1 = EditTypeStat(indelsegi,barpos = "fill")
    # p2.3.2 = EditTypeStat(indelsegf,barpos = "fill")
    p2.3.2 = EditTypeStat(indelseg_del4,barpos = "fill")
    # p2.3.4 = EditTypeStat(indelseg_del4f,barpos = "fill")
    names(p2.3.1) = names(p2.3.2) = names(indelsegi)
    
    #total seg
    {
      totalseg$indeltype = factor(totalseg$indeltype,levels = c("Unedit","Insertion","MultipleDeletion",
                                                        "SingleDeletion"))
      segidden = totalseg %>% group_by(group,id,indeltype) %>%summarise(diversity = length(unique(indel)))
      
      p2.4.1 = ggplot(segidden,aes(x = id, y = diversity,fill = indeltype)) + 
        geom_bar(color = "black",stat = "identity",position = "fill") + facet_grid(~group,scales = "free") + 
        theme_pubr() + scale_x_continuous(breaks = c(unique(totalseg$id))) +
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127"))+ xlab("") +
        ggtitle("total")
      p2.4.1
      
      totalseg$sample_group = unlist(lapply(strsplit(totalseg$sample,split = "_"),"[[",1))
      segidden2 = totalseg %>% group_by(sample,group,id,indeltype) %>%
        summarise(diversity = length(unique(indel)),sample_group = unique(sample_group))
      segidden2 = segidden2 %>% group_by(sample,group,id) %>% mutate(diversity.prop = diversity/ sum(diversity))
      p2.4.2 = ggbarplot(
        segidden2, x = "id", y = "diversity.prop", facet.by = c("group","sample_group"),
        scales = "free",
        add = c("mean_sd"), 
        # color = "indeltype", 
        fill = "indeltype", 
        position = position_dodge(0.8)
      ) + 
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127")) +
        xlab("") + ylab("") 
        # scale_x_continuous(breaks = c(1:8))
      p2.4.2
      
      
      segidden3 = totalseg %>% group_by(sample_group,group,id,indeltype) %>%
        summarise(diversity = length(unique(indel)))
      
      p2.4.3 = ggbarplot(
        segidden3, x = "id", y = "diversity", 
        facet.by = c("group","sample_group"),
        scales = "free",
        # color = "indeltype", 
        fill = "indeltype", 
        position = position_fill()
      ) + 
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127")) +
        xlab("") + ylab("") 
      # scale_x_continuous(breaks = c(1:8))
      p2.4.3
      
      #E11 prog compare
      progsegv1 = datals$v1$E11$result_withcellan[which(datals$v1$E11$result_withcellan$Cell.type %in% 
                                              c("NPBM","NPBL","Rgl1")),]
      progsegv2 = datals$v2$E11$result_withcellan[which(datals$v2$E11$result_withcellan$Cell.type %in% 
                                              c("NPBM","NPBL","Rgl1")),]
      
      # progsegv1 = datals$v1$E11$result_withcellan
      # progsegv2 = datals$v2$E11$result_withcellan
      
      progseg = rbind(progsegv1,progsegv2)
      progseg = merge(totalseg[which(totalseg$sample_group == "E11"),],progseg, by.x = "array",by.y = "pattern")
      
      segidden4 = progseg %>% group_by(sample,group,id,indeltype,Cell.type) %>%
        summarise(diversity = length(unique(indel)))
      segidden4 = segidden4 %>% group_by(sample,group,id,Cell.type) %>% 
        mutate(diversity.prop = diversity / sum(diversity))
      p2.4.4 = ggbarplot(
        segidden4, x = "id", y = "diversity.prop", facet.by = c("group","Cell.type"),
        scales = "free",
        add = c("mean_sd"), 
        # color = "indeltype", 
        fill = "indeltype", 
        position = position_dodge(0.8)
      ) + 
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127")) +
        xlab("") + ylab("") 
      # scale_x_continuous(breaks = c(1:8))
      p2.4.4
      
      segidden5 = progseg %>% group_by(Cell.type,group,id,indeltype) %>%
        summarise(diversity = length(unique(indel)))
      p2.4.5 = ggbarplot(
        segidden5, x = "id", y = "diversity", 
        facet.by = c("group","Cell.type"),
        scales = "free",
        # color = "indeltype", 
        fill = "indeltype", 
        position = position_fill()
      ) + 
        scale_fill_manual(values = c("white","#4391C5","#F69173","#CC2127")) +
        xlab("") + ylab("")  
      # scale_x_continuous(breaks = c(1:8))
      p2.4.5
      
      segidden6 = unique(totalseg[-1])
      segidden6 = segidden6 %>% group_by(sample_group,sample,group) %>%
        summarise(arraynum = length(unique(array)),
               multiarray = length(unique(array[which(indeltype == "MultipleDeletion")])),
               multiarrayprop = multiarray/arraynum,
               indelnum = length(unique(indel)),
               multiindel = length(unique(indel[which(indeltype == "MultipleDeletion")])),
               multiindelprop = multiindel/indelnum
        )
    
      
      p2.4.6 = ggbarplot(
        segidden6, x = "sample_group", y = "multiarrayprop", 
        scales = "free",
        add = c("mean_sd"), 
        # color = "indeltype", 
        fill = "group", 
        position = position_dodge(0.8)
      ) + 
        scale_fill_jama() +
        stat_compare_means(aes(group = group),
                           label = "p.signif",
                           method = "t.test") +
        xlab("") + ylab("Ritio of arrays with multiple site deletion") 

      p2.4.6
      
      
    }
    
    indelsegtypels = list("indelseg"= indelseg_del4, "totalseg" = totalseg,
                          "indelseg_figurels" = p2.3.1,
                          "indelseg_figurels_delv14" = p2.3.2,
                          "indelseg_total_sample" = list("total" = p2.4.1,"dodge_errorbar" = p2.4.2,"fill" = p2.4.3),
                           "indelseg_E11_prog" = list("dodge_errorbar" = p2.4.4,"fill" = p2.4.5),
                          "indelseg_multipleedit_v1v2" = p2.4.6)
    
    # ggexport(p2.3.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_diversity_stat.pdf",width = 6,height = 4)
    # ggexport(p2.3.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_diversity_stat_bl.pdf",width = 6,height = 4)
    # ggexport(p2.3.3,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_diversity_stat_del4.pdf",width = 6,height = 4)
    # ggexport(p2.3.4,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_diversity_stat_del4_bl.pdf",width = 6,height = 4)
    # 
    #
    
    #(3) plot ribbon
    {
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
      segrefv1 = read.table("raw_data/ref/Cutsite CREST/V1.cutSites")
      segrefv1 = segrefv1[-1,]
      scarrefv1 = ReadCutsite(segrefv1)
      
      segrefv2 = read.table("raw_data/ref/Cutsite CREST/V2.cutSites")
      segrefv2 = segrefv2[-1,]
      scarrefv2 = ReadCutsite(segrefv2)
      e11vf = BlacklistFilter2(e11v1$result,e11v2$result,bl$v1_array,bl$v2_array)
      INDELRibbon = function(pattern,bl,scarrefv1){
        edit_stat_e11 = Editstat(pattern)
        ind_allsite_per = NULL
        for (i in 1:nrow(edit_stat_e11)) {
          print(i)
          line = edit_stat_e11[i,]
          line = line[rep(1,line$width),]
          line$pos = unique(line$start):(unique(line$start)+unique(line$width)-1)
          # linerep = line[rep(1:nrow(line),unique(line$clone_size)),]
          ind_allsite_per = rbind(ind_allsite_per,line)
        }
        
        ind_allsite_sum_t = ind_allsite_per %>% 
          group_by(pos,type) %>% summarise(num = sum(clone_size))
        max_editnum = length(pattern)
        ind_allsite_sum_t$pop = ind_allsite_sum_t$num/max_editnum
        
        pt = ggplot()+geom_ribbon(data = scarrefv1,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
          geom_line(data = ind_allsite_sum_t,aes(x = pos,y = pop,color = type),size=1) +
          scale_x_continuous(breaks=c()) +
          xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
        pt
        
        ind_allsite_sum = ind_allsite_per[which(!ind_allsite_per$array %in% bl$Var1),] %>% 
          group_by(pos,type) %>% summarise(num = sum(clone_size))
        max_editnum = length(pattern[which(!pattern %in% bl$Var1)])
        ind_allsite_sum$pop = ind_allsite_sum$num/max_editnum
        
        pf = ggplot()+geom_ribbon(data = scarrefv1,aes(x=scar,ymin=0,ymax=1,fill=type),alpha=0.1) +
          geom_line(data = ind_allsite_sum,aes(x = pos,y = pop,color = type),size=1) +
          scale_x_continuous(breaks=c()) +
          xlab("")+ylab("Edit proportion")+theme_bw()+theme(legend.position = "none")
        pf
        
        return(list("ribbon" = pt, "ribbon_fil" = pf,"ind_allsite_per" = ind_allsite_per))
      }
      
      ribpl = list()
      ribplf = list()
      for (i in 1:length(v1ls)) {
        v1rib = INDELRibbon(v1ls[[i]]$pattern,bl$v1_array,scarrefv1)
        v2rib = INDELRibbon(v2ls[[i]]$pattern,bl$v2_array,scarrefv2)
        ribpl[[i]] = v1rib$ribbon + ggtitle(paste0(names(v1ls)[i]," V1")) + v2rib$ribbon + ggtitle(paste0(names(v1ls)[i]," V2"))
        ribplf[[i]] = v1rib$ribbon_fil + ggtitle(paste0(names(v1ls)[i]," V1"))  + v2rib$ribbon_fil + ggtitle(paste0(names(v1ls)[i]," V2"))
      }
      ribbonls = list("v1" = v1rib,"v2" = v2rib)
      rebuttalst$ribbonls = ribbonls 
      ggexport(ribpl,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_cutsite_stat.pdf",width = 8,height = 4)
      ggexport(ribplf,filename = "final_result/CREST_final_analysis/rebbutal_analysis/segment_indel_cutsite_stat_blfiltered.pdf",width = 8,height = 4)
      
    }
    
  }
  scbasicstls$indelsegtypels = indelsegtypels
  
  # stat text data
  {
    nrow(datals$blacklist$v1_array)
    nrow(datals$blacklist$v2_array)
    blfiltst = NULL
    for (i in 1:length(datals$v1)) {
      samplei = names(datals$v1)[i]
      v1i = datals$v1[[i]]$resultraw
      v2i = datals$v2[[i]]$resultraw
      v1ibnum = v1i[which(v1i$pattern %in% bl$v1_array$Var1),]
      v2ibnum = v2i[which(v2i$pattern %in% bl$v2_array$Var1),]
      cmbi = length(unique(datals$arraydata[[samplei]]$pattern))
      cmbibnum = cmbi - length(unique(datals$arraydata_bl[[samplei]]$pattern))
      blfiltst = rbind(blfiltst,data.frame("sample" = samplei,
                                           "v1" = length(unique(v1i$pattern)),
                                           "v1_bl" = length(unique(v1ibnum$pattern)),
                                           "v2" = length(unique(v2i$pattern)),
                                           "v2_bl" = length(unique(v2ibnum$pattern)),
                                           "combine" = cmbi,
                                           "combine_bl" = cmbibnum))
    }
    blfiltst$v1 = paste0(blfiltst$v1_bl,"/",blfiltst$v1)
    blfiltst$v2 = paste0(blfiltst$v2_bl,"/",blfiltst$v2)
    blfiltst$combine = paste0(blfiltst$combine_bl,"/",blfiltst$combine)
    blfiltst = blfiltst[-c(3,5,7)]
    write.csv(blfiltst,file = "final_result/CREST_final_analysis/rebbutal_analysis/blacklistfilterst.csv",
                row.names = F,quote = F)
    
    
    #edit data stat
    edit_statt1 = Editstat(c(datals$v1$E11$resultraw$pattern,datals$v1$E15$resultraw$pattern,
                             datals$v1$pandaE11$resultraw$pattern,datals$v1$OGN$resultraw$pattern))
    edit_statt2 = Editstat(c(datals$v2$E11$resultraw$pattern,datals$v2$E15$resultraw$pattern,
                             datals$v2$pandaE11$resultraw$pattern,datals$v2$OGN$resultraw$pattern))
    summary(unique(edit_statt1[which(edit_statt1$type == "deletion"),])$width)
    summary(unique(edit_statt1[which(edit_statt1$type == "insertion"),])$width)
    summary(unique(edit_statt1[,c(1,3)])$edit_num)
    table(edit_statt1$type)/nrow(edit_statt1)
    summary(unique(edit_statt2[which(edit_statt2$type == "deletion"),])$width)
    summary(unique(edit_statt2[which(edit_statt2$type == "insertion"),])$width)
    summary(unique(edit_statt2[,c(1,3)])$edit_num)
    table(edit_statt2$type)/nrow(edit_statt2)
    
    #qc cell stat table
    head(datals$statdf)
    qcstat = dcast(datals$statdf,sample+group~process,value.var = "count")
    transst = NULL
    for (i in 1:length(cellls)) {
      transst = rbind(transst,data.frame("sample" = names(cellls)[i],
                                         "transcell" = nrow(cellls[[i]])))
    }
    qcstat = merge(qcstat,transst)
    qcstat[which(qcstat$group == "v1"),"v1_v2_combine"] = qcstat[which(qcstat$group == "v1_v2"),"v1_v2_combine"]
    qcstat[which(qcstat$group == "v2"),"v1_v2_combine"] = qcstat[which(qcstat$group == "v1_v2"),"v1_v2_combine"]
    qcstat[which(qcstat$group == "v1"),"v1_v2_combine_blfil"] = qcstat[which(qcstat$group == "v1_v2"),"v1_v2_combine_blfil"]
    qcstat[which(qcstat$group == "v2"),"v1_v2_combine_blfil"] = qcstat[which(qcstat$group == "v1_v2"),"v1_v2_combine_blfil"]
    qcstat = qcstat[which(qcstat$group != "v1_v2"),]
    head(qcstat)
    write.csv(qcstat,file = "final_result/CREST_final_analysis/rebbutal_analysis/qcstat.csv",
              row.names = F,quote = F)
    qcstat = read.csv("final_result/CREST_final_analysis/rebbutal_analysis/qcstat.csv")
    
    
    #stat supp table
    {
      names(datals$arraydata_bl)
      
      #single cell
      BarcodeStat = function(arrayi){
        arrayi = arrayi[which(!arrayi %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                             "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"))]
        arrayi = arrayi[!is.na(arrayi)]
        indi2 = strsplit(split = "_|&",arrayi)
        arl = unlist(lapply(indi2 , function(x){
          x = unique(x)
          xd = x[grep("D",x)]
          xdn = as.numeric(unlist(lapply(strsplit(xd,"D"),"[[",1)))
          return(sum(xdn))
        }))
        
        indi3 = strsplit(split = "_|&",unique(arrayi))
        arnum = unlist(lapply(indi3, function(x){
          x = unique(x)
          return(length(x[which(x != "NONE")]))
        }))
        
        return(c("meanlen" = mean(arl), "meanedit" = mean(arnum),"medianedit" = median(arnum)))

      }
      
      suptable = NULL
      for (i in 1:12) {
        sampleid = names(datals$arraydata)[i]
        #barcode stat
        bcst.v1raw = nrow(datals$v1[[i]]$result_withcellan)
        bcst.v2raw = nrow(datals$v2[[i]]$result_withcellan)
        bcst.v1_v2 = length(intersect(datals$v1[[i]]$result_withcellan$Cell.BC,
                                      datals$v2[[i]]$result_withcellan$Cell.BC))
        bcst.v1orv2 = nrow(datals$arraydata[[i]])
        bcst.v1orv2f = nrow(datals$arraydata_bl[[i]])
        
        #average UMI
        avumi1 = datals$v1[[i]]$resultraw
        avumi1 = avumi1[which(avumi1$Cell.BC %in% datals$v1[[i]]$result_withcellan$Cell.BC),]
        avumi1 = mean(avumi1$umi_num)
        
        avumi2 = datals$v2[[i]]$resultraw
        avumi2 = avumi2[which(avumi2$Cell.BC %in% datals$v2[[i]]$result_withcellan$Cell.BC),]
        avumi2 = mean(avumi2$umi_num)
        
        #barcode proportion
        v1st = BarcodeStat(datals$arraydata_bl[[i]]$pattern.v1)
        v2st = BarcodeStat(datals$arraydata_bl[[i]]$pattern.v2)
        v1st["meanlen"] = 233 - v1st["meanlen"]
        v2st["meanlen"] = 212 - v2st["meanlen"]
        
        #edit type stat
        editst = unique(datals$figure$scbasicstls$indelsegtypels$indelseg[[sampleid]][-1])
        editst = editst %>% group_by(group,indeltype) %>%
          summarise(arraynum = n())
        editst = editst %>% group_by(group) %>% mutate(prop = arraynum/ sum(arraynum))
        
        #clone number
        clst = datals$arraydata_bl[[i]] %>% group_by(pattern) %>% summarise(cellnum = n())
        clnum = nrow(clst)
        clnumf = nrow(clst[which(clst$cellnum > 1),])
        
        suptable = rbind(suptable, data.frame("V1" = bcst.v1raw , "V2" = bcst.v2raw ,
                                              "both V1 and V2" = bcst.v1_v2,
                                              "V1 or V2" = bcst.v1orv2,
                                              "V1 or V2 fil" = bcst.v1orv2f,
                                              "average length" = paste0(round(v1st["meanlen"],3),"/",
                                                                        round(v2st["meanlen"],3)),
                                              "average UMI" = paste0(round(avumi1,3),"/",
                                                                     round(avumi2,3)),
                                              "mean/median number of edits" = paste0("V1:",round(v1st["meanedit"],3),"/",
                                                                                     round(v1st["medianedit"],3),
                                                                                     ",","V2:",round(v2st["meanedit"],3),"/",
                                                                                     round(v2st["medianedit"],3)),
                                              "single deletion" = paste0(round(unlist(editst[3,4]),3)*100, "/",
                                                                         round(unlist(editst[7,4]),3)*100),
                                              "multi-deletion" = paste0(round(unlist(editst[2,4]),3)*100, "/",
                                                                        round(unlist(editst[6,4]),3)*100),
                                              "insertion" = paste0(round(unlist(editst[1,4]),3)*100, "/",
                                                                   round(unlist(editst[5,4]),3)*100),
                                              "Clone (V1+V2)" = clnum,
                                              "total clone(size > 1)" = clnumf
                                              ))
        
        
      }
      
    }
    
    #bulk st
    {
      # qsave(bulkls,file = "bulk_result/X.bulk_in_vivo/bulk_total_1_5.qs")
      suptablebulk = NULL
      bulkls = qread("bulk_result/X.bulk_in_vivo/bulk_total_1_5.qs")
      indelsegbulk = qread("final_result/CREST_final_analysis/rebbutal_analysis/indelsegbulk.qs")
      bulkls$v1$Dual_E10.5_1
      for (i in 1:length(bulkls$v1)) {
        sampleid = names(bulkls$v1)[i]
        v1i = bulkls$v1[[i]]
        v1i = v1i[which(v1i$reads_pro.main >= 0.5 & v1i$reads_num >= 3),]
        v1i = v1i[which(!v1i$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                         "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
        v2i = bulkls$v2[[i]]
        v2i = v2i[which(v2i$reads_pro.main >= 0.5 & v2i$reads_num >= 3),]
        v2i = v2i[which(!v2i$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                         "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
        
        #average UMI
        avumi1 = mean(v1i$reads_num)
        avumi2 = mean(v2i$reads_num)
        
        #barcode proportion
        v1st = BarcodeStat(v1i$main)
        v2st = BarcodeStat(v2i$main)
        v1st["meanlen"] = 233 - v1st["meanlen"]
        v2st["meanlen"] = 212 - v2st["meanlen"]
        
        #edit type stat
        editst = unique(indelsegbulk[[2]][[sampleid]][-c(2,5)])
        editst = editst %>% group_by(group,indeltype) %>%
          summarise(arraynum = n())
        editst = editst %>% group_by(group) %>% mutate(prop = arraynum/ sum(arraynum))
        
        #clone number
        # clst = v2i %>% group_by(pattern) %>% summarise(cellnum = n())
        # clnum = nrow(clst)
        # clnumf = nrow(clst[which(clst$cellnum > 1),])
        
        suptablebulk = rbind(suptablebulk, data.frame("V1" = NA, "V2" = NA,
                                              "both V1 and V2" = NA,
                                              "V1 or V2" = NA,
                                              "V1 or V2 fil" = NA,
                                              "average length" = paste0(round(v1st["meanlen"],3),"/",
                                                                        round(v2st["meanlen"],3)),
                                              "average UMI" = paste0(round(avumi1,3),"/",
                                                                     round(avumi2,3)),
                                              "mean/median number of edits" = paste0("V1:",round(v1st["meanedit"],3),"/",
                                                                                     round(v1st["medianedit"],3),
                                                                                     ",","V2:",round(v2st["meanedit"],3),"/",
                                                                                     round(v2st["medianedit"],3)),
                                              "single deletion" = paste0(round(unlist(editst[3,4]),3)*100, "/",
                                                                         round(unlist(editst[7,4]),3)*100),
                                              "multi-deletion" = paste0(round(unlist(editst[2,4]),3)*100, "/",
                                                                        round(unlist(editst[6,4]),3)*100),
                                              "insertion" = paste0(round(unlist(editst[1,4]),3)*100, "/",
                                                                        round(unlist(editst[5,4]),3)*100),
                                              "Clone (V1+V2)" = NA,
                                              "total clone(size > 1)" = NA
        ))
        
        
      }
      
    }
    suptabletotal = rbind(suptable,suptablebulk)
    suptabletotal$sample = c(names(datals$arraydata)[1:12],names(bulkls$v1))
    write.csv(suptabletotal,file = "final_result/CREST_final_analysis/rebbutal_analysis/suptable.csv")
  }
  
  dir.create("final_result/CREST_final_analysis/E11_E15_basic_stat")
  qsave(scbasicstls,file = "final_result/CREST_final_analysis/E11_E15_basic_stat/scbasicstls.qs")
  datals$figure$scbasicstls = scbasicstls
  qsave(datals,file = "final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_19_edit.qs")
  
  #stat sc occurence
  StatArrayOccSc = function(v1lsxsc){
    arraydft = NULL
    for (i in 1:length(v1lsxsc)) {
      arraydf = as.data.frame(table(v1lsxsc[[i]]$resultraw$pattern))
      arraydf = arraydf[order(-arraydf$Freq),]
      arraydf$Freq = arraydf$Freq/sum(arraydf$Freq)
      
      colnames(arraydf)[2] = names(v1lsxsc)[i]
      if(i == 1){
        arraydft = arraydf
      }else{
        arraydft = merge(arraydft, arraydf, by = "Var1", all = T)
      }
    }
    
    arraydft[is.na(arraydft)] = 0
    ayoc = as.data.frame(table(rowSums(arraydft[-1]>0)))
    ayoc$norm = ayoc$Freq/sum(ayoc$Freq)
    return(ayoc)
  }
  StatArrayOccScDual = function(arylsxsc){
    arraydft = NULL
    for (i in 1:length(arylsxsc)) {
      aryl = arylsxsc[[i]][which(!is.na(arylsxsc[[i]]$pattern.v1) & !is.na(arylsxsc[[i]]$pattern.v2)),]
      arraydf = as.data.frame(table(paste0(aryl$pattern.v1,"-",
                                           aryl$pattern.v2)))
      arraydf = arraydf[order(-arraydf$Freq),]
      arraydf$Freq = arraydf$Freq/sum(arraydf$Freq)
      
      colnames(arraydf)[2] = names(arylsxsc)[i]
      if(i == 1){
        arraydft = arraydf
      }else{
        arraydft = merge(arraydft, arraydf, by = "Var1", all = T)
      }
    }
    
    arraydft[is.na(arraydft)] = 0
    ayoc = as.data.frame(table(rowSums(arraydft[-1]>0)))
    ayoc$norm = ayoc$Freq/sum(ayoc$Freq)
    return(ayoc)
  }

  v1lsxsc = datals$v1[c(1:9)]
  v1lsxsc[[4]]$resultraw = rbind(v1lsxsc[[4]]$resultraw, datals$v1[[10]]$resultraw)
  v1lsxsc[[5]]$resultraw = rbind(v1lsxsc[[5]]$resultraw, datals$v1[[11]]$resultraw)
  v1lsxsc[[6]]$resultraw = rbind(v1lsxsc[[6]]$resultraw, datals$v1[[12]]$resultraw)
  
  v2lsxsc = datals$v2[c(1:9)]
  v2lsxsc[[4]]$resultraw = rbind(v2lsxsc[[4]]$resultraw, datals$v2[[10]]$resultraw)
  v2lsxsc[[5]]$resultraw = rbind(v2lsxsc[[5]]$resultraw, datals$v2[[11]]$resultraw)
  v2lsxsc[[6]]$resultraw = rbind(v2lsxsc[[6]]$resultraw, datals$v2[[12]]$resultraw)
  
  arylsxsc = datals$arraydata_bl[c(1:9)]
  arylsxsc[[4]] = rbind(datals$arraydata_bl[[4]], datals$arraydata_bl[[10]])
  arylsxsc[[5]] = rbind(datals$arraydata_bl[[5]], datals$arraydata_bl[[11]])
  arylsxsc[[6]] = rbind(datals$arraydata_bl[[6]], datals$arraydata_bl[[12]])
  
  ayocv1 = StatArrayOccSc(v1lsxsc)
  ayocv2 = StatArrayOccSc(v2lsxsc)
  ayoc = StatArrayOccScDual(arylsxsc)
  ayocv1$group = "v1"
  ayocv2$group = "v2"
  ayoc$group = "v1_v2"
  ayocv = rbind(ayocv1,ayocv2,ayoc)
  p1 = ggplot(ayocv,aes(x = Var1, y = norm, fill = group)) + 
    geom_bar(stat = "identity",position = "dodge",width = 0.7) +
    # scale_fill_manual(values=c('red','blue')) +
    scale_fill_jama()+
    theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))
  print(p1)
  ggsave(p1,filename = "final_result/CREST_final_analysis/E11_E15_basic_stat/barcode_occrurency_in_scsample_qc.pdf",width = 7,height = 5)
  write.table(ayocv,file = "final_result/CREST_final_analysis/E11_E15_basic_stat/barcode_occrurency_in__scsample_qc.txt",row.names = F,quote = F,
              sep = "\t")
  
}

#compare transcriptome&lineage
{
  ogne15cmp = list()
  #compare OGN & E15 correrlation
  {
    #combine OGN & E15 heatmap
    {
      ognln = datals$figure$heatmap$OGN
      e15ln = datals$figure$heatmap$E15
      
      ognlnmx = ognln$abspearson@matrix
      e15lnmx = e15ln$abspearson@matrix
      ggvenn(list("ogn" = colnames(ognlnmx),"e15" = colnames(e15lnmx)))
      cmcell = intersect(colnames(ognlnmx),colnames(e15lnmx))
      cmcellorde = orderls$E15[which(orderls$E15 %in% cmcell)]
      cmcellorde = orderls$OGN[which(orderls$OGN %in% cmcell)]
      e15lnmx = e15lnmx[cmcellorde,cmcellorde]
      ognlnmx = ognlnmx[cmcellorde,cmcellorde]
      
      col1 = colorRamp2(c(0, 0.1, 0.15, max(ognlnmx)), c("#08519C","#FEE5D9","#FC9F81", "#A50F15"))
      col2 = colorRamp2(c(0, 0.1, 0.2,max(e15lnmx)), c("#08519C","#FEE5D9","#FC9F81", "#A50F15"))
      
      pdf("final_result/CREST_final_analysis/rebbutal_analysis/ogne15cmpheatmap_ogn.pdf",
          width = 8,height = 6)
      ognht = Heatmap(ognlnmx,heatmap_legend_param = list(title = "abpearson"),
                   col = col1,
                   cluster_rows =  F, cluster_columns = F)
      ognht
      dev.off()
      
      pdf("final_result/CREST_final_analysis/rebbutal_analysis/ogne15cmpheatmap_e15.pdf",
          width = 8,height = 6)
      e15ht = Heatmap(e15lnmx,heatmap_legend_param = list(title = "abpearson"),
                      col = col2,
                      cluster_rows =  F, cluster_columns = F)
      e15ht
      dev.off()
      
      hce15 = hclust(as.dist(1-e15lnmx))
      hcogn = hclust(as.dist(1-ognlnmx))
      orderls$E15
      #heatmap combine
      mxt = ognlnmx
      for (i in 1:nrow(e15lnmx)) {
        mxt[i,i] = NA
        for (j in 1:(i-1)) {
          mxt[i,j] = e15lnmx[i,j]
        }
      }
      mxt[is.na(mxt)] = 0
      mxt[1,1] = 0
      library(circlize)
      col1 = colorRamp2(c(0, 0.1, max(mxt)), c("white", "#FC9F81", "#A50F15"))
      col1 = colorRamp2(c(0, 0.1, max(mxt)), c("#08519C","white", "#A50F15"))
      # col1 = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4))
      pdf("final_result/CREST_final_analysis/rebbutal_analysis/E15_OGN_correlation_cmp_heatmap2_1_9.pdf",
          width = 8, height = 6)
      ht = Heatmap(mxt,heatmap_legend_param = list(title = "abpearson"),
                   col = col1,
                   # col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
                   # border = T, rect_gp = gpar(col = "black"),
                   layer_fun = function(j, i, x, y, w, h, fill) {
                     l = i > j
                     grid.rect(x[l], y[l], w[l], h[l], 
                               gp = gpar(fill = col1(pindex(mxt, i[l], j[l])),col = "black",lwd = 1))
                     l = i < j
                     grid.rect(x[l], y[l], w[l], h[l], 
                               gp = gpar(fill = col1(pindex(mxt, i[l], j[l])),col = "black",lwd = 1))
                   },
                   # show_heatmap_legend = F,
                   cluster_rows =  F, cluster_columns = F)
      ogne15cmp$combine_heatmap_ognlow_e15up = ht
      lgd = Legend(labels = c("OGN(lowertri)", "snapE15.5(upertri)"), title = "border",
                   legend_gp = gpar(fill = c("black", "black")))
      draw(ht, heatmap_legend_list = list(lgd))
      dev.off()
    }
    ognaray = datals$arraydata_bl$OGN
    E15aray = datals$arraydata_bl$E15
    cmcell = intersect(ognaray$Cell.type,E15aray$Cell.type)
    ognaray = ognaray[which(ognaray$Cell.type %in% cmcell),]
    E15aray = E15aray[which(E15aray$Cell.type %in% cmcell),]
    
    MyHeatmapCal2 = function(arraydata,title,order = NULL){
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
                     col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(8)))
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
                     col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(8)),
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
    ogncmpht = MyHeatmapCal2(ognaray,"OGN",orderls$E15)
    e15cmpht = MyHeatmapCal2(E15aray,"E15",orderls$E15)
    e15cmpht$abspearson_withpv
    ogncmpht$abspearson_withpv
    
    #combine OGN & E15 correlation lineage tree
    {
      library(dendextend)
      dend1 <- as.dendrogram (hce15)
      dend2 <- as.dendrogram (hcogn)
      
      pdf("final_result/CREST_final_analysis/rebbutal_analysis/hclust_tree_cmp_1_9.pdf",width = 8,height = 6)
      # plot(hce15,cex = 1,hang=-1,main = "E15")
      # plot(hcogn,cex = 1,hang=-1,main = "OGN")
      lntreecmp = dendlist(dend1, dend2) %>%
        untangle(method = "step1side") %>% # Find the best alignment layout
        tanglegram() 
      ogne15cmp$lineagetreecmp = lntreecmp
      dev.off()
    }
    
    
    e15lnmxl = datals$figure$heatmap_order$E15$cormx
    ognlnmxl = datals$figure$heatmap_order$OGN$cormx
    
    cmplnmxl = merge(e15lnmxl,ognlnmxl,by = c('Var1','Var2'),suffixes = c('.e15','.ogn'))
    cmplnmxl = cmplnmxl[which(cmplnmxl$Var1 != cmplnmxl$Var2),]
    cmplnmxl = cmplnmxl[!duplicated(cmplnmxl[c(3,5)]),]
    
    cmplnmxl$pvalue.e15 = p.adjust(cmplnmxl$pvalue.e15)
    cmplnmxl$pvalue.ogn = p.adjust(cmplnmxl$pvalue.ogn)
    cmplnmxl = cmplnmxl[which(cmplnmxl$pvalue.e15 < 0.05 | cmplnmxl$pvalue.ogn < 0.05),]
    # cmplnmxl$corr.e15 = log2(cmplnmxl$corr.e15);cmplnmxl$corr.ogn = log2(cmplnmxl$corr.ogn)
    p1.0 = ggplot(cmplnmxl,aes(x = log2(corr.e15+0.00001), y = log2(corr.ogn+0.00001))) + geom_point() + 
      annotate("text",x = -10, y = 0,
               label = paste0("r = ", round(cor(cmplnmxl$corr.e15,cmplnmxl$corr.ogn),3),
                              ",t test p = ", round(t.test(cmplnmxl$corr.e15,cmplnmxl$corr.ogn,paired = T)$p.value,3))) +
      # ggrepel::geom_text_repel( data = cmplnmxl[which(cmplnmxl$corr.e15 > 0.1),], color = "red",
      #                           aes(label = paste0(Var1,"-",Var2))) +
      xlab("E15 lineage correlation(log2) of cell clusters") + ylab("OGN lineage correlation(log2) of cell clusters") +
      theme_pubr()
    p1.0
    ggsave(p1.0,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15_OGN_correlation_cmp_1_9.pdf",
           width = 6,height = 6)
    
    
    #PCA analysis
    e15ardf = datals$arraydata_bl$E15
    ognardf = datals$arraydata_bl$OGN
    e15clone = dcast(e15ardf[which(e15ardf$Cell.type %in% cmcell),],pattern~Cell.type)
    ognclone = dcast(ognardf[which(ognardf$Cell.type %in% cmcell),],pattern~Cell.type)
    e15clone$sample = "E15"
    ognclone$sample = "OGN"
    cmpclone = rbind(e15clone,ognclone)
    rownames(cmpclone) = cmpclone$pattern;cmpclone = cmpclone[-1]
    cmpclone = cmpclone[rowSums(cmpclone[-ncol(cmpclone)]) > 1,]
    pc <- prcomp(cmpclone[-ncol(cmpclone)],
                 center = TRUE,
                 scale. = TRUE)
    attributes(pc)
    library(ggfortify)
    library(cluster)
    p1.1.1 = autoplot(pc, data = cmpclone, colour = 'sample',size = 0.1,alpha = 0.5) + 
      scale_color_nejm()+
      theme_pubr() + theme(legend.position = "bottom")
    pcdf = as.data.frame(pc$x)
    pcdf$sample = cmpclone[rownames(pcdf),"sample"]
    
    library(patchwork)
    library(ggpubr)
    p1.1.3 = ggplot(pcdf,aes(x = PC1, color = sample,fill = sample)) + 
      # geom_boxplot(width = 0.5,color = "black") + 
      geom_density(alpha = 0.5) + 
      theme_pubr() +
      scale_fill_nejm()+
      # coord_flip()+
      xlab("") + ylab("") +
      # stat_compare_means(comparisons = list(c("E15","OGN")), 
      #                    label = "p.signif",
      #                    method = "t.test") +
      theme(legend.position = "none",axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # t.test(pcdf[which(pcdf$sample == "E15"),"PC1"],pcdf[which(pcdf$sample == "OGN"),"PC1"],paired = F)
    p1.1.3
    p1.1.4 = ggplot(pcdf,aes(y = PC2, color = sample,fill = sample)) + 
      # geom_boxplot(width = 0.5,color = "black") + 
      theme_pubr() +
      geom_density(alpha = 0.5) + 
      scale_fill_nejm()+
      # stat_compare_means(comparisons = list(c("E15","OGN")), 
      #                    label = "p.signif",
      #                    method = "t.test") +
      xlab("") + ylab("") +
      theme(legend.position = "none",axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    p1.1.4
    p1.1.2 = p1.1.3 + plot_spacer() + p1.1.1 + p1.1.4 + 
      plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    p1.1.2
    ggsave(p1.1.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/E15_OGN_clone_pca_1_9.pdf",
           width = 6,height = 6)
    
    ogne15cmp$lineagepointcmp = p1.0
    ogne15cmp$lineagepcacmp = p1.1.2
    
  }
  
  #compare cell composistion(ref vs snapE15.5)
  {
    cellogn = cellls$OGN
    celle15 = cellls$E15
    cellstogn = as.data.frame(table(cellogn$Cell.type))
    cellste15 = as.data.frame(table(celle15[which(celle15$Cell.type != "En1/Ai9-E15.5-ref"),]$Cell.type))
    cellst = merge(cellste15,cellstogn,by = "Var1")
    cellst$prop.e15 = cellst$Freq.x/sum(cellst$Freq.x)
    cellst$prop.ogn = cellst$Freq.y/sum(cellst$Freq.y)
    cor(cellst$prop.e15,cellst$prop.ogn)
    
    
    celle15$ref = "ref"
    celle15[which(celle15$group != "En1/Ai9-E15.5-ref"), "ref"] = "treat"
    cellcmp = celle15 %>% group_by(ref,ident) %>% summarise(cellnum = n())
    cellcmp = cellcmp %>% group_by(ref) %>% mutate(cellprop = cellnum/sum(cellnum))
    cellcmp = dcast(cellcmp,ident~ref,value.var = "cellprop")
    cor(cellcmp$ref,cellcmp$treat)
    
    cellcmpl = melt(cellcmp[1:3])
    tmp = t.test(cellcmp$ref,cellcmp$treat,paired = T)
    tmp$p.value
    tmp$statistic
    
    p1.2.1 = ggplot(cellcmp,aes(x = ref,y = treat)) + geom_point() +
      # scale_color_manual(values = c("black","red")) +
      # ggrepel::geom_text_repel(data = cellcmp[which(cellcmp$highlight == "Y"),],aes(label = ident)) +
      xlim(c(0,0.17)) + ylim(c(0,0.17)) +
      annotate("text",x = 0.03, y = 0.17,label = paste0("r = ", round(cor(cellcmp$ref,cellcmp$treat),3),
                                                        ", p = ",tmp$p.value)) +
      xlab("En1/Ai9-E15.5-ref different cell proportion") + ylab("snapE15.5 different cell proportion") +
      theme_pubr() + theme(legend.position =  "none")
    p1.2.1
    
    cellcmpl = cellcmpl %>% group_by(ident) %>% mutate(meanprop = mean(value))
    cellcmpl = cellcmpl[order(cellcmpl$meanprop),]
    cellcmpl$ident = factor(cellcmpl$ident,levels = unique(as.character(cellcmpl$ident)))
    p1.2.2 = ggplot(cellcmpl,aes(x = value,y = ident,fill = variable)) + 
      geom_line(aes(group = ident))+
      geom_point(shape = 21,size = 3, colour = "black") + 
      annotate("text",x = 0.12, y = 3 ,label = paste0("r = ", round(cor(cellcmp$ref,cellcmp$treat),3),
                                                      ", p = ",tmp$p.value)) +
      # scale_fill_manual(label = c("En1/Ai9-E15.5-ref","snapE15.5")) +
      scale_fill_nejm(label = c("En1/Ai9-E15.5-ref","snapE15.5"))+
      xlab("cell proportion") + ylab("") +
      theme_bw()
    p1.2.2
    ggsave(p1.2.1,filename = "final_result/CREST_final_analysis/rebbutal_analysis/ref_E15_cellcmp_point.pdf",
           width = 5,height = 4)
    ggsave(p1.2.2,filename = "final_result/CREST_final_analysis/rebbutal_analysis/ref_E15_cellcmp_dumbell.pdf",
           width = 6,height = 6)
    E15_ref_cellcmp = list("point" = p1.2.1,"dumbell" = p1.2.2)
    
  }
  
  
  #compare cell composition and transcriptome
  {
    ptmp = DimPlot(E15trans, group.by='orig.ident',pt.size=0.05,raster = F) + 
      scale_color_manual(values = c("grey","#219ebc","#126782","#023047","#ffb703","#fd9e02","#fb8500"))
    ptmp[[1]]$layers[[1]]$aes_params$alpha =  .5
    ptmp[[1]]
    ggsave(ptmp,filename = "final_result/CREST_final_analysis/annotation_20221207/E15.5_UMAP_orig.pdf",
           width = 8,height = 6)
    #transcriptome compare
    E15trans = qread("final_result/CREST_final_analysis/annotation_20221207/E15_subclustered_renamed_20221207.qs")
    snapE15 = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                        c("snapE15.5-rep1","snapE15.5-rep2","snapE15.5-rep3"))]))
    refE15 = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                       c("En1/Ai9-E15.5-ref"))]))
    snapdf = AverageExpression(snapE15)
    refdf = AverageExpression(refE15)
    # snapdfl = melt(snapdf$RNA)
    snapdfl = melt(snapdf$integrated)
    # refdfl = melt(refdf$RNA)
    refdfl = melt(refdf$integrated)
    cmpdf = merge(snapdfl,refdfl,by = c('Var1','Var2'))
    colnames(cmpdf) = c("gene","Cell.type","snapE15.5","refE15.5")
    cmpdf$Cell.type = as.character(cmpdf$Cell.type)
    cellid = unique(cmpdf$Cell.type)
    p1.3 = list()
    for (i in 1:length(cellid)) {
      cmpdfi = cmpdf[which(cmpdf$Cell.type == cellid[i]),]
      # tmp = cor.test(cmpdfi$snapE15.5, cmpdfi$refE15.5)
      # tmp = t.test(cmpdfi$snapE15.5, cmpdfi$refE15.5, paired = F)
      p1.3.1 = ggplot(cmpdfi,aes(x = log2(snapE15.5+1),y = log2(refE15.5+1))) + geom_point(size = 0.5) + 
        annotate("text",x = log2(max(cmpdfi$snapE15.5)+1)/5, y = log2(max(cmpdfi$refE15.5)+1),
                 label = paste0("r = ", round(cor(cmpdfi$snapE15.5,cmpdfi$refE15.5),3)),
                                # ", p = ",round(tmp$p.value,3)),
                 size = 2) +
        geom_smooth()+
        theme_classic() + xlab(paste0("log2 snap ",cellid[i])) + ylab(paste0("log2 ref ",cellid[i]))
      p1.3[[i]] = p1.3.1
      names(p1.3)[i] = cellid[i]
    }
    library(ggExtra)
    library(gridExtra)
    pdf("final_result/CREST_final_analysis/rebbutal_analysis/ref_E15_transcmp2.pdf",width = 10,height = 8)
    do.call("grid.arrange", c(p1.3, ncol=6))
    dev.off()
    
    
    E15trans = qread("final_result/CREST_final_analysis/annotation_20221207/E15_Organoid_allintegrated_20230108.qs")
    snapE15 = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                        c("snapE15.5-rep1","snapE15.5-rep2","snapE15.5-rep3"))]))
    pandaOGN = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                         c("pandaOGN-rep1","pandaOGN-rep2","pandaOGN-rep3"))]))
    Idents(pandaOGN) = pandaOGN$ClusterName
    Idents(snapE15) = snapE15$ClusterName
    DimPlot(pandaOGN,group.by = "ClusterName")
    DimPlot(snapE15,group.by = "ClusterName")
    snapdf = AverageExpression(snapE15)
    ogndf = AverageExpression(pandaOGN)
    # snapdfl = melt(snapdf$RNA)
    snapdfl = melt(snapdf$integrated)
    # refdfl = melt(refdf$RNA)
    ogndfl = melt(ogndf$integrated)
    cmpdf = merge(snapdfl,ogndfl,by = c('Var1','Var2'))
    colnames(cmpdf) = c("gene","Cell.type","snapE15.5","pandaOGN")
    cmpdf$Cell.type = as.character(cmpdf$Cell.type)
    cellid = unique(cmpdf$Cell.type)
    p1.4 = list()
    for (i in 1:length(cellid)) {
      cmpdfi = cmpdf[which(cmpdf$Cell.type == cellid[i]),]
      p1.4.1 = ggplot(cmpdfi,aes(x = log2(snapE15.5+1),y = log2(pandaOGN+1))) + 
        geom_point(size = 0.5) + 
        geom_smooth() +
        annotate("text",x = max(log2(cmpdfi$snapE15.5+1))/5, y = max(log2(cmpdfi$pandaOGN+1)),
                 label = paste0("r = ", round(cor(cmpdfi$snapE15.5,cmpdfi$pandaOGN),3)),
                 size = 2) +
        theme_classic() + xlab(paste0("log2 snap ",cellid[i])) + ylab(paste0("log2 ogn ",cellid[i]))
      p1.4[[i]] = p1.4.1
      names(p1.4)[i] = cellid[i]
    }
    library(ggExtra)
    library(gridExtra)
    pdf("final_result/CREST_final_analysis/rebbutal_analysis/OGN_E15_transcmp2_1_9.pdf",width = 12,height = 6)
    do.call("grid.arrange", c(p1.4, nrow=3))
    dev.off()
    transcmp = list("E15vsRef" = p1.3,"E15vsOGN" = p1.4)
    
  }
  
  #
  ognrevisels$transcmp = p1.4
}

#NbMP Rgl1 DA deg
{
  E15trans = qread("final_result/CREST_final_analysis/annotation_20221207/E15_Organoid_neu_renamed_integrated_20221211.qs")
  snapE15 = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                      c("snapE15.5-rep1","snapE15.5-rep2","snapE15.5-rep3"))]))
  pandaOGN = subset(E15trans,cells = names(E15trans$orig.ident[which(E15trans$orig.ident %in% 
                                                                       c("pandaOGN-rep1","pandaOGN-rep2","pandaOGN-rep3"))]))
  E11trans = qread("final_result/CREST_final_analysis/annotation_20221207/E11_subclustered_renamed_UMAPswap_20221209.qs")
  
  library(Seurat)
  #NbMP DA & Rgl1 DA
  ognda = qread("final_result/CREST_final_analysis/annotation_20221207/panda_DA_cell.BC_20230109.qs")
  pandaOGN = RenameCells(pandaOGN,new.names = substr(colnames(pandaOGN),1,20))
  npda = subset(pandaOGN,cells = ognda$NP_DA)
  rglda = subset(pandaOGN,cells = ognda$Rgl1_DA)
  npdadf = AverageExpression(npda)
  rgldadf = AverageExpression(rglda)
  npdadfl = melt(npdadf$integrated)
  rgldadfl = melt(rgldadf$integrated)
  cmpdf = merge(npdadfl,rgldadfl,by = c('Var1','Var2'))
  colnames(cmpdf) = c("gene","Cell.type","NPBM","Rgl1")
  p1.5 = ggplot(cmpdf,aes(x = NPBM,y = Rgl1)) + geom_point(size = 0.5) + 
    annotate("text",x = max(cmpdf$NPBM)/5, y = max(cmpdf$Rgl1),
             label = paste0("r = ", round(cor(cmpdf$NPBM,cmpdf$Rgl1),3),
                            ", t test p = ", round(t.test(cmpdf$NPBM,cmpdf$Rgl1,paired = T)$p.value,3)
                            ),
             size = 3) +
    ggtitle("Gene expression correlation") +
    theme_classic() + xlab("NPBM Derived DA") + ylab("Rgl1 Derived DA")
  p1.5
  ggsave(p1.5,filename = "final_result/CREST_final_analysis/rebbutal_analysis/OGN_NPBM-Rgl1_DA_transcmp_point_1_9.pdf",
         width = 5,height = 4)
  
  
  #compare NPBM DA/ none DA
  colnames(E11trans)
  
  npbmda = qread("final_result/CREST_final_analysis/annotation_20221207/panda_NPBM_DA_fate.BC_20230109.qs")
  # npbmda$NProgBM_DA = paste0(substr(npbmda$NProgBM_DA,1,19), as.numeric(substr(npbmda$NProgBM_DA,20,20))-3)
  # npbmda$NProgBM_others = paste0(substr(npbmda$NProgBM_others,1,19), as.numeric(substr(npbmda$NProgBM_others,20,20))-3)
  npda = subset(E11trans,cells = npbmda$NProgBM_DA)
  otherda = subset(E11trans,cells = npbmda$NProgBM_others)
  npdadf = AverageExpression(npda)
  otherdadf = AverageExpression(otherda)
  npdadfl = melt(npdadf$integrated)
  otherdadfl = melt(otherdadf$integrated)
  cmpdf = merge(npdadfl,otherdadfl,by = c('Var1','Var2'))
  colnames(cmpdf) = c("gene","Cell.type","DA_NPBM","Other")
  p1.6 = ggplot(cmpdf,aes(x = DA_NPBM,y = Other)) + geom_point(size = 0.5) + 
    annotate("text",x = max(cmpdf$DA_NPBM)/5, y = max(cmpdf$Other),
             label = paste0("r = ", round(cor(cmpdf$DA_NPBM,cmpdf$Other),3),
                            ", t test p = ", round(t.test(cmpdf$DA_NPBM,cmpdf$Other,paired = T)$p.value,3)),
             size = 3) +
    ggtitle("Gene expression correlation") +
    theme_classic() + xlab("NPBM DA potential") + ylab("Other NPBM")
  p1.6
  ggsave(p1.6,filename = "final_result/CREST_final_analysis/rebbutal_analysis/OGN_NPBM-DA_NPBM-other_transcmp_point_1_9.pdf",
         width = 5,height = 4)
  
  #NPBM none DA DEG
  unique(E11trans$orig.ident)
  tmp = E11trans$orig.ident
  tmp = tmp[which(tmp == "pandaE11.5-rep3")]
  
  npdat = subset(E11trans,cells = c(npbmda$NProgBM_DA,npbmda$NProgBM_others))
  newgroup = c(rep("DA-fate",length(npbmda$NProgBM_DA)),rep("noDA-fate",length(npbmda$NProgBM_others)))
  names(newgroup) = c(npbmda$NProgBM_DA,npbmda$NProgBM_others)
  Idents(npdat) = newgroup
  DimPlot(npdat)
  npdadeg = FindMarkers(npdat,ident.1 = "DA-fate",ident.2 = "noDA-fate",min.pct = 0.1,logfc.threshold = 0.1)
  library(ggplot2)
  pt = VocanoFigure(npdadeg,"DA fate vs non-DA fate NPBM",thres = 0.4)
  ggsave(pt,filename = "final_result/CREST_final_analysis/OGN_analysis/NPBMDA_noDA_DEG_1_9.pdf",width = 7,height = 5)
  
  ognrevisels$DAgenecmp = list("NbMP_Rgl1_DA_cmp" = p1.5, "NPBM_nonepro_DA_cmp" = p1.6,"NPBM_nonepro_DA_DEG" = pt)
  qsave(ognrevisels,file = "final_result/CREST_final_analysis/ognrevisels.qs")
  
  
}


datalsold = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_11_edit.qs")
datals = qread("final_result/CREST_final_analysis/array_recover_data_total_withfigure_12_16.qs")
names(datalsold$figure)
names(datals$figure)
datals$figure$pandaE11_OGNls_f6 = datalsold$figure$pandaE11_OGNls_f6
datals$figure$sister_e11_cmp_fs9 = datalsold$figure$sister_e11_cmp_fs9
datals$figure$bulkstatls_fs2 = datalsold$figure$bulkstatls_fs2

