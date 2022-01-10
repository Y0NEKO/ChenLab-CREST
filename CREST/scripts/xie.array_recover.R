#handle array recover problem of xie's data

setwd("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/")

source("myscript/anascipt/tree_method.R")

#Dual AT E15.5-----------------------
#read data
{
  samplels = list.files("sc_result/dual_E15/",full.names = T)[c(1,4)]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
  names(v1lsx) = list.files("sc_result/dual_E15/")[c(1,4)]
  names(v2lsx) = list.files("sc_result/dual_E15/")[c(1,4)]
  cell = ReadCell("sc_result/dual_E15/celltype_subclustered0901.csv")
  samplels = list.files("sc_result/dual_E15/",full.names = T)[c(1,4)]
}

#filter data test
{
  
  AT = v1lsx$AT
  ATDu = v1lsx$AT_Dual1_2
  Pha = v1lsx$Pha
  
  ATDu$Cell.BC
  
  # comBC = intersect(AT$Cell.BC,Pha$Cell.BC)
  
  ATDuf = ATDu[which(ATDu$umi_pro>0  & ATDu$umi_num >10),]
  ATf = AT[which(AT$reads_pro.umim>0  & AT$umi_num >10),]
  Phaf = Pha[which(Pha$reads_pro.umim>0  & Pha$umi_num >10),]
  
  comBC = Reduce(intersect,list(ATf$Cell.BC,ATDuf$Cell.BC,Phaf$Cell.BC))
  ATDuf = ATDu[which(ATDu$Cell.BC %in% comBC),]
  ATf = AT[which(AT$Cell.BC %in% comBC),]
  Phaf = Pha[which(Pha$Cell.BC %in% comBC),]
  
  ggvenn::ggvenn(list("AT" = paste0(ATf$Cell.BC,"-",ATf$umim),
                      "AT_Dual1_2" = paste0(ATDuf$Cell.BC,"-",ATDuf$umim),
                      "Pha" = paste0(Phaf$Cell.BC,"-",Phaf$umim)),
                 stroke_alpha = 0.5)
  
  comarray = Reduce(intersect,list("AT" = paste0(ATf$Cell.BC,"-",ATf$umim),
                                   "AT_Dual1_2" = paste0(ATDuf$Cell.BC,"-",ATDuf$umim),
                                   "Pha" = paste0(Phaf$Cell.BC,"-",Phaf$umim)))
  
  arstat1 = NULL
  for (i in 1:20) {
    ATDuf = v1lsx$AT_Dual1_2[which(v1lsx$AT_Dual1_2$reads_num >i),]
    ATf = v1lsx$AT[which(v1lsx$AT$reads_num >i),]
    Phaf = v1lsx$Pha[which(v1lsx$Pha$reads_num >i),]
    
    comBC = Reduce(intersect,list(ATf$Cell.BC,ATDuf$Cell.BC,Phaf$Cell.BC))
    ATDuf = ATDuf[which(ATDuf$Cell.BC %in% comBC),]
    ATf = ATf[which(ATf$Cell.BC %in% comBC),]
    Phaf = Phaf[which(Phaf$Cell.BC %in% comBC),]
    comarray = Reduce(intersect,list("AT" = paste0(ATf$Cell.BC,"-",ATf$umim),
                                     "AT_Dual1_2" = paste0(ATDuf$Cell.BC,"-",ATDuf$umim),
                                     "Pha" = paste0(Phaf$Cell.BC,"-",Phaf$umim)))
    arstat1 = rbind(arstat1,
                    data.frame("reads_num" = i,"total" = length(comBC),"intersect" = length(comarray)))
    
  }
  arstat2 = NULL
  for (i in 1:20) {
    ATDuf = v2lsx$AT_Dual1_2[which(v2lsx$AT_Dual1_2$reads_num >i),]
    ATf = v2lsx$AT[which(v2lsx$AT$reads_num >i),]
    Phaf = v2lsx$Pha[which(v2lsx$Pha$reads_num >i),]
    
    comBC = Reduce(intersect,list(ATf$Cell.BC,ATDuf$Cell.BC,Phaf$Cell.BC))
    ATDuf = ATDuf[which(ATDuf$Cell.BC %in% comBC),]
    ATf = ATf[which(ATf$Cell.BC %in% comBC),]
    Phaf = Phaf[which(Phaf$Cell.BC %in% comBC),]
    comarray = Reduce(intersect,list("AT" = paste0(ATf$Cell.BC,"-",ATf$umim),
                                     "AT_Dual1_2" = paste0(ATDuf$Cell.BC,"-",ATDuf$umim),
                                     "Pha" = paste0(Phaf$Cell.BC,"-",Phaf$umim)))
    arstat2 = rbind(arstat2,
                    data.frame("reads_num" = i,"total" = length(comBC),"intersect" = length(comarray)))
    
  }
  # arstat2 = arstat
  arstat1$overpop = round(arstat1$intersect/arstat1$total,2)
  arstat2$overpop = round(arstat2$intersect/arstat2$total,2)
  
  p1.1 = ggplot(arstat1)  + 
    geom_bar(aes(x=reads_num, y=total),stat="identity")+
    geom_line(aes(x=reads_num, y=intersect),stat="identity")+
    geom_text(aes(label=total, x=reads_num, y=total), colour="black",size = 3)+
    geom_text(aes(label=overpop, x=reads_num, y=intersect*0.8), colour="black",size = 3) +
    scale_fill_viridis() + labs(title="V1") +
    theme_bw(base_family = "sans")
  p1.1
  
  p1.2 = ggplot(arstat2)  + 
    geom_bar(aes(x=reads_num, y=total),stat="identity")+
    geom_line(aes(x=reads_num, y=intersect),stat="identity")+
    geom_text(aes(label=total, x=reads_num, y=total), colour="black",size = 3)+
    geom_text(aes(label=overpop, x=reads_num, y=intersect*0.8), colour="black",size = 3) +
    scale_fill_viridis() +labs(title="V2") +
    theme_bw(base_family = "sans")
  p1.2
  
  
  ggsave(p1.1+p1.2,filename = "final_result/xDual/v1_v2_recover_test.pdf",width = 15,height = 4)
  
}

#build tree
{
  prefix = c("v1.","v2.")
  
  
  #测试这样过滤得到的细胞建树情况
  # ATDu1 = v1lsx$AT_Dual1_2[which(v1lsx$AT_Dual1_2$reads_num >10),]
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
    
    p1.0 = ggvenn(rawlist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[1])
    p1.1 = ggvenn(rflist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[2])
    p1.2 = ggvenn(baflist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[3])
    p = ggarrange(p1.0,p1.1,p1.2,nrow = 3)
    
    return(list("result" = result,"rawlist" = rawlist,"rflist" = rflist,"baflist" = baflist,"stat" = stat,"venn_plot" = p))
  }
  v1title = c("E15.5 V1 raw BC recover rate","E15.5 V1 reads num(10) filtered BC recover rate",
              "E15.5 V1 reads num(10) filtered BC-pattern recover rate")
  v1test = FilterRecover(v1lsx,cell,10,v1title)
  v1test$stat
  v1test$venn_plot
  
  v2title = c("E15.5 V2 raw BC recover rate","E15.5 V2 reads num(10) filtered BC recover rate",
              "E15.5 V2 reads num(10) filtered BC-pattern recover rate")
  v2test = FilterRecover(v2lsx,cell,10,v2title)
  v2test$stat
  v2test$venn_plot
  
  ggsave(ggarrange(v1test$venn_plot,v2test$venn_plot,nrow = 1),
         filename = "final_result/xE15.5/E15.5_recover_rate_stat.pdf",
         width = 14,height = 10)
  #
  
  # bl1 = read.csv("sc_result/dual_E15/v1_indel_blacklist.csv",header = T,row.names = 1)
  # colnames(bl1) = "tag"
  # bl2 = read.csv("sc_result/dual_E15/v2_indel_blacklist.csv",header = T,row.names = 1)
  # colnames(bl2) = "tag"
  
  
  
  bl = readRDS("final_result/x.blacklist.rds")
  v1testf = BlacklistFilter(v1test$result,bl$v1_array)
  v2testf = BlacklistFilter(v2test$result,bl$v2_array)
  dataE15 = list(v1testf,v2testf)
  outnameE15 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_blay"
  E15tree = TreeConstruct(dataE15,prefix,outnameE15,cell)
  E15tree = CirclePlot1(E15tree,outnameE15)
  E15tree = MyHeatmap(E15tree,outnameE15)
  saveRDS(E15tree,file = "final_result/xE15.5/E15.5_blay_tree.rds")
  
  v1testf = BlacklistFilter(v1test$result,bl$v1_tag)
  v2testf = BlacklistFilter(v2test$result,bl$v2_tag)
  dataE15 = list(v1testf,v2testf)
  outnameE15 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_bltg"
  E15tree = TreeConstruct(dataE15,prefix,outnameE15,cell)
  E15tree = CirclePlot1(E15tree,outnameE15)
  E15tree = MyHeatmap(E15tree,outnameE15)
  saveRDS(E15tree,file = "final_result/xE15.5/E15.5_bltg_tree.rds")
  
  saveRDS(v1test,file = "final_result/xE15.5/E15.5_v1_result.rds")
  saveRDS(v2test,file = "final_result/xE15.5/E15.5_v2_result.rds")
}



#E15.5 2------------------
{
  samplels = list.files("sc_result/x.TAQ_E15/",full.names = T)[c(1,4,5)]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
  names(v1lsx) = list.files("sc_result/x.TAQ_E15/")[c(1,4,5)]
  names(v2lsx) = list.files("sc_result/x.TAQ_E15/")[c(1,4,5)]
  cell = ReadCell("sc_result/x.TAQ_E15/E15clustername20211126 v2.csv")
  # cellt = cell
  # cell[which(cellt$Cell.type == "Rgl1"),"Cell.type"] = "Rgl3"
  # cell[which(cellt$Cell.type == "Rgl3"),"Cell.type"] = "Rgl1"
  # write.csv(cell,file = "sc_result/dual_E15/celltype_subclustered1117.csv",row.names = F,quote = F)
}
#build tree
{
  prefix = c("v1.","v2.")
  
  
  #测试这样过滤得到的细胞建树情况
  # ATDu1 = v1lsx$AT_Dual1_2[which(v1lsx$AT_Dual1_2$reads_num >10),]
  FilterRecover = function(v1lsx,cell,thres,title){
    
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
    
    p1.0 = ggvenn(rawlist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[1])
    p1.1 = ggvenn(rflist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[2])
    p1.2 = ggvenn(baflist,stroke_alpha = 0.5,stroke_size = 0.5) + ggtitle(title[3])
    p = ggarrange(p1.0,p1.1,p1.2,nrow = 3)
    
    return(list("result" = result,"rawlist" = rawlist,"rflist" = rflist,"baflist" = baflist,"stat" = stat,"venn_plot" = p))
  }
  v1title = c("E15.5 V1 raw BC recover rate","E15.5 V1 reads num(10) filtered BC recover rate",
              "E15.5 V1 reads num(10) filtered BC-pattern recover rate")
  v1test = FilterRecover(v1lsx,cell,10,v1title)
  v1test$stat
  v1test$venn_plot
  
  v2title = c("E15.5 V2 raw BC recover rate","E15.5 V2 reads num(10) filtered BC recover rate",
              "E15.5 V2 reads num(10) filtered BC-pattern recover rate")
  v2test = FilterRecover(v2lsx,cell,10,v2title)
  v2test$stat
  v2test$venn_plot
  
  # v1test = readRDS("final_result/xE15.5/e15_v1qc.rds")
  saveRDS(v1test,file = "final_result/xE15.5/e15_v1qc.rds")
  # v2test = readRDS("final_result/xE15.5/e15_v2qc.rds")
  saveRDS(v2test,file = "final_result/xE15.5/e15_v2qc.rds")
  
  
  ggsave(ggarrange(v1test$venn_plot,v2test$venn_plot,nrow = 1),
         filename = "final_result/xE15.5/E15.5_tap_recover_rate_stat.pdf",
         width = 14,height = 10)
  
  bl = readRDS("final_result/x.blacklist.rds")
  v1testf = BlacklistFilter(v1test$result,bl$v1_array)
  v2testf = BlacklistFilter(v2test$result,bl$v2_array)
  dataE15 = list(v1testf,v2testf)
  outnameE15 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_blay_taq"
  E15tree = TreeConstruct(dataE15,prefix,outnameE15,cell)
  E15tree = CirclePlot1(E15tree,outnameE15)
  E15tree = MyHeatmap(E15tree,outnameE15)
  saveRDS(E15tree,file = "final_result/xE15.5/E15.5_blay_taq_tree.rds")
  
  v1testftg = BlacklistFilter(v1test$result,bl$v1_tag)
  v2testftg = BlacklistFilter(v2test$result,bl$v2_tag)
  dataE15tg = list(v1testftg,v2testftg)
  outnameE15tg = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_bltg_taq"
  E15treetg = TreeConstruct(dataE15tg,prefix,outnameE15tg,cell)
  E15treetg = CirclePlot1(E15treetg,outnameE15tg)
  E15treetg = MyHeatmap(E15treetg,outnameE15tg)
  saveRDS(E15treetg,file = "final_result/xE15.5/E15.5_bltg_taq_tree.rds")
  
  # saveRDS(E15tree,file = "final_result/xE15.5/E15.5_taq_tree.rds")
  # saveRDS(v1test,file = "final_result/xE15.5/E15.5_taq_v1_result.rds")
  # saveRDS(v2test,file = "final_result/xE15.5/E15.5_taq_v2_result.rds")
  
  #单独V1, V2
  bl = readRDS("final_result/x.blacklist.rds")
  v1testf = BlacklistFilter(v1test$result,bl$v1_array)
  v2testf = BlacklistFilter(v2test$result,bl$v2_array)
  
  prefix = c("v1.")
  dataE15v1 = list(v1testf)
  outnameE15v1 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_blay_taq_v1"
  E15treev1 = TreeConstruct(dataE15v1,prefix,outnameE15v1,cell)
  E15treev1 = CirclePlot1(E15treev1,outnameE15v1)
  E15treev1 = MyHeatmap(E15treev1,outnameE15v1)
  saveRDS(E15treev1,file = "final_result/xE15.5/E15.5_blay_taq_v1_tree.rds")
  
  prefix = c("v2.")
  dataE15v2 = list(v2testf)
  outnameE15v2 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE15.5/E15.5_blay_taq_v2"
  E15treev2 = TreeConstruct(dataE15v2,prefix,outnameE15v2,cell)
  E15treev2 = CirclePlot1(E15treev2,outnameE15v2)
  E15treev2 = MyHeatmap(E15treev2,outnameE15v2)
  saveRDS(E15treev2,file = "final_result/xE15.5/E15.5_blay_taq_v2_tree.rds")
  
}



#E11-----------------------------------
#Read data
{
  samplels = list.files("sc_result/X.E11",full.names = T)[c(5:7)]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = list.files("sc_result/X.E11/")[c(5:7)]
  names(v2lsx) = list.files("sc_result/X.E11/")[c(5:7)]
  
  cell11 = ReadCell("sc_result/X.E11/E11clustername20211115.csv")
  
}

#Build Tree
{
  prefix = c("v1.","v2.")
  
  v1title = c("E11 V1 raw BC recover rate","E11 V1 reads num(10) filtered BC recover rate",
              "E11 V1 reads num(10) filtered BC-pattern recover rate")
  v1test11 = FilterRecover(v1lsx,cell11,10,v1title)
  v1test11$stat
  v1test11$venn_plot
  
  v2title = c("E11 V2 raw BC recover rate","E11 V2 reads num(10) filtered BC recover rate",
              "E11 V2 reads num(10) filtered BC-pattern recover rate")
  v2test11 = FilterRecover(v2lsx,cell11,10,v2title)
  v2test11$stat
  v2test11$venn_plot
  saveRDS(v1test11,file = "final_result/xE11/e11_v1qc.rds")
  saveRDS(v2test11,file = "final_result/xE11/e11_v2qc.rds")
  v1test11 = readRDS("final_result/xE11/e11_v1qc.rds")
  v2test11 = readRDS("final_result/xE11/e11_v2qc.rds")
  
  ggsave(ggarrange(v1test11$venn_plot,v2test11$venn_plot,nrow = 1),
         filename = "final_result/xE11/E11_recover_rate_stat.pdf",
         width = 14,height = 10)
  
  #
  bl = readRDS("final_result/x.blacklist.rds")
  v1testf11 = BlacklistFilter(v1test11$result,bl$v1_array)
  v2testf11 = BlacklistFilter(v2test11$result,bl$v2_array)
  dataE11 = list(v1testf11,v2testf11)
  outnameE11 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE11/E11_blay"
  E11tree = TreeConstruct(dataE11,prefix,outnameE11,cell11)
  E11tree = CirclePlot1(E11tree,outnameE11)
  E11tree = MyHeatmap(E11tree,outnameE11)
  
  saveRDS(E11tree,file = "final_result/xE11/E11_blay_tree.rds")
  
  v1testf11 = BlacklistFilter(v1test11$result,bl$v1_tag)
  v2testf11 = BlacklistFilter(v2test11$result,bl$v2_tag)
  dataE11 = list(v1testf11,v2testf11)
  outnameE11 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE11/E11_bltg"
  E11treetg = TreeConstruct(dataE11,prefix,outnameE11,cell11)
  E11treetg = CirclePlot1(E11treetg,outnameE11)
  E11treetg = readRDS("final_result/xE11/E11_bltg_tree.rds")
  E11treetg = MyHeatmap(E11treetg,outnameE11)
  saveRDS(E11tree,file = "final_result/xE11/E11_bltg_tree.rds")
  # 
  # saveRDS(v1test,file = "final_result/xE11/E11_v1_result.rds")
  # saveRDS(v2test,file = "final_result/xE11/E11_v2_result.rds")
  
  
  prefix = c("v1.")
  dataE11v1 = list(v1testf11)
  outnameE11v1 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE11/E11_blay_v1"
  E11treev1 = TreeConstruct(dataE11v1,prefix,outnameE11v1,cell11)
  E11treev1 = CirclePlot1(E11treev1,outnameE11v1)
  E11treev1 = readRDS("final_result/xE11/E11_blay_v1_tree.rds")
  E11treev1 = MyHeatmap(E11treev1,outnameE11v1)
  saveRDS(E11treev1,file = "final_result/xE11/E11_blay_v1_tree.rds")
  
  prefix = c("v2.")
  dataE11v2 = list(v2testf11)
  outnameE11v2 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE11/E11_blay_v2"
  E11treev2 = TreeConstruct(dataE11v2,prefix,outnameE11v2,cell11)
  E11treev2 = CirclePlot1(E11treev2,outnameE11v2)
  E11treev2 = readRDS("final_result/xE11/E11_blay_v2_tree.rds")
  E11treev2 = MyHeatmap(E11treev2,outnameE11v2)
  saveRDS(E11treev2,file = "final_result/xE11/E11_blay_v2_tree.rds")
  
}


#STF E11-----------------------------------
#Read data
{
  samplels = list.files("sc_result/STF-E11",full.names = T)[1:3]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = list.files("sc_result/STF-E11")[1:3]
  names(v2lsx) = list.files("sc_result/STF-E11")[1:3]
  
  cell11 = ReadCell("sc_result/STF-E11/STF-E11_celltype.csv")
  
}

#Build Tree
{
  prefix = c("v1.","v2.")
  
  v1title = c("E11 V1 raw BC recover rate","E11 V1 reads num(10) filtered BC recover rate",
              "E11 V1 reads num(10) filtered BC-pattern recover rate")
  v1test11 = FilterRecover(v1lsx,cell11,10,v1title)
  v1test11$stat
  v1test11$venn_plot
  
  v2title = c("E11 V2 raw BC recover rate","E11 V2 reads num(10) filtered BC recover rate",
              "E11 V2 reads num(10) filtered BC-pattern recover rate")
  v2test11 = FilterRecover(v2lsx,cell11,10,v2title)
  v2test11$stat
  v2test11$venn_plot
  saveRDS(v1test11,file = "final_result/STF_E11/e11_v1qc.rds")
  saveRDS(v2test11,file = "final_result/STF_E11/e11_v2qc.rds")
  v1test11 = readRDS("final_result/STF_E11/e11_v1qc.rds")
  v2test11 = readRDS("final_result/STF_E11/e11_v2qc.rds")
  
  ggsave(ggarrange(v1test11$venn_plot,v2test11$venn_plot,nrow = 1),
         filename = "final_result/STF_E11/E11_recover_rate_stat.pdf",
         width = 14,height = 10)
  
  #
  bl = readRDS("final_result/x.blacklist.rds")
  v1testf11 = BlacklistFilter(v1test11$result,bl$v1_array)
  v2testf11 = BlacklistFilter(v2test11$result,bl$v2_array)
  dataE11 = list(v1testf11,v2testf11)
  outnameE11 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/STF_E11/E11_blay"
  E11tree = TreeConstruct(dataE11,prefix,outnameE11,cell11)
  E11tree = CirclePlot1(E11tree,outnameE11)
  E11tree = MyHeatmap(E11tree,outnameE11)
  
  saveRDS(E11tree,file = "final_result/STF_E11/E11_blay_tree.rds")
  
  v1testf11tg = BlacklistFilter(v1test11$result,bl$v1_tag)
  v2testf11tg = BlacklistFilter(v2test11$result,bl$v2_tag)
  dataE11tg = list(v1testf11tg,v2testf11tg)
  outnameE11 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/STF_E11/E11_bltg"
  E11treetg = TreeConstruct(dataE11tg,prefix,outnameE11,cell11)
  E11treetg = CirclePlot1(E11treetg,outnameE11)
  E11treetg = MyHeatmap(E11treetg,outnameE11)
  saveRDS(E11treetg,file = "final_result/STF_E11/E11_bltg_tree.rds")
  # 
  # saveRDS(v1test,file = "final_result/STF_E11/E11_v1_result.rds")
  # saveRDS(v2test,file = "final_result/STF_E11/E11_v2_result.rds")
  
  
  prefix = c("v1.")
  dataE11v1 = list(v1testf11)
  outnameE11v1 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/STF_E11/E11_blay_v1"
  E11treev1 = TreeConstruct(dataE11v1,prefix,outnameE11v1,cell11)
  E11treev1 = CirclePlot1(E11treev1,outnameE11v1)
  E11treev1 = MyHeatmap(E11treev1,outnameE11v1)
  saveRDS(E11treev1,file = "final_result/STF_E11/E11_blay_v1_tree.rds")
  
  prefix = c("v2.")
  dataE11v2 = list(v2testf11)
  outnameE11v2 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/STF_E11/E11_blay_v2"
  E11treev2 = TreeConstruct(dataE11v2,prefix,outnameE11v2,cell11)
  E11treev2 = CirclePlot1(E11treev2,outnameE11v2)
  # E11treev2 = readRDS("final_result/STF_E11/E11_blay_v2_tree.rds")
  E11treev2 = MyHeatmap(E11treev2,outnameE11v2)
  saveRDS(E11treev2,file = "final_result//E11_blay_v2_tree.rds")
  
}

#Build Tree function
bl = readRDS("final_result/x.blacklist.rds")
XBulidTree = function(v1lsx,v2lsx,cell11,outpath,bl){
  prefix = c("v1.","v2.")
  
  v1title = c("E11 V1 raw BC recover rate","E11 V1 reads num(10) filtered BC recover rate",
              "E11 V1 reads num(10) filtered BC-pattern recover rate")
  v1test11 = FilterRecover(v1lsx,cell11,10,v1title)
  v1test11$stat
  # v1test11$venn_plot
  
  v2title = c("E11 V2 raw BC recover rate","E11 V2 reads num(10) filtered BC recover rate",
              "E11 V2 reads num(10) filtered BC-pattern recover rate")
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
  
  #
  v1testf11 = BlacklistFilter(v1test11$result,bl$v1_array)
  v2testf11 = BlacklistFilter(v2test11$result,bl$v2_array)
  dataE11 = list(v1testf11,v2testf11)
  outnameE11 = paste0(outpath,"_blay")
  E11tree = TreeConstruct(dataE11,prefix,outnameE11,cell11)
  E11tree = CirclePlot1(E11tree,outnameE11)
  E11tree = MyHeatmap(E11tree,outnameE11)
  
  saveRDS(E11tree,file = paste0(outpath,"_blay_tree.rds"))
  
  v1testf11tg = BlacklistFilter(v1test11$result,bl$v1_tag)
  v2testf11tg = BlacklistFilter(v2test11$result,bl$v2_tag)
  dataE11tg = list(v1testf11tg,v2testf11tg)
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

#STF E11 DIV7------------------------------
#Read data
{
  samplels = list.files("sc_result/STF-DIV7",full.names = T)[1:3]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = list.files("sc_result/STF-DIV7")[1:3]
  names(v2lsx) = list.files("sc_result/STF-DIV7")[1:3]
  
  celldiv = ReadCell("sc_result/STF-DIV7/STF-DIV7_celltype.csv")
  
}
outpath = "final_result/STF_DIV7/DIV7"
DIVtree = XBulidTree(v1lsx,v2lsx,celldiv,outpath,bl)



#STF E11 filter doublet 
{
  samplels = list.files("sc_result/STF-E11/",full.names = T)[1:3]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = list.files("sc_result/STF-E11")[1:3]
  names(v2lsx) = list.files("sc_result/STF-E11")[1:3]
  
  cell11 = ReadCell("sc_result/STF-E11/STF_E11_singlet_renamed_celltype.csv")
  
}
outpath = "final_result/STF_E11/singlet_E11"
E11Stree = XBulidTree(v1lsx,v2lsx,cell11,outpath,bl)

#E11 filter doublet
{
  samplels = paste0("sc_result/X.E11/",c("lib1","lib2","lib3"))
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = names(v2lsx) = c("lib1","lib2","lib3")
  
  cell11 = ReadCell("sc_result/X.E11/E11_singlet_renamed_celltype.csv")
  
}
outpath = "final_result/xE11_noDoublet/E11_12_21"
E11tree = XBulidTree(v1lsx,v2lsx,cell11,outpath,bl)

#E15 filter doublet
{
  samplels = paste0("sc_result/x.TAQ_E15/",c("AT","TP","TM"))
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T))
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T))
  names(v1lsx) = names(v2lsx) = c("AT","TP","TM")
  cell = ReadCell("sc_result/x.TAQ_E15/E15clustername20211126 v2.csv")
  doublet = read.csv("sc_result/x.TAQ_E15/E15.5_doublet_BC.csv")$x
  cell = cell[which(!cell$Cell.BC %in% substr(doublet,1,16)),]
  # cellt = cell
  # cell[which(cellt$Cell.type == "Rgl1"),"Cell.type"] = "Rgl3"
  # cell[which(cellt$Cell.type == "Rgl3"),"Cell.type"] = "Rgl1"
  # write.csv(cell,file = "sc_result/dual_E15/celltype_subclustered1117.csv",row.names = F,quote = F)
}
outpath = "final_result/xE15_noDoublet/E15"
E15tree = XBulidTree(v1lsx,v2lsx,cell,outpath,bl)




