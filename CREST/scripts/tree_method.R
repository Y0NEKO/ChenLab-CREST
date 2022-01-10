#try heatmap method 
#S2&S3 compare
library(reshape2)
library(dplyr)
library(ggplot2)
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
library(hrbrthemes)
library('pvclust')
library(ComplexHeatmap)
library(circlepackeR)
library(circlize)
library(igraph)
library(ggraph)

test_data = function()
{
  setwd("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/sc_result/")
  samplels = list.files()[c(9,11,12,15,16)]
  #???ļ?--------
  v1ls = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2ls = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  
  celllist = list.files(path = samplels, "/*_celltype.csv",full.names = T)
  celllist2 = celllist;celllist2[4:5] = c("S2_D51/S2_D51_celltype2.csv","S3_D51/S3_D51_celltype2.csv")
  Cells = lapply(celllist, ReadCell)
  names(Cells) = names(v1ls) = names(v2ls) = samplels
  bl1 = read.table("blacklisttag_select_v1.tsv",header = T,sep = "\t")
  bl2 = read.table("blacklisttag_select_v2.tsv",header = T,sep = "\t")
  
  rawcell = c("7-N_Da","0-Prog_FP","3-Prog_NKX2-2","4-Prog_NKX6-1","1-N_Sero","5-N_Glu","6-N_GABA","9-N_Motor",
              "0-FP_SHH", "2-FP_proli","4-Prog_NKX2-9","8-VLMC",
              "5-N_Gaba","2-N_Sero1","3-N_Sero2","1-N_Glu","7-N_Motor2","4-Prog_LFP","6-N_Motor1",
              "1-Meten_Glu","4-DA","5-N_GABA")
  concell = c("DA","Prog_FP","Prog_NKX2-2","Prog_NKX6-1","Sero","Glu","GABA","Motor","FP_SHH","FP_proli","Prog_NKX2-9","VLMC",
              "GABA","Sero1","Sero2","Glu","Motor2","Motor1","Prog_LFP","Meten_Glu","DA","GABA")
  library(stringi)
  for (i in 1:5) {
    Cells[[i]]$Cell.type = stri_replace_all_fixed(Cells[[i]]$Cell.type,rawcell,concell,vectorize_all = FALSE)
  }
  
  
  #?????߼??????Ľ???��??ģ
  #build tree-----------------------------
  dv13 = v1ls$S3_D51
  dv12 = v1ls$S2_D51
  dv23 = v2ls$S3_D51
  dv22 = v2ls$S2_D51
  cell3 = Cells$S3_D51
  cell2 = Cells$S2_D51
  dv13 = merge(dv13[,c("Cell.BC","umim")],cell3,by="Cell.BC");colnames(dv13)[2] = "pattern"
  dv12 = merge(dv12[,c("Cell.BC","umim")],cell2,by="Cell.BC");colnames(dv12)[2] = "pattern"
  dv23 = merge(dv23[,c("Cell.BC","umim")],cell3,by="Cell.BC");colnames(dv23)[2] = "pattern"
  dv22 = merge(dv22[,c("Cell.BC","umim")],cell2,by="Cell.BC");colnames(dv22)[2] = "pattern"
  
  prefix = c("v1.","v2.")
  data3 = list(dv13,dv23)
  data2 = list(dv12,dv22)
  bl = list(bl1,bl2)
  outname3 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/tree_count/tree_test_6_29/S3_test"
  outname2 = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/tree_count/tree_test_6_29/S2_test"
  # outname = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/tree_count/tree_test_6_29/S3_test"
  s3tree = TreeConstruct(data3,bl,prefix,outname3,cell3)
  s2tree = TreeConstruct(data2,bl,prefix,outname2,cell2)
  s3tree = CirclePlot1(s3tree,outname3)
  s2tree = CirclePlot1(s2tree,outname2)
  s2psd = PseudoAna(s2tree,pseudo)
  s3psd = PseudoAna(s3tree,pseudo)
  s2psdst = data.frame("sample" = "S2",
                       "FC" = s2psd$psdfc$pseado.to / s2psd$psdfc$pseado.from) 
  s3psdst = data.frame("sample" = "S3",
                       "FC" = s3psd$psdfc$pseado.to / s3psd$psdfc$pseado.from) 
  psdst = rbind(s2psdst,s3psdst)
  psdst = psdst[which(psdst$FC>1),]
  ggplot(psdst,aes(x = sample,y = log2(FC))) + geom_boxplot()
  
  
  s2tree = readRDS("../tree_count/s2tree.rds")
  s3tree = readRDS("../tree_count/s3tree.rds")
  
  outname2 = "../tree_count/method_test/method1_s2"
  outname3 = "../tree_count/method_test/method1_s3"
  MyHeatmap(s2tree,outname2)
  MyHeatmap(s3tree,outname3)
  #分析s2 s3 编辑深度差异
  {
    
    
    
    leveldfs2 = StatLevel(s2tree$edge,s2tree$vertice,"s2")
    leveldfs3 = StatLevel(s3tree$edge,s3tree$vertice,"s3")
    
    leveldfr = rbind(leveldfr2,leveldfr3)
    # leveldfs2$sample = "s2"
    # leveldfs3$sample = "s3"
    leveldft = rbind(leveldfs2,leveldfs3)
    # leveldft$level = as.factor(as.character(leveldft$level))
    # leveldft$sizet = log2(leveldft$sizet+1)
    
    # p2.1 = ggboxplot(leveldft,x = "level", y ="sizet", color = "sample") + 
    #   stat_compare_means(aes(group = sample),label = "p.format") + ylab("cell numbers") +
    #   theme_ipsum(base_family = "sans")
    # p2.2 = ggboxplot(leveldft,x = "level", y ="celltypes", color = "sample") + 
    #   stat_compare_means(aes(group = sample),label = "p.format") + ylab("cell type numbers") +
    #   theme_ipsum(base_family = "sans")
    # 
    # 
    # levelar = leveldft %>% group_by(sample,level) %>% summarise(count=n())
    # p2.3 = ggbarplot(levelar,x = "level",y = "count",color="sample",position = position_dodge(0.8)) + 
    #   ylab("array number") +
    #   theme_ipsum(base_family = "sans")
    # p2.3
    # 
    # ggexport(ggarrange(p2.1,p2.2,p2.3,nrow = 1,ncol = 1),
    #          filename = "../tree_count/method_test/method7_s2&s3_level_edit_compare.pdf",height = 5,width = 7)
    # 
    length(unique(leveldfs2$arrays))
    length(unique(leveldfs3$arrays))
    
    head(leveldfs2)
    vertice2 = s2tree$vertice
    vertice3 = s3tree$vertice
    vertice2[which(vertice2$cluster == "NONE"),"size"] = 0
    vertice3[which(vertice3$cluster == "NONE"),"size"] = 0
    head(vertice2)
    length(unique(s2tree$cm$tags))/sum(vertice2$size)
    length(unique(s3tree$cm$tags))/sum(vertice3$size)
    sts23 = data.frame("sample" = c("S2","S3"),
                       "array" = c(length(unique(s2tree$cm$tags)),length(unique(s3tree$cm$tags))),
                       "cell" = c(sum(vertice2$size),sum(vertice3$size)))
    sts23$density = sts23$array/sts23$cell
    sts23
    
    
    # leveldfr = leveldfr[order(leveldfr$norm),]
    # leveldfr$rank = 1:nrow(leveldfr)
    
    # ggplot(leveldfr,aes(x=norm, col=sample)) +
    #   stat_ecdf()
    
    p3.1 = ggplot(leveldfr,aes(x=norm,y = ecdf, col=sample)) +
      geom_point(aes(size = complexity),alpha = 0.5) +
      scale_color_ipsum() +
      theme_bw(base_family = "sans") +
      ylab("ecdf of subcluster size") + xlab("normalize size of subcluster")
    p3.1
    
    p3.2 = ggplot(leveldfr,aes(x=norm,y = complexity, col=sample)) +
      geom_point(alpha = 0.5) + 
      scale_color_ipsum() +
      theme_bw(base_family = "sans") +
      ylab("complexity of subcluster") + xlab("normalize size of subcluster")
    
    sts23 = melt(sts23[1:3])
    p3.3 = ggplot(sts23,aes(x=variable, y = value, fill=sample)) +
      geom_bar(stat = "identity",position = position_dodge(),width = 0.5) + 
      scale_fill_ipsum() +
      theme_bw(base_family = "sans") +
      ylab("") + xlab("")
    p3.3
    
    
    ggexport(ggarrange(p3.1,p3.2,p3.3,nrow = 1,ncol = 1),
             filename = "../tree_count/method_test/s2&s3_compare.pdf",height = 5,width = 7)
    
    
    ggplot(leveldfr,aes(x=rank,y = norm, col=sample)) +
      geom_point(aes(size = complexity),alpha = 0.5) +
      ylab("normalize subcluster size") + xlab("subcluster rank by size")
    
    # ggplot(leveldfr,aes(x=rank,y = norm, col=sample)) +
    #   geom_point(aes(size = complexity),alpha = 0.5) +
    #   ylab("normalize subcluster size") + xlab("subcluster rank by size")
    # scale_fill_distiller(palette = "RdPu")
    
    
    
    leveldfr.e = NULL
    for (i in 1:nrow(leveldfr)) {
      leveldfr.e = rbind(leveldfr.e,leveldfr[rep(i,leveldfr$sizet[i]),])
    }
    
    ggplot(leveldfr.e,aes(x=norm, col=sample)) + geom_histogram(binwidth = 0.001)
    
    ggplot(leveldfr[which(leveldfr$norm>0.03),],aes(x=sample, y=complexity)) + geom_boxplot()
    #计算 complexity
    
    write.csv(leveldft,file = "../final_result/s2_s3/s2_s3_compare_leveldf.csv",row.names = F,quote = F)
    write.csv(leveldfr,file = "../final_result/s2_s3/s2_s3_compare_leveldf_root.csv",row.names = F,quote = F)
    
  }
}

ReadCell = function(cellfile){
  Cells<-read.csv(cellfile,header = T,stringsAsFactors = F)
  colnames(Cells) = c("Cell.BC","Cell.type")
  Cells$Cell.BC = substring(Cells$Cell.BC,1,16)
  return(Cells)
}
CircleGroup = function(edge){
  root = unique(edge[which(edge$from == "node.N0"),"to"])
  groupdf = NULL
  for (i in 1:length(root)) {
    level = 0
    # lt = arrays[1]
    lt = root[i]
    lti = lt
    while (1) {
      lti = edge[which(edge$from %in% lti),"to"]
      if(length(lti) == 0){break}
      groupdf = rbind(groupdf,data.frame("tags" = lti,"group" = lt))
    }
    print(i)
  }
  return(groupdf)
}
CirclePlot1 = function(tree,outname){
  #plot circle
  cm2 = tree$cm2
  edge = tree$edge
  edgegrp = CircleGroup(edge)
  celllist = c("NONE",unique(cm2$celltype)[order(unique(cm2$celltype))])
  tree_cell_count = cm2 %>% group_by(tags, celltype) %>% dplyr::summarise(count = sum(n()))
  
  nodesize = tree_cell_count
  # nodesize = s2tree$count
  # nodesize$size = rowSums(s2tree$count[-1])
  colnames(nodesize) = c("tags","cluster","size")
  
  nodesize = nodesize[which(nodesize$tags %in% unique(c(edge$from,edge$to))),]
  N0 = data.frame("tags" = unique(c(edge$from[!(edge$from %in% nodesize$tags)],
                                    edge$to[!(edge$to %in% nodesize$tags)])), 
                  "cluster" = "NONE","size" = 0)
  vertices = rbind(N0,nodesize)
  # vertices$group = NA
  # vertices[which(vertices$tags %in% edgegrp$group),"group"] = vertices[which(vertices$tags %in% edgegrp$group),"tags"]
  # vertices[which(vertices$tags == "node.N0"),"size"] = sum(nodesize$size)
  # vertices$shortname = substr(vertices$tags,1,10)
  # Rebuild the graph object
  #filter edges
  
  mygraph <- graph_from_data_frame(edge, vertices=vertices)
  ### convert to igraph
  
  collist = c("white","#ffadad","#ffd6a5","#fdffb6",
              "#caffbf","#9bf6ff","#a0c4ff","#bdb2ff",
              "#ffc6ff","#b5838d","#230805")
  if(length(celllist)>length(collist)){collist = viridis_pal(option = "D")(length(celllist)-1); collist = c("white",collist)}
  p1 = ggraph(mygraph, 'circlepack', weight = size) +
    geom_node_circle(aes(fill = factor(cluster,levels = celllist)),size = 0.1) +
    scale_fill_manual(values = collist) +
    coord_fixed() +
    labs(fill='cluster') +
    # scale_x_continuous(limits = c(-500,500)) +
    # scale_y_continuous(limits = c(-500,500)) +
    theme_void() + guides(color=guide_legend(title="celltype"))
  p1
  
  
  p2 = ggraph(mygraph, 'circlepack', weight = log2(size+1)) +
    geom_node_circle(aes(fill = factor(cluster,levels = celllist)),size = 0.1) +
    scale_fill_manual(values = collist) +
    coord_fixed() +
    labs(fill='cluster') +
    theme_void() + guides(color=guide_legend(title="celltype"))
  p2
  
  ggsave(p1,filename = paste0(outname,"_circpack.pdf"))
  ggsave(p2,filename = paste0(outname,"_circpack_size.pdf"))
  tree$vertice = vertices
  tree$figure.circpack = p1
  tree$figure.circpack_log = p2
  return(tree)
}
CirclePlot2 = function(tree,pseudo,outname){
  #plot circle
  celllist = c("NONE",unique(tree$cm2$celltype))
  cm2 = tree$cm2
  edge = tree$edge
  cmt = merge(cm2,pseudo,by.x = "Var1",by.y = "BC")
  tree_cell_count = cmt %>% group_by(tags, celltype) %>% dplyr::summarise(count = sum(Freq),pseudotime = mean(pseudo))
  
  nodesize = tree_cell_count
  # nodesize = s2tree$count
  # nodesize$size = rowSums(s2tree$count[-1])
  colnames(nodesize) = c("tags","cluster","size","pseudotime")
  
  nodesize = nodesize[which(nodesize$tags %in% unique(c(edge$from,edge$to))),]
  N0 = data.frame("tags" = unique(c(edge$from[!(edge$from %in% nodesize$tags)],
                                    edge$to[!(edge$to %in% nodesize$tags)])), 
                  "cluster" = "NONE","size" = 1,"pseudotime" = 0)
  vertices = rbind(N0,nodesize)
  vertices$shortname = substr(vertices$tags,1,10)
  # Rebuild the graph object
  #filter edges
  
  mygraph <- graph_from_data_frame(edge, vertices=vertices)
  ### convert to igraph
  
  
  p = ggraph(mygraph, 'circlepack', weight = size) +
    # p = ggraph(mygraph, 'circlepack') +
    geom_node_circle(aes(fill = pseudotime)) +
    scale_fill_gradient2(low = "#FFFDE4",high = "#005AA7")+
    geom_node_circle(aes(color = factor(cluster,levels = celllist))) +
    scale_color_manual(values = c("black","#ffadad","#ffd6a5","#fdffb6",
                                  "#caffbf","#9bf6ff","#a0c4ff","#bdb2ff","#ffc6ff","#b5838d")) +
    # # geom_node_text( aes(label=shortname, filter=leaf, size=size)) +
    # geom_node_text( aes(label=shortname),size=3) +
    coord_fixed() + labs(fill='cluster') +
    theme_void() + guides(color=guide_legend(title="celltype"),
                          fill=guide_legend(title="pseudotime"))
  p
  
  ggsave(p,filename = paste0(outname,"_pseudo_circpack.pdf"))
  tree$vertice = vertices
  tree$figure.circpack = p
  return(tree)
  
}
BlacklistFilter = function(data,bl){
  df = data
  newdf = NULL
  for (j in 1:nrow(df)) {
    x = df$pattern[j]
    flag = 0
    flag = length(bl[,1][which(bl[,1] == x)])
    if(flag == 0){
      x_str = strsplit(x,split = "_")
      x_pre = x_str[[1]]
      x_pre = unique(x_pre[!(x_pre %in% "NONE")])
      flag = length(intersect(x_pre,bl[,1]))
    }
      
    if(flag==0){
      newdf = rbind(newdf,df[j,])
    }
    
  }
  return(newdf)
  
}
BlacklistFilter2 = function(data,bl){
  df = data
  newdf = NULL
  for (j in 1:nrow(df)) {
    x = df$umim[j]
    x_str = strsplit(x,split = "_")
    x_pre = x_str[[1]]
    x_pre = unique(x_pre[!(x_pre %in% "NONE")])
    x_pre = x_pre[which(!x_pre %in% bl[,1])]
    if(length(x_pre)>0){
      df$umim[j] = paste0(x_pre,collapse =  "_")
      newdf = rbind(newdf,df[j,])
    }
  }
  return(newdf)
  
}
TreeDataInitial = function(data,bl,cell){
  
}

TreeConstruct = function(data,prefix,outname,cell){
  TagStat = function(x) {
    x = as.character(x)
    t = x[2]
    t = unlist(strsplit(t,"_|&"))
    # t = t[!t%in%c("NONE",blt)]
    t = t[!t%in%c("NONE")]
    t = unique(t)
    if(length(t) == 0){
      return(NA)
    }else{
      return(data.frame(Cell.BC = x[1],Tag = t,Cell.type = x[3]))
    }
  }
  
  tag = NULL
  
  for (i in 1:length(data)) {
    # if(length(bl)>0){
    #   tagi = apply(data[[i]],1,TagStat,blt = bl[[i]]$tag)
    # }else{
    #   tagi = apply(data[[i]],1,TagStat,blt = "")
    # }
    tagi = apply(data[[i]],1,TagStat)
    tagi = do.call("rbind",tagi)
    tagi = na.omit(tagi)
    # if(length(bl)>0){
    #   tagi = tagi[which(!tagi$Tag %in% bl[[i]]$tag),]
    # }
    #black list filter
    tagi$Tag = paste(prefix[i], tagi$Tag, sep = "")
    tag = rbind(tag,tagi)
    
  }
  
  
  
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
  Tag_1 = Tag[Tag$Var1 %in% tags_tab[tags_tab$num == 1,1],]
  Tag_1 = Tag_1[order(-Tag_1$Freq),]
  
  Tag_2 = Tag
  Tag_2 = Tag_2[which(Tag_2$Var1 %in% setdiff(Tag_2$Var1,Tag_1$Var1)),]
  Tag_2 = Tag_2[order(-Tag_2$Freq),]
  Tag_2$num = 0
  colnames(Tag_2) = colnames(tags_tab)
  
  tags_tab[tags_tab$num == 1,] = tags_tab[tags_tab$num == 1,][match(as.character(Tag_1$Var1),
                                                                    as.character(tags_tab$tags_paste[tags_tab$num==1])),]
  
  tags_tab = rbind(tags_tab,Tag_2)
  tags_uni = lapply(as.character(tags_tab$tags_paste),function(x){strsplit(x,split = "_")[[1]]})
  #node build
  cluster_stat = function(i){
    x = tags_uni[[i]]
    n = tags_tab$num[i]
    tags_belone = NA
    for(y_ind in which(tags_tab$num<n)){
      y=tags_uni[[y_ind]]
      if(length(intersect(x,y))>0){
        if(length(setdiff(y,x))==0){
          tags_belone = y_ind
          break
        }
        
      }else{
        next
      }
    }
    return(tags_belone)
  }
  
  belons = sapply(which(tags_tab$num > 1),cluster_stat)
  belons[tags_tab$num <= 1]= NA
  
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
  
  # nodes = nodes[!is.na(nodes$N2),] 删去只有一个层的节点（为了看单个节点中包含的细胞关系，不删了）
  # bn = names(table(nodes$N1)[which(table(nodes$N1)==1)])
  # nodes = nodes[!nodes$N1 %in% bn,]
  
  nodes = nodes[which(!(is.na(nodes$N2) & (nodes$N1 %in% as.character(Tag_2$tags_paste)))),]
  
  #save tree figure and rds
  population = as.Node(nodes)
  saveNetwork(diagonalNetwork(ToListExplicit(population, unname = TRUE),
                              margin = 10,fontSize = 8,height = 7000),
              file = paste0(outname,".html"))
  
  saveRDS(population,paste0(outname,".rds"))
  
  #save celltype tab
  cell_tab$tags = tags_paste
  cell_tab$celltype = cell$Cell.type[match(as.character(cell_tab$Var1),cell$Cell.BC)]
  write.csv(cell_tab,paste0(outname,"_cell_tab.csv"),row.names = F,quote = F)
  
  tree = population
  cm = cell_tab
  
  print("data build done")
  print("ploting...")
  
  #plot
  {
    
    #1.trans data.tree to newick with internal node
    td = ToDataFrameNetwork(tree)
    from_node = unique(td$from)
    td[which(td$to %in% from_node),"to"] = paste("node.",td[which(td$to %in% from_node),"to"],sep = "")
    from_node = from_node[which(from_node %in% tags_tab[which(tags_tab$num > 0),1])]
    td = rbind(data.frame("from" = from_node, "to" = from_node),td)
    td$from = paste("node.",td$from,sep = "")
    td = td[which(td$to != "N0"),]
    td_wc = NULL
    for (i in 1:nrow(td)) {
      arrow = td[i,]
      mc = unique(cm[which(arrow$to == cm$tags),3:4])
      td_wc = rbind(td_wc,arrow)
      if(nrow(mc)>0){
        td_wc = rbind(td_wc,data.frame("from" = arrow$to,"to" = paste0(mc$tags,"-",mc$celltype)))
      }
      
    }
    
    root = td_wc[which(td_wc$from == "node.N0"),"to"]
    whitetag = NULL
    for (i in 1:length(root)) {
      line = root[i]
      linenw = td_wc[which(td_wc$from == line),]
      while(1)
      {
        line = c(line,linenw$to)
        if(nrow(linenw) > 1){break}
        linenw = td_wc[which(td_wc$from == linenw$to),]
        if(nrow(linenw)==0){whitetag = c(whitetag,line);break}
      }
      
    }
    td_wc2 = td_wc[which(!(td_wc$from %in% whitetag) & !(td_wc$to %in% whitetag)),]
    
    td_node = FromDataFrameNetwork(td_wc2)
    tree_nwk = ToNewick(td_node)
    write(tree_nwk,paste(outname,"nwk",sep = "."))
    
    
    #2.data build
    #(1)tree_data for plotting tree
    tree_data = read.tree(paste(outname,"nwk",sep = "."))
    tl = tree_data$tip.label
    
    #(2)tree_cell_plot for plotting cell composition 
    cm2 = cm;cm2$tags = paste0(cm2$tags,"-",cm2$celltype)
    tree_cell = cm2 %>% group_by(tags, celltype) %>% dplyr::summarise(count = sum(Freq))
    # tree_cell = dcast(tree_cell,tags~celltype,value.var = "count")
    # tree_cell[is.na(tree_cell)] = 0
    tree_cell = tree_cell %>% group_by(celltype) %>% summarise(tags = tags,count=count,norm = count/sum(count))
    tree_cell = tree_cell[,c(2,1,3,4)]
    
    #get
    tree_cell_count = cm %>% group_by(tags, celltype) %>% dplyr::summarise(count = sum(Freq))
    tree_cell_norm = tree_cell_count %>% group_by(celltype) %>% summarise(tags = tags,count=count,norm = count/sum(count))
    tree_cell_count = dcast(tree_cell_count,tags~celltype,value.var = "count")
    tree_cell_norm = dcast(tree_cell_norm,tags~celltype,value.var = "norm")
    tree_cell_count[is.na(tree_cell_count)] = 0
    tree_cell_norm[is.na(tree_cell_norm)] = 0
    tree_cellm = tree_cell_count
    tdm = ToDataFrameNetwork(tree)
    tdm = tdm[which(tdm$to != "N0"),]
    td_cell = merge(tdm,tree_cellm, by.x = "from",by.y = "tags",all.x = T)
    td_cell = merge(td_cell,tree_cell, by.x = "to",by.y = "tags",suffixes = c(".from",".to"),all.x = T)
    td_cell = td_cell[,c(2,1,3:ncol(td_cell))]

    
    #(3)indel for plotting pattern
    #library(stringr)
    
    indel = data.frame(id=NA,start=NA,end=NA,width=NA,type=NA)
    indel = indel[-1,]
    #build indel data frame
    for(i in c(1:length(tl))){
      x = tl[i]
      x1 = strsplit(x,split = "-")[[1]][1]
      x_str = strsplit(x1,split = "_")
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
    p0 = ggtree(tree_data_sort, layout="circular", branch.length = "none") %<+% tree_cell + 
      aes(color=celltype) + geom_tippoint(aes(size=count), alpha=.6) +
      geom_tiplab(label=tree_cell$tags,offset=.1,size=0.5)
    ggsave(p0,filename = paste0(outname,"_circle.pdf"),height = 15,width = 15)
    
    
    
    p = ggtree(tree_data_sort,size = 0.1, ladderize=F) + geom_tiplab(size = 1) +
      xlim_expand(c(0, 200),panel = "Tree")
    
    pnt = ggtree(tree_data_sort,size = 0.1, ladderize=F) + theme_tree()
    
    p1 = p + geom_facet(panel = "cellnum",data = tree_cell,
                        geom = geom_barh,
                        mapping = aes(x = norm,fill = celltype),stat='identity') +
      #    scale_fill_viridis_d(option="D", name="discrete\nvalue") +
      new_scale_color() +
      geom_facet(panel = "indel pattern", data = indel,
                 geom = geom_segment,
                 mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                 size = 0.3)
    
    p1nt = pnt + geom_facet(panel = "cellnum",data = tree_cell,
                        geom = geom_barh,
                        mapping = aes(x = norm,fill = celltype),stat='identity') +
      #    scale_fill_viridis_d(option="D", name="discrete\nvalue") +
      new_scale_color() +
      geom_facet(panel = "indel pattern", data = indel,
                 geom = geom_segment,
                 mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                 size = 0.3)
    
    p2 = p + geom_facet(panel = "",data = tree_cell,
                        geom = geom_tile,
                        mapping = aes(x = 1,fill = norm,color = celltype)) +
      # scale_fill_viridis_d(option="D", name="discrete\nvalue") +
      scale_fill_gradient2(low = "#FFFDE4",high = "#005AA7")+
      new_scale_color() +
      geom_facet(panel = "indel pattern", data = indel,
                 geom = geom_segment,
                 mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                 size = 0.3)
    
    p2nt = pnt + geom_facet(panel = "",data = tree_cell,
                        geom = geom_tile,
                        mapping = aes(x = 1,fill = norm,color = celltype)) +
         # scale_fill_viridis_d(option="D", name="discrete\nvalue") +
      scale_fill_gradient2(low = "#FFFDE4",high = "#005AA7")+
      new_scale_color() +
      geom_facet(panel = "indel pattern", data = indel,
                 geom = geom_segment,
                 mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                 size = 0.3)
    # p1 = facet_plot(p,panel = "cell composition",data = tree_cell_plot,
    #                 geom = geom_barh,
    #                 mapping = aes(x = value,fill = variable),
    #                 stat="identity") %>%
    #   facet_plot(panel = "indel pattern", data = indel,
    #              geom = geom_segment,
    #              mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
    #              size = 1)
    
    
    #change facet grid
    gt1 = ggplot_gtable(ggplot_build(p1))
    # gt1$widths[7] = 0.5*gt1$widths[7] # in this case it was colmun 7 - reduce the width by a half
    # gt1$widths[9] = 0.5*gt1$widths[9] # in this case it was colmun 7 - reduce the width by a half
    # grid.draw(gt1) # plot with grid draw
    gt1nt = ggplot_gtable(ggplot_build(p1nt))
    
    gt2 = ggplot_gtable(ggplot_build(p2))
    gt2$widths[7] = 0.1*gt2$widths[7] # in this case it was colmun 7 - reduce the width by a half
    grid.draw(gt2) # plot with grid draw
    
    gt2nt = ggplot_gtable(ggplot_build(p2nt))
    gt2nt$widths[7] = 0.1*gt2nt$widths[7] # in this case it was colmun 7 - reduce the width by a half
    grid.draw(gt2nt) # plot with grid draw
    
    ggexport(ggarrange(gt1,gt1nt,gt2,gt2nt,nrow = 1,ncol = 1),filename = paste0(outname,"_norm.pdf"),height = 0.12*length(tl),width = 12,units = "cm",limitsize = F)
    # ggsave(gt,filename = paste0(outname,"_norm2.pdf"),height = 0.12*length(tl),width = 15,units = "cm",limitsize = F)
    
    
  }
  return(list("tree" = tree,"cm" = cm,"cm2" = cm2,"count" = tree_cell_count,
              "norm" = tree_cell_norm,"td" = td_cell,"edge" = td_wc2))
}
PseudoBuild = function(){
  #test for trajectories
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  myPaths <- .libPaths()  
  new <- c('/picb/sysgenomics2/people/liuhengxin/software/anaconda3/envs/r4-base/lib/R/library/')
  myPaths <- c(myPaths, new) 
  .libPaths(myPaths) 
  .libPaths()
  # BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
  #                        'limma', 'S4Vectors', 'SingleCellExperiment',
  #                        'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
  # devtools::install_github('cole-trapnell-lab/leidenbase')
  # devtools::install_github('cole-trapnell-lab/monocle3')
  # remotes::install_github('satijalab/seurat-wrappers')
  
  
  trans = readRDS("../raw_data/S2S3_seurat_1212.rds")
  pdf("../tree_count/method_test/umap_s2&s3.pdf",width = 8,height = 6)
  DimPlot(trans,label = T)
  dev.off()
  trans.cds <- as.cell_data_set(trans)
  trans.cds <- cluster_cells(cds = trans.cds, reduction_method = "UMAP")
  trans.cds <- learn_graph(trans.cds, use_partition = TRUE)
  trans.cds = order_cells(trans.cds, reduction_method = "UMAP")
  pdf("../tree_count/method_test/pseudo_s2&s3.pdf",width = 8,height = 6)
  plot_cells(
    cds = trans.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  dev.off()
  trans <- AddMetaData(
    object = trans,
    metadata = trans.cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "pseudotime"
  )
  pseudo = trans.cds@principal_graph_aux@listData$UMAP$pseudotime
  trans = saveRDS(trans,"../tree_count/method_test/s2&s3_trans.rds")
  
  pseudo = as.data.frame(pseudo)
  pseudo$BC = substr(rownames(pseudo),1,16)
  pseudo = pseudo[which(pseudo$pseudo != Inf),]
}

#分析pseudo与深度的关系
PseudoAna = function(tree,pseudo){
  cmt = merge(tree$cm,pseudo,by.x = "Var1",by.y = "BC")
  #计算深度和pseudo的关系
  vertice = tree$vertice
  edge = tree$edge
  
  StatLevel = function(edge){
    arrays = unique(c(edge$from,edge$to))
    arrays = arrays[which(arrays!="node.N0")]
    leveldf = NULL
    for (i in 1:length(arrays)) {
      level = 0
      # lt = arrays[1]
      lt = arrays[i]
      lti = lt
      while (1) {
        lti = edge[which(edge$to == lti),"from"]
        level = 1 + level
        if(lti == "node.N0"){break}
      }
      print(i)
      leveldf = rbind(leveldf, data.frame("arrays" = lt, "level" = level))
    }
    return(leveldf)
  }
  
  leveldf = StatLevel(tree$edge)
  levp = merge(vertice,leveldf,by.x = "tags",by.y = "arrays")
  levp = levp[which(levp$pseudotime>0),]
  p3.1 = ggplot(levp,aes(x = as.character(level),y = pseudotime,color = cluster)) + geom_point() + theme_bw()
  p3.1
  # ggsave(p3.1,filename = "../tree_count/method_test/method6_pseudo_level_s3.pdf",width = 6,height = 4)
  
  #计算一个父子比值
  tree$td
  psdmean = cmt %>% group_by(tags,celltype) %>% summarise(pseado = mean(pseudo),count = sum(Freq))
  td = tree$td
  td = td[,which(!colnames(td)%in%c("celltype","count","norm"))]
  tmp = merge(td,psdmean,by.y = "tags",by.x = "from",all=T,)
  tmp = merge(tmp,psdmean,by.y = "tags",by.x = "to",all=T,suffixes = c(".from",".to"))
  psdfc = tmp
  psdfc = psdfc[which(!is.na(psdfc$pseado.from) & !is.na(psdfc$pseado.to)),]
  psdfc$count = (psdfc$count.from + psdfc$count.to)/(abs(psdfc$count.to - psdfc$count.from)+1)
  
  p3.2 = ggplot(psdfc,aes(x = pseado.from, y = pseado.to)) + 
    geom_point(aes(pch = celltype.from,color = celltype.to,size = count)) + 
    geom_function(lty = 3,fun = function(x) x) +
    theme_bw()
  p3.2
  # ggsave(p3.2,filename = "../tree_count/method_test/method6_pseudo_parent&child_s3.pdf",width = 8,height = 6)
  return(list("levp" = levp,"psdfc" = psdfc,"figure.depth" = p3.1,"figure.parchi" = p3.2))
  
}

#heatmap normalize,filter,spearman&edu
MyHeatmapNv = function(tree,outname){
  library(amap)
  tree_cell_countm = tree$count;rownames(tree_cell_countm) = tree_cell_countm$tags ;tree_cell_countm = tree_cell_countm[-1]
  # tree_cell_countm = tree_cell_countm[,names(which(colSums(tree_cell_countm)>30))]
  tree_cell_countm = as.matrix(tree_cell_countm)
  tree_cell_countm[which(tree_cell_countm>0)] = 1
  
  tree_cell_norm = tree_cell_countm[which(rowSums(tree_cell_countm>0)>1),]
  tree_cell_norm = log(tree_cell_norm+1)
  tree_cell_norm = t(apply(tree_cell_norm, 1, function(x) {x/sum(x)}))
  
  # add rectangles around groups highly supported by the data
  #filter tree norm
  
  #euclidean distance
  mx1 = as.matrix(dist(t(tree_cell_norm),method = "euclidean"))
  mx1 = 1 - (mx1 - min(mx1))/max(mx1 - min(mx1))
  pdf(paste0(outname,"_euclidean.pdf"),width = 10,height = 8)
  p1 = Heatmap(as.matrix(mx1),heatmap_legend_param = list(title = "euclidean"),
               clustering_distance_rows = "euclidean",
               # clustering_distance_rows = function(x){x = dd1; return(x)},
               clustering_distance_columns = "euclidean",
               # clustering_distance_columns = function(x){x = dd1; return(x)},
               clustering_method_columns = "complete",clustering_method_rows = "complete")
  print(p1)
  dev.off()
  #spearman
  mx2 = cor(tree_cell_norm,method = "spearman")
  mx2[is.na(mx2)] = 0
  pdf(paste0(outname,"_spearman.pdf"),width = 10,height = 8)
  p2 = Heatmap(as.matrix(mx2),heatmap_legend_param = list(title = "spearman"),
               clustering_distance_rows = "spearman",
               # clustering_distance_rows = function(x){x = dd2; return(x)},
               clustering_distance_columns = "spearman",
               # clustering_distance_columns = function(x){x = dd2; return(x)},
               clustering_method_columns = "complete",clustering_method_rows = "complete")
  print(p2)
  dev.off()
  
  dd1 = dist(mx1,method = "euclidean")
  # dd1 = dist(t(tree_cell_norm),method = "euclidean")
  hc = hclust(dd1, method = "complete")
  dd2 = Dist(mx2, method = "spearman")
  # dd2 = Dist(t(tree_cell_norm), method = "spearman")
  hc2 = hclust(dd2, method = "complete")
  
  pdf(paste0(outname,"_pvcluster_tree.pdf"))
  plot(hc, cex = 1,hang=-1,main = "euclidean") # display pv cluster dendogram
  plot(hc2, cex = 1,hang=-1,main = "spearman") # display pv cluster dendogram
  dev.off()
  
  tree$distfigure = list("euc" = p1,"spear" = p2,"euchc" = hc,"spearhc" = hc2)
  return(tree)
}
MyHeatmap = function(tree,outname){
  library(amap)
  tree_cell_countm = tree$count;rownames(tree_cell_countm) = tree_cell_countm$tags;tree_cell_countm = tree_cell_countm[-1]
  tree_cell_countm = as.matrix(tree_cell_countm)
  tree_cell_countm[which(tree_cell_countm>0)] = 1
  
  tree_cell_norm = tree_cell_countm[which(rowSums(tree_cell_countm>0)>1),]
  tree_cell_norm = log(tree_cell_norm+1)
  tree_cell_norm = t(apply(tree_cell_norm, 1, function(x) {x/sum(x)}))
  
  pdf(paste0(outname,"_abspearson.pdf"),width = 10,height = 8)
  
  dd = Dist(t(tree_cell_norm), method = "abspearson")
  mx = as.matrix(1 - dd)
  for (i in 1:nrow(mx)) {
    mx[i,i] = 1
  }
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               # clustering_distance_rows = "spearman",
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               clustering_distance_rows = function(x){x = dd; return(x)},
               # clustering_distance_columns = "spearman",-
               clustering_distance_columns = function(x){x = dd; return(x)},
               clustering_method_columns = "complete",clustering_method_rows = "complete")
  print(p1)
  dev.off()
  
  hc = hclust(dd, method = "complete")
  
  pdf(paste0(outname,"_pvcluster_tree.pdf"))
  plot(hc, cex = 1,hang=-1,main = "abspearson") # display pv cluster dendogram
  dev.off()
  
  tree$distfigure = list("abspearson" = p1,"abspearsonhc" = hc)
  return(tree)
}
#heatmap fisher p value
CalEnrichScore = function(cellmx){
  cmt = melt(cellmx)
  colnames(cmt) = c("tags","celltype","counts")
  Mt = cmt %>% group_by(celltype) %>% summarise(counts = sum(counts))
  Nt = sum(Mt$counts)
  
  
  cellft = NULL
  for (i in 1:nrow(cellmx)) {
    line = cellmx[i,]
    n = sum(line[2:ncol(line)])
    ftl = line
    for (j in 2:ncol(line)) {
      k = unlist(line[j])
      if(k>0){
        M = Mt[which(Mt$celltype == names(k)),]$counts
        N = Nt
        d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
        row.names(d) <- c("In_category", "not_in_category")
        ft = fisher.test(d)
        ftl[j] = -log10(ft$p.value)
      }
    }
    cellft = rbind(cellft,ftl)
  }
  return(cellft)
}
MyHeatmap2 = function(tree,outname){
  cellft = CalEnrichScore(tree$cm)
  rownames(cellft) = cellft$tags; cellft = cellft[-1]
  #filter tree norm
  cellft = cellft[which(rowSums(cellft)>0),]
  
  
  fit2 <- pvclust(cellft, method.hclust="ward",
                  method.dist="euclidean")
  pdf(paste0(outname,"_fisher_pvclust_tree.pdf"))
  plot(fit2, cex = 1) # display pv cluster dendogram
  dev.off()
  # add rectangles around groups highly supported by the data
  
  
  #euclidean distance
  mx1 = as.matrix(dist(t(cellft),method = "euclidean"))
  mx1 = 1 - (mx1 - min(mx1))/max(mx1 - min(mx1))
  pdf(paste0(outname,"_fisher_euclidean.pdf"),width = 10,height = 8)
  p1 = Heatmap(as.matrix(mx1),heatmap_legend_param = list(title = "euclidean"),clustering_distance_rows = "euclidean",
               clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2")
  print(p1)
  dev.off()
  #spearman
  mx2 = cor(cellft,method = "spearman")
  mx2[is.na(mx2)] = 0
  pdf(paste0(outname,"_fisher_spearman.pdf"),width = 10,height = 8)
  p2 = Heatmap(as.matrix(mx2),heatmap_legend_param = list(title = "spearman"),clustering_distance_rows = "spearman",
               clustering_distance_columns = "spearman",clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2")
  print(p2)
  dev.off()
  tree$distmx2 = mx1
  return(tree)
}

StatLevel = function(edge,vertice,sample){
  arrays = unique(c(edge$from,edge$to))
  arrays = arrays[which(arrays!="node.N0")]
  leveldf = NULL
  vertice[which(vertice$cluster == "NONE"),"size"] = 0
  for (i in 1:length(arrays)) {
    level = 0
    size = 0
    # lt = arrays[1]
    lt = arrays[i]
    size = vertice[which(vertice$tags == lt),"size"]
    celltype = length(vertice[which(vertice$tags == lt & vertice$cluster != "NONE"),"cluster"])
    complexity = 0
    lti = lt
    ltj = lt
    ltk = lt
    while (1) {
      lti = edge[which(edge$to == lti),"from"]
      level = 1 + level
      size = size + vertice[which(vertice$tags == lti),"size"]
      if(lti == "node.N0"){break}
    }
    while (1) {
      ltj = edge[which(edge$from %in% ltj),"to"]
      if(length(ltj)>0){
        celltype = celltype  + length(unique(vertice[which(vertice$tags %in% ltj & vertice$cluster != "NONE"),"cluster"]))
        size = size + sum(vertice[which(vertice$tags %in% ltj),"size"])
        complexity = complexity  + length(vertice[which(vertice$tags %in% ltj & vertice$cluster != "NONE"),"tags"])
      }
      else{break}
    }
    
    
    print(i)
    leveldf = rbind(leveldf, data.frame("arrays" = lt, "level" = level, "sizet" = size,
                                        "celltypes" = celltype,"complexity" = complexity,
                                        "sample" = sample))
  }
  result = merge(leveldf,vertice,by.y = "tags",by.x = "arrays")
  return(result)
}

LevelRootExtract = function(s2tree,leveldfs2){
  edge2 = s2tree$edge
  
  root = edge2[which(edge2$from == "node.N0"),]
  edge2 = edge2[which(edge2$from %in% root$to),]
  leveldfr2 = leveldfs2[which(leveldfs2$arrays %in% edge2$from),]
  leveldfr2$norm = leveldfr2$sizet / sum(leveldfr2$sizet)
  
  leveldfr2 = leveldfr2[order(leveldfr2$sizet),]
  
  leveldfr2$ecdf = 0
  for (i in 1:nrow(leveldfr2)) {
    leveldfr2$ecdf[i] = sum(leveldfr2$norm[1:i])
  }
  
  return(leveldfr2)
  
}







