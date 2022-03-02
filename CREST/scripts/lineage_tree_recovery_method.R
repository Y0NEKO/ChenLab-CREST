#lineage tree recovery analysis functions
#22.03.02 by lhx

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
  colnames(nodesize) = c("tags","cluster","size")
  
  nodesize = nodesize[which(nodesize$tags %in% unique(c(edge$from,edge$to))),]
  N0 = data.frame("tags" = unique(c(edge$from[!(edge$from %in% nodesize$tags)],
                                    edge$to[!(edge$to %in% nodesize$tags)])), 
                  "cluster" = "NONE","size" = 0)
  vertices = rbind(N0,nodesize)
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
    data[[i]] = data[[i]][which(data[[i]]$pattern %in% names(table(data[[i]]$pattern))[table(data[[i]]$pattern) > 1]),]
    tagi = apply(data[[i]],1,TagStat)
    tagi = do.call("rbind",tagi)
    tagi = na.omit(tagi)
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

    
    
    #change facet grid
    gt1 = ggplot_gtable(ggplot_build(p1))
    gt1nt = ggplot_gtable(ggplot_build(p1nt))
    gt2 = ggplot_gtable(ggplot_build(p2))
    gt2$widths[7] = 0.1*gt2$widths[7] # in this case it was colmun 7 - reduce the width by a half
    grid.draw(gt2) # plot with grid draw
    
    gt2nt = ggplot_gtable(ggplot_build(p2nt))
    gt2nt$widths[7] = 0.1*gt2nt$widths[7] # in this case it was colmun 7 - reduce the width by a half
    grid.draw(gt2nt) # plot with grid draw
    
    ggexport(ggarrange(gt1,gt1nt,gt2,gt2nt,nrow = 1,ncol = 1),filename = paste0(outname,"_norm.pdf"),height = 0.12*length(tl),width = 12,units = "cm",limitsize = F)    
    
  }
  return(list("tree" = tree,"cm" = cm,"cm2" = cm2,"count" = tree_cell_count,
              "norm" = tree_cell_norm,"td" = td_cell,"edge" = td_wc2))
}


#heatmap normalize,filter,spearman&edu
MyHeatmap = function(tree,outname){
  library(amap)
  tree_cell_countm = tree$count;rownames(tree_cell_countm) = tree_cell_countm$tags;tree_cell_countm = tree_cell_countm[-1]
  tree_cell_countm = as.matrix(tree_cell_countm)
  
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







