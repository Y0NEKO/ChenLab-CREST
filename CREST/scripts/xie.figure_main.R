#The article completion
#21.09.25

setwd("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/")
library(Seurat)
#1. 编辑pattern ribon图 label修改
INDELRibbonPlot = function(del_ranges_df,scarseg,out){
  all_site_stat = function(i,dat){
    return(c(start(dat[i]):end(dat[i])))
  }
  all_site_per_stat = function(i,dat){
    return(length(which(dat==i))/length(dat))
  }
  
  del_allsite = unlist(lapply(1:length(del_ranges_df),all_site_stat,dat=del_ranges_df))
  del_allsite = del_allsite[which(del_allsite >= segref[1,2] & del_allsite <= segref[nrow(segref),3])]
  del_allsite_per = data.frame("Site" = unique(del_allsite),
                               "Freq" = unlist(lapply(unique(del_allsite),all_site_per_stat,dat=del_allsite)))
  siteplus = c(segref[1,2]:segref[7,3])[!c(segref[1,2]:segref[7,3] %in% del_allsite_per$Site)]
  siteplus = data.frame(Site=siteplus,Freq=rep(0,length(siteplus)))
  del_allsite_per = rbind(del_allsite_per,siteplus)
  p = ggplot()+geom_ribbon(data = scarseg,aes(x=scar,ymin=0,ymax=max(del_allsite_per$Freq),fill=type),alpha=0.1) +
    geom_line(data = del_allsite_per,aes(x = Site,y = Freq),size=2,colour="lightskyblue") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks=c()) +
    xlab("")+ylab("Frequency")+theme_bw()
  print(p)
  ggsave(p, filename = out,width = 8,height = 2.5)
  return(p)
}
{
  
  
  
  
}


#2. diversity预测
{
  ExtractScar = function(v1lsx,blst){
    arraydata = NULL
    for (i in 2:length(v1lsx)) {
      arrayline = v1lsx[[i]]
      arrayline = arrayline[which(!arrayline$umim %in% blst$Var1),]
      arrayline = arrayline[which(!arrayline$umim %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      arraydata = c(arraydata,arrayline$umim)
    }
    
    
    statd = NULL
    for (n in seq(1,length(arraydata),100)) {
      #    s = read[sample(1:nrow(read),n),]
      #    pass = s[which(s$quality == "pass"),]
      s = sample(arraydata,n)
      d = length(unique(s))
      line = data.frame("umi" = n,"diversity" = d)
      statd = rbind(statd,line)
    }
    return(statd)
  }
  
  v1fdiv = ExtractScar(v1lsx,bl$v1_array)
  v2fdiv = ExtractScar(v2lsx,bl$v2_array)
  
  
  #模拟v1 + v2
  v3lsxrand = list()
  v3lsxmin = list()
  for (i in 1:length(v1lsx)) {
    line1 = v1lsx[[i]]
    line2 = v2lsx[[i]]
    line1 = line1[which(!line1$umim %in% bl$v1_array$Var1),]
    line1 = line1[which(!line1$umim %in% "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"),]
    line2 = line2[which(!line2$umim %in% bl$v2_array$Var1),]
    line2 = line2[which(!line2$umim %in% "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"),]
    line1$umim = paste0("v1.",line1$umim)
    line2$umim = paste0("v2.",line2$umim)
    if(nrow(line1)>nrow(line2)){linel = line1;lines = line2}else{linel = line2;lines = line1}
    # line1 = line1[which(!line1$umim %in% c(bl$v1_array$Var1,"NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
    # line2 = line2[which(!line2$umim %in% c(bl$v2_array$Var1,"NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
    
    #最少的组合
    indel1 = table(lines$umim);indel1 = indel1[order(indel1)]
    indel2 = table(linel$umim);indel2 = indel2[order(indel2)]
    tmp1 = indel1
    tmp2 = indel2
    cmb = NULL
    for (f1 in 1:length(indel1)) {
      eq1 = indel1[f1]
      eq2 = tmp2[which(tmp2 == eq1)]
      if(length(eq2) != 0){
        eq2 = eq2[1]
        tmp2 = tmp2[which(!names(tmp2) %in% names(eq2))]
        tmp1 = tmp1[which(!names(tmp1) %in% names(eq1))]
        cmb = rbind(cmb,data.frame(v1 = names(eq1),v2 = names(eq2),
                                   v1.count = eq1,v2.count = eq2))
      }
    } 
    
    f2 = 1
    flag = 1
    tmpi = tmp1
    eq1 = NULL
    while (1) {
      eq1 = c(eq1,tmpi[flag])
      eq2 = tmp2[which(tmp2 == sum(eq1))]
      if(length(eq2) != 0){
        eq2 = eq2[1]
        tmp2 = tmp2[which(!names(tmp2) %in% names(eq2))]
        tmp1 = tmp1[which(!names(tmp1) %in% names(eq1))]
        cmb = rbind(cmb,data.frame(v1 = names(eq1),v2 = names(eq2),
                                   v1.count = eq1,v2.count = eq2))
        eq1 = NULL
        f2 = f2 + 1
        flag = f2
      }else{
        flag = flag + 1
        if(flag > length(tmp1)){f2 = f2 + 1; flag = f2; eq1 = NULL}
      }
      if(f2 > length(tmp1)){break}
      
    }
    
    tmpi = tmp1
    f3 = 1
    while(1){
      if(length(tmpi) == 0){break}
      eq1 = tmpi[f3]
      eq2 = tmp2[which(tmp2 > eq1)]
      if(length(eq2) != 0){
        eq2 = eq2[1]
        tmp2[which(names(tmp2) %in% names(eq2))] = tmp2[which(names(tmp2) %in% names(eq2))] - eq1
        tmp1 = tmp1[which(!names(tmp1) %in% names(eq1))]
        cmb = rbind(cmb,data.frame(v1 = names(eq1),v2 = names(eq2),
                                   v1.count = eq1,v2.count = eq2))
        f3 = f3 + 1
      }else{
        if(length(tmp2) == 0){break}
        tmp2 = tmp2[order(tmp2)]
        eq2 = tmp2[length(tmp2)]
        tmp2 = tmp2[which(!names(tmp2) %in% names(eq2))]
        cmb = rbind(cmb,data.frame(v1 = names(eq1),v2 = names(eq2),
                                   v1.count = eq1,v2.count = eq2))
        tmpi[f3] = tmpi[f3] - eq2
      }
      if(f3 > length(tmpi)){break}
    }
    for (j in 1:nrow(cmb)) {
      cmbline = cmb[j,]
      if(cmbline$v1.count < cmbline$v2.count){
        linel[which(linel$umim == cmbline$v2),"umim"][1:cmbline$v1.count] = paste0(cmbline$v1,"-",cmbline$v2)
      }else{
        linel[which(linel$umim == cmbline$v2),"umim"] = paste0(cmbline$v1,"-",cmbline$v2)
      }
      
    }
    v3lsxmin[[i]] = linel
    
    #随机的组合
    if(nrow(line1)>nrow(line2)){linel = line1;lines = line2}else{linel = line2;lines = line1}
    
    cmid = sample(1:nrow(linel),nrow(lines))
    linel$umim[cmid] = paste0(lines$umim,"-",linel$umim[cmid])
    v3lsxrand[[i]] = linel
    names(v3lsxrand)[i] = names(v3lsxmin)[i]= names(v1lsx)[i]
    
  }
  # v3blmin = SelectBlacklist1(v3lsxmin)
  v3fdivmin = ExtractScar(v3lsxmin,v3blmin)

  # v3blrand = SelectBlacklist1(v3lsxrand)
  v3fdivrand = ExtractScar(v3lsxrand,v3blrand)
  
  tail(v3fdivmin)
  tail(v3fdivrand)
  tail(v2fdiv)
  
  
  write.table(v1fdiv, file = "final_result/xBulk/diversity_predict_v1.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v2fdiv, file = "final_result/xBulk/diversity_predict_v2.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v3fdivmin, file = "final_result/xBulk/diversity_predict_v3_min.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v3fdivrand, file = "final_result/xBulk/diversity_predict_v3_rand.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  
  
  v1fun = function(x){(6337*x + 3.842e+06)/(x + 1.485e+04)}
  v2fun = function(x){(1.83e+04*x + 1.91e+07)/(x + 3.494e+04)}
  v3funrand = function(x){(1.774e+05*x + 1.32e+08)/(x + 2.37e+05)}
  v3funmin = function(x){(2.313e+04*x + 2.513e+07)/(x + 3.94e+04)}
  # optimize(v1fun, c(0, 100000000000000000), tol = 0.0001,maximum = T)
  # optimize(v2fun, c(0, 1000000000000000000), tol = 0.0001,maximum = T)
  
  v1fdiv$type = "v1"
  v2fdiv$type = "v2"
  vp = rbind(v1fdiv,v2fdiv)
  options(scipen = 3)
  
  p1 = ggplot(vp,aes(x = umi,y = diversity,color = type)) + xlim(0,500000) + ylim(0,20000) +
    geom_point(size = 0.5) +
    geom_function(aes(colour = "v1"),fun = v1fun) +
    geom_function(aes(colour = "v2"),fun = v2fun) +
    annotate("text", x=100000, y=20000, label= "max(v1) = 6337,max(v2) = 18300") + 
    theme_ipsum(base_size = 18, base_family = "sans") + labs(x = NULL, y = NULL)
  
  p1
  
  v3fdivrand$type = "random"
  v3fdivmin$type = "min"
  vp = rbind(v3fdivrand,v3fdivmin)
  p2 = ggplot(vp,aes(x = umi,y = diversity,color = type)) + xlim(0,5000000) + ylim(0,200000) +
    geom_point(size = 0.5) +
    geom_function(aes(colour = "v1_v2_random"),fun = v3funrand,linetype = "dashed") +
    geom_function(aes(colour = "v1_v2_min"),fun = v3funmin,linetype = "dashed") +
    annotate("text", x=1500000, y=200000, label= "max(random) = 177400,max(min) = 23130") + 
    theme_ipsum(base_size = 18, base_family = "sans") + labs(x = NULL, y = NULL)
  
  p2
  
  ggsave(p1, filename = "final_result/xBulk/Diversity_predict.pdf",width = 8,height = 6)
  ggsave(p2, filename = "final_result/xBulk/Diversity_predict_v1v2_simulation.pdf",width = 10,height = 6)
  
}


#3. bulk top array 图 & 编辑size统计
{
  samplels = list.files("bulk_result/X.bulk_in_vivo/",full.names = T)
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = names(v2lsx) = list.files("bulk_result/X.bulk_in_vivo/")
  bl1 = read.csv("sc_result/dual_E15/v1_indel_blacklist.csv",header = T,row.names = 1)
  colnames(bl1) = "tag"
  bl2 = read.csv("sc_result/dual_E15/v2_indel_blacklist.csv",header = T,row.names = 1)
  colnames(bl2) = "tag"
  bl1 = rbind(bl1,data.frame(tag = "202D+17"))
  #(3)indel for plotting pattern
  #library(stringr)
  
  
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
    
    for (i in 2:length(v1lsx)) {
      
      scardf = as.data.frame(table(v1lsx[[i]]$umim))
      scardf = scardf[which(!scardf$Var1 %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                             "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      scardf = scardf[order(-scardf$Freq),]
      scardf$Freq = scardf$Freq/sum(scardf$Freq)
      
      scardf = scardf[1:num,]
      colnames(scardf)[2] = names(v1lsx)[i]
      if(i == 2){
        scardft = scardf
      }else{
        scardft = merge(scardft, scardf, by = "Var1", all = T)
        
      }
      
    }
    
    
    scardft[is.na(scardft)] = 0
    scardft = scardft[do.call(order, c(scardft[-1], list(decreasing=TRUE))),]
    indel = ScartoIndel(as.character(scardft$Var1))
    index = data.frame(id = scardft$Var1,order = nrow(scardft):1)
    indel = merge(indel,index,by = "id")
    scardfti = melt(scardft)
    colnames(scardfti) = c("id","sample","proportion")
    scardfti = merge(scardfti,index,by = "id")
    
    scardfti$sample = factor(scardfti$sample,levels = colrank)
    p3.0 = ggplot() + geom_segment(data = indel,aes(x = start, xend = end, y = order, yend = order,color=type),
                                   size = 1) +
      theme_void()  + theme(legend.position = "left") + ylab("") + 
      ggplot() + geom_tile(data = scardfti,aes(x = sample,y = order,fill = proportion,height = 1,width = 1),size = 2) + 
      ylab("") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(),axis.title.x = element_blank()) + 
      # scale_fill_viridis(option = "A")
      scale_fill_gradient2(midpoint = 0.001,low = "blue",high = "red")
    return(list("scardf" = scardfti,"indel" = indel,"figure" = p3.0))
  }
  
  
  names(v1lsx) = names(v2lsx) = c("ctrl","E10-1","E10-2","E12","E15-1a","E15-1b",
                                  "E15-2a","E15-2b","E15-3a","E15-3b","E9","P0")
  v1lsx.r = v1lsx[c(1,11,2:10,12)]
  v2lsx.r = v2lsx[c(1,11,2:10,12)]
  colrank = names(v1lsx.r) = names(v2lsx.r) = c("ctrl","E9","E10-1","E10-2","E12","E15-1a","E15-1b",
              "E15-2a","E15-2b","E15-3a","E15-3b","P0")
  scardf1 = ScarPlot(v1lsx.r,20,colrank)
  scardf2 = ScarPlot(v2lsx.r,20,colrank)
  scardf1$figure
  scardf2$figure
  
  # dir.create("final_result/xBulk")
  ggsave(scardf1$figure,filename = "final_result/xBulk/fig2.C.top_array_v1.pdf",width = 10,height = 6)
  ggsave(scardf2$figure,filename = "final_result/xBulk/fig2.D.top_array_v2.pdf",width = 10,height = 6)
  
  #统计每种array在多少个样本中出现
  StatArrayOcc = function(v1lsx){
    arraydft = NULL
    for (i in 2:length(v1lsx)) {
      
      arraydf = as.data.frame(table(v1lsx[[i]]$umim))
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
    arraydft = arraydft[c(-(5:10))]
    
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
    scale_fill_manual(values=c('lightpink1','lightblue2')) +
    theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))
  print(p1)
  ggsave(p1,filename = "final_result/xBulk/fig2.I.barcode_occrurency_in_sample.pdf",width = 7,height = 5)
  write.table(ayocv,file = "final_result/xBulk/barcode_occrurency_in_sample.txt",row.names = F,quote = F,
              sep = "\t")
  
  #Select blacklist
  SelectBlacklist1 = function(v1lsx){
    arraydft = NULL
    for (i in 2:length(v1lsx)) {
      
      arraydf = as.data.frame(table(v1lsx[[i]]$umim))
      arraydf = arraydf[which(!arraydf$umim %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
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
    arraydft$E15_1 = (arraydft$`Dual_E15.5_1a`+ arraydft$`Dual_E15.5_1b`)/2
    arraydft$E15_2 = (arraydft$`Dual_E15.5_2a`+ arraydft$`Dual_E15.5_2b`)/2
    arraydft$E15_3 = (arraydft$`Dual_E15.5_3a`+ arraydft$`Dual_E15.5_3b`)/2
    arraydft = arraydft[c(-(5:10))]
    
    blackls = arraydft[which(rowSums(arraydft[-1]>0.005)>1),]
    blackls = arraydft
    return(blackls)
    
  }
  SelectBlacklist2 = function(v1lsx){
    tagdft = NULL
    for (i in 2:length(v1lsx)) {
      
      tag = NULL
      df = v1lsx[[i]]$umim
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
    tagdft$E15_1 = (tagdft$`E15-1a`+ tagdft$`E15-1b`)/2
    tagdft$E15_2 = (tagdft$`E15-2a`+ tagdft$`E15-2b`)/2
    tagdft$E15_3 = (tagdft$`E15-3a`+ tagdft$`E15-3b`)/2
    tagdft = tagdft[c(-(5:10))]
    head(tagdft)
    
    blackls = tagdft[which(rowSums(tagdft[-1]>0.005)>1),]
    
  }
  
  # v2blay[which(v2blay$Var1 == "108D+17_108D+17_108D+17_108D+17_108D+17_NONE_NONE_NONE"),]
  # tmp = as.data.frame(table(v2test11$result$pattern))
  # tmp[which(tmp$Var1 == "189D+17_189D+17_189D+17_189D+17_189D+17_189D+17_189D+17_189D+17_NONE"),]
  # 
  v1blay = SelectBlacklist1(v1lsx)
  v2blay = SelectBlacklist1(v2lsx)
  v1bltg = SelectBlacklist2(v1lsx)
  v2bltg = SelectBlacklist2(v2lsx)
  
  blt = list("v1_array" = v1blay,"v2_array" = v2blay,
             "v1_tag" = v1bltg,"v2_tag" = v2bltg)
  saveRDS(blt,file = "final_result/x.blacklist.rds")
  blt = readRDS("final_result/x.blacklist.rds")
  FilterBlacklist = function(v1lsx,bl1){
    v1lsx.f = list()
    for (i in 1:length(v1lsx)){
      df = v1lsx[[i]]
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
  
  v1lsx.f = FilterBlacklist(v1lsx,blt$v1_array)
  v2lsx.f = FilterBlacklist(v2lsx,blt$v2_array)
  v1lsx.ftg = FilterBlacklist(v1lsx,blt$v1_tag)
  v2lsx.ftg = FilterBlacklist(v2lsx,blt$v2_tag)
  
  unlist(lapply(v1lsx, nrow))
  unlist(lapply(v1lsx.ftg, nrow))
  unlist(lapply(v1lsx.f, nrow))
  
  unlist(lapply(v2lsx, nrow))
  unlist(lapply(v2lsx.ftg, nrow))
  unlist(lapply(v2lsx.f, nrow))
  
  
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
  ggsave(v1edp,filename = "final_result/xBulk/fig2.E.edit_size_v1.pdf",width = 10,height = 6)
  ggsave(v2edp,filename = "final_result/xBulk/fig2.F.edit_size_v2.pdf",width = 10,height = 6)
  
}


#4. E11，E15 BC,pattern recover rate
#最终过滤在 array_recover脚本中完成
{
  E11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  E15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  
  # saveRDS(E15tree,file = "final_result/xE15.5/E15.5_taq_tree.rds")
  # saveRDS(v1test,file = "final_result/xE15.5/E15.5_taq_v1_result.rds")
  # saveRDS(v2test,file = "final_result/xE15.5/E15.5_taq_v2_result.rds")
  # 
  # saveRDS(E11tree,file = "final_result/xE11/E11_1_tree.rds")
  # saveRDS(v1test,file = "final_result/xE11/E11_v1_result.rds")
  # saveRDS(v2test,file = "final_result/xE11/E11_v2_result.rds")
  
  E15v1 = readRDS("final_result/xE15.5/E15.5_taq_v1_result.rds")
  E15v2 = readRDS("final_result/xE15.5/E15.5_taq_v2_result.rds")
  E11v1 = readRDS("final_result/xE11/E11_v1_result.rds")
  E11v2 = readRDS("final_result/xE11/E11_v2_result.rds")
  cell11 = ReadCell("sc_result/X.E11/E11_celltype_0923.csv")
  cell15 = ReadCell("sc_result/dual_E15/celltype_subclustered0901.csv")
  #统计E11
  v1bc11 = Reduce(intersect,E11v1$rawlist)
  v2bc11 = Reduce(intersect,E11v2$rawlist)
  v1v2bc11 = intersect(v1bc11,v2bc11)
  reBC11 = data.frame("BC" = substr(colnames(E11trans),1,16),
                    "recover" = "unrecover") 
  reBC11[which(reBC11$BC %in% v1bc11),"recover"] = "V1"
  reBC11[which(reBC11$BC %in% v2bc11),"recover"] = "V2"
  reBC11[which(reBC11$BC %in% v1v2bc11),"recover"] = "V1_V2"
  table(reBC11$recover)
  
  
  
  #统计E15
  v1bc15 = Reduce(intersect,E15v1$rawlist[-1])
  v2bc15 = Reduce(intersect,E15v2$rawlist)
  v1v2bc15 = intersect(v1bc15,v2bc15)
  reBC15 = data.frame("BC" = substr(colnames(E15trans),1,16),
                    "recover" = "unrecover") 
  reBC15[which(reBC15$BC %in% v1bc15),"recover"] = "V1"
  reBC15[which(reBC15$BC %in% v2bc15),"recover"] = "V2"
  reBC15[which(reBC15$BC %in% v1v2bc15),"recover"] = "V1_V2"
  table(reBC15$recover)
  
  reBC = rbind(as.data.frame(table(reBC11$recover)),
               as.data.frame(table(reBC15$recover)))
  colnames(reBC) = c("recover","cell.counts")
  reBC$sample = rep(c("E11","E15"),each = 4)
  
  p4.0 = ggplot(reBC[which(reBC$sample == "E11"),],aes(x = sample,y = cell.counts,fill = recover)) + 
    geom_bar(position="fill", stat="identity",width = 0.2) +
    theme_bw()
  p4.0
  
  p4.1 = ggplot(reBC[which(reBC$sample == "E15"),],aes(x = sample,y = cell.counts,fill = recover)) + 
    geom_bar(position="fill", stat="identity",width = 0.2) +
    theme_bw()
  p4.1
  
  
  ggsave(p4.0,filename = "final_result/xE11/E11_recover_rate_unfilter.pdf",width = 4,height = 3)
  ggsave(p4.1,filename = "final_result/xE15.5/E15_recover_rate_unfilter.pdf",width = 4,height = 3)
  
}


#5.  E11, E15 domain 分散度和pseudo time的关系
{
  E11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  E15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  E11tree = readRDS("final_result/xE11/E11_blay.rds")
  E15tree = readRDS("final_result/xE15.5/E15.5_blay_taq_tree.rds")
  
  #Monocle
  {
    library(monocle)
    cds15 = as.CellDataSet(E15trans)
    cds15 <- estimateSizeFactors(cds15)
    cds15 <- estimateDispersions(cds15)
    
    disp_table <- dispersionTable(cds15)
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
    cds15 <- setOrderingFilter(cds15, unsup_clustering_genes$gene_id)
    
    
    cds15 <- reduceDimension(
      cds15,
      max_components = 2,
      method = 'DDRTree')
    # saveRDS(cds15,file = "../final_result/EN1/30&48_monocle_cds15.rds")
    cds15 <- orderCells(cds15)
    
    pdf("final_result/xE15.5/E15.5_pseudotime_sample.pdf")
    plot_cell_trajectory(cds15,color_by = "Pseudotime")
    dev.off()
    
    pseudo15 = data.frame(BC = colnames(cds15),pseudo = cds15$Pseudotime)
    
    pseudo15$BC = substr(pseudo15$BC,1,16)
    pseudo15 = pseudo15[which(pseudo15$pseudo != Inf),]
    
    saveRDS(cds15,file = "final_result/xE15.5/E15.5_monocle_cds15.rds")
    write.csv(pseudo15,file = "final_result/xE15.5/E15.5_monocle_pseudo.csv",row.names = F,quote = F)
    
  }
  
  
  #E11
  {
    cds11 = as.CellDataSet(E11trans)
    cds11 <- estimateSizeFactors(cds11)
    cds11 <- estimateDispersions(cds11)
    
    disp_table <- dispersionTable(cds11)
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
    cds11 <- setOrderingFilter(cds11, unsup_clustering_genes$gene_id)
    
    
    cds11 <- reduceDimension(
      cds11,
      max_components = 2,
      method = 'DDRTree')
    # saveRDS(cds11,file = "../final_result/EN1/30&48_monocle_cds11.rds")
    cds11 <- orderCells(cds11)
    
    pdf("final_result/xE11/E11_pseudotime_sample.pdf")
    plot_cell_trajectory(cds11,color_by = "Pseudotime")
    dev.off()
    
    pseudo11 = data.frame(BC = colnames(cds11),pseudo = cds11$Pseudotime)
    
    pseudo11$BC = substr(pseudo11$BC,1,16)
    pseudo11 = pseudo11[which(pseudo11$pseudo != Inf),]
    
    saveRDS(cds11,file = "final_result/xE11/E11_monocle_cds11.rds")
    write.csv(pseudo11,file = "final_result/xE11/E11_monocle_pseudo.csv",row.names = F,quote = F)
    
  }
  
  cds15 = readRDS(file = "final_result/xE15.5/E15.5_monocle_cds15.rds")
  pseudo15 = read.csv("final_result/xE15.5/E15.5_monocle_pseudo.csv")
  cds11 = readRDS(file = "final_result/xE11/E11_monocle_cds11.rds")
  pseudo11 = read.csv("final_result/xE11/E11_monocle_pseudo.csv")
  
  
  E11psd = PseudoAna(E11tree,pseudo11)
  E15psd = PseudoAna(E15tree,pseudo15)
  
  E11psd$figure.depth + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  E15psd$figure.depth + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(E11psd$figure.depth + theme(axis.text.x = element_text(angle = 45,hjust = 1)),filename = "final_result/xE11/E11_pseudo_level_boxplot.pdf",width = 10, height = 4)
  ggsave(E15psd$figure.depth + theme(axis.text.x = element_text(angle = 45,hjust = 1)),filename = "final_result/xE15.5/E15_pseudo_level_boxplot.pdf",width = 10, height = 4)
  saveRDS(E11psd,"final_result/xE11/E11psd.rds")
  saveRDS(E15psd,"final_result/xE15.5/E15.5psd.rds")
  

  
  #计算分散度和分模块circleplot
  #E15分模块
  domain_anno = read.csv("raw_data/xie/domain annotation.csv",header = T)
  domain_anno$domainid = rep(c(1,2,3,4,5,6,"unknown"),table(domain_anno$domain))
  
  domain_vect = as.vector(domain_anno$domain)
  names(domain_vect) = toupper(domain_anno$ClusterName)
  
  tree = E15tree
  outnameE15 = "final_result/xE15.5/E15.5_taq_domain"
  tree$cm2$celltype = str_replace_all(toupper(tree$cm2$celltype),domain_vect)
  tree = CirclePlot1(tree,outnameE15)
  
  #E15分散度计算
  {
    
    #build root circle
    edge = tree$edge
    
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
    groupdf = CircleGroup(edge)
    
    
    celllist = c("NONE",unique(tree$cm2$celltype))
    
    cm2 = tree$cm2
    cmt = merge(cm2,pseudo15,by.x = "Var1",by.y = "BC")
    cmt = merge(cmt,groupdf,by = "tags")
    root_psd = cmt %>% group_by(group) %>% dplyr::summarise(count = sum(n()),
                                                            pseudo.mean = mean(pseudo),pseudo.sd = sd(pseudo))
    
    #统计分散度
    cmt = merge(cmt,unique(domain_anno[-1]),by.x = "celltype", by.y = "domain")
    root_spr = NULL
    
    for (i in 1:length(unique(cmt$group))) {
      root = unique(cmt$group)[i]
      slice = cmt[which(cmt$group == root),]
      sdm = table(slice$domainid)
      sdm = sdm[order(-sdm)]
      sdm = sdm/sum(sdm)
      H = -sum(sdm*log(sdm))
      
      sdid = names(sdm)
      bias = 0
      if(length(sdid)>1){
        if("unknown" %in% sdid){
          if(length(sdid)>2){
            bias = sd(as.numeric(sdid[which(sdid!="unknown")]))
            bias = bias* (length(sdid)/(length(sdid)-1))
          }else{
            bias = 1
          }
          
        }else{
          bias = sd(as.numeric(sdid))
        }
      }
      Hb = H + bias
      mainModel = names(table(slice$celltype))[which.max(table(slice$celltype))]
      root_spr = rbind(root_spr,data.frame(group = root,Hb = Hb,main.domain = mainModel,
                                           stringsAsFactors = F))
      
    
  }
  
  root_stat = merge(root_psd,root_spr,by = "group")
  
  p5.1 = ggplot(root_stat,aes(x = pseudo.mean, y = Hb,color = main.domain,alpha = 0.7)) + geom_point(aes(size = count)) + 
    scale_size_continuous(range = c(0, 8)) +
    theme_bw()
  p5.1
  
  p5.2 = ggplot(root_stat,aes(x = log(count), y = Hb,color = main.domain,alpha = 0.7)) + geom_point() + 
    theme_bw()
  p5.2
  
  ggboxplot(root_stat$Hb)
  ecdf(root_stat$Hb)
  
  ggplot(root_stat, aes(Hb)) + stat_ecdf(geom = "step")
  ggdensity(root_stat$Hb)
  
  ggsave(p5.1,filename = "final_result/xE15.5/E15_dissemination_pseudo_root.pdf",width = 7,height = 5)
  
  
  }
  
  
  #E15 pseudo time 与 clone作图尝试
  {
    #不使用pseudo time，而是clone成分
    cm = E15tree$cm
    cmtn = merge(cm,pseudo15,by.x = "Var1",by.y = "BC")
    cmtnt = cmtn
    # cmtn[which(cmtnt$celltype == "Rgl1"),"celltype"] = "Rgl3"
    # cmtn[which(cmtnt$celltype == "Rgl3"),"celltype"] = "Rgl1"
    cmtn$domain = str_replace_all(toupper(cmtn$celltype),domain_vect)
    
    head(cmt)
    
    #细胞基本情况
    p1.1 = ggplot(cmt,aes(x = celltype,y = pseudo)) + geom_boxplot(width = 0.7) + theme_bw()
    p1.2 = ggplot(cmt,aes(x = celltype)) + geom_bar(width = 0.7) + theme_bw()
    
    p1.1 + p1.2
    ggsave(p1.1 + p1.2,filename = "final_result/xE15.5/E15_pseudo_cell_box_bar_plot.pdf",width = 8,height = 4)
    
    
    #将pseudo time分为不同阶段，做一个桑基图
    p2.0 = ggplot(cmtn, aes(pseudo)) + stat_ecdf(geom = "step") + geom_vline(xintercept=c(25,35),lty=3,col="black",lwd=0.5)
    p2.0
    ggsave(p2.0,filename = "final_result/xE15.5/E15_pseudo_ecdf.pdf",width = 5,height = 4)
    # cmtn$phase = "0-25"
    # cmtn[which(cmtn$pseudo>25 & cmtn$pseudo <= 35),"phase"] = "25-35"
    # cmtn[which(cmtn$pseudo > 35),"phase"] = "35-"
    
    cmtn$phase = "0-30"
    cmtn[which(cmtn$pseudo > 30),"phase"] = "30-"
    
    head(cmtn)
    
    p2.1 = ggplot(cmtn,aes(x = as.character(phase),fill = domain)) + 
      geom_bar(stat = "count",position = "fill",width = 0.5) + 
      xlab("pseudo time") + theme_bw()
    p2.1
    ggsave(p2.1,filename = "final_result/xE15.5/E15_pseudo_phase_module_bar.pdf",width = 5,height = 5)
    
    Statalluv = function(cmtn){
      
      cmtnp = cmtn %>% group_by(tags,phase,domain) %>% summarise(counts = sum(n()))
      cmtnp = dcast(cmtn,domain~phase)
      cmtnp$order = as.numeric(cmtnp$tags)
      moduleid = c("m1","m2-4","m4","m5","m6","m7","unknown")
      
      
      #统计所有可能性，暂时无视clone size
      cmtnpmx = data.frame("phase1" = rep(moduleid,each = 7),"phase2" = rep(moduleid,7), "Freq" = 0)
      clone = unique(cmtnp$tags)
      for (i in 1:length(clone)) {
        line = cmtnp[which(cmtnp$tags == clone[i]),]
        if(length(unique(line$phase)) > 1){
          line.from = line[which(line$phase == "0-30"),]
          line.to = line[which(line$phase == "30-"),]
          cmtnpmx[which(cmtnpmx$phase1 %in% line.from$domain & 
                          cmtnpmx$phase2 %in% line.to$domain),"Freq"] = cmtnpmx[which(cmtnpmx$phase1 %in% line.from$domain & 
                                                                                        cmtnpmx$phase2 %in% line.to$domain),"Freq"] + 1
        }
      }
      p2.2 = ggplot(data = cmtnpmx,
             aes(axis1 = phase1, axis2 = phase2, 
                 y = Freq)) +
        scale_x_continuous(breaks = 1:2, labels = c("0-30", "30-")) +
        geom_alluvium(aes(fill = phase1)) +
        geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_minimal() + theme_void()
      p2.2
      ggsave(p2.2,filename = "final_result/xE15.5/E15_pseudo_phase_module_stratum_1.pdf",width = 8,height = 5)
      
      
    }
    
    
    clonemd = as.data.frame(table(cmtn[c("tags","domain")]))
    clonemd = dcast(clonemd,tags~domain,value.var = "Freq")
    clonemd.norm = apply(clonemd, 1, ...)
    rownames(clonemd) = clonemd$tags;clonemd = clonemd[-1]
    Heatmap(as.matrix(clonemd),show_row_names = F)
    # p1 = Heatmap(as.matrix(clonemd))
    # print(p1)
    
    #根据pseudotime对celltype分类
    cellps = cmtn %>% group_by(celltype) %>% summarise(pseudo.mean = mean(pseudo),domain = unique(domain))
    border = cellps[which(cellps$celltype == "Nb_M"),]$pseudo.mean
    earlyid = cellps[which(cellps$pseudo.mean <= border),]$celltype
    earlyid = c("NProg_ML_m6","NProg_Foxa2_m5","NProg_Six3","Nb_M",
                "NProg_L","Rgl1","Rgl2_m4","Rgl2_Pax7_m1","Rgl3")
    lateid = setdiff(cellps$celltype,earlyid)
    cmtn$phase = "early"
    cmtn[which(cmtn$celltype %in% lateid),"phase"] = "late"
    
    cmtnp2 = cmtn %>% group_by(tags,phase,celltype) %>% summarise(counts = sum(n()))
      
    #统计所有可能性，暂时无视clone size
    {
      cmtnp2mx = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                            "phase2" = rep(lateid,length(earlyid)), "Freq" = 0)
      clone = unique(cmtnp2$tags)
      for (i in 1:length(clone)) {
        line = cmtnp2[which(cmtnp2$tags == clone[i]),]
        if(length(unique(line$phase)) > 1){
          line.from = line[which(line$phase == "early"),]
          line.to = line[which(line$phase == "late"),]
          cmtnp2mx[which(cmtnp2mx$phase1 %in% line.from$celltype & 
                           cmtnp2mx$phase2 %in% line.to$celltype),"Freq"] = cmtnp2mx[which(cmtnp2mx$phase1 %in% line.from$celltype & 
                                                                                             cmtnp2mx$phase2 %in% line.to$celltype),"Freq"] + 1
        }
      }
      
      is_alluvia_form(cmtnp2mx,weight = "Freq")
      domain_df = unique(cmtn[c("celltype","domain")])
      domain_df = domain_df[order(domain_df$domain),]
      cmtnp2mx_lod = to_lodes_form(cmtnp2mx[which(cmtnp2mx$phase1 != "Rgl1"),],
                                   key = "phase", value = "celltype", id = "alluvium",
                                   axes = 1:2)
      
      cmtnp2mx_lod = merge(cmtnp2mx_lod,domain_df,
                           by = "celltype")
      cmtnp2mx_lod = cmtnp2mx_lod[order(cmtnp2mx_lod$phase,cmtnp2mx_lod$domain),]
      head(cmtnp2mx_lod)
      str(cmtnp2mx_lod)
      cmtnp2mx_lod$celltype = factor(cmtnp2mx_lod$celltype,
                                     levels = unique(as.character(cmtnp2mx_lod$celltype)))
      
      p2.3 = ggplot(cmtnp2mx_lod,
                    aes(x = phase, stratum = celltype, alluvium = alluvium,
                        y = Freq,
                        fill = domain, label = celltype)) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_flow(aes(color = celltype)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      p2.3
      ggsave(p2.3,filename = "final_result/xE15.5/E15_pseudo_phase_module_stratum_celltype_arti_noRgl1.pdf",width = 8,height = 5)
      
    }
      
    #尝试计算每个BC的富集值，利用富集值来做桑基图
    {
      cmtnp3 = cmtn %>% group_by(tags,phase,celltype) %>% 
        summarise(counts = sum(n()),domain = unique(domain),pseudo.mean = mean(pseudo))
      head(cmtnp3)
      
      clone = unique(cmtnp3$tags)
      
      hpsdf = NULL
      for (i in 1:length(clone)) {
        line = cmtnp3[which(cmtnp3$tags == clone[i]),]
        sdm = line$counts
        sdm = sdm[order(-sdm)]
        sdm = sdm/sum(sdm)
        H = -sum(sdm*log(sdm))
        # H = -sum(sdm*log(sdm))/log2(length(sdm))
        psm = mean(line$pseudo.mean)
        psd = sd(line$pseudo.mean)
        size = sum(line$counts)
        celltypenum = nrow(line)
        linedf = data.frame("H" = H,"pseudo.mean" = psm,"pseudo.sd" = psd,"size" = size,
                            "celltypenum" = celltypenum)
        hpsdf =rbind(hpsdf,linedf)
      }
      hpsdf[which(is.nan(hpsdf$H)),"H"] = 0
      hpsdf[which(is.na(hpsdf$pseudo.sd)),"pseudo.sd"] = 0
      head(hpsdf)
      
      ggplot(hpsdf,aes(x = pseudo.mean, y = celltypenum,color = celltypenum,alpha = 0.7)) + 
        geom_point(aes(size = size)) + 
        scale_size_continuous(range = c(0, 8)) +
        theme_bw()
      ggplot(hpsdf,aes(x = celltypenum, y = H,alpha = 0.7)) + 
        geom_point(aes(size = size)) + 
        scale_size_continuous(range = c(0, 8)) +
        theme_bw()
      
      
      #模拟分布，设定一个p值                  
      cellmx.e = dcast(cmtnp3[which(cmtnp3$phase == "early"),c("tags","celltype","counts")],
                       tags~celltype,value.var = "counts",fun.aggregate = sum)
      cellmx.l = dcast(cmtnp3[which(cmtnp3$phase == "late"),c("tags","celltype","counts")],
                       tags~celltype,value.var = "counts",fun.aggregate = sum)
      Me = cmtnp3[which(cmtnp3$phase == "early"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Ne = sum(Me$counts)
      
      Ml = cmtnp3[which(cmtnp3$phase == "late"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Nl = sum(Ml$counts)
      
      cellft.e = NULL
      for (i in 1:nrow(cellmx.e)) {
        line = cellmx.e[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Me[which(Me$celltype == names(k)),]$counts
            N = Ne
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.e = rbind(cellft.e,ftl)
      }
      
      
      cellft.l = NULL
      for (i in 1:nrow(cellmx.l)) {
        line = cellmx.l[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Ml[which(Ml$celltype == names(k)),]$counts
            N = Nl
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.l = rbind(cellft.l,ftl)
      }
      cellft = merge(cellft.e,cellft.l,by = "tags")
      
      #计算cellft
      phasest = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                            "phase2" = rep(lateid,length(earlyid)), "Freq" = 0)
      for (i in 1:nrow(cellft)) {
        line = cellft[i,]
        line.e = unlist(line[,which(colnames(line) %in% earlyid)])
        line.e = line.e[which(line.e>0)]
        line.l = unlist(line[,which(colnames(line) %in% lateid)])
        line.l = line.l[which(line.l>0)]
        for (j in 1:length(line.e)) {
          for (k in 1:length(line.l)) {
            seg = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq"]
            phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq"] = seg + mean(c(line.e[j],line.l[k]))
          }
        }
      }
      
      domain_df = unique(cmtn[c("celltype","domain")])
      domain_df = domain_df[order(domain_df$domain),]
      
      phasestf = phasest[which(phasest$phase1 != "Rgl1"),]
      phasestf = phasest
      gghistogram(phasestf$Freq) + scale_x_continuous(breaks = seq(0,200,10))
      #为什么用20做阈值?
      {
        #尝试做随机分布
        cellmxr.e = NULL
        for (i in 1:nrow(cellmx.e)) {
          line = cellmx.e[i,]
          new = NULL
          res = sum(line[-1])
          for (ie in 2:length(line)) {
            if(res >= 1){
              new = c(new,floor(runif(1,0,res)))
              res = n - sum(new)
            }else{
              new = c(new,0)
            }
          }
          line[sample(2:length(line))] = new
          cellmxr.e = rbind(cellmxr.e,line)
        }
        
        cellmxr.l = NULL
        for (i in 1:nrow(cellmx.l)) {
          line = cellmx.l[i,]
          new = NULL
          res = sum(line[-1])
          for (ie in 2:length(line)) {
            if(res >= 1){
              new = c(new,floor(runif(1,0,res)))
              res = n - sum(new)
            }else{
              new = c(new,0)
            }
          }
          line[sample(2:length(line))] = new
          cellmxr.l = rbind(cellmxr.l,line)
        }
        
        Mer = colSums(cellmxr.e[-1])
        Ner = sum(Mer)
        
        Mlr = colSums(cellmxr.l[-1])
        Nlr = sum(Mlr)
        cellftr.e = NULL
        for (i in 1:nrow(cellmxr.e)) {
          line = cellmxr.e[i,]
          n = sum(line[2:ncol(line)])
          ftl = line
          for (j in 2:ncol(line)) {
            k = unlist(line[j])
            if(k>0){
              M = Mer[which(names(Mer) == names(k))]
              N = Ner
              d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
              row.names(d) <- c("In_category", "not_in_category")
              ft = fisher.test(d)
              ftl[j] = -log10(ft$p.value)
            }
          }
          cellftr.e = rbind(cellftr.e,ftl)
        }
        
        
        cellftr.l = NULL
        for (i in 1:nrow(cellmxr.l)) {
          line = cellmxr.l[i,]
          n = sum(line[2:ncol(line)])
          ftl = line
          for (j in 2:ncol(line)) {
            k = unlist(line[j])
            if(k>0){
              M = Mlr[which(names(Mlr) == names(k))]
              N = Nlr
              d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
              row.names(d) <- c("In_category", "not_in_category")
              ft = fisher.test(d)
              ftl[j] = -log10(ft$p.value)
            }
          }
          cellftr.l = rbind(cellftr.l,ftl)
        }
        
        cellftr = merge(cellftr.e,cellftr.l,by = "tags")
        
        #计算cellft
        phasestr = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                             "phase2" = rep(lateid,length(earlyid)), "Freq" = 0)
        for (i in 1:nrow(cellftr)) {
          line = cellftr[i,]
          line.e = unlist(line[,which(colnames(line) %in% earlyid)])
          line.e = line.e[which(line.e>0)]
          line.l = unlist(line[,which(colnames(line) %in% lateid)])
          line.l = line.l[which(line.l>0)]
          for (j in 1:length(line.e)) {
            for (k in 1:length(line.l)) {
              seg = phasestr[which(phasestr$phase1 == names(line.e)[j] & phasestr$phase2 == names(line.l)[k]),"Freq"]
              phasestr[which(phasestr$phase1 == names(line.e)[j] & phasestr$phase2 == names(line.l)[k]),"Freq"] = seg + mean(c(line.e[j],line.l[k]))
            }
          }
        }
        
        
        
        sum(phasest$Freq)
        sum(phasestr$Freq)
        
        
      }
      
      phasestf = phasestf[which(phasestf$Freq>20),]
      library(ggalluvial)
      phasest_lod = to_lodes_form(phasestf,
                                   key = "phase", value = "celltype", id = "alluvium",
                                   axes = 1:2)
      
      phasest_lod = merge(phasest_lod,domain_df,
                           by = "celltype")
      phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domain),]
      head(phasest_lod)
      phasest_lod$celltype = factor(phasest_lod$celltype,
                                     levels = unique(as.character(phasest_lod$celltype)))
      
      
      p2.4 = ggplot(phasest_lod,
                    aes(x = phase, stratum = celltype, alluvium = alluvium,
                        y = Freq,
                        fill = domain, label = celltype)) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_flow(aes(color = celltype)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      p2.4
      ggsave(p2.4,filename = "final_result/xE15.5/E15_pseudo_phase_module_stratum_celltype_arti_fisher_noRgl1_filter20.pdf",width = 8,height = 5)
      
      
    }
    
    
    # write.csv(cmtn,"final_result/xE15.5/cell_bc_pseudo_domain_alllist.csv",row.names = F,quote = F)
    # write.csv(cmtnp3,"final_result/xE15.5/cell_tags_pseudo_domain_alllist.csv",row.names = F,quote = F)
    # write.csv(cellft,"final_result/xE15.5/cell_fisher_enrich_matrix_phase.csv",row.names = F,quote = F)
    cmtn = read.csv("final_result/xE15.5/cell_bc_pseudo_domain_alllist.csv")
    
    cmtnp3 = read.csv("final_result/xE15.5/cell_tags_pseudo_domain_alllist.csv")
    cellft = read.csv("final_result/xE15.5/cell_fisher_enrich_matrix_phase.csv")
    #细胞聚类
    {
      #算一个total BC 富集值
      cmtnp2 = E15tree$cm %>% group_by(tags,celltype) %>% 
        summarise(counts = sum(n()))
      #算一个total BC 富集值
      cellmx2 = dcast(cmtnp2[,c("tags","celltype","counts")],
                     tags~celltype,value.var = "counts",fun.aggregate = sum)
      Mt = cmtnp2 %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Nt = sum(Mt$counts)
      # cellmx2 = dcast(cmtnp3[,c("tags","celltype","counts")],
      #                  tags~celltype,value.var = "counts",fun.aggregate = sum)
      # Mt = cmtnp3 %>% group_by(celltype) %>% summarise(counts = sum(counts))
      # Nt = sum(Mt$counts)
      
      
      cellft2 = NULL
      for (i in 1:nrow(cellmx2)) {
        line = cellmx2[i,]
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
        cellft2 = rbind(cellft2,ftl)
      }
      cellft2l = melt(cellft2)
      head(cellft2l)
      cellft2l = merge(cmtnp3,cellft2l,by.x = c("tags","celltype"),by.y = c("tags","variable"))
      write.csv(cellft2,file = "final_result/xE15.5/cell_fisher_enrich_matrix.csv",row.names = F,quote = F)
      
      
      #尝试聚类分析
      # data(nutrient,package = "flexclust")
      # head(nutrient)
      
      rownames(cellft2) = cellft2$tags
      cellft2 = cellft2[-1]
      cellft2n = scale(cellft2)
      
      #测试聚类数目
      {
        
        library(cluster) 
        library(factoextra)
        
        pdf("final_result/xE15.5/Kmeans_try.pdf")
        #1
        library(mclust)
        m_clust <- Mclust(as.matrix(cellft2n), G=1:20) #聚类数目从1一直试到20
        summary(m_clust)
        plot(m_clust, "BIC")
        
        #4
        # library(NbClust)
        # set.seed(1234) #因为method选择的是kmeans，所以如果不设定种子，每次跑得结果可能不同
        # nb_clust <- NbClust(cellft2n,  distance = "euclidean",
        #                     min.nc=2, max.nc=15, method = "kmeans",
        #                     index = "alllong", alphaBeale = 0.1)
        # barplot(table(nb_clust$Best.nc[1,]),xlab = "聚类数",ylab = "支持指标数")
        
        
        #2
        library(factoextra)
        library(ggplot2)
        set.seed(1234)
        fviz_nbclust(cellft2n, kmeans, method = "wss") +
          geom_vline(xintercept = 3, linetype = 2)
        
        
        #3
        library(fpc)
        pamk.best <- pamk(cellft2n)
        pamk.best$nc
        library(cluster)
        clusplot(pam(cellft2n, pamk.best$nc))
        
        
        #4
        library(apcluster)
        ap_clust <- apcluster(negDistMat(r=3), cellft2n)
        length(ap_clust@clusters)
        heatmap(ap_clust)
        
        #5
        library(factoextra)
        fviz_nbclust(cellft2n, kmeans, method = "silhouette")
        
        #6
        set.seed(123)
        gap_clust <- clusGap(cellft2n, kmeans, 10, B = 500, verbose = interactive())
        gap_clust
        fviz_gap_stat(gap_clust)
        
        #7
        h_dist <- dist(as.matrix(cellft2n))
        h_clust<-hclust(h_dist)
        plot(h_clust, hang = -1, labels = FALSE)
        rect.hclust(h_clust,6)
        
        dev.off()
        
        
      }
      
      library(ggfortify)
      library(ggplot2)
      fit = kmeans(cellft2n,6)
      fit
      barplot(t(fit$centers),beside = T,xlab = "cluster",ylab="value")
      plot(cellft2n,col=fit$cluster)
      autoplot(kmeans(cellft2n, 6),data=cellft2n,label=TRUE, label.size=0.5, frame=TRUE) + theme_bw()
      tagclu = fit$cluster
      tagclu = data.frame("tags" = names(tagclu),"cluster" = tagclu)
      cellft2l = merge(cellft2l,tagclu,by = "tags")
      
      ggplot(cellft2l,aes(x = celltype,y = value)) + geom_boxplot() +
        facet_wrap(~cluster)
      
      library(ComplexHeatmap)
      p3.1 = Heatmap(as.matrix(cellft2), name = "cell_enrich", row_km = 10,show_row_names = F)
      p3.1
      ggexport(p3.1, filename = "final_result/xE15.5/tag_cluster_heatmap.pdf",width = 8,height = 8)
      
      # library(GGally)
      tmp = table(cellft2l$tags)
      tmp = tmp[which(tmp>1)]
      cellft2f = cellft2[which(rownames(cellft2)%in% names(tmp)),]
      cellmx2f = cellmx2[-1];rownames(cellmx2f) = cellmx2$tags;cellmx2f = cellmx2f[which(rownames(cellmx2f)%in% names(tmp)),]
      
      library(ggfortify)
      library(ggplot2)
      
      fit = kmeans(cellft2f,17)
      fit
      barplot(t(fit$centers),beside = T,xlab = "cluster",ylab="value")
      plot(cellft2,col=fit$cluster)
      autoplot(kmeans(cellft2, 9,algorithm = "MacQueen"),data=cellft2n,label=TRUE, label.size=0.5, frame=TRUE) + theme_bw()
      tagclu = fit$cluster
      tagclu = data.frame("tags" = names(tagclu),"cluster" = tagclu)
      cellft2l = merge(cellft2l,tagclu,by = "tags")
      
      colnames(cellft2f)

      cellorder = c("Glu4","Glu3","Glu1","Rgl2_Pax7_m1","NProg_L","Gaba2_m1","Glu2",
                    "Rgl3","NProg_Six3","Gaba4","Gaba5","Gaba1_m4","Gaba3","Rgl2_m4","Gaba6_Foxa2_m5","NProg_Foxa2_m5",
                    "NProg_ML_m6","Nb_RN",
                    "Rgl1","Nb_Glu","Glu_RN","Glu_M","Nb_M","DA","Nb_DA")
      corder15 = E15tree$distfigure$abspearsonhc$labels[E15tree$distfigure$abspearsonhc$order]
      cellorder = corder15
      
      # split <- factor(as.character(fit$cluster), levels=as.character(c(3,7,11,2,1,5,10,8,6,9,4)))
      n = 30
      
      # cellft2fb = as.matrix(cellft2f);cellft2fb[which(cellft2fb>0)] = 1
      fit = kmeans(cellft2f,n)
      split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
      reorder.hmap <- Heatmap(as.matrix(cellft2fb), split=split, cluster_row_slices = FALSE,
                              show_row_names = F,column_order = cellorder)
      reorder.hmap
      
      ggexport(reorder.hmap, filename = "final_result/xE15.5/tag_cluster_heatmap.pdf",width = 8,height = 8)
      
      cellmx2fn = t(t(as.matrix(cellmx2f))/colSums(cellmx2f))
      n = 30
      fit = kmeans(cellmx2fn,n)
      
      
      #merge cluster
      
      clist = list(c1 = c(1,3,13,17,26),
                   c2 = c(21,23),
                   c3 = c(10,11,12,26),
                   c4 = c(4,24,9),
                   c5 = c(14,15,16,20),
                   c6 = c(6,5,18,27,28),
                   c7 = c(25,19,30),
                   c8 = c(7,8),
                   c9 = c(2),
                   c10 = c(29),
                   c11 = c(22))
      split <- fit$cluster
      for (i in 1:length(clist)) {
        split[which(split %in% clist[[i]])] = paste0("c",i)
      }
      
      split <- factor(as.character(split), levels=as.character(c(paste0("c",1:length(clist)) )))      
      names(split) = names(fit$cluster)
      
      e15_clu = readRDS("final_result/xE15.5/cellcounts_cluster.rds")
      reorder.hmap <- Heatmap(as.matrix(e15_clu$cellcounts), split=e15_clu$cluster, cluster_row_slices = FALSE,
                              col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                              show_row_names = F,column_order = cellorder)
      reorder.hmap
      
      ggexport(reorder.hmap, filename = "final_result/xE15.5/tag_cluster_heatmap_counts_arti.pdf",width = 8,height = 8)
      
      saveRDS(list(cluster = split,fit = fit,cellcounts = cellmx2fn),file = "final_result/xE15.5/cellcounts_cluster.rds")
      
      #cluster show in dimplot
      ccls = readRDS("final_result/xE15.5/cellcounts_cluster.rds")
      
      ccidents = as.data.frame(ccls$cluster)
      ccidents = merge(ccidents,E15tree$cm,by.x = "row.names",by.y = "tags")
      colnames(ccidents) = c("tags","cluster","CellBC","Freq","Celltype")
      ccidents$Celltype = factor(ccidents$Celltype,levels = cellorder)
      ccidents[which(ccidents$CellBC %in% E15tree$cm[E15tree$cm$celltype == "Rgl1","Var1"]),"Celltype"] = "Rgl3"
      ccidents[which(ccidents$CellBC %in% E15tree$cm[E15tree$cm$celltype == "Rgl3","Var1"]),"Celltype"] = "Rgl1"
      write.csv(ccidents,file = "final_result/xE15.5/cluster_BC.csv",row.names = F,quote = F)
      
      p3.2 = ggplot(ccidents,aes(x = cluster,fill = Celltype)) +
        geom_bar(stat = "count",position = "fill",width = 0.7) +
        theme_bw()
      p3.2
      ggsave(p3.2,filename = "final_result/xE15.5/tag_cluster_celltype_bar.pdf",width = 10,height = 6)
      
      p3.3 = list()
      for (i in 1:11) {
        cc = ccidents[which(ccidents$cluster == paste0("c",i)),]
        tarcell = colnames(E15trans)[which(substr(colnames(E15trans),1,16) %in% cc$CellBC)]
        pt = DimPlot(E15trans,label = TRUE, label.size = 3,cells.highlight = tarcell, 
                cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
        p3.3[[i]] = pt
        names(p3.3)[i] = paste0("c",i)
        
      }
      p3.3.0 = ggarrange(plotlist =  p3.3,labels = names(p3.3))
      p3.3
      
      pdf("final_result/xE15.5/tag_cluster_umap.pdf",width = 10,height = 8)
      print(p3.3.0)
      invisible(lapply(p3.3, print))
      dev.off()
      
    }
    
    #pseudo time图
    {
      
      cmtn = cmtn[order(cmtn$domain),]
      cmtn$celltype = factor(as.character(cmtn$celltype),levels = unique(cmtn$celltype))
      p4.1 = ggboxplot(cmtn[which(cmtn$domain != "unknown"),],x = "celltype",y = "pseudo", color = "domain") + theme_bw() + 
        theme(axis.text.x =  element_text(angle = 45,vjust = 0.5))
      p4.1
      
      colnames(cellft2l)[7] = "enrich.degree"
      head(cellft2l)
      cellft2l = cellft2l[which(cellft2l$domain != "unknown"),]
      p4.2 = ggplot(cellft2l,aes(x = pseudo.mean,y = enrich.degree,
                                 size = counts,color = domain)) + 
        geom_point(alpha = 0.5) + theme_bw()
      p4.2
      
      
      p4.3 = ggplot(cellft2l,aes(x = pseudo.mean,y = enrich.degree,
                                 size = counts,color = celltype)) + 
        geom_point(alpha = 0.5) + theme_bw()
      p4.3
      
      p4.4 = ggplot(cellft2l,aes(x = enrich.degree,fill = domain)) +
        geom_density() + scale_fill_viridis(discrete=TRUE) + facet_wrap(~domain) + theme_ipsum(base_family = "sans")
      
      p4.4
      
      p4.1 + p4.2 + p4.3
      
      #计算array的分散度
      {
        # cmtnf = cmtn[which(cmtn$domain != "unknown"),]
        cmtnf = cmtn[which(cmtn$domain != "unknown"),]
        cmtnf$celltype = factor(cmtnf$celltype,levels = unique(cmtnf$celltype))
        cmtnf = cmtnf[order(cmtnf$domain),]
        cmtnf$domainid = factor(as.numeric(factor(cmtnf$domain,levels = unique(cmtnf$domain))))
        
        root_spr  =NULL
        root_spr_mod  =NULL
        for (i in 1:length(unique(cmtnf$tags))) {
          root = unique(cmtnf$tags)[i]
          slice = cmtnf[which(cmtnf$tags == root),]
          pseudom = mean(slice$pseudo)
          
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
            line_m = data.frame(group = root,H = H_m,Hwi = Hwi_m,Hb = Hb,
                              main.domain = mainModel_m,
                              pseudo.mean = pseudom,
                              count = nrow(slice),
                              stringsAsFactors = F)
            cc_m = as.data.frame.list(cc_m)
            line_m = cbind(line_m,cc_m)
            root_spr_mod = rbind(root_spr_mod,line_m)
            
            mainModel = names(table(slice$celltype))[which.max(table(slice$celltype))]
            meanModel = mean(as.numeric(slice$domainid))
            line = data.frame(group = root,H = H_m,Hwi = Hwi_m,Hb = Hb,
                              main.cell = mainModel,
                              mean.model = meanModel,
                              pseudo.mean = pseudom,
                              count = nrow(slice),
                              stringsAsFactors = F)
            cc = table(slice$celltype)
            cc = as.data.frame.list(cc)
            line = cbind(line,cc)
            root_spr = rbind(root_spr,line)
          }
          
        }
        head(root_spr)
        head(root_spr_mod)
        colnames(root_spr_mod)[-(1:7)] = c("m1","m2-4","m4","m5","m6","m7")
        
        library(scatterpie)
        
        root_sprn = root_spr
        root_sprn[-(1:8)] = t(t(as.matrix(root_sprn[-c(1:8)]))/colSums(root_sprn[-(1:8)]))
        
        tem2 <- as.data.frame(spline(root_sprn$mean.model, root_sprn$Hb, n=100))
        p5.1 = ggplot() + 
          geom_scatterpie(data = root_sprn,aes(x=mean.model, y=Hb, group=group,
                                              r = log10(count)/15,alpha = 0.5),color=NA,
                                          cols=colnames(root_sprn)[-(1:8)]) + coord_fixed() +
          geom_smooth(data = tem2, aes(x=x, y=y), se = F, method = 'loess',color = "black") +
          scale_fill_viridis(discrete=TRUE) + scale_y_continuous(breaks = c(0,1,2,3)) + 
          xlab("mean.model") + ylab("Hb") +
          theme_bw()
        p5.1
        
        
        root_spr_modn = root_spr_mod
        root_spr_modn[-(1:7)] = t(t(as.matrix(root_spr_modn[-c(1:7)]))/colSums(root_spr_modn[-(1:7)]))
        
        # tem3 <- as.data.frame(spline(root_spr_modn$pseudo.mean, root_spr_modn$Hb*10, n=100))
        p5.2 = ggplot() + 
          geom_scatterpie(data = root_spr_modn[-14],aes(x=pseudo.mean, y=Hb*10, group=group,
                                              r = log10(count)/2.5,alpha = 0.5),color=NA,
                          cols=colnames(root_spr_modn)[-c(1:7,14)]) + coord_fixed() + 
          scale_fill_viridis(discrete=TRUE) +
          theme_bw()
        p5.2
        
        
        #5.2延伸
        root_spr_modn$max = 0
        for (i in 1:nrow(root_spr_modn)) {
          line = root_spr_modn[i,]
          root_spr_modn[i,"main.domain"] = names(which.max(line[-(1:7)]))
          root_spr_modn$max[i] = max(line[-(1:7)]/sum(line[-(1:7)]))
        }
        
        p5.2.1 = ggplot(root_spr_modn,aes(x = pseudo.mean,fill = main.domain)) +
          geom_density() + scale_fill_viridis(discrete=TRUE) + facet_wrap(~main.domain) + theme_ipsum(base_family = "sans")
        p5.2.1
        
        p5.2.2 = ggplot(root_spr_modn,aes(x = Hb,fill = main.domain)) +
          geom_density() + scale_fill_viridis(discrete=TRUE) + facet_wrap(~main.domain) + theme_ipsum(base_family = "sans")

        p5.2.2
        
        
        p5.2.3 = ggplot(root_spr_modn,aes(x = Hb, y = max,color = main.domain)) +
          geom_point() + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans") + 
          ylab("max_percentage")
        
        p5.2.3
        
        
        ggsave(p5.2.5,filename = "final_result/xE15.5/pseudo_Hindex_domain_Hb_boxplot.pdf",width = 8,height = 6)
        
        p5.2.4 = ggplot(root_spr_modn,aes(x = pseudo.mean, y = max,color = main.domain)) +
          geom_point() + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans") + 
          ylab("max_percentage")
        
        p5.2.4
        
        p5.2.5 = ggplot(root_spr_modn,aes(x = main.domain, y = Hb)) +
          geom_boxplot(width = 0.6) + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans")
        p5.2.5
        
        p5.2.6 = ggplot(root_spr_modn,aes(x = main.domain, y = H)) +
          geom_boxplot(width = 0.6) + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans")
        p5.2.6
        
        p5.2.7 = ggplot(root_spr_modn,aes(x = log2(count), y = Hb,color = main.domain)) +
          geom_point() + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans")
        p5.2.7
        cor(log2(root_spr_modn$count),root_spr_modn$Hb)
        
      }
      
      
      pdf("final_result/xE15.5/pseudo_Hindex_all.pdf",width = 10,height = 6)
      invisible(lapply(list(p4.1,p4.2,p4.3,p4.4,p5.1,p5.2,p5.2.1,p5.2.2,p5.2.3,p5.2.4,p5.2.5,p5.2.6,p5.2.7), print))
      dev.off()
      
      
      
    }
    
  }
  

}


#6. E15分化路线统计
{
  E11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  E15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  E11tree = readRDS("final_result/xE11/E11_blay_tree.rds")
  E15tree = readRDS("final_result/xE15.5/E15.5_blay_taq_tree.rds")
  
  hubc = c("Rgl1","DA","Glu_M","Nb_DA","Nb_Glu","Nb_M")
  cm = E15tree$cm
  hubcm = cm[which(cm$celltype %in% hubc),]
  hubst = as.data.frame(table(hubcm[c("tags","celltype")]))
  # hubst = hubcm %>% group_by(tags,celltype) %>% summarise(counts = sum(n()))

  hubst = dcast(hubst,tags~celltype)
  rownames(hubst) = hubst$tags;hubst = hubst[-1]
  head(hubst)
  
  
  hubstn = t(t(as.matrix(hubst))/colSums(hubst))
  hubstn = hubstn[which(rowSums(hubstn>0)>1),]
  Heatmap(as.matrix(hubstn),row_km = 4,show_row_names = F)
  
  #尝试Upset图
  {
    library(UpSetR) 
    movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
    mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
    
    upset(movies)
    
    hubstn.up = hubstn
    hubstn.up[which(hubstn.up>0)] = T
    hubstn.up[which(hubstn.up==0)] = F
    hubstn.up = as.data.frame(hubstn.up)
    
    upset(hubstn.up,nsets = 6, keep.order = T)
    
    
    head(hubstn)
    nrow(hubstn.up[which(hubstn.up$Nb_M>0 & hubstn.up$Glu_M>0),])
    
    
    
  }
  
  #效果不好
  {
    #试试Apriori算法
    # Loading Libraries
    library(arules)
    library(arulesViz)
    library(RColorBrewer)
    
    #创建数据集
    hubapr = list()
    hubaprc= NULL
    for (i in 1:nrow(hubstn.up)) {
      line = hubstn.up[i,]
      set = colnames(line[,which(line>0)]) 
      set = set[order(set)]
      hubapr[[i]] = set
      names(hubapr)[i] = rownames(line)
      hubaprc = c(hubaprc,paste0(set,collapse = "-"))
    }
    table(hubaprc)
    
    MyTrans = as(hubapr,"transactions")
    summary(MyTrans)
    inspect(MyTrans)
    image(MyTrans)
    
    # import dataset
    data("Groceries")
    
    # using apriori() function
    rules <- apriori(MyTrans,
                     parameter = list(supp = 0.01, conf = 0.5,minlen=2))
    inspect(rules)
    size(rules)
    
    rules.sorted<-sort(x=rules,by="lift",decreasing=TRUE)
    inspect(rules.sorted)
    plot(rules.sorted,method = "grouped")
    plot(rules.sorted,method = "graph")
    
    
    # using inspect() function
    inspect(rules[1:10])
    inspect(head(Groceries, 3))
    
    # using itemFrequencyPlot() function
    arules::itemFrequencyPlot(MyTrans, topN = 20,
                              col = brewer.pal(8, 'Pastel2'),
                              main = 'Relative Item Frequency Plot',
                              type = "relative",
                              ylab = "Item Frequency (Relative)")
  }
  
  #自己构建变量作图
  {

    #试试把两群细胞合并
    merhub = hubst
    merhub$DAt = merhub$DA + merhub$Nb_DA
    merhub$Glut = merhub$Nb_Glu + merhub$Glu_M
    merhub = merhub[-c(1:4)]
    head(merhub)
    
    merhub.up = as.matrix(merhub)
    merhub.up[which(merhub.up>0)] = T
    merhub.up[which(merhub.up==0)] = F
    
    hubstn.up[which(hubstn.up$DA == 0 & hubstn.up$Nb_DA == 0 &
                      hubstn.up$Glu_M > 0 & hubstn.up$Nb_Glu > 0 & hubstn.up$Rgl3 == 0),]
    
    p6.1 = upset(as.data.frame(merhub.up),nsets = 4)
    p6.1
    ggexport(p6.1,filename = "final_result/xE15.5/Nb_direct_stat_upset.pdf",width = 8,height = 6)
    print(1)
    
    Heatmap(as.matrix(merhub),row_km = 4,show_row_names = F)
    #计算enrich
    #算一个total BC 富集值
    Mt = colSums(merhub[3:4])
    Nt = sum(Mt)
    
    merhubf = merhub[which(rowSums(merhub[3:4])>0),1:4]
    merhub.ft = NULL
    for (i in 1:nrow(merhubf)) {
      line = merhubf[i,]
      # line = merhub[which(rownames(merhub) == "v2.2I+98A_v2.54D+17_v2.81D+131"),]
      ftl = line
      if(rowSums(line[3:4]>0)==2){
        n = sum(line[3:4])
        k = unlist(line[3])
        M = Mt[which(names(Mt) == names(k))]
        N = Nt
        d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
        row.names(d) <- c("In_category", "not_in_category")
        ft = fisher.test(d)
        
        if(ft$p.value < 0.1){
          print(line)
          if(k/(n-k) > (M-k)/(N-M-n+k)){
            ftl[3] = 1
            ftl[4] = 0
          }else{
            ftl[3] = 0
            ftl[4] = 1
          }
        }
      }
      
      merhub.ft = rbind(merhub.ft,ftl)
    }
    
    merhubft.up = as.matrix(merhub.ft)
    merhubft.up[which(merhubft.up>0)] = T
    merhubft.up[which(merhubft.up==0)] = F
    p6.2 = upset(as.data.frame(merhubft.up),nsets = 4, keep.order = T,order.by = "freq")
    p6.2
    ggexport(p6.2,filename = "final_result/xE15.5/DA_Glu_direction_stat_upset_ft.pdf",width = 8,height = 6)
    p6.1
    
    merhub.ft[which(rowSums(merhub.ft[3:4]>0)>1),]
    
    #对分类结果作图
    #分组
    group.da = merhub.ft[which(merhub.ft$DAt > 0  & merhub.ft$Glut == 0),]
    group.glu = merhub.ft[which(merhub.ft$DAt == 0 & merhub.ft$Glut > 0),]
    group.all = merhub.ft[which(merhub.ft$DAt > 0 & merhub.ft$Glut > 0),]
    
    
    merhubfl = as.data.frame(t(t(as.matrix(merhubf))/colSums(merhubf)))
    
    merhubfl = cbind(rownames(merhubfl),merhubfl)
    merhubfl = melt(merhubfl)
    colnames(merhubfl) = c("tags","cell","proportion")
    merhubfl = merhubfl %>% group_by(tags) %>% summarise(cell = cell,proportion = proportion,
                                                             abspro = proportion/sum(proportion))
    merhubfl$group = "unclass"
    merhubfl[which(merhubfl$tags %in% rownames(group.da)),"group"] = "DA"
    merhubfl[which(merhubfl$tags %in% rownames(group.glu)),"group"] = "Glu"
    merhubfl[which(merhubfl$tags %in% rownames(group.all)),"group"] = "DA_Glu"
    merhubfl = merhubfl[order(merhubfl$group),]
    merhubfl$tags = factor(merhubfl$tags,levels = unique(merhubfl$tags))
    
    p6.3 = ggplot(merhubfl,aes(x = tags,y = cell,fill = proportion)) + 
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
    p6.3
    ggsave(p6.3,filename = "final_result/xE15.5/DA_Glu_direction_point.pdf",width = 10,height = 3)
    # Heatmap(as.matrix(scale(merhubf)),show_row_names = F,heatmap_width = unit(0.5,"npc"),width = 0.5)
    
    
  }
  
  #根据构建的分组表征clone并差异分析
  {
    E15trans
    
    hubcb = names(E15trans$ClusterName[which(E15trans$ClusterName %in% hubc)])
    hubcb = substr(hubcb,1,16)
    
    trans.da = E15tree$cm[which(E15tree$cm$tags %in% "v2.140D+44_v2.8D+199"),"Var1"]
    trans.da2 = E15tree$cm[which(E15tree$cm$tags %in% "v2.2D+123_v2.2D+179_v2.36D+71_v2.46D+200_v2.8D+14"),"Var1"]
    trans.glu = E15tree$cm[which(E15tree$cm$tags %in% "v2.11D+206_v2.162D+17"),"Var1"]
    trans.all = E15tree$cm[which(E15tree$cm$tags %in% "v2.128D+12_v2.21D+151_v2.5D+204"),"Var1"]
    cell_names.da <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.da]
    cell_names.glu <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.glu]
    cell_names.all <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.all]
    # 
    # 
    # tmp = table(substr(colnames(E15trans),1,16))
    # tmp[which(tmp>1)]
    
    p6.4.0 = DimPlot(E15trans,label = TRUE,label.size = 3) + NoLegend()
    p6.4.1 = DimPlot(E15trans,label = TRUE, label.size = 3,cells.highlight = cell_names.da, 
            cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.2 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.glu, 
            cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.3 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.all, 
            cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    
    p6.4 = p6.4.0 + p6.4.1 + p6.4.2 + p6.4.3
    p6.4
    ggsave(p6.4,filename = "final_result/xE15.5/DA_Glu_direction_umap_example.pdf",width = 10,height = 8)
    #contour plot
    
    
    #差异分析
    trans.da = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.da)),"Var1"]
    trans.glu = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.glu)),"Var1"]
    trans.all = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.all)),"Var1"]
    
    # E15trans.hub = subset(E15trans, idents = c("Rgl1"))
    E15trans.hub = subset(E15trans, idents = c("Nb_M","DA","Nb_DA","Glu_ML"))
    
    tmp = Idents(E15trans.hub)
    tmp = as.character(tmp);names(tmp) = names(Idents(E15trans.hub))
    tmp[which(substr(names(tmp),1,16) %in% trans.da)] = "DA_d"
    tmp[which(substr(names(tmp),1,16) %in% trans.glu)] = "Glu_d"
    tmp[which(substr(names(tmp),1,16) %in% trans.all)] = "All_d"
    tmp = as.factor(tmp)
    Idents(E15trans.hub) = tmp
    
    E15trans.hub = subset(E15trans.hub, idents = c("DA_d","Glu_d","All_d"))
    saveRDS(E15trans.hub,"final_result/xE15.5/E15trans_DA.rds")
    
    DA.markers <- FindMarkers(E15trans.hub, ident.1 = "DA_d",min.pct = 0.25)
    Glu.markers <- FindMarkers(E15trans.hub, ident.1 = "Glu_d", min.pct = 0.25)
    All.markers <- FindMarkers(E15trans.hub, ident.1 = "All_d", min.pct = 0.25)
    head(Glu.markers)
    head(DA.markers)
    head(All.markers)
    head(Glu.markers[order(-Glu.markers$avg_log2FC),])
    head(All.markers[order(-All.markers$avg_log2FC),])
    # Markgene = c("En2", "En1","Pbx1","Nr2f2","Ntm","Tle4","Lhx1os","C1ql1","Nefl","Pitx3")
    Markgene = c("Slc4a4", "Timp4","Fgfbp3","Tafa5","Ecm2","Nap1l5","Tmem184b","Ank","Ocrl")
    p6.5.1 = VlnPlot(E15trans.hub, features = Markgene)
    p6.5.1
    p6.5.2 = FeaturePlot(E15trans.hub, features = Markgene)
    p6.5.2
    #火山图
    VocanoFigure = function(DA.markers,title){
      DA.markers$reg = "NoSig"
      DA.markers[which(DA.markers$p_val<0.05 & DA.markers$avg_log2FC > 1),"reg"] = "Up"
      DA.markers[which(DA.markers$p_val<0.05 & DA.markers$avg_log2FC < -1),"reg"] = "Down"
      DA.markers$gene =rownames(DA.markers)
      DA.markers$reg = factor(DA.markers$reg,levels = c("Down","NoSig","Up"))
      # library(ggrepel)
      p6.5.3 = ggplot(DA.markers,aes(x = avg_log2FC,y = -log10(p_val),color = reg)) + geom_point() + 
        scale_color_manual(values=c("#DC143C","#808080","#00008B"))+#确定点的颜色
        ggrepel::geom_text_repel(
          data = DA.markers[which(DA.markers$p_val<0.05 & abs(DA.markers$avg_log2FC) > 1),],
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
        geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
        geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
      p6.5.3
      return(p6.5.3)
    }
    p6.5.3.1 = VocanoFigure(DA.markers,"DA_unique")
    p6.5.3.2 = VocanoFigure(Glu.markers,"Glu_unique")
    p6.5.3.3 = VocanoFigure(All.markers,"All_unique")
    p6.5.3 = p6.5.3.1 + p6.5.3.2 + p6.5.3.3
    p6.5.3
    
    pdf("final_result/xE15.5/DA_Glu_direction_DEG_analysis_Rgl1.pdf",width = 10,height = 6)
    invisible(lapply(list(p6.5.1,p6.5.2,p6.5.3.1,p6.5.3.2,p6.5.3.3), print))
    dev.off()
    write.csv(DA.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_DA_Rgl1.csv",row.names = T,quote = F)
    write.csv(Glu.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_Glu_Rgl1.csv",row.names = T,quote = F)
    write.csv(All.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_All_Rgl1.csv",row.names = T,quote = F)
    
  }
  
}

 
#7. E11桑基图，并与E15连接
{
  cm = E11tree$cm
  cmtn = merge(cm,pseudo11,by.x = "Var1",by.y = "BC")
  cmtn$domain = "E11"
  # lateid = c("NProg_ML_m6","NProg_Foxa2_m5","NProg_Six3","Nb_M",
  #             "NProg_L","Rgl1","Rgl2_m4","Rgl2_Pax7_m1")
  earlyid = c("NProg_L","Prog_BP","Rgl1")
  lateid = setdiff(unique(cmtn$celltype),earlyid)
  
  cmtn$phase = "late"
  cmtn[which(cmtn$celltype %in% earlyid),"phase"] = "early"
  
  p7.1 = ggboxplot(cmtn[which(cmtn$domain != "unknown"),],x = "celltype",y = "pseudo",color = "phase") + 
    theme_bw() + 
    theme(axis.text.x =  element_text(angle = 45,vjust = 0.5))
  p7.1
  
  cellps = cmtn %>% group_by(celltype) %>% summarise(pseudo.mean = mean(pseudo),domain = unique(domain))
  
  
  cmtnp2 = cmtn %>% group_by(tags,phase,celltype) %>% summarise(counts = sum(n()))
  
  #尝试计算每个BC的富集值，利用富集值来做桑基图
  {
    cmtnp3 = cmtn %>% group_by(tags,phase,celltype) %>% 
      summarise(counts = sum(n()),domain = unique(domain),pseudo.mean = mean(pseudo))
    head(cmtnp3)
    
    clone = unique(cmtnp3$tags)
    
    hpsdf = NULL
    for (i in 1:length(clone)) {
      line = cmtnp3[which(cmtnp3$tags == clone[i]),]
      sdm = line$counts
      sdm = sdm[order(-sdm)]
      sdm = sdm/sum(sdm)
      H = -sum(sdm*log(sdm))
      # H = -sum(sdm*log(sdm))/log2(length(sdm))
      psm = mean(line$pseudo.mean)
      psd = sd(line$pseudo.mean)
      size = sum(line$counts)
      celltypenum = nrow(line)
      linedf = data.frame("H" = H,"pseudo.mean" = psm,"pseudo.sd" = psd,"size" = size,
                          "celltypenum" = celltypenum)
      hpsdf =rbind(hpsdf,linedf)
    }
    hpsdf[which(is.nan(hpsdf$H)),"H"] = 0
    hpsdf[which(is.na(hpsdf$pseudo.sd)),"pseudo.sd"] = 0
    head(hpsdf)
    
    ggplot(hpsdf,aes(x = pseudo.mean, y = celltypenum,color = celltypenum,alpha = 0.7)) + 
      geom_point(aes(size = size)) + 
      scale_size_continuous(range = c(0, 8)) +
      theme_bw()
    ggplot(hpsdf,aes(x = celltypenum, y = H,alpha = 0.7)) + 
      geom_point(aes(size = size)) + 
      scale_size_continuous(range = c(0, 8)) +
      theme_bw()
    
    
    cellmx.e = dcast(cmtnp3[which(cmtnp3$phase == "early"),c("tags","celltype","counts")],
                     tags~celltype,value.var = "counts",fun.aggregate = sum)
    cellmx.l = dcast(cmtnp3[which(cmtnp3$phase == "late"),c("tags","celltype","counts")],
                     tags~celltype,value.var = "counts",fun.aggregate = sum)
    Me = cmtnp3[which(cmtnp3$phase == "early"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
    Ne = sum(Me$counts)
    
    Ml = cmtnp3[which(cmtnp3$phase == "late"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
    Nl = sum(Ml$counts)
    
    cellft.e = NULL
    for (i in 1:nrow(cellmx.e)) {
      line = cellmx.e[i,]
      n = sum(line[2:ncol(line)])
      ftl = line
      for (j in 2:ncol(line)) {
        k = unlist(line[j])
        if(k>0){
          M = Me[which(Me$celltype == names(k)),]$counts
          N = Ne
          d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
          row.names(d) <- c("In_category", "not_in_category")
          ft = fisher.test(d)
          ftl[j] = -log10(ft$p.value)
        }
      }
      cellft.e = rbind(cellft.e,ftl)
    }
    
    
    cellft.l = NULL
    for (i in 1:nrow(cellmx.l)) {
      line = cellmx.l[i,]
      n = sum(line[2:ncol(line)])
      ftl = line
      for (j in 2:ncol(line)) {
        k = unlist(line[j])
        if(k>0){
          M = Ml[which(Ml$celltype == names(k)),]$counts
          N = Nl
          d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
          row.names(d) <- c("In_category", "not_in_category")
          ft = fisher.test(d)
          ftl[j] = -log10(ft$p.value)
        }
      }
      cellft.l = rbind(cellft.l,ftl)
    }
    cellft = merge(cellft.e,cellft.l,by = "tags")
    
    #计算cellft
    phasest = data.frame("phase1" = rep(earlyid,each = length(lateid)),
                         "phase2" = rep(lateid,length(earlyid)), "Freq" = 0)
    for (i in 1:nrow(cellft)) {
      line = cellft[i,]
      line.e = unlist(line[,which(colnames(line) %in% earlyid)])
      line.e = line.e[which(line.e>0)]
      line.l = unlist(line[,which(colnames(line) %in% lateid)])
      line.l = line.l[which(line.l>0)]
      for (j in 1:length(line.e)) {
        for (k in 1:length(line.l)) {
          seg = phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq"]
          phasest[which(phasest$phase1 == names(line.e)[j] & phasest$phase2 == names(line.l)[k]),"Freq"] = seg + mean(c(line.e[j],line.l[k]))
        }
      }
    }
    
    
    phasestf = phasest
    gghistogram(phasestf$Freq) + scale_x_continuous(breaks = seq(0,200,10))
    
    phasestf = phasestf[which(phasestf$Freq>20),]
    
    phasest_lod = to_lodes_form(phasestf,
                                key = "phase", value = "celltype", id = "alluvium",
                                axes = 1:2)
    phasest_lod = phasest_lod[order(phasest_lod$phase),]
    head(phasest_lod)
    phasest_lod$celltype = factor(phasest_lod$celltype,
                                  levels = c("Glu_dorsal","GABA_m5","NProg_GABA_late",
                                             "GABA4","NProg_GABA_early","NProg_Glu_Dorsal",
                                             "NProg_L","Peric","OMTN","Glu_ML","Prog_BP",
                                             "Nb_ML","Nb_M","Rgl1"))
    
    
    p7.4 = ggplot(phasest_lod,
                  aes(x = phase, stratum = celltype, alluvium = alluvium,
                      y = Freq,
                      fill = celltype, label = celltype)) +
      scale_x_discrete(expand = c(.1, .1)) +
      # geom_flow(aes(color = celltype)) +
      geom_flow(aes()) +
      geom_stratum(alpha = .5) +
      geom_text(stat = "stratum", size = 3) + theme_void()
    p7.4
    ggsave(p7.4,filename = "final_result/xE11/E11_pseudo_phase_module_stratum_celltype_arti_fisher_filter20.pdf",width = 8,height = 5)
    
    
  }
  
}


#8. 11.17补图
{
  E11tree = readRDS("final_result/xE11/E11_blay_tree.rds")
  E15tree = readRDS("final_result/xE15.5/E15.5_blay_taq_tree.rds")
  earlyid11 = c("NProgBL","NProgBM","Rgl1")
  earlyid15 = c("NProgAL","Rgl1","Rgl2AL","Rgl2BL","GabaProgBI","GabaProgBL","NbBM1","NbFP","Rgl3")
  corder11 = E11tree$distfigure$spearhc$labels[E11tree$distfigure$spearhc$order]
  corder11 = c('GluAL','NbGluAL','NProgBL','GabaProg0','GabaProgBL','GabaProgBI','NProgBM','OMTN', 'NbBM0','NbBM1','Rgl1','NbFP','Peric')
  corder11 = c("NbBM1","NbBM0","NProgBM","GabaProg0","GabaProgBL","NbGluAL","NProgBL",
               "OMTN","Peric","GabaProgBI","NbFP","Rgl1","GluAL")
  corder15 = E15tree$distfigure$abspearsonhc$labels[E15tree$distfigure$abspearsonhc$order]
  
  #1. 桑基图改善---------------------------------------------------------
  #(1) 单独的E11,E15tree
  E11tree = MyHeatmap(E11tree,"final_result/xE11/E11_heatmap")
  E15tree = MyHeatmap(E15tree,"final_result/xE15.5/E15_heatmap")
  
  #(2) 滤掉widespread的barcode，画弦图
  #封装函数 计算每个BC的富集值
  pseudo15 = read.csv("final_result/xE15.5/E15.5_monocle_pseudo.csv")
  pseudo11 = read.csv("final_result/xE11/E11_monocle_pseudo.csv")
  domain11 = read.csv("sc_result/X.E11/E11pseudoAxis20211126V4.csv")
  domain15 = read.csv("sc_result/x.TAQ_E15/E15pseudoAxis20211126v5.csv")
  colnames(domain11) = colnames(domain15) = c("celltype","domainid","domain")
  
  
  tree = E15tree
  pseudo = pseudo15
  domain = domain15
  earlyid = earlyid15
  
  #计算富集度
  BCEnrichment = function(tree,domain,earlyid,outname){
      cm = tree$cm
      cmtn = merge(cm,domain,by = "celltype")
      lateid = setdiff(unique(cmtn$celltype),earlyid)
      
      cmtn$phase = "late"
      cmtn[which(cmtn$celltype %in% earlyid),"phase"] = "early"
      
      cmtn = cmtn[order(cmtn$phase,cmtn$domain),]
      
      #尝试计算每个BC的富集值，利用富集值来做桑基图
      {
        cmtnp = cmtn %>% group_by(tags,phase,celltype) %>% 
          summarise(counts = sum(n()),domain = unique(domain))
        #
        cellmx.e = dcast(cmtnp[which(cmtnp$phase == "early"),c("tags","celltype","counts")],
                         tags~celltype,value.var = "counts",fun.aggregate = sum)
        cellmx.l = dcast(cmtnp[which(cmtnp$phase == "late"),c("tags","celltype","counts")],
                         tags~celltype,value.var = "counts",fun.aggregate = sum)
        Me = cmtnp[which(cmtnp$phase == "early"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
        Ne = sum(Me$counts)
        
        Ml = cmtnp[which(cmtnp$phase == "late"),] %>% group_by(celltype) %>% summarise(counts = sum(counts))
        Nl = sum(Ml$counts)
        
        cellft.e = NULL
        for (i in 1:nrow(cellmx.e)) {
          line = cellmx.e[i,]
          n = sum(line[2:ncol(line)])
          ftl = line
          for (j in 2:ncol(line)) {
            k = unlist(line[j])
            if(k>0){
              M = Me[which(Me$celltype == names(k)),]$counts
              N = Ne
              d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
              row.names(d) <- c("In_category", "not_in_category")
              ft = fisher.test(d)
              ftl[j] = -log10(ft$p.value)
            }
          }
          cellft.e = rbind(cellft.e,ftl)
        }
        
        
        cellft.l = NULL
        for (i in 1:nrow(cellmx.l)) {
          line = cellmx.l[i,]
          n = sum(line[2:ncol(line)])
          ftl = line
          for (j in 2:ncol(line)) {
            k = unlist(line[j])
            if(k>0){
              M = Ml[which(Ml$celltype == names(k)),]$counts
              N = Nl
              d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
              row.names(d) <- c("In_category", "not_in_category")
              ft = fisher.test(d)
              ftl[j] = -log10(ft$p.value)
            }
          }
          cellft.l = rbind(cellft.l,ftl)
        }
        cellft = merge(cellft.e,cellft.l,by = "tags")
      }
      cellft = cellft[which(rowSums(cellft[-1]>0)< ncol(cellft[-1])/3),]
      #绘制桑基图
      {
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
        
        domain_df = unique(cmtn[c("celltype","domainid","domain")])
        domain_df = domain_df[order(domain_df$domain),]
        phasestf = phasest[which(phasest$Freq.from > 0.001 | phasest$Freq.to > 0.001),]
        # phasestf = phasestf[which(phasestf$Freq>20),]
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
        phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domainid),]
        phasest_lod$celltype = factor(phasest_lod$celltype,
                                      levels = unique(as.character(phasest_lod$celltype)))
        
        
        pf1 = ggplot(phasest_lod,
                      aes(x = phase, stratum = celltype, alluvium = alluvium,
                          y = Freq,
                          fill = domain, label = celltype)) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_flow(aes(color = celltype)) +
          geom_flow(aes()) +
          geom_stratum(alpha = .5) +
          geom_text(stat = "stratum", size = 3) + theme_void()
        pf1
        
        library(circlize)
        library(randomcoloR)
        require(RColorBrewer)
        brewer.pal(9, "Set1")
        
        phasestm = phasest
        phasestm$Freq.from = phasestm$Freq.from/sum(phasestm$Freq.from)
        phasestm$Freq.to = phasestm$Freq.to/sum(phasestm$Freq.to)
        
        phasestm = merge(phasestm,domain_df,
                            by.x = "phase1",by.y = "celltype")
        phasestm = phasestm[which(phasestm$phase1 > 0.001 | phasestm$phase2 > 0.001),]
        gridcol = data.frame(celltype = unique(c(phasestm$phase1,phasestm$phase2)),
                             phase = c(rep(1,length(unique(phasestm$phase1))),
                                       rep(2,length(unique(phasestm$phase2)))))
        gridcol$domain = domain[match(gridcol$celltype,domain$celltype),"domain"]
        domaincol = data.frame(domain = unique(gridcol$domain),
                               col = distinctColorPalette(length(unique(gridcol$domain))))
        gridcol$col = domaincol[match(gridcol$domain,domaincol$domain),"col"]
        gridcol = gridcol[order(-gridcol$phase,gridcol$domain),]
        tmp = gridcol[which(gridcol$phase == 1),]
        gridcol[which(gridcol$phase == 1),] = tmp[order(tmp$domain,decreasing = T),]
        gridcolv = gridcol$col;names(gridcolv) = gridcol$celltype
        
        pdf(paste0(outname,"_chord_plot.pdf"))
        circos.par(start.degree = 90)
        chordDiagram(phasestm[,1:4],annotationTrack =  c("name", "grid"),
                     grid.col = gridcolv,
                     directional = 1, direction.type = c("diffHeight", "arrows"),
                     link.arr.type = "big.arrow",
                     order = unique(gridcol$celltype))
        # Restart circular layout parameters
        circos.clear()
        dev.off()
        
      }
      
      return(list("cellft" = cellft,"phasest" = phasest,"phasest_lod" = phasest_lod,"figure.allu" = pf1))
  }
  
  cellallu15 = BCEnrichment(E15tree,domain15,earlyid15,"final_result/xE15.5/alluv15")
  cellallu11 = BCEnrichment(E11tree,domain11,earlyid11,"final_result/xE11/alluv11")
  saveRDS(cellallu11,file = "final_result/xE11/alluv11.rds")
  saveRDS(cellallu15,file = "final_result/xE15.5/alluv15.rds")
  
  p2.1 = cellallu11$figure.allu
  p2.2 = cellallu15$figure.allu
  ggsave(p2.2,filename = "final_result/xE15.5/alluv15_plot.pdf")
  
  #结合二者的桑基图
  phasest15 = cellallu15$phasest
  phasest11 = cellallu11$phasest
  # tmp11 = phasest11[which(phasest11$phase2 %in% phase2id),] %>% group_by(phase2) %>% summarise(count = sum(Freq.to))
  # tmp11$norm = tmp11$count/ sum(tmp11$count)
  # tmp15 = phasest15[which(phasest15$phase1 %in% phase2id),] %>% group_by(phase1) %>% summarise(count = sum(Freq.to))
  # tmp15$norm = tmp15$count/ sum(tmp15$count)
  # tmp11
  # tmp15
  
  phase1id = unique(cellallu11$phasest$phase1)
  phase2id = c("GabaProgBI","GabaProgBL","NbBM1","NbFP")
  phase3id = unique(cellallu15$phasest$phase2)
  phase3id = phase3id[which(!phase3id %in% c("GabaAL","GluAL1","GluAL2","GluAL3","GluAL4","Rgl3"))]
  phasest = data.frame("phase1" = rep(phase1id,each = length(phase2id)*length(phase3id)),
                       "phase2" = rep(rep(phase2id,each = length(phase3id)),length(phase1id)),
                       "phase3" = rep(phase3id,length(phase1id)*length(phase2id)),
                       "Freq" = 0)
  
  phasest_lod = to_lodes_form(phasest,
                              key = "phase", value = "celltype", id = "alluvium",
                              axes = 1:3)
  domaint = unique(rbind(domain11,domain15))
  phasest_lod = merge(phasest_lod,domaint[c("celltype","domain")],by = "celltype")
  phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domain),]
  phasest_lod$celltype = factor(phasest_lod$celltype,
                                levels = unique(as.character(phasest_lod$celltype)))
  phasest_lod$domain = factor(phasest_lod$domain)
  
  for (i in 1:max(phasest_lod$alluvium)) {
    line = phasest_lod[which(phasest_lod$alluvium == i),]
    p1c = as.character(line[which(line$phase == "phase1"),"celltype"])
    p2c = as.character(line[which(line$phase == "phase2"),"celltype"])
    p3c = as.character(line[which(line$phase == "phase3"),"celltype"])
    freq1 = phasest11[which(phasest11$phase1 == p1c & phasest11$phase2 == p2c),"Freq.from"]
    freq21 = phasest11[which(phasest11$phase1 == p1c & phasest11$phase2 == p2c),"Freq.to"]
    freq22 = phasest15[which(phasest15$phase1 == p2c & phasest15$phase2 == p3c),"Freq.from"]
    freq2 = mean(c(freq21,freq22))
    freq3 = phasest15[which(phasest15$phase1 == p2c & phasest15$phase2 == p3c),"Freq.to"]
    
    phasest_lod[which(phasest_lod$alluvium == i & phasest_lod$phase == "phase1"),"Freq"] = freq1
    phasest_lod[which(phasest_lod$alluvium == i & phasest_lod$phase == "phase2"),"Freq"] = freq2
    phasest_lod[which(phasest_lod$alluvium == i & phasest_lod$phase == "phase3"),"Freq"] = freq3
    
    
  }
  phasest_lod = phasest_lod %>% group_by(phase) %>% summarise(celltype = celltype,Freq=Freq,norm = Freq/sum(Freq),
                                                              domain = domain ,alluvium = alluvium)
  phasest_lod = as.data.frame(phasest_lod)
  
  filterid = unique(phasest_lod[which(phasest_lod$norm>0.001),"alluvium"])
  phasest_lodf = phasest_lod
  phasest_lodf = phasest_lod[which(phasest_lod$alluvium %in% filterid),]
  phasest_lodf = phasest_lodf %>% group_by(phase) %>% summarise(celltype = celltype,Freq=Freq,norm = Freq/sum(Freq),
                                                              domain = domain ,alluvium = alluvium)
  
  p2.3 = ggplot(phasest_lodf,
         aes(x = phase, stratum = celltype, alluvium = alluvium,
             y = norm,
             fill = domain,
             label = celltype)) +
    scale_x_discrete(expand = c(.1, .1)) +
    # geom_flow(aes(color = celltype)) +
    geom_flow(aes()) +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) + theme_void()
  p2.3
  
  pdf("final_result/xtotal/E11_E15_alluvium.pdf",width = 10,height = 6)
  invisible(lapply(list(p2.3), print))
  dev.off()
  write.csv(phasest_lodf,file = "final_result/xtotal/E11_E15_alluvium.csv",row.names = F,quote = F)
  
  #2.e11 circle plot 有问题,已解决-------------------------------
  
  #3.pseudotime disper 图------------------------------------
  PseudoDisper = function(tree,domain){
    library(scatterpie)
    cm = tree$cm
    # cmtn = merge(cm,pseudo,by.x = "Var1",by.y = "BC")
    cmtn = merge(cm,domain,by = "celltype")
    cmtn = cmtn[order(cmtn$domainid),]
    cmtn$celltype = factor(as.character(cmtn$celltype),levels = unique(cmtn$celltype))
    unique(domain$domain)
    
    #计算array的分散度
    # cmtnf = cmtn[which(cmtn$domain != "unknown"),]
    cmtnf = cmtn[which(cmtn$domain != "unknown"),]
    cmtnf$celltype = factor(cmtnf$celltype,levels = unique(cmtnf$celltype))
    cmtnf = cmtnf[order(cmtnf$domainid),]
    cmtnf$domainid = factor(cmtnf$domainid)
    # cmtnf$domainid = factor(as.numeric(factor(cmtnf$domain,levels = unique(cmtnf$domain))))
      
    root_spr  =NULL
    root_spr_mod  =NULL
    for (i in 1:length(unique(cmtnf$tags))) {
        root = unique(cmtnf$tags)[i]
        slice = cmtnf[which(cmtnf$tags == root),]
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
          line_m = data.frame(group = root,H = H_m,Hwi = Hwi_m,Hb = Hb,
                              main.domain = mainModel_m,
                              # pseudo.mean = pseudom,
                              count = nrow(slice),
                              stringsAsFactors = F)
          cc_m = as.data.frame.list(cc_m)
          line_m = cbind(line_m,cc_m)
          root_spr_mod = rbind(root_spr_mod,line_m)
          
          mainModel = names(table(slice$celltype))[which.max(table(slice$celltype))]
          meanModel = mean(as.numeric(as.character(slice$domainid)))
          
          
          #cell H
          cc = table(slice$celltype)
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
          cc = table(slice$celltype)
          cc = as.data.frame.list(cc)
          line = cbind(line,cc)
          root_spr = rbind(root_spr,line)
        }
        
      }
    colnames(root_spr_mod)[-(1:6)] = unique(cmtnf$domain)
    root_sprn = root_spr
    root_sprn[-(1:9)] = t(t(as.matrix(root_sprn[-c(1:9)]))/colSums(root_sprn[-(1:9)]))
      
    # tem2 <- as.data.frame(spline(root_sprn$mean.model, root_sprn$Hb,n = nrow(root_sprn),ties = mean))
    p5.1 = ggplot() + 
        geom_scatterpie(data = root_sprn,aes(x=mean.model, y=Hb, group=group,
                                             r = log10(count)/15,alpha = 0.5),color=NA,
                        cols=colnames(root_sprn)[-(1:9)]) + 
        # geom_histogram(data = root_sprn,aes(mean.model),bins = 100) +
        coord_fixed() +
        # geom_smooth(data = tem2, aes(x=x, y=y), se = F, method = 'loess',color = "black") +
        scale_fill_viridis(discrete=TRUE) + scale_y_continuous(breaks = c(0,1,2,3)) + 
        scale_x_continuous(breaks = c(1,2,3,4,5)) + 
        xlab("mean.model") + ylab("Hb") + 
        theme_bw() + theme(legend.position="bottom")
    p5.1
    
    p5.1.1 = ggplot(root_sprn,aes(x=mean.model)) + 
      # scale_fill_continuous() +
      geom_histogram(bins = 200,fill = "black") + 
      theme_bw() +
      theme(axis.title.x = element_blank())
            # axis.text.x = element_blank(),)
    p5.1.1
    
    p5.1.t = ggarrange(p5.1.1, p5.1, heights = c(0.2, 1),
              ncol = 1, nrow = 2)
    
    root_spr_modn = root_spr_mod
    root_spr_modn[-(1:6)] = t(t(as.matrix(root_spr_modn[-c(1:6)]))/colSums(root_spr_modn[-(1:6)]))
      
    #5.2延伸
    root_spr_modn$max = 0
    for (i in 1:nrow(root_spr_modn)) {
        line = root_spr_modn[i,]
        root_spr_modn[i,"main.domain"] = names(which.max(line[-(1:6)]))
        root_spr_modn$max[i] = max(line[-(1:6)]/sum(line[-(1:6)]))
    }
    p5.2.5 = ggplot(root_spr_modn,aes(x = main.domain, y = Hb)) +
        geom_boxplot(width = 0.6) + scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans")
    p5.2.5
      
    
    return(list("pvs_cell" = root_sprn,"pvs_domain" = root_spr_modn,
                "pvs_cell_count" = root_spr,"pvs_domain_count" = root_spr_mod,
                "figure_pvs" = p5.1.t,"figure_pvs_box" = p5.2.5))
    
  }
  # pvs11 = PseudoDisper(E11tree,pseudo11,domain11)
  pvs11 = PseudoDisper(E11tree,domain11)
  # pvs15 = PseudoDisper(E15tree,pseudo15,domain15)
  pvs15 = PseudoDisper(E15tree,domain15)
  pvs11$figure_pvs
  pvs15$figure_pvs
  pvs11$figure_pvs_box
  ggsave(pvs11$figure_pvs,filename = "final_result/xE11/E11_pseudo_disper_cmb.pdf",width = 8,height = 6)
  ggsave(pvs15$figure_pvs,filename = "final_result/xE15.5/E15_pseudo_disper_cmb.pdf",width = 8,height = 6)
  
  #domain heatmap
  ardm = pvs11$pvs_domain
  ardm = ardm[-c(1:7,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xE11/E11_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               row_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  ardm = pvs15$pvs_domain
  ardm = ardm[-c(1:7,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xE15.5/E15_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
               row_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  
  
  #4.Clone size vs enrich.degree------------------
  pvs11d = pvs11$pvs_cell
  pvs15d = pvs15$pvs_cell
  
  pvsd = pvs15d
  pvsd$cellclass = rowSums(pvsd[-c(1:10)] > 0)
  pvsd$complexity = pvsd$cellclass*(1+pvsd$Hwic)
  
  ggplot(pvsd,aes(x = log(count), y = complexity,color = main.cell,size = Hc)) + 
    geom_point(alpha = 0.5)  +
    scale_color_viridis(discrete=TRUE) + theme_ipsum(base_family = "sans")
  
  
  
  #5.E11 clone 聚类图------------------
  #细胞聚类
  {
    cmtnp = E11tree$cm %>% group_by(tags,celltype) %>% 
      summarise(counts = sum(n()))
    #算一个total BC 富集值
    cellmx = dcast(cmtnp[,c("tags","celltype","counts")],
                    tags~celltype,value.var = "counts",fun.aggregate = sum)
    Mt = cmtnp %>% group_by(celltype) %>% summarise(counts = sum(counts))
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
    cellftl = melt(cellft)
    head(cellftl)
    
    #尝试聚类分析
    rownames(cellft) = cellft$tags;cellft = cellft[-1];cellftn = scale(cellft)
    
    tmp = table(cellftl[which(cellftl$value>0),]$tags)
    tmp = tmp[which(tmp>1)]
    cellftf = cellft[which(rownames(cellft)%in% names(tmp)),]
    cellmxf = cellmx[-1];rownames(cellmxf) = cellmx$tags;cellmxf = cellmxf[which(rownames(cellmxf)%in% names(tmp)),]
    
    library(ggfortify)
    library(ggplot2)
    
    colnames(cellftf)
    
    cellorder = corder11
    
    n = 30
    
    # fit = kmeans(cellftf,n)
    # split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    # reorder.hmap <- Heatmap(as.matrix(cellftf), split=split, cluster_row_slices = FALSE,
    #                         show_row_names = F,column_order = cellorder)
    # reorder.hmap
    
    # ggexport(reorder.hmap, filename = "final_result/xE15.5/tag_cluster_heatmap.pdf",width = 8,height = 8)
    
    cellmxfn = t(t(as.matrix(cellmxf))/colSums(cellmxf))
    fit = kmeans(cellmxfn,n)
    split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    reorder.hmap <- Heatmap(as.matrix(cellmxfn), split=split, cluster_row_slices = FALSE,
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    #merge cluster
    clist = list(c2 = c(15,20,23,29),
                 c7 = c(28,22,18,7,19),
                 c4 = c(1,2,3,11,30),
                 c3 = c(6,8,13,16,21,22,24),
                 c6 = c(17,25,14),
                 c9 = c(10,12),
                 c5 = c(5,9),
                 c8 = c(4,27),
                 c1 = c(26))
                 # c10 = c(29),)
    split <- fit$cluster
    for (i in 1:length(clist)) {
      split[which(split %in% clist[[i]])] = names(clist)[i]
    }
    
    split <- factor(as.character(split), levels=as.character(c(paste0("c",1:length(clist)) )))      
    names(split) = names(fit$cluster)
    # e11_clu = readRDS("final_result/xE11/E11_cellcounts_cluster.rds")
    reorder.hmap <- Heatmap(as.matrix(e11_clu$cellcounts), split=e11_clu$cluster, cluster_row_slices = FALSE,
                            col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    ggexport(reorder.hmap, filename = "final_result/xE11/E11_tag_cluster_heatmap_counts_arti.pdf",width = 8,height = 8)
    
    saveRDS(list(cluster = split,fit = fit,cellcounts = cellmxfn,cluster.fig = reorder.hmap),
            file = "final_result/xE11/E11_cellcounts_cluster.rds")
    
    #cluster show in dimplot
    ccls = readRDS("final_result/xE11/E11_cellcounts_cluster.rds")
    
    ccidents = as.data.frame(ccls$cluster)
    ccidents = merge(ccidents,E11tree$cm,by.x = "row.names",by.y = "tags")
    colnames(ccidents) = c("tags","cluster","CellBC","Freq","Celltype")
    ccidents$Celltype = factor(ccidents$Celltype,levels = cellorder)
    write.csv(ccidents,file = "final_result/xE11/E11_tag_cluster_BC.csv",row.names = F,quote = F)
    
    p3.2 = ggplot(ccidents,aes(x = cluster,fill = Celltype)) +
      geom_bar(stat = "count",position = "fill",width = 0.7) +
      theme_bw()
    p3.2
    ggsave(p3.2,filename = "final_result/xE11/E11_tag_cluster_celltype_bar.pdf",width = 8,height = 6)
    
    p3.3 = list()
    for (i in 1:9) {
      cc = ccidents[which(ccidents$cluster == paste0("c",i)),]
      tarcell = colnames(E11trans)[which(substr(colnames(E11trans),1,16) %in% cc$CellBC)]
      pt = DimPlot(E11trans,label = TRUE, label.size = 3,cells.highlight = tarcell, 
                   cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
      p3.3[[i]] = pt
      names(p3.3)[i] = paste0("c",i)
      
    }
    p3.3.0 = ggarrange(plotlist =  p3.3,labels = names(p3.3))
    p3.3.0
    
    pdf("final_result/xE11/E11_tag_cluster_umap.pdf",width = 10,height = 8)
    print(p3.3.0)
    invisible(lapply(p3.3, print))
    dev.off()
    
  }
  
  
  
  #6.DIV7与E11关系
  {
    Divtrans = readRDS("sc_result/x.E11_DIV7/RetroDIV7_renamed_0904.RDS")
    v1ldiv = read.csv("sc_result/x.E11_DIV7/v1/final_scarform.csv",header = T) 
    v2ldiv = read.csv("sc_result/x.E11_DIV7/v2/final_scarform.csv",header = T) 
    celldiv = ReadCell("sc_result/x.E11_DIV7/DIV7_celltype.csv")
    
    prefix = c("v1.","v2.")
    bl = readRDS("final_result/x.blacklist.rds")
    v1ldiv = merge(v1ldiv[,c("Cell.BC","umim")],celldiv,by="Cell.BC");colnames(v1ldiv)[2] = "pattern"
    v2ldiv = merge(v2ldiv[,c("Cell.BC","umim")],celldiv,by="Cell.BC");colnames(v2ldiv)[2] = "pattern"
    v1ldivf = BlacklistFilter(v1ldiv,bl$v1_array)
    v2ldivf = BlacklistFilter(v2ldiv,bl$v2_array)
    datadiv = list(v1ldivf,v2ldivf)
    outnamediv = "/picb/sysgenomics2/people/liuhengxin/P6_lineartree/final_result/xE11_DIV7/E11_DIV"
    Divtree = TreeConstruct(datadiv,prefix,outnamediv,celldiv)
    Divtree = CirclePlot1(Divtree,outnamediv)
    Divtree = MyHeatmap(Divtree,outnamediv)
    saveRDS(Divtree,file = "final_result/xE11_DIV7/E11_div_tree.rds")
    
    # ggvenn(list("DIV7" = unique(Divtree$cm$tags),"E11" = unique(E11tree$cm$tags)),stroke_alpha = 0.5,stroke_size = 0.5)
    
    corderdiv = Divtree$distfigure$spearhc$labels[Divtree$distfigure$spearhc$order]
    
    count1 = Divtree$count
    count2 = E11tree$count
    count1 = count1[,c("tags",corderdiv)]
    count2 = count2[,c("tags",corder11)]
    
    colnames(count1)[-1] = paste0("DIV7_",colnames(count1)[-1])
    colnames(count2)[-1] = paste0("E11_",colnames(count2)[-1])
    e11divcms = merge(count2,count1,by = "tags")
    
    e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$tags;e11divcms.pd = e11divcms.pd[-1]
    e11divcms.pd = as.matrix(e11divcms.pd)
    e11divcms.pd[which(e11divcms.pd>0)] = 1
    
    split = factor(rep(c("E11","DIV7"),c(ncol(E11tree$count[-1]),ncol(Divtree$count[-1]))),levels = c("E11","DIV7"))
    names(split) = colnames(e11divcms.pd)
    
    e11divcms.pd = e11divcms.pd[do.call(order, as.data.frame(-e11divcms.pd)), ]
    

    tmp = e11divcms.pd[do.call(order, as.data.frame(-e11divcms.pd)), ]
    
    pdf("final_result/xE11_DIV7/E11-DIV7_array_cmp_heatmap.pdf",width = 6,height = 10)
    Heatmap(log2(e11divcms.pd+1),show_row_names = F,
            # col = c("#23174C","#376BB4","white"),
            col = colorRampPalette(brewer.pal(9, "Blues"))(30),
            column_gap = unit(5, "mm"),
            cluster_row_slices = FALSE,row_order = rownames(tmp),
            column_split = split,column_order = colnames(tmp))
    dev.off()
    
    
    #量化分析
    count1 = Divtree$count
    count2 = E11tree$count
    count1 = count1[,c("tags",corderdiv)]
    count2 = count2[,c("tags",corder11)]
    #方法1
    count2f = count2[which(rowSums(count2[-1]>0) == 1),]
    cmt1 = melt(count1)
    cmt2 = melt(count2f)
    e11divcmf = merge(cmt2,cmt1,by = c("tags"),suffixes = c(".e11",".div"))
    e11divcmf = e11divcmf[which(e11divcmf$value.e11>0 & e11divcmf$value.div > 0),]
    ggplot(e11divcmf,aes(y = variable.e11,x = variable.div,
                         size = value.div,color = value.e11)) + 
      geom_point() + theme_bw()
    
    
    #方法2
    tmp = apply(count2[-1], 1, function(x){return(paste0(names(x[which(x>0)]),collapse = "-"))})
    count2f = count2;count2f$collapse = tmp;
    tmp = table(tmp);tmp = tmp[which(tmp>5)];tmp = tmp[order(tmp)]
    count2f = count2f[which(count2f$collapse %in% names(tmp)),]
    count2f$clone.counts = as.vector(tmp[match(count2f$collapse,names(tmp))])
    
    cmt2 = count2f[,c("tags","collapse","clone.counts")]
    e11divcmf = merge(cmt2,cmt1,by = c("tags"),suffixes = c(".e11",".div"))
    e11divcmf = e11divcmf[which(e11divcmf$value > 0),]
    e11divcmf$collapse = factor(e11divcmf$collapse,levels = names(tmp))
    p6.1 = ggplot(e11divcmf,aes(y = collapse,x = variable, color = clone.counts,
                         size = value)) + 
      scale_color_gradient(low = "#30E8BF", high = "#FF8235")+
      geom_point() + theme_bw()
    p6.1
    ggsave(p6.1,filename = "final_result/xE11_DIV7/E11_DIV_collapse_clone_point.pdf",width = 10,height = 6)
    
    #尝试 heatmap
    cmt1 = melt(count1)
    cmt2 = melt(count2)
    #方法1
    e11divcmf.mx = cor(as.matrix(e11divcms[-1]),method = "spearman")
    e11divcmf.mx = e11divcmf.mx[1:(ncol(count2)-1),(ncol(count2)):ncol(e11divcmf.mx)]
    pdf("final_result/xE11_DIV7/E11_DIV_cor_spearman_heatmap.pdf",width = 6,height = 5)
    Heatmap(e11divcmf.mx,
            clustering_method_rows = "complete")
    dev.off()
    
    #方法2
    e11divcmf = merge(cmt2,cmt1,by = c("tags"),suffixes = c(".e11",".div"))
    e11divcmf = e11divcmf[which(e11divcmf$value.e11>0 & e11divcmf$value.div > 0),]
    
    e11divcmf = e11divcmf %>% group_by(tags) %>% summarise(e11.cell = variable.e11,
                                                           div.cell = variable.div,
                                                           e11.size = value.e11/sum(value.e11),
                                                           div.size = value.div/sum(value.div))
    e11divcmf$weight = (e11divcmf$e11.size + e11divcmf$div.size)/2
    e11divcmf.gp = e11divcmf %>% group_by(e11.cell,div.cell) %>% summarise(weight = mean(weight))
    e11divcmf.gp = dcast(e11divcmf.gp,e11.cell~div.cell,value.var = "weight")
    e11divcmf.gp[is.na(e11divcmf.gp)] = 0
    e11divcmf.gp = melt(e11divcmf.gp)
    colnames(e11divcmf.gp) = c("DIV","E11","weight")
    
    ggplot(e11divcmf.gp,aes(x = DIV, y = E11,fill = weight)) + 
      geom_tile() +
      scale_fill_gradient(low = "white",high = "#005AA7") +
      theme(axis.ticks = element_blank(),panel.background = element_blank())
    
    #判断出三块有意义的区域，尝试画一个分块的桑基图
    
    # library(ggsankey)
    # 
    # df <- mtcars %>%
    #   make_long(cyl, vs, am, gear, carb)
    # df <- e11divcmf[c(2,3)] %>%
    #   make_long(e11.cell,div.cell)
    # 
    # 
    # ggplot(df, aes(x = x, 
    #                next_x = next_x, 
    #                node = node, 
    #                next_node = next_node,
    #                fill = factor(node))) +
    #   geom_sankey()
    # 
    # ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), 
    #                label = node)) +
    #   geom_sankey(flow.alpha = .6,
    #               node.color = "gray30") +
    #   geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    #   scale_fill_viridis_d() +
    #   theme_sankey(base_size = 18) +
    #   labs(x = NULL) +
    #   theme(legend.position = "none",
    #         plot.title = element_text(hjust = .5)) +
    #   ggtitle("Car features")
    
    
    
    #尝试计算每个BC的富集值，利用富集值来做桑基图
    {
      cellmx.e = count2
      cellmx.l = count1
      cmtnp.e = melt(cellmx.e)
      cmtnp.l = melt(cellmx.l)
      colnames(cmtnp.e) = c("tags","celltype","counts")
      colnames(cmtnp.l) = c("tags","celltype","counts")
      
      Me = cmtnp.e %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Ne = sum(Me$counts)
      
      Ml = cmtnp.l %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Nl = sum(Ml$counts)
      
      cellft.e = NULL
      for (i in 1:nrow(cellmx.e)) {
        line = cellmx.e[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Me[which(Me$celltype == names(k)),]$counts
            N = Ne
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.e = rbind(cellft.e,ftl)
      }
      
      
      cellft.l = NULL
      for (i in 1:nrow(cellmx.l)) {
        line = cellmx.l[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Ml[which(Ml$celltype == names(k)),]$counts
            N = Nl
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.l = rbind(cellft.l,ftl)
      }
      cellft = merge(cellft.e,cellft.l,by = "tags",suffixes = c(".E11",".DIV"))
      
    }
    cellft = cellft[which(rowSums(cellft[-1]>0) < ncol(cellft[-1])/3),]
    #绘制桑基图
    {
      library(ggalluvial)
      
      colnames(cellft)
      colnames(cellft.e)
      earlyid = colnames(cellft)[-1][1:ncol(cellft.e)]
      lateid = setdiff(colnames(cellft)[-1],earlyid)
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
      
      domain_df = data.frame("celltype" = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL",
                                            "GluAL","GabaProgBI","NProgBM","OMTN","NbBM0","NbBM1",
                                            "Rgl1.E11","NbFP","Peric",
                                            "Gaba3_Nkx2-2","Gaba4_Six3","Glu_Dorsal",
                                            "Glu_RN",
                                            "DA","Rgl2","Glu_M","NbM","Rgl1.DIV",
                                            "GABA_unknown"),
                             "domain" = rep(c("c1","c2","c3","c1","c2","c3","unknown"),c(4,6,3,3,1,5,1)))
      domain_df = domain_df[order(domain_df$domain),]
      
      phasestf = phasest[which(phasest$phase2 != "GABA_unknown"),]
      # phasestf = phasest[which(phasest$phase1 > 0.001 | phasest$phase2 > 0.001),]
      # phasestf = phasestf[which(phasestf$Freq>20),]
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
      phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domain),]
      phasest_lod$celltype = factor(phasest_lod$celltype,
                                    levels = unique(as.character(phasest_lod$celltype)))
      
      
      pf1 = ggplot(phasest_lod,
                   aes(x = phase, stratum = celltype, alluvium = alluvium,
                       y = Freq,
                       fill = domain, label = celltype)) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_flow(aes(color = celltype)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      pf1
      
      ggsave(pf1,filename = "final_result/xE11_DIV7/E11_DIV_class_allu_plot.pdf",width = 8,height = 8)
      
      library(circlize)
      library(randomcoloR)
      require(RColorBrewer)
      brewer.pal(9, "Set1")
      
      phasestm = phasestf
      phasestm$Freq.from = phasestm$Freq.from/sum(phasestm$Freq.from)
      phasestm$Freq.to = phasestm$Freq.to/sum(phasestm$Freq.to)
      
      phasestm = merge(phasestm,domain_df,
                       by.x = "phase1",by.y = "celltype")
      phasestm = phasestm[which(phasestm$phase1 > 0.001 | phasestm$phase2 > 0.001),]
      gridcol = data.frame(celltype = unique(c(phasestm$phase1,phasestm$phase2)),
                           phase = c(rep(1,length(unique(phasestm$phase1))),
                                     rep(2,length(unique(phasestm$phase2)))))
      gridcol$domain = domain_df[match(gridcol$celltype,domain_df$celltype),"domain"]
      domaincol = data.frame(domain = unique(gridcol$domain),
                             col = distinctColorPalette(length(unique(gridcol$domain))))
      gridcol$col = domaincol[match(gridcol$domain,domaincol$domain),"col"]
      gridcol = gridcol[order(-gridcol$phase,gridcol$domain),]
      tmp = gridcol[which(gridcol$phase == 1),]
      gridcol[which(gridcol$phase == 1),] = tmp[order(tmp$domain,decreasing = T),]
      gridcolv = gridcol$col;names(gridcolv) = gridcol$celltype
      
      pdf(paste0("final_result/xE11_DIV7/E11_DIV_class_chord_plot.pdf"))
      circos.par(start.degree = 90)
      chordDiagram(phasestm[,1:4],annotationTrack =  c("name", "grid"),
                   grid.col = gridcolv,
                   directional = 1, direction.type = c("diffHeight", "arrows"),
                   link.arr.type = "big.arrow",
                   order = unique(gridcol$celltype))
      # Restart circular layout parameters
      circos.clear()
      dev.off()
      
    }
    
    
  }
  
  #7.分化诊断模型
  {
    
    #输入数据集
    # write.csv(DA.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_DA_Rgl1.csv",row.names = T,quote = F)
    # write.csv(Glu.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_Glu_Rgl1.csv",row.names = T,quote = F)
    # write.csv(All.markers,"final_result/xE15.5/DA_Glu_direction_DEG_analysis_All_Rgl1.csv",row.names = T,quote = F)
    # saveRDS(E15trans.hub,file = "final_result/xE15.5/E15trans_DA.rds")
    # 
    library(pROC)
    FindVariableFeatures(E15trans.hub, selection.method = "vst", nfeatures = 2000)
    
    
    DA.markers = read.csv("final_result/xE15.5/DA_Glu_direction_DEG_analysis_DA.csv")
    Glu.markers = read.csv("final_result/xE15.5/DA_Glu_direction_DEG_analysis_Glu.csv")
    All.markers = read.csv("final_result/xE15.5/DA_Glu_direction_DEG_analysis_All.csv")
    trans.da = subset(E15trans.hub, idents = "DA_d")
    trans.glu = subset(E15trans.hub, idents = "Glu_d")
    trans.two = subset(E15trans.hub, idents = "All_d")
    damx = as.data.frame(trans.da[["RNA"]]@counts)
    glumx = as.data.frame(trans.glu[["RNA"]]@counts)
    twomx = as.data.frame(trans.two[["RNA"]]@counts)
    DA.markers = DA.markers[which(abs(DA.markers$avg_log2FC) > 0.5),]
    Glu.markers = Glu.markers[which(abs(Glu.markers$avg_log2FC) > 0.5),]
    All.markers = All.markers[which(abs(All.markers$avg_log2FC) > 0.5),]
    
    markers = unique(c(DA.markers$X,Glu.markers$X))
    markers <- head(VariableFeatures(E15trans.hub), 2000)
    nbbmid = cell[which(cell$Cell.type == "NbFP"), "Cell.BC"]
    
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
      dadata = dadata[which(substr(rownames(dadata),1,16) %in% nbbmid),]
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
    
    tmp = colnames(E15trans.hub)[which(substr(colnames(E15trans.hub),1,16)%in% nbbmid)]
    E15trans.NBFP = subset(E15trans.hub,cells = tmp)
    FindVariableFeatures(E15trans.NBFP, selection.method = "vst", nfeatures = 2000)
    
    
    marker1 = head(VariableFeatures(E15trans.NBFP), 256)
    marker2 = unique(c(DA.markers$X,Glu.markers$X))
    marker3 = rownames(E15trans.NBFP)[sample(1:length(rownames(E15trans.NBFP)),256)]
    marker4 = read.delim("raw_data/tf.name.txt",header = F)
    marker4 = marker4$V2
    marker4 = head(VariableFeatures(E15trans.NBFP)[VariableFeatures(E15trans.NBFP) %in% marker4],256)
    
    aucdf = NULL
    for (i in 1:100) {
      roc1 = NeurBuild(marker1)
      roc2 = NeurBuild(marker2)
      roc3 = NeurBuild(marker3)
      roc4 = NeurBuild(marker4)
      aucdf = rbind(aucdf,data.frame(sample = i, markers = c("Highly variable genes(n=256)",
                                                             "Differentially expressed(n=256)",
                                                             "Random genes(n=256)",
                                                             "Transcription factors(n=256)"),
                                     AUC = c(roc1$auc,roc2$auc,roc3$auc,roc4$auc)))
    }
    
    aucdf$markers = factor(aucdf$markers,levels = c("Differentially expressed(n=256)",
                                                    "Highly variable genes(n=256)",
                                                    "Random genes(n=256)",
                                                    "Transcription factors(n=256)"))
    p7 = ggboxplot(aucdf, x = "markers", y = "AUC",
            fill = "markers", palette = "jco",
            width = 0.2) + 
        theme_bw() + theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           legend.position = "none",
                           axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      stat_compare_means(method = "anova", label.y = 1) +
      stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = "Differentially expressed(n=256)")    
    p7
    ggsave(p7,filename = "final_result/xE15.5/DA_Glu_neuralnet_auc_boxplot.pdf",width = 4,height = 4)
      
    aucdf %>% group_by(markers) %>% summarise(mean(AUC))  
    write.csv(aucdf,file = "final_result/xE15.5/DA_Glu_neuralnet_auc.csv",row.names = F,quote = F)
    write.csv(dadata,file = "final_result/xE15.5/DA_Glu_neuralnet_NbFP_dataset.csv",row.names = T,quote = F)
    
    
    plot.roc(modelroc, print.auc=TRUE, auc.polygon=TRUE, add=FALSE,
             grid.col=c("green", "red"), max.auc.polygon=TRUE,asp = NA,
             auc.polygon.col="skyblue", print.thres=TRUE)
    
    
    # min(test.hat) + (max(test.hat) - min(test.hat))*0.6
    test.predict2 = apply(test.hat, 1, function(x){ifelse(x>0.298,1,0)})
    predict.table = table(testdata$type,test.predict2)
    predict.table
    
    model_accuracy <- sum(diag(predict.table))/sum(predict.table)
    model_accuracy
    
    
    pdf("final_result/xtotal/DA_Glu_predictor.pdf",width = 5,height = 4)
    
    dev.off()
    
    
    
  }
  
  #8.bulk indel occurrence histogram
  {
    v1lsx[[5]] = rbind(v1lsx[[5]],v1lsx[[6]])
    v1lsx[[7]] = rbind(v1lsx[[7]],v1lsx[[8]])
    v1lsx[[9]] = rbind(v1lsx[[9]],v1lsx[[10]])
    names(v1lsx)[c(5,7,9)] = c("E15-1","E15-2","E15-3")
    v1lsx = v1lsx[-c(6,8,10)]
    
    v2lsx[[5]] = rbind(v2lsx[[5]],v2lsx[[6]])
    v2lsx[[7]] = rbind(v2lsx[[7]],v2lsx[[8]])
    v2lsx[[9]] = rbind(v2lsx[[9]],v2lsx[[10]])
    names(v2lsx)[c(5,7,9)] = c("E15-1","E15-2","E15-3")
    v2lsx = v2lsx[-c(6,8,10)]
    
    
    id1 = NULL
    id2 = NULL
    for (i in 2:length(v1lsx)) {
      line1 = v1lsx[[i]]
      line2 = v2lsx[[i]]
      line.id1 = unique(unlist(strsplit(unique(line1$umim),split = "_|&")))
      line.id2 = unique(unlist(strsplit(unique(line2$umim),split = "_|&")))
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
    
    p1 = ggplot(idocv,aes(x = Var1, y = norm, fill = group)) + 
      geom_bar(stat = "identity",position = "dodge",width = 0.7) +
      scale_fill_manual(values=c('lightpink1','lightblue2')) +
      theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
      scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))
    print(p1)
    
    ggsave(p1,filename = "final_result/xBulk/fig2.I.indel_occrurency_in_sample.pdf",width = 7,height = 5)
    write.csv(idocv,file = "final_result/xBulk/fig2.I.indel_occrurency_in_sample.csv",row.names = F,quote = F)
    
  }
  
  
  
}

#9. 12.02补图
{
  #加一个v1 + v2 多样性图
  head(v1testf)
  #E15
  scdv = data.frame(sample = rep(c("E11","E15"),c(3,3)),
                    group = rep(c("v1","v2","v1&v2"),2),
                    diversity = c(length(unique(v1testf11$pattern)),
                                  length(unique(v2testf11$pattern)),
                                  length(unique(E11tree$cm$tags)),
                                  length(unique(v1testf$pattern)),
                                  length(unique(v2testf$pattern)),
                                  length(unique(E15tree$cm$tags))))
  
  #查看barcode表达量
  samplels = list.files("sc_result/X.E11",full.names = T)[c(5:7)]
  v1lse11 = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lse11 = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lse11) = list.files("sc_result/X.E11/")[c(5:7)]
  names(v2lse11) = list.files("sc_result/X.E11/")[c(5:7)]
  head(v1lse11$lib1)
  e11rna = as.data.frame(E11trans@assays$RNA@counts)
  e11rna = as.data.frame(t(e11rna))
  e11rna[1:5,1:5]
  e11rna$Cell.BC = substr(rownames(e11rna),1,16) 
  
  
  aycount = merge(v1lse11$lib1[c("Cell.BC","umi_num")],e11rna,by = "Cell.BC")
  head(aycount)
  aystat = NULL
  for (i in 1:1000) {
    lineay = unlist(aycount[i,-1])
    lineay = lineay[which(lineay>0)]
    lineay = lineay[order(-lineay)]
    print(i)
    aystat = c(aystat,which(names(lineay) == "umi_num"))
  }
  gghistogram(aystat,bins = 100) + scale_x_continuous(limits = c(0,100))
  aystat = as.data.frame(aystat)
  ggplot(aystat, aes(log10(aystat))) + stat_ecdf(geom = "step") + 
    geom_vline(xintercept=c(2,2.6),lty=3,col="black",lwd=0.5) + scale_y_continuous(breaks = c(seq(0,1,0.1)))
  
  max(aystat)
  
}

#12.20补图 E11 DIV新数据
{
  corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
               "NbBM1","Rgl1","NbFP_0")
  domain11[12,1] = "NbFP_0"
  #E11 单独重新作图
  cellallu11 = BCEnrichment(E11Stree,domain11,earlyid11,"final_result/STF_E11/singlet_alluv11")
  cellallu11$figure.allu
  saveRDS(cellallu11,file = "final_result/STF_E11/singlet_alluv11.rds")
  ggsave(cellallu11$figure.allu,filename = "final_result/STF_E11/singlet_alluv11_plot.pdf")
  
  
  # cellalludiv = BCEnrichment(DIVtree,domain11,earlyid11,"final_result/STF_DIV7/alluv11")
  # cellallu11$figure.allu
  pvs11 = PseudoDisper(E11Stree,domain11)
  pvs11$figure_pvs
  pvs11$figure_pvs_box
  ggsave(pvs11$figure_pvs,filename = "final_result/STF_E11/singlet_E11_pseudo_disper_cmb.pdf",width = 8,height = 6)
  ggsave(pvs11$figure_pvs_box,filename = "final_result/STF_E11/singlet_E11_pseudo_disper_box.pdf",width = 5,height = 3)
  
  #domain heatmap
  ardm = pvs11$pvs_domain
  ardm = ardm[-c(1:6,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  
  mx = cor(ardm,method = "pearson")
  pdf("final_result/STF_E11/singlet_E11_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               row_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  #E11 clone 聚类图------------------
  #细胞聚类
  {
    cmtnp = E11Stree$cm %>% group_by(tags,celltype) %>% 
      summarise(counts = sum(n()))
    #算一个total BC 富集值
    cellmx = dcast(cmtnp[,c("tags","celltype","counts")],
                   tags~celltype,value.var = "counts",fun.aggregate = sum)
    Mt = cmtnp %>% group_by(celltype) %>% summarise(counts = sum(counts))
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
    cellftl = melt(cellft)
    head(cellftl)
    
    #尝试聚类分析
    rownames(cellft) = cellft$tags;cellft = cellft[-1];cellftn = scale(cellft)
    
    tmp = table(cellftl[which(cellftl$value>0),]$tags)
    tmp = tmp[which(tmp>1)]
    cellftf = cellft[which(rownames(cellft)%in% names(tmp)),]
    cellmxf = cellmx[-1];rownames(cellmxf) = cellmx$tags;cellmxf = cellmxf[which(rownames(cellmxf)%in% names(tmp)),]
    
    library(ggfortify)
    library(ggplot2)
    
    colnames(cellftf)
    
    cellorder = corder11n
    
    n = 15
    
    # fit = kmeans(cellftf,n)
    # split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    # reorder.hmap <- Heatmap(as.matrix(cellftf), split=split, cluster_row_slices = FALSE,
    #                         show_row_names = F,column_order = cellorder)
    # reorder.hmap
    
    # ggexport(reorder.hmap, filename = "final_result/xE15.5/tag_cluster_heatmap.pdf",width = 8,height = 8)
    cellmxf = cellmxf[,colnames(cellmxf) != "Peri"]
    cellmxfn2 = t(t(as.matrix(cellmxf))/colSums(cellmxf))
    cellmxfn = as.matrix(cellmxf)
    cellmxfn[which(cellmxfn>0)] = 1
    fit = kmeans(cellmxfn,n)
    split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    reorder.hmap <- Heatmap(as.matrix(cellmxfn), split=split, cluster_row_slices = FALSE,
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    #merge cluster
    clist = list(c1 = c(2,11),
                 c2 = c(12),
                 c3 = c(14),
                 c4 = c(10),
                 c5 = c(6),
                 c6 = c(13),
                 c7 = c(5),
                 c8 = c(7),
                 c9 = c(1),
                 c10 = c(4),
                 c11 = c(8),
                 c12 = c(9),
                 c13 = c(3),
                 c14 = c(15))
    # c10 = c(29),)
    split <- fit$cluster
    for (i in 1:length(clist)) {
      split[which(split %in% clist[[i]])] = names(clist)[i]
    }
    
    split <- factor(as.character(split), levels=as.character(c(paste0("c",1:length(clist)) )))      
    names(split) = names(fit$cluster)
    # e11_clu = readRDS("final_result/xE11/E11_cellcounts_cluster.rds")
    reorder.hmap <- Heatmap(as.matrix(cellmxfn2), split=split, cluster_row_slices = FALSE,
                            col =  colorRamp2(c(0, 0.01, max(cellmxfn2)), c("#08306B", "white", "red")),
                            name = "normlize enrichment score",
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    ggexport(reorder.hmap, filename = "final_result/STF_E11/singlet_E11_tag_cluster_heatmap_counts_arti.pdf",width = 8,height = 8)
    
    saveRDS(list(cluster = split,fit = fit,cellcounts = cellmxfn2,cluster.fig = reorder.hmap),
            file = "final_result/STF_E11/singlet_E11_cellcounts_cluster.rds")
    
    #cluster show in dimplot
    cclsSE11 = readRDS("final_result/STF_E11/singlet_E11_cellcounts_cluster.rds")
    reorder.hmap <- Heatmap(as.matrix(cclsSE11$cellcounts), split=cclsSE11$cluster, cluster_row_slices = FALSE,
                            col =  colorRamp2(c(min(xm), 0, max(xm)), c("#08306B", "white", "Reds")),
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    
    ccidents = as.data.frame(cclsSE11$cluster)
    ccidents = merge(ccidents,E11Stree$cm,by.x = "row.names",by.y = "tags")
    colnames(ccidents) = c("tags","cluster","CellBC","Freq","Celltype")
    ccidents$Celltype = factor(ccidents$Celltype,levels = cellorder)
    write.csv(ccidents,file = "final_result/STF_E11/singlet_E11_tag_cluster_BC.csv",row.names = F,quote = F)
    
    ccidents = ccidents[which(!is.na(ccidents$Celltype)),]
    p3.2 = ggplot(ccidents,aes(x = cluster,fill = Celltype)) +
      geom_bar(stat = "count",position = "fill",width = 0.7) +
      theme_bw()
    p3.2
    ggsave(p3.2,filename = "final_result/STF_E11/singlet_E11_tag_cluster_celltype_bar.pdf",width = 8,height = 6)
    
    p3.3 = list()
    for (i in 1:length(clist)) {
      cc = ccidents[which(ccidents$cluster == paste0("c",i)),]
      tarcell = colnames(E11trans)[which(substr(colnames(E11trans),1,16) %in% cc$CellBC)]
      pt = DimPlot(E11trans,label = TRUE, label.size = 3,cells.highlight = tarcell, 
                   cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
      p3.3[[i]] = pt
      names(p3.3)[i] = paste0("c",i)
      
    }
    p3.3.0 = ggarrange(plotlist =  p3.3,labels = names(p3.3))
    p3.3.0
    
    pdf("final_result/STF_E11/singlet_E11_tag_cluster_umap.pdf",width = 10,height = 8)
    print(p3.3.0)
    invisible(lapply(p3.3, print))
    dev.off()
    
  }
  
  #DIV7与E11关系
  {
    corderdiv = DIVtree$distfigure$abspearsonhc$labels[DIVtree$distfigure$abspearsonhc$order]
    corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
                  "NbBM1","Rgl1","NbFP_0","Peri")
    
    count1 = DIVtree$count
    count2 = E11Stree$count
    
    count1 = count1[,c("tags",corderdiv)]
    count2 = count2[,c("tags",corder11n)]
    
    colnames(count1)[-1] = paste0("DIV7_",colnames(count1)[-1])
    colnames(count2)[-1] = paste0("E11_",colnames(count2)[-1])
    e11divcms = merge(count2,count1,by = "tags")
    
    e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$tags;e11divcms.pd = e11divcms.pd[-1]
    e11divcms.pd = as.matrix(e11divcms.pd)
    e11divcms.pdt = e11divcms.pd
    e11divcms.pdt[which(e11divcms.pdt>0)] = 1
    
    split = factor(rep(c("E11","DIV7"),c(ncol(E11Stree$count[-1]),ncol(DIVtree$count[-1]))),levels = c("E11","DIV7"))
    names(split) = colnames(e11divcms.pd)
    
    e11divcms.pd = e11divcms.pd[do.call(order, as.data.frame(-e11divcms.pdt)), ]
    
    tmp = e11divcms.pdt[do.call(order, as.data.frame(-e11divcms.pdt)), ]
    # e11divcms.pd = t(t(e11divcms.pd)/colSums(e11divcms.pd))
    # e11divcms.pd[is.nan(e11divcms.pd)] = 0
    
    pdf("final_result/STF_DIV7/singlet_E11-DIV7_array_cmp_heatmap_log2counts.pdf",width = 6,height = 8)
    Heatmap(log2(e11divcms.pd+1),show_row_names = F,
            # col = c("#23174C","#376BB4","white"),
            col = colorRampPalette(brewer.pal(9, "Blues"))(30),
            column_gap = unit(5, "mm"),name = "log2(counts)",
            cluster_row_slices = FALSE,row_order = rownames(tmp),
            column_split = split,column_order = colnames(tmp))
    dev.off()
    
    
    # #量化分析
    # count1 = Divtree$count
    # count2 = E11Stree$count
    # count1 = count1[,c("tags",corderdiv)]
    # count2 = count2[,c("tags",corder11)]
    #方法1
    # count2f = count2[which(rowSums(count2[-1]>0) == 1),]
    # cmt1 = melt(count1)
    # cmt2 = melt(count2f)
    # e11divcmf = merge(cmt2,cmt1,by = c("tags"),suffixes = c(".e11",".div"))
    # e11divcmf = e11divcmf[which(e11divcmf$value.e11>0 & e11divcmf$value.div > 0),]
    # ggplot(e11divcmf,aes(y = variable.e11,x = variable.div,
    #                      size = value.div,color = value.e11)) + 
    #   geom_point() + theme_bw()
    
    
    #方法2
    tmp = apply(count2[-1], 1, function(x){return(paste0(names(x[which(x>0)]),collapse = "-"))})
    count2f = count2;count2f$collapse = tmp;
    tmp = table(tmp);tmp = tmp[which(tmp>5)];tmp = tmp[order(tmp)]
    count2f = count2f[which(count2f$collapse %in% names(tmp)),]
    count2f$clone.counts = as.vector(tmp[match(count2f$collapse,names(tmp))])
    
    cmt2 = count2f[,c("tags","collapse","clone.counts")]
    e11divcmf = merge(cmt2,cmt1,by = c("tags"),suffixes = c(".e11",".div"))
    e11divcmf = e11divcmf[which(e11divcmf$value > 0),]
    e11divcmf$collapse = factor(e11divcmf$collapse,levels = names(tmp))
    p6.1 = ggplot(e11divcmf,aes(y = collapse,x = variable, color = clone.counts,
                                size = value)) + 
      scale_color_gradient(low = "#30E8BF", high = "#FF8235")+
      geom_point() + theme_bw()
    p6.1
    ggsave(p6.1,filename = "final_result/xE11_DIV7/E11_DIV_collapse_clone_point.pdf",width = 10,height = 6)
    
    #尝试 heatmap
    cmt1 = melt(count1)
    cmt2 = melt(count2)
    #方法1
    e11divcmf.mx = cor(as.matrix(e11divcms[-1]),method = "spearman")
    e11divcmf.mx = e11divcmf.mx[1:(ncol(count2)-1),(ncol(count2)):ncol(e11divcmf.mx)]
    e11divcmf.mx[is.na(e11divcmf.mx)] = 0
    pdf("final_result/STF_DIV7/singlet_E11_DIV_cor_spearman_heatmap.pdf",width = 6,height = 5)
    Heatmap(e11divcmf.mx,
            clustering_method_rows = "complete")
    dev.off()
    
    # count1 = count1[which(rowSums(count1[-1]>0) == 1),]
    count2 = count2[which(rowSums(count2[-1]>0) == 1),]
    e11divcms = merge(count2,count1,by = "tags")
    
    
    #尝试计算每个BC的富集值，利用富集值来做桑基图
    {
      cellmx.e = count2
      cellmx.l = count1
      cmtnp.e = melt(cellmx.e)
      cmtnp.l = melt(cellmx.l)
      colnames(cmtnp.e) = c("tags","celltype","counts")
      colnames(cmtnp.l) = c("tags","celltype","counts")
      
      Me = cmtnp.e %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Ne = sum(Me$counts)
      
      Ml = cmtnp.l %>% group_by(celltype) %>% summarise(counts = sum(counts))
      Nl = sum(Ml$counts)
      
      cellft.e = NULL
      for (i in 1:nrow(cellmx.e)) {
        line = cellmx.e[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Me[which(Me$celltype == names(k)),]$counts
            N = Ne
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.e = rbind(cellft.e,ftl)
      }
      
      
      cellft.l = NULL
      for (i in 1:nrow(cellmx.l)) {
        line = cellmx.l[i,]
        n = sum(line[2:ncol(line)])
        ftl = line
        for (j in 2:ncol(line)) {
          k = unlist(line[j])
          if(k>0){
            M = Ml[which(Ml$celltype == names(k)),]$counts
            N = Nl
            d = data.frame(cell.not.interest=c(M-k, N-M-n+k), cell.in.interest=c(k, n-k))
            row.names(d) <- c("In_category", "not_in_category")
            ft = fisher.test(d)
            ftl[j] = -log10(ft$p.value)
          }
        }
        cellft.l = rbind(cellft.l,ftl)
      }
      cellft = merge(cellft.e,cellft.l,by = "tags",suffixes = c(".E11",".DIV"))
      
      }
    cellft = cellft[which(rowSums(cellft[-1]>0) < ncol(cellft[-1])/3),]
    #绘制桑基图
    {
      library(ggalluvial)
      
      colnames(cellft)
      colnames(cellft.e)
      earlyid = colnames(cellft)[2:ncol(cellft.e)]
      lateid = setdiff(colnames(cellft)[-1],earlyid)
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
      
      domain_df = data.frame("celltype" = c("E11_NbGluAL","E11_NProgBL","E11_GabaProg0","E11_GabaProgBL",
                                            "E11_NbFP_0","E11_Peri","E11_Rgl1","E11_OMTN",
                                            "E11_GabaProgBI","E11_NbBM1","E11_NProgBM","E11_NbBM0",
                                            "DIV7_GabaBL1","DIV7_GabaBL4","DIV7_GabaBL2","DIV7_GluAL1","DIV7_GabaUnknown","DIV7_Unknown","DIV7_GabaBI",
                                            "DIV7_Rgl2","DIV7_Rgl3","DIV7_DA","DIV7_NbGlu","DIV7_NbDA","DIV7_NbFP",
                                            "DIV7_GluBM1","DIV7_GluBM2","DIV7_GluFP"),
                             "domain" = rep(c("c1","c2","c3","c1","c2","c3"),c(4,4,4,7,6,3)))
      domain_df = domain_df[order(domain_df$domain),]
      
      phasestf = phasest[which(phasest$phase1 > 0.001 | phasest$phase2 > 0.001),]
      # phasestf = phasestf[which(phasestf$Freq>20),]
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
      phasest_lod = phasest_lod[order(phasest_lod$phase,phasest_lod$domain),]
      phasest_lod$celltype = factor(phasest_lod$celltype,
                                    levels = unique(as.character(phasest_lod$celltype)))
      
      
      pf1 = ggplot(phasest_lod,
                   aes(x = phase, stratum = celltype, alluvium = alluvium,
                       y = Freq,
                       fill = domain, label = celltype)) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_flow(aes(color = celltype)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      pf1
      
      ggsave(pf1,filename = "final_result/STF_DIV7/singlet_E11_DIV_class_allu_plot.pdf",width = 8,height = 8)
      
      library(circlize)
      library(randomcoloR)
      require(RColorBrewer)
      brewer.pal(9, "Set1")
      
      phasestm = phasestf
      phasestm$Freq.from = phasestm$Freq.from/sum(phasestm$Freq.from)
      phasestm$Freq.to = phasestm$Freq.to/sum(phasestm$Freq.to)
      
      phasestm = merge(phasestm,domain_df,
                       by.x = "phase1",by.y = "celltype")
      phasestm = phasestm[which(phasestm$phase1 > 0.001 | phasestm$phase2 > 0.001),]
      gridcol = data.frame(celltype = unique(c(phasestm$phase1,phasestm$phase2)),
                           phase = c(rep(1,length(unique(phasestm$phase1))),
                                     rep(2,length(unique(phasestm$phase2)))))
      gridcol$domain = domain_df[match(gridcol$celltype,domain_df$celltype),"domain"]
      domaincol = data.frame(domain = unique(gridcol$domain),
                             col = distinctColorPalette(length(unique(gridcol$domain))))
      gridcol$col = domaincol[match(gridcol$domain,domaincol$domain),"col"]
      gridcol = gridcol[order(-gridcol$phase,gridcol$domain),]
      tmp = gridcol[which(gridcol$phase == 1),]
      gridcol[which(gridcol$phase == 1),] = tmp[order(tmp$domain,decreasing = T),]
      gridcolv = gridcol$col;names(gridcolv) = gridcol$celltype
      
      pdf(paste0("final_result/STF_DIV7/singlet_E11_DIV_class_chord_plot.pdf"))
      circos.par(start.degree = 90)
      chordDiagram(phasestm[,1:4],annotationTrack =  c("name", "grid"),
                   grid.col = gridcolv,
                   directional = 1, direction.type = c("diffHeight", "arrows"),
                   link.arr.type = "big.arrow",
                   order = unique(gridcol$celltype))
      # Restart circular layout parameters
      circos.clear()
      dev.off()
      
    }
    
    
  }
  
  
  
  #old E11 单独重新作图-----------
  cellallu11 = BCEnrichment(E11tree,domain11,earlyid11,"final_result/xE11/E11_12_21_alluv11")
  cellallu11$figure.allu
  saveRDS(cellallu11,file = "final_result/xE11_noDoublet/E11_12_21_alluv11.rds")
  ggsave(cellallu11$figure.allu,filename = "final_result/xE11_noDoublet/E11_12_21_alluv11_plot.pdf")
  
  
  # cellalludiv = BCEnrichment(DIVtree,domain11,earlyid11,"final_result/STF_DIV7/alluv11")
  # cellallu11$figure.allu
  pvs11 = PseudoDisper(E11tree,domain11)
  pvs11$figure_pvs
  pvs11$figure_pvs_box
  ggsave(pvs11$figure_pvs,filename = "final_result/xE11_noDoublet/E11_12_21_pseudo_disper_cmb.pdf",width = 8,height = 6)
  ggsave(pvs11$figure_pvs_box,filename = "final_result/xE11_noDoublet/E11_12_21_pseudo_disper_box.pdf",width = 5,height = 3)
  
  #domain heatmap
  ardm = pvs11$pvs_domain
  ardm = ardm[-c(1:6,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xE11_noDoublet/E11_12_21_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               row_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  #5.E11 clone 聚类图------------------
  #细胞聚类
  {
    cmtnp = E11tree$cm %>% group_by(tags,celltype) %>% 
      summarise(counts = sum(n()))
    #算一个total BC 富集值
    cellmx = dcast(cmtnp[,c("tags","celltype","counts")],
                   tags~celltype,value.var = "counts",fun.aggregate = sum)
    Mt = cmtnp %>% group_by(celltype) %>% summarise(counts = sum(counts))
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
    cellftl = melt(cellft)
    head(cellftl)
    
    #尝试聚类分析
    rownames(cellft) = cellft$tags;cellft = cellft[-1];cellftn = scale(cellft)
    
    tmp = table(cellftl[which(cellftl$value>0),]$tags)
    tmp = tmp[which(tmp>1)]
    cellftf = cellft[which(rownames(cellft)%in% names(tmp)),]
    cellmxf = cellmx[-1];rownames(cellmxf) = cellmx$tags;cellmxf = cellmxf[which(rownames(cellmxf)%in% names(tmp)),]
    
    library(ggfortify)
    library(ggplot2)
    
    colnames(cellftf)
    
    cellorder = corder11n[corder11n!="Peri"]
    
    n = 15
    
    # fit = kmeans(cellftf,n)
    # split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    # reorder.hmap <- Heatmap(as.matrix(cellftf), split=split, cluster_row_slices = FALSE,
    #                         show_row_names = F,column_order = cellorder)
    # reorder.hmap
    
    # ggexport(reorder.hmap, filename = "final_result/xE15.5/tag_cluster_heatmap.pdf",width = 8,height = 8)
    cellmxf = cellmxf[,colnames(cellmxf) != "Peri"]
    cellmxfn = as.matrix(cellmxf)
    cellmxfn[which(cellmxfn>0)] = 1
    cellmxfn2 = t(t(as.matrix(cellmxf))/colSums(cellmxf))
    fit = kmeans(cellmxfn,n)
    split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    reorder.hmap <- Heatmap(as.matrix(cellmxfn), split=split, cluster_row_slices = FALSE,
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    #merge cluster
    clist = list(c1 = c(14,8,4),
                 c2 = c(6),
                 c3 = c(15),
                 c4 = c(1),
                 c5 = c(3),
                 c6 = c(7),
                 c7 = c(12),
                 c8 = c(13),
                 c9 = c(10),
                 c10 = c(9),
                 c11 = c(5),
                 c12 = c(11),
                 c13 = c(2))
    # c10 = c(29),)
    split <- fit$cluster
    for (i in 1:length(clist)) {
      split[which(split %in% clist[[i]])] = names(clist)[i]
    }
    
    split <- factor(as.character(split), levels=as.character(c(paste0("c",1:length(clist)) )))      
    names(split) = names(fit$cluster)
    # e11_clu = readRDS("final_result/xE11_noDoublet/E11_cellcounts_cluster.rds")
    reorder.hmap <- Heatmap(as.matrix(cellmxfn2), split=split, cluster_row_slices = FALSE,
                            col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    ggexport(reorder.hmap, filename = "final_result/xE11_noDoublet/E11_12_21_tag_cluster_heatmap_counts_arti.pdf",width = 8,height = 8)
    
    saveRDS(list(cluster = split,fit = fit,cellcounts = cellmxfn,cluster.fig = reorder.hmap),
            file = "final_result/xE11_noDoublet/E11_12_21_cellcounts_cluster.rds")
    
    #cluster show in dimplot
    ccls = readRDS("final_result/xE11_noDoublet/E11_12_21_cellcounts_cluster.rds")
    
    ccidents = as.data.frame(ccls$cluster)
    ccidents = merge(ccidents,E11tree$cm,by.x = "row.names",by.y = "tags")
    colnames(ccidents) = c("tags","cluster","CellBC","Freq","Celltype")
    ccidents$Celltype = factor(ccidents$Celltype,levels = cellorder)
    write.csv(ccidents,file = "final_result/xE11_noDoublet/E11_12_21_tag_cluster_BC.csv",row.names = F,quote = F)
    
    ccidents = ccidents[which(!is.na(ccidents$Celltype)),]
    p3.2 = ggplot(ccidents,aes(x = cluster,fill = Celltype)) +
      geom_bar(stat = "count",position = "fill",width = 0.7) +
      theme_bw()
    p3.2
    ggsave(p3.2,filename = "final_result/xE11_noDoublet/E11_12_21_tag_cluster_celltype_bar.pdf",width = 8,height = 6)
    
    p3.3 = list()
    for (i in 1:length(clist)) {
      cc = ccidents[which(ccidents$cluster == paste0("c",i)),]
      tarcell = colnames(E11trans)[which(substr(colnames(E11trans),1,16) %in% cc$CellBC)]
      pt = DimPlot(E11trans,label = TRUE, label.size = 3,cells.highlight = tarcell, 
                   cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
      p3.3[[i]] = pt
      names(p3.3)[i] = paste0("c",i)
      
    }
    p3.3.0 = ggarrange(plotlist =  p3.3,labels = names(p3.3))
    p3.3.0
    
    pdf("final_result/xE11_noDoublet/E11_12_21_tag_cluster_umap.pdf",width = 10,height = 8)
    print(p3.3.0)
    invisible(lapply(p3.3, print))
    dev.off()
    
  }
  
  #封装聚类函数
  cellorder = corder11n
  tree = E11Stree
  outpath = "E11"
  trans = E11trans
  
  CloneCluster = function(tree, trans, cellorder, outpath){
    cmtnp = tree$cm %>% group_by(tags,celltype) %>% 
      summarise(counts = sum(n()))
    
    #calcluate total BC enrichment score
    cellmx = dcast(cmtnp[,c("tags","celltype","counts")],
                   tags~celltype,value.var = "counts",fun.aggregate = sum)
    Mt = cmtnp %>% group_by(celltype) %>% summarise(counts = sum(counts))
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
    cellftl = melt(cellft)
    head(cellftl)
    
    #first kmer cluster
    rownames(cellft) = cellft$tags;cellft = cellft[-1];cellftn = scale(cellft)
    tmp = table(cellftl[which(cellftl$value>0),]$tags)
    tmp = tmp[which(tmp>1)]
    cellftf = cellft[which(rownames(cellft)%in% names(tmp)),]
    cellmxf = cellmx[-1];rownames(cellmxf) = cellmx$tags;cellmxf = cellmxf[which(rownames(cellmxf)%in% names(tmp)),]
    
    library(ggfortify)
    library(ggplot2)
    
    colnames(cellftf)
    
    #clust to 15 segment
    n = 15

    # cellmxf = cellmxf[,colnames(cellmxf) != "Peri"]
    cellmxfn2 = t(t(as.matrix(cellmxf))/colSums(cellmxf))
    cellmxfn = as.matrix(cellmxf)
    cellmxfn[which(cellmxfn>0)] = 1
    fit = kmeans(cellmxfn,n)
    split <- factor(as.character(fit$cluster), levels=as.character(c(1:n)))
    reorder.hmap <- Heatmap(as.matrix(cellmxfn), split=split, cluster_row_slices = FALSE,
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    #merge cluster
    clist = list(c1 = c(2,11),
                 c2 = c(12),
                 c3 = c(14),
                 c4 = c(10),
                 c5 = c(6),
                 c6 = c(13),
                 c7 = c(5),
                 c8 = c(7),
                 c9 = c(1),
                 c10 = c(4),
                 c11 = c(8),
                 c12 = c(9),
                 c13 = c(3),
                 c14 = c(15))
    split <- fit$cluster
    for (i in 1:length(clist)) {
      split[which(split %in% clist[[i]])] = names(clist)[i]
    }
    
    split <- factor(as.character(split), levels=as.character(c(paste0("c",1:length(clist)) )))      
    names(split) = names(fit$cluster)
    # e11_clu = readRDS("final_result/xE11/E11_cellcounts_cluster.rds")
    reorder.hmap <- Heatmap(as.matrix(cellmxfn2), split=split, cluster_row_slices = FALSE,
                            col =  colorRamp2(c(0, 0.01, max(cellmxfn2)), c("#08306B", "white", "red")),
                            name = "normlize enrichment score",
                            show_row_names = F,column_order = cellorder)
    reorder.hmap
    
    ggexport(reorder.hmap, filename = paste0(outpath,"_tag_cluster_heatmap_counts_arti.pdf"),width = 8,height = 8)
    
    saveRDS(list(cluster = split,fit = fit,cellcounts = cellmxfn2,cluster.fig = reorder.hmap),
            file = paste0(outpath,"_cellcounts_cluster.rds"))
    
    #cluster cell composition
    ccls = readRDS(paste0(outpath,"_cellcounts_cluster.rds"))
    
    ccidents = as.data.frame(ccls$cluster)
    ccidents = merge(ccidents,tree$cm,by.x = "row.names",by.y = "tags")
    colnames(ccidents) = c("tags","cluster","CellBC","Freq","Celltype")
    ccidents$Celltype = factor(ccidents$Celltype,levels = cellorder)
    write.csv(ccidents,file = paste0(outpath,"_tag_cluster_BC.csv"),row.names = F,quote = F)
    
    ccidents = ccidents[which(!is.na(ccidents$Celltype)),]
    p3.2 = ggplot(ccidents,aes(x = cluster,fill = Celltype)) +
      geom_bar(stat = "count",position = "fill",width = 0.7) +
      theme_bw()
    p3.2
    ggsave(p3.2,filename = paste0(outpath,"_tag_cluster_celltype_bar.pdf"),width = 8,height = 6)
    
    #cluster show in dimplot
    p3.3 = list()
    for (i in 1:length(clist)) {
      cc = ccidents[which(ccidents$cluster == paste0("c",i)),]
      tarcell = colnames(trans)[which(substr(colnames(E11trans),1,16) %in% cc$CellBC)]
      pt = DimPlot(E11trans,label = TRUE, label.size = 3,cells.highlight = tarcell, 
                   cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
      p3.3[[i]] = pt
      names(p3.3)[i] = paste0("c",i)
      
    }
    p3.3.0 = ggarrange(plotlist =  p3.3,labels = names(p3.3))
    p3.3.0
    
    pdf(paste0(outpath,"_tag_cluster_umap.pdf"),width = 10,height = 8)
    print(p3.3.0)
    invisible(lapply(p3.3, print))
    dev.off()
    
  }
  
  #E15.5 no doublet 作图
  cellallu15 = BCEnrichment(E15tree,domain15,earlyid15,"final_result/xE15_noDoublet/alluv15")
  cellallu15$figure.allu
  saveRDS(cellallu15,file = "final_result/xE15_noDoublet/alluv15.rds")
  ggsave(cellallu11$figure.allu,filename = "final_result/xE15_noDoublet/alluv15_plot.pdf")
  
}

#E11 DIV DA Glu predict
{
  STFtrans = readRDS("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/raw_data/xie/STF_E11_singlet_renamed_20211221.rds")
  e11divcms.pd = as.data.frame(e11divcms.pd)
  Gluclone = e11divcms.pd[which(e11divcms.pd$E11_NbFP_0 > 0 & rowSums(e11divcms.pd[c("DIV7_NbGlu")]) > 0 & rowSums(e11divcms.pd[c("DIV7_DA","DIV7_NbDA")])==0),]
  DAclone = e11divcms.pd[which(e11divcms.pd$E11_NbFP_0 > 0 & rowSums(e11divcms.pd[c("DIV7_DA","DIV7_NbDA")])>0 & rowSums(e11divcms.pd[c("DIV7_NbGlu")]) == 0),]
  rownames(Gluclone)
  rownames(DAclone)
  
  GluBC = E11Stree$cm[which(E11Stree$cm$tags %in% rownames(Gluclone)),"Var1"]
  DABC = E11Stree$cm[which(E11Stree$cm$tags %in% rownames(DAclone)),"Var1"]
  
  tmp = colnames(STFtrans)[which(substr(colnames(STFtrans),1,16)%in% GluBC)]
  Glutrans = subset(STFtrans,cells = tmp)
  
  tmp = colnames(STFtrans)[which(substr(colnames(STFtrans),1,16)%in% DABC)]
  DAtrans = subset(STFtrans,cells = tmp)
  Sdamx = as.data.frame(DAtrans[["RNA"]]@counts)
  Sglumx = as.data.frame(Glutrans[["RNA"]]@counts)
  
  
  #test model
  markers = marker2
  
  
  Sdamxm = Sdamx[rownames(Sdamx) %in% markers,]
  Sglumxm = Sglumx[rownames(Sglumx) %in% markers,]
  Sdamxm = as.data.frame(t(Sdamxm))
  Sglumxm = as.data.frame(t(Sglumxm))
  
  Sdamxmn =  as.data.frame(t(apply(Sdamxm, 1, Normmx)))
  Sglumxmn =  as.data.frame(t(apply(Sglumxm, 1, Normmx)))
  Sdamxmn$type =  0
  Sglumxmn$type = 1
  
  Sdadata = rbind(Sdamxmn,Sglumxmn)
  colnames(Sdadata)[1:(ncol(Sdadata)-1)] = paste0("x",1:(ncol(Sdadata)-1))
  
  
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
  dadata = dadata[which(substr(rownames(dadata),1,16) %in% nbbmid),]
  colnames(dadata)[1:(ncol(dadata)-1)] = paste0("x",1:(ncol(dadata)-1))
  
  traindata = dadata
  
  f <- as.formula(paste0("type ~ ",
                         paste(names(traindata)[1:(ncol(traindata) - 1)],
                               collapse = " + ")))
  library(neuralnet)
  neur <- neuralnet(f, data = traindata, hidden = 10,
                    act.fct = "logistic", linear.output = T,
                    err.fct = "sse", rep = 1)
  saveRDS(neur,file = "final_result/STF_DIV7/neurmodel.rds")
  
  test.hat = compute(neur,Sdadata[-ncol(Sdadata)])$net.result
  
  test.predict = test.hat[,1]
  modelroc = roc(Sdadata$type,test.predict)
  
  test.predict2 = apply(test.hat, 1, function(x){ifelse(x>0,1,0)})
  predict.table = table(Sdadata$type,test.predict2)
  predict.table
  
  model_accuracy <- sum(diag(predict.table))/sum(predict.table)
  model_accuracy
  plot(modelroc)
  
}

