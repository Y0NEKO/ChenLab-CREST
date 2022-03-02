#analysis of main figures
#22/03/02
source("lineage_tree_reconstruction_method.R")

#-------------1. diversity predict-------------
{
  bl = readRDS("final_result/x.blacklist_0.001_2_25.rds")
  samplels = list.files(c("bulk_result/X.bulk_in_vivo/nkx_p1-1/",
                          "bulk_result/X.bulk_in_vivo/nkx_p1-2/",
                          "bulk_result/X.bulk_in_vivo/nkx_p1-3/"),
                        full.names = T)
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = names(v2lsx) = list.files(c("bulk_result/X.bulk_in_vivo/nkx_p1-1/",
                                             "bulk_result/X.bulk_in_vivo/nkx_p1-2/",
                                             "bulk_result/X.bulk_in_vivo/nkx_p1-3/"),
                                           full.names = T)
  
  #nkx venn plot
  nkx1v1 = v1lsx[1:8]
  nkx2v1 = v1lsx[9:16]
  nkx3v1 = v1lsx[17:24]
  
  nkx1v2 = v2lsx[1:8]
  nkx2v2 = v2lsx[9:16]
  nkx3v2 = v2lsx[17:24]
  
  BulkCombine = function(nkx1v1,blv1){
    nkx1v1t = NULL
    for (i in 1:length(nkx1v1)) {
      nkx1v1t = rbind(nkx1v1t,nkx1v1[[i]])
    }
    nkx1v1t = nkx1v1t[which(!nkx1v1t$main %in% blv1$Var1),]
    nkx1v1t = nkx1v1t[which(!nkx1v1t$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                 "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
    nkx1v1t = nkx1v1t[which(nkx1v1t$reads_pro.main >= 0.5 & nkx1v1t$reads_num >= 3),]
    return(nkx1v1t)
  }
  nkx1v1t = BulkCombine(nkx1v1,bl$v1_array)
  nkx2v1t = BulkCombine(nkx2v1,bl$v1_array)
  nkx3v1t = BulkCombine(nkx3v1,bl$v1_array)
  nkx1v2t = BulkCombine(nkx1v2,bl$v2_array)
  nkx2v2t = BulkCombine(nkx2v2,bl$v2_array)
  nkx3v2t = BulkCombine(nkx3v2,bl$v2_array)
  
  
  nkxv1 = list("nkx1v1" = nkx1v1t$main,"nkx2v1" = nkx2v1t$main,"nkx3v1" = nkx3v1t$main)

  nkxv2 = list("nkx1v2" = nkx1v2t$main,"nkx2v2" = nkx2v2t$main,"nkx3v2" = nkx3v2t$main)

  MyVennUndedup = function(nkxv1){
    comn = intersect(nkxv1[[1]],intersect(nkxv1[[2]],nkxv1[[3]]))
    d = list(comn,
             setdiff(intersect(nkxv1[[1]],nkxv1[[2]]),comn),
             setdiff(intersect(nkxv1[[1]],nkxv1[[3]]),comn),
             setdiff(intersect(nkxv1[[2]],nkxv1[[3]]),comn), 
             setdiff(nkxv1[[1]],union(nkxv1[[2]],nkxv1[[3]])),
             setdiff(nkxv1[[2]],union(nkxv1[[1]],nkxv1[[3]])),
             setdiff(nkxv1[[3]],union(nkxv1[[1]],nkxv1[[2]])))
    dn = c(length(nkxv1[[1]][which(nkxv1[[1]] %in% d[[1]])]) + 
             length(nkxv1[[2]][which(nkxv1[[2]] %in% d[[1]])]) +
             length(nkxv1[[3]][which(nkxv1[[3]] %in% d[[1]])]),
           length(nkxv1[[1]][which(nkxv1[[1]] %in% d[[2]])]) + 
             length(nkxv1[[2]][which(nkxv1[[2]] %in% d[[2]])]),
           length(nkxv1[[1]][which(nkxv1[[1]] %in% d[[3]])]) + 
             length(nkxv1[[3]][which(nkxv1[[3]] %in% d[[3]])]),
           length(nkxv1[[2]][which(nkxv1[[2]] %in% d[[4]])]) + 
             length(nkxv1[[3]][which(nkxv1[[3]] %in% d[[4]])]),
           length(nkxv1[[1]][which(nkxv1[[1]] %in% d[[5]])]),
           length(nkxv1[[2]][which(nkxv1[[2]] %in% d[[6]])]),
           length(nkxv1[[3]][which(nkxv1[[3]] %in% d[[7]])])
    )
    ml = data.frame("value" = dn,"percent" = round(dn/sum(dn),3),
                    "nkx1v1" = c(TRUE,  TRUE, TRUE,  FALSE, TRUE, FALSE, FALSE),
                    "nkx2v1" = c(TRUE,  TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                    "nkx3v1" = c(TRUE,  FALSE, TRUE, TRUE, FALSE, FALSE,TRUE))
    colnames(ml) = c("value","percent",names(nkxv1))
    return(ml)
  }
  
  mlv1 = MyVennUndedup(nkxv1)
  mlv2 = MyVennUndedup(nkxv2)
  p2.0.1 = ggvenn(nkxv1,
                  stroke_alpha = 0)
  # detach("package:EnsDb.Mmusculus.v79", unload=TRUE)
  p2.0.2 = ggplot(mlv1) +
    geom_venn(aes(A = nkx1v1, B = nkx2v1, C = nkx3v1,label = paste0(value,"(",percent*100,"%)"))) +
    coord_fixed() +
    theme_void()

  p2.0.1
  p2.0.2
  
  p2.1.1 = ggvenn(nkxv2,
                  stroke_alpha = 0)
  p2.1.2 = ggplot(mlv2) +
    geom_venn(aes(A = nkx1v2, B = nkx2v2, C = nkx3v2,label = paste0(value,"(",percent*100,"%)"))) +
    coord_fixed() +
    theme_void()
  p2.1.1
  p2.1.2
  
  p2.0 = list(p2.0.1,p2.0.2)
  p2.1 = list(p2.1.1,p2.1.2)
  
  pdf("final_result/xdraft/main_figure_2_25/Bulk_Bulk_nkxv1_venn.pdf")
  invisible(lapply(p2.0, print))
  dev.off()
  
  pdf("final_result/xdraft/main_figure_2_25/Bulk_nkxv2_venn.pdf")
  invisible(lapply(p2.1, print))
  dev.off()
  
  
  #(3)indel for plotting pattern
  ExtractScar = function(v1lsx,blst){
    arraydata = NULL
    for (i in 1:length(v1lsx)) {
      arrayline = v1lsx[[i]]
      arrayline = arrayline[which(arrayline$reads_pro.main >= 0.5 & arrayline$reads_num >= 3),]
      arrayline = arrayline[which(!arrayline$main %in% blst$Var1),]
      arrayline = arrayline[which(!arrayline$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      arraydata = c(arraydata,arrayline$main)
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
  ExtractScar2 = function(v1lsx,blst){
    arraydata = NULL
    for (i in 1:length(v1lsx)) {
      arrayline = v1lsx[[i]]
      arrayline = arrayline[which(arrayline$reads_pro.main >= 0.5 & arrayline$reads_num >= 3),]
      arrayline = arrayline[which(!arrayline$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                         "v1.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE-v2.NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      
      arraydata = c(arraydata,arrayline$main)
    }
    
    
    statd = NULL
    for (n in seq(1,length(arraydata),100)) {
      #    s = read[sample(1:nrow(read),n),]
      #    pass = s[which(s$quality == "pass"),]
      arraydata[which(arraydata %in% blst$Var1)] = "pass"
      s = sample(arraydata,n)
      d = length(unique(s[which(s!="pass")]))
      line = data.frame("umi" = n,"diversity" = d)
      statd = rbind(statd,line)
      
      
    }
    return(statd)
  }
  
  
  v1fdiv = ExtractScar(v1lsx,NULL)
  v2fdiv = ExtractScar(v2lsx,NULL)
  
  v3lsxrand = list()
  v3lsxmin = list()
  for (i in 1:length(v1lsx)) {
    line1 = v1lsx[[i]]
    line2 = v2lsx[[i]]
    line1 = line1[which(line1$reads_pro.main >= 0.5 & line1$reads_num >= 3),]
    line1 = line1[which(!line1$main %in% "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"),]
    line2 = line2[which(line2$reads_pro.main >= 0.5 & line2$reads_num >= 3),]
    line2 = line2[which(!line2$main %in% "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE"),]
    line1$main = paste0("v1.",line1$main)
    line2$main = paste0("v2.",line2$main)
    if(nrow(line1)>nrow(line2)){linel = line1;lines = line2}else{linel = line2;lines = line1}
    
    indel1 = table(lines$main);indel1 = indel1[order(indel1)]
    indel2 = table(linel$main);indel2 = indel2[order(indel2)]
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
      if((cmbline$v1 %in% paste0("v1.",bl$v1_array$Var1)) & (cmbline$v2 %in% paste0("v2.",bl$v2_array$Var1))){
        next
      }
      if(cmbline$v1.count < cmbline$v2.count){
        linel[which(linel$main == cmbline$v2),"main"][1:cmbline$v1.count] = paste0(cmbline$v1,"-",cmbline$v2)
      }else{
        linel[which(linel$main == cmbline$v2),"main"] = paste0(cmbline$v1,"-",cmbline$v2)
      }
      
    }
    v3lsxmin[[i]] = linel
    
    if(nrow(line1)>nrow(line2)){
      linel = line1;
      lp = paste0("v1.",bl$v1_array$Var1)
      lines = line2
      sp = paste0("v2.",bl$v2_array$Var1)
    }else{
      linel = line2;
      lp = paste0("v2.",bl$v2_array$Var1)
      lines = line1
      sp = paste0("v1.",bl$v1_array$Var1)}
    
    cmid = sample(1:nrow(linel),nrow(lines))
    linel$short = "NONE"
    linel$short[cmid] = lines$main
    linel$long = linel$main
    linel$main[cmid] = paste0(lines$main,"-",linel$main[cmid])
    linel = linel[which(!(linel$short %in% sp) | !(linel$long %in% lp)),]
    v3lsxrand[[i]] = linel
    names(v3lsxrand)[i] = names(v3lsxmin)[i]= names(v1lsx)[i]
    
  }
  
  ExtractScar3 = function(v1lsx){
    arraydata = NULL
    for (i in 1:length(v1lsx)) {
      arrayline = v1lsx[[i]]
      arraydata = c(arraydata,arrayline$main)
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
  v3fdivmin = ExtractScar3(v3lsxmin)
  v3fdivrand = ExtractScar3(v3lsxrand)
  
  v1fdiv[nrow(v1fdiv),]
  v2fdiv[nrow(v2fdiv),]
  v3fdivrand[nrow(v3fdivrand),]
  v3fdivmin[nrow(v3fdivmin),]
  
  
  write.table(v1fdiv, file = "final_result/xdraft/main_figure_2_25/Bulk_diversity_predict_v1_nobl.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v2fdiv, file = "final_result/xdraft/main_figure_2_25/Bulk_diversity_predict_v2_nobl.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v3fdivmin, file = "final_result/xdraft/main_figure_2_25/Bulk_diversity_predict_v3_min.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  write.table(v3fdivrand, file = "final_result/xdraft/main_figure_2_25/Bulk_diversity_predict_v3_rand.txt", 
              quote = F, row.names = F, col.names = F,sep = "\t")
  
  #The optimal parameters were obtained by fitting the model by matlab
  v1p1 =   2.456e+04;
  v1p2 =    5.39e+07;
  v1q1 =    6.36e+04;
  v2p1 =   4.777e+04;
  v2p2 =   3.563e+08;
  v2q1 =   1.591e+05;
  v3mp1 =   6.661e+04;
  v3mp2 =   3.157e+08;
  v3mq1 =   1.724e+05;
  v3rp1 =   9.224e+05;
  v3rp2 =   1.824e+10;
  v3rq1 =    2.52e+06;
  
  v1fun = function(x){(v1p1*x + v1p2)/(x + v1q1)}
  v2fun = function(x){(v2p1*x + v2p2)/(x + v2q1)}
  v3funrand = function(x){(v3rp1*x + v3rp2)/(x + v3rq1)}
  v3funmin = function(x){(v3mp1*x + v3mp2)/(x + v3mq1)}
  
  v1fdiv$type = "v1"
  v2fdiv$type = "v2"
  vp1 = rbind(v1fdiv,v2fdiv)
  options(scipen = 3)
  
  library(ggplot2)
  library(viridis)
  p1 = ggplot(vp1,aes(x = umi,y = diversity,color = type)) + xlim(0,2000000) + ylim(0,60000) +
    geom_point(size = 0.5) +
    geom_function(aes(colour = "v1"),fun = v1fun,linetype = "dashed") +
    geom_function(aes(colour = "v2"),fun = v2fun,linetype = "dashed") +
    annotate("text", x=800000, y=50000, label= "max(v1) = 24560,max(v2) = 47770") + 
    theme_ipsum(base_size = 18, base_family = "sans") + labs(x = NULL, y = NULL)
  
  p1
  
  v3fdivrand$type = "random"
  v3fdivmin$type = "min"
  vp2 = rbind(v3fdivrand,v3fdivmin)
  p2 = ggplot(vp2,aes(x = umi,y = diversity,color = type)) + xlim(0,100000000) + ylim(0,1000000) +
    geom_point(size = 0.5) +
    geom_function(aes(colour = "random"),fun = v3funrand,linetype = "dashed") +
    geom_function(aes(colour = "min"),fun = v3funmin,linetype = "dashed") +
    annotate("text", x=50000000, y=1000000, label= "max(random) = 922400,max(min) = 66610") + 
    theme_ipsum(base_size = 18, base_family = "sans") + labs(x = NULL, y = NULL)
  
  p2
  
  ggsave(p1, filename = "final_result/xdraft/main_figure_2_25/Bulk_Diversity_predict_nobl.pdf",width = 8,height = 6)
  ggsave(p2, filename = "final_result/xdraft/main_figure_2_25/Bulk_Diversity_predict_v1v2_simulation.pdf",width = 10,height = 6)
  
  
  #cell unique rate predict
  v1lsxt = BulkCombine(v1lsx,NULL)
  v2lsxt = BulkCombine(v2lsx,NULL)
  
  CellUniqueStat = function(v1lsxt,bl1){
    unistv1 = NULL
    for (i in seq(1,nrow(v1lsxt),100)) {
      print(i)
      line = v1lsxt[sample(1:nrow(v1lsxt),i),]
      alrate = table(line$main)
      alratef = alrate[which(!names(alrate) %in% bl1)]
      unicln = length(alrate[which(alrate == 1)])
      uniclnf = length(alratef[which(alratef == 1)])
      uniclr =  unicln / i
      uniclrf =  uniclnf / sum(alratef)
      # rareal = length(alrate[which(alrate/i < 0.01)])
      unistv1 = rbind(unistv1,data.frame("cellcounts" = i,"unicellpos" = uniclr,
                                         "unicellwhitepos" = uniclrf,
                                         "unicellnum" = unicln,
                                         "unicellwhitenum" = unicln))
    }
    scaleFactor <-  max(unistv1$unicellnum) / max(unistv1$unicellpos) 
    p3.1 = ggplot(unistv1, aes(x=cellcounts)) +
      geom_line(aes(y = unicellnum),  col = "#012443") +
      geom_line(aes(y = unicellpos * scaleFactor),  col="#B61919") +
      geom_line(aes(y = unicellwhitepos * scaleFactor),  col="#297F87") +
      scale_y_continuous(name="counts of unique cells", sec.axis=sec_axis(~./scaleFactor, name="proportion of unique cells")) +
      theme_pubr() +
      theme(
        axis.title.y.left=element_text(color="#012443"),
        axis.text.y.left=element_text(color="#012443"),
        axis.ticks.y.left=element_line(color="#012443"),
        axis.line.y.left=element_line(color="#012443"),
        axis.title.y.right=element_text(color="#B61919"),
        axis.text.y.right=element_text(color="#B61919"),
        axis.ticks.y.right=element_line(color="#B61919"),
        axis.line.y.right=element_line(color="#B61919")
      )
    v1unist = list("unist" = unistv1,"fg" = p3.1)
    return(v1unist)
  }
  
  v1unist = CellUniqueStat(v1lsxt,bl$v1_array$Var1)
  v2unist = CellUniqueStat(v2lsxt,bl$v2_array$Var1)
  
  library(ggforce)
  scaleFactorv1 <-  max(v1unist$unist$unicellnum) / max(v1unist$unist$unicellpos) 
  scaleFactorv2 <-  max(v2unist$unist$unicellnum) / max(v2unist$unist$unicellpos) 
  v1unist$fg = ggplot(v1unist$unist, aes(x=cellcounts)) +
    geom_line(aes(y = unicellnum),  col = "#012443") +
    geom_line(aes(y = unicellpos * scaleFactorv1),  col="#B61919") +
    geom_line(aes(y = unicellwhitepos * scaleFactorv1),  col="#297F87") +
    facet_zoom(xlim = c(0,v1unist$unist[which(v1unist$unist$unicellpos < 0.5),"cellcounts"][1]))+
    scale_y_continuous(name="counts of unique cells", sec.axis=sec_axis(~./scaleFactorv1, name="proportion of unique cells")) +
    theme_pubr() +
    theme(
      axis.title.y.left=element_text(color="#012443"),
      axis.text.y.left=element_text(color="#012443"),
      axis.ticks.y.left=element_line(color="#012443"),
      axis.line.y.left=element_line(color="#012443"),
      axis.title.y.right=element_text(color="#B61919"),
      axis.text.y.right=element_text(color="#B61919"),
      axis.ticks.y.right=element_line(color="#B61919"),
      axis.line.y.right=element_line(color="#B61919")
    )
  
  v2unist$fg = ggplot(v2unist$unist, aes(x=cellcounts)) +
    geom_line(aes(y = unicellnum),  col = "#012443") +
    geom_line(aes(y = unicellpos * scaleFactorv2),  col="#B61919") +
    geom_line(aes(y = unicellwhitepos * scaleFactorv2),  col="#297F87") +
    facet_zoom(xlim = c(0,v2unist$unist[which(v2unist$unist$unicellpos < 0.5),"cellcounts"][1]))+
    scale_y_continuous(name="counts of unique cells", sec.axis=sec_axis(~./scaleFactorv2, name="proportion of unique cells")) +
    theme_pubr() +
    theme(
      axis.title.y.left=element_text(color="#012443"),
      axis.text.y.left=element_text(color="#012443"),
      axis.ticks.y.left=element_line(color="#012443"),
      axis.line.y.left=element_line(color="#012443"),
      axis.title.y.right=element_text(color="#B61919"),
      axis.text.y.right=element_text(color="#B61919"),
      axis.ticks.y.right=element_line(color="#B61919"),
      axis.line.y.right=element_line(color="#B61919")
    )
  v1unist$fg
  v2unist$fg
  
  ggsave(v1unist$fg, filename = "final_result/xdraft/main_figure_2_25/Bulk_unicell_v1.pdf",width = 6,height = 6)
  ggsave(v2unist$fg, filename = "final_result/xdraft/main_figure_2_25/Bulk_unicell_v2.pdf",width = 6,height = 6)
  saveRDS(v1unist,file = "final_result/xdraft/main_figure_2_25/Bulk_unicell_v1.rds")
  saveRDS(v2unist,file = "final_result/xdraft/main_figure_2_25/Bulk_unicell_v2.rds")
  
  
  
}

#-------------2. bulk array stat--------------
{
  samplels = list.files("bulk_result/X.bulk_in_vivo/",full.names = T)[c(1:8,9,10,12:13)]
  v1lsx = lapply(paste0(samplels,"/v1/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  v2lsx = lapply(paste0(samplels,"/v2/final_scarform.csv"), function(x) read.csv(x,header = T)) 
  names(v1lsx) = names(v2lsx) = list.files("bulk_result/X.bulk_in_vivo/")[c(1:8,9,10,12:13)]
  
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
      
      scardf = as.data.frame(table(v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]$umim))
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
  
  
  names(v1lsx) = names(v2lsx) = c("ctrl","E10-1","E10-2","E12","E15-1a","E15-1b",
                                  "E15-2a","E15-2b","E15-3a","E15-3c","E9","P0")
  v1lsx.r = v1lsx[c(1,11,2:10,12)]
  v2lsx.r = v2lsx[c(1,11,2:10,12)]
  colrank = names(v1lsx.r) = names(v2lsx.r) = c("ctrl","E9","E10-1","E10-2","E12","E15-1a","E15-1b",
                                                "E15-2a","E15-2b","E15-3a","E15-3c","P0")
  library(stringr)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  scardf1 = ScarPlot(v1lsx.r,20,colrank)
  scardf2 = ScarPlot(v2lsx.r,20,colrank)
  scardf1$figure
  scardf2$figure
  
  ggsave(scardf1$figure,filename = "final_result/xBulk/bulk_top_array_v1_2_17_qc.pdf",width = 10,height = 6)
  ggsave(scardf2$figure,filename = "final_result/xBulk/bulk_top_array_v2_2_17_qc.pdf",width = 10,height = 6)
  
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
    arraydft$E15_3 = (arraydft$`E15-3a`+ arraydft$`E15-3c`)/2
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
  ggsave(p1,filename = "final_result/xBulk/barcode_occrurency_in_sample_2_17_qc.pdf",width = 7,height = 5)
  write.table(ayocv,file = "final_result/xBulk/barcode_occrurency_in_sample_2_17_qc.txt",row.names = F,quote = F,
              sep = "\t")
  
  #bulk indel occurrence histogram
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
  
  p1 = ggplot(idocv,aes(x = Var1, y = norm, fill = group)) + 
    geom_bar(stat = "identity",position = "dodge",width = 0.7) +
    scale_fill_manual(values=c('lightpink1','lightblue2')) +
    theme_bw() + xlab("Frequency of occurrence") + ylab("Proportion") + 
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))
  print(p1)
  
  ggsave(p1,filename = "final_result/xBulk/indel_occrurency_in_sample_2_17_qc.pdf",width = 7,height = 5)
  write.table(idocv,file = "final_result/xBulk/indel_occrurency_in_sample_2_17_qc.txt",row.names = F,
              quote = F,sep = "\t")
  
  
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
  ggsave(v1edp,filename = "final_result/xBulk/edit_size_v1_2_17_qc.pdf",width = 10,height = 6)
  ggsave(v2edp,filename = "final_result/xBulk/edit_size_v2_2_17_qc.pdf",width = 10,height = 6)
  
  
  #edit array proportions
  names(v1lsx)
  EditArrayStat = function(v1lsx){
    propft = NULL
    for (i in 2:length(v1lsx)) {
      arraydf = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
      arraydff = arraydf[which(!arraydf$main %in% c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE",
                                                    "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")),]
      propl = nrow(arraydff) / nrow(arraydf)
      propft = rbind(propft,data.frame("sample" = names(v1lsx)[i],"edit_prop" = propl))
    }
    
    #intact target sites
    
    tarsite = NULL
    for (i in 2:length(v1lsx)) {
      # arraydf = v1lsx[[i]]
      arraydf = v1lsx[[i]][which(v1lsx[[i]]$reads_pro.main >= 0.5 & v1lsx[[i]]$reads_num >= 3),]
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
  write.table(v1propft,file = "final_result/xBulk/eidt_stat_v1_2_17_qc.txt",row.names = F,
              quote = F,sep = "\t")
  write.table(v2propft,file = "final_result/xBulk/eidt_stat_v2_2_17_qc.txt",row.names = F,
              quote = F,sep = "\t")
  
  
  
  
}

#-------------3. cell clusters linkage relationship--------------
{
  E11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  E15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  E11tree = readRDS("final_result/xdraft/tree_blacklistchange_2_24/E11_bl0.001_blay_tree.rds")
  E15tree = readRDS("final_result/xdraft/tree_blacklistchange_2_24/E15_bl0.001_blay_tree.rds")
  E11stree = readRDS("final_result/xdraft/tree_blacklistchange_2_24/E11_STF_bl0.001_blay_tree.rds")
  Divtree = readRDS("final_result/xdraft/tree_blacklistchange_2_24/Div_STF_bl0.001_blay_tree.rds")
  
  domain11 = read.csv("sc_result/X.E11/E11pseudoAxis20211126V4.csv")
  domain15 = read.csv("sc_result/x.TAQ_E15/E15pseudoAxis20211126v5.csv")
  colnames(domain11) = colnames(domain15) = c("celltype","domainid","domain")
  domain11[12,1] = "NbFP_0"
  earlyid11 = c("NProgBL","NProgBM","Rgl1")
  earlyid15 = c("NProgAL","Rgl1","Rgl2AL","Rgl2BL","GabaProgBI","GabaProgBL","NbBM1","NbFP","Rgl3")
  
  BCEnrichment2gen = function(tree,domain,earlyid){
    cm = tree$cm
    cmtn = merge(cm,domain,by = "celltype")
    lateid = setdiff(unique(cmtn$celltype),earlyid)
    cmtn$phase = "late"
    cmtn[which(cmtn$celltype %in% earlyid),"phase"] = "early"
    cmtn = cmtn[order(cmtn$phase,cmtn$domain),]
    {
      cmtnp = cmtn %>% group_by(tags,phase,celltype) %>% 
        summarise(counts = sum(n()),domain = unique(domain))
      #
      cellmx.e = dcast(cmtnp[which(cmtnp$phase == "early"),c("tags","celltype","counts")],
                       tags~celltype,value.var = "counts",fun.aggregate = sum)
      cellmx.l = dcast(cmtnp[which(cmtnp$phase == "late"),c("tags","celltype","counts")],
                       tags~celltype,value.var = "counts",fun.aggregate = sum)
      cellft = merge(cellmx.e,cellmx.l,by = "tags")
      
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
    
    #simulation phasest
    # #method1
    # phasestm = data.frame("phase1" = rep(earlyid,each = length(lateid)),
    #                      "phase2" = rep(lateid,length(earlyid)), "Freq.from" = 0, "Freq.to" = 0)
    # Me = phasest %>% group_by(phase1) %>% summarise(counts = sum(Freq.from))
    # Ne = sum(Me$counts)
    # Ml = phasest %>% group_by(phase2) %>% summarise(counts = sum(Freq.to))
    # Nl = sum(Ml$counts)
    # for (i in 1:nrow(phasestm)) {
    #   lme = Me[which(Me$phase1 == phasestm[i,]$phase1),"counts"] 
    #   lml = Ml[which(Ml$phase2 == phasestm[i,]$phase2),"counts"]
    #   phasestm[i,]$Freq.from = lme*lml/Nl
    #   phasestm[i,]$Freq.to = lml*lme/Ne
    # }
    
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
  BCEnrichment2 = function(allu){
    {
      phasestf = NULL
      phasest = allu$phasest
      phasestmle = allu$phasestmle
      phasestmll = allu$phasestmll
      cmtn = allu$cmtn
      for (i in 1:nrow(phasest)) {
        line = phasest[i,]
        if(line$Freq.from > quantile(phasestmle[[i]],probs = 0.05) | 
           line$Freq.to > quantile(phasestmll[[i]],probs = 0.05)){
          phasestf = rbind(phasestf,line)
        }
        
      }
      
      domain_df = unique(cmtn[c("celltype","domainid","domain")])
      domain_df = domain_df[order(domain_df$domain),]
      # phasestf = phasest[which(phasest$Freq.from > 0.001 | phasest$Freq.to > 0.001),]
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
                       fill = domain,color = celltype,label = celltype)) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_flow(aes(color = celltype)) +
        geom_flow(aes()) +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) + theme_void()
      pf1
      
      
    }
    allu$pf1 = pf1
    return(allu)
  }
  cellallu15 = BCEnrichment2gen(E15tree,domain15,earlyid15)
  cellallu15 = BCEnrichment2(cellallu15)
  cellallu11 = BCEnrichment2gen(E11tree,domain11,earlyid11)
  cellallu11 = BCEnrichment2(cellallu11)
  
  cellallu15$pf1
  cellallu11$pf1
  
  #E11-DIV7
  {
    cellmx.e = E11stree$count
    cellmx.l = Divtree$count
    colnames(cellmx.e)[-1] = paste0(colnames(cellmx.e)[-1],".E11")
    colnames(cellmx.l)[-1] = paste0(colnames(cellmx.l)[-1],".DIV")
    cellft = merge(cellmx.e, cellmx.l, by = "tags")
    earlyid = colnames(cellmx.e)[-1]
    lateid = colnames(cellmx.l)[-1]
    cellft = cellft[which(rowSums(cellft[-1]>0) < ncol(cellft[-1])/3),]
    cellft = cellft[,which(!colnames(cellft) %in% c("OMTN.E11","Unknown.DIV", "GabaUnknown.DIV", 
                                                    "Rgl2.DIV", "Rgl3.DIV"))]
    earlyid = earlyid[which(!earlyid %in% c("OMTN.E11"))]
    lateid = lateid[which(!lateid %in% c("Unknown.DIV", "GabaUnknown.DIV", 
                                         "Rgl2.DIV", "Rgl3.DIV"))]
    
    domain_df = data.frame("celltype" = c("NbGluAL.E11","NProgBL.E11","GabaProg0.E11","GabaProgBL.E11",
                                          "NbFP_0.E11","Peri.E11","Rgl1.E11","OMTN.E11",
                                          "GabaProgBI.E11","NbBM1.E11","NProgBM.E11","NbBM0.E11",
                                          "GabaBL1.DIV","GabaBL4.DIV","GabaBL2.DIV","GluAL1.DIV",
                                          "GabaUnknown.DIV","Unknown.DIV","GabaBI.DIV",
                                          "Rgl2.DIV","Rgl3.DIV","DA.DIV","NbGlu.DIV","NbDA.DIV",
                                          "NbFP.DIV","GluBM1.DIV","GluBM2.DIV","GluFP.DIV"),
                           "domain" = rep(c("c1","c2","c3","c1","c2","c3"),c(4,4,4,7,6,3)))
    domain_df = domain_df[order(domain_df$domain),]
    celllevel = c('NProgBL','NProgBM','Rgl1','NbGluAL','GabaProgBL','GabaProg0','GabaProgBI','NbBM0','NbBM1','NbFP_0',
                  'GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2','GluFP','NbGlu','NbFP','NbDA','DA')
    celllevel1 = c('NProgBL','NProgBM','Rgl1',
                   'GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2','GluFP','NbGlu','NbFP','NbDA','DA')
    celllevel2 = c( 'NbGluAL','GabaProgBL','GabaProg0','GabaProgBI','NbBM0','NbBM1','NbFP_0',
                    'GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2','GluFP','NbGlu','NbFP','NbDA','DA')
    celllevel = paste0(celllevel, rep(c(".E11",".DIV"),c(10,12)))
    celllevel1 = paste0(celllevel1, rep(c(".E11",".DIV"),c(3,12)))
    celllevel2 = paste0(celllevel2, rep(c(".E11",".DIV"),c(7,12)))
    # cellft1 = cellft[,which(colnames(cellft) %in% celllevel1)]
    # cellft2 = cellft[,which(colnames(cellft) %in% celllevel2)]
    
    DIVAllu = function(cellft,celllevel){
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
        if(line$Freq.from > quantile(phasestmle[[i]],probs = 0.05) | 
           line$Freq.to > quantile(phasestmll[[i]],probs = 0.05)){
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
                         fill = domain,color = celltype,label = celltype)) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_flow(aes(color = celltype)) +
          geom_flow(aes()) +
          geom_stratum(alpha = .5) +
          geom_text(stat = "stratum", size = 3) + theme_void()
        pf1
        return(pf1)
        
      }
      phasest1 = phasest[which(phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11')),]
      phasest2 = phasest[which(!phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11')),]
      phasestf1 = phasestf[which(phasestf$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11')),]
      phasestf2 = phasestf[which(!phasestf$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11')),]
      pf1 = PhasePlot(phasestf1,domain_df,celllevel)
      pf2 = PhasePlot(phasestf2,domain_df,celllevel)
      
      which(phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11'))
      
      cellallud71 = list("cellft" = cellft,"phasest" = phasest1,"phasestmle" = phasestmle[which(phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11'))],
                         "phasestmll" = phasestmll[which(phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11'))],"pf1" = pf1)
      cellallud72 = list("cellft" = cellft,"phasest" = phasest2,"phasestmle" = phasestmle[which(!phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11'))],
                         "phasestmll" = phasestmll[which(!phasest$phase1 %in% c('NProgBL.E11','NProgBM.E11','Rgl1.E11'))],"pf1" = pf2)
      cellallud7 = list("prog" = cellallud71,"other" = cellallud72)
      
      return(cellallud7)
      
    }
    
    cellallud7 = DIVAllu(cellft,celllevel)
    
  }
  cellallud7$prog$pf1
  cellallud7$other$pf1
  ggsave(cellallu11$pf1,filename = "final_result/xdraft/main_figure_2_25/E11_alluv_fil0.05.pdf",width = 6,height = 6)
  ggsave(cellallu15$pf1,filename = "final_result/xdraft/main_figure_2_25/E15_alluv_fil0.05.pdf",width = 6,height = 6)
  ggsave(cellallud7$prog$pf1,filename = "final_result/xdraft/main_figure_2_25/DIV7_alluv_fil0.05_prog.pdf",width = 6,height = 6)
  ggsave(cellallud7$other$pf1,filename = "final_result/xdraft/main_figure_2_25/DIV7_alluv_fil0.05_other.pdf",width = 6,height = 6)
  
  #Plot Pvalue
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
    cellallu11$pv = p1.2
    return(cellallu11)
  }
  cellallu11 = CalAlluPvalue(cellallu11)
  cellallu15 = CalAlluPvalue(cellallu15)
  cellallud7$prog = CalAlluPvalue(cellallud7$prog)
  cellallud7$other = CalAlluPvalue(cellallud7$other)
  cellallu11$pv
  cellallu15$pv
  cellallud7$prog$pv
  cellallud7$other$pv
  
  ggsave(cellallu11$pv,filename = "final_result/xdraft/main_figure_2_25/E11_alluv_pvalue.pdf",width = 6,height = 6)
  ggsave(cellallu15$pv,filename = "final_result/xdraft/main_figure_2_25/E15_alluv_pvalue.pdf",width = 6,height = 6)
  ggsave(cellallud7$prog$pv,filename = "final_result/xdraft/main_figure_2_25/DIV7_prog_alluv_pvalue.pdf",width = 6,height = 6)
  ggsave(cellallud7$other$pv,filename = "final_result/xdraft/main_figure_2_25/DIV7_other_alluv_pvalue.pdf",width = 6,height = 6)
  
  saveRDS(cellallu15,file = "final_result/xdraft/main_figure_2_25/E15_alluv.rds")
  saveRDS(cellallu11,file = "final_result/xdraft/main_figure_2_25/E11_alluv.rds")
  saveRDS(cellallud7,file = "final_result/xdraft/main_figure_2_25/DIV7_alluv.rds")
  
  
}

#-------------4. domain dispersion-------------
{
  PseudoDisper = function(tree,domain){
    library(scatterpie)
    cm = tree$cm
    cmtn = merge(cm,domain,by = "celltype")
    cmtn = cmtn[order(cmtn$domainid),]
    cmtn$celltype = factor(as.character(cmtn$celltype),levels = unique(cmtn$celltype))
    unique(domain$domain)
    
    cmtnf = cmtn[which(cmtn$domain != "unknown"),]
    cmtnf$celltype = factor(cmtnf$celltype,levels = unique(cmtnf$celltype))
    cmtnf = cmtnf[order(cmtnf$domainid),]
    cmtnf$domainid = factor(cmtnf$domainid)
    
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
  pvs11 = PseudoDisper(E11tree,domain11)
  pvs15 = PseudoDisper(E15tree,domain15)
  pvs11$figure_pvs
  pvs15$figure_pvs
  pvs11$figure_pvs_box
  ggsave(pvs11$figure_pvs,filename = "final_result/xdraft/main_figure_2_25/E11_domain_disper_cmb.pdf",width = 8,height = 6)
  ggsave(pvs15$figure_pvs,filename = "final_result/xdraft/main_figure_2_25/E15_domain_disper_cmb.pdf",width = 8,height = 6)
  
  #domain heatmap
  ardm = pvs11$pvs_domain
  ardm = ardm[-c(1:6,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xdraft/main_figure_2_25/E11_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
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
  ardm = ardm[-c(1:6,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xdraft/main_figure_2_25/E15_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
               row_order = unique(domain15$domain)[!unique(domain15$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  
  corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
                "NbBM1","Rgl1","NbFP_0")
  domain11[12,1] = "NbFP_0"

  pvss11 = PseudoDisper(E11stree,domain11)
  pvss11$figure_pvs
  pvss11$figure_pvs_box
  ggsave(pvss11$figure_pvs,filename = "final_result/xdraft/main_figure_2_25/E11_STF_domain_disper_cmb.pdf",width = 8,height = 6)
  ggsave(pvss11$figure_pvs_box,filename = "final_result/xdraft/main_figure_2_25/E11_STF_domain_disper_box.pdf",width = 5,height = 3)
  
  #domain heatmap
  ardm = pvss11$pvs_domain
  ardm = ardm[-c(1:6,ncol(ardm))]
  ardm = ardm[which(rowSums(ardm>0)>1 & rowSums(ardm>0.01)<ncol(ardm)),];ardm = as.matrix(ardm)
  
  mx = cor(ardm,method = "pearson")
  pdf("final_result/xdraft/main_figure_2_25/E11_STF_domain_pearson_heatmap.pdf",width = 6,height = 5)
  
  p1 = Heatmap(mx,heatmap_legend_param = list(title = "dist"),
               col = c("#08519C",colorRampPalette(brewer.pal(6, "Reds"))(4)),
               column_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               row_order = unique(domain11$domain)[!unique(domain11$domain)%in%"unknown"],
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", mx[i, j]), x, y, gp = gpar(fontsize = 10))
               })
  print(p1)
  dev.off()
  
  
  
  
}

#-------------5. DIV E11 relationship-------------
{
  #edit heatmap
  corderdiv = c('GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2','GluFP','NbGlu','NbFP','NbDA','DA',
                "GabaUnknown","Rgl2","Rgl3","Unknown")
  corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
                "NbBM1","Rgl1","NbFP_0")
  
  count1 = Divtree$count
  count2 = E11stree$count
  
  count1 = count1[,c("tags",corderdiv)]
  count2 = count2[,c("tags",corder11n)]
  
  colnames(count1)[-1] = paste0("DIV7_",colnames(count1)[-1])
  colnames(count2)[-1] = paste0("E11_",colnames(count2)[-1])
  e11divcms = merge(count2,count1,by = "tags")
  
  e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$tags;e11divcms.pd = e11divcms.pd[-1]
  e11divcms.pd = as.matrix(e11divcms.pd)
  e11divcms.pdt = e11divcms.pd
  e11divcms.pdt[which(e11divcms.pdt>0)] = 1
  
  split = factor(rep(c("E11","DIV7"),c(ncol(E11stree$count[-1]),ncol(Divtree$count[-1]))),levels = c("E11","DIV7"))
  names(split) = colnames(e11divcms.pd)
  
  e11divcms.pd = e11divcms.pd[do.call(order, as.data.frame(-e11divcms.pdt)), ]
  
  tmp = e11divcms.pdt[do.call(order, as.data.frame(-e11divcms.pdt)), ]
  
  pdf("final_result/xdraft/main_figure_2_25/STF_E11-DIV7_array_cmp_heatmap_log2counts.pdf",width = 6,height = 8)
  Heatmap(log2(e11divcms.pd+1),show_row_names = F,
          # col = c("#23174C","#376BB4","white"),
          col = colorRampPalette(brewer.pal(9, "Reds"))(30),
          column_gap = unit(5, "mm"),name = "log2(counts)",
          cluster_row_slices = FALSE,row_order = rownames(tmp),
          column_split = split,column_order = colnames(tmp))
  dev.off()
  
  #cor heatmap
  e11divcmf.mx = cor(as.matrix(e11divcms.pd),method = "pearson")
  e11divcmf.mx = e11divcmf.mx[1:(ncol(count2)-1),(ncol(count2)):ncol(e11divcmf.mx)]
  pdf("final_result/xdraft/main_figure_2_25/STF_E11-DIV7_cor_spearman_heatmap.pdf",width = 6,height = 5)
  Heatmap(e11divcmf.mx,
          clustering_method_rows = "complete")
  dev.off()
  
}

#-------------6. clone cluster-------------
{
  library(ggfortify)
  library(ggplot2)
  library(dplyr)
  cellorder = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
                "NbBM1","Rgl1","NbFP_0","Peri")
  tree = E11tree
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
    
    
    #first kmer cluster
    rownames(cellft) = cellft$tags;cellft = cellft[-1];cellftn = scale(cellft)
    tmp = table(cellftl[which(cellftl$value>0),]$tags)
    tmp = tmp[which(tmp>1)]
    cellftf = cellft[which(rownames(cellft)%in% names(tmp)),]
    cellmxf = cellmx[-1];rownames(cellmxf) = cellmx$tags;cellmxf = cellmxf[which(rownames(cellmxf)%in% names(tmp)),]
    
    #1.clone cluster heatmap------------------
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
    
    #2.cluster cell composition----------------
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
    
    #3.cluster show in dimplot---------------
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
  tree$cm = tree$cm$celltype[which(tree$cm$celltype %in% c("GabaUnknown","Unknown","Rgl2","Rgl3")),]
  
  #clone cluster suplement-------------
  e15_clu = readRDS("final_result/xdraft/main_figure_2_25/E15_cellcounts_cluster_C6_20220225.rds")
  cellordern = c('GluAL4','GluAL3','GluAL1','Rgl2AL','NProgAL','GabaAL','GluAL2',
                 'Rgl3','GabaProgBL','GabaBL2','GabaBL3','GabaBL4','GabaBL1','Rgl2BL',
                 'GabaBI','GabaProgBI','NbBM1','GluBM1',
                 'GluBM2','Rgl1','NbGlu','GluFP','NbFP','DA','NbDA')
  e15cc = E15tree$count
  rownames(e15cc) = e15cc$tags;e15cc = e15cc[-1]
  e15cc = e15cc[which(rowSums(e15cc) > 1),]
  e15cc = as.data.frame(t(t(as.matrix(e15cc))/colSums(e15cc)))
  e15cc = e15cc[,cellordern]
  e15cc = e15cc[do.call(order, -e15cc),]
  adid15 = rownames(e15cc)[which(!rownames(e15cc) %in% names(e15_clu$cluster))]
  tmp = e15cc[which(rownames(e15cc) %in% adid),]
  
  e15add.hmap <- Heatmap(as.matrix(tmp), name = "filtered clones' cell portion",
                         row_order = rownames(tmp),
                         # col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                         show_row_names = F,column_order = cellordern)
  e15add.hmap
  e15addls = list("mx" = tmp,"roworder" = rownames(tmp),"colorder" = cellordern,"fg" = e15add.hmap)
  saveRDS(e15addls,"final_result/xdraft/main_figure_2_25/E15_tag_cluster_heatmap_counts_arti_supplement.rds")
  
  ggexport(e15add.hmap, filename = "final_result/xdraft/main_figure_2_25/E15_tag_cluster_heatmap_counts_arti_supplement.pdf",width = 8,height = 8)
  corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","NbBM0",
                "NbBM1","OMTN","Rgl1","NbFP_0",'Peri')
  e11_clu = readRDS("final_result/xdraft/main_figure_2_25/E11_cellcounts_cluster.rds")
  e11_clu$cluster.fig
  E11tree_peri = readRDS("final_result/xdraft/tree_blacklistchange_2_24/E11_bl0.001_blay_tree.rds")
  E11tree_peri$count = dcast(E11tree_peri$cm,tags~celltype)
  e11cc = E11tree_peri$count
  rownames(e11cc) = e11cc$tags;e11cc = e11cc[-1]
  e11cc = e11cc[which(rowSums(e11cc)>1),]
  e11cc = as.data.frame(t(t(as.matrix(e11cc))/colSums(e11cc)))
  e11cc = e11cc[,corder11n]
  e11ccor = e11cc;e11ccord = apply(e11cc,1,which.max);
  for (i in 1:nrow(e11ccor)) {
    e11ccor[i,] = 0
    e11ccor[i,e11ccord[i]] = e11cc[i,e11ccord[i]]
  }
  e11cc = e11cc[do.call(order, -e11ccor),]
  adid11 = rownames(e11cc)[which(!rownames(e11cc) %in% names(e11_clu$cluster))]
  tmp11 = e11cc[which(rownames(e11cc) %in% adid11),]
  # e15_clu$cellcounts = e15_clu$cellcounts[which(rownames(e15_clu$cellcounts) %in% names(nsplit)),]
  e11add.hmap <- Heatmap(as.matrix(tmp11), name = "filtered clones' cell portion",
                         row_order = rownames(tmp11),
                         # col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
                         show_row_names = F,column_order = corder11n)
  e11add.hmap
  
  
  Heatmap(as.matrix(e11addls$mx), name = "filtered clones' cell portion",
          row_order = e11addls$roworder,
          # col = c("#08306B",colorRampPalette(brewer.pal(9, "Reds"))(10)),
          show_row_names = F,column_order = e11addls$colorder)
  
  e11addls = list("mx" = tmp11,"roworder" = rownames(tmp11),"colorder" = corder11n,"fg" = e11add.hmap)
  saveRDS(e11addls,"final_result/xdraft/main_figure_2_25/E11_tag_cluster_heatmap_counts_arti_supplement.rds")
  
  ggexport(e11add.hmap, filename = "final_result/xdraft/main_figure_2_25/E11_tag_cluster_heatmap_counts_arti_supplement.pdf",width = 8,height = 8)
  
  
  
  
}

#-------------7.0 NbFP direction-------------
{
  e15_clu$cluster.fig
  #DA Glu predict only c7-----------------
  cmhub = E15tree$cm[which(E15tree$cm$tags %in% c(names(e15_clu$cluster[which(e15_clu$cluster %in% c("c6"))]),adid15)),]
  E15trans = readRDS("raw_data/xie/Dual_E15_rerun210904.rds")
  
  hubc = c("Rgl1", "DA", "GluFP", "NbDA", "NbGlu", "NbFP")
  hubcm = cmhub[which(cmhub$celltype %in% hubc),]
  hubst = dcast(hubcm,tags~celltype)
  rownames(hubst) = hubst$tags;hubst = hubst[-1]
  
  merhub = hubst
  merhub$DAt = merhub$DA + merhub$NbDA
  merhub$Glut = merhub$NbGlu + merhub$GluFP
  merhub = merhub[c("NbFP","Rgl1","DAt","Glut")]
  
  merhub.up = as.matrix(merhub)
  merhub.up[which(merhub.up>0)] = T
  merhub.up[which(merhub.up==0)] = F
  
  library(UpSetR)
  p6.1 = upset(as.data.frame(merhub.up),nsets = 4)
  p6.1
  # ggexport(p6.1,filename = "final_result/xE15.5/Nb_direct_stat_upset_2_23.pdf",width = 8,height = 6)
  merhubf = merhub[which(rowSums(merhub[3:4])>0),1:4]
  merhubf.up = as.matrix(merhubf)
  merhubf.up[which(merhubf.up>0)] = T
  merhubf.up[which(merhubf.up==0)] = F
  p6.2 = upset(as.data.frame(merhubf.up),nsets = 4, keep.order = T)
  p6.2
  ggexport(p6.2,filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_stat_upset.pdf",width = 8,height = 6)
  
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
  merhubfl$cell = factor(merhubfl$cell,levels = c("NbFP","Glut","DAt","Rgl1"))
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
  ggsave(p6.3,filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_stat_point.pdf",width = 10,height = 3)

  {
    hubcb = names(E15trans$ClusterName[which(E15trans$ClusterName %in% hubc)])
    hubcb = substr(hubcb,1,16)
    trans.da = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.da)),"Var1"]
    trans.glu = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.glu)),"Var1"]
    trans.all = E15tree$cm[which(E15tree$cm$tags %in% rownames(group.all)),"Var1"]
    cell_names.da <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.da]
    cell_names.glu <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.glu]
    cell_names.all <- colnames(E15trans)[substr(colnames(E15trans),1,16) %in% trans.all]
    
    p6.4.0 = DimPlot(E15trans,label = TRUE,label.size = 3) + NoLegend()
    p6.4.1 = DimPlot(E15trans,label = TRUE, label.size = 3,cells.highlight = cell_names.da, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.2 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.glu, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    p6.4.3 = DimPlot(E15trans,label = TRUE,label.size = 3,cells.highlight = cell_names.all, 
                     cols.highlight = "red", cols = "gray", order = TRUE) + NoLegend()
    
    p6.4 = p6.4.0 + p6.4.1 + p6.4.2 + p6.4.3
    p6.4
    ggsave(p6.4,filename = "final_result/xdraft/E15.5_DA_Glu_direction_umap_example.pdf",width = 10,height = 8)
    
    
    # E15trans.hub = subset(E15trans, idents = c("Rgl1"))
    E15trans.hub = subset(E15trans, idents = c("Nb_M","DA","Nb_DA","Glu_ML"))
    E15trans.hub = subset(E15trans, idents = c("Nb_M"))
    tmp = Idents(E15trans.hub)
    tmp = as.character(tmp);names(tmp) = names(Idents(E15trans.hub))
    tmp[which(substr(names(tmp),1,16) %in% trans.da)] = "DA_d"
    tmp[which(substr(names(tmp),1,16) %in% trans.glu)] = "Glu_d"
    tmp[which(substr(names(tmp),1,16) %in% trans.all)] = "All_d"
    tmp = as.factor(tmp)
    Idents(E15trans.hub) = tmp
    
    E15trans.hub = subset(E15trans.hub, idents = c("DA_d","Glu_d","All_d"))
    
    DA.markers <- FindMarkers(E15trans.hub, ident.1 = "DA_d",min.pct = 0.25)
    Glu.markers <- FindMarkers(E15trans.hub, ident.1 = "Glu_d", min.pct = 0.25)
    All.markers <- FindMarkers(E15trans.hub, ident.1 = "All_d", min.pct = 0.25)
    levels(E15trans.hub) = c("DA_d","All_d","Glu_d")

    write.csv(DA.markers,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_DA_NbFP.csv",row.names = T,quote = F)
    write.csv(Glu.markers,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_Glu_NbFP.csv",row.names = T,quote = F)
    write.csv(All.markers,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_All_NbFP.csv",row.names = T,quote = F)
 
    dagluid = as.data.frame(Idents(E15trans.hub))
    write.csv(dagluid,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_idents_total.csv")
    
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
    
    p6.6.3.0 = dens1 + plot_spacer() + p6.6.3.0 + dens2 + 
      plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    p6.6.3.0
    
    library(ggpubr)
    p6.6.3.1 = ggplot(pcmx,aes(y = PC_1, x = class)) + 
      geom_boxplot(width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = list(c("DA_d","All_d"), c("DA_d","Glu_d"), c("All_d","Glu_d")),
                         label = "p.signif", method = "t.test") +
      theme(legend.position = c(5, 20))
    
    p6.6.3.2 = ggplot(pcmx,aes(y = PC_2, x = class)) + 
      geom_boxplot(width = 0.5) + theme_bw() +
      stat_compare_means(comparisons = list(c("DA_d","All_d"), c("DA_d","Glu_d"), c("All_d","Glu_d")),
                         label = "p.signif", method = "t.test")
    p6.6.1
    ggsave(p6.6.3.0, filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_PCA_onlyNbFP.pdf",width = 8,height = 6)
    p6.6.2
    ggsave(p6.6.2, filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_PC12_gene_onlyNbFP.pdf",width = 7,height = 4)
    p6.6.3.1+p6.6.3.2
    ggsave(p6.6.3.1+p6.6.3.2, 
           filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_PC12_class_boxplot_onlyNbFP.pdf",
           width = 6,height = 4)
    
    saveRDS(E15trans.hub,"final_result/xdraft/main_figure_2_25/E15trans_DA_onlyNbFP.rds")
  }

  
}

#-------------7.1 NbFP direction NW model-------------
{
  {
    
    E15trans.hubdg = subset(E15trans.hub, idents = c("DA_d","Glu_d"))
    E15trans.hubdg = FindVariableFeatures(E15trans.hubdg, selection.method = "vst", nfeatures = 32285)
    DGdegs = FindMarkers(E15trans.hubdg,ident.1 = "DA_d",min.pct = 0.25)
    write.csv(DGdegs,"final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_direction_DEG_analysis_DAvsGlu.csv",row.names = T,quote = F)
    
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
    ggsave(p7,filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_neuralnet_auc_boxplot_onlyNbFP.pdf",width = 4,height = 4)
    aucst = aucdf %>% group_by(markers) %>% summarise(mean(AUC))
    write.csv(aucdf,file = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_neuralnet_auc.csv",row.names = F,quote = F)
    write.csv(aucst,file = "final_result/xdraft/main_figure_2_25/E15.5_DA_Glu_neuralnet_auc_stat.csv",row.names = F,quote = F)
    
    
  }
  
  {
    count1 = Divtree$count
    count2 = E11stree$count
    
    colnames(count1)[-1] = paste0("DIV7_",colnames(count1)[-1])
    colnames(count2)[-1] = paste0("E11_",colnames(count2)[-1])
    e11divcms = merge(count2,count1,by = "tags")
    
    e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$tags;e11divcms.pd = e11divcms.pd[-1]
    
    STFtrans = readRDS("/picb/sysgenomics2/people/liuhengxin/P6_lineartree/raw_data/xie/STF_E11_singlet_renamed_20211221.rds")
    Gluclone = e11divcms.pd[which(e11divcms.pd$E11_NbFP_0 > 0 & rowSums(e11divcms.pd[c("DIV7_NbGlu")]) > 0 & rowSums(e11divcms.pd[c("DIV7_DA","DIV7_NbDA")])==0),]
    DAclone = e11divcms.pd[which(e11divcms.pd$E11_NbFP_0 > 0 & rowSums(e11divcms.pd[c("DIV7_DA","DIV7_NbDA")])>0 & rowSums(e11divcms.pd[c("DIV7_NbGlu")]) == 0),]
    
    GluBC = E11stree$cm[which(E11stree$cm$tags %in% rownames(Gluclone) & E11stree$cm$celltype == "NbFP_0"),"Var1"]
    DABC = E11stree$cm[which(E11stree$cm$tags %in% rownames(DAclone) & E11stree$cm$celltype == "NbFP_0"),"Var1"]
    
    tmp = colnames(STFtrans)[which(substr(colnames(STFtrans),1,16)%in% GluBC)]
    Glutrans = subset(STFtrans,cells = tmp)
    
    tmp = colnames(STFtrans)[which(substr(colnames(STFtrans),1,16)%in% DABC)]
    DAtrans = subset(STFtrans,cells = tmp)
    Sdamx = as.data.frame(DAtrans[["RNA"]]@counts)
    Sglumx = as.data.frame(Glutrans[["RNA"]]@counts)
    
    
    #test model
    markers = marker1
    
    
    #calculate roc
    {
      NeurBuild2 = function(markers){
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
        
        Sdamxm = Sdamx[rownames(Sdamx) %in% markers,]
        Sglumxm = Sglumx[rownames(Sglumx) %in% markers,]
        Sdamxm = as.data.frame(t(Sdamxm))
        Sglumxm = as.data.frame(t(Sglumxm))
        Normmx = function(x){
          return((x-min(x))/(max(x)-min(x)))
        }
        Sdamxmn =  as.data.frame(t(apply(Sdamxm, 1, Normmx)))
        Sglumxmn =  as.data.frame(t(apply(Sglumxm, 1, Normmx)))
        Sdamxmn$type =  0
        Sglumxmn$type = 1
        
        Sdadata = rbind(Sdamxmn,Sglumxmn)
        colnames(Sdadata)[1:(ncol(Sdadata)-1)] = paste0("x",1:(ncol(Sdadata)-1))
        
        testdata = Sdadata
        traindata = dadata
        
        f <- as.formula(paste0("type ~ ",
                               paste(names(traindata)[1:(ncol(traindata) - 1)],
                                     collapse = " + ")))
        library(neuralnet)
        neur <- neuralnet(f, data = traindata, hidden = 10,
                          act.fct = "logistic", linear.output = T,
                          err.fct = "sse", rep = 1)
        
        test.hat = neuralnet::compute(neur,testdata[-ncol(testdata)])$net.result
        test.predict = test.hat[,1]
        modelroc = roc(testdata$type,test.predict)
        
        return(list("neur.model" = neur,"test.hat" = test.hat, "roc" = modelroc))
      }
      
      
      aucdfdiv = NULL
      for (i in 1:100) {
        roc1 = NeurBuild2(marker1)
        roc2 = NeurBuild2(marker2)
        roc3 = NeurBuild2(marker3)
        roc4 = NeurBuild2(marker4)
        aucdfdiv = rbind(aucdfdiv,data.frame(sample = i, markers = c("Differentially expressed",
                                                               "Highly variable genes",
                                                               "Random genes",
                                                               "Transcription factors"),
                                       AUC = c(roc1$roc$auc,roc2$roc$auc,roc3$roc$auc,roc4$roc$auc)))
      }
      
      aucdfdiv$markers = factor(aucdfdiv$markers,levels = c("Differentially expressed",
                                                      "Highly variable genes",
                                                      "Random genes",
                                                      "Transcription factors"))
      p8 = ggboxplot(aucdfdiv, x = "markers", y = "AUC",
                     fill = "markers", palette = "jco",
                     width = 0.2) + 
        theme_bw() + theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           legend.position = "none",
                           axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
        stat_compare_means(method = "anova", label.y = 1) +
        stat_compare_means(label = "p.signif", method = "t.test",
                           ref.group = "Differentially expressed")    
      p8
      ggsave(p8,filename = "final_result/xdraft/main_figure_2_25/STF_E11-DIV_DA_Glu_neuralnet_auc_boxplot_DIV7.pdf",width = 4,height = 4)
      write.csv(aucdf,"final_result/xdraft/main_figure_2_25/STF_E11-DIV_DA_Glu_neuralnet_auc_boxplot_DIV7.csv")
      
      aucdf %>% group_by(markers) %>% summarise(mean(AUC))  
      
    }
  }
  
  #DA portion
  {
    cellcount = E15tree$count[which(E15tree$count$tags %in% cmhub$tags),]
    nbfpp = cellcount[which(cellcount$NbFP > 0 & (cellcount$DA + cellcount$NbDA) > 0 &
                              (cellcount$GluFP + cellcount$NbGlu) == 0),]
    nbfpm = cellcount[which(cellcount$NbFP == 0 & (cellcount$DA + cellcount$NbDA) > 0 &
                              (cellcount$GluFP + cellcount$NbGlu) == 0),]
    nbfpp$DA_percentage = (nbfpp$DA)/(nbfpp$DA + nbfpp$NbDA)
    nbfpm$DA_percentage = (nbfpm$DA)/(nbfpm$DA + nbfpm$NbDA)
    nbfpst = data.frame("DA_percentage" = c(nbfpp$DA_percentage,nbfpm$DA_percentage),
                        "type" = rep(c("NBFP+","NBFP-"),c(nrow(nbfpp),nrow(nbfpm))))
    
    my_comparisons <- list( c("NBFP+", "NBFP-"))
    p4 = ggboxplot(nbfpst, x = "type", y = "DA_percentage",
                   fill = "type", palette = "jco",
                   width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(label = "p.signif", method = "t.test",label.y = 1,comparisons = my_comparisons)    
    p4
    ggsave(p4,filename = "final_result/xdraft/main_figure_2_25/E15.5_DA_percentage_withornot_NbFP.pdf",width = 4,height = 4)
    write.csv(nbfpst,"final_result/xdraft/main_figure_2_25/E15.5_DA_percentage_withornot_NbFP.csv")
    
    nbfpgp = cellcount[which(cellcount$NbFP > 0 & (cellcount$DA + cellcount$NbDA) == 0 &
                               (cellcount$GluFP + cellcount$NbGlu) > 0),]
    nbfpgm = cellcount[which(cellcount$NbFP == 0 & (cellcount$DA + cellcount$NbDA) == 0 &
                               (cellcount$GluFP + cellcount$NbGlu) > 0),]
    nbfpgp$Glu_percentage = (nbfpgp$GluFP)/(nbfpgp$GluFP + nbfpgp$NbGlu)
    nbfpgm$Glu_percentage = (nbfpgm$GluFP)/(nbfpgm$GluFP + nbfpgm$NbGlu)
    nbfpgst = data.frame("Glu_percentage" = c(nbfpgp$Glu_percentage,nbfpgm$Glu_percentage),
                         "type" = rep(c("NBFP+","NBFP-"),c(nrow(nbfpgp),nrow(nbfpgm))))
    
    my_comparisons <- list( c("NBFP+", "NBFP-"))
    p4.1 = ggboxplot(nbfpgst, x = "type", y = "Glu_percentage",
                     fill = "type", palette = "jco",
                     width = 0.2) + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 45,vjust = 0.8,size=8,hjust = 0.8)) +
      xlab("") +
      stat_compare_means(label = "p.signif", method = "t.test",label.y = 1,comparisons = my_comparisons)    
    p4.1
    ggsave(p4.1,filename = "final_result/xdraft/main_figure_2_25/E15.5_Glu_percentage_withornot_NbFP.pdf",width = 4,height = 4)
    write.csv(nbfpgst,"final_result/xdraft/main_figure_2_25/E15.5_Glu_percentage_withornot_NbFP.csv",row.names = F)
    
    
    
  }
  
}

#-------------8. E11 sister clone similarity-------------
{
  
  E11trans = readRDS("raw_data/xie/RetroE11_renamed.RDS")
  
  #build data
  e11cm = E11tree$cm
  cmrid = sample(1:nrow(e11cm),nrow(e11cm)/2)
  e11cmar = e11cm[cmrid,]
  e11cmbr = e11cm[setdiff(1:nrow(e11cm),cmrid),]
  cmtagc = intersect(e11cmar$tags,e11cmbr$tags)
  e11cma = dcast(e11cmar,celltype~tags)
  e11cmb = dcast(e11cmbr,celltype~tags)
  rownames(e11cma) = e11cma$celltype;e11cma = e11cma[-1]
  rownames(e11cmb) = e11cmb$celltype;e11cmb = e11cmb[-1]

  #RNA data
  {
    # e11clrnaa
    cmida = match(cmtagc,colnames(e11cma))
    cmidb = match(cmtagc,colnames(e11cmb))
    smpname = c(colnames(e11cma),colnames(e11cmb))
    
    colnames(e11cma) = paste0("a-",colnames(e11cma))
    colnames(e11cmb) = paste0("b-",colnames(e11cmb))
    cellseur <- CreateSeuratObject(counts = cbind(e11cma,e11cmb), min.cells = 0)
    group = rep(c("A","B"),c(ncol(e11cma),ncol(e11cmb)))
    names(group) = colnames(cbind(e11cma,e11cmb))
    cellseur$group = group
    cellseur <- ScaleData(cellseur)
    cellseur <- RunPCA(cellseur, features = rownames(cellseur))
    cellseur <- FindNeighbors(cellseur, dims = 1:5)
    # cellseur <- FindClusters(cellseur, resolution = 0.8)
    cellseur <- FindClusters(cellseur, resolution = 0.1)
    cellseur <- RunUMAP(cellseur, dims = 1:10)
    
    dimcluster = e11_clu$cluster
    dimcluster = dimcluster[which(names(dimcluster) %in% smpname)]
    dimcluster = dimcluster[smpname]
    clonecluster = as.character(c(dimcluster))
    names(clonecluster) = colnames(cbind(e11cma,e11cmb))
    cellseur$clonecluster = clonecluster
    
    p5.1 = DimPlot(cellseur) | DimPlot(cellseur,split.by = "group",group.by = "clonecluster")
    p5.1
  }
  
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
        scale_color_manual(values = c("#D9534F","#FFAD60")) + 
        labs(color='group') +
        xlab("Cluster distance") + ylab("ecdf of paired clone") +
        theme_bw()
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
      p5 = list(p5.1,p5.2,p5.2.1,p5.3,p5.4,p5.5,p5.6)
      datals = list("funseur"= funseur,"funump" = funump, "distdf" = distdf,
                    "cordf" = cordf, "pairstatl" = pairstatl,"clsispmx" = clsispmx)
      resultls = list("data" = datals,"fg" = p5)
      
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
      smpnamea = colnames(dcast(e11cmar,celltype~tags))[-1]
      smpnameb = colnames(dcast(e11cmbr,celltype~tags))[-1]
      linebca = as.character(e11cmar[which(e11cmar$tags == smpnamea[1]),"Var1"])
      lclrnaa = e11rna[,which(substr(colnames(e11rna),1,16) %in% linebca)]
      if(!is.null(ncol(lclrnaa))){
        lclrnaas = as.data.frame(rowMeans(lclrnaa))
      }else{
        lclrnaas = as.data.frame(lclrnaa)
      }
      e11clrnaa = lclrnaas
      colnames(e11clrnaa)[1] = smpnamea[1]
      
      for (i in 2:length(smpnamea)) {
        print(i)
        linebca = as.character(e11cmar[which(e11cmar$tags == smpnamea[i]),"Var1"])
        lclrnaa = e11rna[,which(substr(colnames(e11rna),1,16) %in% linebca)]
        if(!is.null(ncol(lclrnaa))){
          lclrnaas = as.data.frame(rowMeans(lclrnaa))
        }else{
          lclrnaas = as.data.frame(lclrnaa)
        }
        e11clrnaa = cbind(e11clrnaa,lclrnaas)
        colnames(e11clrnaa)[i] = smpnamea[i]
      }
      
      linebcb = as.character(e11cmbr[which(e11cmbr$tags == smpnameb[1]),"Var1"])
      lclrnab = e11rna[,which(substr(colnames(e11rna),1,16) %in% linebcb)]
      if(!is.null(ncol(lclrnab))){
        lclrnabs = as.data.frame(rowMeans(lclrnab))
      }else{
        lclrnabs = as.data.frame(lclrnab)
      }
      e11clrnab = lclrnabs
      colnames(e11clrnab)[1] = smpnameb[1]
      for (i in 2:length(smpnameb)) {
        print(i)
        linebcb = as.character(e11cmbr[which(e11cmbr$tags == smpnameb[i]),"Var1"])
        lclrnab = e11rna[,which(substr(colnames(e11rna),1,16) %in% linebcb)]
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
    
    e11rna = as.data.frame(GetAssayData(E11trans,assay = "RNA"))
    # cmtagr = intersect(smpnamea,smpnameb)
    e11cmarp = e11cmar[which(e11cmar$celltype %in% c("Rgl1","NProgBL","NProgBM")),]
    e11cmbrp = e11cmbr[which(e11cmbr$celltype %in% c("Rgl1","NProgBL","NProgBM")),]
    cmtagrp = intersect(unique(e11cmarp$tags), unique(e11cmbrp$tags))
    e11cmaro = e11cmar[which(!e11cmar$celltype %in% c("Rgl1","NProgBL","NProgBM")),]
    e11cmbro = e11cmbr[which(!e11cmbr$celltype %in% c("Rgl1","NProgBL","NProgBM")),]
    cmtagro = intersect(unique(e11cmaro$tags), unique(e11cmbro$tags))
    
    transeurls = ClonetoCell(e11cmar,e11cmbr,cmtagc)
    transeurpls = ClonetoCell(e11cmarp,e11cmbrp,cmtagrp)
    transeurols = ClonetoCell(e11cmaro,e11cmbro,cmtagro)
    
    transeurls$seur <- FindClusters(transeurls$seur, resolution = 0.8)
    transeurls$seur <- RunUMAP(transeurls$seur, dims = 1:10)
    
    transeurpls$seur <- FindClusters(transeurpls$seur, resolution = 0.8)
    transeurpls$seur <- RunUMAP(transeurpls$seur, dims = 1:10)
    
    transeurols$seur <- FindClusters(transeurols$seur, resolution = 0.8)
    transeurols$seur <- RunUMAP(transeurols$seur, dims = 1:10)
    
    
    cellcmpls = SisterUmapAna(cellseur, e11cma, e11cmb,
                              cmida,cmidb,cmtagc)
    transcmpls = SisterUmapAna(transeurls$seur,transeurls$cma,transeurls$cmb,
                               transeurls$ida,transeurls$idb,cmtagc)
    transcmppls = SisterUmapAna(transeurpls$seur,transeurpls$cma,transeurpls$cmb,
                                transeurpls$ida,transeurpls$idb,cmtagrp)
    transcmpols = SisterUmapAna(transeurols$seur,transeurols$cma,transeurols$cmb,
                                transeurols$ida,transeurols$idb,cmtagro)
    ggsave(cellcmpls$fg[[1]],filename = "final_result/xdraft/main_figure_2_25/E11_sister_clone_cellcomp_umap.pdf",width = 10,height = 5)
    pdf("final_result/xdraft/main_figure_2_25/E11_sister_clone_cellcomp_cmp.pdf",width = 7,height = 5)
    invisible(lapply(cellcmpls$fg[-1], print))
    dev.off()
    cellcmpls$data = cellcmpls$data[-1]
    transcmpls$data = transcmpls$data[-1]
    transcmppls$data = transcmppls$data[-1]
    transcmpols$data = transcmpols$data[-1]
    saveRDS(cellcmpls$data,"final_result/xdraft/main_figure_2_25/E11_sister_cellcmpls_data.rds")
    
    ggsave(transcmpls$fg[[1]],filename = "final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_umap.pdf",width = 10,height = 5)
    pdf("final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_cmp.pdf",width = 7,height = 5)
    invisible(lapply(transcmpls$fg[-1], print))
    dev.off()
    saveRDS(transcmpls,"final_result/xdraft/main_figure_2_25/E11_sister_transcmpls.rds")
    saveRDS(transcmpls$data,"final_result/xdraft/main_figure_2_25/E11_sister_transcmpls_data.rds")
    
    ggsave(transcmppls$fg[[1]],filename = "final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_umap_prog.pdf",width = 10,height = 5)
    pdf("final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_cmp_prog.pdf",width = 7,height = 5)
    invisible(lapply(transcmppls$fg[-1], print))
    dev.off()
    saveRDS(transcmppls$data,"final_result/xdraft/main_figure_2_25/E11_sister_cellcmpls_prog_data.rds")
    
    ggsave(transcmpols$fg[[1]],filename = "final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_umap_other.pdf",width = 10,height = 5)
    pdf("final_result/xdraft/main_figure_2_25/E11_sister_clone_trans_cmp_other.pdf",width = 7,height = 5)
    invisible(lapply(transcmpols$fg[-1], print))
    dev.off()
    saveRDS(transcmpols$data,"final_result/xdraft/main_figure_2_25/E11_sister_cellcmpls_other_data.rds")
    
  }

  
  
} 

#-------------9. DIV E11 clone state and pie plot DEG-------------
{
  corderdiv = c('GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2','GluFP','NbGlu','NbFP','NbDA','DA',
                "GabaUnknown","Rgl2","Rgl3","Unknown")
  corder11n = c("NbGluAL","NProgBL","GabaProg0","GabaProgBL","GabaProgBI","NProgBM","OMTN","NbBM0",
                "NbBM1","Rgl1","NbFP_0")
  Divtrans = readRDS("raw_data/xie/RetroDIV7_renamed_0904.RDS")
  count1 = Divtree$count
  count2 = E11stree$count
  
  count1 = count1[,c("tags",corderdiv)]
  count2 = count2[,c("tags",corder11n)]
  
  colnames(count1)[-1] = paste0("DIV7_",colnames(count1)[-1])
  colnames(count2)[-1] = paste0("E11_",colnames(count2)[-1])
  e11divcms = merge(count2,count1,by = "tags")
  
  e11divcms.pd = e11divcms;rownames(e11divcms.pd) = e11divcms.pd$tags;e11divcms.pd = e11divcms.pd[-1]
  e11divcms.pdn = as.data.frame(t(t(e11divcms.pd)/colSums(e11divcms.pd)))
  
  
  #DIV sclone cluster
  div_cl = readRDS("final_result/xdraft/main_figure_2_25/div_w_o_unknown_cellcounts_cluster_C9_20220227_V2.rds")
  div_cl$cluster.fig
  div_cl$cluster[which(names(div_cl$cluster) %in% rownames(edrgl1))]
  div_add = Divtree$count[which(!Divtree$count$tags %in% names(div_cl$cluster)),]
  div_addcl = colnames(div_add[-1])[apply(div_add[-1], 1, which.max)]
  maxcount = NULL
  for (i in 1:nrow(div_add)) {
    line = div_add[i,]
    maxcount =c(maxcount, unlist(line[div_addcl[i]]))
  }
  div_addcl = div_addcl[which(maxcount>1)]
  div_add = div_add[which(maxcount>1),]
  
  div_add = data.frame("clone" = div_add$tags, "group" = div_addcl)
  div_add[which(div_add$group %in% c('GluAL1')),"group"] = "uni-AL"
  div_add[which(div_add$group %in% c('GabaBL2','GabaBL4','GabaBL1')),"group"] = "uni-BL"
  div_add[which(div_add$group %in% c('GabaBI')),"group"] = "uni-BI"
  div_add[which(div_add$group %in% c('GluBM1','GluBM2','GluFP')),"group"] = "uni-BM"
  div_add[which(div_add$group %in% c('NbGlu','NbFP','NbDA','DA')),"group"] = "uni-FP"
  div_add = div_add[which(!div_add$group %in% c("GabaUnknown","Rgl2","Rgl3","Unknown")),]
  div_clg = data.frame("clone" = names(div_cl$cluster),"group" = as.character(div_cl$cluster))
  div_clg = rbind(div_clg,div_add)
  div_clg$group = factor(div_clg$group,levels = c(paste0("c",c(1:9)),"uni-AL","uni-BL","uni-BI","uni-BM","uni-FP"))
  
  #build prog data
  # edrprg = e11divcms.pdn[,1:11]
  edrprg = e11divcms.pdn[,c("E11_Rgl1", "E11_NProgBL", "E11_NProgBM")]
  edrprg = edrprg/rowSums(edrprg)
  thes = 1
  edrgl1 = e11divcms.pdn[which(edrprg$E11_Rgl1 >= thes  & edrprg$E11_NProgBL == 0 & edrprg$E11_NProgBM == 0),]
  edbl = e11divcms.pdn[which(edrprg$E11_Rgl1 == 0 & edrprg$E11_NProgBL >= thes & edrprg$E11_NProgBM == 0),]
  edbm = e11divcms.pdn[which(edrprg$E11_Rgl1 == 0 & edrprg$E11_NProgBL == 0 & edrprg$E11_NProgBM >= thes),]
  #merge cell option
  {
    edrgl1$DAt = edrgl1$DIV7_DA + edrgl1$DIV7_NbDA
    edrgl1$Glut = edrgl1$DIV7_GluFP + edrgl1$DIV7_NbGlu
    edrgl1$other = rowSums(edrgl1[,paste0("DIV7_",c('GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI','GluBM1','GluBM2',
                                                    "GabaUnknown","Rgl2","Rgl3","Unknown"))])
    edrgl1 = edrgl1[,c("E11_Rgl1","DIV7_GluAL1","DAt","Glut","other")]
    colnames(edrgl1) = c("Rgl1","GluAL1","DAt","Glut","other")
    
    
    edbl$other = rowSums(edbl[,paste0("DIV7_",c('GluBM1','GluBM2',
                                                "GabaUnknown","Rgl2","Rgl3","Unknown"))])
    edbl = edbl[,c("E11_NProgBL",'DIV7_GluAL1','DIV7_GabaBL1','DIV7_GabaBL2','DIV7_GabaBL4','DIV7_GabaBI',"other")]
    colnames(edbl) = c("NProgBL",'GluAL1','GabaBL1','GabaBL2','GabaBL4','GabaBI',"other")
    
    edbm$DAt = edbm$DIV7_DA + edbm$DIV7_NbDA
    edbm$Glut = edbm$DIV7_GluFP + edbm$DIV7_NbGlu
    edbm$other = rowSums(edbm[,paste0("DIV7_",c('GluAL1','GabaBL1','GabaBL2','GabaBL4',
                                                "GabaUnknown","Rgl2","Rgl3","Unknown"))])
    edbm = edbm[,c("E11_NProgBM",'DIV7_GabaBI','DIV7_GluBM1','DIV7_GluBM2',"DAt","Glut","other")]
    colnames(edbm) = c("NProgBM",'GabaBI','GluBM1','GluBM2',"DAt","Glut","other")
  }
  
  
  #stat
  aimc = "NProgBM"
  PieStat = function(div_clg,edrgl1,aimc){
    edcp = div_clg[which(div_clg$clone %in% rownames(edrgl1)),]
    edcps = E11stree$cm[which(E11stree$cm$celltype == aimc & E11stree$cm$tags %in% edcp$clone),]
    edcps = merge(edcps,div_clg,by.x = "tags",by.y = "clone")
    
    plcps = as.data.frame(table(edcps$group))
    plcps$prop = round(plcps$Freq/sum(plcps$Freq),3)
    plcps = plcps[which(plcps$prop > 0),]
    levels(plcps$Var1) = c(levels(plcps$Var1),"others")
    plcps[which(plcps$prop < 0.05),"Var1"] = "others"
    plcps = plcps %>% group_by(Var1) %>% summarise(prop = sum(prop))
    plcps$group = plcps$Var1
    labs <- paste0(plcps$group, " (", plcps$prop*100, "%)")
    
    p9 = ggdonutchart(plcps, x = "prop", label = rev(labs), color = "white",
                      lab.pos = "in", lab.font = "black",
                      fill = "group") + NoLegend() + ggtitle(paste0(aimc," cellnum = ",nrow(edcps)))
    return(list("data" = edcps[-3], "stat" = plcps,"piefig" = p9))
    
  }
  PieStat2 = function(edrgl1){
    edrgl1cp = div_clg[which(div_clg$clone %in% rownames(edrgl1)),]
    edrgl1cps = table(edrgl1cp$group)
    edrgl1cps$prop = round(edrgl1cps$Freq/sum(edrgl1cps$Freq),3)
    edrgl1cps = edrgl1cps[which(edrgl1cps$prop > 0),]
    levels(edrgl1cps$edrgl1cp) = c(levels(edrgl1cps$edrgl1cp),"others")
    edrgl1cps[which(edrgl1cps$prop < 0.05),"edrgl1cp"] = "others"
    edrgl1cps = edrgl1cps %>%group_by(edrgl1cp) %>% summarise(prop = sum(prop))
    edrgl1cps$group = edrgl1cps$edrgl1cp
    labs <- paste0(edrgl1cps$group, " (", edrgl1cps$prop*100, "%)")
    
    p9 = ggdonutchart(edrgl1cps, x = "prop", label = rev(labs), color = "white",
                      lab.pos = "in", lab.font = "black",
                      fill = "group") + NoLegend()
    return(list("stat" = edrgl1cps,"piefig" = p9))
    
  }
  
  piergl1 = PieStat(div_clg,edrgl1,"Rgl1")
  piebl = PieStat(div_clg,edbl,"NProgBL")
  piebm = PieStat(div_clg,edbm,"NProgBM")
  p9 = piergl1$piefig +
  piebm$piefig +
  piebl$piefig
  p9
  ggsave(p9,filename = "final_result/xdraft/main_figure_2_25/STF_E11_DIV7_prog_clone_cluster_pie.pdf",width = 10,height = 5)
  
  ProgPreDEG = function(piergl1){
    # piergl1$data$CB = colnames(STFtrans)[match(piergl1$data$Var1,substr(colnames(STFtrans),1,16))]
    rglda = Divtree$cm[which(Divtree$cm$tags %in% piergl1$data$tags & Divtree$cm$celltype == "DA"),]
    rglda = merge(rglda,unique(piergl1$data[c(1,4)]),by = "tags",all = F)
    # rglda$CB = colnames(STFtrans)[match(rglda$Var1,substr(colnames(STFtrans),1,16))]
    
    rglga = Divtree$cm[which(Divtree$cm$tags %in% piergl1$data$tags & Divtree$cm$celltype == "GabaBI"),]
    rglga = merge(rglga,unique(piergl1$data[c(1,4)]),by = "tags",all = F)
    # rglga$CB = colnames(STFtrans)[match(rglga$Var1,substr(colnames(STFtrans),1,16))]
    return(list("prog" = piergl1$data,"DA" = rglda[-3],"GabaBI" = rglga[-3]))
  }
  
  rgl1degd = ProgPreDEG(piergl1)
  bmdegd = ProgPreDEG(piebm)
  bldegd = ProgPreDEG(piebl)
  degd = list("Rgl1" = rgl1degd,"BM" = bmdegd,"BL" = bldegd)
  saveRDS(degd,"final_result/xdraft/main_figure_2_25/STF_E11_DIV7_prog_clone_cluster_degpredata.rds")
  head(degd$Rgl1$prog)
  head(degd$Rgl1$DA)
  
  #clone cell type stat
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
      ggtitle(paste0(title," clonesize = ",nrow(edrgl1)))
    p9.2
    return(p9.2)
  }
  colorRampPalette(brewer.pal(9, "Blues"))(10)
  e15cts = dcast(E15tree$cm,tags~celltype)
  e15cts = e15cts[which(rowSums(e15cts[-1]) > 1),]
  e11cts = dcast(E11tree$cm,tags~celltype)
  e11cts = e11cts[which(rowSums(e11cts[-1]) > 1),]
  
  p9.1.1 = CloneCnum(e11cts[-1],"E11")
  p9.1.2 = CloneCnum(e15cts[-1],"E15")
  p9.1 = p9.1.1 + p9.1.2
  p9.1
  ggsave(p9.1,filename = "final_result/xdraft/main_figure_2_25/E11_E15_clone_cluster_celltype_num_stat.pdf",width = 9,height = 6)
  
  p9.2.1 = ClonePie(e11cts[-1],"E11")
  p9.2.2 = ClonePie(e15cts[-1],"E15")
  p9.2 = p9.2.1 + p9.2.2
  p9.2
  ggsave(p9.2,filename = "final_result/xdraft/main_figure_2_25/E11_E15_clone_cluster_celltype_num_stat_pie.pdf",width = 9,height = 6)
  
  
  
  print()
  
  
}


#-------------10. Rgl1/Prog DA trans similarity-------------
{
  rglda = read.csv("final_result/xdraft/DA_from_rgl1.csv")
  pgda = read.csv("final_result/xdraft/DA_from_BM.csv")
  dapg = read.csv("final_result/xdraft/NProgBM_DA.csv")
  otpg = read.csv("final_result/xdraft/NProgBM_others.csv")
  Divtrans = readRDS("raw_data/xie/STF_DIV7_scar_added.rds")
  
  dapg.tr = subset(STFtrans, cells = colnames(STFtrans)[substr(colnames(STFtrans),1,16) %in% dapg$Var1])
  otpg.tr = subset(STFtrans, cells = colnames(STFtrans)[substr(colnames(STFtrans),1,16) %in% otpg$Var1])
  rglda.tr = subset(Divtrans, cells = colnames(Divtrans)[substr(colnames(Divtrans),1,16) %in% rglda$Var1])
  pgda.tr = subset(Divtrans, cells = colnames(Divtrans)[substr(colnames(Divtrans),1,16) %in% pgda$Var1])
  dapg.tr
  dapg.trav = as.data.frame(AverageExpression(dapg.tr,group.by = "orig.ident")$RNA)
  otpg.trav = as.data.frame(AverageExpression(otpg.tr,group.by = "orig.ident")$RNA)
  rglda.trav = as.data.frame(AverageExpression(rglda.tr,group.by = "orig.ident")$RNA)
  pgda.trav = as.data.frame(AverageExpression(pgda.tr,group.by = "orig.ident")$RNA)
  
  pgtr = cbind(dapg.trav,otpg.trav)
  colnames(pgtr) = c("DA","other")
  pgtr = pgtr[which(rowSums(pgtr) > 0),]
  pgtr$fc = log2((pgtr$DA)/(pgtr$other))
  pgtr$gene = rownames(pgtr)
  pgtr[is.infinite(pgtr$fc),"fc"] = 0
  ggdensity(pgtr$fc)
  
  p10.1 = ggscatter(pgtr,x = "DA",y = "other",size = 1) + 
    ggrepel::geom_text_repel(
      data = pgtr[c("Fabp7","Malat1","Shh"),],
      aes(label = gene,color = "red"),
      size = 3,
      show.legend = FALSE ) + 
    geom_text(x = 10,y = 100,label = paste0("p = 1, cor = ",round(cor(pgtr$DA,pgtr$other),3)),
              size = 5)
  p10.1
  t.test(pgtr$DA, pgtr$other, paired = TRUE, alternative = "two.sided")
  
  datr = cbind(rglda.trav,pgda.trav)
  colnames(datr) = c("rgl1","NProgBM")
  datr = datr[which(rowSums(datr) > 0),]
  cor(datr[1],datr[2])
  datr$fc = log2((datr$rgl1)/(datr$NProgBM))
  datr$gene = rownames(datr)
  datr[is.infinite(datr$fc),"fc"] = 0
  ggdensity(datr$fc)
  t.test(datr$rgl1, datr$NProgBM, paired = TRUE, alternative = "two.sided")
  
  
  p10.2 = ggscatter(datr,x = "rgl1",y = "NProgBM",size = 1) + 
    geom_text(x = 10,y = 115,label = paste0("p = 1, cor = ",round(cor(datr$rgl1,datr$NProgBM),3)),
              size = 5)
  p10.2
  ggsave(p10.1, filename = "final_result/xdraft/da_progbm_v_other_progbm_cor.pdf",width = 5,height = 4)
  ggsave(p10.2, filename = "final_result/xdraft/rgl1_da_v_progbm_da_cor.pdf",width = 5,height = 4)
  
    
  
}



