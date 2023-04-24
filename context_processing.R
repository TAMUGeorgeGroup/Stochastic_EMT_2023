#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' This script is written to apply the EMT metric to the context specificity data paper 
#' EMT trajectories will be compared in a pairwise fashion by finding the DTW distance



#' MLR scoring function/ not needed in this script... would need additional file
#' for helper functions to run
#'
#' @param count.matrix 
#' @param gene.names 
#'
#' @return
#' @export
#'
#' @examples
MLR.label.me <- function(count.matrix, gene.names){
  Yhats<-NA
  for(i in 1:dim(count.matrix)[2]){
    rownames(count.matrix)->gene.names
    manipulate.data(as.matrix(count.matrix[,i]), gene.names, 'scRNASeq')->data.lst
    dataset<-data.lst[1][[1]]; labels<-data.lst[2][[1]]
    normalize.data(transformation.data, dataset, labels, express.type="probeset", gene.score.path)->norm.res
    score.EMT(dataset, norm.res, labels)->resres
    resres[[2]]->Yhat.norm
    c(Yhats, Yhat.norm)->Yhats
  }
  return(Yhats)
}

#' KS scoring function was taken from Priyanka's Github and modified
#' 
#' @param exp.mat gene expression matrix 
#' @param genes genes 
#' @param top200 receives top 200 highly variable genes for scoring
#'
#' @return KS score
KS.label.me <- function(exp.mat, genes, topgenes){
  common.sig <- intersect(EMT.sig[,1],genes)
  EMT.exp.idx <- match(common.sig, genes)
  exp.mat[,EMT.exp.idx]->EMT.exp
  topgenes[topgenes %in% EMT.sig[EMT.sig$...2=="Mes",][,1]]->Mes.genes
  topgenes[topgenes %in% EMT.sig[EMT.sig$...2=="Epi",][,1]]->Epi.genes
  Epi.idx <- colnames(EMT.exp[colnames(EMT.exp) %in% Epi.genes])
  Mes.idx <- colnames(EMT.exp[colnames(EMT.exp) %in% Mes.genes])
  
  
  data.frame(matrix(0, nrow <- dim(EMT.exp)[1], ncol <- 6))->ks
  for(i in 1:nrow(EMT.exp)){
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]))$statistic->ks$X1[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]))$p.value->ks$X2[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]), alternative="greater")$statistic->ks$X3[i]
    ks.test(as.numeric(EMT.exp[i,Mes.idx]), as.numeric(EMT.exp[i,Epi.idx]), alternative="greater")$p.value->ks$X4[i]
    ks.test(as.numeric(EMT.exp[i,Epi.idx]), as.numeric(EMT.exp[i,Mes.idx]), alternative="greater")$statistic->ks$X5[i]
    ks.test(as.numeric(EMT.exp[i,Epi.idx]), as.numeric(EMT.exp[i,Mes.idx]), alternative="greater")$p.value->ks$X6[i]
  }
  final.score <- matrix(0,nrow = nrow(ks),ncol = 1)
  for(i in 1:dim(ks)[1]){
    print(i)
    if(ks$X4[i]<0.05){
      final.score[i,1] <- -1*ks$X3[i]
    }else if(ks$X6[i]<0.05){
      final.score[i,1] <- ks$X5[i]
    }else{
      if(ks$X5[i] == max(c(ks$X3[i], ks$X5[i]))){
        final.score[i,1] <- max(c(ks$X3[i], ks$X5[i]))
      }else{
        final.score[i,1] <- -1 * max(c(ks$X3[i], ks$X5[i]))
      }
    }
  }
  return(final.score[,1])
}


scProcess.me <- function(data.input, data, meta.data, genes.names, plt.labels, cutoff){
  
  ####FIRST STEP: Initial filtering
  
  rownames(data)<-gene.names #Set the gene names
  data$V1 <- NULL
  t(data)->data #Transpose data
  colnames(data)<- gene.names
  rownames(data)<-as.character(seq(1, length(rownames(data))))#Change rownames to numbers so we can feed to MAGIC
  
  #Filter cells having at least 10 expressed genes
  keep_cols <- colSums(data > 0) > 10 
  data.filtered.not.norm <- data[,keep_cols]
  
  #library size normalize
  data.filtered <- library.size.normalize(data)
  
  #Find highly variable genes
  as.data.frame(library.size.normalize(data.filtered.not.norm))->data.df
  data.df[colnames(data.df) %in% EMT.genes]->data.df
  
  ######################################################
  #data.df->data.hold
  #data.hold$time<-meta.data$Time
  #data.hold %>% filter(!(time %in% c("8h_rm", "1d_rm", "3d_rm", "7d", "8h")))->data.hold
  #Seurat::CreateSeuratObject(t(data.hold[1:(ncol(data.hold)-1)]))->seur
  #Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=100)->seur
  #top100 <- head(VariableFeatures(seur), 100)
  #write.csv(top100, "/Users/annicenajafi/Desktop/Variable_Genes_UpToDay/EGF/MCF7_EGF_day3.csv")
  ######################################################
  
  ####SECOND STEP: Extracting highly variable genes
  
  Seurat::CreateSeuratObject(t(data.df))->seur
  Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=60)->seur
  top.var.genes <- head(VariableFeatures(seur), cutoff)
  ####Only uncomment if interested in seeing variable genes
  #Seurat::VariableFeaturePlot(seur, cols = c(color.scheme[3], color.scheme[8]))->plot1
  #plot2 <- Seurat::LabelPoints(plot = plot1, points = top.var.genes, repel = TRUE)
  # plot2 <- plot2 + 
  #   labs(x="Average Expression", y="Standardized Variance", tag = "A")+
  #   #ggtitle(bquote(.(data.input$CellLine)~ - TGF ~beta))+
  #   ggtitle(bquote(.(data.input$CellLine)~ - TNF))+
  #   theme(
  #     # Remove panel border
  #     panel.border = element_blank(),  
  #     # Remove panel grid lines
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     # Remove panel background
  #     panel.background = element_blank(),
  #     # Add axis line
  #     axis.line = element_line(colour = "black"),
  #     #legend.position = "none",
  #     plot.title = element_text(hjust = 0.5, size=20),
  #     axis.text = element_text(size = 15),
  #     text = element_text(size=18)
  #   ) 
  #Find top 200 highly variable genes
  Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=200)->seur
  
  ####STEP III: IMPUTATION VIA MAGIC
  #Run MAGIC
  data_MAGIC <- magic(data.filtered, genes=top.var.genes,
                      knn=15)
  #Add time to MAGIC data
  as.data.frame(data_MAGIC)->data.magic
  data.magic$time <- meta.data$Time
  #Set color scheme / hardcoded for the context specificity data
  meta.data$color[meta.data$Time=="0d"]<-color.scheme[1]
  meta.data$color[meta.data$Time=="8h"]<-color.scheme[2]
  meta.data$color[meta.data$Time=="1d"]<-color.scheme[3]
  meta.data$color[meta.data$Time=="3d"]<-color.scheme[4]
  meta.data$color[meta.data$Time=="7d"]<-color.scheme[5]
  meta.data$color[meta.data$Time=="8h_rm"]<-color.scheme[6]
  meta.data$color[meta.data$Time=="1d_rm"]<-color.scheme[7]
  meta.data$color[meta.data$Time=="3d_rm"]<-color.scheme[8]
  
  meta.data$Time[meta.data$Time=="0d"]<-0
  meta.data$Time[meta.data$Time=="8h"]<-0.33
  meta.data$Time[meta.data$Time=="1d"]<-1
  meta.data$Time[meta.data$Time=="3d"]<-3
  meta.data$Time[meta.data$Time=="7d"]<-7
  meta.data$Time[meta.data$Time=="8h_rm"]<-7.33
  meta.data$Time[meta.data$Time=="1d_rm"]<-8
  meta.data$Time[meta.data$Time=="3d_rm"]<-10
  

  ####STEP IV: Run UMAP
  umapped <- umap(data.magic[1:(length(colnames(data.magic))-1)])
  meta.data$Time->my.colors
  my.colors[my.colors=="0"]<-color.scheme[1]
  my.colors[my.colors=="0.33"]<-color.scheme[2]
  my.colors[my.colors=="1"]<-color.scheme[3]
  my.colors[my.colors=="3"]<-color.scheme[4]
  my.colors[my.colors=="7"]<-color.scheme[5]
  my.colors[my.colors=="7.33"]<-color.scheme[6]
  my.colors[my.colors=="8"]<-color.scheme[7]
  my.colors[my.colors=="10"]<-color.scheme[8]
  
  #Perform archetypal analysis
  archetypes(data = umapped$layout, k = 3, verbose = TRUE)->a
  xyplot(a, umapped$layout, chull = chull(umapped$layout),my.colors, xlab=list("UMAPI",cex=1.25),ylab=list("UMAPII",cex=1.25), main=list(bquote(.(data.input$CellLine)~ - TGF~beta), cex=2),scales=list(cex=10,x=list(relation="free")), par.settings = list(strip.background=list(col=color.scheme)))
  
  ####Uncomment if doing archetypal analysis/ otherwise skip
  #lapply(as.list(1:dim(umapped$layout)[1]), function(x) umapped$layout[x[1],])->mydf
  
  #archetypes = generate_arc(arc_coord = mydf, mean=0, sd=1)
  #data= generate_data(archetypes$XC, N_examples=1e4, jiiter=0.04, size=0.9)
  #plot_arc(arc_data=archetypes, data=data, which_dimensions = 1:2)+theme_bw()
  
  #####Min spanning tree finding / uncomment to try MST
  #as.data.frame(umapped$layout)->umap.df
  #colnames(umap.df)<-c("umapI", "umapII")
  #max(ComputeMST(umap.df)$distance)
  
  
  ####STEP V: Run kmeans with three centers
  kmeans(umapped$layout, 3)->model
  
  data.frame(as.data.frame(umapped$layout)$V1, as.data.frame(umapped$layout)$V2, model$cluster)->df.model
  df.model %>% filter(model.cluster==1)->df.model.1
  df.model %>% filter(model.cluster==2)->df.model.2
  df.model %>% filter(model.cluster==3)->df.model.3
  
  
  pca_x <- princomp(df.model[1:2])
  x_cluster <- data.frame(pca_x$scores,model$cluster, meta.data$color, meta.data$Time)
  x_cluster[order(as.numeric(x_cluster$meta.data.Time)),]->x_cluster
  
  model$cluster->meta.data$cluster
  data.magic$cluster<-meta.data$cluster
  
  ####STEP VI: Filter 
  top100 <- head(VariableFeatures(seur), 100)
  ####Labeling the clusters with KS test
  data.df$cluster<-data.magic$cluster
  
  #Run MAGIC
  data_MAGIC <- magic(data.filtered, genes=top100,
                      knn=15)
  as.data.frame(data_MAGIC)->data.magic
  data.magic$time <- meta.data$Time
  data.magic$cluster<-meta.data$cluster
  
  (data.magic %>% filter(cluster==1))[1:(dim(data.magic)[2]-2)]->exp.mat
  #data[which(data.magic$cluster==1),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.1.lab
  
  (data.magic %>% filter(cluster==2))[1:(dim(data.magic)[2]-2)]->exp.mat
  #data[which(data.magic$cluster==2),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.2.lab
  
  (data.magic %>% filter(cluster==3))[1:(dim(data.magic)[2]-2)]->exp.mat
  #data[which(data.magic$cluster==3),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.3.lab
  
  x_cluster$ks.scores<-0
  x_cluster[x_cluster$model.cluster==1,]$ks.scores<-clust.1.lab
  x_cluster[x_cluster$model.cluster==2,]$ks.scores<-clust.2.lab
  x_cluster[x_cluster$model.cluster==3,]$ks.scores<-clust.3.lab
  
  clust.labs <- data.frame(type=c(rep("1", length(clust.1.lab)), rep("2", length(clust.2.lab)), rep("3", length(clust.3.lab))),
                         value = c(clust.1.lab, clust.2.lab, clust.3.lab))
  
  
  clust.labs %>%
    mutate(meanViz = mean(value)) %>%
    group_by(type) %>%
    arrange((meanViz))
  
  
  as.matrix(as.numeric(clust.labs$value))->my.mat
  colnames(my.mat)<-"KS.score"
  rownames(my.mat)<-as.character(seq(1, length(my.mat[,1])))
  cbind(my.mat, as.numeric(clust.labs$type))->my.mat
  colnames(my.mat)<-c("KS.score", "Cluster")
  
  t <- list(
    family = "Serif",
    size = 15,
    color = "Rblack")
  
fig1<-  plot_ly(z=my.mat, type="heatmap",colors = c(color.scheme[1], color.scheme[5], color.scheme[7], color.scheme[8]))%>% layout(title = list(text = paste0(data.input$CellLine," - ", data.input$Factor), y=0.99),
                                              xaxis = list(title = 'KS Score                            Cluster',zeroline = TRUE, showticklabels = FALSE),
                                              yaxis = list(title = 'Cell'), font=t)
fig1e<-  plot_ly(z=my.mat, type="heatmap",colors = c(color.scheme[1], color.scheme[5], color.scheme[7], color.scheme[8]))%>% layout(title = list(text = paste0(data.input$CellLine," - ", data.input$Factor), y=0.99),
                                                                                                                               xaxis = list(title = 'KS Score                            Cluster',zeroline = TRUE, showticklabels = FALSE),
                                                                                                                               yaxis = list(title = 'Cell'), font=t)

  heatmap(my.mat)
  
  ggplot(data = clust.labs[order(clust.labs$cell),], mapping = aes(x = type,
                                         y = factor(cell),
                                         fill = value)) +
    geom_tile() +
    xlab(label = "Sample")
  
  clust.labs |>
    arrange(type) |>
    mutate(dummy = 1) |>
    ggplot(aes(x = factor(dummy),
               y = factor(value),
               fill = value)) +
    #Add rectangle
    geom_rect(aes(xmin=factor(-0.5),
                  xmax=factor(0.5),
                  ymin=0.5,
                  ymax=1.5),
              colour = "black",
              fill = "transparent") +
    geom_tile() +
    # Change breaks for x axis
    #scale_x_discrete(breaks = c(0,1)) +
    xlab(label = "Sample")
  
 
  # mean.1<- mean((clust.labs %>% filter(type=="1"))$value)
  # mean.2<- mean((clust.labs %>% filter(type=="2"))$value)
  # mean.3<- mean((clust.labs %>% filter(type=="3"))$value)
  # 
  # min.dist <- min(abs(mean.1 - mean.2), abs(mean.2 - mean.3), abs(mean.1 - mean.3))
  # 
  # for(i in 1:10){
  #   find.traj(data.df, data.magic,meta.data)->new.res.lst
  #   new.res.lst[[1]]->clust.labs.comp; new.res.lst[[2]]->x_cluster.comp
  #   mean.1<- mean((clust.labs.comp %>% filter(type=="1"))$value)
  #   mean.2<- mean((clust.labs.comp %>% filter(type=="2"))$value)
  #   mean.3<- mean((clust.labs.comp %>% filter(type=="3"))$value)
  #   
  #   min.dist.comp <- min(abs(mean.1 - mean.2), abs(mean.2 - mean.3), abs(mean.1 - mean.3))
  #   
  #   if(min.dist.comp> min.dist){
  #     clust.labs<-clust.labs.comp
  #     x_cluster<-x_cluster.comp
  #   }
  # }
  #uncomment only if interested in heatmaps for ks scores
  fig_1 <- plot_ly(x=x_cluster$Comp.1, y=x_cluster$Comp.2, z=x_cluster$ks.scores, type="scatter3d", mode="markers",
          color=~ x_cluster$meta.data.Time, colors=color.scheme)%>% layout(title = bquote(.(data.input$CellLine)~ - TGF ~beta),
                                                                           xaxis = list(title = 'UMAP I',zeroline = TRUE),yaxis = list(title = 'UMAP II'))
  
  fig_2 <- plot_ly(x=x_cluster$meta.data.Time, y=x_cluster$ks.scores, colorscale = list(c(0, 0.5, 1), c(color.scheme[1], color.scheme[6], color.scheme[8]))) %>%
    add_trace(type='histogram2dcontour')%>% layout(title = bquote(.(data.input$CellLine)~ - TGF ~beta),
                                                   xaxis = list(title = '',zeroline = TRUE),yaxis = list(title = ''))
    subplot(fig_1, fig_2)
    subplot(
    plot_ly(x = as.numeric(x_cluster$meta.data.Time), color = I("black"), type = 'histogram'), 
    plotly_empty(), 
    plot_ly(x = as.numeric(x_cluster$meta.data.Time),y = x_cluster$ks.scores, type = 'histogram2dcontour', showscale = F), 
    plot_ly(y = x_cluster$ks.scores, color = I("black"), type = 'histogram'),
    nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), 
    shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
  )
    
    data_MAGIC.mlr <- magic(data.filtered, genes=c(cutoff, "CLDN7", "VIM", "CDH1"),
                        knn=15)
    as.data.frame(data_MAGIC.mlr)->data.magic.mlr
    data.magic.mlr$time <- meta.data$Time
    data.magic.mlr$cluster<-meta.data$cluster
    fig1.mlr <- plot_ly(data.magic.mlr, x=~CLDN7,y=~VIM, z=~CDH1, color=~as.factor(cluster), colors=c("#513252", "#7A4069", "#CA4E79"))%>% 
      add_markers()%>% layout(title = list(text = paste0(data.input$CellLine," - ", data.input$Factor), y=0.99))
  #clust.labs %>%
  #  mutate(meanViz = mean(value)) %>%
  #  group_by(type) %>%
  #  arrange((meanViz))
  
  colnames(clust.labs)<-c("cluster", "score")
  
  plt.1<- ggplot(clust.labs, aes(x=score, fill=as.factor(cluster))) + 
    geom_density()+ 
    geom_vline(aes(xintercept=mean((clust.labs %>% filter(cluster=="1"))$score)),
               color="#513252", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=mean((clust.labs %>% filter(cluster=="2"))$score)),
               color="#7A4069", linetype="dashed", size=1)+
    geom_vline(aes(xintercept=mean((clust.labs %>% filter(cluster=="3"))$score)),
               color="#CA4E79", linetype="dashed", size=1)+
    facet_grid(cluster ~ .)+
    ggtitle("KS Score Distribution")+
    labs(x="Value", y="Density", tag = plt.labels[1])+
    theme(
      # Remove panel border
      panel.border = element_blank(),  
      # Remove panel grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "black"),
      #legend.position = "none",
      plot.title = element_text(hjust = 0.5, size=20),
      axis.text = element_text(size = 15),
      text = element_text(size=18)
    ) +
    scale_fill_manual(values= c("#513252", "#7A4069", "#CA4E79"))+
    labs(color = "cluster")+
    guides(fill=guide_legend("Clusters"))
    
 
  
  #clust.labs <- data.frame(1:3, c(clust.1.lab, clust.2.lab, clust.3.lab))
  
  clust.labs[order(clust.labs$score),]->clust.labs
  
  #x_cluster %>% filter(x_cluster$model.cluster=="cluster 1")->x1
  #x_cluster %>% filter(x_cluster$model.cluster=="cluster 2")->x2
  #x_cluster %>% filter(x_cluster$model.cluster=="cluster 3")->x3
  
  #x1$clust<-unique(clust.labs$cluster)[1]
  #x2$clust<-clust.labs$cluster[2]
  #x3$clust<-clust.labs$cluster[3]
  
  #x_cluster<-NULL
  #rbind(x1, x2, x3)->x_cluster
  #x_cluster$model.cluster<-as.character(x_cluster$clust)
  x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[1],]$model.cluster<-"Epithelial"
  x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[2],]$model.cluster<-"Hybrid"
  x_cluster[x_cluster$model.cluster==unique(clust.labs$cluster)[3],]$model.cluster<-"Mesenchymal"
  
  x_cluster$meta.data.Time <- factor(x_cluster$meta.data.Time,                 # Relevel group factor
                                     levels = c("0", "0.33", "1", "3", "7", "7.33", "8", "10"))
  
  
  plt.2 <- ggplot(x_cluster, aes(x = Comp.1, y = Comp.2, shape = as.factor(model.cluster), color = meta.data.Time), fill=meta.data$Time) + 
    
    stat_ellipse(aes(x=Comp.1, y=Comp.2, group=model.cluster, color=unique(model.cluster)),type = "norm",geom = "polygon",alpha = 0.9,color='black', fill=NA, segments = 10)+
    #geom_mark_ellipse(aes(x=Comp.1, y=Comp.2, group=model.cluster, color=model.cluster, fill = model.cluster, label = model.cluster)) +
    geom_point() +
    geom_point(aes(x=mean((x_cluster %>% filter(model.cluster=="Mesenchymal"))$Comp.1), y=mean((x_cluster %>% filter(model.cluster=="Mesenchymal"))$Comp.2)), colour=annice.emt.color.scheme.bold[1], shape=10, size=5, stroke = 2)+
    
    geom_point(aes(x=mean((x_cluster %>% filter(model.cluster=="Hybrid"))$Comp.1), y=mean((x_cluster %>% filter(model.cluster=="Hybrid"))$Comp.2)), colour=annice.emt.color.scheme.bold[2], shape=10, size=5, stroke = 2)+
    
    geom_point(aes(x=mean((x_cluster %>% filter(model.cluster=="Epithelial"))$Comp.1), y=mean((x_cluster %>% filter(model.cluster=="Epithelial"))$Comp.2)), colour=annice.emt.color.scheme.bold[3], shape=10, size=5, stroke = 2)+
    labs(x="UMAP I", y="UMAP II", tag = plt.labels[2], color="Time (Days)")+
    #ggtitle(bquote(.(data.input$CellLine)~ - TGF ~beta))+
    ggtitle(bquote(.(data.input$CellLine)~ - TNF))+
    guides(fill="none", shape="none", color = guide_legend(
      override.aes=list(shape = 15, size=8))) +
    scale_color_manual(values=c(color.scheme, "#D96098", "#24A19C", "#BEAEE2"),
                       labels = c("0",
                                  "0.33",
                                  "1",
                                  "3",
                                  "7",
                                  "7.33",
                                  "8",
                                  "10", "Hybrid", "Mesenchymal", "Epithelial"))+
    theme(
      # Remove panel border
      panel.border = element_blank(),  
      # Remove panel grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "black"),
      #legend.position = "none",
      plot.title = element_text(hjust = 0.5, size=20),
      axis.text = element_text(size = 15),
      text = element_text(size=18)
    ) 
  
  
  
  #+
  # stat_ellipse(type = "t",geom = "polygon",alpha = 0.4)
  clusterM.day <- NA
  clusterH.day <- NA
  clusterE.day <- NA
  for(i in 1:length(unique(x_cluster$meta.data.Time))){
    timepoint <- unique(x_cluster$meta.data.Time)[i]
    x_cluster %>% filter(meta.data.Time==timepoint)->data.time
    nrow(data.time %>% filter(model.cluster=="Mesenchymal"))/nrow(data.time)->clust.M.perc
    nrow(data.time %>% filter(model.cluster=="Hybrid"))/nrow(data.time)->clust.H.perc
    nrow(data.time %>% filter(model.cluster=="Epithelial"))/nrow(data.time)->clust.E.perc
    clusterM.day <- c(clusterM.day, clust.M.perc)
    clusterH.day <- c(clusterH.day, clust.H.perc)
    clusterE.day <- c(clusterE.day, clust.E.perc)
  }
  clusterM.day[2:length(clusterM.day)]->clusterM.day
  clusterH.day[2:length(clusterH.day)]->clusterH.day
  clusterE.day[2:length(clusterE.day)]->clusterE.day
  data.frame(clusterM.day, clusterH.day, clusterE.day, as.numeric(levels(unique(x_cluster$meta.data.Time))))->df.cluster
  colnames(df.cluster)<-c("Mesenchymal", "Hybrid", "Epithelial", "time")
  melt(df.cluster, "time")->melted.clust
  write.csv(melted.clust, "/Users/annicenajafi/Downloads/timecourse_data.csv")
  
  plt.3 <- ggplot(melted.clust, aes(time, value, color = variable)) +
    #geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = -Inf, ymax = Inf),
    #          fill="black", 
     #         alpha = .2)+
    #geom_rect(aes(xmin = 7,xmax = Inf,ymin = -Inf, ymax = Inf),
     #         fill="#FFFFFF", 
      #        alpha = .2)+
    #"#e3f6f5"
    #"#F0EBE3"
    geom_rect(aes(xmin = -Inf,xmax = 7,ymin = -Inf, ymax = Inf),
              fill="#e3f6f5", 
              alpha = .2)+
    geom_point(shape=8, size=1, stroke = 1.5) + geom_line(size=0.8)+
    geom_vline(xintercept=7)+
    ggtitle("EMT Trajectory")+
    
    scale_color_manual(values=annice.emt.color.scheme.bold,
                       labels = c("Mesenchymal", "Hybrid", "Epithelial"))+labs(x="Time (Days)", y="Cell Fraction", tag = plt.labels[3], color="Type", shape=8)+
    theme(
      # Remove panel border
      panel.border=element_blank(),  
      #plot.border = element_blank(),
      # Remove panel grid lines
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "black"),
      #legend.position = "none",
      plot.title = element_text(hjust = 0.5, size=20),
      axis.text = element_text(size = 15),
      text = element_text(size=18)
    ) 
  #+  scale_x_continuous(limits =c(0.0,10.0))
  return(plt.1+plt.2+plt.3)
  
}






find.traj <- function(data.df, data.magic,  meta.data){
  #Run UMAP
  umapped.new <- umap(data.magic[1:(length(colnames(data.magic))-2)])
  
  #Run kmeans with three centers
  kmeans(umapped.new$layout, 3)->model.new
  
  data.frame(as.data.frame(umapped.new$layout)$V1, as.data.frame(umapped.new$layout)$V2, model.new$cluster)->df.model.new
  df.model.new %>% filter(model.new.cluster==1)->df.model.1.new
  df.model.new %>% filter(model.new.cluster==2)->df.model.2.new
  df.model.new %>% filter(model.new.cluster==3)->df.model.3.new
  
  
  
  
  
  pca_x.new <- princomp(df.model.new[1:2])
  x_cluster.new <- data.frame(pca_x.new$scores,model$cluster, meta.data$color, meta.data$Time)
  x_cluster.new[order(as.numeric(x_cluster.new$meta.data.Time)),]->x_cluster.new
  
  model.new$cluster->meta.data$cluster
  data.magic$cluster<-meta.data$cluster
  
  ####Labeling the clusters with KS test
  data.df$cluster<-data.magic$cluster
  (data.df %>% filter(cluster==1))[1:(dim(data.df)[2]-1)]->exp.mat
  #data[which(data.magic$cluster==1),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.1.lab
  
  (data.df %>% filter(cluster==2))[1:(dim(data.df)[2]-1)]->exp.mat
  #data[which(data.magic$cluster==2),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.2.lab
  
  (data.df %>% filter(cluster==3))[1:(dim(data.df)[2]-1)]->exp.mat
  #data[which(data.magic$cluster==3),]->exp.mat
  colnames(exp.mat)->genes
  KS.label.me(exp.mat, genes, top100)->clust.3.lab
  
  clust.labs <- data.frame(type=c(rep("1", length(clust.1.lab)), rep("2", length(clust.2.lab)), rep("3", length(clust.3.lab))),
                           value = c(clust.1.lab, clust.2.lab, clust.3.lab))
  
  clust.labs %>%
    mutate(meanViz = mean(value)) %>%
    group_by(type) %>%
    arrange((meanViz))
  return(list(clust.labs, x_cluster))
}












