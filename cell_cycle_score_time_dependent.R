#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' This script is written to cell cycle score time-dependent data and plot circular plots

#Load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(gridExtra)
library(grid)
library(lattice)
library(patchwork)
library(data.table)
library(tidyr)
library(reshape2)
library(todor)
library(tools)
library(berryFunctions)
library(R.matlab)
library(dtw)
library(rmatio)
library(Rmagic)
library(umap)
library(phytools)
library(msigdbr)
library(edgeR)
library(readxl)
library(Seurat)
library(pheatmap)

#Read list of EMT markers
read.csv("/Users/annicenajafi/Downloads/EMTGeneListForAnnice.csv")->>EMT.genes
EMT.genes$name->>EMT.genes
#Read EMT signatures for KS method
EMT.sig <<- data.frame(read_excel("/Users/annicenajafi/Downloads/EM_gene_signature_cellLine_KS.xlsx",
                                  col_names=FALSE))

#Set the working directory
setwd('/Users/annicenajafi/Desktop/context_paper/')
source("cellcyclescoring_functions.R")
#msigdb genes
gene.sets <<- msigdbr(species = "Homo sapiens",
                      category="H") %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gene.sets$gene_symbol->>msig.genes

#Read inputs
read.csv("DataTables/DataTable.csv")->data.inputs
letters <- LETTERS
color.scheme <<- c("#191935", "#1B1B3A", "#2F2246", "#422951", "#693668", "#e3f6f5", "#bae8e8", "#2c698d")
annice.emt.color.scheme <<- c('#CDF0EA', '#F7DBF0', '#BEAEE2') #Mesenchymal, Hybrid, Epithelial
annice.emt.color.scheme.bold <<- c("#24A19C","#D96098",  "#BEAEE2")

read.csv("/Users/annicenajafi/Downloads/g2m_genes.csv")$x->g2m_genes
read.csv("/Users/annicenajafi/Downloads/s_genes.csv")$x->s_genes
append(g2m_genes, s_genes)->cell.cyc.genes


for(q in 0:2){
  data.frame()->phase.df
  matrix(nrow=4,ncol=3)->s.score.mat
  
  for(i in c(1+q, 4+q, 7+q, 10+q)){
    
    data.inputs[i,]->data.input
    data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/", data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_UMI_matrix.csv.gz")))
    data$V1->gene.names #Extract gene names
    meta.data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/",data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_metadata.csv.gz")))
    sort(unique(meta.data$Time))->time.points
    
    meta.data$Time[meta.data$Time=="0d"]<-0
    meta.data$Time[meta.data$Time=="8h"]<-0.33
    meta.data$Time[meta.data$Time=="1d"]<-1
    meta.data$Time[meta.data$Time=="3d"]<-3
    meta.data$Time[meta.data$Time=="7d"]<-7
    meta.data$Time[meta.data$Time=="8h_rm"]<-7.33
    meta.data$Time[meta.data$Time=="1d_rm"]<-8
    meta.data$Time[meta.data$Time=="3d_rm"]<-10
    
    levels(meta.data$Time)<-c(0, 0.33, 1, 3, 7, 7.33, 8, 10)
    
    plt.labels<-c(letters[i*3-2], letters[i*3-1], letters[i*3])
    
    #rbind(phase.df, data.frame(data.input$CellLine, time.points[j], s_phase, g1_phase, g2m_phase))->phase.df
    preprocess_data(data, gene.names)->seurat_phase
    cellcycle.score(seurat_phase, meta,data, g2m_genes, s_genes, time.points)->hold.phase.df
    print(dim(hold.phase.df))
    rbind(phase.df, hold.phase.df)->phase.df
    
  }
  
  
  #pheatmap(s.score.mat, display_numbers = F, color=c("#191935", "#1B1B3A", "#2F2246", "#422951", "#693668", "#e3f6f5", "#bae8e8", "#60B3A2"), cluster_rows=FALSE, cluster_cols=FALSE, main = "G2M Score")
  # Transform data in a tidy format (long format)
  
  colnames(phase.df)<-c("CellLine", "Time", "s_phase", "g1_phase", "g2m_phase")
  levels(phase.df$Time)<-c("0", "0.33", "1", "3", "7", "7.33", "8", "10")
  data <- phase.df %>% gather(key = "observation", value="value", -c(1,2)) 
  data$value<-data$value*100
  
  
  empty_bar <- 2
  nObsType <- nlevels(as.factor(data$observation))
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$CellLine)*nObsType, ncol(data)) )
  colnames(to_add) <- colnames(data)
  to_add$CellLine <- rep(levels(data$CellLine), each=empty_bar*nObsType )
  data <- rbind(data, to_add)
  data <- data %>% arrange(CellLine, Time)
  data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
  
  
  # Get the name and the y position of each label
  label_data <- data %>% group_by(id, Time) %>% summarize(tot=sum(value))
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  
  # prepare a data frame for base lines
  base_data <- data %>% 
    group_by(CellLine) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  
  # Make the plot
  p <- ggplot(data) +      
    
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
    scale_fill_manual(values= c("#513252", "#7A4069", "#CA4E79")) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=0.1, size=0.1 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=0.1, size=0.1 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=0.1, size=0.1 , inherit.aes = FALSE ) +
    #geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    #geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    #ggplot2::annotate("text", x = rep(max(data$id),3), y = c(0, 50, 100), label = c("0", "50", "100") , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
    
    ylim(-150,max(label_data$tot, na.rm=T)) +
    
    theme_minimal() +
    theme(
      #legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size=20)
      #plot.margin = unit(rep(-1,4), "cm") 
    ) +
    #labs(tag = "A", cex=20)+
    coord_polar() +
    
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot, label=Time, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=CellLine), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)+
    #ggtitle(bquote(TGF ~beta))
    ggtitle(bquote(TNF))
  
  if(q==0){
    p <- p + ggtitle(bquote(TGF ~beta))
  }else if(q==1){
    p <- p + ggtitle(bquote(EGF))
  }
  
  
  nam <- paste("plt.", q, sep = "")
  assign(nam, p)
  
}

#Show final plots
plt.0+plt.1+plt.2


#colnames(phase.df)<-c("CellLine", "Time", "s_phase", "g1_phase", "g2m_phase")
#levels(phase.df$Time)<-c("0", "0.33", "1", "3", "7", "7.33", "8", "10")



#ggplot(((phase.df)[order(phase.df$Time),]%>% filter(CellLine=="DU145") %>% filter(Time %in% c(0, 0.33, 1, 3, 7))), aes(x=Time, y=y)) +
#  geom_point() +
#  geom_smooth()



















































