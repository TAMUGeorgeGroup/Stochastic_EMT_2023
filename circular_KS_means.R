#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' This script is written to plot circular plots for the mean KS scores for every cluster

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
library(stringr)

#set working directory
wd<-"/Users/annicenajafi/Desktop/dots"
setwd(wd)

#Iterate over IFs
for(IF in list.files()){
  #Initialize matrix
  setNames(data.frame(matrix(ncol = 6, nrow = 0)),c( "mm", "name", "value", "cellline"))->dots.binded
  #Iterate through files/ there should be files for every cell line in a directory named after cell lines
  for(file in list.files(IF)){
    read.csv(paste0(IF,"/", file))->dots
    keys <- colnames(dots)[!grepl('cutoff....seq.10..100..5.',colnames(dots))]
    dots_dt <- as.data.table(dots)
    setNames(data.frame(matrix(ncol = 4, nrow = 0)),c("cutoff", "mean.E", "mean.H", "mean.M"))->mother.df.cut
    for(cutoff in unique(dots_dt$cutoff....seq.10..100..5.)){
      mean((dots_dt %>% filter(cutoff....seq.10..100..5.==cutoff))$score.E)->mean.E
      mean((dots_dt %>% filter(cutoff....seq.10..100..5.==cutoff))$score.H)->mean.H
      mean((dots_dt %>% filter(cutoff....seq.10..100..5.==cutoff))$score.M)->mean.M
      data.frame(cutoff, mean.E, mean.H, mean.M)->cut.df
      rbind(mother.df.cut, cut.df)->mother.df.cut
    }
    melt(mother.df.cut, "cutoff")->dots_dt
    colnames(dots_dt)<-c("mm", "name", "value")
    str_extract(file, "[^_]+")->cellline
    cbind(dots_dt, cellline)->dots_dt
    rbind(dots.binded, dots_dt)->dots.binded
    
  }
    ############################Make plot
    empty_bar <- 3
    nObsType <- nlevels(as.factor(dots.binded$name))
    to_add <- data.frame( matrix(NA, empty_bar*nlevels(dots.binded$cellline)*nObsType, ncol(dots.binded)) )
    colnames(to_add) <- colnames(dots.binded)
    to_add$xellline <- rep(levels(dots.binded$cellline), each=empty_bar*nObsType )
    dots.binded <- rbind(dots.binded, to_add)
    dots.binded <- dots.binded %>% arrange(cellline, mm)
    dots.binded$id <- rep( seq(1, nrow(dots.binded)/nObsType) , each=nObsType)
    
    
    # Get the name and the y position of each label
    label_data <- dots.binded %>% group_by(id, mm) %>% summarize(tot=sum(value))
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    
    
    # prepare a data frame for base lines
    base_data <- dots.binded %>% 
      group_by(cellline) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    
    
    
    # prepare a data frame for grid (scales)
    grid_data <- base_data
    grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1,]
    
    
    LDP<-c(10, 29,48,67)
    var <- c("A549", "DU145", "MCF7", "OVCA420")
    # Make the plot
    p <- ggplot(dots.binded) +      
      
      # Add the stacked bar
      geom_point(aes(x=as.factor(id), y=value, color=name), stat="identity", alpha=1) +
      #scale_color_manual(values= c("#513252", "#7A4069", "#CA4E79"),labels = c("E", "H", "M")) +
      scale_color_manual(values= c("#BEAEE2", "#D96098","#24A19C"),labels = c("E", "H", "M")) +
      #scale_fill_manual(values= rep(c("#513252", "#7A4069", "#CA4E79"), 76)) +
      geom_hline(aes(linetype = "baseline", yintercept = 0)) +
      #geom_hline(aes(linetype = "baseline", yintercept = 0.7)) +
      geom_vline(xintercept = c(19.5,38.5,57.5,76.5)) +
      ylim(-0.7,1) +
      
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
      #geom_text(data=label_data, aes(x=id, y=tot, label=mm, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
      
      # Add base line information
      #geom_segment(data=base_data, aes(x = start, y = -0.7, xend = end, yend = 0.8), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
      geom_text(data=base_data, aes(x = title, y = 0.98, label=cellline), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)+
      #ggtitle(bquote(TGF ~beta))+
      #ggtitle(bquote(TNF))+
      guides(color=guide_legend("Mean KS Score"))
    
    myplt <- paste0("plt", as.character(which(list.files()==IF)))
    p
    assign(myplt,p)
}


#Attach them together
plt1+ggtitle(bquote(TGF ~beta))+plt2+ggtitle(bquote(EGFF))+plt3+ggtitle(bquote(TNF))












