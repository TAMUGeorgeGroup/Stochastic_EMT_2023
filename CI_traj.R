#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' This script is written to find EMT related trajectories with confidence intervals 
#' for every cutoff (0.95)

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

#Determine cell lines here:
CellLines<<- c("A549","DU145", "MCF7", "OVCA420")
annice.emt.color.scheme <<- c('#CDF0EA', '#F7DBF0', '#BEAEE2') #Mesenchymal, Hybrid, Epithelial
annice.emt.color.scheme.bold <<- c("#24A19C","#D96098",  "#BEAEE2")
#For TGF\beta only in this case
#dir.path <- paste0("/Users/annicenajafi/Desktop/cutoffs/timecourse_data_", CellLine,"_TGFB1_", k, ".csv")
setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("time","variable", "value"))->binded

h<-1
for(CellLine in CellLines){
  
  for(cutoff in seq(10, 100, 5)){
    setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("time","variable", "value"))->binded
  for(k in 1:10){
  #read.csv(paste0(dir.path, "timecourse_data_", CellLine,"_", i, ".csv"))->data
    dir.path <- paste0("/Users/annicenajafi/Desktop/cutoffs/timecourse_data_", CellLine,"_TGFB1_",cutoff, "_", k, ".csv")
    read.csv(dir.path)->data
    rbind(binded, data)->binded
  }  
    ggplot(binded, aes(x=time, y=value, group=variable, color=variable, stroke=1.5), fill=c("#BE79DF", annice.emt.color.scheme.bold[2], annice.emt.color.scheme.bold[1]))+ 
    #ggtitle(bquote(.(CellLine)~ - TGF ~beta))+
    ggtitle(cutoff)+
    geom_rect(aes(xmin = -Inf,xmax = 7,ymin = -Inf, ymax = Inf),
              fill="#DAEAF1", 
              alpha = .2)+
    geom_vline(xintercept=7)+
    stat_summary(geom="ribbon", fun.data=mean_cl_normal, width=0.1, conf.int=0.95, fill = c(rep(annice.emt.color.scheme[3], 8), rep(annice.emt.color.scheme[2], 8), rep(annice.emt.color.scheme[1], 8)))+
    stat_summary(geom="line", fun.y=mean, linetype="dashed", fill=c(annice.emt.color.scheme[3], annice.emt.color.scheme[2], annice.emt.color.scheme[1]))+
    stat_summary(geom="point", fun.y=mean, color=c(rep("#BE79DF", 8), rep(annice.emt.color.scheme.bold[2], 8), rep(annice.emt.color.scheme.bold[1], 8)), shape=8, size=1)+
    scale_color_manual(values=(c("#BE79DF", annice.emt.color.scheme.bold[2], annice.emt.color.scheme.bold[1])),
                       labels = c("Epithelial", "Hybrid", "Mesenchymal"))+labs(x="Time (Days)",
                                                                               y="Cell Fraction", tag = LETTERS[h], color="States", shape=8)+
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
    )  + guides(color=guide_legend(override.aes=list(fill=NA)))+scale_x_continuous(limits=c(0, 10))->plt
  
  
    nam <- paste("plt.", h, sep = "")
    assign(nam, plt)
    print(h)
    h<-h+1
}
#(plt.1+plt.2+plt.3+plt.3+plt.4)/(plt.5+plt.6+plt.7+plt.8+plt.9)/(plt.10+plt.11+plt.12+plt.13+plt.14)/(plt.15+plt.16+plt.17+plt.18+plt.19)

  ggplot() +                      # Draw ggplot2 plot with text only
  annotate("text",
           x = 1,
           y = 1,
           size = 8,
           #label = bquote(.(CellLine)~ - TNF)) + 
           label = bquote(.(CellLine)~ - TGF ~beta)) + 
  theme_void()+ theme(panel.background = element_rect(fill = '#DAEAF1', colour = 'black'), text = element_text(size=28))->plt
grid.arrange(plt, plt.1, plt.2, plt.3, plt.4,
             plt.5,plt.6, plt.7, plt.8, plt.9, plt.10,
             plt.11, plt.12, plt.13, plt.14, plt.15,
             plt.16, plt.17, plt.18, plt.19,
             top =bquote(.(CellLine)~ - TGF ~beta), nrow = 4)


}

# ggplot(binded, aes(x=time, y=value, group=variable, color=variable, stroke=1.5), fill=c("#BE79DF", annice.emt.color.scheme.bold[2], annice.emt.color.scheme.bold[1]))+ 
#   ggtitle(bquote(.(CellLine)~ - TGF ~beta))+
#   
#   geom_rect(aes(xmin = -Inf,xmax = 7,ymin = -Inf, ymax = Inf),
#             fill="#DAEAF1", 
#             alpha = .2)+
#   geom_vline(xintercept=7)+
#   stat_summary(geom="ribbon", fun.data=mean_cl_normal, width=0.1, conf.int=0.95, fill = c(rep(annice.emt.color.scheme[3], 8), rep(annice.emt.color.scheme[2], 8), rep(annice.emt.color.scheme[1], 8)))+
#   stat_summary(geom="line", fun.y=mean, linetype="dashed", fill=c(annice.emt.color.scheme[3], annice.emt.color.scheme[2], annice.emt.color.scheme[1]))+
#   stat_summary(geom="point", fun.y=mean, color=c(rep("#BE79DF", 8), rep(annice.emt.color.scheme.bold[2], 8), rep(annice.emt.color.scheme.bold[1], 8)), shape=8, size=1)+
#   scale_color_manual(values=(c("#BE79DF", annice.emt.color.scheme.bold[2], annice.emt.color.scheme.bold[1])),
#                      labels = c("Epithelial", "Hybrid", "Mesenchymal"))+labs(x="Time (Days)",
#                                                                              y="Cell Fraction", tag = "D", color="States", shape=8)+
#   theme(
#     # Remove panel border
#     panel.border=element_blank(),  
#     #plot.border = element_blank(),
#     # Remove panel grid lines
#     panel.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Add axis line
#     axis.line = element_line(colour = "black"),
#     #legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20),
#     axis.text = element_text(size = 15),
#     text = element_text(size=18)
#   )  + guides(color=guide_legend(override.aes=list(fill=NA)))+scale_x_continuous(limits=c(0, 10))->plt.4
