#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' This script is written to apply the EMT metric to the context specificity data paper
#' EMT trajectories will be compared in a pairwise fashion by finding the DTW distance

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
#Read list of EMT markers from the cancer research paper
read.csv("/Users/annicenajafi/Downloads/EMTGeneListForAnnice.csv")->>EMT.genes
EMT.genes$name->>EMT.genes
#Read EMT signatures for KS method
EMT.sig <<- data.frame(read_excel("/Users/annicenajafi/Downloads/EM_gene_signature_cellLine_KS.xlsx",
                                  col_names=FALSE))

#Set the working directory
setwd('/Users/annicenajafi/Desktop/context_paper/')
#msigdb genes
gene.sets <<- msigdbr(species = "Homo sapiens",
                      category="H") %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gene.sets$gene_symbol->>msig.genes

#Read inputs
read.csv("DataTables/DataTable.csv")->data.inputs
letters <- LETTERS

#Set color scheme
color.scheme <<- c("#191935", "#1B1B3A", "#2F2246", "#422951", "#693668", "#e3f6f5", "#bae8e8", "#2c698d")
annice.emt.color.scheme <<- c('#CDF0EA', '#F7DBF0', '#BEAEE2') #Mesenchymal, Hybrid, Epithelial
annice.emt.color.scheme.bold <<- c("#24A19C","#D96098",  "#BEAEE2")

#Read cell cycle genes
read.csv("/Users/annicenajafi/Downloads/g2m_genes.csv")$x->g2m_genes
read.csv("/Users/annicenajafi/Downloads/s_genes.csv")$x->s_genes
append(g2m_genes, s_genes)->cell.cyc.genes



for(i in 1:length(data.inputs$Type)){
  #Read data input
  data.inputs[i,]->data.input
  data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/", data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_UMI_matrix.csv.gz")))
  data$V1->gene.names #Extract gene names
  meta.data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/",data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_metadata.csv.gz")))
  sort(unique(meta.data$Time))->time.points
  plt.labels<-c(letters[i*3-2], letters[i*3-1], letters[i*3])
  scProcess.me(data.input, data, meta.data, genes.names, plt.labels)
  
}
