#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' contains cell cycle scoring functions


####Cell cycle scoring functions

#' preprocess_data
#'
#' @param data raw count matrix 
#' @param gene.names names of genes
#'
#' @return processed seurat file
preprocess_data<-function(data, gene.names){
  data$V1 <- NULL
  t(data)->data #Transpose data
  colnames(data)<- gene.names
  rownames(data)<-as.character(seq(1, length(rownames(data)))) 
  
  #Filter cells having at least 10 expressed genes
  keep_cols <- colSums(data > 0) > 10 
  data.filtered.not.norm <- data[,keep_cols]
  
  #Find highly variable genes
  as.data.frame(data.filtered.not.norm)->data.df
  
  Seurat::CreateSeuratObject(t(data.df))->seur
  NormalizeData(seur)->seur
  seurat_phase <- CellCycleScoring(seur, 
                                   g2m.features = g2m_genes, 
                                   s.features = s_genes)
  
  return(seurat_phase)
}


#' cellcycle.score
#'
#' @param seurat_phase preprocessed data
#' @param meta metadata
#' @param data raw count matrix
#' @param g2m_genes list of g2m genes
#' @param s_genes list of s phase genes
#' @param time.points timepoints
#'
#' @return dataframe with cell cycle fractions for every timepoint
cellcycle.score <- function(seurat_phase, meta,data, g2m_genes, s_genes, time.points){
  
  data.frame()->phase.df #instantiate dataframe
  data.frame(seurat_phase@meta.data$Phase, meta.data$Time)->cell.cyc
  unique(meta.data$Time)->time.points
  #Loop through df and 
  for(j in 1:length(time.points)){
    
    s_phase <- nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]) %>%
                      filter(seurat_phase.meta.data.Phase=="S"))/nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]))
    
    g1_phase <- nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]) %>%
                       filter(seurat_phase.meta.data.Phase=="G1"))/nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]))
    
    g2m_phase <- nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]) %>%
                        filter(seurat_phase.meta.data.Phase=="G2M"))/nrow(cell.cyc %>% filter(meta.data.Time==time.points[j]))
    
    rbind(phase.df, data.frame(data.input$CellLine, time.points[j], s_phase, g1_phase, g2m_phase))->phase.df
    
  }
  return(phase.df)
}
