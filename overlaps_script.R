#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2022
#' The following script plots a circos plot showing the intersecting highly variable genes
#' between three EMT induction factors and four cell lines. 


matrix(ncol=12, nrow=12)->mat

for (i in 1:12){
  for(j in 1:12){
    
    #row.val <- eval(parse(text = paste0("top200.", i)))[1:100]
    #col.val <- eval(parse(text = paste0("top200.", j)))[1:100]
    #row.val[which(row.val %in% EMT.sig$...1)]->row.val
    #col.val[which(row.val %in% EMT.sig$...1)]->col.val
    data.inputs[i,]->data.input
    data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/", data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_UMI_matrix.csv.gz")))
    #data$V1->gene.names #Extract gene names
    meta.data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/",data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_metadata.csv.gz")))
    #rownames(data)->gene.names #Set the gene names
    data$V1 -> gene.names
    rownames(data)<-gene.names
    data$V1 <- NULL
    t(data)->data #Transpose data
    colnames(data)<- gene.names
    rownames(data)<-as.character(seq(1, length(rownames(data))))#Change rownames to numbers so we can feed to MAGIC
    
    #Filter cells having at least 10 expressed genes
    keep_cols <- colSums(data > 0) > 10 
    data.filtered.not.norm <- data[,keep_cols]
    
    #library size normalize
    data.filtered <- library.size.normalize(data)
    
    #Run MAGIC
    #data_MAGIC <- magic(data.filtered, genes=EMT.genes,
    #                   knn=15)
    
    #Find highly variable genes
    as.data.frame(library.size.normalize(data.filtered.not.norm))->data.df
    data.df[colnames(data.df) %in% EMT.genes]->data.df
    Seurat::CreateSeuratObject(t(data.df))->seur
    Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=100)->seur
    top100.i <- head(VariableFeatures(seur), 100)
    
    
    data.inputs[j,]->data.input
    data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/", data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_UMI_matrix.csv.gz")))
    #data$V1->gene.names #Extract gene names
    meta.data <- fread(file.path(paste0("/Users/annicenajafi/Desktop/context_paper/TimeCourse/transcription_factors/",data.input$Factor,"/",data.input$CellLine,  "/GSE147405_",data.input$CellLine,"_",data.input$Factor, "_TimeCourse_metadata.csv.gz")))
    #rownames(data)->gene.names #Set the gene names
    data$V1 -> gene.names
    rownames(data)<-gene.names
    data$V1 <- NULL
    t(data)->data #Transpose data
    colnames(data)<- gene.names
    rownames(data)<-as.character(seq(1, length(rownames(data))))#Change rownames to numbers so we can feed to MAGIC
    
    #Filter cells having at least 10 expressed genes
    keep_cols <- colSums(data > 0) > 10 
    data.filtered.not.norm <- data[,keep_cols]
    
    #library size normalize
    data.filtered <- library.size.normalize(data)
    
    #Run MAGIC
    #data_MAGIC <- magic(data.filtered, genes=EMT.genes,
    #                   knn=15)
    
    #Find highly variable genes
    as.data.frame(library.size.normalize(data.filtered.not.norm))->data.df
    data.df[colnames(data.df) %in% EMT.genes]->data.df
    Seurat::CreateSeuratObject(t(data.df))->seur
    Seurat::FindVariableFeatures(seur, selection.method='vst', nfeatures=100)->seur
    top100.j <- head(VariableFeatures(seur), 100)
    
    #length(intersect(row.val, col.val))->mat.val
    length(intersect(top100.i, top100.j))->mat.val
    mat[i,j]<-mat.val
  }
}
diag(mat)<-NA
colnames(mat)<-c("A549_TGFBeta", "A549_EGF","A549_TNF",
                 "DU145_TGFBeta", "DU145_EGF","DU145_TNF",
                 "MCF7_TGFBeta", "MCF7_EGF","MCF7_TNF",
                 "OVCA420_TGFBeta", "OVCA420_EGF","OVCA420_TNF")
rownames(mat)<-c("A549_TGFBeta", "A549_EGF","A549_TNF",
                 "DU145_TGFBeta", "DU145_EGF","DU145_TNF",
                 "MCF7_TGFBeta", "MCF7_EGF","MCF7_TNF",
                 "OVCA420_TGFBeta", "OVCA420_EGF","OVCA420_TNF")
chordDiagram(mat, grid.col = c(blues9[3], blues9[4], blues9[5], 
                               "#FFC4C4", "#EE6983", "#850E35", 
                               "#D2D79F", "#90B77D", "#42855B", 
                               "#423F3E", "#2B2B2B", "#171010"))

chordDiagram(mat, grid.col = c(blues9[3], blues9[4], blues9[5], 
                                                               "#FFC4C4", "#EE6983", "#850E35", 
                                                               "#D2D79F", "#90B77D", "#42855B", 
                                                               "#423F3E", "#2B2B2B", "#171010"), annotationTrack = c("name", "grid"), link.sort = TRUE)

chordDiagram(mat, grid.col = c(blues9[3], blues9[4], blues9[5], 
                               "#FFC4C4", "#EE6983", "#850E35", 
                               "#D2D79F", "#90B77D", "#42855B", 
                               "#423F3E", "#2B2B2B", "#171010"), annotationTrack = c("name", "grid"), link.sort = TRUE, directional = 1, diffHeight = mm_h(3))
title("Intersection of Top 100 Variable Genes", cex = 15)

chordDiagram(data.mat, grid.col = c(blues9[3], blues9[4], blues9[5], blues9[6],
                               "#FFC4C4", "#EE6983", "#850E35", "#400719",
                               "#D2D79F", "#90B77D", "#42855B", 
                               "#423F3E", 
                               #"#2B2B2B", "#171010"
                               ), annotationTrack = c("name", "grid"), link.sort = TRUE, directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow")






chordDiagram(data.mat, grid.col = c(blues9[3], blues9[4], blues9[5], blues9[6],
                                    "#FFC4C4", "#EE6983", "#850E35", "#400719",
                                    "#D2D79F", "#90B77D", "#42855B", 
                                    "#295238"), annotationTrack = c("name", "grid"), link.sort = TRUE, directional = 1, direction.type = c("diffHeight", "arrows"),
link.arr.type = "big.arrow")
