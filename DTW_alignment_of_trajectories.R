#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Fall 2022
#' The following program calculates the DTW distances between the flow cytometry data
#' of Jia et al and the context specific data of Cook and depicts them in a heatmap

library(dtw)
library(plotly)

####Inputs
file.path("/Users/annicenajafi/Desktop/fitting_markovs/mean/timecourse_data_")->mean.val.dir

file.path("/Users/annicenajafi/Desktop/fitting_markovs/cutoffs/timecourse_data_")->input.dir

c(7.33, 8, 10)->MET.range

data.input<-data.inputs[1,]

cutoff.range<-seq(10, 100, 5) 
#Read the flow cytometry data
read.csv("/Users/annicenajafi/Desktop/fitting_markovs/markov_gold.csv")->gold.std
gold.std[1:7,]->gold.std
#Set the working directory
setwd('/Users/annicenajafi/Desktop/context_paper/')
#Read the data table
read.csv("DataTables/DataTable.csv")->data.inputs
#States
states <- c("Epithelial", "Hybrid", "Mesenchymal")

state.no<-length(states)

replicate.no <- 10 #How many runs of the algorithm

TGFBeta.cases<-c(1, 4, 7, 10)

#Loop through the TGFBeta cases
for(i in TGFBeta.cases){
  
  data.inputs[i,]->data.input #Read the input file
  #Make csv files containing the mean value
  populate.dtw.files(mean.val.dir, data.input, input.dir, cutoff.range, replicate.no)
  #Find the DTW distance
  DTW.calculate(mean.val.dir, data.input, state.no, cutoff.range, MET.range)->sum.mats
  #Make heatmap
  plot_ly(z=matmat, type="heatmap",colors = c(color.scheme[1], color.scheme[5], color.scheme[7], color.scheme[8]))%>% layout(title = list(text = paste0(data.input$CellLine," - ", data.input$Factor), y=0.99),
                                                                                                                           xaxis = list(title = '     E              H              M          Total',zeroline = TRUE, showticklabels = FALSE),
                                                                                                                           yaxis = list(title = 'Cutoff', showticklabels = FALSE, nticks=29), font=t)
}


####Functions
#' Title
#'
#' @param mean.val.dirthe prefix to the mean value dir
#' @param data.input the data input info should be like:
#'  CellLine Factor       Type metadata count       Name
#'     MCF7  TGFB1 TimeCourse metadata   UMI TGFB1_MCF7
#' @param input.dir directory to read files from
#' @param cutoff.range #range of cutoffs example: seq(10, 100, 5)
#' @param replicate.no how many times you ran the algorithm
#'
#' @return does not return - saves everything in the mean val dir as csv files
populate.dtw.files <- function(mean.val.dir, data.input, input.dir, cutoff.range, replicate.no){
  
  for(cutoff in cutoff.range){
    
    df <- data.frame(matrix(ncol = 3, nrow = 0))
    
    for(k in 1:replicate.no){
      
      read.csv(paste0(input.dir, data.input$CellLine, "_",
                      data.input$Factor, "_", cutoff, "_", k, ".csv"))->df.bind
      
      rbind(df, df.bind)->df
      
    }
    
    counter<-1
    
    df.hold <- setNames(data.frame(matrix(ncol=3, nrow=0)), c("time", "X", "value"))
    #iterate through states
    for(state in states){
      #Set the counter
      counter+1->counter
      #Find mean for every state
      df[df$variable==state,] %>% 
        
        group_by(time) %>%
        
        summarise(across(where(is.numeric), mean))->df.in
      
      df.in$variable<-state
      #Assign variables
      nam <- paste("df.", counter, sep = "")
      
      assign(nam, df.in)
      
      rbind(df.hold, df.in)->df.hold
    }
    #Save worj as csv file
    write.csv(df.hold, paste0(mean.val.dir, data.input$CellLine, "_",
                              data.input$Factor, "_", cutoff, ".csv"))
  }
}

#' Title
#'
#' @param mean.val.dir prefix to the mean value dir 
#' @param data.input the data input info should be like: 
#' #'  CellLine Factor       Type metadata count       Name
#'     MCF7  TGFB1 TimeCourse metadata   UMI TGFB1_MCF7
#' @param state.no Number of states
#' @param cutoff.range range of cutoffs for highly variable genes 
#' @param MET.range range to exclude example c(3, 4, 5)

DTW.calculate <- function(mean.val.dir, data.input, state.no, cutoff.range, MET.range){
  #initialize matrix
  matmat<-matrix(nrow=length(cutoff.range), ncol=(state.no+1))
  h<-1
  for(cutoff in cutoff.range){
    
    df<- read.csv(paste0(mean.val.dir, data.input$CellLine, "_",
                         data.input$Factor, "_", cutoff, ".csv"))
    #Exclude data from MET range
    df %>% filter(!(time %in% MET.range))->df
    #iterate over states and find DTW
    for(state in unique(df$variable)){
      
      df %>% filter(variable==state)->df.state
      
      dtw(df.state$value, gold.std[state])$distance->dtw.d
      
      matmat[h, which(unique(df$variable)==state)]<-dtw.d
    }
    h+1->h
  }
  rownames(matmat)<-cutoff.range
  #Save in matrix
  sum.mats<-NA
  for(state in state.no){
    
    matmat[,state.no]+sum.mats->sum.mats
    
  }
  #Find total DTW distance
  sum.mats->matmat[,(state.no+1)]
  #matmat[,1]+matmat[,2]+matmat[,3]->matmat[,4]
  return(sum.mats)
}


###A549 --> cutoff 45
###DU145 --> cutoff 40
###MCF7 --> cutoff 50
###OVCA420 --> cutoff 55




