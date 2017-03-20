#String operations
library(stringr)

#Reading and writing tables
library(readr)
library(openxlsx)

#For plotting
library(ggplot2)

#For microarray stuff
library(Biobase)
library(matrixStats)
library(lumi)

#Longitudinal analysis
library(ebdbNet)
library(BoolNet)

#Plotting
library(Cairo)
library(WGCNA)
library(heatmap.plus)
enableWGCNAThreads()
library(igraph)
library(TeachingDemos)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

lumi.collapse <- readRDS.gz("../WGCNA/save/lumi.collapse.rda") 
expr.tgdox <- exprs(lumi.collapse[,lumi.collapse$Gt.Treatment == "Tg.DOX"])
saveRDS.gz(expr.tgdox, "./save/expr.tgdox.rda")


