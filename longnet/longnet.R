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
library(abind)
library(lumi)
library(lumiHumanAll.db)
library(annotate)

#Network Analysis
library(minerva)
library(minet)
library(parallel)
library(moments)
library(Rgraphviz)
library(flashClust)
library(WGCNA)
enableWGCNAThreads()

#Functional programming
library(magrittr)
library(purrr)
library(functional)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

library(Cairo)

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

intensities.means <- readRDS.gz("../dtw/save/intensities.means.rda")
intensities.ica.genes <- readRDS.gz("../ica/save/intensities.ica.genes.rda")
dtw.doxnd <- readRDS.gz("../dtw/save/dtw.doxnd.rda")
dtw.tgwt <- readRDS.gz("../dtw/save/dtw.tgwt.rda")
dtw.control <- readRDS.gz("../dtw/save/dtw.control.rda")
fdata <- readRDS.gz("../dtw/save/fdata.rda")

intensities.tgdox <- intensities.means$Tg.DOX %>% data.frame
ica.tgdox <- intensities.ica.genes$Tg.DOX %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) 
dtw.tgdox <- unique(c(dtw.doxnd[1:2000,]$Symbol, dtw.tgwt[1:2000,]$Symbol)) 
all.tgdox <- unique(c(ica.tgdox, dtw.tgdox)) %>% str_replace("\\-", "\\.")
saveRDS.gz(ica.tgdox, "./save/ica.tgdox.rda")

rownames(intensities.tgdox) %<>% str_replace("\\-", "\\.")
mi.expr.df <- slice(intensities.tgdox, match(all.tgdox, rownames(intensities.tgdox))) 
rownames(mi.expr.df) <- all.tgdox

mi.mrnet <- minet(t(mi.expr.df), method = "mrnet", estimator = "mi.shrink", disc = "equalfreq")
#mine.mat <- mine(t(mi.expr.df), alpha = 1, n.cores = 6)
#mine.mrnet <- mrnet(mine.mat$MIC)

TOM.matrix <- TOMsimilarity(mi.mrnet, verbose = 5)
dissimilarity.TOM <- 1 - TOM.matrix

CairoPDF("scalefree", height = 6, width = 6)
scaleFreePlot(TOM.matrix) #Meh
dev.off()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")

CairoPDF(file = "genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 50
#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF("gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(t(mi.expr.df), colors = dynamic.colors, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - cor(ME.genes)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

module.expr.out <- data.frame(Symbol = rownames(mi.expr.df), module.color = dynamic.colors)
module.symbols <- split(module.expr.out, factor(module.expr.out$module.color))
write.xlsx(module.expr.out, "./modules.out.xlsx")

source("../../code/GO/enrichr.R")
stringdb.submit <- function(module.color, module.df)
{
    module.genes <- module.df[[module.color]]
    get.stringdb(module.genes, module.color, species.id = 10090)
}
map(names(module.symbols), stringdb.submit, module.symbols)
turquoise.stringdb <- get.stringdb(module.symbols$turquoise, "turquoise", species.id = 10090, edge.threshold = 700)

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- full.df[[index]]
    dir.create(file.path("./enrichr", index), recursive = TRUE, showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
trap1 <- map(names(module.symbols), enrichr.submit, module.symbols, enrichr.terms, FALSE)
#test <- enrichr.submit("black", module.symbols, enrichr.terms, FALSE)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort
