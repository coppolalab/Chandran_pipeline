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
library(limma)

#Longitudinal analysis
library(betr)

#Plotting
library(Cairo)
library(WGCNA)
library(heatmap.plus)
library(flashClust)
enableWGCNAThreads()

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

lumi.import <- readRDS.gz("../WGCNA/save/lumi.collapse.rda") 

lumi.tgdox <- lumi.import[,lumi.import$Gt.Treatment == "Tg.DOX"]
tgdox.out <- betr(lumi.tgdox, timepoint = factor(lumi.tgdox$TimePoint), replicate = factor(lumi.tgdox$Replicate), twoCondition = FALSE)
saveRDS.gz(tgdox.out, "./save/tgdox.out.rda")

tgdox.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = tgdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(tgdox.genes, "./tgdox.out.xlsx")

lumi.tgnd <- lumi.import[,lumi.import$Gt.Treatment == "Tg.ND"]
tgnd.out <- betr(lumi.tgnd, timepoint = factor(lumi.tgnd$TimePoint), replicate = factor(lumi.tgnd$Replicate), twoCondition = FALSE)
saveRDS.gz(tgnd.out, "./save/tgnd.out.rda")

tgnd.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = tgnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(tgnd.genes, "./tgnd.out.xlsx")

lumi.wtdox <- lumi.import[,lumi.import$Gt.Treatment == "WT.DOX" & lumi.import$Replicate < 4]
wtdox.out <- betr(lumi.wtdox, timepoint = factor(lumi.wtdox$TimePoint), replicate = factor(lumi.wtdox$Replicate), twoCondition = FALSE)
saveRDS.gz(wtdox.out, "./save/wtdox.out.rda")

wtdox.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = wtdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(wtdox.genes, "./wtdox.out.xlsx")

lumi.doxnd <- lumi.import[,grepl("Tg.DOX$|Tg.ND", lumi.import$Gt.Treatment)]
doxnd.out <- betr(lumi.doxnd, cond = factor(lumi.doxnd$Gt.Treatment), timepoint = factor(lumi.doxnd$TimePoint), replicate = factor(lumi.doxnd$Replicate), twoCondition = TRUE)
saveRDS.gz(doxnd.out, "./save/doxnd.out.rda")

doxnd.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = doxnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(doxnd.genes, "./doxnd.out.xlsx")

lumi.tgwt <- lumi.import[,grepl("Tg.DOX$|WT.DOX", lumi.import$Gt.Treatment)]
lumi.tgwt <- lumi.tgwt[,!(lumi.tgwt$Gt.Treatment == "WT.DOX" & lumi.tgwt$Replicate == 4)]
tgwt.out <- betr(lumi.tgwt, cond = factor(lumi.tgwt$Gt.Treatment), timepoint = factor(lumi.tgwt$TimePoint), replicate = factor(lumi.tgwt$Replicate), twoCondition = TRUE)
saveRDS.gz(tgwt.out, "./save/tgwt.out.rda")

tgwt.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = tgwt.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(tgwt.genes, "./tgwt.out.xlsx")

lumi.control <- lumi.import[,grepl("Tg.ND|WT.DOX", lumi.import$Gt.Treatment)]
lumi.control <- lumi.control[,!(lumi.control$Gt.Treatment == "WT.DOX" & lumi.control$Replicate == 4)]
control.out <- betr(lumi.control, cond = factor(lumi.control$Gt.Treatment), timepoint = factor(lumi.control$TimePoint), replicate = factor(lumi.control$Replicate), twoCondition = TRUE)
saveRDS.gz(control.out, "./save/control.out.rda")

control.genes <- data.frame(Symbol = rownames(expr.collapse), Probability = control.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(control.genes, "./control.out.xlsx")

patient.genes <- readRDS.gz("../../FRDA\\ project/betr/save/patient.out.rda")
carrier.genes <- readRDS.gz("../../FRDA\\ project/betr/save/carrier.out.rda")
control.genes <- readRDS.gz("../../FRDA\\ project/betr/save/control.out.rda")
pca.genes <- read.xlsx("../../FRDA project/betr/pca.out.xlsx")
pco.genes <- read.xlsx("../../FRDA project/betr/pco.out.xlsx")
cc.genes <- read.xlsx("../../FRDA project/betr/cc.out.xlsx")

human.genes <- names(pca.genes) %>% toupper
mouse.genes <- names(tgdox.out) %>% toupper  
mouse.overlap <- mouse.genes[mouse.genes %in% human.genes]

cutoff <- 500
doxnd.ranked <- doxnd.genes[1:cutoff,]$Symbol %>% toupper #Special exception
tgwt.ranked <- tgwt.genes[1:cutoff,]$Symbol %>% toupper 
pca.ranked <- pca.genes[1:cutoff,]$Symbol
pco.ranked <- pco.genes[1:cutoff,]$Symbol
doxnd.pca <- doxnd.ranked[doxnd.ranked %in% pca.ranked] %>% sort
doxnd.pco <- doxnd.ranked[doxnd.ranked %in% pco.ranked] %>% sort
tgwt.pca <- tgwt.ranked[tgwt.ranked %in% pca.ranked] %>% sort
tgwt.pco <- tgwt.ranked[tgwt.ranked %in% pco.ranked] %>% sort

doxnd.pca.hyper <- phyper(length(doxnd.pca), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
doxnd.pco.hyper <- phyper(length(doxnd.pco), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
tgwt.pca.hyper <- phyper(length(tgwt.pca), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
tgwt.pco.hyper <- phyper(length(tgwt.pco), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)

hyper.list <- c(doxnd.pca.hyper, doxnd.pco.hyper, tgwt.pca.hyper, tgwt.pco.hyper)

#DRG
drg.import <- readRDS.gz("../dtw/save/lumi.drg.rda")
drg.tgdox <- drg.import[,drg.import$Gt.Treatment == "Tg.DOX"]
drg.tgdox.out <- betr(drg.tgdox, timepoint = factor(drg.tgdox$TimePoint), replicate = factor(drg.tgdox$Replicate), twoCondition = FALSE)
saveRDS.gz(drg.tgdox.out, "./save/drg.tgdox.out.rda")

drg.tgdox.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.tgdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.tgdox.genes, "./drg.tgdox.out.xlsx")

drg.tgnd <- drg.import[,drg.import$Gt.Treatment == "Tg.ND"]
drg.tgnd.out <- betr(drg.tgnd, timepoint = factor(drg.tgnd$TimePoint), replicate = factor(drg.tgnd$Replicate), twoCondition = FALSE)
saveRDS.gz(drg.tgnd.out, "./save/drg.tgnd.out.rda")

drg.tgnd.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.tgnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.tgnd.genes, "./drg.tgnd.out.xlsx")

drg.wtdox <- drg.import[,drg.import$Gt.Treatment == "WT.DOX" & drg.import$Replicate < 4]
drg.wtdox.out <- betr(drg.wtdox, timepoint = factor(drg.wtdox$TimePoint), replicate = factor(drg.wtdox$Replicate), twoCondition = FALSE)
saveRDS.gz(drg.wtdox.out, "./save/drg.wtdox.out.rda")

drg.wtdox.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.wtdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.wtdox.genes, "./drg.wtdox.out.xlsx")

drg.doxnd <- drg.import[,grepl("Tg.DOX$|Tg.ND", drg.import$Gt.Treatment)]
drg.doxnd.out <- betr(drg.doxnd, cond = factor(drg.doxnd$Gt.Treatment), timepoint = factor(drg.doxnd$TimePoint), replicate = factor(drg.doxnd$Replicate), twoCondition = TRUE)
saveRDS.gz(drg.doxnd.out, "./save/drg.doxnd.out.rda")

drg.doxnd.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.doxnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.doxnd.genes, "./drg.doxnd.out.xlsx")

drg.tgwt <- drg.import[,grepl("Tg.DOX$|WT.DOX", drg.import$Gt.Treatment)]
drg.tgwt <- drg.tgwt[,!(drg.tgwt$Gt.Treatment == "WT.DOX" & drg.tgwt$Replicate == 4)]
drg.tgwt.out <- betr(drg.tgwt, cond = factor(drg.tgwt$Gt.Treatment), timepoint = factor(drg.tgwt$TimePoint), replicate = factor(drg.tgwt$Replicate), twoCondition = TRUE)
saveRDS.gz(drg.tgwt.out, "./save/drg.tgwt.out.rda")

drg.tgwt.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.tgwt.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.tgwt.genes, "./drg.tgwt.out.xlsx")

drg.control <- drg.import[,grepl("Tg.ND|WT.DOX", drg.import$Gt.Treatment)]
drg.control <- drg.control[,!(drg.control$Gt.Treatment == "WT.DOX" & drg.control$Replicate == 4)]
drg.control.out <- betr(drg.control, cond = factor(drg.control$Gt.Treatment), timepoint = factor(drg.control$TimePoint), replicate = factor(drg.control$Replicate), twoCondition = TRUE)
saveRDS.gz(drg.control.out, "./save/drg.control.out.rda")

drg.control.genes <- data.frame(Symbol = featureNames(drg.import), Probability = drg.control.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(drg.control.genes, "./drg.control.out.xlsx")

cutoff <- 500
drg.doxnd.ranked <- drg.doxnd.genes[1:cutoff,]$Symbol %>% toupper #Special exception
drg.tgwt.ranked <- drg.tgwt.genes[1:cutoff,]$Symbol %>% toupper 
drg.pca.ranked <- pca.genes[1:cutoff,]$Symbol
drg.pco.ranked <- pco.genes[1:cutoff,]$Symbol
drg.doxnd.pca <- drg.doxnd.ranked[drg.doxnd.ranked %in% drg.pca.ranked] %>% sort
drg.doxnd.pco <- drg.doxnd.ranked[drg.doxnd.ranked %in% drg.pco.ranked] %>% sort
drg.tgwt.pca <- drg.tgwt.ranked[drg.tgwt.ranked %in% drg.pca.ranked] %>% sort
drg.tgwt.pco <- drg.tgwt.ranked[drg.tgwt.ranked %in% drg.pco.ranked] %>% sort

drg.doxnd.pca.hyper <- phyper(length(drg.doxnd.pca), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
drg.doxnd.pco.hyper <- phyper(length(drg.doxnd.pco), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
drg.tgwt.pca.hyper <- phyper(length(drg.tgwt.pca), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
drg.tgwt.pco.hyper <- phyper(length(drg.tgwt.pco), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)

drg.hyper.list <- c(drg.doxnd.pca.hyper, drg.doxnd.pco.hyper, drg.tgwt.pca.hyper, drg.tgwt.pco.hyper)

#Cerebellum
cerebellum.import <- readRDS.gz("../dtw/save/lumi.cerebellum.rda")
cerebellum.tgdox <- cerebellum.import[,cerebellum.import$Gt.Treatment == "Tg.DOX"]
cerebellum.tgdox.out <- betr(cerebellum.tgdox, timepoint = factor(cerebellum.tgdox$TimePoint), replicate = factor(cerebellum.tgdox$Replicate), twoCondition = FALSE)
saveRDS.gz(cerebellum.tgdox.out, "./save/cerebellum.tgdox.out.rda")

cerebellum.tgdox.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.tgdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.tgdox.genes, "./cerebellum.tgdox.out.xlsx")

cerebellum.tgnd <- cerebellum.import[,cerebellum.import$Gt.Treatment == "Tg.ND"]
cerebellum.tgnd.out <- betr(cerebellum.tgnd, timepoint = factor(cerebellum.tgnd$TimePoint), replicate = factor(cerebellum.tgnd$Replicate), twoCondition = FALSE)
saveRDS.gz(cerebellum.tgnd.out, "./save/cerebellum.tgnd.out.rda")

cerebellum.tgnd.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.tgnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.tgnd.genes, "./cerebellum.tgnd.out.xlsx")

#Oh noes
cerebellum.wtdox <- cerebellum.import[,cerebellum.import$Gt.Treatment == "WT.DOX" & cerebellum.import$Replicate < 4]
cerebellum.wtdox.out <- betr(cerebellum.wtdox, timepoint = factor(cerebellum.wtdox$TimePoint), replicate = factor(cerebellum.wtdox$Replicate), twoCondition = FALSE)
saveRDS.gz(cerebellum.wtdox.out, "./save/cerebellum.wtdox.out.rda")

cerebellum.wtdox.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.wtdox.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.wtdox.genes, "./cerebellum.wtdox.out.xlsx")

cerebellum.doxnd <- cerebellum.import[,grepl("Tg.DOX$|Tg.ND", cerebellum.import$Gt.Treatment)]
cerebellum.doxnd.out <- betr(cerebellum.doxnd, cond = factor(cerebellum.doxnd$Gt.Treatment), timepoint = factor(cerebellum.doxnd$TimePoint), replicate = factor(cerebellum.doxnd$Replicate), twoCondition = TRUE)
saveRDS.gz(cerebellum.doxnd.out, "./save/cerebellum.doxnd.out.rda")

cerebellum.doxnd.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.doxnd.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.doxnd.genes, "./cerebellum.doxnd.out.xlsx")

#Oh noes
cerebellum.tgwt <- cerebellum.import[,grepl("Tg.DOX$|WT.DOX", cerebellum.import$Gt.Treatment)]
cerebellum.tgwt <- cerebellum.tgwt[,!(cerebellum.tgwt$Gt.Treatment == "WT.DOX" & cerebellum.tgwt$Replicate == 4)]
cerebellum.tgwt.out <- betr(cerebellum.tgwt, cond = factor(cerebellum.tgwt$Gt.Treatment), timepoint = factor(cerebellum.tgwt$TimePoint), replicate = factor(cerebellum.tgwt$Replicate), twoCondition = TRUE)
saveRDS.gz(cerebellum.tgwt.out, "./save/cerebellum.tgwt.out.rda")

cerebellum.tgwt.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.tgwt.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.tgwt.genes, "./cerebellum.tgwt.out.xlsx")

#Oh noes
cerebellum.control <- cerebellum.import[,grepl("Tg.ND|WT.DOX", cerebellum.import$Gt.Treatment)]
cerebellum.control <- cerebellum.control[,!(cerebellum.control$Gt.Treatment == "WT.DOX" & cerebellum.control$Replicate == 4)]
cerebellum.control.out <- betr(cerebellum.control, cond = factor(cerebellum.control$Gt.Treatment), timepoint = factor(cerebellum.control$TimePoint), replicate = factor(cerebellum.control$Replicate), twoCondition = TRUE)
saveRDS.gz(cerebellum.control.out, "./save/cerebellum.control.out.rda")

cerebellum.control.genes <- data.frame(Symbol = featureNames(cerebellum.import), Probability = cerebellum.control.out) %>% arrange(desc(Probability)) %>% slice(1:500)
write.xlsx(cerebellum.control.genes, "./cerebellum.control.out.xlsx")

cerebellum.doxnd.ranked <- cerebellum.doxnd.genes[1:cutoff,]$Symbol %>% toupper #Special exception
cerebellum.pca.ranked <- pca.genes[1:cutoff,]$Symbol
cerebellum.pco.ranked <- pco.genes[1:cutoff,]$Symbol
cerebellum.doxnd.pca <- cerebellum.doxnd.ranked[cerebellum.doxnd.ranked %in% cerebellum.pca.ranked] %>% sort
cerebellum.doxnd.pco <- cerebellum.doxnd.ranked[cerebellum.doxnd.ranked %in% cerebellum.pco.ranked] %>% sort

cerebellum.doxnd.pca.hyper <- phyper(length(cerebellum.doxnd.pca), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
cerebellum.doxnd.pco.hyper <- phyper(length(cerebellum.doxnd.pco), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)

cerebellum.hyper.list <- c(cerebellum.doxnd.pca.hyper, cerebellum.doxnd.pco.hyper)

drg.heart.doxnd <- drg.doxnd.ranked[drg.doxnd.ranked %in% doxnd.ranked]
cerebellum.heart.doxnd <- cerebellum.doxnd.ranked[cerebellum.doxnd.ranked %in% doxnd.ranked]
cerebellum.drg.doxnd <- cerebellum.doxnd.ranked[cerebellum.doxnd.ranked %in% drg.doxnd.ranked]

drg.heart.tgwt <- drg.tgwt.ranked[drg.tgwt.ranked %in% tgwt.ranked]

drg.heart.doxnd.hyper <- phyper(length(drg.heart.doxnd), cutoff, nrow(lumi.import) - cutoff, cutoff, lower.tail = F)
cerebellum.heart.doxnd.hyper <- phyper(length(cerebellum.heart.doxnd), cutoff, nrow(lumi.import) - cutoff, cutoff, lower.tail = F)
cerebellum.drg.doxnd.hyper <- phyper(length(cerebellum.drg.doxnd), cutoff, nrow(lumi.import) - cutoff, cutoff, lower.tail = F)

drg.heart.tgwt.hyper <- phyper(length(drg.heart.tgwt), cutoff, nrow(lumi.import) - cutoff, cutoff, lower.tail = F)

overlap.list <- c(drg.heart.doxnd.hyper, cerebellum.heart.doxnd.hyper, cerebellum.drg.doxnd.hyper, drg.heart.tgwt.hyper)
