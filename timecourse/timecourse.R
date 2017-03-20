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
library(annotate)
library(sva)
library(peer)
library(limma)

#Longitudinal analysis
library(timecourse)

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
tgdox.reps <- rep(ncol(lumi.tgdox)/5, nrow(lumi.tgdox))
tgdox.repgrp <- rep(1:4, (ncol(lumi.tgdox)/4))
tgdox.out <- mb.long(lumi.tgdox, method = "1D", type = "robust", times = 5, reps = tgdox.reps, rep.grp = tgdox.repgrp, time.grp = factor(lumi.tgdox$TimePoint))
saveRDS.gz(tgdox.out, "./save/tgdox.out.rda")

tgdox.genes <- data.frame(Symbol = rownames(lumi.tgdox), Hotelling = tgdox.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(tgdox.genes, "./tgdox.out.xlsx")

lumi.tgnd <- lumi.import[,lumi.import$Gt.Treatment == "Tg.ND"]
tgnd.reps <- rep(ncol(lumi.tgnd)/5, nrow(lumi.tgnd))
tgnd.repgrp <- rep(1:4, (ncol(lumi.tgnd)/4))
tgnd.out <- mb.long(lumi.tgnd, method = "1D", type = "robust", times = 5, reps = tgnd.reps, rep.grp = tgnd.repgrp, time.grp = factor(lumi.tgnd$TimePoint))
saveRDS.gz(tgnd.out, "./save/tgnd.out.rda")

tgnd.genes <- data.frame(Symbol = rownames(lumi.tgnd), Hotelling = tgnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(tgnd.genes, "./tgnd.out.xlsx")

lumi.wtdox <- lumi.import[,lumi.import$Gt.Treatment == "WT.DOX" & lumi.import$Replicate < 4]
wtdox.reps <- rep(ncol(lumi.wtdox)/5, nrow(lumi.wtdox))
wtdox.repgrp <- rep(1:3, (ncol(lumi.wtdox)/3))
wtdox.out <- mb.long(lumi.wtdox, method = "1D", type = "robust", times = 5, reps = wtdox.reps, rep.grp = wtdox.repgrp, time.grp = factor(lumi.wtdox$TimePoint))
saveRDS.gz(wtdox.out, "./save/wtdox.out.rda")

wtdox.genes <- data.frame(Symbol = rownames(lumi.wtdox), Hotelling = wtdox.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(wtdox.genes, "./wtdox.out.xlsx")

lumi.doxnd <- lumi.import[,grepl("Tg.DOX$|Tg.ND", lumi.import$Gt.Treatment)]
doxnd.reps <- cbind(rep(length(which(lumi.doxnd$Gt.Treatment == "Tg.DOX"))/5, nrow(lumi.doxnd)), rep(length(which(lumi.doxnd$Gt.Treatment == "Tg.ND"))/5, nrow(lumi.doxnd)))
doxnd.out <- mb.long(lumi.doxnd, method = "2D", type = "robust", times = 5, reps = doxnd.reps, rep.grp = lumi.doxnd$Replicate, time.grp = factor(lumi.doxnd$TimePoint), condition.grp = factor(lumi.doxnd$Gt.Treatment))
saveRDS.gz(doxnd.out, "./save/doxnd.out.rda")

doxnd.genes <- data.frame(Symbol = rownames(lumi.doxnd), Hotelling = doxnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(doxnd.genes, "./doxnd.out.xlsx")

lumi.tgwt <- lumi.import[,grepl("Tg.DOX$|WT.DOX", lumi.import$Gt.Treatment)]
lumi.tgwt <- lumi.tgwt[,!(lumi.tgwt$Gt.Treatment == "WT.DOX" & lumi.tgwt$Replicate == 4)]
tgwt.reps <- cbind(rep(length(which(lumi.tgwt$Gt.Treatment == "Tg.DOX"))/5, nrow(lumi.tgwt)), rep(length(which(lumi.tgwt$Gt.Treatment == "WT.DOX"))/5, nrow(lumi.tgwt)))
tgwt.out <- mb.long(lumi.tgwt, method = "2D", type = "robust", times = 5, reps = tgwt.reps, rep.grp = lumi.tgwt$Replicate, time.grp = factor(lumi.tgwt$TimePoint), condition.grp = factor(lumi.tgwt$Gt.Treatment))
saveRDS.gz(tgwt.out, "./save/tgwt.out.rda")

tgwt.genes <- data.frame(Symbol = rownames(lumi.tgwt), Hotelling = tgwt.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(tgwt.genes, "./tgwt.out.xlsx")

lumi.control <- lumi.import[,grepl("WT.DOX|Tg.ND", lumi.import$Gt.Treatment)]
lumi.control <- lumi.control[,!(lumi.control$Gt.Treatment == "WT.DOX" & lumi.control$Replicate == 4)]
control.reps <- cbind(rep(length(which(lumi.control$Gt.Treatment == "Tg.ND"))/5, nrow(lumi.control)), rep(length(which(lumi.control$Gt.Treatment == "WT.DOX"))/5, nrow(lumi.control)))
control.out <- mb.long(lumi.control, method = "2D", type = "robust", times = 5, reps = control.reps, rep.grp = lumi.control$Replicate, time.grp = factor(lumi.control$TimePoint), condition.grp = factor(lumi.control$Gt.Treatment))
saveRDS.gz(control.out, "./save/control.out.rda")

control.genes <- data.frame(Symbol = rownames(lumi.control), Hotelling = control.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(control.genes, "./control.out.xlsx")

patient.genes <- readRDS.gz("../../FRDA\\ project/betr/save/patient.out.rda")
#carrier.genes <- readRDS.gz("../../FRDA\\ project/betr/save/carrier.out.rda")
#control.genes <- readRDS.gz("../../FRDA\\ project/betr/save/control.out.rda")
pca.genes <- read.xlsx("../../FRDA project/timecourse/pca.out.xlsx")
pco.genes <- read.xlsx("../../FRDA project/timecourse/pco.out.xlsx")
cc.genes <- read.xlsx("../../FRDA project/timecourse/cc.out.xlsx")

human.genes <- names(patient.genes) %>% toupper
mouse.genes <- rownames(lumi.tgdox) %>% toupper  
mouse.overlap <- mouse.genes[mouse.genes %in% human.genes]

cutoff <- 200
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
drg.tgdox.reps <- rep(ncol(drg.tgdox)/5, nrow(drg.tgdox))
drg.tgdox.repgrp <- rep(1:4, (ncol(drg.tgdox)/4))
drg.tgdox.out <- mb.long(drg.tgdox, method = "1D", type = "robust", times = 5, reps = drg.tgdox.reps, rep.grp = drg.tgdox.repgrp, time.grp = factor(drg.tgdox$TimePoint))
saveRDS.gz(drg.tgdox.out, "./save/drg.tgdox.out.rda")

drg.tgdox.genes <- data.frame(Symbol = rownames(drg.tgdox), Hotelling = drg.tgdox.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.tgdox.genes, "./drg.tgdox.out.xlsx")

drg.tgnd <- drg.import[,drg.import$Gt.Treatment == "Tg.ND"]
drg.tgnd.reps <- rep(ncol(drg.tgnd)/5, nrow(drg.tgnd))
drg.tgnd.repgrp <- rep(1:4, (ncol(drg.tgnd)/4))
drg.tgnd.out <- mb.long(drg.tgnd, method = "1D", type = "robust", times = 5, reps = drg.tgnd.reps, rep.grp = drg.tgnd.repgrp, time.grp = factor(drg.tgnd$TimePoint))
saveRDS.gz(drg.tgnd.out, "./save/drg.tgnd.out.rda")

drg.tgnd.genes <- data.frame(Symbol = rownames(drg.tgnd), Hotelling = drg.tgnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.tgnd.genes, "./drg.tgnd.out.xlsx")

drg.wtdox <- drg.import[,drg.import$Gt.Treatment == "WT.DOX" & drg.import$Replicate < 4 ]
drg.wtdox.reps <- rep(ncol(drg.wtdox)/5, nrow(drg.wtdox))
drg.wtdox.repgrp <- rep(1:3, (ncol(drg.wtdox)/3))
drg.wtdox.out <- mb.long(drg.wtdox, method = "1D", type = "robust", times = 5, reps = drg.wtdox.reps, rep.grp = drg.wtdox.repgrp, time.grp = factor(drg.wtdox$TimePoint))
saveRDS.gz(drg.wtdox.out, "./save/drg.wtdox.out.rda")

drg.wtdox.genes <- data.frame(Symbol = rownames(drg.wtdox), Hotelling = drg.wtdox.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.wtdox.genes, "./drg.wtdox.out.xlsx")

drg.doxnd <- drg.import[,grepl("Tg.DOX$|Tg.ND", drg.import$Gt.Treatment)]
drg.doxnd.reps <- cbind(rep(length(which(drg.doxnd$Gt.Treatment == "Tg.DOX"))/5, nrow(drg.doxnd)), rep(length(which(drg.doxnd$Gt.Treatment == "Tg.ND"))/5, nrow(drg.doxnd)))
drg.doxnd.out <- mb.long(drg.doxnd, method = "2D", type = "robust", times = 5, reps = drg.doxnd.reps, rep.grp = drg.doxnd$Replicate, time.grp = factor(drg.doxnd$TimePoint), condition.grp = factor(drg.doxnd$Gt.Treatment))
saveRDS.gz(drg.doxnd.out, "./save/drg.doxnd.out.rda")

drg.doxnd.genes <- data.frame(Symbol = rownames(drg.doxnd), Hotelling = drg.doxnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.doxnd.genes, "./drg.doxnd.out.xlsx")

drg.tgwt <- drg.import[,grepl("Tg.DOX$|WT.DOX", drg.import$Gt.Treatment)]
drg.tgwt <- drg.tgwt[,!(drg.tgwt$Gt.Treatment == "WT.DOX" & drg.tgwt$Replicate == 4)]
drg.tgwt.reps <- cbind(rep(length(which(drg.tgwt$Gt.Treatment == "Tg.DOX"))/5, nrow(drg.tgwt)), rep(length(which(drg.tgwt$Gt.Treatment == "WT.DOX"))/5, nrow(drg.tgwt)))
drg.tgwt.out <- mb.long(drg.tgwt, method = "2D", type = "robust", times = 5, reps = drg.tgwt.reps, rep.grp = drg.tgwt$Replicate, time.grp = factor(drg.tgwt$TimePoint), condition.grp = factor(drg.tgwt$Gt.Treatment))
saveRDS.gz(drg.tgwt.out, "./save/drg.tgwt.out.rda")

drg.tgwt.genes <- data.frame(Symbol = rownames(drg.tgwt), Hotelling = drg.tgwt.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.tgwt.genes, "./drg.tgwt.out.xlsx")

drg.control <- drg.import[,grepl("WT.DOX|Tg.ND", drg.import$Gt.Treatment)]
drg.control <- drg.control[,!(drg.control$Gt.Treatment == "WT.DOX" & drg.control$Replicate == 4)]
drg.control.reps <- cbind(rep(length(which(drg.control$Gt.Treatment == "Tg.ND"))/5, nrow(drg.control)), rep(length(which(drg.control$Gt.Treatment == "WT.DOX"))/5, nrow(drg.control)))
drg.control.out <- mb.long(drg.control, method = "2D", type = "robust", times = 5, reps = drg.control.reps, rep.grp = drg.control$Replicate, time.grp = factor(drg.control$TimePoint), condition.grp = factor(drg.control$Gt.Treatment))
saveRDS.gz(control.out, "./save/control.out.rda")

drg.control.genes <- data.frame(Symbol = rownames(drg.control), Hotelling = drg.control.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(drg.control.genes, "./drg.control.out.xlsx")

cutoff <- 200
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
cerebellum.tgdox.reps <- rep(ncol(cerebellum.tgdox)/5, nrow(cerebellum.tgdox))
cerebellum.tgdox.repgrp <- rep(1:4, (ncol(cerebellum.tgdox)/4))
cerebellum.tgdox.out <- mb.long(cerebellum.tgdox, method = "1D", type = "robust", times = 5, reps = cerebellum.tgdox.reps, rep.grp = cerebellum.tgdox.repgrp, time.grp = factor(cerebellum.tgdox$TimePoint))
saveRDS.gz(cerebellum.tgdox.out, "./save/cerebellum.tgdox.out.rda")

cerebellum.tgdox.genes <- data.frame(Symbol = rownames(cerebellum.tgdox), Hotelling = cerebellum.tgdox.out$HotellingT2) %>% arrange(desc(Hotelling)) #%>% slice(1:500)
write.xlsx(cerebellum.tgdox.genes, "./cerebellum.tgdox.out.xlsx")

cerebellum.tgnd <- cerebellum.import[,cerebellum.import$Gt.Treatment == "Tg.ND"]
cerebellum.tgnd.reps <- rep(ncol(cerebellum.tgnd)/5, nrow(cerebellum.tgnd))
cerebellum.tgnd.repgrp <- rep(1:4, (ncol(cerebellum.tgnd)/4))
cerebellum.tgnd.out <- mb.long(cerebellum.tgnd, method = "1D", type = "robust", times = 5, reps = cerebellum.tgnd.reps, rep.grp = cerebellum.tgnd.repgrp, time.grp = factor(cerebellum.tgnd$TimePoint))
saveRDS.gz(cerebellum.tgnd.out, "./save/cerebellum.tgnd.out.rda")

cerebellum.tgnd.genes <- data.frame(Symbol = rownames(cerebellum.tgnd), Hotelling = cerebellum.tgnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(cerebellum.tgnd.genes, "./cerebellum.tgnd.out.xlsx")

#Problemas!
cerebellum.wtdox <- cerebellum.import[,cerebellum.import$Gt.Treatment == "WT.DOX" ]
cerebellum.wtdox.reps <- rep(ncol(cerebellum.wtdox)/5, nrow(cerebellum.wtdox))
cerebellum.wtdox.repgrp <- rep(1:3, (ncol(cerebellum.wtdox)/3))
cerebellum.wtdox.out <- mb.long(cerebellum.wtdox, method = "1D", type = "robust", times = 5, reps = cerebellum.wtdox.reps, rep.grp = cerebellum.wtdox.repgrp, time.grp = factor(cerebellum.wtdox$TimePoint))
saveRDS.gz(cerebellum.wtdox.out, "./save/cerebellum.wtdox.out.rda")

cerebellum.wtdox.genes <- data.frame(Symbol = rownames(cerebellum.wtdox), Hotelling = cerebellum.wtdox.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(cerebellum.wtdox.genes, "./cerebellum.wtdox.out.xlsx")

cerebellum.doxnd <- cerebellum.import[,grepl("Tg.DOX$|Tg.ND", cerebellum.import$Gt.Treatment)]
cerebellum.doxnd.reps <- cbind(rep(length(which(cerebellum.doxnd$Gt.Treatment == "Tg.DOX"))/5, nrow(cerebellum.doxnd)), rep(length(which(cerebellum.doxnd$Gt.Treatment == "Tg.ND"))/5, nrow(cerebellum.doxnd)))
cerebellum.doxnd.out <- mb.long(cerebellum.doxnd, method = "2D", type = "robust", times = 5, reps = cerebellum.doxnd.reps, rep.grp = cerebellum.doxnd$Replicate, time.grp = factor(cerebellum.doxnd$TimePoint), condition.grp = factor(cerebellum.doxnd$Gt.Treatment))
saveRDS.gz(cerebellum.doxnd.out, "./save/cerebellum.doxnd.out.rda")

cerebellum.doxnd.genes <- data.frame(Symbol = rownames(cerebellum.doxnd), Hotelling = cerebellum.doxnd.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(cerebellum.doxnd.genes, "./cerebellum.doxnd.out.xlsx")

cerebellum.tgwt <- cerebellum.import[,grepl("Tg.DOX$|WT.DOX", cerebellum.import$Gt.Treatment)]
cerebellum.tgwt <- cerebellum.tgwt[,!(cerebellum.tgwt$Gt.Treatment == "WT.DOX" & cerebellum.tgwt$Replicate == 4)]
cerebellum.tgwt.reps <- cbind(rep(length(which(cerebellum.tgwt$Gt.Treatment == "Tg.DOX"))/5, nrow(cerebellum.tgwt)), rep(length(which(cerebellum.tgwt$Gt.Treatment == "WT.DOX"))/5, nrow(cerebellum.tgwt)))
cerebellum.tgwt.out <- mb.long(cerebellum.tgwt, method = "2D", type = "robust", times = 5, reps = cerebellum.tgwt.reps, rep.grp = cerebellum.tgwt$Replicate, time.grp = factor(cerebellum.tgwt$TimePoint), condition.grp = factor(cerebellum.tgwt$Gt.Treatment))
saveRDS.gz(cerebellum.tgwt.out, "./save/cerebellum.tgwt.out.rda")

cerebellum.tgwt.genes <- data.frame(Symbol = rownames(cerebellum.tgwt), Hotelling = cerebellum.tgwt.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(cerebellum.tgwt.genes, "./cerebellum.tgwt.out.xlsx")

cerebellum.control <- cerebellum.import[,grepl("WT.DOX|Tg.ND", cerebellum.import$Gt.Treatment)]
cerebellum.control <- cerebellum.control[,!(cerebellum.control$Gt.Treatment == "WT.DOX" & cerebellum.control$Replicate == 4)]
cerebellum.control.reps <- cbind(rep(length(which(cerebellum.control$Gt.Treatment == "Tg.ND"))/5, nrow(cerebellum.control)), rep(length(which(cerebellum.control$Gt.Treatment == "WT.DOX"))/5, nrow(cerebellum.control)))
cerebellum.control.out <- mb.long(cerebellum.control, method = "2D", type = "robust", times = 5, reps = cerebellum.control.reps, rep.grp = cerebellum.control$Replicate, time.grp = factor(cerebellum.control$TimePoint), condition.grp = factor(cerebellum.control$Gt.Treatment))
saveRDS.gz(control.out, "./save/control.out.rda")

cerebellum.control.genes <- data.frame(Symbol = rownames(cerebellum.control), Hotelling = cerebellum.control.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% slice(1:500)
write.xlsx(cerebellum.control.genes, "./cerebellum.control.out.xlsx")

cutoff <- 200
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
