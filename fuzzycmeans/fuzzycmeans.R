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
library(sva)
library(peer)
library(limma)

#Longitudinal analysis
library(Mfuzz)

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

csel.repeat <- function(eset, m, crange, max.runs, csel.mat = matrix(), count.runs = 0)
{
    csel.out <- cselection(eset, m = m, crange = crange, repeats = 5, visu = FALSE)
    if (count.runs == 0)
    {
        csel.return <- csel.out
    }
    else
    {
        csel.return <- rbind(csel.mat, csel.out)
    }
    count.runs <- count.runs + 1
    print(count.runs)
    if (count.runs < max.runs)
    {
        csel.repeat(eset, m, crange, max.runs, csel.return, count.runs)
    }
    else
    {
        return(csel.return)
    }
}

dmin.repeat <- function(eset, m, crange, max.runs, dmin.mat = vector(), count.runs = 0)
{
    dmin.out <- Dmin(eset, m = m, crange = crange, repeats = 5, visu = FALSE)
    if (count.runs == 0)
    {
        dmin.return <- dmin.out
    }
    else
    {
        dmin.return <- rbind(dmin.mat, dmin.out)
    }
    count.runs <- count.runs + 1
    if (count.runs < max.runs)
    {
        dmin.repeat(eset, m, crange, max.runs, dmin.return, count.runs)
    }
    else
    {
        return(dmin.return)
    }
}

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

intensities.mean <- readRDS.gz("../dtw/save/intensities.means.rda") 
intensities.tgdox <- intensities.mean$Tg.DOX

tgdox.betr.genes <- read.xlsx("../betr/tgdox.out.xlsx") %>% select(Symbol)
tgdox.timecourse.genes <- read.xlsx("../timecourse/tgdox.out.xlsx") %>% select(Symbol)
longnet.tgdox <- readRDS.gz("../longnet/save/ica.tgdox.rda")
ica.tgdox <- readRDS.gz("../longnet/save/ica.tgdox.rda")

tgdox.genes <- c(unlist(tgdox.betr.genes), unlist(tgdox.timecourse.genes)) %>% unique %>% match.exact
#tgdox.genes <- match.exact(ica.tgdox)

tgdox.cluster <- intensities.tgdox[grepl(tgdox.genes, rownames(intensities.tgdox)),] 
tgdox.cluster <- intensities.tgdox[match(longnet.tgdox, rownames(intensities.tgdox)),] 
tgdox.eset <- ExpressionSet(assayData = as.matrix(tgdox.cluster)) %>% standardise

m.estimate <- mestimate(tgdox.eset) 

csel.runs <- csel.repeat(tgdox.eset, m = m.estimate, crange = 4:20, max.runs = 5)
csel.ratio <- (4:20) - colMeans(csel.runs) 

dmin.runs <- dmin.repeat(tgdox.eset, m = m.estimate, crange = seq(4, 40, 4), max.runs = 5)

seed.mfuzz <- function(eset, c.num, m, mfuzz.list = list(), iter.count = 0)
{
   if (iter.count < 250) 
   {
        mfuzz.new <- mfuzz(eset, c = c.num, m = m)
        iter.count <- iter.count + 1
        print(iter.count)
        mfuzz.add <- mfuzz.new$membership
        colnames(mfuzz.add) <- paste("X", 1:c.num, sep = "")
        mfuzz.list[[iter.count]] <- mfuzz.add
        seed.mfuzz(eset = eset, c.num = c.num, m = m, mfuzz.list = mfuzz.list, iter.count = iter.count)
   }
   else
   {
        return(mfuzz.list)
   }
}

cluster.tgdox7 <- seed.mfuzz(eset = tgdox.eset, c.num = 4, m = m.estimate)
median.tgdox7 <- melt(cluster.tgdox7) %>% dcast(Var1 ~ Var2, median)

cluster.tgdox13 <- seed.mfuzz(eset = tgdox.eset, c.num = 13, m = m.estimate)
median.tgdox13 <- melt(cluster.tgdox13) %>% dcast(Var1 ~ Var2, median)

cluster.tgdox14 <- seed.mfuzz(eset = tgdox.eset, c.num = 14, m = m.estimate)
median.tgdox14 <- melt(cluster.tgdox14) %>% dcast(Var1 ~ Var2, median)

cluster.tgdox <- mfuzz(tgdox.eset, c = 7, m = m.estimate)
mfuzz.plot(tgdox.eset, cl = cluster.tgdox, mfrow = c(4,4), time.labels = 1:4)
#cluster.filter <- (cluster.patient$membership > 0.4) %>% apply(1, any)
#cluster.filtered <- cluster.patient$membership[cluster.filter,]
#patient.members <- apply(cluster.patient$membership, 1, which.max)
patient.df <- data.frame(Symbol = names(cluster.patient$cluster), Cluster = cluster.patient$cluster)
write.xlsx(patient.df, "./patient.cmeans.xlsx")

partcoef.tgdox <- partcoef(tgdox.eset)

permut.timeseries <- function()
