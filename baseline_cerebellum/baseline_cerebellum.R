#Some utils
library(R.utils)

#For DE analysis
library(limma)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)
library(annotate)
library(lumi)
library(WGCNA) #for fastcor
library(siggenes)

#For batch correction and PEER
library(sva)

#Data arrangement
library(dplyr)
library(parallel)

#Functional programming
library(magrittr)
library(purrr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Reading and writing tables
library(readr)
library(openxlsx)

saveRDS.gz <- function(object, file, threads=parallel::detectCores()) 
{
    con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
    saveRDS(object, file = con)
    close(con)
}

readRDS.gz <- function(file, threads=parallel::detectCores()) 
{
    con <- pipe(paste0("pigz -d -c -p",threads," ",file))
    object <- readRDS(file = con)
    close(con)
    return(object)
}

#Boxplot function
gen.boxplot <- function(filename, lumi.object, colorscheme, maintext, ylabtext)
{
    #dataset %<>% t %>% data.frame
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Heatmap bar function
gen.heatmapbars <- function(batch.colorscheme, genotype.colorscheme, treatment.colorscheme, tissue.colorscheme, targetset)
{
    batch.heatmap <- data.frame("Batch" = seq(1:length(batch.colorscheme)), "Batch.Color" = batch.colorscheme)
    genotype.heatmap <- data.frame("Genotype" = levels(factor(targetset$Gt)), "Genotype.Color" = genotype.colorscheme)
    treatment.heatmap <- data.frame("Treatment" = levels(targetset$Treatment), "Treatment.Color" = treatment.colorscheme)
    tissue.heatmap <- data.frame("Tissue" = levels(factor(targetset$Tissue)), "Tissue.Color" = tissue.colorscheme)
    colorscheme <- data.frame("Batch" = targetset$Batch, "Genotype" = targetset$Gt, "Treatment" = targetset$Treatment, "Tissue" = targetset$Tissue) %>% join(batch.heatmap) %>% join(genotype.heatmap) %>% join(treatment.heatmap) %>% join(tissue.heatmap)
    colorscheme <- as.matrix(select(colorscheme, Batch.Color, Genotype.Color, Treatment.Color, Tissue.Color))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, lumi.object, maintitle)
{
    intensities1.cor <- corFast(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = cbind(lumi.object$Batch.Color, lumi.object$Genotype.Color, lumi.object$Treatment.Color, lumi.object$Tissue.Color), scale = "none", cexCol = 0.17, cexRow = 0.17, main = maintitle)
    dev.off()
}

#IAC detection of outliers vix 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = corFast(dataset, use = "p")
    cluster1 = flashClust::hclust(as.dist(1 - IAC))
    CairoPDF(filename, width = 13, height = 10)
    plot(cluster1, main = paste(maintitle, " (no = ", dim(IAC)[2], ")"))
    dev.off()
    return(IAC)
}

#Create plot of standard deviations of all interarray correlations.  
gen.sdplot <- function(filename, dataset, maintitle)
{
    meanIAC <- apply(dataset, 2, mean)
    sdCorr <- sd(meanIAC)
    numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
    numbersd.plot <- data.frame(Sample.Num = 1:ncol(dataset), Z.score = numbersd, Sample.Status = colnames(dataset))

    p <- ggplot(numbersd.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Status) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(numbersd)
}

#Median absolute deviation standardization function
standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes
gen.topgenes <- function(filename, lumi.object, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    dataset <- exprs(lumi.object)
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- flashClust::hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", ColSideColors = cbind(lumi.object$Batch.Color, lumi.object$Genotype.Color, lumi.object$Treatment.Color, lumi.object$Tissue.Color), cexRow = 0.08, cexCol = 0.08)
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$External.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) 
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Run statistical cutoff tests
gen.decide <- function(test, fit.object, write.results)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] #Eliminate the row for no change in expression
    mysum[2,] <- -(mysum[2,])
    colnames(mysum) <- c("Tg.DOX_vs_Tg.ND", "Tg.DOX_vs_WT.DOX", "Tg.DOXR_vs_Tg.DOX")
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot, width.plot = 6, height.plot = 7)
{
    decide.plot$variable <- str_replace_all(decide.plot$variable, "_", " ")
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = max(value) + 110, label = value), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = min(value) - 110, label = abs(value)), hjust = 1.3, position = position_dodge(width = 1))
    if (length(unique(decide.plot$Test)) > 1)
    {
        p <- p + facet_grid(Num + Test ~ .) 
        #p <- p + ggtitle("Threshold Selection")
    }
    #else
    #{
        #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    #}
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))# + ylab("Differentially Expressed Genes")
    CairoPDF(filename, width = width.plot, height = height.plot)
    print(p)
    dev.off()
}

gen.pval.hist <- function(filename, fit.pvals)
{
    colnames(fit.pvals) <- c("Tg.DOX_vs_Tg.ND", "Tg.DOX_vs_WT.DOX", "Tg.DOXR_vs_Tg.DOX")
    fit.pvals.plot <- melt(fit.pvals)
    fit.pvals.plot$Contrasts <- str_replace_all(fit.pvals.plot$Contrasts, "_", " ")
    p <- ggplot(fit.pvals.plot, aes(x = value)) + geom_histogram(binwidth = 1/80) + facet_grid(. ~ Contrasts)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("P-value distribution across contrasts") + theme(axis.title.y = element_blank()) + xlab("p-value")
    CairoPDF(filename, height = 7, width = 21)
    print(p)
    dev.off()
}

gen.venndiagram <- function(filename, results)
{
    sumV <- colSums(summary(results)[-2,])
    v <- paste(c("Tg.DOX vs Tg.ND", "Tg.DOX vs WT.DOX", "Tg.DOXR vs Tg.DOX"), " (", sumV, ")", sep = "")
    CairoPDF(filename, width = 6, height = 6)
    vennDiagram(results, names = v, main = "", include = c("up", "down"), counts.col=c(2,3), cex = 0.8)
    dev.off()
}

gen.ratios <- function(lumi.final)
{
    all.tgdox <- exprs(lumi.final[,lumi.final$Gt.Treatment == "Tg.DOX"])
    all.tgdox.means <- rowMeans(all.tgdox)
    all.tgnd <- exprs(lumi.final[,lumi.final$Gt.Treatment == "Tg.ND"])
    all.tgnd.means <- rowMeans(all.tgnd)
    all.wtdox <- exprs(lumi.final[,lumi.final$Gt.Treatment == "WT.DOX"])
    all.wtdox.means <- rowMeans(all.wtdox)
    all.tgdoxr <- exprs(lumi.final[,lumi.final$Gt.Treatment == "Tg.DOXR"])

    coef.tgdox.tgnd <- all.tgdox - all.tgnd.means
    colnames(coef.tgdox.tgnd) <- paste(colnames(coef.tgdox.tgnd), ".doxnd", sep = "")
    coef.tgdox.wtdox <- all.tgdox - all.wtdox.means
    colnames(coef.tgdox.wtdox) <- paste(colnames(coef.tgdox.wtdox), ".tgwt", sep = "")
    coef.tgdoxr.tgdox <- all.tgdoxr - all.tgdox.means
    colnames(coef.tgdoxr.tgdox) <- paste(colnames(coef.tgdoxr.tgdox), ".rescue", sep = "")

    coefs <- data.frame("Symbol" = rownames(coef.tgdox.tgnd), coef.tgdox.tgnd, coef.tgdox.wtdox, coef.tgdoxr.tgdox)
    all.samples <- data.frame("Symbol" = rownames(lumi.final), exprs(lumi.final))
    colnames(all.samples)[2:ncol(all.samples)] <- paste(colnames(all.samples[2:ncol(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(coefs, all.samples)
    return(ratio.exp)
}

gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(Tg.DOX_vs_Tg.ND = Tg.DOX - Tg.ND, Tg.DOX_vs_WT.DOX = Tg.DOX - WT.DOX, Tg.DOXR_vs_Tg.DOX = Tg.DOXR - Tg.DOX, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
    
    top.object.doxnd <- topTable(fitb, coef = 1, n = Inf)
    saveRDS.gz(top.object.doxnd, "./save/top.object.doxnd.rda")
    top.object.tgwt <- topTable(fitb, coef = 2, n = Inf) 
    saveRDS.gz(top.object.tgwt, "./save/top.object.tgwt.rda")
    top.object.rescue <- topTable(fitb, coef = 3, n = Inf) 
    saveRDS.gz(top.object.rescue, "./save/top.object.rescue.rda")

    return(fitb)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:7, widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    
    fitsel.ratio.all <- merge(treat.de, ratio.exp)
    fitsel.return.all <- select(fitsel.ratio.all, Symbol, contains("Coef."), contains("p.value."), F, F.p.value, contains("Res."), contains("t.time1."), A, matches("doxnd|tgwt|rescue|expr")) %>% arrange(desc(F))

    anovalist <- select(fitsel.return.all, contains("Res.")) %>% apply(1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef.") %>% which
    colnames(fitsel.return)[coef.cols] <- c("Coef.Tg.Dox vs. Tg.ND", "Coef.Tg.Dox vs. Wt.Dox", "Coef.Tg.DoxR vs. Tg.Dox")
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, "_time1.xlsx", sep = ""))

    #Carrier - Control
    doxnd.abs <- as.numeric(fitsel.return$"Coef.Tg.Dox vs. Tg.ND") %>% abs
    fitsel.doxnd <- select(fitsel.return, Symbol, matches("Coef.Tg.Dox vs. Tg.ND"), matches("Coef.Tg.Dox vs. Wt.Dox"), matches("Coef.Tg.DoxR vs. Tg.Dox")) %>% mutate(doxnd.abs) %>% arrange(desc(doxnd.abs)) %>% select(-doxnd.abs)
    gen.small.workbook(fitsel.doxnd, paste("./significant_geneList_", suffix, "_time1_doxnd.xlsx", sep = ""))

    #Patient - Control
    tgwt.abs <- as.numeric(fitsel.return$"Coef.Tg.Dox vs. Wt.Dox") %>% abs
    fitsel.tgwt <- select(fitsel.return, Symbol, matches("Coef.Tg.Dox vs. Wt.Dox"), matches("Coef.Tg.Dox vs. Tg.ND"), matches("Coef.Tg.DoxR vs. Tg.Dox")) %>% mutate(tgwt.abs) %>% arrange(desc(tgwt.abs)) %>% select(-tgwt.abs)
    gen.small.workbook(fitsel.tgwt, paste("./significant_geneList_", suffix, "_time1_tgwt.xlsx", sep = ""))

    #Patient - Carrier
    rescue.abs <- as.numeric(fitsel.return$"Coef.Tg.DoxR vs. Tg.Dox") %>% abs
    fitsel.rescue <- select(fitsel.return, Symbol, matches("Coef.Tg.DoxR vs. Tg.Dox"), matches("Coef.Tg.Dox vs. Wt.Dox"), matches("Coef.Tg.Dox vs. Tg.ND")) %>% mutate(rescue.abs) %>% arrange(desc(rescue.abs)) %>% select(-rescue.abs)
    gen.small.workbook(fitsel.rescue, paste("./significant_geneList_", suffix, "_time1_rescue.xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_time1_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:4, widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}
    
#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.tgnd <- filter(dataset, Res.Tg.DOX_vs_Tg.ND != 0) %>% select(matches("doxnd"))
    plot.wtdox <- filter(dataset, Res.Tg.DOX_vs_WT.DOX != 0) %>% select(matches("tgwt"))
    plot.rescue <- filter(dataset, Res.Tg.DOXR_vs_Tg.DOX != 0) %>% select(matches("rescue"))
    gen.anova.heatmap(paste("./anova_heatmap_tg.dox_vs_tg.nd", suffix, sep = "_"), plot.tgnd, "Tg.DOX vs. Tg.ND") 
    gen.anova.heatmap(paste("./anova_heatmap_tg.dox_vs_wt.dox", suffix, sep = "_"), plot.wtdox, "Tg.DOX vs. WT.DOX")
    gen.anova.heatmap(paste("./anova_heatmap_tg.doxr_vs_tg.dox", suffix, sep = "_"), plot.rescue, "Tg.DOXR vs. Tg.DOX")
    return(list(plot.tgnd, plot.wtdox, plot.rescue))
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 8, height = 8)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", labCol = "", keysize = 0.9)
    dev.off()
}

#GO functions
enrichr.submit <- function(colname, dataset, enrichr.terms, subdir)
{
    filter.cond <- paste(paste(colname, '==', '1'), paste(colname, "==", '-1'), sep = "|")
    dataset.submit <- filter_(dataset, filter.cond)
    #write.xlsx(dataset.down, file.path("./enrichr", subdir, paste(colname.formatted, "down.xlsx", sep = "_")))

    dir.create(file.path("./enrichr", subdir), showWarnings = FALSE)
    colname.formatted <- str_replace_all(colname, "time1\\.", "") %>% str_replace("Res\\.", "")
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset.submit, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]

    names(enrichr.data) <- enrichr.names

    map(names(enrichr.data), enrichr.wkbk, enrichr.data, colname, subdir, "enrichr")
}

enrichr.wkbk <- function(database, full.df, colname, subdir, direction)
{
    dataset <- full.df[[database]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    colname.formatted <- str_replace_all(colname, "time1\\.", "") %>% str_replace("Res\\.", "")
    
    dir.create(file.path("./enrichr", subdir, colname.formatted, direction), recursive = TRUE)
    filename = paste(file.path("./enrichr", subdir, colname.formatted, direction, database), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

lumi.vst <- readRDS.gz("../baseline/save/lumi.vst.rda")

#Cerebellum
lumi.cerebellum <- lumi.vst[,lumi.vst$Tissue == "Cerebellum"]
lumi.cerebellum.norm <- lumiN(lumi.cerebellum, method = "rsn") #Normalize with robust spline regression
lumi.cerebellum.qual <- lumiQ(lumi.cerebellum.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.cerebellum.cutoff <- detectionCall(lumi.cerebellum.qual) #Get the count of probes which passed the detection threshold per sample
lumi.cerebellum.expr <- lumi.cerebellum.qual[which(lumi.cerebellum.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.cerebellum <- getSYMBOL(rownames(lumi.cerebellum.expr), 'lumiMouseAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.cerebellum.annot <- lumi.cerebellum.expr[!symbols.lumi.cerebellum,] #Drop any probe which is not annotated
saveRDS.gz(lumi.cerebellum.annot, file = "lumi.cerebellum.annot.rda")

batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green")
gen.boxplot("cerebellum_intensity_norm.jpg", lumi.cerebellum.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
gen.heatmap("cerebellum_heatmap", lumi.cerebellum.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.cerebellum <- gen.IACcluster("cerebellum_IAC", exprs(lumi.cerebellum.annot), "Outlier removed")
sd.raw.norm.cerebellum <- gen.sdplot("cerebellum_sd", IAC.norm.cerebellum, "Outlier removed")
saveRDS.gz(IAC.norm.cerebellum, file = "./save/IAC.norm.cerebellum.rda")
saveRDS.gz(sd.raw.norm.cerebellum, file = "./save/sd.raw.norm.cerebellum.rda")

lumi.cerebellum.mads <- apply(exprs(lumi.cerebellum.annot), 1, mad)
lumi.cerebellum.ordered <- order(lumi.cerebellum.mads, decreasing = TRUE)
top500.dist.norm.cerebellum <- gen.topgenes("cerebellum_heatmap_500", lumi.cerebellum.annot, "Clustering Based on the Top 500 Most Variable Genes", lumi.cerebellum.ordered, 500)
top1000.dist.norm.cerebellum <- gen.topgenes("cerebellum_heatmap_1000", lumi.cerebellum.annot, "Clustering Based on the Top 1000 Most Variable Genes", lumi.cerebellum.ordered, 1000)

cm1.cerebellum <- cmdscale(top1000.dist.norm.cerebellum, eig = TRUE)
gen.pca("cerebellum_mds_genotype", cm1.cerebellum, phenoData(lumi.cerebellum.annot), genotype.colors, "Gt")
gen.pca("cerebellum_mds_batch", cm1.cerebellum, phenoData(lumi.cerebellum.annot), batch.colors, "Batch")
gen.pca("cerebellum_mds_treatment", cm1.cerebellum, phenoData(lumi.cerebellum.annot), treatment.colors, "Treatment")

qcsum.cerebellum <- lumi.cerebellum.annot@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum.cerebellum) %<>% str_replace("\\.0\\.01\\.", "")
qcsum.cerebellum$Sample.Name <- rownames(qcsum.cerebellum)
qcsum.cerebellum$RIN <- lumi.cerebellum.annot$RIN
qcsum.cerebellum$Sample.Num <- lumi.cerebellum.annot$TimePoint
arrange(qcsum.cerebellum, distance.to.sample.mean)

lumi.cerebellum.rmone <- lumi.cerebellum.annot[,!lumi.cerebellum.annot$Batch == 15]

#Use ComBat for batch effect correction
lumi.cerebellum.rmone$Gt.Treatment <- paste(lumi.cerebellum.rmone$Gt, lumi.cerebellum.rmone$Treatment, sep = ".")

model.gt.treatment <- model.matrix( ~ 0 + factor(lumi.cerebellum.rmone$Gt.Treatment) ) #Make dummy variables out of Sex
colnames(model.gt.treatment) <- c("Tg.DOX", "Tg.DOXR", "Tg.ND", "WT.DOX") 
model.gt.treatment.reduce <- model.gt.treatment[,-4] #Remove one column (like done with Status)

model.sex <- model.matrix( ~ 0 + factor(lumi.cerebellum.rmone$Sex) ) #Make dummy variables out of Sex
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1] #Remove one column (like done with Status)

model.combat <- data.frame(model.gt.treatment.reduce, Male = model.sex.reduce, Weight = as.numeric(lumi.cerebellum.rmone$Body.weight.at.Diss), RIN = lumi.cerebellum.rmone$RNA_RIN, Age = lumi.cerebellum.rmone$Age.at.diss)  #Assemble covariate matrix with both categorical and continuous covariates
model.combat$RIN <- model.combat$RIN - min(model.combat$RIN)

expr.combat <- ComBat(dat = exprs(lumi.cerebellum.rmone), batch = factor(lumi.cerebellum.rmone$Batch), mod = model.combat) #Run ComBat
lumi.combat <- lumi.cerebellum.rmone #Create a new lumi object as a copy of lumi.rmreps.annot
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities
saveRDS.gz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
gen.boxplot("intensity_combat.jpg", lumi.combat, batch.colors, "Batch-corrected intensity", "Intensity") #Replot intensities to make sure they look okay
gen.heatmap("combat_heatmap", lumi.combat, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.combat <- gen.IACcluster("combat_IAC", exprs(lumi.combat), "Outlier removed")
sd.raw.norm.combat <- gen.sdplot("combat_sd", IAC.norm.combat, "Outlier removed")
saveRDS.gz(IAC.norm.combat, file = "./save/IAC.norm.combat.rda")
saveRDS.gz(sd.raw.norm.combat, file = "./save/sd.raw.norm.combat.rda")

lumi.combat.mads <- apply(exprs(lumi.combat), 1, mad)
lumi.combat.ordered <- order(lumi.combat.mads, decreasing = TRUE)
top500.dist.norm.combat <- gen.topgenes("combat_heatmap_500", lumi.combat, "Clustering Based on the Top 500 Most Variable Genes", lumi.combat.ordered, 500)
top1000.dist.norm.combat <- gen.topgenes("combat_heatmap_1000", lumi.combat, "Clustering Based on the Top 1000 Most Variable Genes", lumi.combat.ordered, 1000)

cm1.combat <- cmdscale(top1000.dist.norm.combat, eig = TRUE)
gen.pca("combat_mds_genotype_treatment", cm1.combat, phenoData(lumi.combat), genotype.colors, "Gt.Treatment")
gen.pca("combat_mds_batch", cm1.combat, phenoData(lumi.combat), batch.colors, "Batch")

#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
model.cov <- select(model.combat, Male, Weight, RIN, Age) #Create model matrix of covariates to be removed
model.design <- model.matrix( ~ 0 + factor(lumi.combat$Gt.Treatment) )
colnames(model.design) %<>% str_split("\\)") %>% map_chr(`[`, 2) 
model.design.reduce <- model.design[,-4]
export.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.cov, design = model.design) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
export.lumi <- lumi.combat #Make a copy of lumi object
exprs(export.lumi) <- export.expr #Transfer cleaned expression values into new lumi object
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")

gen.boxplot("baseline_intensity_corrected.jpg", export.lumi, batch.colors, "Covariate-corrected intensity", "Intensity") #Replot intensities to make sure they look okay
gen.heatmap("cov_heatmap", export.lumi, "")
IAC.norm.cov <- gen.IACcluster("cov_IAC", exprs(export.lumi), "Outlier removed")
sd.raw.norm.cov <- gen.sdplot("cov_sd", IAC.norm.cov, "")
saveRDS.gz(IAC.norm.cov, file = "./save/IAC.norm.cov.rda")
saveRDS.gz(sd.raw.norm.cov, file = "./save/sd.raw.norm.cov.rda")

export.lumi.mads <- apply(exprs(export.lumi), 1, mad)
export.lumi.ordered <- order(export.lumi.mads, decreasing = TRUE)
top500.dist.norm.cov <- gen.topgenes("cov_heatmap_500", export.lumi, "", export.lumi.ordered, 500)
top1000.dist.norm.cov <- gen.topgenes("cov_heatmap_1000", export.lumi, "Clustering Based on the Top 1000 Most Variable Genes", export.lumi.ordered, 1000)

cm1.cov <- cmdscale(top1000.dist.norm.cov, eig = TRUE)
gen.pca("cov_mds_genotype_treatment", cm1.cov, phenoData(export.lumi), genotype.colors, "Gt.Treatment")
gen.pca("cov_mds_batch", cm1.cov, phenoData(export.lumi), batch.colors, "Batch")

#Collapse the data by symbol
lumi.combat <- lumi.combat[!is.na(getSYMBOL(featureNames(lumi.combat), 'lumiMouseAll.db')),] #Remove the probes with missing symbols (this isn't supposed to be necessary!)
expr.collapse <- collapseRows(exprs(lumi.combat), factor(fdata$SYMBOL), rownames(lumi.combat), method = "function", methodFunction = colMeans)$datETcollapsed #collapseRows by symbol
lumi.collapse <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.combat))
saveRDS.gz(expr.collapse, file = "./save/expr.collapse.rda")

#lumi.time4 <- lumi.collapse[,lumi.collapse$TimePoint == "T4"] 
#lumi.time5 <- lumi.collapse[,lumi.collapse$TimePoint == "T5"]

#model.time4 <- model.matrix( ~ 0 + factor(lumi.time4$Gt.Treatment) )
#colnames(model.time4) <- c("Tg.DOX", "Tg.DOXR", "Tg.ND", "WT.DOX")
#model.sex4 <- model.matrix( ~ 0 + factor(lumi.time4$Sex) )
#colnames(model.sex4) <- c("Female", "Male")
#model.sex4.reduce <- model.sex4[,-1]
#model.cov4 <- cbind(Male = model.sex4.reduce, Weight = lumi.time4$Body.weight.at.Diss, RIN = lumi.time4$RNA_RIN, Age = lumi.time4$Age.at.diss)

#model.time5 <- model.matrix( ~ 0 + factor(lumi.time5$Gt.Treatment) )
#colnames(model.time5) <- c("Tg.DOX", "Tg.DOXR", "Tg.ND", "WT.DOX")
#model.sex5 <- model.matrix( ~ 0 + factor(lumi.time5$Sex) )
#colnames(model.sex5) <- c("Female", "Male")
#model.sex5.reduce <- model.sex5[,-1]
#model.cov5 <- cbind(Male = model.sex5.reduce, Weight = lumi.time5$Body.weight.at.Diss, RIN = lumi.time5$RNA_RIN, Age = lumi.time5$Age.at.diss)

#model.limma4 <- data.frame(model.time4, model.cov4) #Make covariate matrix for limma
#saveRDS.gz(model.limma4, file = "./save/model.limma4.rda")

#model.limma5 <- data.frame(model.time5, model.cov5) #Make covariate matrix for limma
#saveRDS.gz(model.limma5, file = "./save/model.limma5.rda")

##Fit limma linear model 
#fit.object4 <- gen.fit(exprs(lumi.time4), model.limma4)
#saveRDS.gz(fit.object4, file = "./save/fit.object4.rda")
#fit.object5 <- gen.fit(exprs(lumi.time5), model.limma5)
#saveRDS.gz(fit.object5, file = "./save/fit.object5.rda")

#ebam.doxnd <- limma2ebam(fit.object5, coef = 1)
#ebam.doxnd.df <- data.frame(Z.score = ebam.doxnd@z, Posterior = ebam.doxnd@posterior)
#ebam.doxnd.df$Significant <- ebam.doxnd.df$Posterior > 0.9
#ebam.doxnd.df$Symbol <- rownames(ebam.doxnd.df)
#saveRDS.gz(ebam.doxnd.df, "ebam.doxnd.df.rda")

#ebam.tgwt <- limma2ebam(fit.object5, coef = 2)
#ebam.tgwt.df <- data.frame(Z.score = ebam.tgwt@z, Posterior = ebam.tgwt@posterior)
#ebam.tgwt.df$Significant <- ebam.tgwt.df$Posterior > 0.9
#ebam.tgwt.df$Symbol <- rownames(ebam.tgwt.df)
#saveRDS.gz(ebam.tgwt.df, "ebam.tgwt.df.rda")

#ebam.rescue <- limma2ebam(fit.object5, coef = 3)
#ebam.rescue.df <- data.frame(Z.score = ebam.rescue@z, Posterior = ebam.rescue@posterior)
#ebam.rescue.df$Significant <- ebam.rescue.df$Posterior > 0.9
#ebam.rescue.df$Symbol <- rownames(ebam.rescue.df)
#saveRDS.gz(ebam.rescue.df, "ebam.rescue.df.rda")

lumi.doxnd <- lumi.combat[,grepl("Tg.DOX$|Tg.ND", lumi.combat$Gt.Treatment)]
lumi.tgwt <- lumi.combat[,grepl("Tg.DOX$|WT.DOX", lumi.combat$Gt.Treatment)]
lumi.rescue <- lumi.combat[,grepl("Tg.DOXR|Tg.DOX", lumi.combat$Gt.Treatment)]

model.doxnd <- model.matrix(~ Sex + Age.at.diss + Body.weight.at.Diss + RNA_RIN, data = pData(lumi.doxnd)) #Create model matrix of covariates to be removed
rmcov.doxnd.expr <- removeBatchEffect(exprs(lumi.doxnd), covariates = model.doxnd[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.doxnd <- lumi.doxnd #Make a copy of lumi object
exprs(rmcov.doxnd) <- rmcov.doxnd.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
doxnd.collapse.expr <- collapseRows(exprs(rmcov.doxnd), getSYMBOL(featureNames(rmcov.doxnd), 'lumiMouseAll.db'), rownames(rmcov.doxnd)) #collapseRows by symbol
doxnd.collapse <- ExpressionSet(assayData = doxnd.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.doxnd))
saveRDS.gz(doxnd.collapse, file = "./save/doxnd.collapse.rda")

model.tgwt <- model.matrix(~ Sex + Age.at.diss + Body.weight.at.Diss + RNA_RIN, data = pData(lumi.tgwt)) #Create model matrix of covariates to be removed
rmcov.tgwt.expr <- removeBatchEffect(exprs(lumi.tgwt), covariates = model.tgwt[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.tgwt <- lumi.tgwt #Make a copy of lumi object
exprs(rmcov.tgwt) <- rmcov.tgwt.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
tgwt.collapse.expr <- collapseRows(exprs(rmcov.tgwt), getSYMBOL(featureNames(rmcov.tgwt), 'lumiMouseAll.db'), rownames(rmcov.tgwt)) #collapseRows by symbol
tgwt.collapse <- ExpressionSet(assayData = tgwt.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.tgwt))
saveRDS.gz(tgwt.collapse, file = "./save/tgwt.collapse.rda")

model.rescue <- model.matrix(~ Sex + Age.at.diss + Body.weight.at.Diss + RNA_RIN, data = pData(lumi.rescue)) #Create model matrix of covariates to be removed
rmcov.rescue.expr <- removeBatchEffect(exprs(lumi.rescue), covariates = model.rescue[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.rescue <- lumi.rescue #Make a copy of lumi object
exprs(rmcov.rescue) <- rmcov.rescue.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
rescue.collapse.expr <- collapseRows(exprs(rmcov.rescue), getSYMBOL(featureNames(rmcov.rescue), 'lumiMouseAll.db'), rownames(rmcov.rescue)) #collapseRows by symbol
rescue.collapse <- ExpressionSet(assayData = rescue.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.rescue))
saveRDS.gz(rescue.collapse, file = "./save/rescue.collapse.rda")

doxnd.collapse$Gt.Treatment %<>% factor(levels = c("Tg.ND", "Tg.DOX"))
a0.doxnd <- find.a0(doxnd.collapse, as.integer(doxnd.collapse$Gt.Treatment), B = 1000, rand = 12345)
ebam.doxnd <- ebam(a0.doxnd)
ebam.doxnd.df <- data.frame(Z.score = ebam.doxnd@z, Posterior = ebam.doxnd@posterior)
ebam.doxnd.df$Significant <- ebam.doxnd.df$Posterior > 0.9
ebam.doxnd.df$Symbol <- rownames(ebam.doxnd.df)
saveRDS.gz(ebam.doxnd.df, "ebam.doxnd.df.rda")

tgwt.collapse$Gt.Treatment %<>% factor(levels = c("WT.DOX", "Tg.DOX"))
a0.tgwt <- find.a0(tgwt.collapse, as.integer(tgwt.collapse$Gt.Treatment), B = 1000, rand = 12345)
ebam.tgwt <- ebam(a0.tgwt)
ebam.tgwt.df <- data.frame(Z.score = ebam.tgwt@z, Posterior = ebam.tgwt@posterior)
ebam.tgwt.df$Significant <- ebam.tgwt.df$Posterior > 0.9
ebam.tgwt.df$Symbol <- rownames(ebam.tgwt.df)
saveRDS.gz(ebam.tgwt.df, "ebam.tgwt.df.rda")

rescue.collapse$Gt.Treatment %<>% factor(levels = c("Tg.DOX", "Tg.DOXR"))
a0.rescue <- find.a0(rescue.collapse, as.integer(rescue.collapse$Gt.Treatment), B = 1000, rand = 12345)
ebam.rescue <- ebam(a0.rescue)
ebam.rescue.df <- data.frame(Z.score = ebam.rescue@z, Posterior = ebam.rescue@posterior)
ebam.rescue.df$Significant <- ebam.rescue.df$Posterior > 0.9
ebam.rescue.df$Symbol <- rownames(ebam.rescue.df)
saveRDS.gz(ebam.rescue.df, "ebam.rescue.df.rda")

#Calculate ratios for use in tables
ratio.exp4 <- gen.ratios(lumi.time4)
saveRDS.gz(ratio.exp4, file = "./save/ratio.exp4.rda")
ratio.exp5 <- gen.ratios(lumi.time5)
saveRDS.gz(ratio.exp5, file = "./save/ratio.exp5.rda")

#Generate statisical cutoff
decide <- list(c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot4 <- ldply(decide, gen.decide, fit.object4, FALSE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_4", decide.plot4, 10, 10) #Plot different significance cutoffs
decide.plot5 <- ldply(decide, gen.decide, fit.object5, FALSE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection_5", decide.plot5, 10, 10) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), fit.object5, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

#Make tables
de.object <- read_tsv("./fit_none.tsv") #Read in unadjusted fit object
rownames(de.object) <- featureNames(lumi.time5)

fit.selection <- gen.tables(de.object, lumi.time5, ratio.exp5, "pLess005") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection, file = "./save/fit.selection.rda")

cutoff <- 300
all.human.genes <- read_csv("../../FRDA project/baseline_lumi/complete_genelist_time1_fdrLess01.csv")
all.human.reduce <- select(all.human.genes, Symbol:A)
rm(all.human.genes)

all.human.pca.up <- filter(all.human.reduce, Coef.time1.patient_vs_time1.carrier >= 0) %>% arrange(p.value.time1.patient_vs_time1.carrier) %>% slice(1:cutoff)
all.human.pca.down <- filter(all.human.reduce, Coef.time1.patient_vs_time1.carrier < 0) %>% arrange(p.value.time1.patient_vs_time1.carrier) %>% slice(1:cutoff)
all.human.pco.up <- filter(all.human.reduce, Coef.time1.patient_vs_time1.control >= 0) %>% arrange(p.value.time1.patient_vs_time1.control) %>% slice(1:cutoff)
all.human.pco.down <- filter(all.human.reduce, Coef.time1.patient_vs_time1.control < 0) %>% arrange(p.value.time1.patient_vs_time1.control) %>% slice(1:cutoff)
saveRDS.gz(all.human.pca.up, './save/all.human.pca.up.rda')
saveRDS.gz(all.human.pca.down, './save/all.human.pca.down.rda')
saveRDS.gz(all.human.pco.up, './save/all.human.pco.up.rda')
saveRDS.gz(all.human.pco.down, './save/all.human.pco.down.rda')

#all.mouse.genes <- read_csv("./complete_genelist_time1_pLess001.csv")
#all.mouse.reduce <- select(all.mouse.genes, Symbol:A)
all.mouse.doxnd.up <- filter(all.mouse.reduce, Coef.Tg.DOX_vs_Tg.ND >= 0) %>% arrange(p.value.Tg.DOX_vs_Tg.ND) %>% slice(1:cutoff)
all.mouse.doxnd.down <- filter(all.mouse.reduce, Coef.Tg.DOX_vs_Tg.ND < 0) %>% arrange(p.value.Tg.DOX_vs_Tg.ND) %>% slice(1:cutoff)
all.mouse.tgwt.up <- filter(all.mouse.reduce, Coef.Tg.DOX_vs_WT.DOX >= 0) %>% arrange(p.value.Tg.DOX_vs_WT.DOX) %>% slice(1:cutoff)
all.mouse.tgwt.down <- filter(all.mouse.reduce, Coef.Tg.DOX_vs_WT.DOX < 0) %>% arrange(p.value.Tg.DOX_vs_WT.DOX) %>% slice(1:cutoff)
saveRDS.gz(all.mouse.doxnd.up, './save/all.mouse.doxnd.up.rda')
saveRDS.gz(all.mouse.doxnd.down, './save/all.mouse.doxnd.down.rda')
saveRDS.gz(all.mouse.tgwt.up, './save/all.mouse.tgwt.up.rda')
saveRDS.gz(all.mouse.tgwt.down, './save/all.mouse.tgwt.down.rda')

de.mouse.doxnd.up <- all.mouse.doxnd.up$Symbol %>% toupper
de.mouse.doxnd.down <- all.mouse.doxnd.down$Symbol %>% toupper
de.mouse.tgwt.up <- all.mouse.tgwt.up$Symbol %>% toupper
de.mouse.tgwt.down <- all.mouse.tgwt.down$Symbol %>% toupper
saveRDS.gz(de.mouse.doxnd.up, "./save/de.mouse.doxnd.up.rda")
saveRDS.gz(de.mouse.doxnd.down, "./save/de.mouse.doxnd.down.rda")
saveRDS.gz(de.mouse.tgwt.up, "./save/de.mouse.tgwt.up.rda")
saveRDS.gz(de.mouse.tgwt.down, "./save/de.mouse.tgwt.down.rda")

de.human.pca.up <- all.human.pca.up$Symbol %>% toupper
de.human.pca.down <- all.human.pca.down$Symbol %>% toupper
de.human.pco.up <- all.human.pco.up$Symbol %>% toupper
de.human.pco.down <- all.human.pco.down$Symbol %>% toupper
mouse.overlap <- readRDS.gz("../dtw/save/mouse.overlap.rda")

de.doxnd.pca.up <- de.mouse.doxnd.up[de.mouse.doxnd.up %in% de.human.pca.up]
de.doxnd.pca.down <- de.mouse.doxnd.down[de.mouse.doxnd.down %in% de.human.pca.down]
de.doxnd.pco.up <- de.mouse.doxnd.up[de.mouse.doxnd.up %in% de.human.pco.up]
de.doxnd.pco.down <- de.mouse.doxnd.down[de.mouse.doxnd.down %in% de.human.pco.down]
de.tgwt.pca.up <- de.mouse.tgwt.up[de.mouse.tgwt.up %in% de.human.pca.up]
de.tgwt.pca.down <- de.mouse.tgwt.down[de.mouse.tgwt.down %in% de.human.pca.down]
de.tgwt.pco.up <- de.mouse.tgwt.up[de.mouse.tgwt.up %in% de.human.pco.up]
de.tgwt.pco.down <- de.mouse.tgwt.down[de.mouse.tgwt.down %in% de.human.pco.down]

hyper.doxnd.pca.up <- phyper(length(de.doxnd.pca.up), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.doxnd.pca.down <- phyper(length(de.doxnd.pca.down), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.doxnd.pco.up <- phyper(length(de.doxnd.pco.up), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.doxnd.pco.down <- phyper(length(de.doxnd.pco.down), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.tgwt.pca.up <- phyper(length(de.tgwt.pca.up), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.tgwt.pca.down <- phyper(length(de.tgwt.pca.down), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.tgwt.pco.up <- phyper(length(de.tgwt.pco.up), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)
hyper.tgwt.pco.down <- phyper(length(de.tgwt.pco.down), cutoff, length(mouse.overlap) - cutoff, cutoff, lower.tail = F)

hyper.list <- c(hyper.doxnd.pca.up, hyper.doxnd.pca.down, hyper.doxnd.pco.up, hyper.doxnd.pco.down, hyper.tgwt.pca.up, hyper.tgwt.pca.down, hyper.tgwt.pco.up, hyper.tgwt.pco.down)

#P-value histogram
gen.pval.hist("./hist_pvalue", fit.object5$p.value) 

#Venn diagram
results5 <- select(fit.selection, contains("Res.")) #%>% apply(2, factor)
gen.venndiagram("./venn", results5) 

#Anova heatmaps
clust.none <- gen.anova(fit.selection, "none")

#Code for getting size of objects in memory
objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.nofdr <- select(fit.selection, Symbol, Res.time1.carrier_vs_time1.control, Res.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.carrier)

comparison.cols <- names(enrichr.nofdr[-(1:2)])
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 

trap1 <- map(comparison.cols, enrichr.submit, enrichr.fdr, enrichr.terms, "fdr")
trap2 <- map(comparison.cols, enrichr.submit, enrichr.nofdr, enrichr.terms, "nofdr")

