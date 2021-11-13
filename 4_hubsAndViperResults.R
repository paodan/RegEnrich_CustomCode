
library(funcTools)
library(RegEnrich)
library(microbenchmark)

library(RegEnrich)
library(doParallel)
library(foreach)
library(ggplot2)
library(DOSE)
library(WGCNA)
library(igraph)

data(Lyme_GSE63085)
log2FPKM = log2(Lyme_GSE63085$FPKM + 1)
sampleInfo = Lyme_GSE63085$sampleInfo

data(TFs)

patientID = factor(sampleInfo$patientID, levels = unique(sampleInfo$patientID))
week = factor(sampleInfo$week, unique(sampleInfo$week))

pData = data.frame(patientID, week, row.names = rownames(sampleInfo))

# Define design matrix, contrast and coef
design = model.matrix(~0 + patientID + week, data = pData)
contrast = c(rep(0, ncol(design) - 1), 1)
coef = NULL

# Differential expression analysis to get p values and log2 fold changes.
baseCut = 10^-3  # minimum average expression cutoff
de = DEA(expr = log2FPKM, pData = pData, method = "limma", 
         designCont = list(design, contrast, coef), minMeanExpr = baseCut)
# P-value and log2 fold change
pFC = de$pFC

# Volcano plot
plot(pFC$logFC, -log10(pFC$p), xlim = c(-3, 3))
abline(h = -log(0.05))

log2FPKMhi = log2FPKM[rowMeans(log2FPKM) >= baseCut, , drop = FALSE]

regulators = unique(TFs$TF_name)

# Co-expression network
coeNet = COEN(expr = log2FPKMhi, reg = regulators, softPower = 14)

save(coeNet, file = "../rData/coeNet.Rdata")

topNet5 = topNet(network = coeNet$weightHi, percent = 5, 
                 dirrected = FALSE, reg = regulators)

p = setNames(pFC[, 2], rownames(pFC))
# Fisher exact test
resFET = regFET(object = topNet5, namedScores = p, 
                namedScoresCutoffs = 0.05, 
                minSize = 5, maxSize = 5000)
regScoreFET = rankScore(resEnrich = resFET, pFC = pFC)

regReg = regTarNet(rankScores = regScoreFET, tpNet = topNet5, 
                   nodeSizeScale = 2, cutoffFC = c(-0.3, 0.3), 
                   interactPlot = FALSE, fontSizeScale = 1, seed = 123)
rownames(regReg$topRankScores) = regReg$topRankScores$reg
topReg = sortDataframe(regReg$topRankScores, "score", decreasing = T)

####################### heatmap #############################
topScores = sortDataframe(regScoreFET$score, "score", decreasing = T)



################## Viper ##################
library(ARACNe)
netFile = "../rData/ARACNe_network_Lyme_GSE63085.txt"
aracneNet= read.csv(netFile, header = FALSE, sep = "\t", quote = "")
file = "../rData/regulatorOrderOfRegEnrichAndViper.Rdata"
load(file)
ARACNe_viper = as.character(ARACNe_viper)



##################### hubs #############################
# https://sna.stanford.edu/lab.php?l=4
# ####### hubness ##########
central_viper = hubness(aracneNet, reg = unique(as.character(aracneNet$V1)))
rownames(central_viper) = central_viper$reg
central_viper$viper = 0
central_viper[ARACNe_viper, "viper"] = 1
save(central_viper, file = "../rData/central_viper.rData")

# all degree and all closeness
g = topHubVenn(hubness = central_viper, colID = c(7,8, 9), topN = 50)
plotSave("../fig/Hubs_comparison_allDegree_allCloseness_viper_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)


# out degree and out closeness
g = topHubVenn(hubness = central_viper, colID = c(3,5, 9), topN = 50)
plotSave("../fig/Hubs_comparison_outDegree_outCloseness_viper_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)

# in degree and in closeness
g = topHubVenn(hubness = central_viper, colID = c(2,4, 9), topN = 50)
plotSave("../fig/Hubs_comparison_inDegree_inCloseness_viper_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)

# in degree, in closeness, Betweenness, and viper
g = topHubVenn(hubness = central_viper, colID = c(2,4,6, 9), topN = 50)


################# compare hubs between networks (ARACNE and WGCNA) #################
load("../rData/central_coeNetTop5.rData")
central_viper$reg = as.character(central_viper$reg)
central_coeNetTop5$reg = as.character(central_coeNetTop5$reg)
reg0 = unique(c(central_viper$reg, central_coeNetTop5$reg))

central = merge(central_viper, central_coeNetTop5, by = "reg", all = T)
colnames(central) = gsub("\\.x$", "_viper",
                         gsub("\\.y$", "_WGCNA", colnames(central)))
g = topHubVenn(hubness = central, colID = c(7,15), topN = 50)
plotSave("../fig/Hubs_comparison_allDegree_ARACNEvsWGCNA_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)

g = topHubVenn(hubness = central, colID = c(8,16), topN = 50)
plotSave("../fig/Hubs_comparison_allCloseness_ARACNEvsWGCNA_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)

g = topHubVenn(hubness = central, colID = c(9,17), topN = 50)
plotSave("../fig/Hubs_comparison_viperVsRegEnrichWGCNA_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)
