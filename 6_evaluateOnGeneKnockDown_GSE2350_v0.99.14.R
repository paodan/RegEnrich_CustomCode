# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)
commandArgs()

source("0_functionsInThePaper.R")

### GPL8300 platform, [HG_U95Av2] Affymetrix Human Genome U95 Version 2 Array
# Total number of rows: 12625
#### Ref32_GSE2350_RAW_BCL6
geoAcession = "GSE2350"
if (file.exists("../rData/GSE2350.Rdata")){
  print("load data")
  load("../rData/GSE2350.Rdata")
} else {
  GSE2350 = getGEO("GSE2350")
  GSE2350 = GSE2350$`GSE2350-GPL8300_series_matrix.txt.gz`
  save(GSE2350, file = "../rData/GSE2350.Rdata")
}


gse = GSE2350[, tail(colnames(GSE2350), 8)] 
tmp = pData(gse)
tmp$group = factor(rep(c("NT", "BCL6"), each = 4), levels = c("NT", "BCL6"))
pData(gse) = tmp

phenoGSE2350 = pData(gse)
dataGSE2350 = log10(exprs(gse))
geneName = gse@featureData@data$`Gene Symbol`
geneNameValid = geneName[geneName != ""]
dataGSE2350 = dataGSE2350[geneName != "", ]
geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = " /// ")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE2350 = dataGSE2350[geneNameValidNoDupID,]
rownames(dataGSE2350) = geneNameValidNoMul[geneNameValidNoDupID]
dataGSE2350 = scale(dataGSE2350)
dataGSE2350TXT = data.frame(gene = rownames(dataGSE2350), dataGSE2350)
write.table(dataGSE2350TXT, file = "../rawData/dataGSE2350.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)



# viper method on data
netFile = "../rawData/ARACNe_network_GSE2350.txt"
if(!file.exists(netFile)){
  expr = "../rawData/dataGSE2350.txt"
  tfs = "../rData/tfs.txt"
  outputFolder = "../rData/GSE2350_ARACNe/"
  dir.create(outputFolder)
  pval = 10^-8
  aracne(expr = expr, tfs = tfs, 
         outputFolder = outputFolder, 
         threads = 1, calculateThreshold = T)
  
  nets = lapply(setNames(1:100, paste0("b", 1:100)), function(x) {
    aracne(expr = expr, tfs = tfs,pvalue = pval,
           outputFolder = outputFolder, seed = x)})
  net = aracne(outputFolder = outputFolder, consolidate = TRUE)
  if (ncol(net) == 4) net = net[, -4]
  colnames(net) = c("tf", "target", "mi")
  write.table(net, file = netFile, quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
}

eset = ExpressionSet(assayData=dataGSE2350, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE2350))

try(silent = TRUE, expr = {
  regul <- aracne2regulon(afile = netFile, eset = eset, 
                          format = "3col", verbose = FALSE)
  
  # # BCL6 vs NT
  # signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "Fox")
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "BCL6")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "Fox", per = 1000,
                         repos = TRUE, verbose = FALSE)
  
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  summary(mrs)
})



# # RegEnrich method on data
# # second batch
data(TFs)
regulators = unique(TFs$TF_name) 
design = model.matrix(~group, data = phenoGSE2350)

contrast = c(0, 1) # NT vs BCL6 

# Burkitt lymphoma cell line 
cell = "BurkittLymphoma"


library(BiocParallel)
bpparam = register(MulticoreParam(2))
object = RegenrichSet(expr = dataGSE2350, colData = phenoGSE2350, method = "limma", 
                      minMeanExpr = -20, 
                      design = design, contrast = contrast, 
                      reg = regulators, 
                      BPPARAM = bpparam,
                      networkConstruction = 'COEN',
                      softPower = NULL, RsquaredCut = 0.8,
                      maxSize = 15000,
                      enrichTest = 'FET')
print(dim(object))

# FET
object = object %>% regenrich_diffExpr() %>% 
  regenrich_network() %>% 
  regenrich_enrich() %>% 
  regenrich_rankScore()

score = results_score(object)
write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_BCL6_FET.tsv"), sep = "\t")
save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_BCL6_FET.Rdata"))

# GSEA
set.seed(1234)
object = object %>% regenrich_enrich(enrichTest = 'GSEA', nperm = 10000) %>% 
  regenrich_rankScore()

score = results_score(object)
write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_BCL6_GSEA.tsv"), sep = "\t")
save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_BCL6_GSEA.Rdata"))



if(FALSE){
  rQsub("./", "6_evaluateOnGeneKnockDown_GSE2350_v0.99.14.R", 
        jobName = "Job_06GSE2350", threaded = 2, memoryG = "30", 
        rTimeHour = 20, logFile = "logFile_06GSE2350.log", 
        email = "w.tao-2@umcutrecht.nl", 
        preCMD = "echo \"Rscript ", param1 = 1)
}
