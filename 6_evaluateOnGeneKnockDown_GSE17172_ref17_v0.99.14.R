# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)
commandArgs()

source("0_functionsInThePaper.R")

# GPL8300 platform, [HG_U95Av2] Affymetrix Human Genome U95 Version 2 Array
# Total number of rows: 12625
#### Ref17_GSE17172_RAW_FOXM1_MYB
geoAcession = "GSE17172"
if (file.exists("../rData/gse_GSE17172.Rdata")){
  print("load data")
  load("../rData/gse_GSE17172.Rdata")
} else {
  gse = getGEO(filename = "../rawData/GeneKnockDownDataset/Ref17_GSE17172_series_matrix_FOXM1_MYB.txt")
  save(gse, file = "../rData/gse_GSE17172.Rdata")
}
phenoGSE17172 = pData(gse)

tmp = strsplit2(phenoGSE17172$title, "_")[,2]
phenoGSE17172$group = sub("\\d$", "", strsplit2(tmp, " ")[,1])
phenoGSE17172$normMethod = gsub("[()]", "", strsplit2(tmp, " ")[,2])
phenoGSE17172 = subset(phenoGSE17172,normMethod == "MAS5")

# Only the data of MAS5 method
idMAS5 = rownames(subset(phenoGSE17172,normMethod == "MAS5"))
gseMAS5 = gse[, idMAS5]
dataGSE17172 = exprs(gseMAS5)
geneName = gseMAS5@featureData@data$`Gene Symbol`
geneNameValid = geneName[geneName != ""]
dataGSE17172 = dataGSE17172[geneName != "", ]
geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = " /// ")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE17172 = dataGSE17172[geneNameValidNoDupID,]
rownames(dataGSE17172) = geneNameValidNoMul[geneNameValidNoDupID]
dataGSE17172 = scale(dataGSE17172)
dataGSE17172TXT = data.frame(gene = rownames(dataGSE17172), dataGSE17172)
write.table(dataGSE17172TXT, file = "../rawData/dataGSE17172.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# viper method on data
netFile = "../rawData/ARACNe_network_GSE17172.txt"
if(!file.exists(netFile)){
  expr = "../rawData/dataGSE17172.txt"
  tfs = "../rData/tfs.txt"
  outputFolder = "../rData/GSE17172_ARACNe/"
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

eset = ExpressionSet(assayData=dataGSE17172, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE17172))

try(silent = TRUE, expr = {
  regul <- aracne2regulon(afile = netFile, eset = eset, 
                          format = "3col", verbose = FALSE)
  
  # FOX vs NT
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "Fox")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "Fox", per = 1000,
                         repos = TRUE, verbose = FALSE)
  
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  summary(mrs)
  save(mrs, file = paste0("../rData/mrs_", geoAcession, "FOX.Rdata"))
  
})



# RegEnrich method on data
cell = "ST486"
gene = "FOXM1" #"Fox"
data(TFs)
regulators = unique(TFs$TF_name)
phenoGSE17172$group = factor(phenoGSE17172$group, 
                             levels = c("NT", "Fox", "MYB"))
design = model.matrix(~1 + group, data = phenoGSE17172)
# contrast = c(-1, 0, 1) # NT vs MYB
contrast = c(0, 1, 0) # NT vs Fox (******* This works ******* !!!!!!)


library(BiocParallel)
bpparam = register(MulticoreParam(2))
object = RegenrichSet(expr = dataGSE17172, colData = phenoGSE17172, method = "limma", 
                      minMeanExpr = -20, 
                      design = design, contrast = contrast, 
                      reg = regulators, 
                      BPPARAM = bpparam,
                      networkConstruction = 'COEN',
                      softPower = NULL, #RsquaredCut = 0.8,
                      maxSize = 15000,
                      enrichTest = 'FET')
print(dim(object))

# FET
object = object %>% regenrich_diffExpr() %>% 
  regenrich_network() %>% 
  regenrich_enrich() %>% 
  regenrich_rankScore()

score = results_score(object)
write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", gene, "_FET.tsv"), sep = "\t")
save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", gene, "_FET.Rdata"))

# GSEA
set.seed(1234)
object = object %>% regenrich_enrich(enrichTest = 'GSEA', nperm = 10000) %>% 
  regenrich_rankScore()

score = results_score(object)
write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", gene, "_GSEA.tsv"), sep = "\t")
save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", gene, "_GSEA.Rdata"))

network_user = results_topNet(object)



if(FALSE){
  rQsub("./", "6_evaluateOnGeneKnockDown_GSE17172_ref17_v0.99.14.R", 
        jobName = "Job_06GSE17172", threaded = 2, memoryG = "60", 
        rTimeHour = 20, logFile = "logFile_06GSE17172.log", 
        email = "w.tao-2@umcutrecht.nl", 
        preCMD = "echo \"Rscript ", param1 = 1)
}
