# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)

source("0_functionsInThePaper.R")


# GPL6947 platform, Illumina HumanHT-12 V3.0 expression beadchip
# Total number of rows: 49576
#### Ref18_GSE19114_series_matrix_STAT3_BHLHB2_FosL2_RunX1_CEBPb

geoAcession = "GSE19114"
if (file.exists("../rData/gseGSE19114.Rdata")){
  print("load data")
  load("../rData/gseGSE19114.Rdata")
} else {
  gseGSE19114 = getGEO("GSE19114")
  gseGSE19114 = gseGSE19114$GSE19114_series_matrix.txt.gz
  tmp = pData(gseGSE19114)
  tmp$cell = tmp$`cell type/line:ch1`
  kdGene = apply(strsplit2(tmp$`transduction:ch1`, " ")[, c(1, 6)], 
                 1, paste0, collapse = "_")
  tmp$group = gsub("NON-TARGET", "NT", toupper(sub("_$", "", kdGene)))
  pData(gseGSE19114) = tmp
  save(gseGSE19114, file = "../rData/gseGSE19114.Rdata")
}


## The third batch, the last 12 samples
gse = gseGSE19114[, tail(colnames(gseGSE19114), 12)] 
phenoGSE19114 = pData(gse)

dataGSE19114 = exprs(gse)
geneName = gse@featureData@data$ILMN_Gene
geneNameValid = geneName[geneName != ""]
dataGSE19114 = dataGSE19114[geneName != "", ]

geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = "/")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE19114 = dataGSE19114[geneNameValidNoDupID,]
rownames(dataGSE19114) = geneNameValidNoMul[geneNameValidNoDupID]
# dataGSE19114 = scale(dataGSE19114)
dataGSE19114TXT = data.frame(gene = rownames(dataGSE19114), dataGSE19114)
write.table(dataGSE19114TXT, file = "../rawData/dataGSE19114_3.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# viper method on data
netFile = "../rawData/ARACNe_network_GSE19114_3.txt"
if (!file.exists(netFile)){
  expr = "../rawData/dataGSE19114_3.txt"
  tfs = "../rData/tfs.txt"
  outputFolder = "../rData/GSE19114_ARACNe_3/"
  dir.create(outputFolder)
  pval = 10^-8
  aracne(expr = expr, tfs = tfs, pvalue = pval,
         outputFolder = outputFolder, 
         threads = 1, calculateThreshold = T)
  
  n = 100
  nets = lapply(setNames(1:n, paste0("b", 1:n)), function(x) {
    aracne(expr = expr, tfs = tfs, pvalue = pval, threads = 10,
           outputFolder = outputFolder, seed = x)})
  net = aracne(outputFolder = outputFolder, consolidate = TRUE)
  if (ncol(net) == 4) net = net[, -4]
  colnames(net) = c("tf", "target", "mi")
  write.table(net, file = netFile, quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
}


eset = ExpressionSet(assayData=dataGSE19114, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE19114))
regul <- aracne2regulon(afile = netFile, eset = eset, 
                        format = "3col", verbose = FALSE)

mrsFile = "../rData/mrs_GSE19114_ARACNe_3_STAT3.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # "STAT3" knockout (not found)
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "STAT3", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(mrs$regulon))
sortDataframe(res, "p.value")



mrsFile = "../rData/mrs_GSE19114_ARACNe_3_CEBPB.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # "CEBPB" knockout (not found)
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "CEBPB")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "CEBPB", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(mrs$regulon))
sortDataframe(res, "p.value")


mrsFile = "../rData/mrs_GSE19114_ARACNe_3_STAT3_CEBPB.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # "STAT3_CEBPB" knockout (not found)
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3_CEBPB")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "STAT3_CEBPB", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(mrs$regulon))
sortDataframe(res, "p.value")



# # RegEnrich method on data
cell = "SNB19"
# # third batch
data(TFs)
regulators = unique(TFs$TF_name) 
phenoGSE19114$group = factor(phenoGSE19114$group,
                             levels = c("NT", "STAT3", "CEBPB", "STAT3_CEBPB")) # third batch
design = model.matrix(~0 + group, data = phenoGSE19114)

contrast = list(STAT3 = c(-1, 1, 0, 0), # NT vs STAT3 (ranked No 14 *****)
                CEBPB = c(-1, 0, 1, 0), # NT vs CEBPB
                STAT3_CEBPB = c(-1, 0, 0, 1)) # NT vs STAT3_CEBPB

library(BiocParallel)
bpparam = register(MulticoreParam(2))
object = RegenrichSet(expr = dataGSE19114, colData = phenoGSE19114, method = "limma", 
                      minMeanExpr = -20, 
                      design = design, contrast = contrast[[1]], 
                      reg = regulators, 
                      BPPARAM = bpparam,
                      networkConstruction = 'COEN',
                      softPower = 2, 
                      maxSize = 15000,
                      enrichTest = 'FET')
print(dim(object))

# FET
object = object %>% regenrich_diffExpr() %>% 
  regenrich_network() 
network_user = results_topNet(object)

for(gene in c("STAT3", "CEBPB", "STAT3_CEBPB")){
  object = object %>% regenrich_diffExpr(contrast = contrast[[gene]])
  regenrich_network(object) = network_user
  # FET
  object = object %>% regenrich_enrich(enrichTest = 'FET') %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_thirdBatch_cell_", cell, "_", gene, "_FET.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_thirdBatch_cell_", cell, "_", gene, "_FET.Rdata"))
  
  # GSEA
  set.seed(1234)
  object = object %>% regenrich_enrich(enrichTest = 'GSEA', nperm = 10000) %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_thirdBatch_cell_", cell, "_", gene, "_GSEA.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_thirdBatch_cell_", cell, "_", gene, "_GSEA.Rdata"))
}


if(FALSE){
  rQsub("./", "6_evaluateOnGeneKnockDown_GSE19114_thirdBatch_v0.99.14.R", 
        jobName = "Job_06_GSE19114_thirdBatch", threaded = 2, memoryG = "60", 
        rTimeHour = 40, logFile = "logFile_06_GSE19114_thirdBatch.log", 
        email = "w.tao-2@umcutrecht.nl", 
        preCMD = "echo \"Rscript ", param1 = 1)
}





