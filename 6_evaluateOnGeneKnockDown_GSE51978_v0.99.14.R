# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)
commandArgs()

source("0_functionsInThePaper.R")

# GPL570 platform, [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# Total number of rows: 54675
### GSE51978, CHAF1A knockdown, two time points (day 5 and 10)
geoAcession = "GSE51978"
if (file.exists("../rData/GSE51978.Rdata")){
  print("load data")
  load("../rData/GSE51978.Rdata")
} else {
  GSE51978 = getGEO("GSE51978")
  GSE51978 = GSE51978$GSE51978_series_matrix.txt.gz
  
  tmp = pData(GSE51978)
  tmp$cell = "neuroblastoma"
  tmp$group = factor(rep(c("D00", "D05", "D10"), each = 3), 
                     levels = c("D00", "D05", "D10"))
  pData(GSE51978) = tmp
  save(GSE51978, file = "../rData/GSE51978.Rdata")
}

gse = GSE51978

phenoGSE51978 = pData(gse)
dataGSE51978 = log10(exprs(gse))
geneName = gse@featureData@data$`Gene Symbol`
geneNameValid = geneName[geneName != ""]
dataGSE51978 = dataGSE51978[geneName != "", ]

geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = "/")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE51978 = dataGSE51978[geneNameValidNoDupID,]
rownames(dataGSE51978) = geneNameValidNoMul[geneNameValidNoDupID]
dataGSE51978TXT = data.frame(gene = rownames(dataGSE51978), dataGSE51978)
write.table(dataGSE51978TXT, file = "../rawData/dataGSE51978.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)


# viper method on data (Failed to construct network)
netFile = "../rawData/ARACNe_network_GSE51978.txt"
if (!file.exists(netFile)){
  expr = "../rawData/dataGSE51978.txt"
  tfs = "../rData/tfs.txt"
  outputFolder = "../rData/GSE51978_ARACNe/"
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

eset = ExpressionSet(assayData=dataGSE51978, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE51978))


# netFile is empty *****
try(silent = TRUE, expr = {
  regul <- aracne2regulon(afile = netFile, eset = eset, 
                          format = "3col", verbose = FALSE)
  
  mrsFile = "../rData/mrs_GSE51978_ARACNe_CHAF1A.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "CHAF1A" knockout (Failed)
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "CHAF1A")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                           group2 = "CHAF1A", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(mrs$regulon))
  sortDataframe(res, "p.value")
})




### RegEnrich
library(RegEnrich)
data(TFs)
regulators = unique(TFs$TF_name)
design = model.matrix(~group, data = phenoGSE51978)

cell = "IMR32"


##--------Day 5--------
{
  time = "day5"
  contrast = c(0, 1, 0) # NT vs CHAF1A on day 5
  
  library(BiocParallel)
  bpparam = register(MulticoreParam(2))
  object = RegenrichSet(expr = dataGSE51978, colData = phenoGSE51978, method = "limma", 
                        minMeanExpr = -20, 
                        design = design, contrast = contrast, 
                        reg = regulators, 
                        BPPARAM = bpparam,
                        networkConstruction = 'COEN',
                        softPower = NULL,
                        maxSize = 15000,
                        enrichTest = 'FET')
  print(dim(object))
  
  # FET
  object = object %>% regenrich_diffExpr() %>% 
    regenrich_network() %>% 
    regenrich_enrich() %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_FET.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_FET.Rdata"))
  
  # GSEA
  set.seed(1234)
  object = object %>% regenrich_enrich(enrichTest = 'GSEA', nperm = 10000) %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_GSEA.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_GSEA.Rdata"))
}
network_user = results_topNet(object)


##--------Day 10--------
{
  time = "day10"
  contrast = c(0, 0, 1) # NT vs CHAF1A on day 10
  
  object = object %>% regenrich_diffExpr(contrast = contrast) 
  
  # FET
  regenrich_network(object) = network_user
  object = object %>% regenrich_enrich(enrichTest = 'FET') %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_FET.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_FET.Rdata"))
  
  
  # GSEA
  set.seed(1234)
  object = object %>% regenrich_enrich(enrichTest = 'GSEA', nperm = 10000) %>% 
    regenrich_rankScore()
  
  score = results_score(object)
  write.table(score, paste0("../rData/score_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_GSEA.tsv"), sep = "\t")
  save(object, file = paste0("../rData/object_06_", geoAcession, "_cell_", cell, "_", time, "_CHAF1A_GSEA.Rdata"))
}



if(FALSE){
  rQsub("./", "6_evaluateOnGeneKnockDown_GSE51978_v0.99.14.R", 
        jobName = "Job_06GSE51978", threaded = 2, memoryG = "60", 
        rTimeHour = 20, logFile = "logFile_06GSE51978.log", 
        email = "w.tao-2@umcutrecht.nl", 
        preCMD = "echo \"Rscript ", param1 = 1)
}
