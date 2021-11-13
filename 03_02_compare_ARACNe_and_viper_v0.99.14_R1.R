
# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)

library(funcTools)
library(RegEnrich)
library(microbenchmark)

library(doParallel)
library(foreach)
library(ggplot2)
library(DOSE)
library(WGCNA)
library(igraph)

library(viper)
library(bcellViper)


load("../rData/02_02_Lyme_GSE63085_full.Rdata")
log2FPKM_full = log2(Lyme_GSE63085_full$FPKM_full + 1)
sampleInfo_full = Lyme_GSE63085_full$sampleInfo_full
# sampleInfo$week = factor(sampleInfo$week, levels = c(0, 3))

data(TFs)

library(ARACNe)
expr = "../rData/03_02_log2FPKM_matrix.txt"
tfs = "../rData/03_02_tfs.txt"

log2FPKM1 = data.frame(gene = rownames(log2FPKM_full), log2FPKM_full)
write.table(log2FPKM1, file = expr, sep = "\t", 
            row.names = FALSE, quote = F)
write.table(unique(TFs[,2,drop = FALSE]), file = tfs, quote = F, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
netFile = "../rData/03_02_ARACNe_network_Lyme_GSE63085.txt"
if (!file.exists(netFile)){
  folder = "../rData/03_02_Lyme_GSE63085_ARACNe_outputFolder/"
  if(!dir.exists(folder)) dir.create(folder)
  # # Run on HPC
  # jobnames = paste0("ARA_", 1:100)
  # lognames = paste0("logARA_", 1:100)
  # rQsub(path = "/data/Github/RegEnrich_shared/script",
  #       rFile = "3_qsubArray.R", jobName = jobnames, threaded = 1,
  #       memoryG = 10, rTimeHour = 1, logFile = lognames,
  #       param1 = 1:100)
  
  t = aracne(expr = expr, tfs = tfs, outputFolder = folder,
             calculateThreshold = T)
  for(ar1 in 1:100){
    cat(ar1, "\n")
    write.table(ar1, file = "./03_02_log_aracne_bootstrap.log", append = T)
    aracne(expr = expr, tfs = tfs, outputFolder = folder, seed = ar1)
  }
  
  net = aracne(outputFolder = folder, consolidate = TRUE)
  if (ncol(net) == 4) net = net[, -4]
  colnames(net) = c("tf", "target", "mi")
  # netFile = "../rData/ARACNe_network_Lyme_GSE63085.txt"
  write.table(net, file = netFile, quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
}


# make ExpressionSet: eset
load("../rData/Lyme_GSE63085.rda")
exprs = as.matrix(read.csv("../rData/03_01_log2FPKM_matrix.txt", 
                           sep = "\t", quote = "", row.names = 1))
pData <- Lyme_GSE63085$sampleInfo
metadata <- data.frame(labelDescription=
                         c("GEO accession","Sample title",
                           "Sample source","Organism", "patient ID", 
                           "Week of treatment", "Disease state", 
                           "Description", "Molecule"),
                       row.names=colnames(pData))
phenoData <- new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)
eset = ExpressionSet(assayData=exprs, phenoData=phenoData)

# 
# 
# file = "../rData/03_02_mrs.Rdata"
# if(file.exists(file)){
#   load(file)
# } else {
#   regul <- aracne2regulon(afile = netFile, eset = eset, 
#                           format = "3col", verbose = FALSE)
#   
#   signature0 <- rowTtest(eset, pheno = "week", group1 = 0, group2 = 3)
#   
#   signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
#                   sign(signature0$statistic))[, 1]
#   nullmodel <- ttestNull(x = eset, pheno = "week", group1 = 0, 
#                          group2 = 3, per = 1000,
#                          repos = TRUE, verbose = FALSE)
#   
#   mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
#   mrs_summary = viper:::summary.msviper(mrs, 50) # summary(mrs)
#   save(mrs, mrs_summary, file = file)
# }



##### Using paired t test to perform vipe.
topReg50_viperFile = "../rData/03_02_3topReg50_viper.Rdata"
if (file.exists(topReg50_viperFile)){
  load(topReg50_viperFile)
} else {
  source("0_functionsInThePaper_2_testing_viper.R")
  group1 = 1:26
  group2 = 27:52
  signature1 = f1(exprs, group1, group2, paired = T)
  nullModelFile = "../rData/03_02_nullmodel1_compare_viper_pairedTtest.Rdata"
  if (file.exists(nullModelFile)){
    load(nullModelFile)
  } else {
    set.seed(1)
    nullmodel1 = f2(exprs, group1, group2, repos = T, per = 1000, paired = T)
    save(nullmodel1, file = nullModelFile)
  }
  # mrs1 <- msviper(signature, regul, nullmodel1, verbose = FALSE)
  mrs1 <- msviper(signature1, regul, nullmodel1, verbose = FALSE)
  topReg50_viper = sortDataframe(summary(mrs1, 50), "p.value")
  head(topReg50_viper)
  save(topReg50_viper, file = topReg50_viperFile)
}



source("0_functionsInThePaper.R")
load("../rData/pFC.Rdata")
load("../rData/log2FPKMhi.Rdata")
# regulatorKey = c("ZFP36L1", "ZFPM1", "CNOT7")
regulatorKey = c("ZFP36L1", "HIC1", "MAF1") # first 3 key regulators using the Revised ARACNe network
for (mi in regulatorKey){
  tmp2 = regul[[mi]]$likelihood
  tarLikelihood = setNames(tmp2, names(regul[[mi]]$tfmode))
  topNetTmp = list(net = list(names(tarLikelihood)))
  names(topNetTmp$net) = mi
  
  p = plotRegTarExpr(reg = mi, expr = log2FPKMhi, pFC = pFC,
                     topNet = topNetTmp, 
                     n = Inf,tarCol = alpha("black", 0.1),
                     xlab = "Samples",
                     ylab = "Z-scores", p_threshold = 1)
  # print(p)
  ggsave(paste0("../fig/03_02_viper_", mi, "_and_targes_expression.svg"), 
         plot = p, width = 7, height = 5, dpi = 300)
  ggsave(paste0("../fig/03_02_viper_", mi, "_and_targes_expression.png"), 
         plot = p, width = 7, height = 5, dpi = 300)
}




file = "../rData/03_02_regulatorOrderOfRegEnrichAndViper.Rdata"
if (file.exists(file)){
  load(file)
} else {
  load("../rData/02_01_object.Rdata")
  
  if(file.exists("../rData/03_01_object_COEN_GSEA_GRN_GSEA_GRN_FET.Rdata")){
    load("../rData/03_01_object_COEN_GSEA_GRN_GSEA_GRN_FET.Rdata")
  } else {
    # Enrichment analysis by GSEA for COEN network
    object_COEN_GSEA = regenrich_enrich(object, enrichTest = "GSEA")
    (object_COEN_GSEA = regenrich_rankScore(object_COEN_GSEA))
    
    # Enrichment analysis by GSEA for GRN network
    object_GRN_GSEA = regenrich_network(object, networkConstruction = "GRN")
    object_GRN_GSEA = regenrich_enrich(object_GRN_GSEA, enrichTest = "GSEA")
    (object_GRN_GSEA = regenrich_rankScore(object_GRN_GSEA))
    
    # Enrichment analysis by FET for GRN network
    object_GRN_FET = regenrich_enrich(object_GRN_GSEA, enrichTest = "FET")
    (object_GRN_FET = regenrich_rankScore(object_GRN_FET))
    
    
    save(object_COEN_GSEA, object_GRN_GSEA, object_GRN_FET, 
         file = "../rData/03_01_object_COEN_GSEA_GRN_GSEA_GRN_FET.Rdata")
  }
  
  topReg50_COEN_FET = data.frame(results_score(object))
  topReg50_COEN_SEA = data.frame(results_score(object_COEN_GSEA))
  topReg50_GRN_FET = data.frame(results_score(object_GRN_FET))
  topReg50_GRN_SEA = data.frame(results_score(object_GRN_GSEA))
  
  COEN_FET = topReg50_COEN_FET$reg[1:50]
  COEN_SEA = topReg50_COEN_SEA$reg[1:50]
  GRN_FET = topReg50_GRN_FET$reg[1:50]
  GRN_SEA = topReg50_GRN_SEA$reg[1:50]
  ARACNe_viper = topReg50_viper$Regulon
  
  save(COEN_FET, COEN_SEA, GRN_FET, GRN_SEA, ARACNe_viper, file = file)
}


### compare the orders of regulators obtained by different methods
top50Regs = data.frame(COEN_FET, COEN_SEA, GRN_FET, GRN_SEA, ARACNe_viper)
for(mi in 1:ncol(top50Regs)){
  for (ni in mi:ncol(top50Regs)){
    mtd = c(colnames(top50Regs[mi]), colnames(top50Regs[ni]))
    if(ni == mi) next() else {
      p = RegEnrich:::plotOrders(top50Regs[[mi]], top50Regs[[ni]]) +
        theme_Publication()+
        xlab(NULL)+ 
        scale_x_discrete(labels=c("top50Regs[[mi]]" = mtd[1], 
                                  "top50Regs[[ni]]" = mtd[2]))
      p
      for(figFile in c(paste0("../fig/03_02_compareRegOrder_", 
                              mtd[1], "_and_", mtd[2], ".svg"), 
                       paste0("../fig/03_02_compareRegOrder_", 
                              mtd[1], "_and_", mtd[2], ".png"))){
        print(figFile)
        ggsave(figFile, plot = p, width = 3, height = 4.5, dpi = 300)
      }
    }
  }
}


### plot venn diagram of common regulators by different methods.
library(VennDiagram)
mtd = colnames(top50Regs)
colors = setNames(c("#CDB4FF", "#A1D4FF", "#FFE5C3", 
                    "#FFA5C3", "#A5FFC3"), mtd)

cmp_groups = list(allMethod = mtd, 
                  RegEnirchMethod = mtd[1:4],
                  COEN_FET_viper = mtd[c(1,5)],
                  COEN_SEA_viper = mtd[c(2,5)],
                  GRN_FET_viper = mtd[c(3,5)],
                  GRN_SEA_viper = mtd[c(4,5)])
for(i in names(cmp_groups)){
  mi = cmp_groups[[i]]
  g = venn.diagram(x = as.list(top50Regs[mi]), 
                   filename = NULL, 
                   fill = colors[mi],
                   # cat.col = colors,
                   alpha = 1, cex = 1.8, 
                   cat.cex = 2, margin = 0.1)
  # grid.newpage(); grid.draw(g)
  plotSave(filename = paste0("../fig/03_02_vennDiagram_commonTop50Regs_", i, ".svg"),
           Plot = g, width = 6, height = 6, dpi = 300)
  plotSave(filename = paste0("../fig/03_02_vennDiagram_commonTop50Regs_", i, ".png"),
           Plot = g, width = 6, height = 6, dpi = 300)
}

