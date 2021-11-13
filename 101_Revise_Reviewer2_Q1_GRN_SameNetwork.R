
if(FALSE){
  options(width = 300)
  rm(list = ls())
  setwd("/hpc/dla_lti/wtao/3_RegEnrich_paper/script/")
  library(funcTools)
  rFile = "101_Revise_Reviewer2_Q1_GRN.R"
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = rFile, threaded = 10, 
         memoryG = 40, rTimeHour = 24, 
         logFile = paste0("../rData/101_Revise_Reviewer2_Q1_GRN.log"), 
         email = "weiyangtao1513@gmail.com", when2Email = "aes",
         preCMD = "Rscript ", 
         param1 = 1)
  qstat2("all")
  
}

library(RegEnrich)

{
  nTF = 500
  nTF_top = 20
  nSample = 100
  
  # Generating expression data and phenotype data
  load("../rData/100_expr_pData_simulate2.Rdata")
  TFs = grep("TF", rownames(expr), value = T)
  
  fileName = "../rData/101_object_rank_100samples_GRN_3.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  load("../rData/101_object_rank_2_GRN.Rdata")
  object100 = object
  network_user = results_topNet(object)
  save(object, file = fileName)
  }
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_100samples_GRN_3_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  set.seed(1234)
  object = regenrich_enrich(object, enrichTest = "GSEA")
  # Regulators ranking
  object = regenrich_rankScore(object)
  save(object, file = fileName)
  }
  
  ## 50 samples
  fileName = "../rData/101_object_rank_50samples_GRN_3.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  load("../rData/101_object_rank_50samples_2_GRN.Rdata")
  
  regenrich_network(object50) = network_user
  # Enrichment analysis by Fisher's exact test (FET)
  object50 = regenrich_enrich(object50)
  # Regulators ranking
  object50 = regenrich_rankScore(object50)
  grep("top", object50@resScore$reg)
  save(object50, file = fileName)
  }
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_50samples_GRN_3_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  set.seed(1234)
  object50 = regenrich_enrich(object50, enrichTest = "GSEA")
  # Regulators ranking
  object50 = regenrich_rankScore(object50)
  save(object50, file = fileName)
  }
  
  
  ## 20 samples
  fileName = "../rData/101_object_rank_20samples_GRN_3.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  load("../rData/101_object_rank_20samples_2_GRN.Rdata")
  
  regenrich_network(object20) = network_user
  # Enrichment analysis by Fisher's exact test (FET)
  object20 = regenrich_enrich(object20)
  # Regulators ranking
  object20 = regenrich_rankScore(object20)
  grep("top", object20@resScore$reg)
  save(object20, file = "../rData/101_object_rank_20samples_GRN_3.Rdata")
  }
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_20samples_GRN_3_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  set.seed(1234)
  object20 = regenrich_enrich(object20, enrichTest = "GSEA")
  # Regulators ranking
  object20 = regenrich_rankScore(object20)
  save(object20, file = fileName)
  }
  
  
  ## 10 samples
  fileName = "../rData/101_object_rank_10samples_GRN_3.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  id = c(1:5, seq(51, 55))
  design = model.matrix(~group, data = pData[id,])
  # Initializing a 'RegenrichSet' object
  object10 = RegenrichSet(expr = expr[, id],
                          colData = pData[id,],
                          method = 'limma', minMeanExpr = 0,
                          design = design,
                          contrast = c(0, 1), reg = TFs,
                          networkConstruction = 'GRN',
                          enrichTest = 'FET')
  object10 = regenrich_diffExpr(object10)
  
  regenrich_network(object10) = network_user
  # Enrichment analysis by Fisher's exact test (FET)
  object10 = regenrich_enrich(object10)
  # Regulators ranking
  object10 = regenrich_rankScore(object10)
  save(object10, file = fileName)
  }
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_10samples_GRN_3_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  set.seed(1234)
  object10 = regenrich_enrich(object10, enrichTest = "GSEA")
  # Regulators ranking
  object10 = regenrich_rankScore(object10)
  save(object10, file = fileName)
  }
  
  ## 6 samples
  fileName = "../rData/101_object_rank_6samples_GRN_3.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  id = c(1:3, seq(51, 53))
  design = model.matrix(~group, data = pData[id,])
  # Initializing a 'RegenrichSet' object
  object10 = RegenrichSet(expr = expr[, id],
                          colData = pData[id,],
                          method = 'limma', minMeanExpr = 0,
                          design = design,
                          contrast = c(0, 1), reg = TFs,
                          networkConstruction = 'GRN',
                          enrichTest = 'FET')
  object10 = regenrich_diffExpr(object10)
  
  regenrich_network(object10) = network_user
  # Enrichment analysis by Fisher's exact test (FET)
  object10 = regenrich_enrich(object10)
  # Regulators ranking
  object10 = regenrich_rankScore(object10)
  save(object10, file = fileName)
  }
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_6samples_GRN_3_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
  set.seed(1234)
  object6 = regenrich_enrich(object6, enrichTest = "GSEA")
  # Regulators ranking
  object6 = regenrich_rankScore(object6)
  save(object6, file = fileName)
  }
  
}



### Plot
fileNames = c("../rData/101_object_rank_100samples_GRN_3.Rdata", 
              "../rData/101_object_rank_50samples_GRN_3.Rdata", 
              "../rData/101_object_rank_20samples_GRN_3.Rdata", 
              "../rData/101_object_rank_10samples_GRN_3.Rdata", 
              "../rData/101_object_rank_6samples_GRN_3.Rdata",
              "../rData/101_object_rank_100samples_GRN_3_GSEA.Rdata", 
              "../rData/101_object_rank_50samples_GRN_3_GSEA.Rdata", 
              "../rData/101_object_rank_20samples_GRN_3_GSEA.Rdata", 
              "../rData/101_object_rank_10samples_GRN_3_GSEA.Rdata", 
              "../rData/101_object_rank_6samples_GRN_3_GSEA.Rdata")
library(ggplot2)
library(fgsea)
p = list()
for(mi in fileNames){
  data = get(load(mi))
  sudoPathway = list(topReg = grep("top", data@resScore$reg, value = T))
  sudoRanks = setNames((length(data@resScore$reg):1)/500 - 0.5, data@resScore$reg)
  enrich = fgsea(sudoPathway, sudoRanks)
  
  nn = substr(strSplit(mi, "\\.R")[,1], 21, 100)
  nn = gsub("sam", " sam", nn)
  p[[nn]] = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(strSplit(nn, "_")[,2], " (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
}

for(nn in names(p)){
  print(nn)
  figFile = paste0("../fig/101_keyTF_", nn, ".png")
  ggsave(filename = figFile, p[[nn]], width = 3, height = 1.8)
}


