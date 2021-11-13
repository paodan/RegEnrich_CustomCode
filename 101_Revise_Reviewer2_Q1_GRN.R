
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
  fileName = "../rData/100_expr_pData_simulate.Rdata"
  load(fileName)
  
  
  library(BiocParallel)
  register(MulticoreParam(10)) # Number of cores
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    design = model.matrix(~group, data = pData)
    # Initializing a 'RegenrichSet' object
    object = RegenrichSet(expr = expr,
                          colData = pData,
                          method = 'limma', minMeanExpr = 0,
                          design = design,
                          contrast = c(0, 1), reg = TFs,
                          networkConstruction = 'GRN', 
                          enrichTest = 'FET')
    
    object = regenrich_diffExpr(object)
    
    # Network inference using 'GRN' method
    object = regenrich_network(object)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object = regenrich_enrich(object)
    
    # Regulators ranking
    object = regenrich_rankScore(object)
    save(object, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object@resScore)
  grep("top", object@resScore$reg)
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_50samples_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object50 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'GRN',
                            enrichTest = 'FET')
    
    object50 = regenrich_diffExpr(object50)
    
    # Network inference using 'COEN' method
    object50 = regenrich_network(object50)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object50 = regenrich_enrich(object50)
    
    # Regulators ranking
    object50 = regenrich_rankScore(object50)
    save(object50, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object50@resScore)
  grep("top", object50@resScore$reg)
  
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_20samples_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(10), seq(51, 60))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object20 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'GRN',
                            enrichTest = 'FET')
    
    object20 = regenrich_diffExpr(object20)
    
    # Network inference using 'COEN' method
    object20 = regenrich_network(object20)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object20 = regenrich_enrich(object20)
    
    # Regulators ranking
    object20 = regenrich_rankScore(object20)
    save(object20, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object20@resScore)
  grep("top", object20@resScore$reg)
  
  
  # # RegEnrich analysis
  # fileName = "../rData/101_object_rank_10samples_GRN.Rdata"
  # if(file.exists(fileName)){
  #   cat("loading", fileName, "\n")
  #   load(fileName)
  # } else {
  #   cat(fileName, "\n")
  #   TFs = grep("TF", rownames(expr), value = T)
  #   
  #   # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
  #   id = c(seq(5), seq(51, 55))
  #   design = model.matrix(~group, data = pData[id,])
  #   # Initializing a 'RegenrichSet' object
  #   object10 = RegenrichSet(expr = expr[, id],
  #                           colData = pData[id,],
  #                           method = 'limma', minMeanExpr = 0,
  #                           design = design,
  #                           contrast = c(0, 1), reg = TFs,
  #                           networkConstruction = 'GRN', softPower = 6,
  #                           enrichTest = 'FET')
  #   
  #   object10 = regenrich_diffExpr(object10)
  #   
  #   # Network inference using 'COEN' method
  #   object10 = regenrich_network(object10)
  #   
  #   # Enrichment analysis by Fisher's exact test (FET)
  #   object10 = regenrich_enrich(object10)
  #   
  #   # Regulators ranking
  #   object10 = regenrich_rankScore(object10)
  #   save(object10, file = fileName)
  # }
  # 
  # # Top Regs
  # as.data.frame(object10@resScore)
  # grep("top", object10@resScore$reg)
}

######## Change the fold change of TF and redo RegEnrich #########
{
  fileName = "../rData/100_expr_pData_simulate2.Rdata"
  load(fileName)
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_2_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    design = model.matrix(~group, data = pData)
    # Initializing a 'RegenrichSet' object
    object = RegenrichSet(expr = expr,
                          colData = pData,
                          method = 'limma', minMeanExpr = 0,
                          design = design,
                          contrast = c(0, 1), reg = TFs,
                          networkConstruction = 'GRN',
                          enrichTest = 'FET')
    
    object = regenrich_diffExpr(object)
    
    # Network inference using 'COEN' method
    object = regenrich_network(object)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object = regenrich_enrich(object)
    
    # Regulators ranking
    object = regenrich_rankScore(object)
    save(object, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object@resScore)
  grep("top", object@resScore$reg)
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_100samples_2_GRN_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
    set.seed(1234)
    object = regenrich_enrich(object, enrichTest = "GSEA")
    # Regulators ranking
    object = regenrich_rankScore(object)
    save(object, file = fileName)
  }
  
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_50samples_2_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object50 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'GRN',
                            enrichTest = 'FET')
    
    object50 = regenrich_diffExpr(object50)
    
    # Network inference using 'COEN' method
    object50 = regenrich_network(object50)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object50 = regenrich_enrich(object50)
    
    # Regulators ranking
    object50 = regenrich_rankScore(object50)
    save(object50, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object50@resScore)
  grep("top", object50@resScore$reg)
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_50samples_2_GRN_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
    set.seed(1234)
    object50 = regenrich_enrich(object50, enrichTest = "GSEA")
    # Regulators ranking
    object50 = regenrich_rankScore(object50)
    save(object50, file = fileName)
  }
  
  
  # RegEnrich analysis
  fileName = "../rData/101_object_rank_20samples_2_GRN.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(10), seq(51, 60))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object20 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'GRN',
                            enrichTest = 'FET')
    
    object20 = regenrich_diffExpr(object20)
    
    # Network inference using 'COEN' method
    object20 = regenrich_network(object20)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object20 = regenrich_enrich(object20)
    
    # Regulators ranking
    object20 = regenrich_rankScore(object20)
    save(object20, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object20@resScore)
  grep("top", object20@resScore$reg)
  
  # Enrichment analysis by GSEA
  fileName = "../rData/101_object_rank_20samples_2_GRN_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
    set.seed(1234)
    object20 = regenrich_enrich(object20, enrichTest = "GSEA")
    # Regulators ranking
    object20 = regenrich_rankScore(object20)
    save(object20, file = fileName)
  }

}





fileNames = c("../rData/101_object_rank_GRN.Rdata", 
              "../rData/101_object_rank_50samples_GRN.Rdata", 
              "../rData/101_object_rank_20samples_GRN.Rdata", 
              "../rData/101_object_rank_2_GRN.Rdata", 
              "../rData/101_object_rank_50samples_2_GRN.Rdata", 
              "../rData/101_object_rank_20samples_2_GRN.Rdata", 
              "../rData/101_object_rank_100samples_2_GRN_GSEA.Rdata", 
              "../rData/101_object_rank_50samples_2_GRN_GSEA.Rdata", 
              "../rData/101_object_rank_20samples_2_GRN_GSEA.Rdata")
library(ggplot2)
p = list()
for(mi in fileNames){
  data = get(load(mi))
  sudoPathway = list(topReg = grep("top", data@resScore$reg, value = T))
  sudoRanks = setNames((length(data@resScore$reg):1)/500 - 0.5, data@resScore$reg)
  enrich = fgsea(sudoPathway, sudoRanks)
  
  nn = substr(strSplit(mi, "\\.R")[,1], 21, 100)
  if(nn == "rank_GRN"){
    nn = "rank_100samples_GRN"
  } else if(nn == "rank_2_GRN"){
    nn = "rank_100samples_2_GRN"
  }
  nn = gsub("sam", " sam", nn)
  p[[nn]] = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(strSplit(nn, "_")[,2], " (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
}

for(nn in names(p)){
  print(nn)
  figFile = paste0("../fig/100_keyTF_", nn, ".png")
  ggsave(filename = figFile, p[[nn]], width = 3, height = 1.8)
}


