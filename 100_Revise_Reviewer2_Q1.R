
if(FALSE){
  options(width = 300)
  rm(list = ls())
  setwd("/hpc/dla_lti/wtao/3_RegEnrich_paper/script/")
  library(funcTools)
  rFile = "100_Revise_Reviewer2_Q1.R"
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = rFile, threaded = 1, 
         memoryG = 40, rTimeHour = 24, 
         logFile = paste0("../rData/100_Revise_Reviewer2_Q1.log"), 
         email = "weiyangtao1513@gmail.com", when2Email = "aes",
         preCMD = "Rscript ", 
         param1 = 1)
  qstat2("all")
  
}

library(RegEnrich)

{
  # Generating expression data and phenotype data
  nTF = 500
  nTF_top = 20
  nSample = 100
  fileName = "../rData/100_expr_pData_simulate.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    
    # Top TF (differential TF)
    set.seed(1234) 
    expr_TF_top = do.call("rbind", lapply(seq(nTF_top), function(x){
      geneName = paste0('G', x)
      y1 = rnorm(n = nSample/2, mean = 10, sd = 5)
      y2 = rnorm(n = nSample/2, mean = 11 + x/10, sd = 5)
      matrix(c(y1, y2), nrow = 1, 
             dimnames = list(geneName, paste0('Samp', seq(nSample))))
    }))
    
    # Other TF  (non-differential TF)
    expr_TF_rest = do.call("rbind", lapply(seq(nTF_top+1, nTF), function(x){
      geneName = paste0('G', x)
      y1 = rnorm(n = nSample, mean = 15, sd = 5)
      matrix(y1, nrow = 1, 
             dimnames = list(geneName, paste0('Samp', seq(nSample))))
    }))
    
    # Targets of top TF (high-correlation)
    expr_Tar_top = do.call("rbind", apply(expr_TF_top, 1, function(x){
      nTar = sample(3:50, size = 1)
      do.call("rbind", lapply(1:nTar, function(y){
        x + rnorm(length(x), mean = mean(x)/5, sd = 3)
      }))
    }))
    
    # Targets of Other TF (low-correlation)
    expr_Tar_rest = do.call("rbind", apply(expr_TF_rest, 1, function(x){
      nTar = sample(3:50, size = 1)
      do.call("rbind", lapply(1:nTar, function(y){
        x + rnorm(length(x), mean = mean(x)/5, sd = 6)
      }))
    }))
    
    expr = rbind(expr_TF_top, expr_TF_rest, expr_Tar_top, expr_Tar_rest)
    rownames(expr) = paste0(rep(c("TF_top", "TF", "Tar_top", "Tar"), 
                                c(nTF_top, nTF-nTF_top, nrow(expr_Tar_top), nrow(expr_Tar_rest))), 
                            seq(nrow(expr)))
    
    pData = data.frame(sampleName = colnames(expr), 
                       group = rep(c("A", "B"), each = nSample/2), 
                       row.names = colnames(expr))
    
    save(expr, pData, file = fileName)
  }
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank.Rdata"
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
                          networkConstruction = 'COEN',
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
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_50samples.Rdata"
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
                            networkConstruction = 'COEN',
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
  fileName = "../rData/100_object_rank_20samples.Rdata"
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
                            networkConstruction = 'COEN',
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
  
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_10samples.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(5), seq(51, 55))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object10 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'COEN', softPower = 6,
                            enrichTest = 'FET')
    
    object10 = regenrich_diffExpr(object10)
    
    # Network inference using 'COEN' method
    object10 = regenrich_network(object10)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object10 = regenrich_enrich(object10)
    
    # Regulators ranking
    object10 = regenrich_rankScore(object10)
    save(object10, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object10@resScore)
  grep("top", object10@resScore$reg)
  
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_6samples.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(3), seq(51, 53))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object6 = RegenrichSet(expr = expr[, id],
                           colData = pData[id,],
                           method = 'limma', minMeanExpr = 0,
                           design = design,
                           contrast = c(0, 1), reg = TFs,
                           networkConstruction = 'COEN', softPower = 6,
                           enrichTest = 'FET')
    
    object6 = regenrich_diffExpr(object6)
    
    # Network inference using 'COEN' method
    object6 = regenrich_network(object6)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object6 = regenrich_enrich(object6)
    
    # Regulators ranking
    object6 = regenrich_rankScore(object6)
    save(object6, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object6@resScore)
  grep("top", object6@resScore$reg)
}

######## Change the fold change of TF and redo RegEnrich #########  ******* To be finished
{
  nTF = 500
  nTF_top = 20
  nSample = 100
  
  fileName = "../rData/100_expr_pData_simulate2.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    
    # Top TF (differential TF)
    set.seed(1234) 
    expr_TF_top = do.call("rbind", lapply(seq(nTF_top), function(x){
      geneName = paste0('G', x)
      y1 = rnorm(n = nSample/2, mean = 10, sd = 5)
      y2 = rnorm(n = nSample/2, mean = 13 + x/10, sd = 5)
      matrix(c(y1, y2), nrow = 1, 
             dimnames = list(geneName, paste0('Samp', seq(nSample))))
    }))
    
    # Other TF  (non-differential TF)
    expr_TF_rest = do.call("rbind", lapply(seq(nTF_top+1, nTF), function(x){
      geneName = paste0('G', x)
      y1 = rnorm(n = nSample, mean = 15, sd = 5)
      matrix(y1, nrow = 1, 
             dimnames = list(geneName, paste0('Samp', seq(nSample))))
    }))
    
    # Targets of top TF (high-correlation)
    expr_Tar_top = do.call("rbind", apply(expr_TF_top, 1, function(x){
      nTar = sample(3:50, size = 1)
      do.call("rbind", lapply(1:nTar, function(y){
        x + rnorm(length(x), mean = mean(x)/5, sd = 3)
      }))
    }))
    
    # Targets of Other TF (low-correlation)
    expr_Tar_rest = do.call("rbind", apply(expr_TF_rest, 1, function(x){
      nTar = sample(3:50, size = 1)
      do.call("rbind", lapply(1:nTar, function(y){
        x + rnorm(length(x), mean = mean(x)/5, sd = 6)
      }))
    }))
    
    expr = rbind(expr_TF_top, expr_TF_rest, expr_Tar_top, expr_Tar_rest)
    rownames(expr) = paste0(rep(c("TF_top", "TF", "Tar_top", "Tar"), 
                                c(nTF_top, nTF-nTF_top, nrow(expr_Tar_top), nrow(expr_Tar_rest))), 
                            seq(nrow(expr)))
    
    pData = data.frame(sampleName = colnames(expr), 
                       group = rep(c("A", "B"), each = nSample/2), 
                       row.names = colnames(expr))
    
    save(expr, pData, file = fileName)
  }
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_2.Rdata"
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
                          networkConstruction = 'COEN',
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
  fileName = "../rData/100_object_rank_100samples_2_GSEA.Rdata"
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
  fileName = "../rData/100_object_rank_50samples_2.Rdata"
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
                            networkConstruction = 'COEN',
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
  fileName = "../rData/100_object_rank_50samples_2_GSEA.Rdata"
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
  fileName = "../rData/100_object_rank_20samples_2.Rdata"
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
                            networkConstruction = 'COEN',
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
  fileName = "../rData/100_object_rank_20samples_2_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
    set.seed(1234)
    object20 = regenrich_enrich(object20, enrichTest = "GSEA")
    # Regulators ranking
    object20 = regenrich_rankScore(object20)
    save(object20, file = fileName)
  }
  
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_10samples_2.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(5), seq(51, 55))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object10 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'COEN', softPower = 6,
                            enrichTest = 'FET')
    
    object10 = regenrich_diffExpr(object10)
    
    # Network inference using 'COEN' method
    object10 = regenrich_network(object10)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object10 = regenrich_enrich(object10)
    
    # Regulators ranking
    object10 = regenrich_rankScore(object10)
    save(object10, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object10@resScore)
  grep("top", object10@resScore$reg)
  
  # Enrichment analysis by GSEA
  fileName = "../rData/100_object_rank_10samples_2_GSEA.Rdata"
  if(file.exists(fileName)){
    load(fileName)
  } else {
    set.seed(1234)
    object10 = regenrich_enrich(object10, enrichTest = "GSEA")
    # Regulators ranking
    object10 = regenrich_rankScore(object10)
    save(object10, file = fileName)
  }
  
  
  # RegEnrich analysis
  fileName = "../rData/100_object_rank_6samples_2.Rdata"
  if(file.exists(fileName)){
    cat("loading", fileName, "\n")
    load(fileName)
  } else {
    cat(fileName, "\n")
    TFs = grep("TF", rownames(expr), value = T)
    
    # id = c(seq(nSample/4), seq(nSample/2 + 1, nSample*3/4))
    id = c(seq(3), seq(51, 53))
    design = model.matrix(~group, data = pData[id,])
    # Initializing a 'RegenrichSet' object
    object6 = RegenrichSet(expr = expr[, id],
                            colData = pData[id,],
                            method = 'limma', minMeanExpr = 0,
                            design = design,
                            contrast = c(0, 1), reg = TFs,
                            networkConstruction = 'COEN', softPower = 6,
                            enrichTest = 'FET')
    
    object6 = regenrich_diffExpr(object6)
    
    # Network inference using 'COEN' method
    object6 = regenrich_network(object6)
    
    # Enrichment analysis by Fisher's exact test (FET)
    object6 = regenrich_enrich(object6)
    
    # Regulators ranking
    object6 = regenrich_rankScore(object6)
    save(object6, file = fileName)
  }
  
  # Top Regs
  as.data.frame(object6@resScore)
  grep("top", object6@resScore$reg)
  
  # Enrichment analysis by GSEA
  fileName = "../rData/100_object_rank_6samples_2_GSEA.Rdata"
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



fileNames = c("../rData/100_object_rank.Rdata", 
              "../rData/100_object_rank_50samples.Rdata", 
              "../rData/100_object_rank_20samples.Rdata", 
              "../rData/100_object_rank_10samples.Rdata", 
              "../rData/100_object_rank_6samples.Rdata", 
              "../rData/100_object_rank_2.Rdata", 
              "../rData/100_object_rank_50samples_2.Rdata", 
              "../rData/100_object_rank_20samples_2.Rdata", 
              "../rData/100_object_rank_10samples_2.Rdata",
              "../rData/100_object_rank_6samples_2.Rdata", 
              "../rData/100_object_rank_100samples_2_GSEA.Rdata", 
              "../rData/100_object_rank_50samples_2_GSEA.Rdata", 
              "../rData/100_object_rank_20samples_2_GSEA.Rdata", 
              "../rData/100_object_rank_10samples_2_GSEA.Rdata",
              "../rData/100_object_rank_6samples_2_GSEA.Rdata")
library(ggplot2)
library(fgsea)
p = list()
for(mi in fileNames){
  data = get(load(mi))
  sudoPathway = list(topReg = grep("top", data@resScore$reg, value = T))
  sudoRanks = setNames((length(data@resScore$reg):1)/500 - 0.5, data@resScore$reg)
  enrich = fgsea(sudoPathway, sudoRanks)
  
  nn = substr(strSplit(mi, "\\.R")[,1], 21, 100)
  if(nn == "rank"){
    nn = "rank_100samples"
  } else if(nn == "rank_2"){
    nn = "rank_100samples_2"
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

