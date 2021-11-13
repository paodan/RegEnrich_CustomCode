
if(FALSE){
  options(width = 300)
  rm(list = ls())
  setwd("/hpc/dla_lti/wtao/3_RegEnrich_paper/script/")
  library(funcTools)
  rFile = "102_Revise_Reviewer2_Q1_viper.R"
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = rFile, threaded = 1, 
         memoryG = 40, rTimeHour = 24, 
         logFile = paste0("../rData/102_Revise_Reviewer2_Q1_viper.log"), 
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
  
  
  library(ARACNe)
  exprFile = "../rData/102_expr.txt"
  tfs = "../rData/102_tfs.txt"
  
  log2FPKM1 = data.frame(gene = rownames(expr), expr)
  write.table(log2FPKM1, file = exprFile, sep = "\t", 
              row.names = FALSE, quote = F)
  write.table(grep("TF", rownames(expr), value = T), file = tfs, quote = F, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  netFile = "../rData/102_ARACNe_network_simulate.txt"
  if (!file.exists(netFile)){
    set.seed(123)
    folder = "../rData/102_ARACNe_network_simulate_outputFolder/"
    if(!dir.exists(folder)) dir.create(folder)
    
    t = aracne(expr = exprFile, tfs = tfs, outputFolder = folder,
               calculateThreshold = T)
    for(ar1 in 1:100){
      cat(ar1, "\n")
      write.table(ar1, file = "./102_log_aracne_bootstrap.log", append = T)
      aracne(expr = exprFile, tfs = tfs, outputFolder = folder, seed = ar1)
    }
    
    net = aracne(outputFolder = folder, consolidate = TRUE)
    if (ncol(net) == 4) net = net[, -4]
    colnames(net) = c("tf", "target", "mi")
    # netFile = "../rData/ARACNe_network_Lyme_GSE63085.txt"
    write.table(net, file = netFile, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  phenoData <- new("AnnotatedDataFrame", data=pData)
  eset = ExpressionSet(assayData = expr, phenoData=phenoData)
  
  library(viper)
  
  file = "../rData/102_mrs.Rdata"
  if(file.exists(file)){
    load(file)
  } else {
    regul <- aracne2regulon(afile = netFile, eset = eset,
                            format = "3col", verbose = FALSE)
    
    # signature0 <- rowTtest(eset, pheno = "group", group1 = "A", group2 = "B")
    signature0 <- rowTtest(eset, pheno = "group", group1 = "B", group2 = "A")

    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    set.seed(123)
    # nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "A",
    #                        group2 = "B", per = 1000,
    #                        repos = TRUE, verbose = FALSE)
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "B",
                           group2 = "A", per = 1000,
                           repos = TRUE, verbose = FALSE)

    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    mrs_summary = viper:::summary.msviper(mrs, 50) # summary(mrs)
    save(mrs, mrs_summary, regul, file = file)
  }
  

  ## Down-sampling
  # regul <- aracne2regulon(afile = netFile, eset = eset,
  #                         format = "3col", verbose = FALSE)
  for(mi in c(100, 50, 20, 10, 6)){
    file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), ".Rdata")
    if(file.exists(file)){
      load(file)
    } else {
      id = c(1:(mi/2), 51:(mi/2 + 50))
      # regul <- aracne2regulon(afile = netFile, eset = eset[,id],
      #                         format = "3col", verbose = FALSE)
      
      # signature0 <- rowTtest(eset[,id], pheno = "group", group1 = "A", group2 = "B")
      signature0 <- rowTtest(eset[,id], pheno = "group", group1 = "B", group2 = "A")
      
      signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                      sign(signature0$statistic))[, 1]
      set.seed(123)
      # nullmodel <- ttestNull(x = eset[,id], pheno = "group", group1 = "A",
      #                        group2 = "B", per = 1000,
      #                        repos = TRUE, verbose = FALSE)
      nullmodel <- ttestNull(x = eset[,id], pheno = "group", group1 = "B",
                             group2 = "A", per = 1000,
                             repos = TRUE, verbose = FALSE)
      
      mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
      mrs_summary = viper:::summary.msviper(mrs, 500) # summary(mrs)
      save(mrs, mrs_summary, file = file)
    }
  }
}


## Plot 
library(fgsea)
library(ggplot2)
for(mi in c(100, 50, 20, 10, 6)){
  file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), ".Rdata")
  load(file)
  mrs_summaryAll = viper:::summary.msviper(mrs, 500) 
  
  sudoPathway = list(topReg = grep("top", mrs_summaryAll$Regulon, value = T))
  sudoRanks = setNames((length(mrs_summaryAll$Regulon):1)/500 - 0.5, mrs_summaryAll$Regulon)
  enrich = fgsea(sudoPathway, sudoRanks)
  p = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(mi, " samples (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  figFile = paste0("../fig/102_keyTF_", mi, "Samples_viper.png")
  ggsave(filename = figFile, p, width = 3, height = 1.8)
}


{
  nTF = 500
  nTF_top = 20
  nSample = 100
  
  # Generating expression data and phenotype data
  fileName = "../rData/100_expr_pData_simulate2.Rdata"
  load(fileName)
  
  
  library(ARACNe)
  exprFile = "../rData/102_expr_2.txt"
  tfs = "../rData/102_tfs.txt"
  
  log2FPKM1 = data.frame(gene = rownames(expr), expr)
  write.table(log2FPKM1, file = exprFile, sep = "\t", 
              row.names = FALSE, quote = F)
  write.table(grep("TF", rownames(expr), value = T), file = tfs, quote = F, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  netFile = "../rData/102_ARACNe_network_simulate_2.txt"
  if (!file.exists(netFile)){
    set.seed(123)
    folder = "../rData/102_ARACNe_network_simulate_outputFolder_2/"
    if(!dir.exists(folder)) dir.create(folder)
    
    t = aracne(expr = exprFile, tfs = tfs, outputFolder = folder,
               calculateThreshold = T)
    for(ar1 in 1:100){
      cat(ar1, "\n")
      write.table(ar1, file = "./102_log_aracne_bootstrap_2.log", append = T)
      aracne(expr = exprFile, tfs = tfs, outputFolder = folder, seed = ar1)
    }
    
    net = aracne(outputFolder = folder, consolidate = TRUE)
    if (ncol(net) == 4) net = net[, -4]
    colnames(net) = c("tf", "target", "mi")
    # netFile = "../rData/ARACNe_network_Lyme_GSE63085.txt"
    write.table(net, file = netFile, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  phenoData <- new("AnnotatedDataFrame", data=pData)
  eset = ExpressionSet(assayData = expr, phenoData=phenoData)
  
  library(viper)
  
  file = "../rData/102_mrs_2.Rdata"
  if(file.exists(file)){
    load(file)
  } else {
    regul <- aracne2regulon(afile = netFile, eset = eset,
                            format = "3col", verbose = FALSE)
    
    # signature0 <- rowTtest(eset, pheno = "group", group1 = "A", group2 = "B")
    signature0 <- rowTtest(eset, pheno = "group", group1 = "B", group2 = "A")
    
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    set.seed(123)
    # nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "A",
    #                        group2 = "B", per = 1000,
    #                        repos = TRUE, verbose = FALSE)
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "B",
                           group2 = "A", per = 1000,
                           repos = TRUE, verbose = FALSE)
    
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    mrs_summary = viper:::summary.msviper(mrs, 50) # summary(mrs)
    save(mrs, mrs_summary, regul, file = file)
  }
  
  
  ## Down-sampling
  # regul <- aracne2regulon(afile = netFile, eset = eset,
  #                         format = "3col", verbose = FALSE)
  for(mi in c(100, 50, 20, 10, 6)){
    file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), "_2.Rdata")
    if(file.exists(file)){
      load(file)
    } else {
      id = c(1:(mi/2), 51:(mi/2 + 50))
      # regul <- aracne2regulon(afile = netFile, eset = eset[,id],
      #                         format = "3col", verbose = FALSE)
      
      # signature0 <- rowTtest(eset[,id], pheno = "group", group1 = "A", group2 = "B")
      signature0 <- rowTtest(eset[,id], pheno = "group", group1 = "B", group2 = "A")
      
      signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                      sign(signature0$statistic))[, 1]
      set.seed(123)
      # nullmodel <- ttestNull(x = eset[,id], pheno = "group", group1 = "A",
      #                        group2 = "B", per = 1000,
      #                        repos = TRUE, verbose = FALSE)
      nullmodel <- ttestNull(x = eset[,id], pheno = "group", group1 = "B",
                             group2 = "A", per = 1000,
                             repos = TRUE, verbose = FALSE)
      
      mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
      mrs_summary = viper:::summary.msviper(mrs, 500) # summary(mrs)
      save(mrs, mrs_summary, file = file)
    }
  }
}

## Plot 
for(mi in c(100, 50, 20, 10, 6)){
  file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), "_2.Rdata")
  load(file)
  mrs_summaryAll = viper:::summary.msviper(mrs, 500) 
  
  sudoPathway = list(topReg = grep("top", mrs_summaryAll$Regulon, value = T))
  sudoRanks = setNames((length(mrs_summaryAll$Regulon):1)/500 - 0.5, mrs_summaryAll$Regulon)
  enrich = fgsea(sudoPathway, sudoRanks)
  p = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(mi, " samples (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  figFile = paste0("../fig/102_keyTF_", mi, "Samples_viper_2.png")
  ggsave(filename = figFile, p, width = 3, height = 1.8)
}
