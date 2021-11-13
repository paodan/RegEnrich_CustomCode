
if(FALSE){
  options(width = 300)
  rm(list = ls())
  setwd("/hpc/dla_lti/wtao/3_RegEnrich_paper/script/")
  library(funcTools)
  rFile = "102_Revise_Reviewer2_Q1_viper_diffNetwork.R"
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = rFile, threaded = 1, 
         memoryG = 40, rTimeHour = 24, 
         logFile = paste0("../rData/102_Revise_Reviewer2_Q1_viper_diffNetwork.log"), 
         email = "weiyangtao1513@gmail.com", when2Email = "aes",
         preCMD = "Rscript ", 
         param1 = 1)
  qstat2("all")[,1:12]
  
}

{
  nTF = 500
  nTF_top = 20
  nSample = 100
  
  # Generating expression data and phenotype data
  fileName = "../rData/100_expr_pData_simulate2.Rdata"
  load(fileName)
  
  
  library(ARACNe)
  library(viper)
  
  ## Down-sampling
  for(mi in c(100, 50, 20)){ # when mi = 10 or 6, ARACN3 fails
    file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), "samples_3.Rdata")
    if(file.exists(file)){
      load(file)
    } else {
      id = c(1:(mi/2), 51:(mi/2 + 50))
      
      exprFile = paste0("../rData/102_expr_2_", mi, "samples.txt")
      tfs = "../rData/102_tfs.txt"
      
      log2FPKM1 = data.frame(gene = rownames(expr[,id]), expr[,id])
      write.table(log2FPKM1, file = exprFile, sep = "\t", 
                  row.names = FALSE, quote = F)
      write.table(grep("TF", rownames(expr[,id]), value = T), file = tfs, quote = F, 
                  sep = "\t", row.names = FALSE, col.names = FALSE)
      
      
      netFile = paste0("../rData/102_ARACNe_network_simulate_2_", mi, "samples.txt")
      if (!file.exists(netFile)){
        set.seed(123)
        folder = paste0("../rData/102_ARACNe_network_simulate_outputFolder_2_", mi, "samples/")
        if(!dir.exists(folder)) dir.create(folder)
        
        t = aracne(expr = exprFile, tfs = tfs, outputFolder = folder,
                   calculateThreshold = T)
        for(ar1 in 1:100){
          cat(ar1, "\n")
          write.table(ar1, file = paste0("./102_log_aracne_bootstrap_2_", mi, "samples.log"), append = T, col.names = F)
          aracne(expr = exprFile, tfs = tfs, outputFolder = folder, seed = ar1)
        }
        
        net = aracne(outputFolder = folder, consolidate = TRUE)
        if (ncol(net) == 4) net = net[, -4]
        colnames(net) = c("tf", "target", "mi")
        # netFile = "../rData/ARACNe_network_Lyme_GSE63085.txt"
        write.table(net, file = netFile, quote = FALSE, 
                    sep = "\t", row.names = FALSE, col.names = FALSE)
      }
      
      phenoData <- new("AnnotatedDataFrame", data=pData[id,])
      eset = ExpressionSet(assayData = expr[,id], phenoData=phenoData)
      
      
      # viper
      regul <- aracne2regulon(afile = netFile, eset = eset,
                              format = "3col", verbose = FALSE)
      
      signature0 <- rowTtest(eset, pheno = "group", group1 = "B", group2 = "A")
      
      signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                      sign(signature0$statistic))[, 1]
      set.seed(123)
      nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "B",
                             group2 = "A", per = 1000,
                             repos = TRUE, verbose = FALSE)
      
      mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
      mrs_summary = viper:::summary.msviper(mrs, 500) # summary(mrs)
      
      save(mrs, mrs_summary, regul, file = file)
    }
  }
}

## Plot 
library(ggplot2)
library(fgsea)
# 10, 6 fails
for(mi in c(100, 50, 20)){
  file = paste0("../rData/102_mrs_", funcTools::formatN(mi, 2), "samples_3.Rdata")
  load(file)
  mrs_summaryAll = viper:::summary.msviper(mrs, 500) 
  
  sudoPathway = list(topReg = grep("top", mrs_summaryAll$Regulon, value = T))
  sudoRanks = setNames((length(mrs_summaryAll$Regulon):1)/500 - 0.5, mrs_summaryAll$Regulon)
  enrich = fgsea(sudoPathway, sudoRanks)
  p = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(mi, " samples (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
  # print(p)
  figFile = paste0("../fig/102_keyTF_", mi, "Samples_viper_3.png")
  ggsave(filename = figFile, p, width = 3, height = 1.8)
}


