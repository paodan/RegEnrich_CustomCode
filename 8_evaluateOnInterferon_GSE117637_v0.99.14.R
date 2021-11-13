
if(FALSE){
    library(funcTools)
    rQsub2("./", "8_evaluateOnInterferon_GSE117637_v0.99.14.R", 
          jobName = "Job_08", threaded = 2, memoryG = "60", 
          rTimeHour = 20, logFile = "logFile_08.log", 
          email = "w.tao-2@umcutrecht.nl", 
          preCMD = "Rscript ", param1 = 1)
    qstat2()
}

# setwd("/hpc/dla_lti/wtao/")
rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)
commandArgs()

source("0_functionsInThePaper.R")

# GPL20301 platform, Illumina HiSeq 4000 (Homo sapiens)

############ Life-threatening influenza pneumonitis in a child with inherited IRF9 deficiency
folder = "../rawData/Interferon/GSE117637/"
geoAcession = "GSE117637"
{ 
  dp = read.csv(paste0(folder, "GSE117637_series_matrix.txt"), sep = "\t", skip = 56, header = F)
  dp = as.data.frame(t(dp), stringsAsFactors = F)[, c(1,2,8)]
  dp = dp[-1,]
  
  colnames(dp) = c("title", "accession", "cell")#, "diseaseStatus", "stimulation", "ID1", "ID2")
  dp$cell[dp$cell == "B-cell"] = "Bcell"
  tmp = strSplit(dp$title, "_")
  dp$donor = paste0(tmp[,1], ".", ("repeat" == tmp[,2]) + 1)
  for(mi in 1:3){
    dp$donor[dp$donor == paste0("Control",mi,".1")] = paste0("Control.",mi)
  }
  dp$stimulation = c("non", "IFNa2b")[(tmp[,2] == "IFNa2b" | tmp[,3] == "IFNa2b") + 1]
  dp$stimulation = factor(dp$stimulation, 
                          levels = c("non","IFNa2b"))
  rownames(dp) = paste(stringr:::str_to_lower(dp$cell), dp$donor, dp$stimulation, sep = "_")
  
  # library(funcTools)
  # dp[,4] = strSplit(dp[,4], ": ")[,2]
  # dp[,5] = strSplit(dp[,5], ": ")[,2]
  # dp$ID = apply(dp[,c("ID1", "ID2")], 1, function(x) {
  #   if (x[2] == "") x[1] else x[2]})
  # rownames(dp) = dp$title
}

{ 
  expr0 = read.csv(paste0(folder, "GSE117637_GEO_normalized_log.csv"), sep = ",", stringsAsFactors = F)
  sampleName = strSplit(substr(colnames(expr0)[3:38], 15, 100), "\\.")[,1]
  sampleName2 = gsub("IRF9_", "IRF9.1_", 
                     gsub("STAT1_", "STAT1.1_", 
                          gsub("STAT2_", "STAT2.1_", 
                               gsub("_repeat", ".2", 
                                    gsub("ifna2b|ifn2ab", "IFNa2b", 
                                         gsub("ns", "non", 
                                              gsub("[Cc]ontrol", "Control.", 
                                                   sampleName)))))))
  expr = subset(expr0, !duplicated(expr0$EnsemblID))[, -2]
  rownames(expr) = expr[,1]
  expr = expr[,-1]
  colnames(expr) = sampleName2
  expr = expr[, rownames(dp)]
  
  expr = expr[rowMeans(expr) > 1,]
  dim(expr)
  # gene IDs
  geneIDs = subset(expr0, !duplicated(expr0$EnsemblID))[, 1:2]
  rownames(geneIDs) = geneIDs[,1]
}

# eSet
library(Biobase)
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF))
{
  edata = ExpressionSet(assayData = as.matrix(expr), 
                        phenoData = AnnotatedDataFrame(data = dp))
  exprs(edata)
  pData(edata)
  dim(edata)
}
# Healthy donor

######## cells
cells = c("Bcell", "Fibroblast")
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    file = paste0("../rData/object_08_", geoAcession, "_", cell, ".Rdata")
    if(file.exists(file)){
        cat("load ", file)
        load(file)
    } else {
        
        # IFNa2
        edata_1 = edata[, pData(edata)$cell == cell]
        edata_2 = edata_1[, grep("Control", pData(edata_1)$donor)]
        
        # 
        expr = exprs(edata_2)
        pdata = pData(edata_2)
        
        softPower = if(cell == "Fibroblast") 18 else NULL
        design = ~ donor + stimulation
        reduced = ~ donor
        
        
        ### RegEnrich
        library(BiocParallel)
        bpparam = register(MulticoreParam(2))
        object = RegenrichSet(expr, colData = pdata, method = "LRT_LM", 
                              minMeanExpr = 1, 
                              design = design, reduced = reduced, 
                              reg = regulators, 
                              BPPARAM = bpparam,
                              networkConstruction = 'COEN',
                              softPower = softPower,
                              enrichTest = 'FET')
        print(dim(object))
        
        object = object %>% regenrich_diffExpr() %>% 
            regenrich_network() %>% 
            regenrich_enrich() %>% 
            regenrich_rankScore()
        
        save(object, file = file)
    }
    score = results_score(object)
    score$reg = TFs[score$reg, 2]
    write.table(score, paste0("../rData/score_08_", geoAcession, "_", 
                              cell, "FET.tsv"), sep = "\t")
    
    ## GSEA
    file2 = paste0("../rData/object_08_", geoAcession, "_", 
                   cell, "_GSEA.Rdata")
    if(file.exists(file2)){
        cat("load", file2)
        load(file2)
    } else {
        object = object %>% 
            regenrich_enrich(enrichTest = "GSEA") %>% 
            regenrich_rankScore()
        save(object, file = file2)
    }
    score = results_score(object)
    score$reg = TFs[score$reg, 2]
    write.table(score, paste0("../rData/score_08_", geoAcession, "_", 
                              cell, "_GSEA.tsv"), sep = "\t")
    
}





