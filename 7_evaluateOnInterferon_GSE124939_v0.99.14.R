

if(FALSE){
    library(funcTools)
    rQsub2("./", "7_evaluateOnInterferon_GSE124939_v0.99.14.R", 
          jobName = "Job_07", threaded = 2, memoryG = "50", 
          rTimeHour = 20, logFile = "logFile_07.log", 
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



############ Determinants in Cutaneous Lupus Keratinocytes Reveal Key Mechanistic Hypersensitive IFN Responses in Lupus
folder = "../rawData/Interferon/GSE124939/"
geoAcession = "GSE124939"
{ dp = read.csv(paste0(folder, "GSE124939_series_matrix.txt"), sep = "\t", skip = 25, header = F)
  dp = as.data.frame(t(dp), stringsAsFactors = F)[, c(1,2,8,10,11,17, 18)]
  dp = dp[-1,]
  colnames(dp) = c("title", "accession", "cell", "diseaseStatus", "stimulation", "ID1", "ID2")
  library(funcTools)
  dp[,4] = strSplit(dp[,4], ": ")[,2]
  dp[,5] = strSplit(dp[,5], ": ")[,2]
  dp$ID = apply(dp[,c("ID1", "ID2")], 1, function(x) {
    if (x[2] == "") x[1] else x[2]})
  rownames(dp) = dp$title
  dp$stimulation = factor(dp$stimulation, 
                          levels = c("control","IFNa2","IFNa2_IL18", "IL18","IFNa6","IFNb","IFNg"))
}

{ expr = read.csv(paste0(folder, "GSE124939_countmatrix.txt"), sep = "\t")
  expr = expr[,paste0("X", dp$ID)]
  colnames(expr) = dp$title
  expr = expr[(rowMeans(expr)>5),]
}

# eSet
library(Biobase)
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
{
  edata = ExpressionSet(assayData = as.matrix(expr), 
                        phenoData = AnnotatedDataFrame(data = dp))
  exprs(edata)
  pData(edata)
  dim(edata)
}
# Healthy donor

######## IFNs
IFNs = c("IFNb","IFNg", "IFNa2", "IFNa6")
np = c("1p", "5p", "3p")[2]
cell = "keratinocytes"
for (IFN in IFNs){
    file = paste0("../rData/object_07_", geoAcession, "_", 
                  cell, "_", IFN, ".Rdata")
    if(file.exists(file)){
        cat("load ", file, "\n")
        load(file) 
    } else {
        
        # IFNa2
        edata_2 = edata[, pData(edata)$diseaseStatus == "normal"]
        edata_2 = edata_2[, pData(edata_2)$stimulation %in% c("control", IFN)]
        
        # 
        expr = exprs(edata_2)
        pdata = pData(edata_2)
        
        design = ~stimulation
        reduced = ~1
        
        softPower = if(IFN == IFNs[1]) 6 else  8
        
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
                              maxSize = 15000,
                              enrichTest = 'FET')
        print(dim(object))
        
        object = object %>% regenrich_diffExpr() %>% 
            regenrich_network() %>% 
            regenrich_enrich() %>% 
            regenrich_rankScore()
        
        save(object, file = file)
    }
    
    score = results_score(object)
    write.table(score, paste0("../rData/score_07_", geoAcession, "_", 
                              cell, "_", IFN, "_FET.tsv"), sep = "\t")
    
    ## COEN + GSEA
    file2 = paste0("../rData/object_07_", geoAcession, "_", 
                   cell, "_", IFN, "_GSEA.Rdata")
    if (file.exists(file2)){
        cat("load ", file2, "\n")
        load(file2)
    } else {
        object = object %>% 
            regenrich_enrich(enrichTest = "GSEA") %>% 
            regenrich_rankScore()
        save(object, file = file2)
    }
    score = results_score(object)
    write.table(score, paste0("../rData/score_07_", geoAcession, "_", 
                              cell, "_", IFN, "_GSEA.tsv"), sep = "\t")
    
}






