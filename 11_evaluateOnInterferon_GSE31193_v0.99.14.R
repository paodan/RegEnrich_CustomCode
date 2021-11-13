
if(FALSE){
    library(funcTools)
    rQsub2("./", "11_evaluateOnInterferon_GSE31193_v0.99.14.R", 
           jobName = "Job_11", threaded = 2, memoryG = "50", 
           rTimeHour = 20, logFile = "logFile_11.log", 
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
# GPL570 platform, [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# Total number of rows: 54675


############ HCV infection induces a unique hepatic innate immune response associated with robust production of type III interferons. (IFN-α)
geoAcession = "GSE31193"
if(file.exists("../rData/11_gds.Rdata")){
    cat("load ", "../rData/11_gds.Rdata", "\n")
    load("../rData/11_gds.Rdata")
} else {
    folder = paste0("../rawData/Interferon/", geoAcession, "/")
    library(GEOquery)
    exprFile = normalizePath(dir(folder, pattern = "*.txt", full.names = T))
    print(exprFile)
    gds <- getGEO(filename=exprFile)
    save(gds, file = "../rData/11_gds.Rdata")
}

{
    expr = exprs(gds)
    pd0 = pData(gds)
    pd = pd0[, c("title", "geo_accession", "tissue:ch1", "time:ch1", "agent:ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "time", "group")
    pd$cell = "Hepatocytes" # "Primary Human Hepatocytes"
    pd$time = c(`n/a` = 0, `6` = 6, `24` = 24)[pd$time]
    pd$group = factor(pd$group, levels = c("none", "IFN", "IL28B"))
    pData(gds) = pd
    
    edata = gds[, pd$group %in% c("none", "IFN")]
}


# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("Hepatocytes")
times = c(6, 24)
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    for (time in times){
        file = paste0("../rData/object_11_", geoAcession, "_", cell, 
                      "_time", time, ".Rdata")
        if(file.exists(file)){
            cat("load ", file, "\n")
            load(file)
        } else {
            # IFNa
            edata_2 = edata[, pData(edata)$cell == cell & pData(edata)$time %in% c(0, time)]
            pData(edata_2)$group = droplevels(pData(edata_2)$group)
            
            # 
            fdata = fData(edata_2)
            gene0 = strSplit(fdata$`Gene Symbol`, " /// ")[,1]
            expr = exprs(edata_2)[!duplicated(gene0),]
            fdata = fdata[!duplicated(gene0), ]
            # expr = preprocessCore::normalize.quantiles(expr)
            rownames(expr) = gene0[!duplicated(gene0)]
            # expr = expr[rowMeans(expr) > 3, ]
            pdata = pData(edata_2)
            
            ### RegEnrich
            design = ~1+group # The same as design = ~ group
            contrast = c(0,1)
            
            softPower = if(cell == "SKOV3") 20 else NULL
            
            
            library(BiocParallel)
            bpparam = register(MulticoreParam(2))
            object = RegenrichSet(expr, colData = pdata, method = "limma", 
                                  minMeanExpr = 3, 
                                  design = design, contrast = contrast, 
                                  reg = regulators, 
                                  BPPARAM = bpparam,
                                  networkConstruction = 'COEN',
                                  softPower = softPower,
                                  enrichTest = 'FET'
            )
            print(dim(object))
            
            object = object %>% regenrich_diffExpr() %>% 
                regenrich_network() %>% 
                regenrich_enrich(maxSize = 15000) %>% 
                regenrich_rankScore()
            
            collectGarbage()
            save(object, file = file)
        }
        score = results_score(object)
        write.table(score, paste0("../rData/score_11_", geoAcession, "_", cell, 
                                  "_time", time, "_FET.tsv"), sep = "\t")
        
        file2 = paste0("../rData/object_11_", geoAcession, "_", cell, 
               "_time", time, "_GSEA.Rdata")
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
        write.table(score, paste0("../rData/score_11_", geoAcession, "_", cell, 
                                  "_time", time, "_GSEA.tsv"), sep = "\t")
    }
}
