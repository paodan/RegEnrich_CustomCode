
if(FALSE){
    library(funcTools)
    rQsub2("./", "15_evaluateOnInterferon_GSE11886_v0.99.14.R", 
           jobName = "Job_15", threaded = 2, memoryG = "50", 
           rTimeHour = 20, logFile = "logFile_15.log", 
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

############ Gene expression analysis of macrophages derived from ankylosing spondylitis 
# patients reveals interferon-gamma dysregulation.  (IFNg)
geoAcession = "GSE11886"
{
    folder = paste0("../rawData/Interferon/", geoAcession, "/")
    library(GEOquery)
    exprFile = normalizePath(dir(folder, pattern = "*.txt", full.names = T))
    print(exprFile)
    gds <- getGEO(filename=exprFile)
}
{
    expr = exprs(gds)
    pd0 = pData(gds)
    pd = pd0[, c("title", "geo_accession", "source_name_ch1", "characteristics_ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "group")
    pd$cell = "pbMacrophages"
    pd$disease = strSplit(pd$group, " ")[,1]
    pd$group = factor(strSplit(strSplit(pd$group, ", ")[,2], " ")[,1], 
                      levels = c("Untreated", "IFN"))
    pd$time = 24
    pData(gds) = pd
    
    edata = gds[, pd$disease %in% c("Control") & pd$group %in% c("Untreated", "IFN")]
    
    dim(edata)
}
# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("pbMacrophages")
times = c(24)
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    for (time in times){
        # IFNg
        file = paste0("../rData/object_15_", geoAcession, "_", cell, "_time_", time,
                      ".Rdata")
        if(file.exists(file)){
            cat("load ", file, "\n")
            load(file) 
        } else {
            
            edata_2 = edata[, pData(edata)$cell == cell & pData(edata)$time %in% c(time)]
            pData(edata_2)$group = droplevels(pData(edata_2)$group)
            
            # 
            fdata = fData(edata_2)
            gene0 = strSplit(fdata$`Gene Symbol`, " /// ")[,1]
            expr = exprs(edata_2)[!duplicated(gene0),]
            fdata = fdata[!duplicated(gene0), ]
            rownames(expr) = gene0[!duplicated(gene0)]
            pdata = pData(edata_2)
            
            
            design = ~1+group
            contrast = c(0,1)
            
            ### RegEnrich
            
            library(BiocParallel)
            bpparam = register(MulticoreParam(2))
            object = RegenrichSet(expr, colData = pdata, method = "limma", 
                                  minMeanExpr = -10, 
                                  design = design, contrast = contrast, 
                                  reg = regulators, 
                                  BPPARAM = bpparam,
                                  networkConstruction = 'COEN',
                                  softPower = NULL,
                                  maxSize = 15000,
                                  enrichTest = 'FET')
            print(dim(object))
            
            object = object %>% regenrich_diffExpr() %>% 
                regenrich_network() %>% 
                regenrich_enrich() %>% 
                regenrich_rankScore()
            
            collectGarbage()
            
            save(object, file = file)
        }
        
        score = results_score(object)
        write.table(score, paste0("../rData/score_15_", geoAcession, "_", cell, 
                                  "_time_", time, "_FET.tsv"), sep = "\t")
        
        ## COEN + GSEA
        file2 = paste0("../rData/object_15_", geoAcession, "_", cell, "_time_", time,
                       "_GSEA.Rdata")
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
        print(subset(as.data.frame(score), reg == "CIITA"))
        print(subset(as.data.frame(score), reg == "ETV7"))
        write.table(score, paste0("../rData/score_15_", geoAcession, "_", cell, 
                                  "_time_", time, "_GSEA.tsv"), sep = "\t")
        
    }
}






