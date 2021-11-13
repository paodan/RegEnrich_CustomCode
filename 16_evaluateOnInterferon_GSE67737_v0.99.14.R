
if(FALSE){
    library(funcTools)
    rQsub2("./", "16_evaluateOnInterferon_GSE67737_v0.99.14.R", 
           jobName = "Job_16", threaded = 2, memoryG = "50", 
           rTimeHour = 20, logFile = "logFile_16.log", 
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
# GPL10558 Platform, Illumina HumanHT-12 V4.0 expression beadchip
# Total number of rows: 48107

############ Gene expression analysis of macrophages derived from ankylosing spondylitis 
# patients reveals interferon-gamma dysregulation.  (IFNg)
geoAcession = "GSE67737"
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
    pd = pd0[, c("title", "geo_accession", "cell type:ch1", "treatment:ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "group")
    pd$disease = strSplit(strSplit(pd$cell, ", ")[,2], " ")[,1]
    pd$cell = "fibroblast"
    pd$group = factor(pd$group, levels = c("none", "IFN-alpha", "IFN-beta", "IFN-gamma"))
    pd$time = 10
    pData(gds) = pd
    
    edata = gds[, pd$disease %in% c("contol")]
    
    dim(edata)
}
# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("fibroblast")
# times = c(10)
np = c("1p", "5p", "3p")[2]
groups = c("IFN-alpha", "IFN-beta", "IFN-gamma")
for (cell in cells){
    file = paste0("../rData/object_16_", geoAcession, "_", cell, "_time_", groups[1],
                  ".Rdata")
    if(file.exists(file)){
        cat("load ", file, "\n")
        load(file) 
    } else {
    ## total involved expr
    edata_1 = edata[, pData(edata)$cell == cell & pData(edata)$group %in% c("none", groups)]
    pData(edata_1)$group = droplevels(pData(edata_1)$group)
    
    # 
    fdata = fData(edata_1)
    gene0 = strSplit(fdata$Symbol, " /// ")[,1]
    expr1 = exprs(edata_1)[!duplicated(gene0),]
    fdata = fdata[!duplicated(gene0), ]
    expr1 = preprocessCore::normalize.quantiles(expr1)
    rownames(expr1) = gene0[!duplicated(gene0)]
    pdata1 = pData(edata_1)
    
    design = model.matrix(~1+group, pdata1)
    contrast = as.numeric(colnames(design) == paste0("group", groups[1]))
    
    
    library(BiocParallel)
    bpparam = register(MulticoreParam(2))
    object = RegenrichSet(expr1, colData = pdata1, method = "limma", 
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
    write.table(score, paste0("../rData/score_16_", geoAcession, "_", cell, "_time_", groups[1],
                              "_FET.tsv"), sep = "\t")
    
    
    ## COEN + GSEA
    file2 = paste0("../rData/object_16_", geoAcession, "_", cell, "_time_", groups[1],
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
    cat(file2, "\n")
    print(score)
    print(subset(as.data.frame(score), reg == "CIITA"))
    print(subset(as.data.frame(score), reg == "ETV7"))
    write.table(score, paste0("../rData/score_16_", geoAcession, "_", cell, "_time_", groups[1],
                              "_GSEA.tsv"), sep = "\t")
    
    
    network_user = results_topNet(object)
    
    for (time in groups[-1]){
        file = paste0("../rData/object_16_", geoAcession, "_", cell, "_time_", time,
                      ".Rdata")
        if(file.exists(file)){
            cat("load ", file, "\n")
            load(file) 
        } else {
        contrast = as.numeric(colnames(design) == paste0("group", time))
        object = regenrich_diffExpr(object, contrast = contrast)
        
        regenrich_network(object) = network_user
        
        object = object %>% regenrich_enrich() %>% 
            regenrich_rankScore()
        
        
        collectGarbage()
        
        save(object, file = file)
        }
        score = results_score(object)
        write.table(score, paste0("../rData/score_16_", geoAcession, "_", cell, "_time_", time,
                                  "_FET.tsv"), sep = "\t")
        
        
        ## COEN + GSEA
        file2 = paste0("../rData/object_16_", geoAcession, "_", cell, "_time_", time,
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
        cat(file2, "\n")
        print(score)
        print(subset(as.data.frame(score), reg == "CIITA"))
        print(subset(as.data.frame(score), reg == "ETV7"))
        write.table(score, paste0("../rData/score_16_", geoAcession, "_", cell, "_time_", time,
                                  "_GSEA.tsv"), sep = "\t")
    }
}
