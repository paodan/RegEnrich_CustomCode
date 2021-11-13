
if(FALSE){
    library(funcTools)
    rQsub2("./", "12_evaluateOnInterferon_GSE79077_v0.99.14.R", 
           jobName = "Job_12", threaded = 2, memoryG = "50", 
           rTimeHour = 20, logFile = "logFile_12.log", 
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
# GPL17077 platform, Agilent-039494 SurePrint G3 Human GE v2 8x60K Microarray 039381 (Probe Name version)
# Total number of rows: 50737

############ Imatinib Triggers Phagolysosome Acidification and Antimicrobial Activity against 
# Mycobacterium bovis Bacille Calmette-Guérin in Glucocorticoid-Treated Human Macrophages. (IFNg)
geoAcession = "GSE79077"
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
    pd = pd0[, c("title", "geo_accession", "cell type:ch1", "source_name_ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "group")
    pd$cell = "moMacrophages" # Monocyte-derived macrophage
    pd$group = as.factor(c(media = "ctrl", dexa = "dexa", ifng = "IFN", 
                           `ifng + dexa` = "IFN_dexa")[strSplit(pd$group, "_")[,1]])
    pData(gds) = pd
    
    edata = gds[, pd$group %in% c("ctrl", "IFN")]
    
    dim(edata)
}

# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("moMacrophages")
np = c("1p", "5p", "3p")[2]

for (cell in cells){
    # IFNg
    file = paste0("../rData/object_12_", geoAcession, "_", cell, ".Rdata")
    if(file.exists(file)){
        cat("load", file, "\n")
        load(file)
    } else {
        edata_2 = edata[, pData(edata)$cell == cell]
        pData(edata_2)$group = droplevels(pData(edata_2)$group)
        
        # 
        fdata = fData(edata_2)
        gene0 = fdata$GENE_SYMBOL
        expr = exprs(edata_2)[!duplicated(gene0),]
        fdata = fdata[!duplicated(gene0), ]
        rownames(expr) = gene0[!duplicated(gene0)]
        expr = expr[rowMeans(expr) > -6, ]
        dim(expr)
        tmp = dimnames(expr)
        expr = preprocessCore::normalize.quantiles(expr)
        dimnames(expr) = tmp
        pdata = pData(edata_2)
        
        design = ~1+group
        contrast = c(0,1)
        
        ### RegEnrich
        library(BiocParallel)
        bpparam = register(MulticoreParam(2))
        object = RegenrichSet(expr, colData = pdata, method = "limma", 
                              minMeanExpr = NULL, 
                              design = design, contrast = contrast, 
                              reg = regulators, 
                              BPPARAM = bpparam,
                              networkConstruction = 'COEN',
                              # softPower = softPower,
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
    write.table(score, paste0("../rData/score_12_", geoAcession, "_", cell, 
                              "_FET.tsv"), sep = "\t")
    
    
    ## COEN + GSEA
    file2 = paste0("../rData/object_12_", geoAcession, "_", cell, 
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
    write.table(score, paste0("../rData/score_12_", geoAcession, "_", cell, 
                              "_GSEA.tsv"), sep = "\t")
}






