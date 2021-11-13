

if(FALSE){
    library(funcTools)
    rQsub2("./", "13_evaluateOnInterferon_GSE36537_v0.99.14.R", 
           jobName = "Job_13", threaded = 2, memoryG = "50", 
           rTimeHour = 20, logFile = "logFile_13.log", 
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
# GPL6480 platform, Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)
# Total number of rows: 41108

############ Monocyte responses in the context of Q fever: from a static polarized model 
# to a kinetic model of activation (IFNg)
geoAcession = "GSE36537"
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
    pd = pd0[, c("title", "geo_accession", "cell type:ch1", "time:ch1", "treatment:ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "time", "group")
    pd$cell = c(`Monocytes-Derived Macrophages` = "moMacrophages",
                Monocyte = "Monocytes")[pd$cell] 
    pd$group = factor(pd$group, levels = c("NS", "IFNg", "IL4"))
    pd$time = as.numeric(strSplit(pd$time, "h")[,1])
    pData(gds) = pd
    
    edata = gds[, pd$group %in% c("NS", "IFNg")]
    
    dim(edata)
}
# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("Monocytes", "moMacrophages")# Monocyte-derived macrophages
times = c(6, 18)
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    if (cell == "moMacrophages"){
        times = 18
    }
    for (time in times){
        # IFNg
        file = paste0("../rData/object_13_", geoAcession, "_", cell, "_time_", time,
                      ".Rdata")
        if(file.exists(file)){
            cat("load", file, "\n")
            load(file)
        } else {
            edata_2 = edata[, pData(edata)$cell == cell & pData(edata)$time %in% c(time)]
            pData(edata_2)$group = droplevels(pData(edata_2)$group)
            
            # 
            fdata = fData(edata_2)
            gene0 = fdata$GENE_SYMBOL
            expr = exprs(edata_2)[!duplicated(gene0),]
            fdata = fdata[!duplicated(gene0), ]
            expr = preprocessCore::normalize.quantiles(expr)
            rownames(expr) = gene0[!duplicated(gene0)]
            # expr = expr[rowMeans(expr) > 3, ]
            pdata = pData(edata_2)
            
            
            softPower = if(cell == "SKOV3") 20 else NULL
            collectGarbage()
            
            design = ~1+group
            contrast = c(0,1)
            
            ### RegEnrich
            library(BiocParallel)
            bpparam = register(MulticoreParam(2))
            object = RegenrichSet(expr, colData = pdata, method = "limma", 
                                  minMeanExpr = 1, 
                                  design = design, contrast = contrast, 
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
            
            collectGarbage()
            save(object, file = file)
        }
        score = results_score(object)
        write.table(score, paste0("../rData/score_13_", geoAcession, "_", cell, "_time_", time,
                                  "_FET.tsv"), sep = "\t")
        
        
        ## COEN + GSEA
        file2 = paste0("../rData/object_13_", geoAcession, "_", cell, "_time_", time,
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
        cat("cell_", cell, "_time_", time, "_GSEA\n")
        print(score)
        print(subset(as.data.frame(score), reg == "CIITA"))
        print(subset(as.data.frame(score), reg == "ETV7"))
        write.table(score, paste0("../rData/score_13_", geoAcession, "_", cell, "_time_", time,
                                  "_GSEA.tsv"), sep = "\t")
    }
}






