

if(FALSE){
    library(funcTools)
    rQsub2("./", "14_evaluateOnInterferon_GSE12109_v0.99.14.R", 
           jobName = "Job_14", threaded = 2, memoryG = "60", 
           rTimeHour = 20, logFile = "logFile_14.log", 
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
# GPL571 platform, [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array
# Total number of rows: 22277

############ Th17 cytokines interleukin (IL)-17 and IL-22 modulate distinct inflammatory 
# and keratinocyte-response pathways (IFNg)
geoAcession = "GSE12109"
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
    pd = pd0[, c("title", "geo_accession", "source_name_ch1", "Cytokine:ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "group")
    pd$cell = "keratinocytes"
    pd$group = factor(pd$group, levels = c("control", "IFNg", "IL17", "IL22"))
    pd$time = 24
    pData(gds) = pd
    
    edata = gds[, pd$group %in% c("control", "IFNg")]
    
    dim(edata)
}
# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("keratinocytes")
times = c(24)
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    for (time in times){
        # IFNg
        file = paste0("../rData/object_14_", geoAcession, "_", cell, 
                      "_time_", time, ".Rdata")
        if(file.exists(file)){
            cat("load", file)
            load(file)
        } else {
            edata_2 = edata[, pData(edata)$cell == cell & pData(edata)$time %in% c(time)]
            pData(edata_2)$group = droplevels(pData(edata_2)$group)
            
            # 
            fdata = fData(edata_2)
            gene0 = strSplit(fdata$`Gene Symbol`, " /// ")[,1]
            expr = exprs(edata_2)[!duplicated(gene0),]
            fdata = fdata[!duplicated(gene0), ]
            expr = preprocessCore::normalize.quantiles(expr)
            rownames(expr) = gene0[!duplicated(gene0)]
            # expr = expr[rowMeans(expr) > 3, ]
            pdata = pData(edata_2)
            
            
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
                                  softPower = 18,
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
        write.table(score, paste0("../rData/score_14_", geoAcession, "_", cell, 
                                  "_time_", time, "_FET.tsv"), sep = "\t")
        
        ## COEN + GSEA
        file2 = paste0("../rData/object_14_", geoAcession, "_", cell, "_time_", time,
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
        write.table(score, paste0("../rData/score_14_", geoAcession, "_", cell, 
                                  "_time_", time, "_GSEA.tsv"), sep = "\t")
    }
}


