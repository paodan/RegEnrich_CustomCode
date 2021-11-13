if(FALSE){
    library(funcTools)
    rQsub2("./", "9_evaluateOnInterferon_GSE30536_v0.99.14.R", 
           jobName = "Job_09", threaded = 2, memoryG = "60", 
           rTimeHour = 20, logFile = "logFile_09.log", 
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

############ TRAF6 and IRF7 control HIV replication in macrophages (IFN alpha 2)
geoAcession = "GSE30536"
folder = paste0("../rawData/Interferon/", geoAcession, "/")

library(GEOquery)
if(file.exists("../rData/9_gds.Rdata")){
    cat("load ../rData/9_gds.Rdata", "\n")
    load("../rData/9_gds.Rdata")
} else {
    exprFile = normalizePath(dir(folder, pattern = "*.txt", full.names = T))
    print(exprFile)
    gds <- getGEO(filename=exprFile)
    save(gds, file = "../rData/9_gds.Rdata")
}

{
    expr = exprs(gds)
    pd0 = pData(gds)
    pd = pd0[, c("title", "geo_accession", "characteristics_ch1.1", 
                 "description.1", "hiv infection:ch1", "time:ch1")]
    colnames(pd) = c("title", "geo_accession", "ifn",
                     "description", "hiv_infection", "time")
    pd$ifn = strSplit(pd$ifn, ": ")[,2]
    pd$time = as.numeric(strSplit(pd$time, " ")[,1])
    pd$group = as.factor(strSplit(pd$description, "_")[,3])
    pd$cell = "macrophage"
    pData(gds) = pd
    
    edata = gds[, pd$hiv_infection == "none"]
}


# eSet
library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF_name))
# Healthy donor

######## cells
cells = c("macrophage")
np = c("1p", "5p", "3p")[2]
for (cell in cells){
    file = paste0("../rData/object_09_", geoAcession, "_", cell, ".Rdata")
    if(file.exists(file)){
        cat("load ", file, "\n")
        load(file)
    } else {
        # IFNa2
        edata_2 = edata
        
        # 
        fdata = fData(edata_2)
        gene0 = strSplit(fdata$`Gene Symbol`, " /// ")[,1]
        expr = exprs(edata_2)[!duplicated(gene0),]
        fdata = fdata[!duplicated(gene0), ]
        # expr = expr[rowMeans(expr) > 4,]
        # expr = t(scale(t(expr)))
        expr = preprocessCore::normalize.quantiles(expr)
        rownames(expr) = gene0[!duplicated(gene0)]
        pdata = pData(edata_2)
        
        softPower = NULL
        design = ~ time + group
        reduced = ~ time
        
        ### RegEnrich
        library(BiocParallel)
        bpparam = register(MulticoreParam(2))
        object = RegenrichSet(expr, colData = pdata, method = "LRT_LM", 
                              minMeanExpr = NULL, 
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
    write.table(score, paste0("../rData/score_09_", geoAcession, "_", 
                              cell, "_FET.tsv"), sep = "\t")
    
    file2 = paste0("../rData/object_09_", geoAcession, "_", cell, "_GSEA.Rdata")
    if(file.exists(file2)){
        cat("load ", file2)
        load(file2)
    } else {
        object = object %>% 
            regenrich_enrich(enrichTest = "GSEA") %>% 
            regenrich_rankScore()
        save(object, file = file2)
    }
    score = results_score(object)
    write.table(score, paste0("../rData/score_09_", geoAcession, "_", 
                              cell, "_GSEA.tsv"), sep = "\t")
    
}





