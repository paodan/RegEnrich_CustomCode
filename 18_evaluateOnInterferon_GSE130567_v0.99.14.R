

if(FALSE){
    library(funcTools)
    rQsub2("./", "18_evaluateOnInterferon_GSE130567_v0.99.14.R", 
           jobName = "Job_18", threaded = 2, memoryG = "80", 
           rTimeHour = 20, logFile = "logFile_18.log", 
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
# GPL20301 	Illumina HiSeq 4000 (Homo sapiens)


############ IFN-g selectively suppresses a subset of TLR4-activated genes and 
############ enhancers to potentiate M1-like macrophage polarization
geoAcession = "GSE130567"
destFolder = paste0("../rawData/Interferon/", geoAcession)
library(GEOquery)

# pheno data
{
    exprFile = normalizePath(dir(destFolder, pattern = "_series_matrix.txt", full.names = T))
    print(exprFile)
    gds <- getGEO(filename=exprFile)
    pd0 = pData(gds)
    pd = pd0[, c("title", "geo_accession", "cell type:ch1", "cultured in:ch1", "donor:ch1", "treatment:ch1")]
    colnames(pd) = c("title", "geo_accession", "cell", "medium", "donor", "treatment")
    pd$group = paste0(pd$medium, "_", pd$treatment)
}
# download read file
if (FALSE){
    destFile = paste0(destFolder, "/", geoAcession, "_RAW.tar")
    fileURL = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", geoAcession, "&format=file")
    download.file(fileURL, destfile = destFile)
    untar(destFile, exdir = destFolder)
    file.remove(destFile)
    gzfiles = dir(destFolder, pattern = ".gz", full.names = T)
    for (mi in gzfiles){
        if(!file.exists(substr(mi, 1, nchar(mi) - 3))){
            gunzip(mi, remove = T)
        } else {
            file.remove(mi)
        }
    }
}

## read data
datafile = paste0("../rData/data_", geoAcession, ".Rdata")
if(file.exists(datafile)){
    cat("loading data ...\n")
    load(datafile)
} else {
    # read tab files
    tabFiles = dir(destFolder, pattern = ".tab", full.names = T)
    data_tab = read.csv(tabFiles[1], header = F, sep = "\t", row.names = 1)[,1,drop = F]
    colnames(data_tab) = strSplit(strSplit(tabFiles[1], "/")[,5], split = "_")[,1]
    for (mi in tabFiles[-1]){
        colMi = strSplit(strSplit(mi, "/")[,5], split = "_")[,1]
        dataMi = read.csv(mi, header = F, sep = "\t", row.names = 1)[,1]
        if (all(rownames(dataMi) == rownames(data_tab))){
            data_tab[[colMi]] = dataMi
        } else {
            stop("Row names are not identical", mi)
        }
    }
    data_tab = data_tab[-grep("^N|_PAR", rownames(data_tab)), ,drop = FALSE]
    
    # read txt files
    txtFiles = dir(destFolder, pattern = "[^_series_matrix].txt", full.names = T)
    data_txt = read.csv(txtFiles[1], header = F, sep = "\t", row.names = 1)[,1,drop = F]
    colnames(data_txt) = strSplit(strSplit(txtFiles[1], "/")[,5], split = "_")[,1]
    for (mi in txtFiles[-1]){
        colMi = strSplit(strSplit(mi, "/")[,5], split = "_")[,1]
        dataMi = read.csv(mi, header = F, sep = "\t", row.names = 1)[,1]
        if (all(rownames(dataMi) == rownames(data_txt))){
            data_txt[[colMi]] = dataMi
        } else {
            stop("Row names are not identical", mi)
        }
    }
    data_txt = data_txt[-grep("^N|_PAR", rownames(data_txt)), ,drop = FALSE]
    save(data_tab, data_txt, file = datafile)
}

all(subset(pd, medium == "M-CSF")$geo_accession == colnames(data_tab))
all(subset(pd, medium == "GM-CSF")$geo_accession == colnames(data_txt))

pd_tab = subset(pd, medium == "M-CSF")
pd_txt = subset(pd, medium == "GM-CSF") # IFG-48h is IFNg-48h

### download fastq file from SRA database and perform alignment
if (FALSE){
    srrIDsFile = paste0("../rData/", geoAcession, "_srrIDs.Rdata")
    if(file.exists(srrIDsFile)){
        load(srrIDsFile)
    } else {
        # control and IFN treated samples
        GSM_to_download = rownames(subset(rbind(pd_tab, pd_txt), treatment %in% c("NT", "IFNG-3h", "IFG-48h")))
        srr = searchSrrID("PRJNA540657")
        srrSampleName = subset(srr, SampleName %in% GSM_to_download)[, c("Run", "SampleName")]
        srrIDs = subset(srr, SampleName %in% GSM_to_download)$Run
        save(srrIDs, srrSampleName, file = srrIDsFile)
    }
    
    OutDir = "../rawData/Interferon/GSE130567/SRR_fastq"
    if(!dir.exists(OutDir)){
        dir.create(OutDir, recursive = TRUE)
    }
    if(!all(srrIDs %in% dir(OutDir))){
        srrIDs1 = srrIDs[!paste0(srrIDs, "_1.fastq.gz") %in% dir(OutDir)]
        if(length(srrIDs1)>0){
            downloadSrr(srrIDs = srrIDs1, OutDir = OutDir, multipleDownload = 6)
        }
    }
    ### run alignment R script on HPC ###
}




########### to be finished ###########
## geneID and geneName
fileName = paste0("../rData/geneName_", geoAcession, ".Rdata")
if(file.exists(fileName)){
    cat("load(fileName)")
    load(fileName)
} else {
    geneID_M_CSF = rownames(data_tab) # ref_gene_id
    geneID_GM_CSF = rownames(data_txt)
    geneName_M_CSF = ENSGID2GeneName2019(strSplit(geneID_M_CSF, "\\.")[,1], sameOrder = T)
    rownames(geneName_M_CSF) = geneID_M_CSF
    geneName_GM_CSF = ENSGID2GeneName2019(strSplit(geneID_GM_CSF, "\\.")[,1], sameOrder = T)
    rownames(geneName_GM_CSF) = geneID_GM_CSF
    save(geneName_M_CSF, geneName_GM_CSF, file = fileName)
}

## ExpressionSet
{
rownames(subset(pd_tab, treatment %in% c("NT", "IFNG-3h", "IFG-48h")))
id1 = colnames(data_tab) %in% rownames(subset(pd_tab, treatment %in% c("NT", "IFNG-3h", "IFG-48h")))
gds_M_CSF = ExpressionSet(assayData = as.matrix(data_tab[, id1]),
                          phenoData=new("AnnotatedDataFrame",data=pd_tab[id1,], 
                                        varMetadata=data.frame(labelDescription = colnames(pd_tab[id1,]),
                                                               row.names = colnames(pd_tab[id1,]))))

id2 = colnames(data_txt) %in% rownames(subset(pd_txt, treatment %in% c("NT", "IFNG-3h", "IFG-48h")))
gds_GM_CSF = ExpressionSet(assayData = as.matrix(data_txt[,id2]),
                           phenoData=new("AnnotatedDataFrame",data=pd_txt[id2,], 
                                         varMetadata=data.frame(labelDescription = colnames(pd_txt[id2,]),
                                                                row.names = colnames(pd_txt[id2,]))))
}
groups = c("M_CSF", "GM_CSF")
cells = "monocytes"


library(RegEnrich)
data(TFs)
regulators = c(unique(TFs$TF))

for (cell in cells){
    for (time in groups){
        file = paste0("../rData/object_18_", geoAcession, "_", cell, "_time_", time,
                      ".Rdata")
        if(file.exists(file)){
            cat("load ", file, "\n")
            load(file) 
        } else {
            
            if (time == "M_CSF"){
                gds = gds_M_CSF
                minCnt = 1
            } else if (time == "GM_CSF"){
                gds = gds_GM_CSF
                minCnt = 5
            } else {
                gds = NULL
                minCnt = NULL
            }
            
            expr = exprs(gds)
            rownames(expr) = strSplit(rownames(expr), "\\.")[,1]
            pdata = pData(gds)
            pdata$treatment = gsub("-", ".", pdata$treatment)
            
            pdata$donor = as.factor(pdata$donor)
            pdata$treatment = factor(pdata$treatment, 
                                     levels = c("NT", setdiff(unique(pdata$treatment), "NT")))
            
            # DEG
            design = ~ donor + treatment
            reduced = ~ donor
            
            ### RegEnrich
            
            library(BiocParallel)
            bpparam = register(MulticoreParam(2))
            object = RegenrichSet(expr, colData = pdata, method = "LRT_DESeq2", 
                                  minMeanExpr = minCnt, 
                                  design = design, reduced = reduced, 
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
            save(object, file)
        }
        score = results_score(object)
        score$reg = TFs[score$reg, 2]
        write.table(score, paste0("../rData/score_18_", geoAcession, "_", cell, "_time_", time,
                                  "_FET.tsv"), sep = "\t")
        
        ## COEN + GSEA
        file2 = paste0("../rData/object_18_", geoAcession, "_", cell, "_time_", time,
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
        score$reg = TFs[score$reg, 2]
        cat(file2, "\n")
        print(score)
        write.table(score, paste0("../rData/score_18_", geoAcession, "_", cell, "_time_", time,
                                  "_GSEA.tsv"), sep = "\t")
    }
}

