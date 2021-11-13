
library(funcTools)
library(RegEnrich)
library(microbenchmark)

library(RegEnrich)
library(doParallel)
library(foreach)
library(ggplot2)
library(DOSE)
library(WGCNA)
library(igraph)
options(max.print = 100)

load("../rData/Lyme_GSE63085.rda")

### read the full dataset
if(file.exists("../rData/02_02_Lyme_GSE63085_full.Rdata")){
  load("../rData/02_02_Lyme_GSE63085_full.Rdata")
} else {
  # Download GSE63085 gene expression file
  url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63085/suppl/GSE63085_RAW.tar'
  tarFile = "../rawData/GSE63085_RAW-3.tar"
  tarFolder = "../rawData/GSE63085_RAW-3"
  download.file(url, tarFile, method = "wget")
  untar(tarfile = tarFile, exdir = tarFolder)
  
  fileNames = dir(tarFolder, full.names = T)
  fileNames0 = strSplit(dir(tarFolder), "_")[,1]
  FPKM_full = do.call("cbind", 
                      lapply(setNames(fileNames, fileNames0), function(mi){
                        y = read.csv(mi, sep = "\t")[, c("tracking_id", "FPKM")]
                        # keep the first row if duplicate
                        y = y[!duplicated(y[,1]),]
                        y = y[match(rownames(Lyme_GSE63085$FPKM), y$tracking_id),]
                        return(setNames(y[,2], y[,1]))
                      }))
  
  library(GEOquery)
  sampleInfo_full = pData(getGEO(GEO = "GSE63085")[[1]])
  
  Lyme_GSE63085_full = list(FPKM_full = FPKM_full, sampleInfo_full = sampleInfo_full)
  save(Lyme_GSE63085_full, file = "../rData/02_02_Lyme_GSE63085_full.Rdata")
}
