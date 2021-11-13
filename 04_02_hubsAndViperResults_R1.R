rm(list = ls())
gc(reset = TRUE)
options(max.print = 100)

library(funcTools)
library(RegEnrich)
library(microbenchmark)

library(doParallel)
library(foreach)
library(ggplot2)
library(DOSE)
library(WGCNA)
library(igraph)

load("../rData/Lyme_GSE63085.rda")
log2FPKM = log2(Lyme_GSE63085$FPKM + 1)
sampleInfo = Lyme_GSE63085$sampleInfo
sampleInfo$week = factor(sampleInfo$week, levels = c(0, 3))

data(TFs)

################## Viper ##################
library(ARACNe)
netFile = "../rData/03_02_ARACNe_network_Lyme_GSE63085.txt"
aracneNet= read.csv(netFile, header = FALSE, sep = "\t", quote = "")
file = "../rData/03_02_regulatorOrderOfRegEnrichAndViper.Rdata"
load(file)
ARACNe_viper = as.character(ARACNe_viper)



##################### hubs #############################
# https://sna.stanford.edu/lab.php?l=4
# ####### hubness ##########
source("0_functionsInThePaper_2_testing_viper.R")
source("0_functionsInThePaper.R")
if(file.exists("../rData/04_02_central_viper.Rdata")){
  load("../rData/04_02_central_viper.Rdata")
} else {
  central_viper = hubness(aracneNet, reg = unique(as.character(aracneNet$V1)))
  rownames(central_viper) = central_viper$reg
  central_viper$viper = 0
  central_viper[ARACNe_viper, "viper"] = 1
  save(central_viper, file = "../rData/04_02_central_viper.Rdata")
}


# out degree and out closeness
g = topHubVenn(hubness = central_viper, colID = c(3,5, 9), topN = 50)
plotSave("../fig/04_02_Hubs_comparison_outDegree_outCloseness_viper_top50.svg", 
         Plot = g, width = 5, height = 5, dpi = 300)
