
if(!require("aracne.networks")){
  BiocManager::install("aracne.networks")
}
library(aracne.networks)
help(package = "aracne.networks")

source("0_functionsInThePaper.R")


library(RegEnrich)
data(TFs)


## Human Glioblastoma context-specific ARACNe interactome
# SNB19, BTIC 
data(regulongbm)

### Map gene entrez ID to gene name
fileName = "../rData/103_01_regul.Rdata"
if(file.exists(fileName)){
  cat("load", fileName, "\t")
  load(fileName)
} else {
  cat(fileName)
  regul = standardizeRegulon(regulongbm, filterTF = unique(TFs$TF_name))
  save(regul, file = fileName)
}


### Load GSE19114
load("../rData/gseGSE19114.Rdata")


if(TRUE){
  # The second batch, BTIC cells
  gse = gseGSE19114[, rownames(subset(pData(gseGSE19114), cell == "BTIC"))] 
  phenoGSE19114 = pData(gse)
  
  dataGSE19114 = exprs(gse)
  geneName = gse@featureData@data$ILMN_Gene
  geneNameValid = geneName[geneName != ""]
  dataGSE19114 = dataGSE19114[geneName != "", ]
  
  geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = "/")[,1]
  geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
  dataGSE19114 = dataGSE19114[geneNameValidNoDupID,]
  rownames(dataGSE19114) = geneNameValidNoMul[geneNameValidNoDupID]
  
  
  eset = ExpressionSet(assayData=dataGSE19114, 
                       phenoData=new("AnnotatedDataFrame", data=phenoGSE19114))
  
  ############# gene STAT3 ***************
  mrsFile = "../rData/103_01_mrs_GSE19114_ARACNe_2_STAT3.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "STAT3" knockout (ranked on number 11 by two sided method or 7 by one sided)
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                           group2 = "STAT3", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(regul))
  
  which(res$Regulon == "STAT3") ## 54
  
  

  ############# gene CEBPB ***************
  mrsFile = "../rData/103_01_mrs_GSE19114_ARACNe_2_CEBPB.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "CEBPB" knockout, not in the result
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "CEBPB")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT",
                           group2 = "CEBPB", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(regul))
  which(res$Regulon == "CEBPB") ## 454
  
  
  ############# gene STAT3 + CEBPB ***************
  mrsFile = "../rData/103_01_mrs_GSE19114_ARACNe_2_STAT3_CEBPB.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "STAT3" and "CEBPB" knockout, STAT3 ranked on 9 by two side or 6 by one side, CEBPB not found
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3_CEBPB")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                           group2 = "STAT3_CEBPB", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(regul))
  which(res$Regulon == "STAT3") ## 18
  which(res$Regulon == "CEBPB") ## 451
}
#


if(TRUE){
  
  ## The third batch, the last 12 samples, SNB19 cell ****
  gse = gseGSE19114[, tail(colnames(gseGSE19114), 12)] 
  phenoGSE19114 = pData(gse)
  
  dataGSE19114 = exprs(gse)
  geneName = gse@featureData@data$ILMN_Gene
  geneNameValid = geneName[geneName != ""]
  dataGSE19114 = dataGSE19114[geneName != "", ]
  
  geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = "/")[,1]
  geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
  dataGSE19114 = dataGSE19114[geneNameValidNoDupID,]
  rownames(dataGSE19114) = geneNameValidNoMul[geneNameValidNoDupID]
  
  
  eset = ExpressionSet(assayData=dataGSE19114, 
                       phenoData=new("AnnotatedDataFrame", data=phenoGSE19114))
  
  
  ########## STAT3 gene *******
  mrsFile = "../rData/103_01_mrs_GSE19114_ARACNe_3_STAT3.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "STAT3" knockout (not found)
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                           group2 = "STAT3", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(regul))
  
  which(res$Regulon == "STAT3") ## 63
  
  
  
  ###### CEBPB gene *****
  mrsFile = "../rData/103_03_mrs_GSE19114_ARACNe_3_CEBPB.Rdata"
  if (file.exists(mrsFile)){
    load(mrsFile)
    cat("loading mrsFile.\n")
  } else {
    # "CEBPB" knockout (not found)
    signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "CEBPB")
    signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                    sign(signature0$statistic))[, 1]
    nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                           group2 = "CEBPB", per = 1000,
                           repos = TRUE, verbose = FALSE)
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    save(mrs, file = mrsFile)
  }
  res = summary(mrs, length(regul))
  which(res$Regulon == "CEBPB") ## 1133
}


###### STAT3 + CEBPB genes ********
mrsFile = "../rData/103_01_mrs_GSE19114_ARACNe_3_STAT3_CEBPB.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # "STAT3_CEBPB" knockout (not found)
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "STAT3_CEBPB")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "STAT3_CEBPB", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(regul))
which(res$Regulon == "STAT3") ## 33
which(res$Regulon == "CEBPB") ## 762





## Human Acute Myeloid Leukemia context-specific ARACNe interactome
# ST486, Burkitt lymphoma cell line
source("0_functionsInThePaper.R")
data(TFs)
data(regulonlaml)
if(file.exists("103_01_regulBL.Rdata")){
  cat("load 103_01_regulBL.Rdata")
  load("103_01_regulBL.Rdata")
} else {
  cat("103_01_regulBL.Rdata")
  regulBL = standardizeRegulon(regulonlaml, filterTF = unique(TFs$TF_name))
  save(regulBL, file = "103_01_regulBL.Rdata")
}


##########  GSE2350, Burkitt lymphoma cell line, BCL6 ****
load("../rData/GSE2350.Rdata")

gse = GSE2350[, tail(colnames(GSE2350), 8)] 
tmp = pData(gse)
tmp$group = factor(rep(c("NT", "BCL6"), each = 4), levels = c("NT", "BCL6"))
pData(gse) = tmp

phenoGSE2350 = pData(gse)
dataGSE2350 = log10(exprs(gse))
geneName = gse@featureData@data$`Gene Symbol`
# geneName = gse@featureData@data$Symbol
geneNameValid = geneName[geneName != ""]
dataGSE2350 = dataGSE2350[geneName != "", ]
geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = " /// ")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE2350 = dataGSE2350[geneNameValidNoDupID,]
rownames(dataGSE2350) = geneNameValidNoMul[geneNameValidNoDupID]
dataGSE2350 = scale(dataGSE2350)


eset = ExpressionSet(assayData=dataGSE2350, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE2350))

mrsFile = "../rData/103_01_mrs_GSE2350_ARACNe_BCL6.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # BCL6 vs NT
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "BCL6")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "BCL6", per = 1000,
                         repos = TRUE, verbose = FALSE)
  
  mrs <- msviper(signature, regulBL, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(regulBL))
which(res$Regulon == "BCL6") ## 603
#


#### GSE17172
load("../rData/gse_GSE17172.Rdata")
phenoGSE17172 = pData(gse)

tmp = strsplit2(phenoGSE17172$title, "_")[,2]
phenoGSE17172$group = sub("\\d$", "", strsplit2(tmp, " ")[,1])
phenoGSE17172$normMethod = gsub("[()]", "", strsplit2(tmp, " ")[,2])
phenoGSE17172 = subset(phenoGSE17172,normMethod == "MAS5")

# Only the data of MAS5 method
idMAS5 = rownames(subset(phenoGSE17172,normMethod == "MAS5"))
gseMAS5 = gse[, idMAS5]
dataGSE17172 = exprs(gseMAS5)
geneName = gseMAS5@featureData@data$`Gene Symbol`
geneNameValid = geneName[geneName != ""]
dataGSE17172 = dataGSE17172[geneName != "", ]
geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = " /// ")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE17172 = dataGSE17172[geneNameValidNoDupID,]
rownames(dataGSE17172) = geneNameValidNoMul[geneNameValidNoDupID]
dataGSE17172 = scale(dataGSE17172)

eset = ExpressionSet(assayData=dataGSE17172, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE17172))

mrsFile = "../rData/103_01_mrs_GSE17172_ARACNe_Fox.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  # FOX vs NT
  signature0 <- rowTtest(eset, pheno = "group", group1 = "NT", group2 = "Fox")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "NT", 
                         group2 = "Fox", per = 1000,
                         repos = TRUE, verbose = FALSE)
  
  mrs <- msviper(signature, regulBL, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(regulBL))
which(res$Regulon == "FOXM1") ## 13






## IMR32 cell line, CHAF1A gene
load("../rData/GSE51978.Rdata")

gse = GSE51978

phenoGSE51978 = pData(gse)
dataGSE51978 = log10(exprs(gse))
geneName = gse@featureData@data$`Gene Symbol`
geneNameValid = geneName[geneName != ""]
dataGSE51978 = dataGSE51978[geneName != "", ]

geneNameValidNoMul = strsplit2(geneName[geneName != ""], split = "/")[,1]
geneNameValidNoDupID = !duplicated(geneNameValidNoMul)
dataGSE51978 = dataGSE51978[geneNameValidNoDupID,]
rownames(dataGSE51978) = geneNameValidNoMul[geneNameValidNoDupID]


eset = ExpressionSet(assayData=dataGSE51978, 
                     phenoData=new("AnnotatedDataFrame", data=phenoGSE51978))



#   
mrsFile = "../rData/103_01_mrs_GSE51978_ARACNe_CHAF1A_D05.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  signature0 <- rowTtest(eset, pheno = "group", group1 = "D00", group2 = "D05")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "D00", 
                         group2 = "D05", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regulBL, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(regulBL))
which(res$Regulon == "CHAF1A") ## 1319


mrsFile = "../rData/103_01_mrs_GSE51978_ARACNe_CHAF1A_D10.Rdata"
if (file.exists(mrsFile)){
  load(mrsFile)
  cat("loading mrsFile.\n")
} else {
  signature0 <- rowTtest(eset, pheno = "group", group1 = "D00", group2 = "D10")
  signature <- (qnorm(signature0$p.value/2, lower.tail = FALSE) *
                  sign(signature0$statistic))[, 1]
  nullmodel <- ttestNull(x = eset, pheno = "group", group1 = "D00", 
                         group2 = "D10", per = 1000,
                         repos = TRUE, verbose = FALSE)
  mrs <- msviper(signature, regulBL, nullmodel, verbose = FALSE)
  save(mrs, file = mrsFile)
}
res = summary(mrs, length(regulBL))
which(res$Regulon == "CHAF1A") ## 760





