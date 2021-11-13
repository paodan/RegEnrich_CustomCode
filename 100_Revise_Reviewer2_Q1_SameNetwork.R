nTF = 500
nTF_top = 20
nSample = 100

# expr, pData
load("../rData/100_expr_pData_simulate2.Rdata")

# object
load("../rData/100_object_rank_2.Rdata")
network_user = results_topNet(object)
save(object, file = "../rData/100_object_rank_100samples_3.Rdata")


# Enrichment analysis by GSEA
set.seed(1234)
object = regenrich_enrich(object, enrichTest = "GSEA")
# Regulators ranking
object = regenrich_rankScore(object)
save(object, file = "../rData/100_object_rank_100samples_3_GSEA.Rdata")


## 50
load("../rData/100_object_rank_50samples_2.Rdata")
object50_2 = object50
grep("top", object50_2@resScore$reg)

regenrich_network(object50) = network_user
# Enrichment analysis by Fisher's exact test (FET)
object50 = regenrich_enrich(object50)
# Regulators ranking
object50 = regenrich_rankScore(object50)
grep("top", object50@resScore$reg)
save(object50, file = "../rData/100_object_rank_50samples_3.Rdata")

# Enrichment analysis by GSEA
set.seed(1234)
object50 = regenrich_enrich(object50, enrichTest = "GSEA")
# Regulators ranking
object50 = regenrich_rankScore(object50)
save(object50, file = "../rData/100_object_rank_50samples_3_GSEA.Rdata")


## 20
load("../rData/100_object_rank_20samples_2.Rdata")
object20_2 = object20
grep("top", object20_2@resScore$reg)

regenrich_network(object20) = network_user
# Enrichment analysis by Fisher's exact test (FET)
object20 = regenrich_enrich(object20)
# Regulators ranking
object20 = regenrich_rankScore(object20)
grep("top", object20@resScore$reg)
save(object20, file = "../rData/100_object_rank_20samples_3.Rdata")

# Enrichment analysis by GSEA
set.seed(1234)
object20 = regenrich_enrich(object20, enrichTest = "GSEA")
# Regulators ranking
object20 = regenrich_rankScore(object20)
save(object20, file = "../rData/100_object_rank_20samples_3_GSEA.Rdata")


## 10
load("../rData/100_object_rank_10samples_2.Rdata")
object10_2 = object10
grep("top", object10_2@resScore$reg)

regenrich_network(object10) = network_user
# Enrichment analysis by Fisher's exact test (FET)
object10 = regenrich_enrich(object10)
# Regulators ranking
object10 = regenrich_rankScore(object10)
grep("top", object10@resScore$reg)
save(object10, file = "../rData/100_object_rank_10samples_3.Rdata")

# Enrichment analysis by GSEA
set.seed(1234)
object10 = regenrich_enrich(object10, enrichTest = "GSEA")
# Regulators ranking
object10 = regenrich_rankScore(object10)
save(object10, file = "../rData/100_object_rank_10samples_3_GSEA.Rdata")


## 6
load("../rData/100_object_rank_6samples_2.Rdata")
object6_2 = object6
grep("top", object6_2@resScore$reg)

regenrich_network(object6) = network_user
# Enrichment analysis by Fisher's exact test (FET)
object6 = regenrich_enrich(object6)
# Regulators ranking
object6 = regenrich_rankScore(object6)
grep("top", object6@resScore$reg)
save(object6, file = "../rData/100_object_rank_6samples_3.Rdata")

# Enrichment analysis by GSEA
set.seed(1234)
object6 = regenrich_enrich(object6, enrichTest = "GSEA")
# Regulators ranking
object6 = regenrich_rankScore(object6)
save(object6, file = "../rData/100_object_rank_6samples_3_GSEA.Rdata")


### Plot COEN + FET
fileNames = c("../rData/100_object_rank_100samples_3.Rdata", 
              "../rData/100_object_rank_50samples_3.Rdata", 
              "../rData/100_object_rank_20samples_3.Rdata", 
              "../rData/100_object_rank_10samples_3.Rdata",
              "../rData/100_object_rank_6samples_3.Rdata")
library(ggplot2)
library(fgsea)
p = list()
for(mi in fileNames){
  data = get(load(mi))
  sudoPathway = list(topReg = grep("top", data@resScore$reg, value = T))
  sudoRanks = setNames((length(data@resScore$reg):1)/500 - 0.5, data@resScore$reg)
  enrich = fgsea(sudoPathway, sudoRanks)
  
  nn = substr(strSplit(mi, "\\.R")[,1], 21, 100)
  if(nn == "rank"){
    nn = "rank_100samples"
  } else if(nn == "rank_2"){
    nn = "rank_100samples_2"
  }
  nn = gsub("sam", " sam", nn)
  p[[nn]] = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(strSplit(nn, "_")[,2], " (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
}

for(nn in names(p)){
  print(nn)
  figFile = paste0("../fig/100_keyTF_", nn, ".png")
  ggsave(filename = figFile, p[[nn]], width = 3, height = 1.8)
}



### Plot COEN + GSEA
fileNames = c("../rData/100_object_rank_100samples_3_GSEA.Rdata", 
              "../rData/100_object_rank_50samples_3_GSEA.Rdata", 
              "../rData/100_object_rank_20samples_3_GSEA.Rdata", 
              "../rData/100_object_rank_10samples_3_GSEA.Rdata",
              "../rData/100_object_rank_6samples_3_GSEA.Rdata")
library(ggplot2)
library(fgsea)
p = list()
for(mi in fileNames){
  data = get(load(mi))
  sudoPathway = list(topReg = grep("top", data@resScore$reg, value = T))
  sudoRanks = setNames((length(data@resScore$reg):1)/500 - 0.5, data@resScore$reg)
  enrich = fgsea(sudoPathway, sudoRanks)
  
  nn = substr(strSplit(mi, "\\.R")[,1], 21, 100)
  if(nn == "rank"){
    nn = "rank_100samples"
  } else if(nn == "rank_2"){
    nn = "rank_100samples_2"
  }
  nn = gsub("sam", " sam", nn)
  p[[nn]] = plotEnrichment(sudoPathway[[1]], stats = sudoRanks, ) +
    ggtitle(paste0(strSplit(nn, "_")[,2], " (p = ", signif(enrich$pval, 2), ")"))+
    theme(plot.title = element_text(hjust = 0.5))
}

for(nn in names(p)){
  print(nn)
  figFile = paste0("../fig/100_keyTF_", nn, ".png")
  ggsave(filename = figFile, p[[nn]], width = 3, height = 1.8)
}
