
library(funcTools)
library(reshape2)
library(RegEnrich)
library(microbenchmark)

library(doParallel)
library(foreach)
library(ggplot2)
library(DOSE)
library(WGCNA)
library(igraph)
library(GGally)
library(ARACNe)


# library(VennDiagram)
library(viper)
# library(bcellViper)
library(DESeq2)

library(affy)
library(GEOquery)
library(limma)


#' plot the time for each process
#' @param usedTime load("../rData/usedTimes_compareDifferentNetworkConstructionMethods2.Rdata")
#' @param oldNames the original process names.
#' @param newNames the new names for the processes, 
#' this name must match the oldNames.
compareTimePlot = function(usedTime, 
                           oldNames = c("coeNet", "grn1", "grn2"),
                           newNames = c("WGCNA", "parallel_RF", "non_parallel_RF")){
  # load("../rData/usedTimes_compareDifferentNetworkConstructionMethods2.Rdata")
  x = lapply(usedTime, function(x){ 
    x = as.data.frame(x, row.names = as.character(x[,1]))
    y = x[oldNames, 2, drop = F]
    rownames(y) = newNames
    y
  })
  # nm = setNames(newNames, oldNames)
  y = do.call(cbind, x)/10^9 # seconds
  colnames(y) = names(x)
  y$Method = rownames(y)
  z = melt(y, id.vars = "Method")
  z$Method = factor(z$Method, levels = newNames)
  g = ggplot(z, aes(x = variable, y = value, fill = Method)) + 
    geom_bar(stat = "identity", position = "dodge2") + 
    theme_bw()+ 
    # scale_y_log10(breaks = pretty(z$value, n = 6))+
    scale_y_log10(breaks = 10^pretty(log10(range(z$value)), n = 6))+
    ylab("Seconds")+
    xlab("Numbers of genes (Ng) and samples (Ns)")+
    theme_Publication(base_size = 12, x_angle = 90) + 
    scale_fill_Publication()
  return(g)
}


#' calculate the hubness based on different measurement
#' @param netEdge the data frame of network edge information.
#' @param reg the node need to evaluate
#' @return a data frame of degrees (in, out and all), 
#' closeness (in, out and all), and betweeness.
hubness = function(netEdge = topNet5$netTable,
                   reg = unique(topNet5$netTable[[1]])){
  # igraph object
  coeNetTop5 = graph.data.frame(netEdge, directed = TRUE)
  
  # degree
  cat("Calculating degree\n")
  inDegree = degree(coeNetTop5, v = reg, mode = "in")
  outDegree = degree(coeNetTop5, v = reg, mode = "out")
  allDegree = degree(coeNetTop5, v = reg, mode = "all")
  
  ####### centrality (closeness) ########
  cat("Calculating closeness\n")
  inCloseness = closeness(coeNetTop5, vids = reg, mode='in')
  outCloseness = closeness(coeNetTop5, vids = reg, mode='out')
  allCloseness = closeness(coeNetTop5, vids = reg, mode='all')
  
  ####### centrality (betweenness) ########
  cat("Calculating betweenness\n")
  Betweenness = betweenness(coeNetTop5, v = reg)
  
  ####### centrality (eigenvector centrality) ########
  # ev_obj_coeNetTop5 = evcent(coeNetTop5)
  # eigenvector = ev_obj_coeNetTop5$vector[reg]
  central_coeNetTop5 = data.frame(reg, inDegree, outDegree,
                                  inCloseness, outCloseness, 
                                  Betweenness, allDegree, allCloseness)
  return(central_coeNetTop5)
}



#' plot the overlaps of hub genes
topHubVenn = function(hubness, colID = 2:3, topN = 100, topReg = NULL,
                      fill = c("blue", "orange3", "green", "red", "magenta")){
  if (is.null(topReg)){
    nm = rownames(hubness)
    topReg = lapply(hubness[, colID], function(x){
      nm[order(x, decreasing = TRUE)[1:topN]]
    })
  }
  len = length(topReg)
  if (len == 2){
    catPos = c(340, 20)
    catDist = 0.05
  } else if (len == 3){
    # catPos = c(30, 180, 330)
    catPos = c(-20, 20, 180)
    catDist = c(0.1, 0.1, 0.1)
  } else if (len == 4){
    catPos = c(340, 20, 0, 0)
    catDist = c(0.22, 0.22, 0.125, 0.1)
  } else if (len == 5){
    catPos = c(0, 330, 210, 150, 0)
    catDist = c(0.25, 0.25, 0.25, 0.25, 0.25)
  } else {
    stop("Unproper plotID!")
  }
  # all degree and all closeness
  x = venn.diagram(x = topReg, 
                   filename = NULL, 
                   fill = fill[1:len], 
                   cat.col = fill[1:len],
                   alpha = 0.4, cex = 1.8, 
                   cat.pos = catPos, cat.dist = catDist, 
                   cat.cex = 2, margin = 0.1)
  grid.newpage()
  grid.draw(x)
  invisible(x)
}


#' plot regulator and targets expression
#' @param reg regulator
#' @param expr a matrix or data frame, the gene expression data.
#' @param pFC a data frame of gene, p, and log2FC.
#' @param topNet a `topNetwork` object, which shows the sub 
#' network (top edges) of co-expression network.
#' @param n the maximun number of targets to plot.
#' @param tarCol the color of the line for the targets of "reg".
#' @param regCol the color of the line for the "reg".
#' @param xlab x label of plot.
#' @param ylab y label of plot.
#' @param ... other parameters in ggplot function.
#' @return a ggplot object.
#' @export
#' @examples 
#' \dontrun{
#' plotRegTarExpr(reg = "CREG1", expr = log2FPKMhi, pFC = pFC, 
#'                topNet = topNet5, n = 1000, scale = TRUE, 
#'                tarCol = alpha("black", 0.05),
#'                regCol = "#ffaa00", 
#'                xlab = "Samples",
#'                ylab = "Z-scores")
#' }
plotRegTarExpr = function(reg, expr, pFC, topNet, n = 100, 
                          scale = TRUE, 
                          tarCol = alpha("black", 0.1), 
                          regCol = "#ffaa00", 
                          xlab = "Samples", ylab = "Z-scores", 
                          p_threshold = 0.05, ...){
  
  pReg = subset(pFC[topNet$net[[reg]],], p < p_threshold)
  gNames = rownames(sortDataframe(pReg, "p"))[1: min(n, nrow(pReg))]
  
  expr_regTar = expr[c(reg, gNames),]
  
  if (scale){
    expr_regTar = t(scale(t(expr_regTar)))
  }
  
  expr_regTar2 = melt(data.frame(expr_regTar, 
                                 gene = factor(rownames(expr_regTar), 
                                               rownames(expr_regTar))),
                      id.vars = "gene")
  expr_regTar2 = sortDataframe(expr_regTar2, c("gene", "variable"))
  p = ggplot(data = expr_regTar2, 
             aes(x = variable, y = value, 
                 group = gene, color = I(tarCol)), ...) +
    geom_line(show.legend = F) +
    geom_line(aes(color = I(regCol)), 
              data = subset(expr_regTar2, gene == reg)) +
    xlab(xlab) + ylab(ylab) + ggtitle(reg)+
    theme_Publication(x_angle = 90)
  return(p)
}

getRegTarExpr = function(reg, expr, pFC, topNet, n = 100, 
                          scale = TRUE, 
                          tarCol = alpha("black", 0.1), 
                          regCol = "#ffaa00", 
                          xlab = "Samples", ylab = "Z-scores", 
                          p_threshold = 0.05, ...){
  
  pReg = subset(pFC[topNet$net[[reg]],], p < p_threshold)
  gNames = rownames(sortDataframe(pReg, "p"))[1: min(n, nrow(pReg))]
  
  expr_regTar = expr[c(reg, gNames),]
  
  if (scale){
    expr_regTar = t(scale(t(expr_regTar)))
  }
  
  expr_regTar2 = data.frame(expr_regTar, 
                            gene = factor(rownames(expr_regTar), 
                                          rownames(expr_regTar)))
  return(expr_regTar2)
}

#' Fisher exact test for two sets in Venn diagram
#' @param x0 the number of intersection of two sets
#' @param x1 the number of elements only in set 1
#' @param x2 the number of elements only in set 2
#' @param x the number of global elements, including set 1, set 2 
#' and their complmentary set.
#' @examples {
#' @dontrun {
#' viperAlone = 30
#' RegEnrichAlone = 30
#' overlap = 30
#' total = 1200
#' VennFisherTest(overlap, viperAlone, RegEnrichAlone, total)
#' }
#' }
VennFisherTest <- function(x0, x1, x2, x) {
  fisher.test(matrix(c(x0,x1,x2,x-x1-x2-x0), ncol = 2), 
              alternative = "greater")
}


library(biomaRt)
ENSGID2GeneName = function(ENSGID, oneByOne = T, sameOrder = F){
  host38.79 = "mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = 'ensembl_gene_id'
  attributes0 = c('ensembl_gene_id','external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = ENSGID, mart = mart)
  if (!oneByOne) {
    uniqID = unique(ID$ensembl_gene_id)
    geneName = sapply(uniqID, function(x){paste(ID[(ID[,1] %in% x), 2], collapse = ", ")})
    
    ID = data.frame(ENSGID = ENSGID, 
                    geneName = geneName[ENSGID])
  }
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[ENSGID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}


ENSGID2GeneName37 = function(ENSGID, oneByOne = T, sameOrder = F){
  host38.79 = "grch37.ensembl.org/"#"mar2015.archive.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = 'ensembl_gene_id'
  attributes0 = c('ensembl_gene_id','external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = ENSGID, mart = mart)
  if (!oneByOne) {
    uniqID = unique(ID$ensembl_gene_id)
    geneName = sapply(uniqID, function(x){paste(ID[(ID[,1] %in% x), 2], collapse = ", ")})
    
    ID = data.frame(ENSGID = ENSGID, 
                    geneName = geneName[ENSGID])
  }
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[ENSGID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}


ENSGID2GeneName2019 = function(ENSGID, oneByOne = T, sameOrder = F){
  host38.79 = "apr2019.archive.ensembl.org" # "grch37.ensembl.org/"#
  # host38.79 = "http://www.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = 'ensembl_gene_id'
  attributes0 = c('ensembl_gene_id','external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = ENSGID, mart = mart)
  if (!oneByOne) {
    uniqID = unique(ID$ensembl_gene_id)
    geneName = sapply(uniqID, function(x){paste(ID[(ID[,1] %in% x), 2], collapse = ", ")})
    
    ID = data.frame(ENSGID = ENSGID, 
                    geneName = geneName[ENSGID])
  }
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[ENSGID,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}


Entrez2GeneName2019 = function(Entrez, oneByOne = T, sameOrder = F){
  host38.79 = "apr2019.archive.ensembl.org" # "grch37.ensembl.org/"#
  # host38.79 = "http://www.ensembl.org"
  biomart0 = "ENSEMBL_MART_ENSEMBL"
  dataset0 = "hsapiens_gene_ensembl"
  mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
  filters0 = 'entrezgene'
  attributes0 = c('entrezgene','external_gene_name')
  ID = getBM(filters = filters0, attributes = attributes0, values = Entrez, mart = mart)
  if (!oneByOne) {
    uniqID = unique(ID$entrezgene)
    geneName = sapply(setNames(uniqID, uniqID), function(x){paste(ID[(ID[,1] %in% x), 2], collapse = ", ")})
    
    ID = data.frame(Entrez = Entrez, 
                    geneName = geneName[Entrez])
  }
  if(sameOrder) {
    rownames(ID) = ID[,1]
    ID = ID[Entrez,]
    rownames(ID) = 1:nrow(ID)
  }
  return(ID)
}



RefSeq2GeneName = function(RefSeq, oneByOne = T, sameOrder = F){
    host38.79 = "mar2015.archive.ensembl.org"
    biomart0 = "ENSEMBL_MART_ENSEMBL"
    dataset0 = "hsapiens_gene_ensembl"
    mart = useMart(biomart = biomart0, dataset = dataset0, host = host38.79)
    filters0 = 'refseq_mrna'
    attributes0 = c('refseq_mrna','external_gene_name')
    ID = getBM(filters = filters0, attributes = attributes0, values = RefSeq, mart = mart)
    if (!oneByOne) {
        uniqID = unique(ID$refseq_mrna)
        geneName = sapply(uniqID, function(x){paste(ID[(ID[,1] %in% x), 2], collapse = ", ")})
        
        ID = data.frame(RefSeq = RefSeq, 
                        geneName = geneName[RefSeq])
    }
    if(sameOrder) {
        rownames(ID) = ID[,1]
        ID = ID[RefSeq,]
        rownames(ID) = 1:nrow(ID)
    }
    return(ID)
}


library(gdata)

TFfamily = list(ETS = c("ELF1", "ELF2", "ELF4", # ELF
                        "GABPA", # ELG
                        "ERG", "FLI1", "FEV", # ERG
                        "ERF", "ETV3", # ERF
                        "ELF3", "ELF5", "ESE3", # ESE
                        "ETS1", "ETS2", # ETS
                        "SPDEF", # PDEF
                        "ETV4", "ETV5", "ETV1", # PEA3
                        "ETV2", # ER71
                        "SPI1", "SPIB", "SPIC", # SPI
                        "ELK1", "ELK4", "ELK3", # TCF
                        "ETV6", "ETV7"), # TEL
                IRF = paste0("IRF", 1:9),
                STAT = paste0("STAT", c(1:4, "5A", "5B", 6)),
                NFKB = c("NFKB1", "NFKB2", "RELA", "RELB", "REL"))



### plot heatmap of RegEnrich speed and memory usage
library(pheatmap)
#' @param speedMemoData is the data that read from the files generaged 
#' by 19_01 and 19_02 R file. 
#' @param what is the column to plot as the colored cells in the heatmap.
#' The default is "maxMem", which means memory usage.
#' @param log10Transform logical. Whether or not to log10 transform 
#' the `what` column in speedMemoData.
#' @param nThread how many threads (CPUs) was used to generate the 
#' speedMemoData.
#' @param ... other parameters in pheatmap function.
plotSpeedMemory = function(speedMemoData, what = c("maxMem", "time"), 
                           log10Transform = FALSE, nThread = "1", ...){
    what = match.arg(what)
    nThread = as.character(nThread)
    hm_mem = tapply(1:nrow(speedMemoData), speedMemoData[["nThread"]], FUN = function(i){
        x = if(log10Transform) log10(speedMemoData[i, what]) else speedMemoData[i, what]
        matrix(x, nrow = length(unique(speedMemoData[["nGenes"]])), 
               byrow = T, dimnames = list(unique(speedMemoData[["nGenes"]]), 
                                          unique(speedMemoData[["nSamples"]])))
    })
    pm = pheatmap(hm_mem[[nThread]][nrow(hm_mem[[nThread]]):1,], 
                  cluster_rows = F, cluster_cols = F, ...)
    return(invisible(pm))
}

plotSpeedMemory2 = function(speedMemoData, what = c("maxMem", "time"), 
                           log10Transform = FALSE, nThread = "1", 
                           tableFileName = NULL, ...){
    what = match.arg(what)
    nThread = as.character(nThread)
    hm_mem = tapply(1:nrow(speedMemoData), speedMemoData[["nThread"]], FUN = function(i){
        x = if(log10Transform) log10(speedMemoData[i, what]) else speedMemoData[i, what]
        matrix(x, nrow = length(unique(speedMemoData[["nGenes"]])), 
               byrow = T, dimnames = list(unique(speedMemoData[["nGenes"]]), 
                                          unique(speedMemoData[["nSamples"]])))
    })
    if(is.null(tableFileName)){
        print(hm_mem[[nThread]][nrow(hm_mem[[nThread]]):1,])
    } else {
        print(tableFileName)
        write.csv(hm_mem[[nThread]][nrow(hm_mem[[nThread]]):1,], file = tableFileName)
    }
    pm = pheatmap(hm_mem[[nThread]][nrow(hm_mem[[nThread]]):1,], 
                  cluster_rows = F, cluster_cols = F, ...)
    return(invisible(pm))
}


hpcCPUinfo = function(){
    x = system("cat /proc/cpuinfo", T)
    y = gsub("\t", "", removeSpace(x))
    head(y)
    id2 = which(y == "")
    id1 = c(1, head(id2, -1)+1)
    
    di = diff(id2)
    nrow = max(di)
    z = strSplit(y, ":")
    
    if (length(unique(di)) == 1){
        res = matrix(z[,2], nrow = nrow, 
                     dimnames = list(unique(z[,1]), NULL))
    } else {
        res = matrix(NA, nrow = nrow, ncol = 1, 
                     dimnames = list(unique(z[,1]), NULL))
        for (mi in seq_along(id1)){
            res = cbind(res, z[id1[mi]:id2[mi], 2])
        }
    }
    res = t(res)
    return(res)
}


standardizeRegulon = function(regulon, filterTF = unique(TFs$TF_name)){
  
  ### Map gene entrez ID to gene name
  allID = unique(c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)))))
  idMap = Entrez2GeneName2019(allID, oneByOne = T)
  idMap = subset(idMap, !duplicated(idMap$entrezgene))
  rownames(idMap) = idMap$entrezgene
  
  # Map tf name 
  tf = names(regulon)
  tf2 = idMap[tf, "external_gene_name"]
  notNaID = which(!is.na(tf2))
  regul0 = regulon[notNaID]
  names(regul0) = tf2[notNaID]
  
  # Map target name
  regul1 = lapply(regul0, function(x){
    tar = names(x$tfmode)
    tar2 = idMap[tar, "external_gene_name"]
    notNaID = which(!is.na(tar2))
    
    tfmode = x$tfmode[notNaID]
    names(tfmode) = tar2[notNaID]
    lh = x$likelihood[notNaID]
    
    list(tfmode = tfmode, likelihood = lh)
  })
  
  
  # filtered regulon
  regul = regul1[names(regul1) %in% filterTF]
  return(regul)
}