

########## Speed and memory usage of RegEnrich
library(funcTools)
source("./0_functionsInThePaper.R")

# Combine data for both networks
if(file.exists("../rData/19_03_RegEnrichSpeedMemory_R1.Rdata")){
    load("../rData/19_03_RegEnrichSpeedMemory_R1.Rdata")
} else {
    speedMem = list()
    for(network in c("COEN", "GRN")){
        file =  paste0("../rData/RegEnrichSpeed/Speed_nGenes_nSamples_nThread_", 
                       network, ".tsv")
        
        figFolder = "../fig/"
        speedCOEN0 = read.csv(file, sep = "\t")
        speedCOEN0 = subset(speedCOEN0, nGenes != "nGenes")
        for(mi in colnames(speedCOEN0)[1:6]){
            speedCOEN0[[mi]] = as.numeric(as.character(speedCOEN0[[mi]]))
        }
        speedCOEN0 = sortDataframe(speedCOEN0, colnames(speedCOEN0)[c(1, 3, 2)])
        speedCOEN0 = speedCOEN0[!duplicated(speedCOEN0[,1:3]),]
        speedMem[[network]] = speedCOEN0
    }
    speedCOEN = do.call("rbind", speedMem)
    save(speedCOEN, file = "../rData/19_03_RegEnrichSpeedMemory_R1.Rdata")
}

# Different color scale for each panel
for(network in c("COEN", "GRN")){
    ## Speed of COEN
    if(!file.exists(paste0("../fig/19_03_heatmap_", network, "_time_nThread4.png"))){
        ## use memory to plot heatmap
        library(RColorBrewer)
        colorBoard = rev(brewer.pal(n = 7, name ="RdYlBu"))
        {
            # plot the same scale
            breaksList = seq(min(subset(speedCOEN, networkConstruction == network)$maxMem), 
                             max(subset(speedCOEN, networkConstruction == network)$maxMem) + 499, 
                             by = 500)
            color = colorRampPalette(colorBoard)(length(breaksList))
            for(nT in as.character(c(1, 4, 8, 16))){
                pmem1 = plotSpeedMemory2(subset(speedCOEN, networkConstruction == network), 
                                         what = "maxMem", nThread = nT, 
                                         main = paste0("Memory Usage (Mb), NCore = ", nT),
                                         color = color, breaks = breaksList, 
                                         cellwidth = 25, cellheight = 15, 
                                         tableFileName = paste0("../rData/memoryUsage_",network, "_NCore", nT, ".csv"))
                plotSave(filename = paste0(figFolder, "19_03_heatmap_", network, 
                                           "_mem_nThread", nT, ".svg"), 
                         Plot = pmem1, width = 4, height = 8, dpi = 300)
                plotSave(filename = paste0(figFolder, "19_03_heatmap_", network, 
                                           "_mem_nThread", nT, ".png"), 
                         Plot = pmem1, width = 4, height = 8, dpi = 300)
            }
        }
        
        ## use time to plot heatmap
        colorBoard = rev(brewer.pal(n = 7, name ="PiYG"))
        {
            # breaksList = seq(1.5, 4.5, by = 0.05)
            breaksList = seq(min(log10(subset(speedCOEN, networkConstruction == network)$time)), 
                             max(log10(subset(speedCOEN, networkConstruction == network)$time)), 
                             by = 0.05)
            color = colorRampPalette(colorBoard)(length(breaksList))
            # nThread = 1, 4, 8, 16
            for(nT in as.character(c(1, 4, 8, 16))){
                ptime1 = plotSpeedMemory2(subset(speedCOEN, networkConstruction == network), 
                                          what = "time", log10Transform = T, 
                                          nThread = nT, main = paste0("log10(Time) (s), NCore = ", nT),
                                          color = color, breaks = breaksList, 
                                          cellwidth = 25, cellheight = 15,
                                          tableFileName = paste0("../rData/timeUsage_",network, "_NCore", nT, ".csv"))
                plotSave(filename = paste0(figFolder, "19_03_heatmap_", network, 
                                           "_time_nThread", nT, ".svg"), 
                         Plot = ptime1, width = 4, height = 8, dpi = 300)
                plotSave(filename = paste0(figFolder, "19_03_heatmap_", network, 
                                           "_time_nThread", nT, ".png"), 
                         Plot = ptime1, width = 4, height = 8, dpi = 300)
            }
        }
    }
}



# Same color scale for each panel
for(nThreads in c(1, 4, 8, 16)){
    {
        ## use memory to plot heatmap
        library(RColorBrewer)
        colorBoard = rev(brewer.pal(n = 7, name ="RdYlBu"))
        {
            speedCOEN1 = subset(speedCOEN, nThread == nThreads)
            speedCOEN1$nSamples = paste0(speedCOEN1$networkConstruction, speedCOEN1$nSamples)
            speedCOEN1 = sortDataframe(speedCOEN1, by = "nGenes")
            speedCOEN1$maxMem = speedCOEN1$maxMem/1000
            speedCOEN1$nGenes = speedCOEN1$nGenes/1000
            
            # plot the same scale
            breaksList = seq(min(speedCOEN1$maxMem), 
                             max(speedCOEN1$maxMem), 
                             by = 5)
            color = colorRampPalette(colorBoard)(length(breaksList))
            # Memory
            subset(speedCOEN, nThread == nThreads)
            pmem1 = plotSpeedMemory(speedCOEN1, 
                                    what = "maxMem", nThread = nThreads, 
                                    main = paste0("Memory Usage (Mb), NCore = ", nThreads),
                                    color = color, breaks = breaksList, 
                                    cellwidth = 15, cellheight = 9)
            plotSave(filename = paste0("../fig/19_03_heatmap_bothNetworks_memory_nThread", 
                                       nThreads, ".svg"), 
                     Plot = pmem1, width = 6, height = 8, dpi = 300)
            plotSave(filename = paste0("../fig/19_03_heatmap_bothNetworks_memory_nThread", 
                                       nThreads, ".png"), 
                     Plot = pmem1, width = 6, height = 8, dpi = 300)
        }
        
        ## use time to plot heatmap
        colorBoard = rev(brewer.pal(n = 7, name ="PiYG"))
        {
            # breaksList = seq(1.5, 4.5, by = 0.05)
            breaksList = seq(min(log10(speedCOEN1$time)), 
                             max(log10(speedCOEN1$time)), 
                             by = 0.05)
            color = colorRampPalette(colorBoard)(length(breaksList))
            # nThread = 1
            ptime1 = plotSpeedMemory(speedCOEN1, 
                                     what = "time", log10Transform = T, 
                                     nThread = nThreads, 
                                     main = paste0("log10(Time) (s), NCore = ", nThreads),
                                     color = color, breaks = breaksList, 
                                     cellwidth = 15, cellheight = 9)
            plotSave(filename = paste0("../fig/19_03_heatmap_bothNetworks_time_nThread", 
                                       nThreads, ".svg"), 
                     Plot = ptime1, width = 6, height = 8, dpi = 300)
            plotSave(filename = paste0("../fig/19_03_heatmap_bothNetworks_time_nThread", 
                                       nThreads, ".png"), 
                     Plot = ptime1, width = 6, height = 8, dpi = 300)
        }
    }
}





