## Evaluation of speed and memory
if(FALSE){
  library(funcTools)
  param1 = c(2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000,
             12000, 15000, 20000, 25000, 30000, 35000, 40000)
  param2 = c(50, 100, 200)
  
  nthread = c(8, 16)
  
  params = expand.grid(param1, param2, nthread)
  colnames(params) = c("nGenes", "nSamples", "nthread")
  
  memory = (c(4, 6, 6, 6, 6, 6, 6, 8, 10,
              10, 10, 12, 12, 14, 14))+15
  timeH0 = c(2, 2, 3, 3, 4, 5, 7, 7, 8, 
             10, 15, 18, 20, 22, 25) + 2
  timeH = c(timeH0,    2 * timeH0,  3 * timeH0, 
            ceiling(timeH0/2),  timeH0,      2*timeH0)/2
  timeH = ceiling(timeH)
  
  params = data.frame(params, memory, timeH)
  
  # params = params[c(1, 7, 46, 50, 65),]
  
  params$jobNames = paste0("TimeMemo_GRN_nthread_", 
                           params$nthread, "_", 
                           params[[1]], "_", params[[2]])
  
  
  # Parameters
  params = sortDataframe(params, "nGenes")
  
  logPath = "./19_02_log_TimeMemo/"
  if(!dir.exists(logPath)) dir.create(logPath)
  rFile = thisFile()$file
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = params$jobNames, threaded = params$nthread, 
         memoryG = params$memory, rTimeHour = params$timeH, 
         logFile = paste0(logPath, params$jobNames, ".log"), 
         email = "weiyangtao1513@gmail.com", when2Email = "aes",
         preCMD = "Rscript ", 
         param1 = params[[1]], param2 = params[[2]], param3 = params[[3]])
  qstat2("all")
}

message("Running the job.\n")

args = commandArgs(trailingOnly = TRUE)
nGenes = as.numeric(args[1])
nSamples = as.numeric(args[2])
nThread = as.numeric(args[3])


library(RegEnrich)

# constructing a RegenrichSet object
colData = data.frame(patientID = paste0('Sample_', seq(nSamples)),
                     groups = rep(c('0', '1'), each = nSamples/2),
                     row.names = paste0('Sample_', seq(nSamples)), 
                     stringsAsFactors = TRUE)
design = model.matrix(~groups, data = colData)
contrast = c(0, 1)
set.seed(123)
expr = matrix(rnorm(n=nGenes*nSamples), ncol=nSamples,
              dimnames = list(paste0('gene', seq(nGenes)), rownames(colData)))

# tmp = plotSoftPower(expr)
library(BiocParallel)
register(MulticoreParam(nThread)) # Number of cores


object = RegenrichSet(expr = expr,
                      colData = colData,
                      method = 'limma', 
                      design = design, contrast = contrast, 
                      nbTrees = 1000,
                      networkConstruction = 'GRN', minR = 0.0001,
                      enrichTest = 'FET',
                      reg = paste0('gene', seq(1572)))

## RegEnrich analysis
time = system.time(
  object <- object %>% regenrich_diffExpr() %>% 
    RegEnrich:::.regenrich_network() %>% 
    RegEnrich:::.regenrich_enrich() %>% 
    regenrich_rankScore()
)

# maximum memmory
maxMem = as.numeric(funcTools::mem_used2(maximum = T))/10^6

# How fast is the compute node
timeCtrl = system.time(
  {
    set.seed(123)
    tmpM1 = matrix(rnorm(100000), ncol = 1000)
    tmpM2 = matrix(rnorm(100000), ncol = 1000)
    for(mi in 1:500){
      cor(tmpM1, tmpM2)
    }
  }
)
# 
# res2Save = data.frame(nGenes, nSamples, nThread,
#                       method = "limma", 
#                       networkConstruction = "GRN", 
#                       enrichTest = "FET",
#                       time = time[3], 
#                       row.names = "1")

res2Save = data.frame(nGenes, nSamples, nThread,
                      time = time[3], 
                      maxMem = maxMem,
                      timeCtrl = timeCtrl[3],
                      method = "limma", 
                      networkConstruction = "GRN", 
                      enrichTest = "FET",
                      nodeName = Sys.getenv("SLURMD_NODENAME"),
                      row.names = "1")

folder = "../rData/RegEnrichSpeed"
if(!dir.exists("../rData/RegEnrichSpeed")){
  dir.create(folder)
}
# write.table(res2Save, file = paste0(folder, "/Speed_nGenes_nSamples_nThread_GRN.tsv"), 
#             append = TRUE, sep = "\t", col.names = FALSE)

writeFile = paste0(folder, "/Speed_nGenes_nSamples_nThread_GRN.tsv")
write.table(res2Save, file = writeFile, append = TRUE, sep = "\t", col.names=!file.exists(writeFile), row.names = FALSE)
