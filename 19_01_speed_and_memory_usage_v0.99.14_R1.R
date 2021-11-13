## Evaluation of speed and memory
if(FALSE){
  library(funcTools)
  # number of genes
  param1 = c(2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000,
             12000, 15000, 20000, 25000, 30000, 35000, 40000)
  # number of samples
  param2 = c(10, 20, 50, 100)
  
  # number of threads
  param3 = c(8, 16)
  memory = rep(rep(c(4, 6, 8, 10, 12, 14, 16, 18, 20,
                     24, 40, 50, 75, 100, 125), 4), 2)
  timeH = rep(rep(c(1, 1, 2, 2, 3, 4, 6, 9, 12, 
                    15, 20, 25, 30, 35, 35), 4), 2)
  params = expand.grid(param1, param2, param3)
  
  params = cbind(params, memory, timeH)
  
  # params = params[c(1, 10, 50, 65),]
  
  params$jobNames = paste0("TimeMemo_nthread_", params[[3]], "_", 
                           params[[1]], "_", params[[2]])
  
  
  # Parameters
  params = sortDataframe(params, "Var1")
  
  logPath = "./log_TimeMemo/"
  if(!dir.exists(logPath)) dir.create(logPath)
  rFile = thisFile()$file
  
  rQsub2(path = "./", 
         rFile = rFile, 
         jobName = params$jobNames, threaded = params[[3]], 
         memoryG = params$memory, rTimeHour = params$timeH, 
         logFile = paste0(logPath, params$jobNames, ".log"), 
         email = "w.tao-2@umcutrecht.nl", when2Email = "aes",
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

tmp = plotSoftPower(expr)
library(BiocParallel)
register(MulticoreParam(nThread)) # Number of cores

object = RegenrichSet(expr = expr,
                      colData = colData,
                      method = 'limma', 
                      design = design, contrast = contrast,
                      networkConstruction = 'COEN', softPower = 6,
                      enrichTest = 'FET',
                      reg = paste0('gene', seq(1572)))

## RegEnrich analysis
time = system.time(
  object <- object %>% regenrich_diffExpr() %>% 
    regenrich_network() %>% 
    regenrich_enrich() %>% 
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

res2Save = data.frame(nGenes, nSamples, nThread,
                      time = time[3], 
                      maxMem = maxMem,
                      timeCtrl = timeCtrl[3],
                      method = "limma", 
                      networkConstruction = "COEN", 
                      enrichTest = "FET",
                      nodeName = Sys.getenv("SLURMD_NODENAME"),
                      row.names = "1")
folder = "../rData/RegEnrichSpeed"
if(!dir.exists("../rData/RegEnrichSpeed")){
  dir.create(folder)
}
writeFile = paste0(folder, "/Speed_nGenes_nSamples_nThread_COEN.tsv")
write.table(res2Save, file = writeFile, append = TRUE, sep = "\t", col.names=!file.exists(writeFile), row.names = FALSE)

