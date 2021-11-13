ttestNull1 = function (x, ...) 
{
  .local <- function (x, pheno, group1, group2, per = 1000, 
                      repos = TRUE, seed = 1, verbose = TRUE) 
  {
    pos1 <- pData(x)[[pheno]] %in% group1
    pos2 <- pData(x)[[pheno]] %in% group2
    if (length(pos1) == 0) 
      stop(paste(pheno, " is not present in the ExpressionSet object", 
                 sep = ""), call. = FALSE)
    if (length(which(pos1)) == 0) 
      stop(paste(group1, " was nor found in ", pheno, 
                 sep = ""), call. = FALSE)
    if (length(which(pos2)) == 0) 
      stop(paste(group2, " was nor found in ", pheno, 
                 sep = ""), call. = FALSE)
    ttestNull0(exprs(x)[, pos1], exprs(x)[, pos2], per = per, 
              repos = repos, seed = seed, verbose = verbose)
  }
  .local(x, ...)
}

rowVars = viper:::rowVars

ttestNull0 = function (x, ...) 
{
  .local <- function (x, y, per = 1000, repos = TRUE, seed = 1, 
                      cores = 1, verbose = TRUE) 
  {
    # if (seed > 0) 
    #   set.seed(round(seed))
    pb <- NULL
    if (verbose) {
      message(date(), "\nComputing the null model distribution by ", 
              per, " permutations.")
    }
    if (cores == 1){
      if (verbose) 
        pb <- txtProgressBar(max = per, style = 3)
      res <- sapply(1:per, function(i, x, y, repos, pb, 
                                    verbose) {
        if (verbose) 
          setTxtProgressBar(pb, i)
        expset <- cbind(x, y)
        repeat {
          sorder <- sample(ncol(expset), replace = repos)
          if (length(unique(sorder[1:ncol(x)])) > 1 & 
              length(unique(sorder[-(1:ncol(x))])) > 1) 
            break
          if (verbose) 
            message("-", appendLF = FALSE)
        }
        x1 <- filterColMatrix(expset, sorder[1:ncol(x)])
        y1 <- filterColMatrix(expset, sorder[-(1:ncol(x))])
        largo <- rowSums(!is.na(x1))
        largoy <- rowSums(!is.na(y1))
        t <- ((rowMeans(x1, na.rm = TRUE) - rowMeans(y1, 
                                                     na.rm = TRUE))/sqrt(((largo - 1) * rowVars(x1) + 
                                                                            (largoy - 1) * rowVars(y1))/(largo + largoy - 
                                                                                                           2))/sqrt(1/largo + 1/largoy))[, 1]
        t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), 
                   lower.tail = FALSE) * sign(t)
        names(t) <- rownames(x)
        return(t)
      }, x = x, y = y, repos = repos, pb = pb, verbose = verbose)
    }
    colnames(res) <- 1:per
    if (verbose) 
      message("\n", date())
    return(res)
  }
  .local(x, ...)
}



# Do t test for each row 
f1 = function(exprs, group1, group2, x1 = NULL, y1 = NULL, paired = F){
  if (is.null(x1)){
    t_p = sapply(setNames(1:nrow(exprs), rownames(exprs)), function(x){
      tmp = t.test(x = exprs[x,group1], y = exprs[x,group2], 
                   paired = paired, var.equal = TRUE)
      return(c(statistics = tmp$statistic, p.value = tmp$p.value))
    })
  } else {
    t_p = sapply(setNames(1:nrow(x1), rownames(x1)), function(x){
      tmp = t.test(x = x1[x,], y = y1[x,], paired = paired, var.equal = TRUE)
      return(c(statistics = tmp$statistic, p.value = tmp$p.value))
    })
  }
  t_p = t(t_p)
  # change the p.values of t.test to the quantiles of standard normal distribution 
  signature2 <- (qnorm(t_p[,2]/2, lower.tail = FALSE) * sign(t_p[,1]))
  
  return(signature2)
}


# Do permulation first and do t test for each row 
f2 = function(exprs, group1, group2, repos = T, per = 100, paired = F){
  res = sapply(1:per, function(i){
    print(i)
    repeat {
      sorder <- sample(ncol(exprs), replace = repos)
      if (length(unique(sorder[1:length(group1)])) > 1 & 
          length(unique(sorder[-(1:length(group1))])) > 1) {
        break
      }
    }
    x1 <- filterColMatrix(exprs, sorder[1:length(group1)])
    y1 <- filterColMatrix(exprs, sorder[-(1:length(group1))])
    f1(x1 = x1, y1 = y1)
  })
  colnames(res) = 1:ncol(res)
  return(res)
}

