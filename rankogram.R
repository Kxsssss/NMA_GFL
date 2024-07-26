library(meta)
library(netmeta)
library(metarep)

# Function to ensure rankings of all simulations from joint distribution of relative effects
# will include all levels (i.e., 1:ntrts)
# X is a ntrts x nsim matrix of the treatments ranks for each simulation (e.g., rnk in ranksampling_plus)
# ntrts is the number of treatments in the network (e.g., x$n)
rank.tab.list <- function(X, ntrts){
  
  tt <- list()
  for(i in 1:dim(X)[1]){
    tt[[i]] <- factor(X[i,], levels = 1:ntrts)
  }
  
  #names(tt) <- 1:ntrts
  
  lapply(tt, table)
  
}

# Updated ranksamping() to solve problems of same level ranking 
# and rank of only two reatments.
ranksampling_plus <- function(x, nsim,
                              pooled = "random", small.values = "desirable",
                              keep.sample = FALSE) {
  chkclass(x, "netmeta")
  pooled <- setchar(pooled, c("common", "random"))
  small.values <- setsv(small.values)
  chklogical(keep.sample)
  ##
  if (pooled == "common") {
    TE.pooled <- x$TE.common
    Cov.pooled <- x$Cov.common
  }
  else {
    TE.pooled <- x$TE.random
    Cov.pooled <- x$Cov.random
  }
  ##
  if (small.values == "desirable")
    theta <- TE.pooled[, 1]
  else
    theta <- TE.pooled[1, ]
  ##
  compMatrix <- matrix(0, nrow = nrow(Cov.pooled), ncol = length(x$trts))
  rownames(compMatrix) <- rownames(Cov.pooled)
  colnames(compMatrix) <- x$trts
  ##
  allcomps <- compsplit(rownames(compMatrix), x$sep.trts)
  for (i in seq_len(nrow(Cov.pooled)))
    compMatrix[i, allcomps[[i]]] <- 1
  ##
  var.theta <- as.vector(ginv(compMatrix) %*% diag(Cov.pooled))
  ##
  sample <- mvtnorm::rmvnorm(nsim, theta, diag(var.theta))
  rownames(sample) <- seq_len(nrow(sample))
  colnames(sample) <- x$trts
  ##
  ## Ranks
  ##
  rnk <- apply(sample, 1, rank, ties.method = "random")
  ##
  ## Rankogram
  ##
  tab <- rank.tab.list(rnk, ntrts = x$n)
  names(tab) = x$trts
  ##
  if (is.list(tab)) {
    rankogram <- matrix(0, nrow = x$n, ncol = x$n)
    rownames(rankogram) <- names(tab)
    colnames(rankogram) <- seq_len(x$n)
    ##
    for (i in names(tab))
      rankogram[i, names(tab[[i]])] <- tab[[i]][names(tab[[i]])]
  }   else
    rankogram <- t(as.data.frame(tab))
  ##
  ## Cumulative ranks
  ##
  cumrank <- t(apply(rankogram, 1, cumsum))
  ##
  ## SUCRAs
  ##
  if( dim(TE.pooled)[2] == 2 ){
    ranking <- cumrank[, -x$n]
  } else {
    ranking <- apply(cumrank[, -x$n], 1, sum) / (x$n - 1)
  }
  ##
  ## Return results
  ##
  res <- list(ranking = ranking / nsim,
              rankogram = rankogram / nsim,
              cumrank = cumrank / nsim,
              ##
              nsim = nsim,
              pooled = pooled,
              small.values = small.values,
              keep.sample = keep.sample,
              ##
              compMatrix = compMatrix)
  ##
  if (keep.sample)
    res[["sample"]] <- sample
  ##
  res
}

# Updated rankogram() by using function ranksampling_plus().
rankogram_plus <- function(x, nsim = 1000,
                           common = x$common, random = x$random,
                           small.values = x$small.values,
                           cumulative.rankprob = FALSE,
                           nchar.trts = x$nchar.trts,
                           warn.deprecated = gs("warn.deprecated"),
                           ...) {
  
  ##
  ##
  ## (1) Check for netmeta object and upgrade object
  ##
  ##
  chkclass(x, "netmeta")
  x <- updateversion(x)
  ##
  is.installed.package("mvtnorm")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(nsim, min = 1, length = 1)
  ##
  small.values <- setsv(small.values)
  ##
  chklogical(cumulative.rankprob)
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  chknumeric(nchar.trts, length = 1)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Resampling to calculate ranking probabilites and SUCRAs
  ##
  ##
  ranking.common  <- ranking.matrix.common  <- cumrank.matrix.common  <- NULL
  ranking.random <- ranking.matrix.random <- rank.cum.random <- NULL
  ##
  if (common) {
    res.f <- ranksampling_plus(x, nsim, "common", small.values)
    ##
    ranking.common <- res.f$ranking
    ranking.matrix.common <- res.f$rankogram
    cumrank.matrix.common <- res.f$cumrank
  }
  ##
  if (random) {
    res.r <- ranksampling_plus(x, nsim, "random", small.values)
    ##
    ranking.random <- res.r$ranking
    ranking.matrix.random <- res.r$rankogram
    rank.cum.random <- res.r$cumrank
  }
  
  
  ##
  ##
  ## (4) Create rankogram object
  ##
  ##
  res <- list(ranking.common = ranking.common,
              ranking.matrix.common = ranking.matrix.common,
              cumrank.matrix.common = cumrank.matrix.common,
              ##
              ranking.random = ranking.random,
              ranking.matrix.random = ranking.matrix.random,
              cumrank.matrix.random = rank.cum.random,
              ##
              nsim = nsim,
              ##
              common = common,
              random = random,
              small.values = small.values,
              cumulative.rankprob = cumulative.rankprob,
              ##
              nchar.trts = nchar.trts,
              x = x,
              version = packageDescription("netmeta")$Version
  )
  ##
  ## Backward compatibility
  ##
  res$fixed <- common
  ##
  res$ranking.fixed <- ranking.common
  res$ranking.matrix.fixed <- ranking.matrix.common
  res$cumrank.matrix.fixed <- cumrank.matrix.common
  ##
  class(res) <- "rankogram"
  
  res
}





#' @rdname rankogram
#' @method print rankogram
#' @export


print.rankogram <- function(x,
                            common = x$common,
                            random = x$random,
                            cumulative.rankprob = x$cumulative.rankprob,
                            nchar.trts = x$nchar.trts,
                            digits = gs("digits.prop"),
                            legend = TRUE,
                            warn.deprecated = gs("warn.deprecated"),
                            ...) {
  
  ##
  ##
  ## (1) Check for rankogram object and upgrade object
  ##
  ##
  chkclass(x, "rankogram")
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(cumulative.rankprob)
  ##
  chknumeric(nchar.trts, length = 1)
  ##
  chknumeric(digits, length = 1)
  chklogical(legend)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random", warn.deprecated)
  chklogical(random)
  
  
  ##
  ##
  ## (3) Print results
  ##
  ##
  if (common | random)
    cat(if (cumulative.rankprob)
      "Cumulative ranking probabilities" else "Rankogram",
      " (based on ", x$nsim, " simulation",
      if (x$nsim > 1) "s", ")\n\n",
      sep = "")
  ##
  if (common) {
    if (cumulative.rankprob)
      rank.common <- x$cumrank.matrix.common
    else
      rank.common <- x$ranking.matrix.common
    rownames(rank.common) <- treats(rank.common, nchar.trts)
    ##
    cat("Common effects model: \n\n")
    prmatrix(formatN(rank.common, digits), quote = FALSE, right = TRUE, ...)
    if (random)
      cat("\n")
  }
  ##
  if (random) {
    if (cumulative.rankprob)
      rank.random <- x$cumrank.matrix.random
    else
      rank.random <- x$ranking.matrix.random
    rownames(rank.random) <-
      treats(rank.random, nchar.trts)
    ##
    cat("Random effects model: \n\n")
    prmatrix(formatN(rank.random, digits), quote = FALSE, right = TRUE, ...)
  }
  ##
  ## Add legend with abbreviated treatment labels
  ##
  if ((common | random) & legend) {
    if (common)
      trts <- rownames(x$ranking.matrix.common)
    else if (random)
      trts <- rownames(x$ranking.matrix.random)
    ##
    legendabbr(trts, treats(trts, nchar.trts), TRUE)
  }
  
  
  invisible()
}


