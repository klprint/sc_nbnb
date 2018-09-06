####
# 09.08.2018
# Kevin Leiss
# Updated NB fitting and convenience functions
####


filter_genes = function(counts, min_cells = 20){
  counts = counts[rowSums(counts) > min_cells, ]
  return(counts)
}


anscombe_transform = function(counts){
  apply(counts, 2, function(umi) {
    ans <- sqrt(umi+3/8) - sqrt(3/8)
    ans/sqrt(sum(ans^2)) })
}


calculate_size_factor = function(counts){
  Matrix::colSums(counts) / mean(Matrix::colSums(counts))
}

disp <- function(k) (var(k) - mean(k)) / mean(k) / mean(k) # var = mu + disp * mu^2


fitNB <- function(k, sf, initialVals = c(round(mean(k), 10), round(disp(k), 10))) {
  k = as.numeric(k)
  o = optim(
    initialVals,  # mu, alpha
    function(x) {
      -sum( lgamma( k + 1/x[2] ) - lgamma( 1/x[2] ) - lgamma( k+1 ) - ( k + 1/x[2] ) * log( 1 + sf * x[2] * x[1] ) + k * log( sf * x[2] * x[1] ) )
    },

    function(x) c(
      -sum( ( k - sf * x[1] ) / ( x[1] + sf * x[2] * x[1]^2 ) ),
      -sum( ( x[2] * ( k - sf * x[1] ) / ( 1 + sf * x[2] * x[1] ) + log( 1 + sf * x[2] * x[1] ) - digamma( k + 1/x[2] ) + digamma( 1/x[2] ) ) / x[2]^2 ) ),
    hessian = TRUE,
    lower = c( 1e-10, 1e-10 ),
    method = "L-BFGS-B" )
  c( mean = o$par[1],
     disp = o$par[2],
     SE_m = 1 / sqrt(o$hessian[1,1]),
     SE_d = 1 / o$hessian[2,2] ) }


library(purrr) # for possibly
safe_fitNB <- possibly(fitNB, otherwise = c(mean=NA, disp=NA, SE_m = NA, SE_d = NA))


stack_disp_table = function(dispList){
  dtable = NULL
  for(i in 1:length(dispList)){
    tmp = dispList[[i]]
    dtable = rbind(dtable, tmp)
  }

  return(dtable)
}


calculate_disp_table = function(counts, sf = calculate_size_factor(counts)){
  require(pbmcapply)


  disp.table = pbmclapply(rownames(counts), function(g){
    NBparams = safe_fitNB(counts[g, ], sf)
    data.frame(gene = g, t(NBparams))
  })

  return(stack_disp_table(disp.table))
}


inspect_disptable = function(disptable){
  signif_dispersion = (disptable$disp / disptable$SE_d) > 2

  noNA = colSums(apply(disptable, 1, is.na)) == 0

  disptable$noNA = noNA
  disptable$signif_dispersion = signif_dispersion


  disp_plot = ggplot(disptable, aes(x = mean, y = disp)) +
    geom_point(size = 0.25, alpha = 0.2, aes(color = noNA & signif_dispersion))+
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(base = 10, size = 0.2)

  return(list("signif_dispersion" = signif_dispersion,
              "noNA" = noNA,
              "plot" = disp_plot))

}

trainNB <- function(countmatrix, isPositive, isNegative, sf) {
  ## FYI: faster for dense matrices
  pbmclapply(rownames(countmatrix), function(g) {
    if(sum(countmatrix[g, ]) == 0) return(
      data.frame(gene = NA, meanPos=NA, dispPos=NA, meanNeg=NA, dispNeg=NA))

    posCounts <- countmatrix[g, isPositive]
    negCounts <- countmatrix[g, isNegative]
    # if one group has 100% zeros this overemphasizes the genes, so remedy:
    if(sum(posCounts) == 0) posCounts <- c(1, rep(0, sum(isPositive)-1))
    if(sum(negCounts) == 0) negCounts <- c(1, rep(0, sum(isNegative)-1))

    pos =  safe_fitNB(posCounts, sf[isPositive])
    neg =  safe_fitNB(negCounts, sf[isNegative])
    data.frame(  gene = g,
                 meanPos = pos["mean"], dispPos = pos["disp"],
                 meanNeg = neg["mean"], dispNeg = neg["disp"],
                 stringsAsFactors = F, row.names = NULL
    )
  })
}


NBNB <- function(countmatrix, dispersionTable, sf = NULL) {
  if(is.null(sf)) {
    sf <- colSums(countmatrix)
    names(sf) <- colnames(countmatrix)
  }
  sapply(colnames(countmatrix), function(cell) {

    sum(dnbinom(x  = countmatrix[dispersionTable$gene, cell],
                mu = dispersionTable$meanPos * sf[cell],
                size = 1 / dispersionTable$dispPos, log=T), na.rm = T) +
      -
      sum(dnbinom(x  = countmatrix[dispersionTable$gene, cell],
                  mu = dispersionTable$meanNeg * sf[cell],
                  size = 1 / dispersionTable$dispNeg, log=T), na.rm = T)
  })


}


