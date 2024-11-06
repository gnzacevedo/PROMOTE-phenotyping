calculateNRS <- function(x, features){
  require(CATALYST)
  features <- CATALYST:::.get_features(x, features)
  y <- assay(x, 'exprs')
  y <- y[features, ]
  cs_by_s <- split(seq_len(ncol(x)), x$sample_id)
  nrs <- lapply(cs_by_s, function(cs) CATALYST:::.nrs(y[, cs, drop = FALSE]))
  
  #df <- as.data.frame(do.call(rbind, vectorList))
  nrs.df <- do.call(rbind, nrs) %>% as.data.frame()
  return(nrs.df)
}
