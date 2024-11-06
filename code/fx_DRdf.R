#extracts reduced dimensionality coordinates from SingleCellExperiment, extracting selected md columns alongside
DRdf <- function(x, dr, md=NULL){
  kids <- names(cluster_codes(x))
  u <- c(rownames(x), colnames(colData(x)), kids)
  stopifnot(length(SingleCellExperiment:::reducedDimNames(x)) != 0)
  if (is.null(dr)) {
    dr <- reducedDimNames(x)[1]
  }
  stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
  xy <- SingleCellExperiment:::reducedDim(x, dr)
  df <- data.frame(colData(x), xy)
  
  if(!is.null(md)){
    for(mdcol in md){
      if (mdcol %in% rownames(x)) {
        df[[mdcol]] <- assay(x, "exprs")[mdcol, ]
      }
      
      if (mdcol %in% kids) {
        df[[mdcol]] <- cluster_ids(x, mdcol)
      }
    }
  }
  return(df)
}
