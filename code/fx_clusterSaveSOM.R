# Main FlowSOM function####
clusterSaveSOM <- function (x, features = "type", xdim = 10, ydim = 10, maxK = 20,
                            som.rlen=20,
                            verbose = TRUE, seed = 1, marker.importance=NULL) {
  stopifnot(is(x, "SingleCellExperiment"))
  stopifnot(is.logical(verbose), length(verbose) == 1, vapply(list(xdim, 
                                                                   ydim, maxK, seed), function(arg) is.numeric(arg) && 
                                                                length(arg) == 1, logical(1)))
  features <- CATALYST:::.get_features(x, features)
  if (is.null(marker_classes(x))) {
    rowData(x)$marker_class <- factor(c("state", "type")[as.numeric(rownames(x) %in% 
                                                                      features) + 1], levels = c("type", "state", "none"))
  }
  rowData(x)$used_for_clustering <- rownames(x) %in% features
  if (verbose) 
    message("o running FlowSOM clustering...")
  fsom <- FlowSOM::ReadInput(flowFrame(t(assay(x, "exprs"))))
  som <- FlowSOM::BuildSOM(fsom, colsToUse = features, silent = TRUE, 
                           xdim = xdim, ydim = ydim, importance=marker.importance, rlen=som.rlen)
  som <- FlowSOM::BuildMST(som)
  if (verbose) 
    message("o running ConsensusClusterPlus metaclustering...")
  pdf(NULL)
  mc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(t(som$map$codes), 
                                                                    maxK = maxK, reps = 100, distance = "euclidean", seed = seed, 
                                                                    plot = NULL))
  dev.off()
  k <- xdim * ydim
  mcs <- seq_len(maxK)[-1]
  codes <- data.frame(seq_len(k), purrr::map(mc[-1], "consensusClass"))
  codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
  colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", 
                                                    mcs))
  x$cluster_id <- factor(som$map$mapping[, 1])
  metadata(x)$cluster_codes <- codes
  metadata(x)$SOM_codes <- som$map$codes
  metadata(x)$delta_area <- CATALYST:::.plot_delta_area(mc)
  metadata(x)$flowSOM_object <- som
  return(x)
}