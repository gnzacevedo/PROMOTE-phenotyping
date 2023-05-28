#Modded function from package Spectre, required to overcome error when creating a data table from a SCE without column names for each cell
create.dt.mod <- function (dat, from = NULL) 
{
  require("data.table")
  if (class(dat)[1] == "Seurat") {
    object.type <- "Seurat"
  }
  if (class(dat)[1] == "SingleCellExperiment") {
    object.type <- "SingleCellExperiment"
  }
  if (class(dat)[1] == "flowFrame") {
    object.type <- "flowFrame"
  }
  if (!is.null(from)) {
    object.type <- from
  }
  if (!exists("object.type")) {
    stop("Could not determine the object type. Currently only Seurat objects are supported. You can also try manually specifying the object type using the 'from' argument (e.g. from = 'Seurat'")
  }
  if (object.type == "Seurat") {
    message(object.type, " detected")
    require("dplyr")
    require("Seurat")
    require("patchwork")
    a <- GetAssayData(object = dat)
    assays <- names(dat@assays)
    dim.reds <- names(dat@reductions)
    var.features <- VariableFeatures(dat)
    if (length(var.features) > 0) {
      var.features.top10 <- head(VariableFeatures(dat), 
                                 10)
    }
    geneNames <- a@Dimnames[[1]]
    cellNames <- a@Dimnames[[2]]
    res <- as.data.table(cellNames)
    names(res) <- "cellNames"
    if (!is.null(dat@meta.data)) {
      col.meta <- dat@meta.data
      col.meta <- as.data.table(col.meta)
      meta.cols <- names(col.meta)
      res <- cbind(res, col.meta)
    }
    for (i in assays) {
      types <- vector()
      if (ncol(dat@assays[[i]]@counts) > 0) {
        types <- c(types, "counts")
        x1 <- GetAssayData(object = dat, assay = i, 
                           slot = "counts")
        x1 <- as.data.table(x1)
        x1 <- data.table::transpose(x1)
        names(x1) <- geneNames
        names(x1) <- paste0(names(x1), "_", i, "_", 
                            "counts")
        res <- cbind(res, x1)
      }
      if (ncol(dat@assays[[i]]@data) > 0) {
        types <- c(types, "data")
        x2 <- GetAssayData(object = dat, assay = i, 
                           slot = "data")
        x2 <- as.data.table(x2)
        x2 <- data.table::transpose(x2)
        names(x2) <- geneNames
        names(x2) <- paste0(names(x2), "_", i, "_", 
                            "data")
        res <- cbind(res, x2)
      }
      if (ncol(dat@assays[[i]]@scale.data) > 0) {
        types <- c(types, "scale.data")
        x3 <- GetAssayData(object = dat, assay = i, 
                           slot = "scale.data")
        x3 <- as.data.table(x3)
        x3 <- data.table::transpose(x3)
        names(x3) <- geneNames
        names(x3) <- paste0(names(x3), "_", i, "_", 
                            "scale.data")
        res <- cbind(res, x3)
      }
      rm(i)
      rm(x1)
      rm(x2)
      rm(x3)
    }
    for (i in dim.reds) {
      tmp <- dat@reductions[[i]]@cell.embeddings
      tmp <- as.data.table(tmp)
      names(tmp) <- paste0(i, "_", names(tmp))
      res <- cbind(res, tmp)
    }
    final.res <- list()
    final.res$data.table <- res
    final.res$geneNames <- geneNames
    final.res$cellNames <- cellNames
    final.res$meta.data <- meta.cols
    if (length(var.features) > 0) {
      final.res$var.features <- var.features
      final.res$var.features.top10 <- var.features.top10
    }
    final.res$assays <- paste0("_", assays)
    final.res$slots <- paste0("_", types)
    final.res$dim.reds <- paste0(dim.reds, "_")
    message(paste0("Converted a ", object.type, " object into a data.table stored in a list"))
    return(final.res)
  }
  if (object.type == "SingleCellExperiment") {
    message(object.type, " detected")
    require("SingleCellExperiment")
    geneNames <- rownames(dat)
    geneNames
    cellNames <- colnames(dat)
    # GA edits #
    if(is.null(cellNames)){
      cellNames <- paste0('cell', 1:ncol(dat))
    }
    ##
    cellNames
    assays <- names(dat@assays@data)
    assays
    dim.reds <- names(reducedDims(dat))
    dim.reds
    res <- as.data.table(cellNames)
    names(res) <- "cellNames"
    message("-- Adding metadata")
    if (!is.null(dat@colData)) {
      col.meta <- dat@colData
      col.meta <- as.data.table(col.meta)
      meta.cols <- names(col.meta)
      res <- cbind(res, col.meta)
      res
    }
    message("-- Adding assay data")
    for (i in assays) {
      tmp <- dat@assays@data[[i]]
      tmp <- as.matrix(tmp)
      tmp <- as.data.table(tmp)
      tmp <- data.table::transpose(tmp)
      names(tmp) <- geneNames
      names(tmp) <- paste0(names(tmp), "_", i)
      res <- cbind(res, tmp)
      rm(i)
      rm(tmp)
    }
    message("-- Adding DimRed data")
    for (i in dim.reds) {
      tmp <- reducedDims(dat)[[i]]
      tmp <- as.data.table(tmp)
      new.names <- paste0(i, "_", c(1:length(names(tmp))))
      names(tmp) <- new.names
      res <- cbind(res, tmp)
      rm(i)
      rm(tmp)
    }
    message("-- Finalising")
    final.res <- list()
    final.res$data.table <- res
    final.res$geneNames <- geneNames
    final.res$cellNames <- cellNames
    final.res$meta.data <- meta.cols
    final.res$assays <- paste0("_", assays)
    final.res$dim.reds <- paste0(dim.reds, "_")
    message(paste0("Converted a ", object.type, " object into a data.table stored in a list"))
    return(final.res)
  }
  if (object.type == "flowFrame") {
    message(object.type, " detected")
    require("flowCore")
    require("data.table")
    res <- exprs(dat)
    res <- res[1:nrow(res), 1:ncol(res)]
    res <- as.data.table(res)
    for (i in names(res)) {
      if (!any(is.na(as.numeric(as.character(res[[i]]))))) {
        res[[i]] <- as.numeric(res[[i]])
      }
    }
    final.res <- list()
    final.res$data.table <- res
    final.res$parameters <- dat@parameters
    final.res$description <- dat@description
    message(paste0("Converted a ", object.type, " object into a data.table stored in a list"))
    return(final.res)
  }
}
