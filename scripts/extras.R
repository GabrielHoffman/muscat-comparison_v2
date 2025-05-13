
getExprGeneNames <- function(sceObj, assays = assayNames(sceObj), min.cells = 5, min.count = 5, min.samples = 4, min.prop = .4, min.total.count = 15, normalize.method = "TMM"){

  # checks
  stopifnot(is(sceObj, "SingleCellExperiment"))

  # check colnames of SCE
  if (is.null(colnames(sceObj))) {
    stop("colnames(sceObj) is NULL.  Column names are needed for internal filtering")
  }

  # extract metadata shared across assays
  data_constant <- droplevels(as.data.frame(colData(sceObj)))

  # check if assays are valid
  if (any(!assays %in% assayNames(sceObj))) {
    idx <- which(!assays %in% assayNames(sceObj))
    txt <- paste("Assays are not found in dataset:", paste(head(assays[idx]), collapse = ", "))
    stop(txt)
  }

  # extract cell counts
  n.cells_full <- cellCounts(sceObj)

  # extract all unique colnames
  colNamesAll <- unique(unlist(lapply(assayNames(sceObj), function(x) colnames(assay(sceObj, x)))))

  # check for colnames missing cell counts
  if (any(!colNamesAll %in% rownames(n.cells_full))) {
    stop("Cell counts could not be extracted.\n  Check that colnames(sceObj) or rownames(colData(sceObj))\n  have not been modified manually after running aggregateToPseudoBulk()")
  }

   # for each assay
  resList <- lapply(assays, function(k) {
    startTime <- Sys.time()

    y <- as.matrix(assay(sceObj, k))

    # cell counts
    n.cells <- n.cells_full[colnames(y), k, drop = FALSE]

    # merge data_constant (data constant for all cell types)
    # with metadata(sceObj)$aggr_means (data that varies)
    data <- dreamlet:::merge_metadata(
      data_constant,
      dreamlet:::get_metadata_aggr_means(sceObj),
      k,
      metadata(sceObj)$agg_pars$by
    )

    y = y[, rownames(data), drop = FALSE]

    # samples to include of they have enough observed cells
    include <- (n.cells >= min.cells)

    # if no samples are retained
    if (sum(include) == 0) {
      return(NULL)
    }

    # subset expression and data
    y <- y[, include, drop = FALSE]
    data <- droplevels(data[include, , drop = FALSE])

    # if there are too few remaining samples
    if (nrow(data) < min.samples | nrow(y) == 0) {
      return(NULL)
    }

    # Get count data and normalize
    y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
    y <- calcNormFactors(y, method = normalize.method)

    # get samples with enough cells
    # filter genes
    # design: model.matrix( subbars(formula), data)
    # Design often includes batch and donor, which are very small
    #   this causes too many genes to be retained
    keep <- suppressWarnings(
    filterByExpr(y, 
      min.count = min.count, 
      min.prop = min.prop, 
      min.total.count = min.total.count)
    )
    names(keep)[keep]
  })
  names(resList) <- assays
  resList
}



pbWeights <- function(sce, sample_id, cluster_id, geneList = NULL, method = c("delta", "ncells"), shrink = TRUE, prior.count = 0.5, maxRatio = 20, h5adBlockSizes = 1e9, details = FALSE, verbose = TRUE) {
  method <- match.arg(method)

  # check for NA values
  if (any(is.na(sce[[sample_id]]))) {
    stop("NA values are not allowed in sample_id column")
  }
  if (any(is.na(sce[[cluster_id]]))) {
    stop("NA values are not allowed in cluster_id column")
  }

  if (method == "ncells") {
    W.list <- .pbWeights_ncells(sce, sample_id, cluster_id)
  } else {
    # delta approximation

    # update block size for reading h5ad file from disk
    tmp <- getAutoBlockSize()
    suppressMessages(setAutoBlockSize(h5adBlockSizes))
    on.exit(suppressMessages(setAutoBlockSize(tmp)))

    # compute variances
    var.lst <- getVarList(sce, sample_id, cluster_id, geneList, shrink, prior.count, details = details, verbose = verbose)

    # for each cell type
    W.list <- lapply( names(var.lst), function(id) {
      v.mat = var.lst[[id]]
      # regularize reciprocal with offset
      # get offset tau as the max of x
      #  to give a maximum ratio of maxRatio
      # for each cell type
      ret <- t(apply(v.mat, 1, function(x) {
        tau <- get_offset(x, maxRatio)
        1 / (x + tau)
      }))

      if( details ) attr(ret, "details") <- attr(var.lst[[id]], "details") 
      ret
    })
    names(W.list) <- names(var.lst)
  }

  W.list
}

