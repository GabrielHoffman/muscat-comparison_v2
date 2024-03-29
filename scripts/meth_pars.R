# this determines which types of methods to include
# names(ids) <- ids <- c("pb", "ad", "scdd", "mast", "mm")
names(ids) <- ids <- c("pb", "mast")

# aggregation-based ------------------------------------------------------------
pb <- dplyr::bind_rows(
    expand.grid(
        stringsAsFactors = FALSE,
        assay = "counts", fun = "sum", scale = FALSE, 
        method = c("edgeR", "limma-voom"),
        treat = c(TRUE, FALSE)
    ),
    expand.grid(
        stringsAsFactors = FALSE,
        assay = "counts", fun = "sum", scale = FALSE, 
        # , "dreamlet_delta", "dreamlet_ncells"
        method = c("dreamlet_deltaW", "dreamlet_none", "DESeq2"),
        treat = c(FALSE)),
    # expand.grid(
    #     stringsAsFactors = FALSE,
    #     assay = "counts", fun = "sum", scale = FALSE, 
    #     method = "dreamlet_no_cell_weights",
    #     treat = c(FALSE)),
    expand.grid(
        stringsAsFactors = FALSE, scale = FALSE,
        assay = c("logcounts", "vstresiduals"),
        fun = "mean", method = "limma-trend"
    ),
    data.frame(
        stringsAsFactors = FALSE, scale = TRUE,
        assay = "cpm", fun = "sum", method = "edgeR",
        treat = FALSE
    )    
)
pb$treat[is.na(pb$treat)] <- FALSE
pb$id <- with(pb, sprintf("%s%s.%s.%s%s", 
    method, ifelse(treat, "-treat", ""), 
    fun, ifelse(scale, "scale", ""), assay))

# mixed-models -----------------------------------------------------------------
mm <- data.frame(
    stringsAsFactors = FALSE,
    method = c("dream", "dream2", "nbinom", "vst"),
    vst = c("", "", "", "sctransform"),
    ddf = "Satterthwaite")
mm$id <- with(mm, paste0("MM-", method))

# Anderson-Darling -------------------------------------------------------------
ad <-  expand.grid(
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstresiduals"),
    var = c("sample_id", "group_id"))
ad$id <- with(ad, sprintf("AD-%s.%s", gsub("(.).*", "\\1id", var), assay))

# scDD -------------------------------------------------------------------------
scdd <- data.frame(
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstresiduals"))
scdd$id <- with(scdd, sprintf("scDD.%s", assay))

# MAST -------------------------------------------------------------------------
mast <- data.frame(
    stringsAsFactors = FALSE,
    assay = "logcounts",
    id = "MAST.logcounts")

# write method IDs to .csv -----------------------------------------------------
for (id in ids) {
    ps <- get(id)
    assign(id, split(ps, ps$id))
}

pars <- sapply(ids, get)
typs <- rep.int(ids, sapply(pars, length))
pars <- purrr::flatten(pars)
mids <- names(pars)

mids_df <- data.frame(
    row.names = NULL,
    stringsAsFactors = FALSE,
    id = mids, type = typs)

write.csv(mids_df, config$mids)

# write method parameters to .json ---------------------------------------------
fns <- paste0(mids, ".json")
fns <- paste0(config$meth_pars, fns)
names(fns) <- mids

for (id in mids) {
    new <- as.list(pars[[id]])
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!identical(old, new))
            write(jsonlite::toJSON(new), fns[id])
    } else {
        write(jsonlite::toJSON(new), fns[id])
    }
}
