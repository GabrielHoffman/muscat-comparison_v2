.typ_cols <- c(pb = "grey50", mm = "violet", ad = "lightgreen", mast = "orange", scdd = "steelblue")
.typ_labs <- c(pb = "PB", mm = "MM", ad = "AD", mast = "MAST", scdd = "scDD")

.cat_cols <- c("royalblue", "cornflowerblue", "red3", "tomato", "orange", "gold")
names(.cat_cols) <- c("ee", "ep", "de", "dp", "dm", "db")

# colors = c(dreamlet = "#E41A1C", 
#       'limma-voom' = "#377EB8", 
#       edgeR = "#4DAF4A", 
#       DESeq2 = "#984EA3", 
#       glmer = "#FF7F00", 
#       MAST = "#f2e718")


.meth_cols <- rev(c(    
    "dreamlet_deltaW.sum.counts" = "#E41A1C",
    # "dreamlet_delta.sum.counts" = "#E41A1C",
    # "dreamlet_ncells.sum.counts" = "#f76b5c",
    "dreamlet_none.sum.counts" = "#d1696a",

    "limma-voom.sum.counts" = "#377EB8",
    "limma-trend.mean.logcounts"    = "#95bcdb",
    "limma-trend.mean.vstresiduals" = "#cfeaff",

    "DESeq2.sum.counts" = "#984EA3",    

    "edgeR.sum.counts" = "#4DAF4A",
    "edgeR.sum.scalecpm" = "#9cf099",

    #"scDD.logcounts"    = "#0056B2",
    #"scDD.vstresiduals" = "#009EF6",     
    
    "MAST.logcounts"    = "#ed68b8"))
    
#     "MM-dream"  = "#005E5C",
#     "MM-dream2" = "#00ABA6",
#     "MM-nbinom" = "#00E3DD",
#     "MM-vst"    = "#95EEE8",
    
#     "AD-gid.logcounts"    = "#9A41B3",
#     "AD-gid.vstresiduals" = "#FFA9FF",
#     "AD-sid.logcounts"    = "#E56D4B",
#     "AD-sid.vstresiduals" = "#FBB6A2")

#cols <- .meth_cols
#hist(seq_along(cols), breaks = c(seq_along(cols) - 0.5, length(cols) + 0.5), col = cols)

i <- grep("edgeR|limma-v", names(.meth_cols))
.treat_mids <- names(.meth_cols)[i]
.treat_mids <- c(.treat_mids, gsub("(.*?)(\\.)(.*)", "\\1-treat.\\3", .treat_mids))[
    c(sapply(seq_along(.treat_mids), function(i) c(i, i + length(.treat_mids))))]
.treat_cols <- rep(.meth_cols[i], each = 2)
names(.treat_cols) <- .treat_mids

.read_res <- function(fns, include = "all", slot = "tbl") {
    if (include == "treat") {
        fns <- fns[matrixStats::rowAnys(sapply(.treat_mids, grepl, fns))]
        mids <- .treat_mids
    } else {
        fns <- fns[!grepl("treat", fns)]
        mids <- names(.meth_cols)
    }
    res <- map(lapply(fns, readRDS), slot)
    rmv <- vapply(res, function(u) 
        is.null(u) | inherits(u, "error"), 
        logical(1))
    res <- res[!rmv]
    if (slot == "tbl")
        res <- map(res, mutate_if, is.factor, as.character) %>% 
            bind_rows %>% mutate_if(is.character, as.factor) %>% 
            mutate_at("category", factor, levels = muscat:::cats) %>% 
            mutate_at("mid", factor, levels = mids) %>% 
            mutate_if(is.factor, droplevels)
    return(res)
}

.filter_matrix <- function(m, n = 100) {
    while (any(m < n)) {
        # get candidate rows/cols for removal
        i <- m < n
        r <- apply(i, 1, any)
        c <- apply(i, 2, any)
        # get smallest row/col
        rs <- rowSums(m)
        cs <- colSums(m)
        r <- which(r)[which.min(rs[r])]
        c <- which(c)[which.min(cs[c])]
        # priorities removal of rows over cols
        if (rs[r] <= cs[c]) {
            m <- m[-r, , drop = FALSE]
        } else {
            m <- m[, -c, drop = FALSE]
        }
        if (any(dim(m) == 1)) 
            break
    }
    return(m)
}

.update_sce <- function(sce) {
    # update colData
    cd <- as.data.frame(colData(sce))
    cd <- mutate_if(cd, is.factor, droplevels) 
    colData(sce) <- DataFrame(cd, row.names = colnames(sce))
    # update metadata
    ei <- metadata(sce)$experiment_info
    ei <- ei[ei$sample_id %in% levels(sce$sample_id), ]
    ei <- mutate_if(ei, is.factor, droplevels)
    metadata(sce)$experiment_info <- ei
    return(sce)
}

.filter_sce <- function(sce, kids, sids) {
    cs1 <- sce$cluster_id %in% kids
    cs2 <- sce$sample_id %in% sids
    sce <- sce[, cs1 & cs2]
    sce <- .update_sce(sce)
    return(sce)
}

.update_cd <- function(sce) {
    colData(sce) <- as.data.frame(colData(sce)) %>% 
        mutate_if(is.factor, droplevels) %>% 
        DataFrame(row.names = colnames(sce))
    return(sce)
}

.prettify <- function(theme = NULL, ...) {
    if (is.null(theme)) theme <- "classic"
    base <- paste0("theme_", theme)
    base <- getFromNamespace(base, "ggplot2")
    base(base_size = 8) + theme(
        aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(2, "mm"),
        strip.background = element_rect(fill = NA),
        plot.margin = unit(rep(1, 4), "mm"),
        panel.spacing = unit(0, "mm"),
        legend.margin = margin(0,0,1,0,"mm"),
        ...)}

.plot_perf_points <- function(df, 
    include = "all", color_by = "method", facet = "splitval") {
    df <- filter(ungroup(df), TPR + FDR != 0)
    df$treat <- grepl("treat", df$method)
    if (any(rmv <- table(df$splitval) < 2))
        df <- filter(df, !splitval %in% levels(df$splitval)[rmv])
    df$method <- factor(
        gsub("-treat", "", df$method), 
        levels = levels(df$method))
    p <- ggplot(filter(df, FDR + TPR != 0),
        aes_string(x = "FDR", y = "TPR", col = color_by)) +
        facet_wrap(facet, labeller = labeller(.multi_line = FALSE)) +
        geom_vline(size = 0.2, lty = 2, aes(xintercept = thr)) + 
        geom_point(size = 1, alpha = 0.8) + 
        geom_line(aes(lty = treat), size = 0.4, alpha = 0.4, show.legend = (include == "treat")) +
        scale_color_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
        scale_x_sqrt(expand = c(0, 0.05),
            breaks = c(c(0.01, switch(include == "treat", 0.05, NULL), 0.1), seq(0.2, 1, 0.2)), 
            limits = c(0, ifelse(include == "treat", 0.05, 1)),
            labels = function(x) format(x, drop0trailing = TRUE)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
        .prettify(theme = "bw", legend.position = "bottom", legend.box.just = "center") + 
        guides(
            lty = guide_legend(ncol = 1, keywidth = unit(4, "mm"), 
                override.aes = list(alpha = 1, lty = c(1, 3))),
            col = guide_legend(order = 1, nrow = 4, 
                override.aes = list(alpha = 1, size = 2)))
    suppressWarnings(suppressMessages(p))
}
