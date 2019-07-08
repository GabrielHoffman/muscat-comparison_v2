source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

x <- snakemake@wildcards$x
#fns <- list.files("/users/helena/dropbox/portmac/results/kang", "ds10_ns;", full.names = TRUE)
res <- lapply(snakemake@input$res, readRDS) %>% map("tbl")
rmv <- vapply(res, inherits, what = "error", logical(1))
res <- map(res[!rmv], mutate_if, is.factor, as.character) %>% 
    bind_rows %>% setDT %>% split(by = "j", flatten = FALSE) %>% 
    map(group_by, mid) %>% map(function(u) 
        set_names(group_split(u), group_keys(u)[[1]]))

cd <- lapply(seq_along(res), function(i) {
    map(res[[i]][[1]], select, c(x, "is_de")) %>% 
        bind_rows
    truth <- lapply(c(x, "is_de"), map, .x = res[[i]]) %>% 
        map(unlist) %>% set_names(c(x, "is_de")) %>% 
        data.frame(row.names = NULL, check.names = FALSE)
    pvals <- lapply(c("p_val", "p_adj.loc"), map, .x = res[[i]]) %>% 
        map(data.frame, check.names = FALSE)
    dfs <- c(list(truth), pvals)
    names(dfs) <- c("truth", "pval", "padj")
    do.call(COBRAData, dfs)
})

perf <- lapply(cd, calculate_performance,
    binary_truth = "is_de", 
    aspects = c("fdrtpr", "fdrtprcurve"),
    splv = x, maxsplit = Inf)

df <- map(perf, "fdrtpr") %>% 
    bind_rows(.id = "j") %>% 
    select(splitval, thr, method, TPR, FDR) %>% 
    dplyr::filter(splitval != "overall") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%  
    mutate_at("splitval", function(u) {
        u <- gsub(paste0(x, ":"), "", u)
        v <- sort(unique(as.numeric(u)))
        factor(u, levels = v)
    }) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)