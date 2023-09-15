suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

# wcs = list(x = "s", did = "kang")
# args <- list(
#     res = list.files("results", sprintf("%s,de10_n%s,", wcs$did, wcs$x), full.names = TRUE),
#     ggp = file.path("plots", paste0(wcs$did, sprintf("-perf_by_n%s.rds", wcs$x))),
#     fig = file.path("plots", paste0(wcs$did, sprintf("-perf_by_n%s.pdf", wcs$x))))

res <- .read_res(args$res) %>% 
    dplyr::mutate(id = sprintf(
        "i%sj%sc%ss%s%s%s", 
        i, j, c, s, gene, cluster_id)) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = "i", flatten = FALSE) %>% 
    map(group_by, mid) %>% map(function(u) 
        setNames(group_split(u), group_keys(u)[[1]]))

# some methods may fail for too low number of cells / replicates;
# the below chunk fills in missing results to match in dimension
for (i in seq_along(res)) {
    u <- res[[i]][[1]]
    any_missing <- which(sapply(res[[i]], nrow) != nrow(u))
    for (j in any_missing) {
        v <- res[[i]][[j]]
        m <- match(setdiff(u$id, v$id), u$id)
        filler <- mutate_if(u[m, ], is.numeric, 
            function(u) replace(u, TRUE, NA))
        v <- rbind(v, filler)
        v <- v[match(u$id, v$id), ]
        res[[i]][[j]] <- v
    }
}

cd <- lapply(seq_along(res), function(i) {
    truth <- res[[i]][[1]][, c("id", "is_de", wcs$x)] %>% 
        data.frame(row.names = NULL, check.names = FALSE)
    ps <- lapply(c("p_val", "p_adj.loc"), map, .x = res[[i]]) %>% 
        map(data.frame, check.names = FALSE)
    dfs <- c(list(truth), ps)
    names(dfs) <- c("truth", "pval", "padj")
    do.call(COBRAData, dfs)
})

perf <- lapply(cd, calculate_performance, 
    aspects = "fdrtpr", binary_truth = "is_de", 
    splv = wcs$x, maxsplit = Inf)

df <- map(perf, "fdrtpr") %>% 
    bind_rows(.id = "j") %>% 
    dplyr::select(splitval, thr, method, TPR, FDR) %>% 
    dplyr::filter(splitval != "overall") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%  
    mutate_at("splitval", function(u) {
        u <- gsub(paste0(wcs$x, ":"), "", u)
        v <- sort(unique(as.numeric(u)))
        factor(u, levels = v)
    }) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)


saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)




# Precision-Recall curves
#########################

library(tidyverse)

thresholds = c(1e-18, 1e-15, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3)
thresholds = sort(c(thresholds, seq(5e-3, 1, length.out=500)))

perf2 <- lapply(cd, calculate_performance, 
    aspects = "fdrtpr", binary_truth = "is_de", 
    splv = wcs$x, maxsplit = Inf, thrs=thresholds)

df <- map(perf2, "fdrtpr") %>% 
    bind_rows(.id = "j") %>%
    mutate(Precision = TP / (TP+FP),
        Recall = TP / (TP+FN)) %>%
    mutate(F1 = 2*Precision*Recall / (Precision + Recall)) %>%
    dplyr::select(splitval, thr, method, TPR, FDR, Precision, Recall, F1) %>% 
    dplyr::filter(splitval != "overall") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%  
    mutate_at("splitval", function(u) {
        u <- gsub(paste0(wcs$x, ":"), "", u)
        v <- sort(unique(as.numeric(u)))
        factor(u, levels = v)
    }) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR", "Precision", "Recall", "F1"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

plot_PR = function(df,
    include = "all", color_by = "method", facet = "splitval") {
    df <- filter(ungroup(df), TPR + FDR != 0)
    df$treat <- grepl("treat", df$method)
    if (any(rmv <- table(df$splitval) < 2))
        df <- filter(df, !splitval %in% levels(df$splitval)[rmv])
    df$method <- factor(
        gsub("-treat", "", df$method),
        levels = levels(df$method))
    p <- ggplot(filter(df, FDR + TPR != 0),
        aes_string(x = "Recall", y = "Precision", col = color_by)) +
        facet_wrap(facet, labeller = labeller(.multi_line = FALSE), nrow=1) +
        # geom_vline(size = 0.2, lty = 2, aes(xintercept = thr)) +
        # geom_point(size = 1, alpha = 0.8) +
        geom_line(aes(lty = treat), size = 1, alpha = 0.7, show.legend = (include == "treat")) +
        scale_color_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
        .prettify(theme = "bw", legend.position = "bottom", legend.box.just = "center") +
        guides(
            lty = guide_legend(ncol = 1, keywidth = unit(4, "mm"),
                override.aes = list(alpha = 1, lty = c(1, 3))),
            col = guide_legend(order = 1, nrow = 4,
                override.aes = list(alpha = 1, size = 2)))
    suppressWarnings(suppressMessages(p))
}


p <- plot_PR(df)

file = gsub("\\.pdf", "\\_PR.pdf", args$fig)

ggsave(file, p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)

# AUPR
#---------

library(DescTools)

df_aupr = df %>%
    # filter(splitval==5, method=="limma-voom.sum.counts") %>%
    mutate(Precision = replace_na(Precision, 1)) %>%
    select(splitval, method, Precision, Recall) %>%
    group_by(splitval, method) %>%
    select(-thr) %>% 
    unique %>%
    summarize(AUPR = AUC(Recall, Precision))

plot_AUPR = function(df,
    include = "all", color_by = "method", facet = "splitval") {

    df_aupr %>% 
        ggplot(aes_string("method", "AUPR", fill = color_by)) +
        geom_bar(stat="identity") +
        scale_color_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
        facet_wrap(facet, labeller = labeller(.multi_line = FALSE), nrow=1) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
        scale_fill_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
        .prettify(theme = "bw", legend.position = "bottom", legend.box.just = "center", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        guides(
            lty = guide_legend(ncol = 1, keywidth = unit(4, "mm"),
                override.aes = list(alpha = 1, lty = c(1, 3))),
            col = guide_legend(order = 1, nrow = 4,
                override.aes = list(alpha = 1, size = 2)))
}


p <- plot_AUPR(df_aupr)

file = gsub("\\.pdf", "\\_AUPR.pdf", args$fig)

ggsave(file, p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)

# Max F1 score
##############

plot_maxF1 = function(df,
    include = "all", color_by = "method", facet = "splitval") {
    df <- filter(ungroup(df), TPR + FDR != 0)
    df$treat <- grepl("treat", df$method)
    if (any(rmv <- table(df$splitval) < 2))
        df <- filter(df, !splitval %in% levels(df$splitval)[rmv])
    df$method <- factor(
        gsub("-treat", "", df$method),
        levels = levels(df$method))

    p <- df %>% 
        filter(FDR + TPR != 0) %>%
        group_by(splitval, method) %>%
        summarize(F1max = max(F1, na.rm=TRUE)) %>%
        ggplot(aes_string("method", "F1max", fill = color_by)) +
            geom_bar(stat="identity") +
            scale_color_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
            facet_wrap(facet, labeller = labeller(.multi_line = FALSE), nrow=1) +
            scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
            scale_fill_manual(NULL, values = switch(include, treat = .treat_cols, .meth_cols)) +
            .prettify(theme = "bw", legend.position = "bottom", legend.box.just = "center", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            guides(
                lty = guide_legend(ncol = 1, keywidth = unit(4, "mm"),
                    override.aes = list(alpha = 1, lty = c(1, 3))),
                col = guide_legend(order = 1, nrow = 4,
                    override.aes = list(alpha = 1, size = 2)))

    suppressWarnings(suppressMessages(p))
}


p2 <- plot_maxF1(df)



file = gsub("\\.pdf", "\\_F1max.pdf", args$fig)

ggsave(file, p2,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)








