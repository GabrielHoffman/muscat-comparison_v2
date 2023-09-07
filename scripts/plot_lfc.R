suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(purrr)
})

#args <- list(res = list.files("results", "kang,d[a-z][0-9]+,", full.names = TRUE))

df <- .read_res(args$res) %>%
    dplyr::filter(!(is.na(sim_lfc) | is.na(est_lfc))) %>% 
    mutate_at("mid", droplevels) %>% 
    mutate_at("sid", factor, 
        levels = paste0(c("de", "dp", "dm", "db"), "10"),
        labels = c("DE", "DP", "DM", "DB")) 

# dowsample to 2k cluster-gene combinations 
# per method, simulation & gene type
set.seed(29)
sub <- ungroup(sample_n(group_by(df, mid, sid, is_de), 2e3))
sub <- sub[sample(nrow(sub)), ]

p <- ggplot(sub, aes(x = sim_lfc, y = est_lfc, col = as.logical(is_de))) +
    facet_grid(vars(sid), vars(mid)) +
    geom_abline(size = 0.1, slope = 1, intercept = 0)  +
    geom_point_rast(size = 1, alpha = 0.4, raster.dpi = 100) +
    scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "tomato")) +
    guides(color = guide_legend("differential", override.aes = list(size = 3, alpha = 1))) +
    scale_x_continuous(limits = c(-6,6), breaks = seq(-4,4,4), expand = c(0,0)) +
    scale_y_continuous(limits = c(-6,6), breaks = seq(-4,4,4), expand = c(0,0)) +
    labs(x = "simulated logFC", y = "estimated logFC") +
    .prettify("bw") + theme(
        legend.position = "bottom",
        strip.text.x = element_text(size = 3),
        strip.text.y = element_text(size = 6))

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 10, units = "cm",
    dpi = 300, useDingbats = FALSE)


# correlation by group
fig = df %>% 
    group_by(mid, sid) %>% 
    summarize(cor = cor(sim_lfc, est_lfc)) %>% 
    ggplot(aes(mid, cor, fill=mid)) +
        geom_bar(stat="identity") + 
        facet_grid(~ sid)  + 
        theme_classic() +
        theme(aspect.ratio=1, legend.position="none") + 
        scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels= c('0', '0.25', '0.5', '0.75', '1'))  + 
        ylab("Correlation between true and estimated logFC") + 
        xlab("Method") + 
        scale_fill_manual(values = .meth_cols) +
        coord_flip() 

file = gsub(".pdf", "_cor.pdf", args$fig)

ggsave(file=file, fig,
    width = 30, height = 8, units = "cm",
    dpi = 300, useDingbats = FALSE)




