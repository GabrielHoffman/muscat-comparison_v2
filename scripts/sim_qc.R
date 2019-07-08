# load packages
suppressMessages({
    library(countsimQC)
    library(DESeq2)
    library(SummarizedExperiment)
})

# load data
devtools::load_all('../muscat')
sce <- readRDS('data/raw_data/magl.rds')
probs <- list(
    table(sce$cluster_id) / ncol(sce),
    table(sce$sample_id)  / ncol(sce),
    NULL)

# simulate data
set.seed(2903)
sim <- simData(sce,
    n_genes = nrow(sce),
    n_cells = ncol(sce), 
    p_dd = diag(6)[1, ],
    probs = probs)

dds_sim <- DESeqDataSetFromMatrix(
    countData = counts(sim),
    colData = colData(sim),
    design = ~ sample_id)

dds_sce <- DESeqDataSetFromMatrix(
    countData = counts(sce),
    colData = colData(sce),
    design = ~ sample_id)

dds_list <- list(
    sim = dds_sim,
    data = dds_sce)

countsimQCReport(
    ddsList = dds_list,
    outputFile = '/users/helena/Desktop/countsimQC-LPS.html',
    outputDir = ".",
    outputFormat = "html_document",
    maxNForCorr = 200,
    maxNForDisp = 500,
    forceOverwrite = TRUE)