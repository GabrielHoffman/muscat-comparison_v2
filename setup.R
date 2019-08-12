config <- yaml::read_yaml("config.yaml")
x <- c("sim_pars", "run_pars", "meth_pars", "sim_data", "results", "figures")
x <- c(unlist(config[x]), 
    file.path(config$results, config$dids), 
    file.path(config$figures, config$dids))
for (dir in x)
    if (!dir.exists(dir) & !isTRUE(grep("\\.", dir) == 1))
        dir.create(dir, recursive = TRUE)

scripts <- file.path(config$scripts,
    paste0(c("sim_pars", "run_pars", "meth_pars"), ".R"))
sapply(scripts, source)