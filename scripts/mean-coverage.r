library(rtracklayer)
library(dplyr)
library(tidyr)
modules::import('klmr/functional/lambda')

samples = readr::read_tsv('raw/samples.tsv')

load_coverage = function (filenames) {
    f = function (filename)
        import.bedGraph(filename) %>%
        coverage(weight = 'score')
    lapply(paste0(filenames, '.gc'), f)
}

roll_mean_cov = function (genomecov, window) {
    f = function (cov)
        as.vector(cov) %>%
        zoo::rollapply(window, mean)
    lapply(genomecov, f)
}

mean_coverage = function (coverages) {
    lapply(seq_along(coverages[[1]]),
           i -> rowMeans(sapply(coverages, `[[`, i))) %>%
    list()
}

normalize = function (data, input) {
    # TODO: For each window on each chromosome, divide data by input.
    f = x ~ y -> {
        r = log(x) - log(y)
        r[is.infinite(r) | is.nan(r)] = max(r, na.rm = TRUE)
        exp(r)
    }
    Map(f, data, input)
}

mean_cov = samples %>%
    mutate(Coverage = load_coverage(File)) %>%
    mutate(Coverage = lapply(Coverage, roll_mean_cov, window = 50)) %>%
    group_by(Condition, Factor) %>%
    summarize(Coverage = mean_coverage(Coverage)) %>%
    spread(Factor, Coverage) %>%
    mutate(Norm = normalize(ChIP, input))
