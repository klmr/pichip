import_package('rlang', attach = TRUE)

mutate_when = function (.data, .filter, ...) {
    dots = dots_definitions(...)$dots
    rows = eval_tidy(enquo(.filter), .data)

    .data[rows, names(dots)] =
        lapply(dots, eval_tidy, data = .data[rows, , drop = FALSE])
    .data
}
