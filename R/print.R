#' Print a Dynamic Measurement Invariance result
#'
#' @param x An object of class `MgDynamic`.
#' @param showPlot Whether to create plots or not.
#' @param ... Additional arguments, currently unused.
#'
#' @export
print.MgDynamic <- function(x, showPlot = TRUE, ...) {
    cat("MgDynamic result\n")
    cat("Estimator:", x$input$est, "\n")
    cat("Factors:", x$input$Factors, "\n")

    if (!is.null(x$input$Group)) {
        cat("Group variable:", x$input$Group, "\n")
    }

    if (!is.null(x$input$Inv)) {
        inv_label <- switch(as.character(x$input$Inv),
            "1" = "Metric",
            "2" = "Metric and scalar",
            "3" = "Metric, scalar and strict",
            as.character(x$input$Inv)
        )
        cat("Invariance level:", inv_label, "\n")
    }

    if (!is.null(x$input$Reps)) {
        cat("Simulation repetitions:", x$input$Reps, "\n")
    }

    cat("\nAvailable components:\n")
    components <- c(
        "fit_indices", "differences", "decision", "dmacs",
        "cutoffs", "parameter_tables", "plot", "outputs", "download"
    )
    available <- components[!vapply(x[components], is.null, logical(1))]
    if (length(available) > 0) {
        cat(" ", paste(available, collapse = ", "), "\n", sep = "")
    } else {
        cat(" none\n")
    }

    if (!is.null(x$fit_indices)) {
        cat("\nFit indices:\n")
        print(x$fit_indices)
        cat("\n\u0394Fit indices:\n ")
        print(diff(x$fit_indices))
    }

    if (!is.null(x$cutoffs)) {
        cat("\nDynamic cutoffs:\n")
        print(x$cutoffs)
    }

    if (!(is.null(x$cutoffs) | is.null(x$fit_indices))) {
        cat("\nComparision\n")
        x$cutoffs |>
            as.data.frame() |>
            tibble::rownames_to_column(var = "measure") |>
            mutate(measure = stringr::str_replace(measure, "\u0394", "")) -> dataCutoff

        diff(x$fit_indices) |>
            t() |>
            as.data.frame() |>
            tibble::rownames_to_column(var = "measure") |>
            dplyr::filter(measure %in% dataCutoff$measure) -> dataFit

        dataCutoff |>
            dplyr::filter(measure %in% dataFit$measure) -> dataCutoff

        print(dataFit |>
            select(measure) |>
            dplyr::bind_cols(
                dataFit |>
                    select(-measure) |>
                    dplyr::mutate(across(
                        dplyr::everything(),
                        ~ dplyr::if_else(.x > dataCutoff[[dplyr::cur_column()]], "above", "below")
                    ))
            ))
    }

    if (showPlot) {
        semPlot::semPaths(x$fit, residuals = FALSE, intercepts = FALSE, thresholds = FALSE)
    }

    invisible(x)
}
