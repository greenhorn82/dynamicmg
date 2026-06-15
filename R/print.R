#' Print a Dynamic Measurement Invariance result
#'
#' @param x An object of class `MgDynamic`.
#' @param ... Additional arguments, currently unused.
#'
#' @export
print.MgDynamic <- function(x, ...) {
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
    }

    if (!is.null(x$cutoffs)) {
        cat("\nDynamic cutoffs:\n")
        print(x$cutoffs)
    }

    invisible(x)
}
