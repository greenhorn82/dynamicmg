example_dmi_data <- function(n = 300) {
    set.seed(123)

    Faktor1 <- rnorm(n, mean = 0, sd = 1)
    Faktor2 <- rnorm(n, mean = 0, sd = 1)

    Item1_F1 <- 0.8 * Faktor1 + rnorm(n, 0, 0.5)
    Item2_F1 <- 0.7 * Faktor1 + rnorm(n, 0, 0.5)
    Item3_F1 <- 0.9 * Faktor1 + rnorm(n, 0, 0.5)

    Item1_F2 <- 0.8 * Faktor2 + rnorm(n, 0, 0.5)
    Item2_F2 <- 0.7 * Faktor2 + rnorm(n, 0, 0.5)
    Item3_F2 <- 0.9 * Faktor2 + rnorm(n, 0, 0.5)

    F3 <- rep(c(0, 1), each = n / 2)

    df <- data.frame(
        Item1_F1,
        Item2_F1,
        Item3_F1,
        Item1_F2,
        Item2_F2,
        Item3_F2,
        F3
    )
}


test_that("Normal run - uncorrelated factors", {
    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")

    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }

    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Reps = 1000,
        FacCor = FALSE,
        Inv = 3
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run1.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})

test_that("Normal run - uncorrelated factors", {
    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")
    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }


    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Reps = 1000,
        FacCor = FALSE,
        Inv = 3
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run1.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})

test_that("Fast run - correlated factors", {
    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }

    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")


    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Inv = 3,
        Reps = 100
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run2.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})

test_that("Fast run - correlated factors + correllated error terms", {
    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }

    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")


    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Inv = 3,
        Reps = 100,
        ResCov = list(c("Item1_F1", "Item1_F2"))
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run3.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})

test_that("Fast run - correlated factors + two correllated error terms", {
    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }

    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")


    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Inv = 3,
        Reps = 100,
        ResCov = list(
            c("Item1_F1", "Item1_F2"),
            c("Item2_F1", "Item3_F2")
        )
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run4.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})


test_that("Three Group model", {
    tol <- if (Sys.info()[["sysname"]] == "Windows") {
        1e-6
    } else {
        .1
    }

    testthat::skip_if_not_installed("lavaan")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("tidyr")
    testthat::skip_if_not_installed("purrr")
    testthat::skip_if_not_installed("MASS")
    testthat::skip_if_not_installed("GenOrd")
    testthat::skip_if_not_installed("semPlot")


    df <- read.table(test_path("testdata", "example1.csv"), header = TRUE, sep = ",")
    df$F3 <- rep(c(0, 1, 3), each = 100)
    result <- calcMG(
        data = df,
        loadings = list(
            c("Item1_F1", "Item2_F1", "Item3_F1"),
            c("Item1_F2", "Item2_F2", "Item3_F2")
        ),
        Group = "F3",
        Inv = 3,
        Reps = 100,
        ResCov = list(
            c("Item1_F1", "Item1_F2"),
            c("Item2_F1", "Item3_F2")
        )
    )

    expect_s3_class(result, "MgDynamic")
    expect_equal(result$input$Factors, 2L)
    expect_equal(result$input$Group, "F3")
    expect_true(is.character(result$model))
    expect_named(result, c(
        "input", "call", "model", "fit", "fit_indices", "differences",
        "decision", "dmacs", "cutoffs", "parameter_tables", "plot",
        "outputs", "download"
    ))
    env <- new.env(parent = emptyenv())
    appResults <- load(test_path("testdata", "run4.rdata"), envir = env)
    expect_equal(result$cutoffs, env$Results, tolerance = tol)
})


test_that("calling inner function without input leads to error", {
    df <- example_dmi_data()

    expect_error(calcMgInner(df), "No input given")
})
