.as_calc_mg_yes_no <- function(x) {
    if (is.logical(x)) {
        return(if (isTRUE(x)) "Yes" else "No")
    }

    if (is.numeric(x)) {
        return(if (x == 1) "Yes" else "No")
    }

    if (tolower(x) %in% c("yes", "y", "true", "1")) {
        return("Yes")
    }

    "No"
}


.build_calc_mg_input <- function(loadings = NULL,
                                 GenFac = FALSE, FacCor = TRUE,
                                 est = "ML", Group = NULL,
                                 Inv = 1, Scale = "N", Reps = 500,
                                 ResCov = NULL,
                                 missing = NULL,
                                 ...) {
    is_non_empty <- function(x) {
        !is.null(x) && length(x) > 0
    }


    if (!is.null(loadings)) {
        if (!is.list(loadings)) {
            loadings <- list(loadings)
        }
        factor_values <- loadings
    } else {
        stop("A list with which items loads on which factor is requiered.", call. = FALSE)
    }

    if (length(factor_values) > 8) {
        stop("A maximum of eight factors is supported.", call. = FALSE)
    }

    input <- list()
    for (i in seq_len(8)) {
        input[[paste0("Factor", i)]] <- if (i <= length(factor_values)) {
            factor_values[[i]]
        } else {
            vector()
        }
    }

    input$Factors <- length(factor_values)
    if (is.na(input$Factors) || input$Factors < 1 || input$Factors > 8) {
        stop("Factors must be an integer between 1 and 8.", call. = FALSE)
    }

    input$GenFac <- if (isTRUE(GenFac) || identical(GenFac, 1)) 1 else 0
    input$FacCor <- .as_calc_mg_yes_no(FacCor)
    input$est <- est
    input$Group <- Group
    input$Inv <- as.integer(Inv)
    input$Scale <- match.arg(Scale, c("N", "L"))
    input$Reps <- Reps
    input$missing <- missing

    rc_args <- ResCov
    rc_values <- rc_args[vapply(rc_args, is_non_empty, logical(1))]


    input$ResCov <- length(rc_values)
    if (is.na(input$ResCov) || input$ResCov < 0) {
        stop("ResCov must be a non-negative integer.", call. = FALSE)
    }
    if (input$ResCov > length(rc_values)) {
        stop("Provide one residual covariance vector for each requested ResCov.", call. = FALSE)
    }

    for (i in seq_len(input$ResCov)) {
        input[[paste0("RC", i)]] <- rc_values[[i]]
    }

    extra <- list(...)
    if (length(extra) > 0) {
        input[names(extra)] <- extra
    }

    input
}


.new_MgDynamic <- function(input, model, fit,
                           fit_indices = NULL,
                           differences = NULL,
                           decision = NULL,
                           dmacs = NULL,
                           cutoffs = NULL,
                           parameter_tables = NULL,
                           plot = NULL,
                           outputs = NULL,
                           download = NULL,
                           call = NULL) {
    structure(
        list(
            input = input,
            call = call,
            model = model,
            fit = fit,
            fit_indices = fit_indices,
            differences = differences,
            decision = decision,
            dmacs = dmacs,
            cutoffs = cutoffs,
            parameter_tables = parameter_tables,
            plot = plot,
            outputs = outputs,
            download = download
        ),
        class = "MgDynamic"
    )
}


#' Calculate Dynamic Measurement Invariance Cutoffs
#'
#' `calcMG()` estimates a two-group confirmatory factor analysis model and
#' simulates model-specific cutoff values for measurement invariance tests.
#' The function was adapted from the Dynamic Measurement Invariance Shiny app
#' and returns the results as an object of class `MgDynamic`.
#'
#' @param data A `data.frame` containing the observed item variables and the
#'   grouping variable.
#' @param loadings A list describing the measurement model. Each element must be
#'   a character vector with the item names loading on one latent factor. Up to
#'   eight factors are supported.
#' @param GenFac Logical. If `TRUE`, a general higher-order factor is added for
#'   models with more than two factors. Defaults to `FALSE`.
#' @param FacCor Logical. If `TRUE`, latent factors are allowed to correlate. If
#'   `FALSE`, factors are estimated as orthogonal. Defaults to `TRUE`.
#' @param est Character string specifying the lavaan estimator. Currently used
#'   values are `"ML"` and `"MLR"`.
#' @param Group Character string naming the grouping variable in `data`.
#' @param Inv Integer indicating the highest invariance level to evaluate:
#'   `1` for metric, `2` for metric and scalar, and `3` for metric, scalar, and
#'   strict invariance.
#' @param Scale Character string indicating the response scale. Use `"N"` for
#'   continuous normally distributed indicators and `"L"` for Likert/ordinal
#'   indicators treated as continuous.
#' @param Reps Number of simulation replications used to estimate the dynamic
#'   cutoff values. Larger values give more stable cutoffs but increase runtime.
#' @param ResCov Optional list of residual covariances. Each element should be a
#'   character vector of length two naming the two item residuals to correlate,
#'   for example `list(c("Item1_F1", "Item1_F2"))`.
#' @param missing Optional missing-data indicator retained for compatibility with
#'   the original Shiny app workflow.
#' @param ... Additional values passed to the internal input list.
#'
#' @return An object of class `MgDynamic`, a list containing the generated
#'   lavaan model syntax, fitted lavaan object, model fit indices, fit-index
#'   differences, D-MACS values, simulated cutoffs, parameter tables, plot object,
#'   and the original input settings.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 300
#'
#' Faktor1 <- rnorm(n)
#' Faktor2 <- rnorm(n)
#'
#' df <- data.frame(
#'     Item1_F1 = 0.8 * Faktor1 + rnorm(n, 0, 0.5),
#'     Item2_F1 = 0.7 * Faktor1 + rnorm(n, 0, 0.5),
#'     Item3_F1 = 0.9 * Faktor1 + rnorm(n, 0, 0.5),
#'     Item1_F2 = 0.8 * Faktor2 + rnorm(n, 0, 0.5),
#'     Item2_F2 = 0.7 * Faktor2 + rnorm(n, 0, 0.5),
#'     Item3_F2 = 0.9 * Faktor2 + rnorm(n, 0, 0.5),
#'     F3 = rep(c(0, 1), each = n / 2)
#' )
#'
#' result <- calcMG(
#'     data = df,
#'     loadings = list(
#'         c("Item1_F1", "Item2_F1", "Item3_F1"),
#'         c("Item1_F2", "Item2_F2", "Item3_F2")
#'     ),
#'     Group = "F3",
#'     Inv = 3,
#'     Reps = 100
#' )
#'
#' print(result)
#'
#' result_with_residual_covariance <- calcMG(
#'     data = df,
#'     loadings = list(
#'         c("Item1_F1", "Item2_F1", "Item3_F1"),
#'         c("Item1_F2", "Item2_F2", "Item3_F2")
#'     ),
#'     Group = "F3",
#'     Inv = 3,
#'     Reps = 100,
#'     ResCov = list(c("Item1_F1", "Item1_F2"))
#' )
#' }
#'
#' @importFrom lavaan anova cfa fitMeasures fitmeasures lavInspect lavNames parameterEstimates partable lav_partable_npar
#' @importFrom dplyr arrange desc filter group_by mutate pull reframe select
#' @importFrom tidyr nest unite
#' @importFrom purrr map map_dfr pmap_dfr
#' @importFrom magrittr %>%
#' @importFrom MASS mvrnorm
#' @importFrom GenOrd ordcont ordsample
#' @importFrom semPlot semPaths
#' @importFrom rmarkdown html_document render
#' @export

calcMG <- function(data, loadings = NULL,
                   GenFac = FALSE, FacCor = TRUE,
                   est = "ML", Group = NULL,
                   Inv = 1, Scale = "N", Reps = 500,
                   ResCov = NULL,
                   missing = NULL,
                   ...) {
    input <- .build_calc_mg_input(
        loadings = loadings,
        GenFac = GenFac, FacCor = FacCor,
        est = est, Group = Group,
        Inv = Inv, Scale = Scale, Reps = Reps,
        ResCov = ResCov,
        missing = missing,
        ...
    )


    result <- calcMgInner(data, input)
    result$call <- match.call()
    result
}


calcMgInner <- function(data, input) {
    if (missing(input) || is.null(input)) {
        stop("No input given. Either use calcMg with parameters or give a
             list with an input object")
    }

    # create list of items on each factor
    # 8 are listed because its the maximum currently allowed
    # easier to list out all 8 that to adaptively change the suffix after a dollar sign
    l <- list(
        input$Factor1, input$Factor2,
        input$Factor3, input$Factor4,
        input$Factor5, input$Factor6,
        input$Factor7, input$Factor8
    )
    output <- list()
    # l <- loadings()
    # left side
    lhs <- list()
    # right side
    rhs <- list()
    # combined
    line <- list()

    if (input$Factors < 3) {
        # loop over factors
        for (m in 1:input$Factors) {
            # left side is factor name
            lhs[[m]] <- paste0("f", m, "=~")
            # begin right side with first item
            rhs[[m]] <- l[[m]][1]

            # loop over number of items
            for (i in 1:(length(l[[m]]) - 1)) {
                # plus sign between item names
                rhs[[m]] <- paste(rhs[[m]], "+", l[[m]][i + 1])
            }

            # combine left and right hand side
            line[[m]] <- paste(lhs[[m]], rhs[[m]])
        }

        ### add residual covariances here###
        if (input$ResCov > 0) {
            for (r in 1:input$ResCov) {
                line[[input$Factors + r]] <- paste(eval(parse(text = paste0("input$RC", r, "[1]"))), "~~", eval(parse(text = paste0("input$RC", r, "[2]"))))
            }
        }
    }

    if (input$Factors > 2) {
        # loop over factors
        for (m in 1:input$Factors) {
            # left side is factor name
            lhs[[m]] <- paste0("f", m, "=~")
            # begin right side with first item
            rhs[[m]] <- l[[m]][1]

            # loop over number of items
            for (i in 1:(length(l[[m]]) - 1)) {
                # plus sign between item names
                rhs[[m]] <- paste(rhs[[m]], "+", l[[m]][i + 1])
            }

            # combine left and right hand side
            line[[m]] <- paste(lhs[[m]], rhs[[m]])
        }

        # if hierarchical factor is present, add line to statement has load general factor 'g' on all factors
        if (as.numeric(input$GenFac == 1)) {
            line[[input$Factors + 1]] <- paste0("gen", "=~", "f1")
            for (m in 2:input$Factors) {
                line[[input$Factors + 1]] <- paste0(line[[input$Factors + 1]], " + ", "f", m)
            }
        }

        ### add residual covariances here###
        if (input$ResCov > 0) {
            for (r in 1:input$ResCov) {
                line[[input$Factors + as.numeric(input$GenFac) + r]] <- paste(eval(parse(text = paste0("input$RC", r, "[1]"))), "~~", eval(parse(text = paste0("input$RC", r, "[2]"))))
            }
        }
    }

    # unlist to put into one model statement
    model <- unlist(line)


    # T/F indicator for correlated factors in lavaan
    if (input$Factors == 1) {
        COR <- TRUE
    } else {
        COR <- ifelse(input$FacCor == "Yes", FALSE, TRUE)
    }

    if (input$est %in% c("ML", "MLR")) {
        FIML <- "ml"
    } else {
        FIML <- "listwise"
    }

    # Fit model in lavaan
    a <- lavaan::cfa(
        data = data, model = model,
        estimator = input$est, orthogonal = COR, std.lv = TRUE, missing = FIML
    )
    lav0 <- lavaan::partable(a)
    df0 <- a@test$standard$df
    imp0 <- a@implied$cov[[1]]
    n0 <- a@Data@nobs[[1]]


    #######################################
    ###### ML vs MLR for fit measures ######
    #######################################
    if (input$est == "ML") {
        a1 <- lavaan::cfa(data = data, model = model, estimator = input$est, group = input$Group, orthogonal = COR, std.lv = TRUE)
        indm1 <- t(as.data.frame(lavaan::fitmeasures(a1, c("rmsea", "srmr", "cfi", "mfi", "chisq", "df", "pvalue"))))
        colnames(indm1) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
        rownames(indm1) <- ("Configural")
        chi0 <- a@test$standard$stat
    }

    if (input$est == "MLR") {
        a1 <- lavaan::cfa(data = data, model = model, estimator = input$est, group = input$Group, orthogonal = COR, std.lv = TRUE)
        indm1 <- t(as.data.frame(lavaan::fitmeasures(a1, c("rmsea.scaled", "srmr", "cfi.scaled", "mfi", "chisq.scaled", "df.scaled", "pvalue.scaled"))))
        # lavaan doesn't have scaled version of mfi, so calculate it manually
        indm1[4] <- exp(-.5 * ((indm1[5] - indm1[6]) / (n0 - 1)))
        colnames(indm1) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
        rownames(indm1) <- ("Configural")
        chi0 <- a@test$yuan.bentler.mplus$stat
    }

    ###########################################################################
    #### code for computing d_macs ####
    ## from Dueber (2020) ##
    # loading package was a mess because
    # of dependcies that are not needed if reading in a lavaan object #
    ############################################################################

    colSD <- function(x, ...) {
        apply(X = x, MARGIN = 2, FUN = sd, ...)
    }

    expected_value <- function(Lambda, Nu, Eta, Thresh = NULL, Theta = NULL, categorical = FALSE) {
        if (categorical) {
            ## Graded Response model with probit link.

            ## Let's give ourselves a y* score
            Mu <- Nu + Lambda * Eta

            ## categories go from 0 to number of thresholds
            max <- length(Thresh)
            ## Make a max+1 category that is impossible to attain
            Thresh[max + 1] <- Inf

            ## Start expected value at 0 and incremented with expectation value from each category
            expected <- 0
            for (i in 1:max) {
                expected <- expected + i * (pnorm(Mu - Thresh[i], sd = sqrt(Theta)) -
                    pnorm(Mu - Thresh[i + 1], sd = sqrt(Theta)))
            }
            expected
        } else {
            ## The linear continuous world is so easy!!
            Nu + Lambda * Eta
        }
    }


    item_dmacs <- function(LambdaR, LambdaF,
                           NuR, NuF,
                           MeanF, VarF, SD,
                           ThreshR = NULL, ThreshF = NULL,
                           ThetaR = NULL, ThetaF = NULL,
                           categorical = FALSE) {
        # Use Thresholds as a check for categorical-ness
        if (!is.null(ThreshR)) {
            categorical <- TRUE
            ## If threshold vectors do not have the same length, throw an error
            if (length(ThreshR) != length(ThreshF)) stop("Item must have same number of thresholds in both reference and focal group")
        }

        ## If item does not load on factor, return NA
        if (LambdaR == 0) {
            return(NA)
        }

        ## Create a function for the integrand using the expected value function expected_value
        ## The change of variable math replaces f_F(Xi) with dnorm(z)

        integrand <- function(z, LambdaR, LambdaF,
                              NuR, NuF,
                              ThreshR, ThreshF,
                              ThetaR, ThetaF,
                              MeanF, VarF, categorical) {
            (expected_value(LambdaF, NuF, MeanF + z * sqrt(VarF), ThreshF, ThetaF, categorical) -
                expected_value(LambdaR, NuR, MeanF + z * sqrt(VarF), ThreshR, ThetaR, categorical))^2 * dnorm(z)
        }

        ## Now, sum it to get the integral, and compute the effect size. Stepsize is in z units, not theta units!!
        sqrt(integrate(integrand, -Inf, Inf,
            LambdaR, LambdaF,
            NuR, NuF,
            ThreshR, ThreshF,
            ThetaR, ThetaF,
            MeanF, VarF,
            categorical = categorical
        )$value) / SD
    }

    dmacs_summary <- function(LambdaList, NuList,
                              MeanList, VarList, SDList,
                              Groups = NULL, RefGroup = 1,
                              ThreshList = NULL, ThetaList = NULL,
                              categorical = FALSE) {
        ## See if we need to get group names, and if we do, try to grab them from the names of LambdaList. Otherwise, just number the groups
        if (is.null(Groups)) {
            if (is.null(names(LambdaList))) {
                Groups <- c(1:length(LambdaList))
            } else {
                Groups <- names(LambdaList)
            }
        }

        ## If RefGroup is a string, lets turn it into an index
        if (is.character(RefGroup)) {
            RefGroup <- match(RefGroup, Groups)
        }

        # The categorical and continuous cases are different from each other
        if (categorical) { # now we are categorical
            ## if only two groups, then call DIF effect summary single right away, else iterate over the focal groups
            if (length(Groups) == 2) {
                dmacs_summary_single(
                    LambdaF = LambdaList[-RefGroup][[1]],
                    NuF = NuList[-RefGroup][[1]],
                    ThreshF = ThreshList[-RefGroup][[1]],
                    ThetaF = ThetaList[-RefGroup][[1]],
                    MeanF = MeanList[-RefGroup][[1]],
                    VarF = VarList[-RefGroup][[1]],
                    SD = SDList[-RefGroup][[1]],
                    LambdaR = LambdaList[[RefGroup]],
                    NuR = NuList[[RefGroup]],
                    ThreshR = ThreshList[[RefGroup]],
                    ThetaR = ThetaList[[RefGroup]],
                    categorical = categorical
                )
            } else {
                mapply(dmacs_summary_single,
                    LambdaF = LambdaList[-RefGroup],
                    NuF = NuList[-RefGroup],
                    ThreshF = ThreshList[-RefGroup],
                    ThetaF = ThetaList[-RefGroup],
                    MeanF = MeanList[-RefGroup],
                    VarF = VarList[-RefGroup],
                    SD = SDList[-RefGroup],
                    MoreArgs = list(
                        LambdaR = LambdaList[[RefGroup]],
                        NuR = NuList[[RefGroup]],
                        ThreshR = ThreshList[[RefGroup]],
                        ThetaR = ThetaList[[RefGroup]],
                        categorical = categorical
                    ),
                    SIMPLIFY = FALSE
                )
            }
        } else { # Continuous indicators
            ## if only two groups, then call DIF effect summary single right away, else iterate over the focal groups
            if (length(Groups) == 2) {
                dmacs_summary_single(
                    LambdaF = LambdaList[-RefGroup][[1]],
                    NuF = NuList[-RefGroup][[1]],
                    MeanF = MeanList[-RefGroup][[1]],
                    VarF = VarList[-RefGroup][[1]],
                    SD = SDList[-RefGroup][[1]],
                    LambdaR = LambdaList[[RefGroup]],
                    NuR = NuList[[RefGroup]],
                    categorical = categorical
                )
            } else {
                mapply(dmacs_summary_single,
                    LambdaF = LambdaList[-RefGroup],
                    NuF = NuList[-RefGroup],
                    MeanF = MeanList[-RefGroup],
                    VarF = VarList[-RefGroup],
                    SD = SDList[-RefGroup],
                    MoreArgs = list(
                        LambdaR = LambdaList[[RefGroup]],
                        NuR = NuList[[RefGroup]],
                        categorical = categorical
                    ),
                    SIMPLIFY = FALSE
                )
            }
        }
    }

    dmacs_summary_single <- function(LambdaR, LambdaF,
                                     NuR, NuF,
                                     MeanF, VarF, SD,
                                     ThreshR = NULL, ThreshF = NULL,
                                     ThetaR = NULL, ThetaF = NULL,
                                     categorical = FALSE) {
        ## Categorical and continuous work a bit differently from each other
        if (categorical) { # Now we are categorical
            categorical <- TRUE
            if (!is.list(ThreshR)) stop("Thresholds must be in a list indexed by item. The thresholds for each item should be a vector")

            ## If unidimensional, then things are straightforward, otherwise not so much!!
            if (ncol(LambdaR) == 1) {
                DMACS <- mapply(
                    item_dmacs,
                    LambdaR, LambdaF,
                    NuR, NuF,
                    MeanF, VarF, SD,
                    ThreshR, ThreshF,
                    ThetaR, ThetaF,
                    categorical
                )
                names(DMACS) <- rownames(LambdaR)

                ItemDeltaMean <- mapply(
                    delta_mean_item,
                    LambdaR, LambdaF,
                    NuR, NuF,
                    MeanF, VarF,
                    ThreshR, ThreshF,
                    ThetaR, ThetaF,
                    categorical
                )
                names(ItemDeltaMean) <- rownames(LambdaR)

                MeanDiff <- sum(ItemDeltaMean, na.rm = TRUE)
                names(MeanDiff) <- colnames(LambdaR)

                list(DMACS = DMACS, ItemDeltaMean = ItemDeltaMean, MeanDiff = MeanDiff)
            } else {
                ## Need to give MeanF and VarF (which are vectors indexed by factor) the same structure as LambdaR (an array indexed by itemsxfactors)
                MeanF <- as.vector(MeanF)
                MeanF <- matrix(rep(MeanF, nrow(LambdaR)), nrow = nrow(LambdaR), byrow = TRUE)
                VarF <- as.vector(VarF)
                VarF <- matrix(rep(VarF, nrow(LambdaR)), nrow = nrow(LambdaR), byrow = TRUE)

                DMACS <- as.data.frame(matrix(
                    mapply(
                        item_dmacs,
                        LambdaR, LambdaF,
                        NuR, NuF,
                        MeanF, VarF, SD,
                        ThreshR, ThreshF,
                        ThetaR, ThetaF,
                        categorical
                    ),
                    nrow = nrow(LambdaR)
                ))
                colnames(DMACS) <- colnames(LambdaR)
                rownames(DMACS) <- rownames(LambdaR)


                ## ItemDeltaMean has the same possible issues as DMACS
                ItemDeltaMean <- as.data.frame(matrix(
                    mapply(
                        delta_mean_item,
                        LambdaR, LambdaF,
                        NuR, NuF,
                        MeanF, VarF,
                        ThreshR, ThreshF,
                        ThetaR, ThetaF,
                        categorical
                    ),
                    nrow = nrow(LambdaR)
                ))
                colnames(ItemDeltaMean) <- colnames(LambdaR)
                rownames(ItemDeltaMean) <- rownames(LambdaR)

                MeanDiff <- colSums(ItemDeltaMean, na.rm = TRUE)

                list(DMACS = DMACS, ItemDeltaMean = ItemDeltaMean, MeanDiff = MeanDiff)
            }
        } else { # Now we are continuous
            ## If unidimensional, then things are straightforward, otherwise not so much!!
            if (ncol(LambdaR) == 1) {
                DMACS <- mapply(item_dmacs, LambdaR, LambdaF,
                    NuR, NuF,
                    MeanF, VarF, SD,
                    categorical = FALSE
                )
                names(DMACS) <- rownames(LambdaR)

                ItemDeltaMean <- mapply(delta_mean_item, LambdaR, LambdaF,
                    NuR, NuF,
                    MeanF, VarF,
                    categorical = FALSE
                )
                names(ItemDeltaMean) <- rownames(LambdaR)

                MeanDiff <- sum(ItemDeltaMean, na.rm = TRUE)
                names(MeanDiff) <- colnames(LambdaR)

                VarDiff <- delta_var(LambdaR, LambdaF, VarF)
                names(VarDiff) <- colnames(LambdaR)
                list(DMACS = DMACS, ItemDeltaMean = ItemDeltaMean, MeanDiff = MeanDiff, VarDiff = VarDiff)
            } else {
                ## Need to give MeanF and VarF (which are vectors indexed by factor) the same structure as LambdaR (an array indexed by itemsxfactors)
                MeanF <- as.vector(MeanF)
                MeanF <- matrix(rep(MeanF, nrow(LambdaR)), nrow = nrow(LambdaR), byrow = TRUE)
                VarF <- as.vector(VarF)
                VarF <- matrix(rep(VarF, nrow(LambdaR)), nrow = nrow(LambdaR), byrow = TRUE)

                DMACS <- as.data.frame(matrix(
                    mapply(item_dmacs,
                        LambdaR, LambdaF,
                        NuR, NuF,
                        MeanF, VarF, SD,
                        categorical = FALSE
                    ),
                    nrow = nrow(LambdaR)
                ))
                colnames(DMACS) <- colnames(LambdaR)
                rownames(DMACS) <- rownames(LambdaR)


                ## ItemDeltaMean has the same possible issues as DMACS
                ItemDeltaMean <- as.data.frame(matrix(
                    mapply(delta_mean_item,
                        LambdaR, LambdaF,
                        NuR, NuF,
                        MeanF, VarF,
                        categorical = FALSE
                    ),
                    nrow = nrow(LambdaR)
                ))
                colnames(ItemDeltaMean) <- colnames(LambdaR)
                rownames(ItemDeltaMean) <- rownames(LambdaR)

                MeanDiff <- colSums(ItemDeltaMean, na.rm = TRUE)

                ## delta_var needs to be redesigned for multidimensional models, so let's leave it off for now
                # VarDiff <- delta_var(LambdaR, LambdaF, VarF)


                list(DMACS = DMACS, ItemDeltaMean = ItemDeltaMean, MeanDiff = MeanDiff) # , VarDiff = VarDiff)
            }
        }
    }

    delta_mean_item <- function(LambdaR, LambdaF,
                                NuR, NuF,
                                MeanF, VarF,
                                ThreshR = NULL, ThreshF = NULL,
                                ThetaR = NULL, ThetaF = NULL,
                                categorical = FALSE) {
        # Use Thresholds as a check for categorical-ness
        if (categorical) {
            categorical <- TRUE
            ## If threshold vectors do not have the same length, throw an error
            if (length(ThreshR) != length(ThreshF)) stop("Item must have same number of thresholds in both reference and focal group.")
        }

        ## If item does not load on factor, return NA
        if (LambdaR == 0) {
            return(NA)
        }

        ## Create a function for the integrand using the expected value function expected_value
        ## The change of variable math replaces f_F(Xi) with dnorm(z)
        integrand <- function(z, LambdaR, LambdaF,
                              NuR, NuF,
                              ThreshR, ThreshF,
                              ThetaR, ThetaF,
                              MeanF, VarF, categorical) {
            (expected_value(LambdaF, NuF, MeanF + z * sqrt(VarF), ThreshF, ThetaF, categorical) -
                expected_value(LambdaR, NuR, MeanF + z * sqrt(VarF), ThreshR, ThetaR, categorical)) * dnorm(z)
        }
        ## Now, integrate
        integrate(integrand, -Inf, Inf,
            LambdaR, LambdaF,
            NuR, NuF,
            ThreshR, ThreshF,
            ThetaR, ThetaF,
            MeanF, VarF,
            categorical = categorical
        )$value
    }

    delta_var <- function(LambdaR, LambdaF, VarF, categorical = FALSE) {
        if (categorical) {
            warning("At this time, delta variance can only be computed for linear models, not for categorical ones")
            return(NULL)
        }
        delta_cov_mat <- matrix(nrow = length(LambdaR), ncol = length(LambdaR))
        ## I know for loops are supposed to be bad, but this is SO CLEAN!
        for (i in 1:length(LambdaR)) {
            for (j in 1:length(LambdaR)) {
                delta_cov_mat[i, j] <- LambdaR[[j]] * (LambdaF[[i]] - LambdaR[[i]]) * VarF
                +LambdaR[[i]] * (LambdaF[[j]] - LambdaR[[j]]) * VarF
                +(LambdaF[[i]] - LambdaR[[i]]) * (LambdaF[[j]] - LambdaR[[j]]) * VarF
            }
        }
        sum(delta_cov_mat)
    }

    lavaan_dmacs <- function(fit, RefGroup = 1, dtype = "pooled", MEtype = "Group") {
        if (grepl("ong", MEtype, fixed = TRUE)) { # Long, Longitudinal, long, longitudinal
            ## Groups are time-points. We ignore correlated residuals!

            # Make a vector of factor names
            Groups <- colnames(lavaan::lavInspect(fit, what = "est")$lambda)

            # If RefTime is a name, turn it into an index
            if (is.character(RefGroup)) {
                RefGroup <- match(RefGroup, Groups)
            }

            # Store the estimates and the data, because I am about to reference them a LOT of times
            FitEst <- lavaan::lavInspect(fit, "est")
            FitData <- lavaan::lavInspect(fit, "data")

            ## factor loadings, item intercepts, factor means, and factor variances are easy
            LambdaList <- lapply(Groups, function(x) {
                Lambdas <- FitEst$lambda[FitEst$lambda[, x] != 0, x]
                matrix(Lambdas, ncol = 1, dimnames = list(names(Lambdas)))
            })
            names(LambdaList) <- Groups

            NuList <- lapply(1:length(Groups), function(x) {
                if (is.null(FitEst$nu)) {
                    # fill in zeros, because if intercepts are not in the model, they are automatically zero
                    rep(0, length(rownames(LambdaList[[x]])))
                } else {
                    FitEst$nu[rownames(LambdaList[[x]]), ]
                }
            })
            names(NuList) <- Groups

            MeanList <- lapply(Groups, function(x) {
                if (is.null(FitEst$alpha)) {
                    # If factor mean is not mentioned in the model, it must be zero!
                    0
                } else {
                    FitEst$alpha[x, 1]
                }
            })
            names(MeanList) <- Groups

            VarList <- lapply(Groups, function(x) {
                FitEst$psi[x, x]
            })
            names(VarList) <- Groups

            ## compute the sds for use in Equation 3 of Nye and Drasgow (2011)
            if (dtype == "pooled") {
                refsd <- colSD(FitData[, rownames(LambdaList[[RefGroup]])], na.rm = TRUE)
                refn <- colSums(!is.na(FitData[, rownames(LambdaList[[RefGroup]])]))
                SDList <- lapply(1:length(Groups), function(x) {
                    focsd <- colSD(FitData[, rownames(LambdaList[[x]])], na.rm = TRUE)
                    focn <- colSums(!is.na(FitData[, rownames(LambdaList[[x]])]))
                    ((focn - 1) * focsd + (refn - 1) * refsd) / ((focn - 1) + (refn - 1))
                })
            } else if (dtype == "glass") { ## Glass says to always use the SD of the reference group
                SDs <- colSD(FitData[, rownames(LambdaList[[RefGroup]])], na.rm = TRUE)
                SDList <- lapply(1:length(Groups), function(x) {
                    SDs
                })
                names(SDList) <- Groups
            } else {
                stop("Only \"pooled\" and \"glass\" SD types are supported")
            }

            ## Check to see if we are using categorical or linear variables, because Thresh and Theta only apply to categorical
            if (length(lavaan::lavNames(fit, type = "ov.ord")) == 0) {
                categorical <- FALSE
                ThreshList <- NULL
                ThetaList <- NULL
            } else {
                categorical <- TRUE

                ## Make a list of thresholds indexed by group
                ThreshList <- lapply(1:length(Groups), function(x) {
                    # Fetch indicator names so we can grepl them
                    ItemNames <- rownames(LambdaList[[x]])

                    # Return a list index by items
                    lapply(ItemNames, function(y) {
                        # now we need to fetch the thresholds for this item.
                        FitEst$tau[grepl(paste0(y, "\\|"), rownames(FitEst$tau))]
                    })
                })

                ## make a list of residual variances indexed by group
                ThetaList <- lapply(1:length(Groups), function(x) {
                    diag(FitEst$theta)[rownames(LambdaList[[x]])]
                })
            }
        } else {
            # Now we are doing multi-group measurement equivalence testing
            Groups <- names(lavaan::lavInspect(fit, "est"))

            ## If RefGroup is a string, turn it into an index
            if (is.character(RefGroup)) {
                RefGroup <- match(RefGroup, Groups)
            } else {
                warning(paste("It is recommended that you provide the name of the reference group as a string; see ?lavaan_dmacs. The reference group being used is", Groups[RefGroup]))
            }

            ## factor loadings, item intercepts, factor means, and factor variances are easy
            LambdaList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                x$lambda
            })
            NuList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                x$nu
            })
            MeanList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                x$alpha
            })
            VarList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                diag(x$psi)
            })


            ## compute the sds for use in Equation 3 of Nye and Drasgow (2011)
            if (dtype == "pooled") {
                refsd <- colSD(lavaan::lavInspect(fit, "data")[[RefGroup]], na.rm = TRUE)
                refn <- colSums(!is.na(lavaan::lavInspect(fit, "data")[[RefGroup]]))
                SDList <- lapply(lavaan::lavInspect(fit, "data"), function(x) {
                    focsd <- colSD(x, na.rm = TRUE)
                    focn <- colSums(!is.na(x))
                    ((focn - 1) * focsd + (refn - 1) * refsd) / ((focn - 1) + (refn - 1))
                })
            } else if (dtype == "glass") { ## Glass says to always use the SD of the reference group
                SDs <- colSD(lavaan::lavInspect(fit, "data")[[RefGroup]], na.rm = TRUE)
                SDList <- lapply(1:length(Groups), function(x) {
                    SDs
                })
                names(SDList) <- Groups
            } else {
                stop("Only \"pooled\" and \"glass\" SD types are supported")
            }


            ## Check to see if we are using categorical or linear variables, because Thresh works differently in those cases
            if (length(lavaan::lavNames(fit, type = "ov.ord")) == 0) {
                categorical <- FALSE
            } else {
                categorical <- TRUE

                ## Need the item names so we can grepl them
                ItemNames <- rownames(lavaan::lavInspect(fit, "est")[[1]]$lambda)

                ## I don't know why I am not doing this as nested for loops!! Nesting lapply inside of lapply is awful
                ThreshList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                    ## This next line makes a LIST indexed by item, which ensures that the mapply in DIF_effect_summary_single iterates over the thresholds properly
                    lapply(
                        ItemNames,
                        ## The funny paste0 is in case one item name is an extension of another item name (e.g., item10 vs item1)
                        function(iname, threshlist) {
                            threshlist[grepl(paste0(iname, "\\|"), rownames(threshlist))]
                        },
                        x$tau
                    )
                })

                # Now we need to get the thetas, too!!
                ThetaList <- lapply(lavaan::lavInspect(fit, "est"), function(x) {
                    diag(x$theta)
                })
            }
        }


        Results <- dmacs_summary(
            LambdaList, NuList,
            MeanList, VarList, SDList,
            Groups, RefGroup,
            ThreshList, ThetaList,
            categorical
        )


        ## Note to self - we may need to insert some names here!!

        Results
    }
    ###############################
    # DMACS code over #
    ###############################

    suppressWarnings(D <- lavaan_dmacs(a1))
    if (ncol(as.data.frame(D$DMACS)) == 1) {
        D_new <- as.data.frame(D$DMACS)
        colnames(D_new) <- "f1"
    }

    if (ncol(as.data.frame(D$DMACS)) > 1) {
        D_new <- D$DMACS
    }

    output$Dmacs <- D_new


    ##########################################################################
    ## ADD ROWS TO TABLE THAT TAKE DIFFERENCE BETWEEN CURRENT AND PREVIOUS ROW#
    ##########################################################################

    if (input$Inv > 0) {
        fit2 <- lavaan::cfa(data = data, model = model, estimator = input$est, orthogonal = COR, std.lv = TRUE, group = input$Group, group.equal = "loadings")

        q <- lavaan::parameterEstimates(fit2)

        names(q)[names(q) == "est"] <- "Estimate"
        names(q)[names(q) == "se"] <- "SE"
        names(q)[names(q) == "z"] <- "Z"
        names(q)[names(q) == "pvalue"] <- "p"
        names(q)[names(q) == "ci.lower"] <- "95% CI Lower Limit"
        names(q)[names(q) == "ci.upper"] <- "95% CI Upper Limit"

        q1 <- q %>% dplyr::filter(group == 1)
        q2 <- q %>% dplyr::filter(group == 2)

        q1load <- q1 %>%
            dplyr::filter(op == "=~")
        q1load <- q1load[, -c(2, 4:6)]
        names(q1load)[names(q1load) == "lhs"] <- "Factor"
        names(q1load)[names(q1load) == "rhs"] <- "Item"

        q1int <- q1 %>%
            dplyr::filter(op == "~1")
        q1int <- q1int[, -c(2:6)]
        names(q1int)[names(q1int) == "lhs"] <- "Variable"

        q1var <- q1 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        q1var <- q1var[, -c(1, 2, 4:6)]
        names(q1var)[names(q1var) == "rhs"] <- "Variable"

        q1cov <- q1 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        q1cov <- q1cov[, -c(2, 4:6)]
        names(q1cov)[names(q1cov) == "lhs"] <- "Variable 1"
        names(q1cov)[names(q1cov) == "rhs"] <- "Variable 2"

        output$Loadm2g1 <- q1load
        output$Intm2g1 <- q1int

        output$Varm2g1 <- q1var

        output$Corrm2g1 <- q1cov


        # group 2
        q2load <- q2 %>%
            dplyr::filter(op == "=~")
        q2load <- q2load[, -c(2, 4:6)]
        names(q2load)[names(q2load) == "lhs"] <- "Factor"
        names(q2load)[names(q2load) == "rhs"] <- "Item"

        q2int <- q2 %>%
            dplyr::filter(op == "~1")
        q2int <- q2int[, -c(2:6)]
        names(q2int)[names(q2int) == "lhs"] <- "Variable"

        q2var <- q2 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        q2var <- q2var[, -c(1, 2, 4:6)]
        names(q2var)[names(q2var) == "rhs"] <- "Variable"

        q2cov <- q2 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        q2cov <- q2cov[, -c(2, 4:6)]
        names(q2cov)[names(q2cov) == "lhs"] <- "Variable 1"
        names(q2cov)[names(q2cov) == "rhs"] <- "Variable 2"

        output$Loadm2g2 <- q2load

        output$Intm2g2 <- q2int

        output$Varm2g2 <- q2var

        output$Corrm2g2 <- q2cov


        if (input$est == "ML") {
            indm2 <- t(as.data.frame(lavaan::fitmeasures(fit2, c("rmsea", "srmr", "cfi", "mfi", "chisq", "df", "pvalue"))))
            colnames(indm2) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
            SBdif1 <- lavaan::anova(a1, fit2)
        }

        if (input$est == "MLR") {
            indm2 <- t(as.data.frame(lavaan::fitmeasures(fit2, c("rmsea.scaled", "srmr", "cfi.scaled", "mfi", "chisq.scaled", "df.scaled", "pvalue.scaled"))))
            # lavaan doesn't have scaled version of mfi, so calculate it manually
            indm2[4] <- exp(-.5 * ((indm2[5] - indm2[6]) / (n0 - 1)))
            # calcuate satorra-bentler delta X2
            SBdif1 <- lavaan::anova(a1, fit2)
        }

        colnames(indm2) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
        rownames(indm2) <- ("Metric")
        Ind <- rbind(indm1, indm2)

        Dif <- t(as.data.frame(indm2[, 1:4] - indm1[, 1:4]))
        DifD1_num <- SBdif1[2, 5] - (indm2[, 6] - indm1[, 6])
        DifD1_num <- ifelse(DifD1_num < 0, 0, DifD1_num)
        DifD1 <- sqrt(2) * sqrt(DifD1_num / ((indm2[, 6] - indm1[, 6]) * (n0 - 2)))

        Dif <- cbind(Dif, DifD1)

        del <- "\u0394"
        colnames(Dif) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        rownames(Dif) <- c("Metric")

        Dif <- t(Dif)

        output$Fit <- Ind
        output$Del <- Dif
    }

    if (input$Inv > 1) {
        fit3 <- lavaan::cfa(data = data, model = model, estimator = input$est, orthogonal = COR, std.lv = TRUE, group = input$Group, group.equal = c("loadings", "intercepts"))

        qq <- lavaan::parameterEstimates(fit3)

        names(qq)[names(qq) == "est"] <- "Estimate"
        names(qq)[names(qq) == "se"] <- "SE"
        names(qq)[names(qq) == "z"] <- "Z"
        names(qq)[names(qq) == "pvalue"] <- "p"
        names(qq)[names(qq) == "ci.lower"] <- "95% CI Lower Limit"
        names(qq)[names(qq) == "ci.upper"] <- "95% CI Upper Limit"

        qq1 <- qq %>% dplyr::filter(group == 1)
        qq2 <- qq %>% dplyr::filter(group == 2)

        qq1load <- qq1 %>%
            dplyr::filter(op == "=~")
        qq1load <- qq1load[, -c(2, 4:6)]
        names(qq1load)[names(qq1load) == "lhs"] <- "Factor"
        names(qq1load)[names(qq1load) == "rhs"] <- "Item"

        qq1int <- qq1 %>%
            dplyr::filter(op == "~1")
        qq1int <- qq1int[, -c(2:6)]
        names(qq1int)[names(qq1int) == "lhs"] <- "Variable"

        qq1var <- qq1 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        qq1var <- qq1var[, -c(1, 2, 4:6)]
        names(qq1var)[names(qq1var) == "rhs"] <- "Variable"

        qq1cov <- qq1 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        qq1cov <- qq1cov[, -c(2, 4:6)]
        names(qq1cov)[names(qq1cov) == "lhs"] <- "Variable 1"
        names(qq1cov)[names(qq1cov) == "rhs"] <- "Variable 2"

        output$Loadm3g1 <-
            qq1load

        output$Intm3g1 <-
            qq1int

        output$Varm3g1 <-
            qq1var

        output$Corrm3g1 <- qq1cov


        # group 2
        qq2load <- qq2 %>%
            dplyr::filter(op == "=~")
        qq2load <- qq2load[, -c(2, 4:6)]
        names(qq2load)[names(qq2load) == "lhs"] <- "Factor"
        names(qq2load)[names(qq2load) == "rhs"] <- "Item"

        qq2int <- qq2 %>%
            dplyr::filter(op == "~1")
        qq2int <- qq2int[, -c(2:6)]
        names(qq2int)[names(qq2int) == "lhs"] <- "Variable"

        qq2var <- qq2 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        qq2var <- qq2var[, -c(1, 2, 4:6)]
        names(qq2var)[names(qq2var) == "rhs"] <- "Variable"

        qq2cov <- qq2 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        qq2cov <- qq2cov[, -c(2, 4:6)]
        names(qq2cov)[names(qq2cov) == "lhs"] <- "Variable 1"
        names(qq2cov)[names(qq2cov) == "rhs"] <- "Variable 2"

        output$Loadm3g2 <- qq2load
        output$Intm3g2 <- qq2int
        output$Varm3g2 <- qq2var
        output$Corrm3g2 <- qq2cov


        if (input$est == "ML") {
            indm3 <- t(as.data.frame(lavaan::fitmeasures(fit3, c("rmsea", "srmr", "cfi", "mfi", "chisq", "df", "pvalue"))))
            colnames(indm3) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
            SBdif2 <- lavaan::anova(fit2, fit3)
        }

        if (input$est == "MLR") {
            indm3 <- t(as.data.frame(lavaan::fitmeasures(fit3, c("rmsea.scaled", "srmr", "cfi.scaled", "mfi", "chisq.scaled", "df.scaled", "pvalue.scaled"))))
            # lavaan doesn't have scaled version of mfi, so calculate it manually
            indm3[4] <- exp(-.5 * ((indm3[5] - indm3[6]) / (n0 - 1)))
            SBdif2 <- lavaan::anova(fit2, fit3)
        }
        colnames(indm3) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
        rownames(indm3) <- ("Scalar")
        Ind <- rbind(indm1, indm2, indm3)

        Dif1 <- t(as.data.frame(indm2[, 1:4] - indm1[, 1:4]))
        DifD1_num <- SBdif1[2, 5] - (indm2[, 6] - indm1[, 6])
        DifD1_num <- ifelse(DifD1_num < 0, 0, DifD1_num)
        DifD1 <- sqrt(2) * sqrt(DifD1_num / ((indm2[, 6] - indm1[, 6]) * (n0 - 2)))
        Dif1 <- cbind(Dif1, DifD1)

        Dif2 <- t(as.data.frame(indm3[, 1:4] - indm2[, 1:4]))
        DifD2_num <- SBdif2[2, 5] - (indm3[, 6] - indm2[, 6])
        DifD2_num <- ifelse(DifD2_num < 0, 0, DifD2_num)
        DifD2 <- sqrt(2) * sqrt(DifD2_num / ((indm3[, 6] - indm2[, 6]) * (n0 - 2)))
        Dif2 <- cbind(Dif2, DifD2)

        Dif <- rbind(Dif1, Dif2)
        del <- "\u0394"
        colnames(Dif) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        rownames(Dif) <- c("Metric", "Scalar")

        Dif <- t(Dif)
        output$Fit <- Ind
        output$Del <- Dif
    }

    if (input$Inv > 2) {
        fit4 <- lavaan::cfa(data = data, model = model, estimator = input$est, orthogonal = COR, std.lv = TRUE, group = input$Group, group.equal = c("loadings", "intercepts", "residuals"))

        qqq <- lavaan::parameterEstimates(fit4)

        names(qqq)[names(qqq) == "est"] <- "Estimate"
        names(qqq)[names(qqq) == "se"] <- "SE"
        names(qqq)[names(qqq) == "z"] <- "Z"
        names(qqq)[names(qqq) == "pvalue"] <- "p"
        names(qqq)[names(qqq) == "ci.lower"] <- "95% CI Lower Limit"
        names(qqq)[names(qqq) == "ci.upper"] <- "95% CI Upper Limit"

        qqq1 <- qqq %>% dplyr::filter(group == 1)
        qqq2 <- qqq %>% dplyr::filter(group == 2)

        qqq1load <- qqq1 %>%
            dplyr::filter(op == "=~")
        qqq1load <- qqq1load[, -c(2, 4:6)]
        names(qqq1load)[names(qqq1load) == "lhs"] <- "Factor"
        names(qqq1load)[names(qqq1load) == "rhs"] <- "Item"

        qqq1int <- qqq1 %>%
            dplyr::filter(op == "~1")
        qqq1int <- qqq1int[, -c(2:6)]
        names(qqq1int)[names(qqq1int) == "lhs"] <- "Variable"

        qqq1var <- qqq1 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        qqq1var <- qqq1var[, -c(1, 2, 4:6)]
        names(qqq1var)[names(qqq1var) == "rhs"] <- "Variable"

        qqq1cov <- qqq1 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        qqq1cov <- qqq1cov[, -c(2, 4:6)]
        names(qqq1cov)[names(qqq1cov) == "lhs"] <- "Variable 1"
        names(qqq1cov)[names(qqq1cov) == "rhs"] <- "Variable 2"

        output$Loadm4g1 <- qqq1load

        output$Intm4g1 <- qqq1int

        output$Varm4g1 <- qqq1var

        output$Corrm4g1 <- qqq1cov


        # group 2
        qqq2load <- qqq2 %>%
            dplyr::filter(op == "=~")
        qqq2load <- qqq2load[, -c(2, 4:6)]
        names(qqq2load)[names(qqq2load) == "lhs"] <- "Factor"
        names(qqq2load)[names(qqq2load) == "rhs"] <- "Item"

        qqq2int <- qqq2 %>%
            dplyr::filter(op == "~1")
        qqq2int <- qqq2int[, -c(2:6)]
        names(qqq2int)[names(qqq2int) == "lhs"] <- "Variable"

        qqq2var <- qqq2 %>%
            dplyr::filter(lhs == rhs & op == "~~")
        qqq2var <- qqq2var[, -c(1, 2, 4:6)]
        names(qqq2var)[names(qqq2var) == "rhs"] <- "Variable"

        qqq2cov <- qqq2 %>%
            dplyr::filter(lhs != rhs & op == "~~")
        qqq2cov <- qqq2cov[, -c(2, 4:6)]
        names(qqq2cov)[names(qqq2cov) == "lhs"] <- "Variable 1"
        names(qqq2cov)[names(qqq2cov) == "rhs"] <- "Variable 2"

        output$Loadm4g2 <- qqq2load

        output$Intm4g2 <- qqq2int

        output$Varm4g2 <- qqq2var

        output$Corrm4g2 <- qqq2cov


        if (input$est == "ML") {
            indm4 <- t(as.data.frame(lavaan::fitmeasures(fit4, c("rmsea", "srmr", "cfi", "mfi", "chisq", "df", "pvalue"))))
            colnames(indm2) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
            SBdif3 <- lavaan::anova(fit3, fit4)
        }

        if (input$est == "MLR") {
            indm4 <- t(as.data.frame(lavaan::fitmeasures(fit4, c("rmsea.scaled", "srmr", "cfi.scaled", "mfi", "chisq.scaled", "df.scaled", "pvalue.scaled"))))
            # lavaan doesn't have scaled version of mfi, so calculate it manually
            indm4[4] <- exp(-.5 * ((indm4[5] - indm4[6]) / (n0 - 1)))
            SBdif3 <- lavaan::anova(fit3, fit4)
        }

        colnames(indm4) <- c("RMSEA", "SRMR", "CFI", "McDonald Non-Centrality", "Chi-Square", "df", "p-value")
        rownames(indm4) <- ("Strict")
        Ind <- rbind(indm1, indm2, indm3, indm4)

        Dif1 <- t(as.data.frame(indm2[, 1:4] - indm1[, 1:4]))
        DifD1_num <- SBdif1[2, 5] - (indm2[, 6] - indm1[, 6])
        DifD1_num <- ifelse(DifD1_num < 0, 0, DifD1_num)
        DifD1 <- sqrt(2) * sqrt(DifD1_num / ((indm2[, 6] - indm1[, 6]) * (n0 - 2)))
        Dif1 <- cbind(Dif1, DifD1)

        Dif2 <- t(as.data.frame(indm3[, 1:4] - indm2[, 1:4]))
        DifD2_num <- SBdif2[2, 5] - (indm3[, 6] - indm2[, 6])
        DifD2_num <- ifelse(DifD2_num < 0, 0, DifD2_num)
        DifD2 <- sqrt(2) * sqrt(DifD2_num / ((indm3[, 6] - indm2[, 6]) * (n0 - 2)))
        Dif2 <- cbind(Dif2, DifD2)

        Dif3 <- t(as.data.frame(indm4[, 1:4] - indm3[, 1:4]))
        DifD3_num <- SBdif3[2, 5] - (indm4[, 6] - indm3[, 6])
        DifD3_num <- ifelse(DifD3_num < 0, 0, DifD3_num)
        DifD3 <- sqrt(2) * sqrt(DifD3_num / ((indm4[, 6] - indm3[, 6]) * (n0 - 2)))
        Dif3 <- cbind(Dif3, DifD3)

        Dif <- rbind(Dif1, Dif2, Dif3)
        del <- "\u0394"
        colnames(Dif) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        rownames(Dif) <- c("Metric", "Scalar", "Strict")

        Dif <- t(Dif)
        output$Fit <- Ind
        output$Del <- Dif
    }


    sig.star <- function(lav = NULL, dat = NULL) {
        ss_mod_load <- suppressMessages(lav %>%
            dplyr::filter(lhs != rhs) %>%
            dplyr::group_by(lhs, op) %>%
            dplyr::filter(op != "~1") %>%
            dplyr::filter(op != "|") %>%
            dplyr::select(lhs, op, rhs, est) %>%
            # dplyr::mutate(est=round(est,digits=8)) %>%
            dplyr::reframe(rhs = paste(est, "*", rhs, collapse = " + ")) %>%
            dplyr::arrange(desc(op)) %>%
            tidyr::unite("mod", lhs, op, rhs, sep = "") %>%
            dplyr::pull(mod))

        ss_mod_res <- suppressMessages(lav %>%
            dplyr::filter(lhs == rhs) %>%
            dplyr::select(lhs, op, rhs, est) %>%
            # dplyr::mutate(est=round(est,digits=8)) %>%
            dplyr::reframe(rhs = paste(lhs, "~~", est, "*", rhs)) %>%
            dplyr::pull())

        # Collapse into one string because my other functions expect that
        mod <- base::paste(c(ss_mod_load, ss_mod_res), sep = "", collapse = "\n")

        # gdat<-dplyr::filter(dat,dat[,noquote(paste0(input$Group))]==grp)
        b <- lavaan::cfa(data = dat, model = mod)

        return(b@implied$cov[[1]])
    }

    sig.star0 <- sig.star(lav0, dat = data)


    ## helper functions for simulation

    duplication.matrix <- function(n = 1) {
        if ((n < 1) | (round(n) != n)) {
            stop("n must be a positive integer")
        }
        d <- matrix(0, n * n, n * (n + 1) / 2)
        count <- 0
        for (j in 1:n) {
            d[(j - 1) * n + j, count + j] <- 1
            if (j < n) {
                for (i in (j + 1):n) {
                    d[(j - 1) * n + i, count + i] <- 1
                    d[(i - 1) * n + j, count + i] <- 1
                }
            }
            count <- count + n - j
        }
        return(d)
    }

    vech <- function(x) {
        if (!is.square.matrix(x)) {
            stop("argument x is not a square numeric matrix")
        }
        return(t(t(x[!upper.tri(x)])))
    }

    is.square.matrix <- function(x) {
        if (!is.matrix(x)) {
            stop("argument x is not a matrix")
        }
        return(nrow(x) == ncol(x))
    }

    matrix.trace <- function(x) {
        if (!is.square.matrix(x)) {
            stop("argument x is not a square matrix")
        }
        return(sum(diag(x)))
    }

    # lav = parameter table, n = sample size in group G, imp= implied covariance matrix in Group G

    datgencov <- function(lav = NULL, n = NULL, imp = NULL, dat = NULL, df = NULL) {
        # preliminary definitions
        Sigma.gamma0 <- imp # model-implied from fitted model
        p <- nrow(Sigma.gamma0) # number of manifest variables
        q <- lav_partable_npar(lav) # number of model parameters
        p.star <- p * (p + 1) / 2 # number of non-duplicated entries

        discrep <- chi0 / n # desired disscrepancy is df / N
        delta <- discrep

        D.mat <- duplication.matrix(p)
        D <- t(D.mat) %*% D.mat
        Sigma.gamma0 <- sig.star(lav, dat)
        W <- Sigma.gamma0
        W.inv <- solve(W)

        # matrix manipulation functions

        # setup tracking matrices
        h <- 1e-8
        Sigma.deriv <- array(NA, c(p, p, q))
        B <- matrix(NA, p.star, q)

        # loop through parameter table and perturb each estimate
        for (i in 1:q) {
            lav1 <- lav
            lav1[i, 14] <- lav1[i, 14] + h # should this be ordered or does it not matter because the parameters are fixed anyway?
            # [order(-lav$free),]
            Sigma.gamma <- sig.star(lav1, dat)
            Sigma.deriv[, , i] <- (Sigma.gamma - Sigma.gamma0) * (1 / h)
            B[, i] <- (-1) * D %*% vech(W.inv %*% Sigma.deriv[, , i] %*% W.inv)
        }

        # randomly draw errors
        set.seed(101492)
        y <- matrix(rnorm(p.star), p.star, 1)
        B.qr <- qr(B)
        e.tilt <- qr.resid(B.qr, y)

        E1 <- matrix(0, p, p)
        index <- 1
        for (i2 in 1:p) {
            for (i1 in i2:p) {
                E1[i1, i2] <- e.tilt[index, 1]
                index <- index + 1
            }
        }

        E2 <- matrix(0, p, p)
        index <- 1
        for (i1 in 1:p) {
            for (i2 in i1:p) {
                E2[i1, i2] <- e.tilt[index, 1]
                index <- index + 1
            }
        }

        E.tilt <- E1 + E2 - diag(diag(E1))

        # rescale errors to have magnitude that will give desired result
        G <- W.inv %*% E.tilt
        get.kappa <- function(kappa, G, I, delta) {
            target <- abs(kappa * matrix.trace(G) - log(det(I + kappa * G)) - delta)
            return(target)
        }

        kappa0 <- sqrt(2 * delta / matrix.trace(G %*% G))
        I <- diag(p)
        res.kappa <- suppressWarnings(nlm(get.kappa, kappa0, G = G, I = I, delta = delta))
        kappa <- res.kappa$estimate
        iter <- res.kappa$iterations

        kappa <- as.numeric(kappa)
        E <- kappa * E.tilt

        # return matrix with perturbations that make chi-square = df
        Sigma.star <- Sigma.gamma0 + E

        return(Sigma.star)
    }

    if (input$Scale == "N") {
        imp <- imp0
    }

    if (input$Scale == "L") {
        varname <- c(unlist(a@Data@ov.names[[1]]))
        p <- list()
        for (i in 1:length(varname)) {
            xx <- as.data.frame(table(data[varname[1]])) %>% mutate(cum_sum = cumsum(Freq))
            xx1 <- xx %>% mutate(cum_prop = cum_sum / cum_sum[nrow(xx)])
            xx2 <- t(xx1[, 4])
            xx3 <- xx2[, 1:ncol(xx2) - 1]
            xx4 <- c(xx3)
            p[[i]] <- xx4
        }
        sd <- diag(imp0)
        R <- cov2cor(imp0)
        RR <- GenOrd::ordcont(marginal = p, Sigma = R)
        Rd <- RR$SigmaC
        imp <- cor2cov(Rd, sqrt(sd))
    }

    # if(input$Scale=="N"){
    #    imp=imp0
    #  }
    #
    #  if(input$Scale=="L"){
    #    sd<-diag(imp0)
    #    R<-cov2cor(imp0)
    #    p<-list()
    # #
    # #   #########################################################################
    # #   #How to select columns in the data that correspond to items in the model#
    # #   #########################################################################
    #     varname<-c(unlist(a@Data@ov.names[[1]]))
    #
    #   for (i in 1:length(varname)){
    #      xx<-as.data.frame(table(data()[,varname[i]]))%>% mutate(cum_sum = cumsum(Freq))
    #      xx1<-xx %>% mutate(cum_prop=cum_sum/cum_sum[nrow(xx)])
    #      xx2<-t(xx1[,4])
    #      xx3<-xx2[,1:ncol(xx2)-1]
    #      xx4<-c(xx3)
    #
    #      p[[i]]<-xx4
    #    }
    # #   ###################################################
    #
    #    a<-GenOrd::ordcont(marginal=p,Sigma=R)
    #    Rd<-a$SigmaC
    #    imp<-cor2cov(Rd,sqrt(sd))
    #    imp=imp0
    #  }

    Sigma.star0 <- datgencov(lav = lav0, n = n0, imp = imp, dat = data, df = df0)

    ##########################################################
    # add feature to get proportion of sample in each group  #
    # Use this to simulate data with correct balance         #
    ##########################################################

    ### Also deal with:
    # Check factor correlation: auto.cov.lv.x VS orthogonal
    # reproduce missing data
    # add option for residual covariance invariance?
    # Big one is to adapt for Likert


    true_fit_MI <- function(model, reps, n0, Scale = input$Scale) {
        # Can make this faster by only doing it once
        # Would need to change table. Not sure what would happen to plot.
        # Already did this

        # Number of reps
        r <- reps

        # Set Seed
        set.seed(8675)

        if (Scale == "N") {
            datax <- as.data.frame(MASS::mvrnorm(n = n0 * r, mu = rep(0, nrow(Sigma.star0)), Sigma.star0))
            colnames(datax) <- c(unlist(a@Data@ov.names[[1]]))
        }

        if (Scale == "L") {
            SigR <- cov2cor(Sigma.star0)
            datax <- as.data.frame(GenOrd::ordsample(n = n0 * r, Sigma = SigR, marginal = p))
            colnames(datax) <- c(unlist(a@Data@ov.names[[1]]))
        }

        prop <- table(data[, input$Group])[1] / length(data[, input$Group])

        g <- rbinom(n0 * r, 1, prop)
        datax1 <- cbind(datax, g)

        rep <- base::rep(1:r, n0)
        data_true <- base::cbind(datax1, rep)

        # Group and list
        true_data <- data_true %>%
            dplyr::group_by(rep) %>%
            tidyr::nest() %>%
            base::as.list()

        if (input$est == "ML") {
            ind <- c("srmr", "rmsea", "cfi", "mfi", "chisq", "df")
        }

        if (input$est == "MLR") {
            ind <- c("srmr", "rmsea.scaled", "cfi.scaled", "mfi", "chisq.scaled", "df.scaled")
        }

        # Run 500 cfa

        ## configural Model
        m1 <- purrr::map(true_data$data, function(x) {
            lavaan::cfa(
                model = model, group = "g", estimator = input$est, orthogonal = COR,
                data = x,
                std.lv = TRUE,
                check.gradient = FALSE,
                check.post = FALSE,
                check.vcov = FALSE,
                control = list(rel.tol = .001)
            )
        })
        m1_fit <- purrr::map_dfr(m1, ~ lavaan::fitMeasures(., ind)) %>%
            # `colnames<-`(c("SRMR_M1","RMSEA_M1","CFI_M1", "McD_M1","Chisq_M1"))
            `colnames<-`(c("SRMR_M1", "RMSEA_M1", "CFI_M1", "McD_M1", "Chisq_M1", "df_M1"))

        # mean(m1_fit$Chisq_M1)

        summary(m1[[1]])
        # metric models
        m2 <- purrr::map(true_data$data, function(x) {
            lavaan::cfa(
                model = model, group = "g", group.equal = "loadings", estimator = input$est, orthogonal = COR,
                data = x,
                std.lv = TRUE,
                check.gradient = FALSE,
                check.post = FALSE,
                check.vcov = FALSE,
                control = list(rel.tol = .001)
            )
        })
        m2_fit <- purrr::map_dfr(m2, ~ lavaan::fitMeasures(., ind)) %>%
            `colnames<-`(c("SRMR_M2", "RMSEA_M2", "CFI_M2", "McD_M2", "Chisq_M2", "df_M2"))
        dif1 <- purrr::pmap_dfr(list(m1, m2), lavaan::anova) %>%
            dplyr::select("Chisq diff") %>%
            na.omit() %>%
            `colnames<-`(c("x2dif_M2"))
        m2_fit <- cbind(m2_fit, dif1)

        if (input$Inv > 1) {
            # scalar models
            m3 <- purrr::map(true_data$data, function(x) {
                lavaan::cfa(
                    model = model, group = "g", group.equal = c("loadings", "intercepts"), estimator = input$est, orthogonal = COR,
                    data = x,
                    std.lv = TRUE,
                    check.gradient = FALSE,
                    check.post = FALSE,
                    check.vcov = FALSE,
                    control = list(rel.tol = .001)
                )
            })
            m3_fit <- purrr::map_dfr(m3, ~ lavaan::fitMeasures(., ind)) %>%
                `colnames<-`(c("SRMR_M3", "RMSEA_M3", "CFI_M3", "McD_M3", "Chisq_M3", "df_M3"))
            dif2 <- purrr::pmap_dfr(list(m2, m3), lavaan::anova) %>%
                dplyr::select("Chisq diff") %>%
                na.omit() %>%
                `colnames<-`(c("x2dif_M3"))
            m3_fit <- cbind(m3_fit, dif2)
        }

        if (input$Inv == 3) {
            # Strict models
            m4 <- purrr::map(true_data$data, function(x) {
                lavaan::cfa(
                    model = model, group = "g", group.equal = c("loadings", "intercepts", "residuals"), estimator = input$est, orthogonal = COR,
                    data = x,
                    std.lv = TRUE,
                    check.gradient = FALSE,
                    check.post = FALSE,
                    check.vcov = FALSE,
                    control = list(rel.tol = .001)
                )
            })
            m4_fit <- purrr::map_dfr(m4, ~ lavaan::fitMeasures(., ind)) %>%
                `colnames<-`(c("SRMR_M4", "RMSEA_M4", "CFI_M4", "McD_M4", "Chisq_M4", "df_M4"))
            dif3 <- purrr::pmap_dfr(list(m3, m4), lavaan::anova) %>%
                dplyr::select("Chisq diff") %>%
                na.omit() %>%
                `colnames<-`(c("x2dif_M4"))
            m4_fit <- cbind(m4_fit, dif3)
        }

        if (input$Inv == 1) {
            ## might have to change comparison

            m <- cbind(m1_fit, m2_fit) %>%
                mutate(n0 = n0) %>%
                mutate(RMSEA_metric = RMSEA_M2 - RMSEA_M1) %>%
                mutate(SRMR_metric = SRMR_M2 - SRMR_M1) %>%
                mutate(CFI_metric = CFI_M2 - CFI_M1) %>%
                mutate(McD_metric = McD_M2 - McD_M1) %>%
                mutate(RMSEAD_metric_num = x2dif_M2 - (df_M2 - df_M1)) %>%
                mutate(RMSEAD_metric_num = ifelse(RMSEAD_metric_num < 0, 0, RMSEAD_metric_num)) %>%
                mutate(RMSEAD_metric = sqrt(2) * sqrt(RMSEAD_metric_num / ((df_M2 - df_M1) * (n0 - 2))))


            res <- apply(m[, c(15:18, 20)], 2, quantile, probs = c(.01, .99))

            # signs are different, so save group by lower-is-better and higher-is-better
            pos <- res[2, c(1:2, 5)]
            neg <- res[1, 3:4]
            all_res <- as.data.frame(c(pos, neg))

            # one-column per fit index
            metric <- all_res[c(1, 2, 4, 5, 3), ]

            # table with more intuitive labels
            table <- cbind(metric)
            colnames(table) <- c("Metric")
            del <- "\u0394"
            rownames(table) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        }

        if (input$Inv == 2) {
            ## might have to change comparison
            ## Should it be 2 vs. 4 and not 3 vs 4?
            m <- cbind(m1_fit, m2_fit, m3_fit) %>%
                mutate(no = n0) %>%
                mutate(RMSEA_metric = RMSEA_M2 - RMSEA_M1) %>%
                mutate(RMSEA_Scalar = RMSEA_M3 - RMSEA_M2) %>%
                mutate(SRMR_metric = SRMR_M2 - SRMR_M1) %>%
                mutate(SRMR_Scalar = SRMR_M3 - SRMR_M2) %>%
                mutate(CFI_metric = CFI_M2 - CFI_M1) %>%
                mutate(CFI_Scalar = CFI_M3 - CFI_M2) %>%
                mutate(McD_metric = McD_M2 - McD_M1) %>%
                mutate(McD_Scalar = McD_M3 - McD_M2) %>%
                mutate(RMSEAD_metric_num = x2dif_M2 - (df_M2 - df_M1)) %>%
                mutate(RMSEAD_Scalar_num = x2dif_M3 - (df_M3 - df_M2)) %>%
                mutate(RMSEAD_metric_num = ifelse(RMSEAD_metric_num < 0, 0, RMSEAD_metric_num)) %>%
                mutate(RMSEAD_Scalar_num = ifelse(RMSEAD_Scalar_num < 0, 0, RMSEAD_Scalar_num)) %>%
                mutate(RMSEAD_metric = sqrt(2) * sqrt(RMSEAD_metric_num / ((df_M2 - df_M1) * (n0 - 2)))) %>%
                mutate(RMSEAD_Scalar = sqrt(2) * sqrt(RMSEAD_Scalar_num / ((df_M3 - df_M2) * (n0 - 2))))

            res <- apply(m[, c(22:29, 32:33)], 2, quantile, probs = c(.01, .99))

            # signs are different, so save group by lower-is-better and higher-is-better
            pos <- res[2, c(1:4, 9:10)]
            neg <- res[1, 5:8]
            all_res <- as.data.frame(c(pos, neg))

            # one-column per fit index
            metric <- all_res[c(1, 3, 7, 9, 5), ]
            scalar <- all_res[c(2, 4, 8, 10, 6), ]

            # table with more intuitive labels
            table <- cbind(metric, scalar)
            colnames(table) <- c("Metric", "Scalar")
            del <- "\u0394"
            rownames(table) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        }

        if (input$Inv == 3) {
            ## might have to change comparison
            ## Should it be 2 vs. 4 and not 3 vs 4?
            m <- cbind(m1_fit, m2_fit, m3_fit, m4_fit) %>%
                mutate(n0 = n0) %>%
                mutate(RMSEA_metric = RMSEA_M2 - RMSEA_M1) %>%
                mutate(RMSEA_Scalar = RMSEA_M3 - RMSEA_M2) %>%
                mutate(RMSEA_Strict = RMSEA_M4 - RMSEA_M3) %>%
                mutate(SRMR_metric = SRMR_M2 - SRMR_M1) %>%
                mutate(SRMR_Scalar = SRMR_M3 - SRMR_M2) %>%
                mutate(SRMR_Strict = SRMR_M4 - SRMR_M3) %>%
                mutate(CFI_metric = CFI_M2 - CFI_M1) %>%
                mutate(CFI_Scalar = CFI_M3 - CFI_M2) %>%
                mutate(CFI_Strict = CFI_M4 - CFI_M3) %>%
                mutate(McD_metric = McD_M2 - McD_M1) %>%
                mutate(McD_Scalar = McD_M3 - McD_M2) %>%
                mutate(McD_Strict = McD_M4 - McD_M3) %>%
                mutate(RMSEAD_metric_num = x2dif_M2 - (df_M2 - df_M1)) %>%
                mutate(RMSEAD_Scalar_num = x2dif_M3 - (df_M3 - df_M2)) %>%
                mutate(RMSEAD_Strict_num = x2dif_M4 - (df_M4 - df_M3)) %>%
                mutate(RMSEAD_metric_num = ifelse(RMSEAD_metric_num < 0, 0, RMSEAD_metric_num)) %>%
                mutate(RMSEAD_Scalar_num = ifelse(RMSEAD_Scalar_num < 0, 0, RMSEAD_Scalar_num)) %>%
                mutate(RMSEAD_Strict_num = ifelse(RMSEAD_Strict_num < 0, 0, RMSEAD_Strict_num)) %>%
                mutate(RMSEAD_metric = sqrt(2) * sqrt(RMSEAD_metric_num / ((df_M2 - df_M1) * (n0 - 2)))) %>%
                mutate(RMSEAD_Scalar = sqrt(2) * sqrt(RMSEAD_Scalar_num / ((df_M3 - df_M2) * (n0 - 2)))) %>%
                mutate(RMSEAD_Strict = sqrt(2) * sqrt(RMSEAD_Strict_num / ((df_M4 - df_M3) * (n0 - 2))))

            res <- apply(m[, c(29:40, 44:46)], 2, quantile, probs = c(.01, .99))

            # signs are different, so save group by lower-is-better and higher-is-better
            pos <- res[2, c(1:6, 13:15)]
            neg <- res[1, 7:12]
            all_res <- as.data.frame(c(pos, neg))

            # one-column per fit index
            metric <- all_res[c(1, 4, 10, 13, 7), ]
            scalar <- all_res[c(2, 5, 11, 14, 8), ]
            strict <- all_res[c(3, 6, 12, 15, 9), ]

            # table with more intuitive labels
            table <- cbind(metric, scalar, strict)
            colnames(table) <- c("Metric", "Scalar", "Strict")
            del <- "\u0394"
            rownames(table) <- c(paste0(del, "RMSEA"), paste0(del, "SRMR"), paste0(del, "CFI"), paste0(del, "McDonald Non-Centrality"), "RMSEA<sub>D</sub>")
        }

        # set.seed(NULL)

        return(round(table, 4))
    }
    Results <- true_fit_MI(model = model, reps = as.numeric(input$Reps), n0 = n0)


    output$DFI <- Results
    # I did not found any hint for decision.    #o utut$Dec <- Decision


    #  PD <- semPlot::semPaths(a, residuals = FALSE, intercepts = FALSE, thresholds = FALSE)
    #  output$PD <- semPlot::semPaths(a, residuals = FALSE, intercepts = FALSE, thresholds = FALSE)

    if (input$Inv == 1) {
        MI <- "Metric"
    }

    if (input$Inv == 2) {
        MI <- "Metric and Scalar"
    }

    if (input$Inv == 3) {
        MI <- "Metric, Scalar, and Strict"
    }

    parameter_tables <- list(
        metric_group_1 = list(
            loadings = if (exists("q1load")) q1load else NULL,
            intercepts = if (exists("q1int")) q1int else NULL,
            variances = if (exists("q1var")) q1var else NULL,
            covariances = if (exists("q1cov")) q1cov else NULL
        ),
        metric_group_2 = list(
            loadings = if (exists("q2load")) q2load else NULL,
            intercepts = if (exists("q2int")) q2int else NULL,
            variances = if (exists("q2var")) q2var else NULL,
            covariances = if (exists("q2cov")) q2cov else NULL
        ),
        scalar_group_1 = list(
            loadings = if (exists("qq1load")) qq1load else NULL,
            intercepts = if (exists("qq1int")) qq1int else NULL,
            variances = if (exists("qq1var")) qq1var else NULL,
            covariances = if (exists("qq1cov")) qq1cov else NULL
        ),
        scalar_group_2 = list(
            loadings = if (exists("qq2load")) qq2load else NULL,
            intercepts = if (exists("qq2int")) qq2int else NULL,
            variances = if (exists("qq2var")) qq2var else NULL,
            covariances = if (exists("qq2cov")) qq2cov else NULL
        ),
        strict_group_1 = list(
            loadings = if (exists("qqq1load")) qqq1load else NULL,
            intercepts = if (exists("qqq1int")) qqq1int else NULL,
            variances = if (exists("qqq1var")) qqq1var else NULL,
            covariances = if (exists("qqq1cov")) qqq1cov else NULL
        ),
        strict_group_2 = list(
            loadings = if (exists("qqq2load")) qqq2load else NULL,
            intercepts = if (exists("qqq2int")) qqq2int else NULL,
            variances = if (exists("qqq2var")) qqq2var else NULL,
            covariances = if (exists("qqq2cov")) qqq2cov else NULL
        )
    )
    parameter_tables <- lapply(parameter_tables, function(x) {
        x[!vapply(x, is.null, logical(1))]
    })
    parameter_tables <- parameter_tables[vapply(parameter_tables, length, integer(1)) > 0]


    .new_MgDynamic(
        input = input,
        model = model,
        fit = a,
        fit_indices = if (exists("Ind")) Ind else NULL,
        differences = if (exists("Dif")) Dif else NULL,
        dmacs = if (exists("D_new")) D_new else NULL,
        cutoffs = if (exists("Results")) Results else NULL,
        parameter_tables = parameter_tables,
        plot = if (exists("PD")) PD else NULL,
        outputs = output,
    )
}
