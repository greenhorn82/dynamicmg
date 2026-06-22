# dynamicMg

`dynamicMg` provides an R function interface for dynamic measurement invariance cutoffs for dependent groups.

 <!-- badges: start -->
  [![R-CMD-check](https://github.com/greenhorn82/dynamicMg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/greenhorn82/dynamicMg/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

## Installation

You can install the development version from GitHub with `devtools`:

```r
install.packages("devtools")
devtools::install_github("greenhorn82/dynamicmg")
```

## Background

This package is a port of the code from [dmcneish18/DG-DMI](https://github.com/dmcneish18/DG-DMI) into an R function, so that DG-DMI results can be integrated directly into Quarto or Markdown workflows.

The theoretical basis is dynamic measurement invariance cutoffs, especially McNeish’s work on DMI cutoffs for fit index differences: [Dynamic measurement invariance cutoffs for two-group fit index differences](https://psycnet.apa.org/record/2026-27788-001). DG-DMI extends this idea to dependent-group settings such as longitudinal or dyadic data.

## Example

The example below follows the structure used in the package tests: two factors, three indicators per factor, and a grouping variable named `F3`.

```r
library(dynamicmg)

set.seed(123)

n <- 300
Faktor1 <- rnorm(n)
Faktor2 <- rnorm(n)

df <- data.frame(
  Item1_F1 = 0.8 * Faktor1 + rnorm(n, 0, 0.5),
  Item2_F1 = 0.7 * Faktor1 + rnorm(n, 0, 0.5),
  Item3_F1 = 0.9 * Faktor1 + rnorm(n, 0, 0.5),
  Item1_F2 = 0.8 * Faktor2 + rnorm(n, 0, 0.5),
  Item2_F2 = 0.7 * Faktor2 + rnorm(n, 0, 0.5),
  Item3_F2 = 0.9 * Faktor2 + rnorm(n, 0, 0.5),
  F3 = rep(c(0, 1), each = n / 2)
)

result <- calcMG(
  data = df,
  loadings = list(
    c("Item1_F1", "Item2_F1", "Item3_F1"),
    c("Item1_F2", "Item2_F2", "Item3_F2")
  ),
  Group = "F3",
  Reps = 100,
  FacCor = FALSE,
  Inv = 3
)

result
```

For final analyses, use a larger number of replications, for example `Reps = 1000` or higher.
