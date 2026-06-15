library(devtools)
set.seed(123)

# Anzahl der Personen / Fälle
n <- 300

# Zwei latente Faktoren simulieren
Faktor1 <- rnorm(n, mean = 0, sd = 1)
Faktor2 <- rnorm(n, mean = 0, sd = 1)

# Items für Faktor 1
Item1_F1 <- 0.8 * Faktor1 + rnorm(n, 0, 0.5)
Item2_F1 <- 0.7 * Faktor1 + rnorm(n, 0, 0.5)
Item3_F1 <- 0.9 * Faktor1 + rnorm(n, 0, 0.5)

# Items für Faktor 2
Item1_F2 <- 0.8 * Faktor2 + rnorm(n, 0, 0.5)
Item2_F2 <- 0.7 * Faktor2 + rnorm(n, 0, 0.5)
Item3_F2 <- 0.9 * Faktor2 + rnorm(n, 0, 0.5)
F3 <- rep(c(0, 1), each = n / 2)

# Dataframe erstellen
df <- data.frame(
    Item1_F1,
    Item2_F1,
    Item3_F1,
    Item1_F2,
    Item2_F2,
    Item3_F2,
    F3
)

calcMG(data = df, loadings = list(
    c("Item1_F1", "Item2_F1", "Item3_F1"),
    c("Item1_F2", "Item2_F2", "Item3_F2")),
    Group = "F3")
