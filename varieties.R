# The Varieties of Conformal Prediction - California Housing

options(java.parameters = "-Xmx8g")

library(tidyverse)
library(ranger)
library(bartMachine)

###

plot_intervals <- function(y, y_hat_tst, lower, upper, max_n = 100, method, color = "black") {
    coverage <- mean(lower <= tst$median_house_value & tst$median_house_value <= upper)
    avg_width <- mean(upper - lower)
 
    tibble(id = seq_along(y), y, y_hat_tst, lower, upper) %>%
        filter(id <= max_n) %>%
        ggplot(aes(x = id)) +
            geom_errorbar(aes(ymin = lower, ymax = upper), color = paste0("dark ", color)) +
            geom_point(aes(y = y_hat_tst), color = "blue", size = 1) +
            geom_point(aes(y = y), color = "red", size = 1) +
            scale_y_continuous(labels = scales::dollar_format(), limits = c(0, 650000)) +
            labs(x = "Test sample unit", y = "Price", title = method) +
            annotate("text", x = 60, y = 625000, hjust = 0, label = sprintf("Coverage = %.1f%%",  100 * coverage)) +
            annotate("text", x = 60, y = 550000, hjust = 0, label = paste0("Average width = $", format(avg_width, digits = 2, big.mark = ",", scientific = FALSE))) +
            theme_bw()
}

###

# Pace and Barry. "Sparse spatial autoregressions"
# Statistics and Probability Letters 33, no. 3 (1997): 291-297

db <- read.csv("california.csv", stringsAsFactors = TRUE)

glimpse(db)

# skimr::skim(db)

seed <- 42

set.seed(seed)

ind <- sample(1:3, size = nrow(db), prob = c(0.80, 0.10, 0.10), replace = TRUE)

trn <- db[ind == 1, ]
cal <- db[ind == 2, ]
tst <- db[ind == 3, ]

alpha <- 1 - 0.9

### Linear regression - Standard Score

lin_reg <- lm(median_house_value ~ ., data = trn)

y_hat_cal <- predict(lin_reg, newdata = cal)
R <- abs(cal$median_house_value - y_hat_cal)
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

y_hat_tst <- predict(lin_reg, newdata = tst)
sqrt(mean((tst$median_house_value - y_hat_tst)^2)) # test error
lower <- pmax(0, y_hat_tst - r_hat)
upper <- y_hat_tst + r_hat

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "SCP - Linear Regression", color = "blue")

### Random Forest - Standard Score

rf <- ranger(median_house_value ~  ., data = trn)

y_hat_cal <- predict(rf, data = cal)$predictions
R <- abs(cal$median_house_value - predict(rf, data = cal)$predictions)
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

y_hat_tst <- predict(rf, data = tst)$predictions
sqrt(mean((tst$median_house_value - y_hat_tst)^2))
lower <- pmax(0, y_hat_tst - r_hat)
upper <- y_hat_tst + r_hat

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "SCP - Random Forest", color = "green")

### Random Forest - Locally Weighted

y_hat_trn <- predict(rf, data = trn)$predictions

trn2 <- trn %>%
    mutate(delta = abs(median_house_value - y_hat_trn)) %>%
    select(-median_house_value)

rf2 <- ranger(delta ~ ., data = trn2)

rho_hat_cal <- predict(rf2, data = cal)$predictions

R <- abs(cal$median_house_value - y_hat_cal) / rho_hat_cal
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

rho_hat_tst <- predict(rf2, data = tst)$predictions

lower <- y_hat_tst - r_hat * rho_hat_tst
upper <- y_hat_tst + r_hat * rho_hat_tst

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "SCP Locally Weighted - Random Forest", color = "red")

### Conformalized Quantile Regression (CQR)

rf <- ranger(median_house_value ~ ., data = trn, quantreg = TRUE, num.trees = 10^3)

alpha_low <- alpha / 2
alpha_high <- 1 - alpha /2

q_hat_cal <- predict(rf, data = cal, type = "quantiles", quantiles = c(alpha_low, alpha_high))$predictions
E <- pmax(q_hat_cal[, 1] - cal$median_house_value, cal$median_house_value - q_hat_cal[, 2])
E_hat <- sort(E)[ceiling((1 - alpha)*(nrow(cal) + 1))]
q_hat_tst <- predict(rf, data = tst, type = "quantiles", quantiles = c(alpha_low, alpha_high))$predictions

lower <- q_hat_tst[, 1] - E_hat
upper <- q_hat_tst[, 2] + E_hat

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "Conformalized Quantile Regression", color = "orange")

### CQR v2

q_hat_cal <- predict(rf, data = cal, type = "quantiles", quantiles = c(0.25, 0.5, 0.75))$predictions
y_hat_cal <- q_hat_cal[, 2]
sig_hat_cal <- q_hat_cal[, 3] - q_hat_cal[, 1]

R <- abs(cal$median_house_value - y_hat_cal) / sig_hat_cal
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

q_hat_tst <- predict(rf, data = tst, type = "quantiles", quantiles = c(0.25, 0.5, 0.75))$predictions
y_hat_tst <- q_hat_tst[, 2]
sig_hat_tst <- q_hat_tst[, 3] - q_hat_tst[, 1]

lower <- y_hat_tst - r_hat * sig_hat_tst
upper <- y_hat_tst + r_hat * sig_hat_tst

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "CQR v2", color = "green")

### BART

X_trn <- trn[, -9]
y_trn <- trn$median_house_value
X_tst <- tst[, -9]
y_tst <- tst$median_house_value

bart <- bartMachine(X_trn, y_trn, seed = seed,
                    num_burn_in = 250,
                    num_iterations_after_burn_in = 1000)

post_tst <- bart_machine_get_posterior(bart, X_tst)

sqrt(mean((tst$median_house_value - post_tst$y_hat)^2))

# posterior predictive quantiles
PPQ_tst <- t(apply(post_tst$y_hat_posterior_samples, 1, function(.x) quantile(.x, prob = c(alpha / 2, 1 - alpha / 2))))

plot_intervals(tst$median_house_value, y_hat_tst, PPQ_tst[, 1], PPQ_tst[, 2],
               method = "BART", color = "blue")

### Conformalized BART

X_cal <- cal[, -9]
y_cal <- cal$median_house_value

post_cal <- bart_machine_get_posterior(bart, X_cal)
PPQ_cal <- t(apply(post_cal$y_hat_posterior_samples, 1, function(.x) quantile(.x, prob = c(alpha / 2, 1 - alpha / 2))))

R <- pmax(PPQ_cal[, 1] - cal$median_house_value, cal$median_house_value - PPQ_cal[, 2])
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

lower <- PPQ_tst[, 1] - r_hat
upper <- PPQ_tst[, 2] + r_hat

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "Conformalized BART", color = "red")

### Conformalized BART v2

y_hat_cal <- post_cal$y_hat
sig_hat_cal <- apply(post_cal$y_hat_posterior_samples, 1, sd)

R <- abs(cal$median_house_value - y_hat_cal) / sig_hat_cal
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

y_hat_tst <- post_tst$y_hat
sig_hat_tst <- apply(post_tst$y_hat_posterior_samples, 1, sd)

lower <- y_hat_tst - r_hat * sig_hat_tst
upper <- y_hat_tst + r_hat * sig_hat_tst

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "Conformalized BART v2", color = "green")

### Out-of-bag Conformal Prediction

rf1 <- ranger(median_house_value ~  ., data = trn, keep.inbag = TRUE, num.trees = 3e3)

y_hat_trn <- predict(rf1, data = trn, predict.all = TRUE)$predictions
out1 <- matrix(unlist(lapply(rf1$inbag.count, \(x) x == 0)), nrow = rf1$num.trees, byrow = TRUE)
mu_hat_trn <- sapply(1:nrow(trn), \(i) mean(y_hat_trn[i, out1[, i]]))

trn2 <- trn %>%
    mutate(delta = abs(median_house_value - mu_hat_trn)) %>%
    select(-median_house_value)

rf2 <- ranger(delta ~ ., data = trn2, keep.inbag = TRUE, num.trees = 3*10^3)

delta_hat_trn <- predict(rf2, data = trn2, predict.all = TRUE)$predictions
out2 <- matrix(unlist(lapply(rf2$inbag.count, \(x) x == 0)), nrow = rf2$num.trees, byrow = TRUE)
sig_hat_trn <- sapply(1:nrow(trn), \(i) mean(delta_hat_trn[i, out2[, i]]))

R <- abs(trn$median_house_value - mu_hat_trn) / sig_hat_trn
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(trn) + 1))]

mu_hat_tst <- predict(rf1, data = tst)$predictions
sig_hat_tst <- predict(rf2, data = tst)$predictions

lower <- mu_hat_tst - r_hat * sig_hat_tst
upper <- mu_hat_tst + r_hat * sig_hat_tst

plot_intervals(tst$median_house_value, mu_hat_tst, lower, upper,
               method = "Out-of-bag Conformal Prediction", color = "cyan")
