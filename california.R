# Split Conformal Prediction - California Housing

library(tidyverse)
library(ranger)

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

set.seed(42)

ind <- sample(1:3, size = nrow(db), prob = c(0.80, 0.10, 0.10), replace = TRUE)

trn <- db[ind == 1, ]
cal <- db[ind == 2, ]
tst <- db[ind == 3, ]

alpha <- 0.1 # nominal miscoverage level

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
               method = "Linear Regression", color = "blue")

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
               method = "Random Forest", color = "green")

### Random Forest - Locally Weighted

y_hat_trn <- predict(rf, data = trn)$predictions

trn_res <- trn %>%
    mutate(delta = abs(median_house_value - y_hat_trn)) %>%
    select(-median_house_value)

rf_res <- ranger(delta ~ ., data = trn_res)

sig_hat_cal <- predict(rf_res, data = cal)$predictions

R <- abs(cal$median_house_value - y_hat_cal) / sig_hat_cal
r_hat <- sort(R)[ceiling((1 - alpha)*(nrow(cal) + 1))]

sig_hat_tst <- predict(rf_res, data = tst)$predictions

lower <- y_hat_tst - r_hat * sig_hat_tst
upper <- y_hat_tst + r_hat * sig_hat_tst

plot_intervals(tst$median_house_value, y_hat_tst, lower, upper,
               method = "Random Forest (Locally Weighted)", color = "red")
