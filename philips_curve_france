# Clear environment
rm(list = ls())
cat("\014")
graphics.off()

# Packages
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(mFilter)
library(forecast)
library(gmm)
library(minpack.lm)
library(maxLik)


# 1) Load data ----
gdp_raw  <- read.xlsx("gdp.xlsx", sheet = "Table", colNames = FALSE)
dates    <- as.character(gdp_raw[1, 2:ncol(gdp_raw)])
gdp_real <- as.numeric(gdp_raw[2, 2:ncol(gdp_raw)])
gdp_nom  <- as.numeric(gdp_raw[3, 2:ncol(gdp_raw)])

ir_raw <- read.xlsx("interest_rates.xlsx", sheet = "Table", colNames = FALSE)
ir_10y <- as.numeric(ir_raw[2, 2:ncol(ir_raw)])
ir_3m  <- as.numeric(ir_raw[3, 2:ncol(ir_raw)])

ulc_raw <- read.xlsx("unit_labour_cost.xlsx", colNames = FALSE)
ulc     <- as.numeric(ulc_raw[2:nrow(ulc_raw), 2])

# Remove last NA column (184 observations: 1980Q1-2025Q4)
dates    <- dates[1:184]
gdp_real <- gdp_real[1:184]
gdp_nom  <- gdp_nom[1:184]
ir_10y   <- ir_10y[1:184]
ir_3m    <- ir_3m[1:184]
ulc      <- ulc[1:184]


# 2) Create ts objects and model variables ----
ts_data <- ts(data.frame(gdp_real, gdp_nom, ir_10y, ir_3m, ulc),
              start = c(1980, 1), frequency = 4)

# GDP deflator and quarterly inflation
gdp_def <- ts_data[, "gdp_nom"] / ts_data[, "gdp_real"] * 100
pi_def  <- diff(log(gdp_def))

# Output gap - HP filter on log real GDP
gdp_hp <- hpfilter(log(ts_data[, "gdp_real"]), freq = 1600, type = "lambda", drift = FALSE)
gap    <- gdp_hp$cycle

# Cyclical component of ULC - proxy for real marginal costs (see GG1999, eq. 13-14)
ulc_hp <- hpfilter(log(ts_data[, "ulc"]), freq = 1600, type = "lambda", drift = FALSE)
ulc_c  <- ulc_hp$cycle

# Interest rate spread and wage inflation
spread <- (ts_data[, "ir_10y"] - ts_data[, "ir_3m"]) / 100
pi_w   <- diff(log(ts_data[, "ulc"]))


# 3) Build ts object, restrict to 1980-2019 (pre-COVID distortions) ----
ts_y <- ts.union(pi_def, gap, ulc_c, spread, pi_w)
ts_y <- window(ts_y, start = c(1980, 1), end = c(2019, 4))


# 4) Visualise model series and cross-correlograms ----
p_series <- autoplot(ts_y, facets = TRUE, linewidth = 1) +
  labs(title = "France: NKPC model variables (1980:1-2019:4)",
       caption = "Source: OECD") +
  theme_bw()
print(p_series)
ggsave("01_model_series.png", plot = p_series, width = 10, height = 8, dpi = 300)

p_ccf_gap <- ggCcf(ts_y[, "gap"], ts_y[, "pi_def"], lag.max = 20) +
  labs(title = "Cross-correlation: gap(t) vs. pi_def(t+k)") +
  theme_bw()
print(p_ccf_gap)
ggsave("02_ccf_gap.png", plot = p_ccf_gap, width = 8, height = 5, dpi = 300)

p_ccf_ulc <- ggCcf(ts_y[, "ulc_c"], ts_y[, "pi_def"], lag.max = 20) +
  labs(title = "Cross-correlation: ulc_c(t) vs. pi_def(t+k)") +
  theme_bw()
print(p_ccf_ulc)
ggsave("03_ccf_ulc.png", plot = p_ccf_ulc, width = 8, height = 5, dpi = 300)


# 5) Prepare data and instruments for GMM ----
yy <- ts_y[, "pi_def"]

Z_1 <- stats::lag(ts_y[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -1)
Z_2 <- stats::lag(ts_y[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -2)
Z_3 <- stats::lag(ts_y[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -3)
Z_4 <- stats::lag(ts_y[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -4)
ZZ  <- ts.union(Z_1, Z_2, Z_3, Z_4)
k_Z <- dim(ZZ)[2]

# Data matrix for ULC models (k_X = 2)
XX_ulc    <- ts.union(ts_y[, "ulc_c"], stats::lag(ts_y[, "pi_def"], k = 1))
k_X       <- dim(XX_ulc)[2]
data_ulc  <- na.omit(ts.union(yy, XX_ulc, ZZ))
x_mat_ulc <- matrix(data_ulc, ncol = (1 + k_X + k_Z))

# Data matrix for gap models (k_X = 2)
XX_gap    <- ts.union(ts_y[, "gap"], stats::lag(ts_y[, "pi_def"], k = 1))
data_gap  <- na.omit(ts.union(yy, XX_gap, ZZ))
x_mat_gap <- matrix(data_gap, ncol = (1 + k_X + k_Z))

# Data matrix for hybrid models (k_X = 3)
XX_hyb    <- ts.union(ts_y[, "ulc_c"],
                      stats::lag(ts_y[, "pi_def"], k =  1),
                      stats::lag(ts_y[, "pi_def"], k = -1))
k_X_hyb   <- dim(XX_hyb)[2]
data_hyb  <- na.omit(ts.union(yy, XX_hyb, ZZ))
x_mat_hyb <- matrix(data_hyb, ncol = (1 + k_X_hyb + k_Z))


# 6) Moment condition functions ----

# Reduced-form NKPC: pi_t = lambda*x_t + beta*E_t{pi_{t+1}}
my_g1 <- function(my_par, x) {
  lambda <- my_par[1]
  beta   <- my_par[2]
  e      <- x[, 1] - lambda * x[, 2] - beta * x[, 3]
  Z      <- x[, (1 + k_X + 1):(1 + k_X + k_Z)]
  as.numeric(e) * Z
}

# Structural NKPC, normalisation (18) from GG1999
my_g2 <- function(my_par, x) {
  theta <- my_par[1]
  beta  <- my_par[2]
  e     <- theta * x[, 1] - (1 - theta) * (1 - beta * theta) * x[, 2] - theta * beta * x[, 3]
  Z     <- x[, (1 + k_X + 1):(1 + k_X + k_Z)]
  as.numeric(e) * Z
}

# Structural NKPC, normalisation (19) from GG1999
my_g3 <- function(my_par, x) {
  theta <- my_par[1]
  beta  <- my_par[2]
  e     <- x[, 1] - (1/theta) * (1 - theta) * (1 - beta * theta) * x[, 2] - beta * x[, 3]
  Z     <- x[, (1 + k_X + 1):(1 + k_X + k_Z)]
  as.numeric(e) * Z
}

# Hybrid NKPC, reduced form (eq. 24 GG1999)
my_g4 <- function(my_par, x) {
  lambda  <- my_par[1]
  gamma_f <- my_par[2]
  gamma_b <- my_par[3]
  e       <- x[, 1] - lambda * x[, 2] - gamma_f * x[, 3] - gamma_b * x[, 4]
  Z       <- x[, (1 + k_X_hyb + 1):(1 + k_X_hyb + k_Z)]
  as.numeric(e) * Z
}

# Hybrid NKPC, structural form (eq. 27+28 GG1999)
my_g5 <- function(my_par, x) {
  theta   <- my_par[1]
  beta    <- my_par[2]
  omega   <- my_par[3]
  phi     <- theta + omega * (1 - theta * (1 - beta))
  lambda  <- (1 - omega) * (1 - theta) * (1 - beta * theta) * phi^(-1)
  gamma_f <- beta * theta * phi^(-1)
  gamma_b <- omega * phi^(-1)
  e       <- x[, 1] - lambda * x[, 2] - gamma_f * x[, 3] - gamma_b * x[, 4]
  Z       <- x[, (1 + k_X_hyb + 1):(1 + k_X_hyb + k_Z)]
  as.numeric(e) * Z
}

# Restricted reduced-form NKPC: beta = 1, normalisation (19)
my_g3_b1 <- function(my_par, x) {
  theta <- my_par[1]
  beta  <- 1
  e     <- x[, 1] - (1/theta) * (1 - theta) * (1 - beta * theta) * x[, 2] - beta * x[, 3]
  Z     <- x[, (1 + k_X + 1):(1 + k_X + k_Z)]
  as.numeric(e) * Z
}

# Restricted hybrid NKPC: beta = 1 => gamma_f + gamma_b = 1
my_g4_b1 <- function(my_par, x) {
  lambda  <- my_par[1]
  gamma_f <- my_par[2]
  gamma_b <- 1 - gamma_f
  e       <- x[, 1] - lambda * x[, 2] - gamma_f * x[, 3] - gamma_b * x[, 4]
  Z       <- x[, (1 + k_X_hyb + 1):(1 + k_X_hyb + k_Z)]
  as.numeric(e) * Z
}


# 7) GMM estimates - reduced form ----

# Model 1: pi_t = lambda*ulc_t + beta*E_t{pi_{t+1}}
Model_1 <- gmm(my_g1, x = x_mat_ulc, t0 = c(0.02, 0.95), type = "twoStep", vcov = "HAC")
summary(Model_1)

# Model 2: pi_t = lambda*gap_t + beta*E_t{pi_{t+1}}
Model_2 <- gmm(my_g1, x = x_mat_gap, t0 = c(0.02, 0.95), type = "twoStep", vcov = "HAC")
summary(Model_2)


# 8) GMM estimates - structural form ----

# Model 3: normalisation (18)
Model_3 <- gmm(my_g2, x = x_mat_ulc, t0 = c(0.75, 0.95), type = "twoStep", vcov = "HAC")
summary(Model_3)

theta_3  <- coefficients(Model_3)[1]
beta_3   <- coefficients(Model_3)[2]
lambda_3 <- (1 - theta_3) * (1 - beta_3 * theta_3) / theta_3
cat("Model 3 | theta =", round(theta_3, 4), "| beta =", round(beta_3, 4),
    "| lambda =", round(lambda_3, 4), "| price fixation duration =", round(1/(1-theta_3), 2), "Q\n")

x_str   <- x_mat_ulc[, 1:(1 + k_X)]
resid_3 <- theta_3 * x_str[, 1] - (1 - theta_3) * (1 - beta_3 * theta_3) * x_str[, 2] - theta_3 * beta_3 * x_str[, 3]

p_acf_3 <- ggAcf(resid_3, lag.max = 20) +
  labs(title = "Residual correlogram - Model 3 (normalisation 18)") + theme_bw()
print(p_acf_3)
ggsave("04_acf_model3.png", plot = p_acf_3, width = 8, height = 5, dpi = 300)

# Model 4: normalisation (19)
Model_4 <- gmm(my_g3, x = x_mat_ulc, t0 = c(0.75, 0.95), type = "twoStep", vcov = "HAC")
summary(Model_4)

theta_4  <- coefficients(Model_4)[1]
beta_4   <- coefficients(Model_4)[2]
lambda_4 <- (1 - theta_4) * (1 - beta_4 * theta_4) / theta_4
cat("Model 4 | theta =", round(theta_4, 4), "| beta =", round(beta_4, 4),
    "| lambda =", round(lambda_4, 4), "| price fixation duration =", round(1/(1-theta_4), 2), "Q\n")

resid_4 <- theta_4 * x_str[, 1] - (1 - theta_4) * (1 - beta_4 * theta_4) * x_str[, 2] - theta_4 * beta_4 * x_str[, 3]

p_acf_4 <- ggAcf(resid_4, lag.max = 20) +
  labs(title = "Residual correlogram - Model 4 (normalisation 19)") + theme_bw()
print(p_acf_4)
ggsave("05_acf_model4.png", plot = p_acf_4, width = 8, height = 5, dpi = 300)


# 9) GMM estimates - hybrid NKPC ----

# Model 5: reduced-form hybrid
Model_5 <- gmm(my_g4, x = x_mat_hyb, t0 = c(0.02, 0.6, 0.3), type = "twoStep", vcov = "HAC")
summary(Model_5)

lambda_5  <- coefficients(Model_5)[1]
gamma_f_5 <- coefficients(Model_5)[2]
gamma_b_5 <- coefficients(Model_5)[3]
cat("Model 5 | lambda =", round(lambda_5, 4), "| gamma_f =", round(gamma_f_5, 4),
    "| gamma_b =", round(gamma_b_5, 4), "\n")

x_hyb   <- x_mat_hyb[, 1:(1 + k_X_hyb)]
resid_5 <- x_hyb[, 1] - lambda_5 * x_hyb[, 2] - gamma_f_5 * x_hyb[, 3] - gamma_b_5 * x_hyb[, 4]

p_acf_5 <- ggAcf(resid_5, lag.max = 20) +
  labs(title = "Residual correlogram - Model 5 (hybrid NKPC, reduced form)") + theme_bw()
print(p_acf_5)
ggsave("06_acf_model5.png", plot = p_acf_5, width = 8, height = 5, dpi = 300)

# Model 6: structural hybrid form
Model_6 <- gmm(my_g5, x = x_mat_hyb, t0 = c(0.75, 0.95, 0.3), type = "iterative", vcov = "HAC")
summary(Model_6)

theta_6   <- coefficients(Model_6)[1]
beta_6    <- coefficients(Model_6)[2]
omega_6   <- coefficients(Model_6)[3]
phi_6     <- theta_6 + omega_6 * (1 - theta_6 * (1 - beta_6))
lambda_6  <- (1 - omega_6) * (1 - theta_6) * (1 - beta_6 * theta_6) * phi_6^(-1)
gamma_f_6 <- beta_6 * theta_6 * phi_6^(-1)
gamma_b_6 <- omega_6 * phi_6^(-1)
cat("Model 6 | theta =", round(theta_6, 4), "| beta =", round(beta_6, 4),
    "| omega =", round(omega_6, 4), "\n")
cat("        | lambda =", round(lambda_6, 4), "| gamma_f =", round(gamma_f_6, 4),
    "| gamma_b =", round(gamma_b_6, 4), "\n")

resid_6 <- x_hyb[, 1] - lambda_6 * x_hyb[, 2] - gamma_f_6 * x_hyb[, 3] - gamma_b_6 * x_hyb[, 4]

p_acf_6 <- ggAcf(resid_6, lag.max = 20) +
  labs(title = "Residual correlogram - Model 6 (hybrid NKPC, structural form)") + theme_bw()
print(p_acf_6)
ggsave("07_acf_model6.png", plot = p_acf_6, width = 8, height = 5, dpi = 300)


# 10) Test for omitted dynamics - Ljung-Box test of residuals ----
cat("\n--- Ljung-Box tests ---\n")
Box.test(resid_3, lag = 8,  type = "Ljung-Box")
Box.test(resid_3, lag = 12, type = "Ljung-Box")
Box.test(resid_4, lag = 8,  type = "Ljung-Box")
Box.test(resid_4, lag = 12, type = "Ljung-Box")
Box.test(resid_5, lag = 8,  type = "Ljung-Box")
Box.test(resid_5, lag = 12, type = "Ljung-Box")
Box.test(resid_6, lag = 8,  type = "Ljung-Box")
Box.test(resid_6, lag = 12, type = "Ljung-Box")

# Robustness check: narrower sample 1980-2000 (pre-euro and pre-GFC)
ts_y_pre2000  <- window(ts_y, start = c(1980, 1), end = c(2000, 4))
yy_pre2000    <- ts_y_pre2000[, "pi_def"]
Z_1p <- stats::lag(ts_y_pre2000[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -1)
Z_2p <- stats::lag(ts_y_pre2000[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -2)
Z_3p <- stats::lag(ts_y_pre2000[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -3)
Z_4p <- stats::lag(ts_y_pre2000[, c("pi_def", "gap", "ulc_c", "spread", "pi_w")], k = -4)
ZZ_p          <- ts.union(Z_1p, Z_2p, Z_3p, Z_4p)
XX_pre2000    <- ts.union(ts_y_pre2000[, "ulc_c"], stats::lag(ts_y_pre2000[, "pi_def"], k = 1))
data_pre2000  <- na.omit(ts.union(yy_pre2000, XX_pre2000, ZZ_p))
x_mat_pre2000 <- matrix(data_pre2000, ncol = (1 + k_X + dim(ZZ_p)[2]))

Model_3_pre2000 <- gmm(my_g2, x = x_mat_pre2000, t0 = c(0.75, 0.95), type = "twoStep", vcov = "HAC")
summary(Model_3_pre2000)

theta_p  <- coefficients(Model_3_pre2000)[1]
beta_p   <- coefficients(Model_3_pre2000)[2]
lambda_p <- (1 - theta_p) * (1 - beta_p * theta_p) / theta_p
cat("Model 3 (1980-2000) | theta =", round(theta_p, 4), "| beta =", round(beta_p, 4),
    "| lambda =", round(lambda_p, 4), "| price fixation duration =", round(1/(1-theta_p), 2), "Q\n")

x_str_p  <- x_mat_pre2000[, 1:(1 + k_X)]
resid_p  <- theta_p * x_str_p[, 1] - (1 - theta_p) * (1 - beta_p * theta_p) * x_str_p[, 2] - theta_p * beta_p * x_str_p[, 3]
Box.test(resid_p, lag = 8,  type = "Ljung-Box")
Box.test(resid_p, lag = 12, type = "Ljung-Box")


# 11) Restricted models - beta = 1 (long-run vertical Phillips curve test) ----

# Structural NKPC with beta = 1
Model_4_b1 <- gmm(my_g3_b1, x = x_mat_ulc, t0 = c(0.75), type = "twoStep", vcov = "HAC")
summary(Model_4_b1)

theta_b1  <- coefficients(Model_4_b1)[1]
lambda_b1 <- (1 - theta_b1)^2 / theta_b1
cat("Model 4 (beta=1) | theta =", round(theta_b1, 4), "| lambda =", round(lambda_b1, 4),
    "| price fixation duration =", round(1/(1-theta_b1), 2), "Q\n")

# Hybrid NKPC with beta = 1 (gamma_f + gamma_b = 1)
Model_5_b1 <- gmm(my_g4_b1, x = x_mat_hyb, t0 = c(0.02, 0.6), type = "twoStep")
summary(Model_5_b1)

lambda_5b1  <- coefficients(Model_5_b1)[1]
gamma_f_5b1 <- coefficients(Model_5_b1)[2]
gamma_b_5b1 <- 1 - gamma_f_5b1
cat("Model 5 (beta=1) | lambda =", round(lambda_5b1, 4), "| gamma_f =", round(gamma_f_5b1, 4),
    "| gamma_b =", round(gamma_b_5b1, 4), "\n")


# 12) NLS estimation ----
nls_data <- as.data.frame(x_mat_ulc)
colnames(nls_data)[1:3] <- c("pi", "ulc", "pi_f")

Model_NLS <- nlsLM(
  pi ~ (1/theta) * (1 - theta) * (1 - beta * theta) * ulc + beta * pi_f,
  data  = nls_data,
  start = list(theta = 0.75, beta = 0.95)
)
summary(Model_NLS)

theta_nls  <- coef(Model_NLS)["theta"]
beta_nls   <- coef(Model_NLS)["beta"]
lambda_nls <- (1 - theta_nls) * (1 - beta_nls * theta_nls) / theta_nls
cat("Lambda (unrestricted NLS):", round(lambda_nls, 4), "\n")

logLik(Model_NLS)

fixation_duration_nls <- 1 / (1 - theta_nls)
cat("Price fixation duration (NLS unrestricted):", round(fixation_duration_nls, 2), "quarters\n")

# NLS restricted model (beta = 1)
Model_NLS_b1 <- nlsLM(
  pi ~ (1/theta) * (1 - theta) * (1 - theta) * ulc + pi_f,
  data  = nls_data,
  start = list(theta = 0.75)
)
summary(Model_NLS_b1)

theta_nls_b1  <- coef(Model_NLS_b1)["theta"]
lambda_nls_b1 <- (1 - theta_nls_b1)^2 / theta_nls_b1
cat("Lambda (restricted NLS):", round(lambda_nls_b1, 4), "\n")

logLik(Model_NLS_b1)

fixation_duration_nls_b1 <- 1 / (1 - theta_nls_b1)
cat("Price fixation duration (NLS restricted):", round(fixation_duration_nls_b1, 2), "quarters\n")

# LR test of long-run verticality for NLS
RSS_unres <- sum(residuals(Model_NLS)^2)
RSS_res   <- sum(residuals(Model_NLS_b1)^2)
LR_NLS    <- nrow(nls_data) * log(RSS_res / RSS_unres)
p_LR_NLS  <- pchisq(LR_NLS, df = 1, lower.tail = FALSE)
cat("LR test (NLS) | stat =", round(LR_NLS, 4), "| p =", round(p_LR_NLS, 4), "\n")


# 13) ML estimation - normal distribution ----
sig_start <- sd(residuals(Model_NLS))

loglik_norm <- function(my_par) {
  theta <- my_par[1]
  beta  <- my_par[2]
  sig   <- my_par[3]
  e     <- nls_data$pi - (1/theta) * (1 - theta) * (1 - beta * theta) * nls_data$ulc - beta * nls_data$pi_f
  -0.5 * log(2 * pi) - 0.5 * log(sig^2) - 0.5 * e^2 / sig^2
}

Model_ML_norm <- maxLik(loglik_norm,
                        start  = c(theta = 0.788, beta = 0.931, sig = sig_start),
                        method = "BHHH")
summary(Model_ML_norm)

loglik_ml <- logLik(Model_ML_norm)

theta_ml <- coef(Model_ML_norm)["theta"]
beta_ml  <- coef(Model_ML_norm)["beta"]

lambda_ml  <- (1 - theta_ml) * (1 - beta_ml * theta_ml) / theta_ml
fixace_ml  <- 1 / (1 - theta_ml)

cat("--- Unrestricted ML estimate ---\n")
cat("Log-Likelihood :", round(loglik_ml, 4), "\n")
cat("Lambda         :", round(lambda_ml, 4), "\n")
cat("Price fixation :", round(fixace_ml, 2), "quarters\n\n")

# ML restricted model (beta = 1)
loglik_norm_b1 <- function(my_par) {
  theta <- my_par[1]
  sig   <- my_par[2]
  e     <- nls_data$pi - (1/theta) * (1 - theta) * (1 - theta) * nls_data$ulc - nls_data$pi_f
  -0.5 * log(2 * pi) - 0.5 * log(sig^2) - 0.5 * e^2 / sig^2
}

Model_ML_norm_b1 <- maxLik(loglik_norm_b1,
                           start  = c(theta = 0.788, sig = sig_start),
                           method = "BHHH")
summary(Model_ML_norm_b1)

LR_norm   <- -2 * (logLik(Model_ML_norm_b1) - logLik(Model_ML_norm))
p_LR_norm <- pchisq(LR_norm, df = 1, lower.tail = FALSE)
cat("LR test (ML normal) | stat =", round(LR_norm, 4), "| p =", round(p_LR_norm, 4), "\n")

loglik_ml_b1 <- logLik(Model_ML_norm_b1)
theta_ml_b1  <- coef(Model_ML_norm_b1)["theta"]
lambda_ml_b1 <- (1 - theta_ml_b1)^2 / theta_ml_b1
fixace_ml_b1 <- 1 / (1 - theta_ml_b1)

cat("--- Restricted ML estimate (beta = 1) ---\n")
cat("Log-Likelihood :", round(loglik_ml_b1, 4), "\n")
cat("Lambda         :", round(lambda_ml_b1, 4), "\n")
cat("Price fixation :", round(fixace_ml_b1, 2), "quarters\n")


# 14) ML estimation - t-distribution ----
loglik_t <- function(my_par) {
  theta <- my_par[1]
  beta  <- my_par[2]
  sig   <- my_par[3]
  nu    <- my_par[4]
  e     <- nls_data$pi - (1/theta) * (1 - theta) * (1 - beta * theta) * nls_data$ulc - beta * nls_data$pi_f
  dt(e / sig, df = nu, log = TRUE) - log(sig)
}

Model_ML_t <- maxLik(loglik_t,
                     start  = c(theta = 0.788, beta = 0.931, sig = sig_start, nu = 5),
                     method = "BHHH")
summary(Model_ML_t)

LR_norm_vs_t   <- -2 * (logLik(Model_ML_norm) - logLik(Model_ML_t))
p_LR_norm_vs_t <- pchisq(LR_norm_vs_t, df = 1, lower.tail = FALSE)

cat("\n--- Normality test (Normal vs. t-distribution) ---\n")
cat("LR test | stat =", round(LR_norm_vs_t, 4), "| p =", round(p_LR_norm_vs_t, 4), "\n")

loglik_ml_t <- logLik(Model_ML_t)
theta_ml_t  <- coef(Model_ML_t)["theta"]
beta_ml_t   <- coef(Model_ML_t)["beta"]
lambda_ml_t <- (1 - theta_ml_t) * (1 - beta_ml_t * theta_ml_t) / theta_ml_t
fixace_ml_t <- 1 / (1 - theta_ml_t)

cat("--- Unrestricted ML estimate (t-distribution) ---\n")
cat("Log-Likelihood :", round(loglik_ml_t, 4), "\n")
cat("Lambda         :", round(lambda_ml_t, 4), "\n")
cat("Price fixation :", round(fixace_ml_t, 2), "quarters\n\n")

# ML t restricted model (beta = 1)
loglik_t_b1 <- function(my_par) {
  theta <- my_par[1]
  sig   <- my_par[2]
  nu    <- my_par[3]
  e     <- nls_data$pi - (1/theta) * (1 - theta) * (1 - theta) * nls_data$ulc - nls_data$pi_f
  dt(e / sig, df = nu, log = TRUE) - log(sig)
}

Model_ML_t_b1 <- maxLik(loglik_t_b1,
                        start  = c(theta = 0.808, sig = sig_start, nu = 5),
                        method = "BHHH")
summary(Model_ML_t_b1)

LR_t   <- -2 * (logLik(Model_ML_t_b1) - logLik(Model_ML_t))
p_LR_t <- pchisq(LR_t, df = 1, lower.tail = FALSE)
cat("LR test (ML t) | stat =", round(LR_t, 4), "| p =", round(p_LR_t, 4), "\n")

loglik_ml_t_b1 <- logLik(Model_ML_t_b1)
theta_ml_t_b1  <- coef(Model_ML_t_b1)["theta"]
lambda_ml_t_b1 <- (1 - theta_ml_t_b1)^2 / theta_ml_t_b1
fixace_ml_t_b1 <- 1 / (1 - theta_ml_t_b1)

cat("--- Restricted ML estimate (t-distribution, beta = 1) ---\n")
cat("Log-Likelihood :", round(loglik_ml_t_b1, 4), "\n")
cat("Lambda         :", round(lambda_ml_t_b1, 4), "\n")
cat("Price fixation :", round(fixace_ml_t_b1, 2), "quarters\n")
