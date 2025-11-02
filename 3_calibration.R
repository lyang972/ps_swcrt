# 3_calibration.R
# Bayesian calibration for rho* and duration-specific lower bounds for lambda.

stopifnot(exists("CONFIG_PATHS"), exists("CONFIG_STAN_AUX"))
stopifnot(file.exists(CONFIG_PATHS$PROCESSED_DATA))
load(CONFIG_PATHS$PROCESSED_DATA)

K_X <- stan_data$K_X; K_C <- stan_data$K_C; D_max <- stan_data$D_max
X_names <- colnames(X_matrix)

options(mc.cores = CONFIG_STAN_AUX$CORES)
rstan_options(auto_write = TRUE)

# Build lagged variables from treated periods
data_calib <- data_model %>%
  dplyr::group_by(i_idx) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(M_lag1 = dplyr::lag(M_std),
                Y_lag1 = dplyr::lag(Y),
                Z_lag1 = dplyr::lag(Z, default = 0)) %>%
  dplyr::ungroup()

data_treated_lags <- data_calib %>%
  dplyr::filter(Z == 1 & Duration > 0) %>%
  dplyr::filter(!is.na(M_lag1) & !is.na(M_std) & !is.na(Y_lag1) & !is.na(Y))
stopifnot(nrow(data_treated_lags) > 0, D_max >= 1)

i_idx_treated <- data_treated_lags$i_idx
i_idx_treated <- i_idx_treated[i_idx_treated <= nrow(X_matrix)]
dt <- data_treated_lags[data_treated_lags$i_idx %in% i_idx_treated, ]

df_calib <- data.frame(
  M_pre = dt$M_lag1, M_post = dt$M_std,
  Y_pre = dt$Y_lag1, Y_post = dt$Y,
  Province = dt$Province, Cluster_ID = dt$j_idx, Duration = dt$Duration
)
if (K_X > 0) {
  X_treated <- X_matrix[i_idx_treated, , drop = FALSE]; colnames(X_treated) <- X_names
  df_calib <- cbind(df_calib, X_treated)
}
X_calib <- if (K_X > 0) as.matrix(df_calib[, X_names, drop=FALSE]) else matrix(0, nrow(df_calib), 0)
if (K_C > 0) {
  C_calib_full <- model.matrix(~ Province, data = df_calib)
  C_calib <- C_calib_full[, -1, drop = FALSE]
} else C_calib <- matrix(0, nrow(df_calib), 0)

# ---- Rho* calibration (duration‑independent) ----
C_levels <- Provinces_List; N_C_levels <- length(C_levels)
Var_X_given_C_list <- list()
if (K_X > 0) {
  for (p in C_levels) {
    idx <- Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == p]
    Xp <- X_matrix[idx, , drop = FALSE]
    V <- if (nrow(Xp) > K_X + 5) tryCatch(cov(Xp), error = function(e) matrix(0, K_X, K_X)) else matrix(0, K_X, K_X)
    if (any(is.na(V))) V <- matrix(0, K_X, K_X)
    Var_X_given_C_list[[as.character(p)]] <- V
  }
} else {
  for (p in C_levels) Var_X_given_C_list[[as.character(p)]] <- matrix(0, 0, 0)
}
Var_X_given_C_array <- array(dim = c(N_C_levels, K_X, K_X))
if (N_C_levels > 0) for (i in 1:N_C_levels) Var_X_given_C_array[i, , ] <- Var_X_given_C_list[[as.character(C_levels[i])]]

stan_data_rho <- list(
  N = nrow(df_calib), D_max = D_max, K_X = K_X, K_C = K_C,
  M = cbind(df_calib$M_pre, df_calib$M_post),
  X = X_calib, C = C_calib, Duration_idx = df_calib$Duration,
  N_C_levels = N_C_levels, Var_X_given_C = Var_X_given_C_array
)
stan_file_rho <- CONFIG_PATHS$STAN_AUX_RHO; stopifnot(file.exists(stan_file_rho))
model_rho <- stan_model(file = stan_file_rho)
fit_rho <- sampling(model_rho, data = stan_data_rho,
                    iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS,
                    control = list(adapt_delta = CONFIG_STAN_AUX$ADAPT_DELTA, max_treedepth = CONFIG_STAN_AUX$MAX_TREEDEPTH))
posterior_rho_star_c <- rstan::extract(fit_rho)$rho_star_c
posterior_rho_star <- apply(posterior_rho_star_c, 1, mean, na.rm = TRUE)

# ---- Lambda lower bound calibration (duration‑specific) ----
stan_file_lambda <- CONFIG_PATHS$STAN_AUX_LAMBDA; stopifnot(file.exists(stan_file_lambda))
model_lambda <- stan_model(file = stan_file_lambda)

base_lambda <- list(N = nrow(df_calib), D_max = D_max, K_X = K_X, K_C = K_C,
                    X = X_calib, C = C_calib, Duration_idx = df_calib$Duration)

# Model for Y_post -> theta_2(d)
dat_Ypost <- base_lambda; dat_Ypost$Y <- df_calib$Y_post; dat_Ypost$M_post <- df_calib$M_post; dat_Ypost$M_pre <- df_calib$M_pre
fit_lambda_Ypost <- sampling(model_lambda, data = dat_Ypost,
                             iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)
posterior_theta_2_d <- rstan::extract(fit_lambda_Ypost)$Beta_M_pre

# Model for Y_pre -> zeta_1(d)
dat_Ypre <- base_lambda; dat_Ypre$Y <- df_calib$Y_pre; dat_Ypre$M_post <- df_calib$M_post; dat_Ypre$M_pre <- df_calib$M_pre
fit_lambda_Ypre <- sampling(model_lambda, data = dat_Ypre,
                            iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)
posterior_zeta_1_d <- rstan::extract(fit_lambda_Ypre)$Beta_M_post

calibration_results_bayesian <- list(
  Posteriors = list(
    rho_star = posterior_rho_star,
    lambda_z_lower = posterior_theta_2_d,
    lambda_minus_z_lower = posterior_zeta_1_d
  )
)
saveRDS(calibration_results_bayesian, file = CONFIG_PATHS$CALIBRATION_RESULTS)
