# 3_calibration.R
# Bayesian calibration for duration-independent rho* and duration-dependent lambda bounds.

if (!"dplyr" %in% loadedNamespaces()) library(dplyr)
if (!"rstan" %in% loadedNamespaces()) library(rstan)
if (!"tidyr" %in% loadedNamespaces()) library(tidyr)
if (!"Matrix" %in% loadedNamespaces()) library(Matrix)

# Load config and processed data
load("Analysis_Config.RData")
load(CONFIG_PATHS$PROCESSED_DATA)

if (!exists("K_X")) K_X = stan_data$K_X
if (!exists("K_C")) K_C = stan_data$K_C
if (!exists("D_max")) D_max = stan_data$D_max
X_names = colnames(X_matrix)

options(mc.cores = CONFIG_STAN_AUX$CORES)
rstan_options(auto_write = TRUE)

# ---------- 1) Prepare data for calibration ----------
data_calibration = data_model %>%
  group_by(i_idx) %>%
  arrange(time) %>%
  mutate(
    M_lag1 = lag(M_std),
    Y_lag1 = lag(Y)
  ) %>%
  ungroup()

data_treated_lags = data_calibration %>%
  filter(Z == 1, Duration > 0) %>%
  filter(!is.na(M_lag1), !is.na(M_std), !is.na(Y_lag1), !is.na(Y))
stopifnot(nrow(data_treated_lags) > 0, D_max >= 1)

i_idx_treated = data_treated_lags$i_idx
i_idx_treated = i_idx_treated[i_idx_treated <= nrow(X_matrix)]
stopifnot(length(i_idx_treated) > 0)

X_treated = X_matrix[i_idx_treated, , drop = FALSE]
data_treated_lags_filtered = data_treated_lags[data_treated_lags$i_idx %in% i_idx_treated, ]

df_calib_input = data.frame(
  M_pre  = data_treated_lags_filtered$M_lag1,
  M_post = data_treated_lags_filtered$M_std,
  Y_pre  = data_treated_lags_filtered$Y_lag1,
  Y_post = data_treated_lags_filtered$Y,
  Province  = data_treated_lags_filtered$Province,
  Cluster_ID = data_treated_lags_filtered$j_idx,
  Duration   = data_treated_lags_filtered$Duration
)

if (K_X > 0) {
  colnames(X_treated) <- X_names
  df_calib_input = cbind(df_calib_input, X_treated)
}
X_calib = if (K_X > 0) as.matrix(df_calib_input[, X_names, drop=FALSE]) else matrix(0, nrow(df_calib_input), 0)

if (K_C > 0) {
  C_calib_full = model.matrix(~ Province, data = df_calib_input)
  C_calib = C_calib_full[, -1, drop = FALSE]
} else {
  C_calib = matrix(0, nrow(df_calib_input), 0)
}

# ---------- 2) Calibrate rho* via aux_model_rho.stan ----------
# Build Var(X|C) by province
C_levels = Provinces_List
N_C_levels = length(C_levels)
Var_X_given_C_list = vector("list", N_C_levels)

if (K_X > 0) {
  for (i in seq_along(C_levels)) {
    p = C_levels[i]
    idx_p = Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == p]
    Xp = X_matrix[idx_p, , drop = FALSE]
    V = if (nrow(Xp) > K_X + 5) tryCatch(cov(Xp), error = function(e) matrix(0, K_X, K_X)) else matrix(0, K_X, K_X)
    if (any(!is.finite(V))) V[,] = 0
    Var_X_given_C_list[[i]] = V
  }
} else {
  for (i in seq_len(N_C_levels)) Var_X_given_C_list[[i]] = matrix(0, 0, 0)
}

Var_X_given_C_array = array(0, dim = c(N_C_levels, K_X, K_X))
if (N_C_levels > 0) for (i in 1:N_C_levels) Var_X_given_C_array[i, , ] = Var_X_given_C_list[[i]]

stan_data_rho = list(
  N = nrow(df_calib_input),
  D_max = D_max,
  K_X = K_X,
  K_C = K_C,
  M = cbind(df_calib_input$M_pre, df_calib_input$M_post),
  X = X_calib,
  C = C_calib,
  Duration_idx = df_calib_input$Duration,
  N_C_levels = N_C_levels,
  Var_X_given_C = Var_X_given_C_array
)

model_rho = stan_model(file = CONFIG_PATHS$STAN_AUX_RHO)
fit_rho = sampling(model_rho, data = stan_data_rho,
                   iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS,
                   control = list(adapt_delta = CONFIG_STAN_AUX$ADAPT_DELTA, max_treedepth = CONFIG_STAN_AUX$MAX_TREEDEPTH))

posterior_rho_star_c = rstan::extract(fit_rho)$rho_star_c   # [S, N_C_levels]
posterior_rho_star   = apply(posterior_rho_star_c, 1, mean, na.rm = TRUE)

# ---------- 3) Calibrate lambda bounds via aux_model_lambda.stan ----------
model_lambda = stan_model(file = CONFIG_PATHS$STAN_AUX_LAMBDA)

stan_data_lambda_base = list(
  N = nrow(df_calib_input),
  D_max = D_max,
  K_X = K_X,
  K_C = K_C,
  X = X_calib,
  C = C_calib,
  Duration_idx = df_calib_input$Duration
)

# Y_post -> theta_2(d) = Beta_M_pre
stan_data_lambda_Ypost = stan_data_lambda_base
stan_data_lambda_Ypost$Y = df_calib_input$Y_post
stan_data_lambda_Ypost$M_post = df_calib_input$M_post
stan_data_lambda_Ypost$M_pre  = df_calib_input$M_pre

fit_lambda_Ypost = sampling(model_lambda, data = stan_data_lambda_Ypost,
                            iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)
posterior_theta_2_d = rstan::extract(fit_lambda_Ypost)$Beta_M_pre  # [S, D_max]

# Y_pre  -> zeta_1(d) = Beta_M_post
stan_data_lambda_Ypre = stan_data_lambda_base
stan_data_lambda_Ypre$Y = df_calib_input$Y_pre
stan_data_lambda_Ypre$M_post = df_calib_input$M_post
stan_data_lambda_Ypre$M_pre  = df_calib_input$M_pre

fit_lambda_Ypre = sampling(model_lambda, data = stan_data_lambda_Ypre,
                           iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)
posterior_zeta_1_d = rstan::extract(fit_lambda_Ypre)$Beta_M_post  # [S, D_max]

# ---------- 4) Save ----------
calibration_results_bayesian = list(
  Posteriors = list(
    rho_star = posterior_rho_star,
    lambda_z_lower = posterior_theta_2_d,
    lambda_minus_z_lower = posterior_zeta_1_d
  )
)
saveRDS(calibration_results_bayesian, file = CONFIG_PATHS$CALIBRATION_RESULTS)

# Console summary
rho_mean = mean(posterior_rho_star)
theta_mean = apply(posterior_theta_2_d, 2, mean)
zeta_mean  = apply(posterior_zeta_1_d, 2, mean)
summary_df = data.frame(
  Duration = 1:D_max,
  Lambda_z_LB_Mean = theta_mean,
  Lambda_minus_z_LB_Mean = zeta_mean,
  Rho_star_Mean = rho_mean
)
print(round(summary_df, 6))
