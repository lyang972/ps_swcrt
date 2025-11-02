# 3_calibration.R
# Bayesian calibration for Rho* and Lambda bounds (Eq 7-9).

# --- 0. Load Configuration and Data ---
# Assumes configuration and data are loaded by 0_driver.R
if (!exists("CONFIG_PATHS") || !exists("CONFIG_STAN_AUX") || !exists("stan_data")) {
  message("Configuration or processed data not found. Ensure 0_driver.R is executed.")
  if (sys.nframe() > 0) return(invisible(NULL))
}

# Ensure necessary variables are accessible (if loaded into a specific environment by driver)
if (!exists("K_X")) K_X = stan_data$K_X
if (!exists("K_C")) K_C = stan_data$K_C
if (!exists("D_max")) D_max = stan_data$D_max
X_names = colnames(X_matrix)

# Configure Stan for auxiliary models
options(mc.cores = CONFIG_STAN_AUX$CORES)
rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. Prepare Data for Calibration
# ==============================================================================

# Create lagged variables
data_calibration = data_model %>%
  group_by(i_idx) %>%
  arrange(time) %>%
  mutate(
    M_lag1 = lag(M_std),
    Y_lag1 = lag(Y),
    Z_lag1 = lag(Z, default = 0)
  ) %>%
  ungroup()

# Filter for treated periods with available lags
data_treated_lags = data_calibration %>%
  dplyr::filter(Z == 1 & Duration > 0) %>%
  dplyr::filter(!is.na(M_lag1) & !is.na(M_std) & !is.na(Y_lag1) & !is.na(Y))

if (nrow(data_treated_lags) == 0 || D_max < 1) {
  stop("Insufficient data for calibration.")
}

# Prepare input dataframe (df_calib_input)
i_idx_treated = data_treated_lags$i_idx
# Basic index safety
i_idx_treated = i_idx_treated[i_idx_treated <= nrow(X_matrix)]
data_treated_lags_filtered = data_treated_lags[data_treated_lags$i_idx %in% i_idx_treated, ]

df_calib_input = data.frame(
  M_pre = data_treated_lags_filtered$M_lag1,
  M_post = data_treated_lags_filtered$M_std,
  Y_pre = data_treated_lags_filtered$Y_lag1,
  Y_post = data_treated_lags_filtered$Y,
  Province = data_treated_lags_filtered$Province,
  Cluster_ID = data_treated_lags_filtered$j_idx,
  Duration = data_treated_lags_filtered$Duration
)

# Append X covariates
if (K_X > 0) {
  X_treated = X_matrix[i_idx_treated, , drop = FALSE]
  colnames(X_treated) <- X_names
  df_calib_input = cbind(df_calib_input, X_treated)
}

# Create matrices for Stan input
X_calib = if (K_X > 0) as.matrix(df_calib_input[, X_names, drop=FALSE]) else matrix(0, nrow(df_calib_input), 0)

# Handle C covariates (Province)
if (K_C > 0) {
  C_calib_full = model.matrix(~ Province, data = df_calib_input)
  C_calib = C_calib_full[, -1, drop = FALSE] # K-1 dummies
} else {
  C_calib = matrix(0, nrow(df_calib_input), 0)
}


# ==============================================================================
# 2. Calibration of Rho (ρ*) - Bayesian Multivariate Regression (Eq 7-8)
# ==============================================================================
message("Starting Bayesian calibration for Rho (ρ*)...")

# --- 2.1 Calculate Empirical Var(X|C) ---
# Required input for aux_model_rho.stan

C_levels = Provinces_List
N_C_levels = length(C_levels)
Var_X_given_C_list = list()

# Estimate Var(X|C) using the full dataset (X_matrix)
if (K_X > 0) {
  for (p in C_levels) {
    i_indices_p = Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == p]
    X_subset_p = X_matrix[i_indices_p, , drop = FALSE]
    
    # Basic check for sufficient data
    if (nrow(X_subset_p) > K_X + 5) {
      V_X_C = tryCatch(cov(X_subset_p), error = function(e) matrix(0, K_X, K_X))
    } else {
      V_X_C = matrix(0, K_X, K_X)
    }
    # Handle potential NA/NaN from cov()
    if (any(is.na(V_X_C))) V_X_C = matrix(0, K_X, K_X)
    
    Var_X_given_C_list[[as.character(p)]] = V_X_C
  }
} else {
  # Dummy structure if K_X=0
  for (p in C_levels) Var_X_given_C_list[[as.character(p)]] = matrix(0, 0, 0)
}

# Convert list to array for Stan [N_C_levels, K_X, K_X]
Var_X_given_C_array = array(dim = c(N_C_levels, K_X, K_X))
if (N_C_levels > 0) {
  for (i in 1:N_C_levels) {
    Var_X_given_C_array[i, , ] = Var_X_given_C_list[[as.character(C_levels[i])]]
  }
}


# --- 2.2 Prepare Stan Data for Rho Model ---
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

# --- 2.3 Compile and Run Rho Model ---
stan_file_rho = CONFIG_PATHS$STAN_AUX_RHO
if (!file.exists(stan_file_rho)) stop(paste("File not found:", stan_file_rho))

message("Compiling and running aux_model_rho.stan...")
model_rho = stan_model(file = stan_file_rho)

fit_rho = sampling(model_rho, data = stan_data_rho,
                   iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS,
                   control = list(adapt_delta = CONFIG_STAN_AUX$ADAPT_DELTA, max_treedepth = CONFIG_STAN_AUX$MAX_TREEDEPTH))

# Extract posterior samples of rho_star_c (Marginalized correlation conditional on C)
posterior_rho_star_c = rstan::extract(fit_rho)$rho_star_c
# Calculate the average rho* across C levels (Section 4.3).
posterior_rho_star = apply(posterior_rho_star_c, 1, mean, na.rm = TRUE)

message("Calibration for Rho finished.")

# ==============================================================================
# 3. Calibration of Lambda (λ) Lower Bounds - Bayesian Logistic Regression (Eq 9)
# ==============================================================================
message("Starting Bayesian calibration for Lambda (λ) bounds...")

# Compile Lambda Model
stan_file_lambda = CONFIG_PATHS$STAN_AUX_LAMBDA
if (!file.exists(stan_file_lambda)) stop(paste("File not found:", stan_file_lambda))

message("Compiling aux_model_lambda.stan...")
model_lambda = stan_model(file = stan_file_lambda)

# --- 3.1 Prepare Stan Data Base ---
stan_data_lambda_base = list(
  N = nrow(df_calib_input),
  D_max = D_max,
  K_X = K_X,
  K_C = K_C,
  X = X_calib,
  C = C_calib,
  Duration_idx = df_calib_input$Duration
)

# --- 3.2 Model Y_post (for theta_2(d)) ---
stan_data_lambda_Ypost = stan_data_lambda_base
stan_data_lambda_Ypost$Y = df_calib_input$Y_post
stan_data_lambda_Ypost$M_post = df_calib_input$M_post
stan_data_lambda_Ypost$M_pre = df_calib_input$M_pre

message("Running Stan sampling (Y_post -> theta_2(d))...")
fit_lambda_Ypost = sampling(model_lambda, data = stan_data_lambda_Ypost,
                            iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)

# Extract theta_2(d) (coefficient of M_pre, Beta_M_pre in Stan)
posterior_theta_2_d = rstan::extract(fit_lambda_Ypost)$Beta_M_pre

# --- 3.3 Model Y_pre (for zeta_1(d)) ---
stan_data_lambda_Ypre = stan_data_lambda_base
stan_data_lambda_Ypre$Y = df_calib_input$Y_pre
stan_data_lambda_Ypre$M_post = df_calib_input$M_post
stan_data_lambda_Ypre$M_pre = df_calib_input$M_pre

message("Running Stan sampling (Y_pre -> zeta_1(d))...")
fit_lambda_Ypre = sampling(model_lambda, data = stan_data_lambda_Ypre,
                           iter = CONFIG_STAN_AUX$ITER, warmup = CONFIG_STAN_AUX$WARMUP, chains = CONFIG_STAN_AUX$CHAINS)

# Extract zeta_1(d) (coefficient of M_post, Beta_M_post in Stan)
posterior_zeta_1_d = rstan::extract(fit_lambda_Ypre)$Beta_M_post

message("Calibration for Lambda finished.")

# ==============================================================================
# 4. Store Results
# ==============================================================================

# Store the posterior distributions for stochastic calibration
calibration_results_bayesian = list(
  Posteriors = list(
    rho_star = posterior_rho_star,                   # rho*
    lambda_z_lower = posterior_theta_2_d,            # lambda(z) lower = theta_2(d)
    lambda_minus_z_lower = posterior_zeta_1_d        # lambda(-z) lower = zeta_1(d)
  )
)

saveRDS(calibration_results_bayesian, file = CONFIG_PATHS$CALIBRATION_RESULTS)
message(paste("Bayesian calibration results saved to", CONFIG_PATHS$CALIBRATION_RESULTS))

# --- Summary Output (Posterior Means) ---
cat("\n================ Bayesian Calibration Summary (Posterior Means) =================\n")
if (length(posterior_rho_star) > 0) {
  rho_mean = mean(posterior_rho_star, na.rm = TRUE)
  theta_2_means = apply(posterior_theta_2_d, 2, mean, na.rm = TRUE)
  zeta_1_means = apply(posterior_zeta_1_d, 2, mean, na.rm = TRUE)
  
  bounds_df = data.frame(
    Duration_d = 1:D_max,
    Theta_2_d = theta_2_means,
    Zeta_1_d = zeta_1_means,
    Rho_star = rho_mean
  )
  colnames(bounds_df) = c("Duration (d)", "Theta_2(d) Mean", "Zeta_1(d) Mean", "Rho* Mean")
  print(round(bounds_df, 6))
}
cat("=================================================================================\n")