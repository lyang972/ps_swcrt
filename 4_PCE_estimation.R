# 4_PCE_estimation.R
# PCE Estimation with Amplitude Scaling (k) for Lambda sensitivity.

# ==============================================================================
# 0. Load Data and Setup
# ==============================================================================

# Assumes configuration and data are loaded by 0_driver.R
if (!exists("CONFIG_PATHS") || !exists("CONFIG_PCE") || !exists("D_max")) {
  message("Configuration or processed data not found. Ensure 0_driver.R is executed.")
  if (sys.nframe() > 0) return(invisible(NULL))
}

# Load Main Model Posteriors
if (!file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES)) stop("Posterior samples (main model) not found.")
posterior_samples = readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)

# Load Bayesian Calibration Results
if (!file.exists(CONFIG_PATHS$CALIBRATION_RESULTS)) stop("Bayesian calibration results not found.")
calibration_results_bayesian = readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)

# Ensure essential constants are defined
if (!exists("K_X") || !exists("K_C") || !exists("T_max") || !exists("sn_sd")) {
  stop("Essential constants (K_X, K_C, T_max, sn_sd) missing from environment.")
}

# --- Load Sensitivity Configuration ---
DELTA_LIST = CONFIG_PCE$DELTA_LIST
RHO_SCENARIOS = CONFIG_PCE$RHO_SCENARIOS
LAMBDA_K_SCALING = CONFIG_PCE$LAMBDA_K_SCALING

# --- Load Bayesian Calibration Posteriors ---
# These posteriors (theta_2(d), zeta_1(d)) serve as the calibrated Lambda distributions (k=1.0).

if (!is.null(calibration_results_bayesian$Posteriors)) {
  Posteriors_Aux = calibration_results_bayesian$Posteriors
  
  rho_star_posterior = Posteriors_Aux$rho_star
  # Loaded as 'lower' from 3_calibration.R, representing theta_2(d) and zeta_1(d)
  lambda_z_lower_posterior = as.matrix(Posteriors_Aux$lambda_z_lower)
  lambda_minus_z_lower_posterior = as.matrix(Posteriors_Aux$lambda_minus_z_lower)
  
  # Filter out invalid draws (NA or unstable Rho)
  valid_draws_rho = !is.na(rho_star_posterior) & (rho_star_posterior > -0.999) & (rho_star_posterior < 0.999)
  
  # Ensure dimensions match before combining validity checks
  S_aux_raw = length(valid_draws_rho)
  if (nrow(lambda_z_lower_posterior) == S_aux_raw && nrow(lambda_minus_z_lower_posterior) == S_aux_raw) {
    valid_draws_lambda_z = rowSums(is.na(lambda_z_lower_posterior)) == 0
    valid_draws_lambda_minus_z = rowSums(is.na(lambda_minus_z_lower_posterior)) == 0
    valid_draws = valid_draws_rho & valid_draws_lambda_z & valid_draws_lambda_minus_z
  } else {
    stop("Dimension mismatch between Rho and Lambda calibration posteriors.")
  }
  
  if (sum(valid_draws) == 0) {
    stop("No valid posterior samples from Bayesian calibration.")
  }
  
  # Apply filtering
  rho_star_posterior = rho_star_posterior[valid_draws]
  lambda_z_lower_posterior = lambda_z_lower_posterior[valid_draws, , drop = FALSE]
  lambda_minus_z_lower_posterior = lambda_minus_z_lower_posterior[valid_draws, , drop = FALSE]
  
  S_aux = length(rho_star_posterior)
  message(sprintf("Loaded %d valid posterior samples from Bayesian calibration.", S_aux))
  
} else {
  stop("Bayesian calibration results structure invalid.")
}


# ==============================================================================
# 1. Setup Parallel Backend and Utility Functions
# ==============================================================================

# Setup parallel backend
n_cores_to_use = if (exists("CONFIG_STAN") && !is.null(CONFIG_STAN$CORES)) {
  CONFIG_STAN$CORES
} else {
  max(1, parallel::detectCores())
}

cl = tryCatch(parallel::makeCluster(n_cores_to_use, type = "SOCK"), error = function(e) NULL)

if (!is.null(cl)) {
  doSNOW::registerDoSNOW(cl)
  message(sprintf("Parallel backend registered with %d cores.", n_cores_to_use))
} else {
  foreach::registerDoSEQ()
  message("Running sequentially.")
}


# Utility Functions
expit = function(x) {
  p <- 1 / (1 + exp(-x));
  p[x > 709] <- 1; p[x < -709] <- 0 # Handle overflow
  return(p)
}

# Gaussian-Hermite Quadrature Setup
N_gh = 20
gh = statmod::gauss.quad(N_gh, kind = "hermite")
gh_nodes = gh$nodes
gh_weights = gh$weights / sqrt(pi)

calculate_covariate_effect = function(covariates, coefficients) {
  # Simplified calculation assuming correct dimensions
  if (length(coefficients) > 0 && length(covariates) == length(coefficients)) {
    return(sum(as.vector(covariates) * as.vector(coefficients)))
  }
  return(0)
}

# ==============================================================================
# 2. Parameter Mapping Functions
# ==============================================================================

# Simplified parameter extraction for the subsetted posterior
map_params = function(s, posterior, D_max_val, K_X_val, K_C_val) {
  
  D_max <- as.integer(D_max_val)
  K_X <- as.integer(K_X_val)
  K_C <- as.integer(K_C_val)
  
  params = list()
  
  # Helper to extract the s-th draw from the subsetted object
  extract_s = function(param_obj, index) {
    if (is.null(param_obj)) return(NULL)
    dims = dim(param_obj)
    if (is.null(dims) || length(dims) == 1) return(param_obj[index])
    if (length(dims) == 2) return(param_obj[index, ])
    if (length(dims) == 3) return(param_obj[index, , ])
    return(NULL)
  }
  
  # Extract parameters
  params$eta1 = extract_s(posterior$eta1, s)
  params$eta2 = extract_s(posterior$eta2, s)
  
  if (D_max > 0) {
    params$gamma = extract_s(posterior$gamma, s)
    params$beta = extract_s(posterior$beta, s)
  } else {
    params$gamma = numeric(0); params$beta = numeric(0)
  }
  
  if (K_X > 0) {
    params$omega1 = extract_s(posterior$omega1, s)
    params$psi1 = extract_s(posterior$psi1, s)
  } else {
    params$omega1 = numeric(0); params$psi1 = numeric(0)
  }
  
  if (K_C > 0) {
    params$omega2 = extract_s(posterior$omega2, s)
    params$psi2 = extract_s(posterior$psi2, s)
  } else {
    params$omega2 = numeric(0); params$psi2 = numeric(0)
  }
  
  params$psi3 = extract_s(posterior$psi3, s)
  params$psi4 = extract_s(posterior$psi4, s)
  
  Sigma_alpha = extract_s(posterior$Sigma_alpha, s)
  Sigma_phi = extract_s(posterior$Sigma_phi, s)
  sigma_epsilon = extract_s(posterior$sigma_epsilon, s)
  
  # Derived variance components (Eq 4-5 implementation details)
  # Basic check for essential components before calculation
  if (is.null(Sigma_alpha) || is.null(Sigma_phi) || is.null(sigma_epsilon)) return(NULL)
  
  params$Var_M = Sigma_alpha[1,1] + Sigma_phi[1,1] + sigma_epsilon^2
  params$Cov_M_RE2 = Sigma_alpha[1,2] + Sigma_phi[1,2]
  params$Var_RE2 = Sigma_alpha[2,2] + Sigma_phi[2,2]
  
  # Ensure positive definiteness for numerical stability
  if (is.na(params$Var_M) || params$Var_M < 1e-9) params$Var_M = 1e-9
  
  params$Var_RE2_given_M = params$Var_RE2 - (params$Cov_M_RE2^2 / params$Var_M)
  
  if (is.na(params$Var_RE2_given_M) || params$Var_RE2_given_M < 1e-9) params$Var_RE2_given_M = 1e-9
  
  params$rho_sensitivity = NA # To be set later
  
  return(params)
}

get_history_effects = function(t, d_z, params) {
  effects = list(E_M_time=0, E_M_duration=0, E_Y_time=0, E_Y_duration=0)
  t_idx = t - 1 # Adjust index (t=2 maps to index 1)
  
  if (t_idx >= 1) {
    if (t_idx <= length(params$eta1)) effects$E_M_time = params$eta1[t_idx]
    if (t_idx <= length(params$eta2)) effects$E_Y_time = params$eta2[t_idx]
  }
  
  if (d_z > 0) {
    if (d_z <= length(params$gamma)) effects$E_M_duration = params$gamma[d_z]
    if (d_z <= length(params$beta)) effects$E_Y_duration = params$beta[d_z]
  }
  
  return(effects)
}

# ==============================================================================
# 3. Core Calculation Functions (Methodology Implementation)
# ==============================================================================

# Calculates E[Y(z)|M(z), X, C] using GH Quadrature (Section 4.2)
calc_E_Y_given_M_Z_X_C = function(t, d_z, m, params, covariates_X, covariates_C) {
  history_effects = get_history_effects(t, d_z, params)
  X_eff = calculate_covariate_effect(covariates_X, params$psi1)
  C_eff = calculate_covariate_effect(covariates_C, params$psi2)
  
  # Mean of Y (linear predictor)
  mu_Y = history_effects$E_Y_time + history_effects$E_Y_duration + X_eff + C_eff + m * params$psi3
  if (d_z > 0) {
    mu_Y = mu_Y + m * params$psi4
  }
  
  # Mean of M (for conditioning)
  X_eff_M = calculate_covariate_effect(covariates_X, params$omega1)
  C_eff_M = calculate_covariate_effect(covariates_C, params$omega2)
  E_M_z = history_effects$E_M_time + history_effects$E_M_duration + X_eff_M + C_eff_M
  
  # Conditional expectation of RE2 given M
  E_RE2_given_M = (params$Cov_M_RE2 / params$Var_M) * (m - E_M_z)
  
  # GH Quadrature integration
  sd_RE2_given_M = sqrt(params$Var_RE2_given_M)
  nodes_transformed = E_RE2_given_M + sqrt(2) * sd_RE2_given_M * gh_nodes
  
  integrand_values = expit(mu_Y + nodes_transformed)
  E_Y_given_M = sum(gh_weights * integrand_values)
  
  return(E_Y_given_M)
}

# Solves the convolution equation (Eq 6) for Delta using Newton-Raphson
solve_Delta = function(t, d_z, d_z_star, m_z, E_Y_obs, params, covariates_X, covariates_C, lambda_diff) {
  
  # Setup for P(M(z*)|M(z)) calculation
  X_eff_M = calculate_covariate_effect(covariates_X, params$omega1)
  C_eff_M = calculate_covariate_effect(covariates_C, params$omega2)
  
  history_z = get_history_effects(t, d_z, params)
  E_M_z = history_z$E_M_time + history_z$E_M_duration + X_eff_M + C_eff_M
  
  history_z_star = get_history_effects(t, d_z_star, params)
  E_M_z_star = history_z_star$E_M_time + history_z_star$E_M_duration + X_eff_M + C_eff_M
  
  rho = params$rho_sensitivity
  Var_M = params$Var_M
  
  # Conditional distribution P(M(z*)|M(z)) (Assumption 5)
  E_M_z_star_given_z = E_M_z_star + rho * (m_z - E_M_z)
  Var_M_z_star_given_z = (1 - rho^2) * Var_M
  
  if (Var_M_z_star_given_z < 1e-9) Var_M_z_star_given_z = 1e-9
  SD_M_z_star_given_z = sqrt(Var_M_z_star_given_z)
  
  # GH Quadrature nodes for M(z*)
  nodes_M_star = E_M_z_star_given_z + sqrt(2) * SD_M_z_star_given_z * gh_nodes
  
  # Newton-Raphson
  MAX_ITER = 50; TOL = 1e-7
  
  # Initialize Delta (logit of observed E[Y|M])
  E_Y_obs_clipped = pmin(1 - 1e-9, pmax(1e-9, E_Y_obs))
  Delta_current = log(E_Y_obs_clipped / (1 - E_Y_obs_clipped))
  
  for (iter in 1:MAX_ITER) {
    # Calculate the integral (f(Delta))
    linear_predictor = Delta_current + lambda_diff * nodes_M_star
    integrand_values = expit(linear_predictor)
    Integral_f = sum(gh_weights * integrand_values)
    
    # Objective function: -E_Y_obs + Integral_f = 0
    f_Delta = -E_Y_obs + Integral_f
    
    if (abs(f_Delta) < TOL) return(Delta_current)
    
    # Calculate the derivative (f'(Delta))
    derivative_integrand = integrand_values * (1 - integrand_values)
    f_prime_Delta = sum(gh_weights * derivative_integrand)
    
    if (abs(f_prime_Delta) < 1e-9) break
    
    # Update Delta
    Delta_current = Delta_current - f_Delta / f_prime_Delta
  }
  
  return(Delta_current) # Return last estimate if max iterations reached
}

# Calculates PCE conditional on X, C using Monte Carlo integration (Theorem 1)
calculate_PCE_XC = function(t, d_z, d_z_star, I, params, covariates_X, covariates_C, lambda_z, lambda_minus_z, N_mc) {
  
  # Setup for P(M(z), M(z*)) calculation
  X_eff_M = calculate_covariate_effect(covariates_X, params$omega1)
  C_eff_M = calculate_covariate_effect(covariates_C, params$omega2)
  
  history_z = get_history_effects(t, d_z, params)
  E_M_z = history_z$E_M_time + history_z$E_M_duration + X_eff_M + C_eff_M
  
  history_z_star = get_history_effects(t, d_z_star, params)
  E_M_z_star = history_z_star$E_M_time + history_z_star$E_M_duration + X_eff_M + C_eff_M
  
  rho = params$rho_sensitivity
  Var_M = params$Var_M
  
  # Define the joint distribution P(M(z), M(z*)) (Assumption 5)
  Mean_Vector = c(E_M_z, E_M_z_star)
  Cov_Matrix = matrix(c(Var_M, rho*Var_M, rho*Var_M, Var_M), nrow=2)
  
  # Monte Carlo Integration: Sample from the joint distribution
  M_samples = tryCatch({
    mvtnorm::rmvnorm(N_mc, mean = Mean_Vector, sigma = Cov_Matrix)
  }, error = function(e) { return(NULL) })
  
  if (is.null(M_samples)) return(NA)
  
  m_z_samples = M_samples[, 1]
  m_z_star_samples = M_samples[, 2]
  
  # Identify samples within the stratum I
  M_diff_samples = m_z_samples - m_z_star_samples
  in_stratum = (M_diff_samples >= I[1]) & (M_diff_samples < I[2])
  
  N_in_stratum = sum(in_stratum)
  if (N_in_stratum == 0) return(0)
  
  m_z_subset = m_z_samples[in_stratum]
  m_z_star_subset = m_z_star_samples[in_stratum]
  
  # Calculate Delta for each sample in the stratum
  Delta_z = numeric(N_in_stratum)
  Delta_z_star = numeric(N_in_stratum)
  
  for (k in 1:N_in_stratum) {
    # Step 1: Calculate E[Y(z)|M(z)] (observed marginal)
    E_Y_obs_z = calc_E_Y_given_M_Z_X_C(t, d_z, m_z_subset[k], params, covariates_X, covariates_C)
    E_Y_obs_z_star = calc_E_Y_given_M_Z_X_C(t, d_z_star, m_z_star_subset[k], params, covariates_X, covariates_C)
    
    # Step 2: Solve for Delta(m, z)
    Delta_z[k] = solve_Delta(t, d_z, d_z_star, m_z_subset[k], E_Y_obs_z, params, covariates_X, covariates_C, lambda_z)
    Delta_z_star[k] = solve_Delta(t, d_z_star, d_z, m_z_star_subset[k], E_Y_obs_z_star, params, covariates_X, covariates_C, lambda_minus_z)
  }
  
  # Step 3: Calculate E[Y(z)|M(z), M(z*)] using MSM (Assumption 6)
  E_Y_z = expit(Delta_z + lambda_z * m_z_star_subset)
  E_Y_z_star = expit(Delta_z_star + lambda_minus_z * m_z_subset)
  
  # Approximate PCE by averaging the difference
  PCE_approx = mean(E_Y_z - E_Y_z_star, na.rm = TRUE)
  
  return(PCE_approx)
}

# Calculates the probability of belonging to stratum I
calculate_Stratum_Probability = function(I, d_z, d_z_star, params) {
  rho = params$rho_sensitivity
  Var_M = params$Var_M
  
  # Determine the difference in duration effects (gamma)
  D_max_params = length(params$gamma)
  gamma_dz = if (d_z > 0 && d_z <= D_max_params) params$gamma[d_z] else 0
  gamma_dz_star = if (d_z_star > 0 && d_z_star <= D_max_params) params$gamma[d_z_star] else 0
  
  # Distribution of M(z) - M(z*)
  Mean_diff = gamma_dz - gamma_dz_star
  Var_diff = 2 * (1 - rho) * Var_M
  
  if (Var_diff < 1e-9) {
    return(ifelse(Mean_diff >= I[1] && Mean_diff < I[2], 1.0, 0.0))
  }
  
  SD_diff = sqrt(Var_diff)
  
  # Calculate probability P(I[1] <= Diff < I[2])
  Prob = pnorm(I[2], mean = Mean_diff, sd = SD_diff) - pnorm(I[1], mean = Mean_diff, sd = SD_diff)
  return(Prob)
}


# ==============================================================================
# 4. Main Execution Loop (Algorithm 1 Implementation)
# ==============================================================================

# --- Define Estimands (Duration L vs 0) ---
estimands_list = list()
DURATIONS_TO_EVALUATE = CONFIG_PCE$DURATIONS

for (L in DURATIONS_TO_EVALUATE) {
  if (L > 0 && L <= D_max) {
    min_time = L + 1
    if (min_time <= T_max) {
      for (t in min_time:T_max) {
        Type_label = ifelse(L == 1, "stPCE (d=1)", paste0("ltPCE (d=", L, ")"))
        est_name = paste0("PCE_d", L, "_t", t)
        estimands_list[[est_name]] = list(t=t, d_z=L, d_z_star=0, Duration=L, Type=Type_label)
      }
    }
  }
}

# Setup for computation
N_posterior = length(posterior_samples$psi3)
S_subset = min(CONFIG_PCE$S_SUBSET, N_posterior)
N_individuals_sample_target = CONFIG_PCE$N_INDIV_SAMPLE
N_MC_Integration = CONFIG_PCE$N_MC_INTEGRATION

set.seed(123)
idx_S = sample(1:N_posterior, S_subset)

# Prepare Province-Specific Individual Indices
Province_Individual_Map = list()
for (P in Provinces_List) {
  idx_P = Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == P]
  if (length(idx_P) > 0) Province_Individual_Map[[as.character(P)]] = idx_P
}

# --- Subset posterior samples before parallelization (Memory Optimization) ---
message("Subsetting posterior samples...")

# Function to subset draws across different dimensions
subset_posterior_draws = function(param, indices) {
  if (is.null(param)) return(NULL)
  dims = dim(param)
  if (is.matrix(param)) return(param[indices, , drop = FALSE])
  if (is.array(param) && length(dims) == 3) return(param[indices, , , drop = FALSE])
  if (is.numeric(param) && length(param) >= max(indices)) return(param[indices])
  return(NULL)
}

posterior_samples_subset = lapply(posterior_samples, subset_posterior_draws, indices = idx_S)
posterior_samples_subset = Filter(Negate(is.null), posterior_samples_subset)

rm(posterior_samples)
gc()

message(sprintf("Starting PCE calculation (S=%d, N_MC=%d)...", S_subset, N_MC_Integration))

# Define constants for stratum definition
LARGE_NUM = 1e9

# Main loop (Parallelized over posterior samples 's')
all_results = foreach::foreach(s_idx = 1:S_subset, .combine = rbind,
                               .packages = c("statmod", "mvtnorm", "dplyr")
) %dopar% {
  
  # Set seed inside worker for reproducibility
  set.seed(s_idx + 12345)
  
  # Use tryCatch for robustness in parallel execution
  tryCatch({
    
    s = s_idx # Index relative to the subset
    original_draw_s = idx_S[s_idx]
    
    # 1. Extract parameters from Main Model
    params = map_params(s, posterior_samples_subset, D_max, K_X, K_C)
    if (is.null(params)) return(NULL)
    
    results_s = list()
    idx = 1
    
    # 2. Sample Calibrated Parameters from Auxiliary Models
    # Sample one realization from the Bayesian auxiliary posteriors (S_aux)
    s_aux_idx = sample(1:S_aux, 1)
    
    # Retrieve the sampled duration-independent Rho*
    rho_star_s = rho_star_posterior[s_aux_idx]
    
    # ==========================================================================
    # Loop over Estimands
    # ==========================================================================
    
    for (est_name in names(estimands_list)) {
      est = estimands_list[[est_name]]
      t = est$t; d_z = est$d_z; d_z_star = est$d_z_star
      d_z_idx = d_z
      
      # 2.1 Retrieve Duration-Specific Sampled Calibrated Lambda (Baseline k=1.0)
      # These correspond to the sampled theta_2(d) and zeta_1(d).
      
      if (d_z_idx > 0 && d_z_idx <= ncol(lambda_z_lower_posterior)) {
        L_z_calib_s = lambda_z_lower_posterior[s_aux_idx, d_z_idx]
        L_minus_z_calib_s = lambda_minus_z_lower_posterior[s_aux_idx, d_z_idx]
      } else {
        L_z_calib_s = NA; L_minus_z_calib_s = NA
      }
      
      # Fallback robustness
      if (is.na(L_z_calib_s)) L_z_calib_s = 0
      if (is.na(L_minus_z_calib_s)) L_minus_z_calib_s = 0
      
      # 3. Loop over Rho Scenarios
      for (CURRENT_RHO_SCENARIO in RHO_SCENARIOS) {
        
        # Determine Rho value
        if (CURRENT_RHO_SCENARIO == "Calibrated_Stochastic") {
          current_rho_value = rho_star_s
        } else if (startsWith(CURRENT_RHO_SCENARIO, "Fixed_")) {
          fixed_rho = as.numeric(sub("Fixed_", "", CURRENT_RHO_SCENARIO))
          if (!is.na(fixed_rho)) current_rho_value = fixed_rho else next
        } else {
          next
        }
        
        params$rho_sensitivity = current_rho_value
        
        # 4. Loop over Delta Values
        for (CURRENT_DELTA in DELTA_LIST) {
          
          # Define Strata (Standardized scale)
          delta_std = CURRENT_DELTA / sn_sd
          
          I1 = c(-delta_std, delta_std)        # Dissociative
          I2 = c(-LARGE_NUM, -delta_std)       # Associative Negative
          I3 = c(delta_std, LARGE_NUM)         # Associative Positive
          strata_list = list(PCE1_Dissociative=I1, PCE2_Associative_Neg=I2, PCE3_Associative_Pos=I3)
          
          # 5. Loop over Lambda K Scaling Factors
          for (CURRENT_K in LAMBDA_K_SCALING) {
            
            # Apply amplitude scaling: lambda_sensitivity = k * lambda_calibrated
            lambda_z = CURRENT_K * L_z_calib_s
            lambda_minus_z = CURRENT_K * L_minus_z_calib_s
            
            # 6. PCE Computation
            # Loop over strata
            for (stratum_name in names(strata_list)) {
              stratum_I = strata_list[[stratum_name]]
              
              # Calculate Stratum Probability
              current_stratum_prob = calculate_Stratum_Probability(stratum_I, d_z, d_z_star, params)
              if (is.na(current_stratum_prob)) next
              
              # --- 6.1 Marginalize over X, Conditional on C (Province) ---
              for (P in names(Province_Individual_Map)) {
                
                # Sample individuals within the province
                idx_P = Province_Individual_Map[[P]]
                N_P = length(idx_P)
                if (N_P > N_individuals_sample_target) {
                  idx_sample = sample(idx_P, N_individuals_sample_target)
                  N_sample = N_individuals_sample_target
                } else {
                  idx_sample = idx_P
                  N_sample = N_P
                }
                
                PCE_i_list = numeric(N_sample)
                
                # Loop over sampled individuals
                for (k in 1:N_sample) {
                  i = idx_sample[k]
                  
                  # Get individual covariates X_i
                  if (K_X > 0 && i <= nrow(X_matrix)) {
                    covariates_X = X_matrix[i, ]
                  } else {
                    covariates_X = numeric(0)
                  }
                  
                  # Get cluster covariates C_j
                  j = Individual_to_Cluster_Map$j_idx[Individual_to_Cluster_Map$i_idx == i][1]
                  
                  if (K_C > 0 && !is.na(j) && j <= nrow(C_matrix)) {
                    covariates_C = C_matrix[j, ]
                  } else {
                    covariates_C = numeric(0)
                  }
                  
                  # Calculate PCE(X_i, C_j)
                  PCE_i = calculate_PCE_XC(t, d_z, d_z_star, stratum_I, params,
                                           covariates_X, covariates_C,
                                           lambda_z, lambda_minus_z, N_MC_Integration)
                  
                  PCE_i_list[k] = PCE_i
                }
                
                # Average over individuals (Marginalize X conditional on C)
                PCE_conditional_C = mean(PCE_i_list, na.rm = TRUE)
                
                # Store result
                results_s[[idx]] = data.frame(
                  Draw=original_draw_s, Estimand=est_name, Stratum=stratum_name,
                  Rho_Scenario=CURRENT_RHO_SCENARIO,
                  Province=P, PCE=PCE_conditional_C,
                  Delta=CURRENT_DELTA,
                  Lambda_K=CURRENT_K,
                  Stratum_Probability=current_stratum_prob,
                  Lambda_z=lambda_z, Lambda_minus_z=lambda_minus_z
                )
                idx = idx + 1
              } # End Province
            } # End Stratum
          } # End Lambda K Scaling
        } # End Delta
      } # End Rho Scenario
    } # End Estimands
    
    # Combine results for this iteration 's'
    do.call(rbind, results_s)
    
  }, error = function(e) {
    # Handle errors during the iteration
    cat(sprintf("Error during iteration %d: %s\n", s_idx, e$message))
    return(NULL)
  })
}

# Stop parallel cluster
if (!is.null(cl)) {
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
}

message("PCE calculation finished.")

# ==============================================================================
# 5. Summarize Results
# ==============================================================================

PCE_summary <- data.frame()
Probability_summary <- data.frame()
Lambda_summary <- data.frame()

if (!is.null(all_results) && nrow(all_results) > 0) {
  
  # Prepare estimand definitions for joining
  estimands_definitions = dplyr::bind_rows(estimands_list, .id="Estimand") %>%
    dplyr::select(Estimand, Time = t, Duration, Type)
  
  # Calculate posterior summaries for PCE
  PCE_summary = all_results %>%
    dplyr::filter(!is.na(PCE) & !is.nan(PCE)) %>%
    dplyr::left_join(estimands_definitions, by = "Estimand") %>%
    # Group by all defining characteristics
    dplyr::group_by(Estimand, Stratum, Time, Duration, Type, Province, Delta, Lambda_K, Rho_Scenario) %>%
    dplyr::summarise(
      Mean = mean(PCE),
      Median = stats::median(PCE),
      Q2.5 = stats::quantile(PCE, 0.025),
      Q97.5 = stats::quantile(PCE, 0.975),
      .groups = 'drop'
    )
  
  # Calculate posterior summaries for Stratum Probabilities
  Probability_summary = all_results %>%
    dplyr::filter(!is.na(Stratum_Probability)) %>%
    dplyr::left_join(estimands_definitions %>% select(Estimand, Duration, Type) %>% distinct(), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Type, Stratum, Delta, Rho_Scenario, Stratum_Probability) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Duration, Type, Stratum, Delta, Rho_Scenario) %>%
    dplyr::summarise(
      Prob_Mean = mean(Stratum_Probability, na.rm = TRUE),
      Prob_Median = stats::median(Stratum_Probability, na.rm = TRUE),
      Prob_Q2.5 = stats::quantile(Stratum_Probability, 0.025, na.rm = TRUE),
      Prob_Q97.5 = stats::quantile(Stratum_Probability, 0.975, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Summarize Lambda distributions
  Lambda_summary = all_results %>%
    dplyr::left_join(estimands_definitions %>% select(Estimand, Duration) %>% distinct(), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Lambda_K, Lambda_z, Lambda_minus_z) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(cols = c(Lambda_z, Lambda_minus_z), names_to = "Lambda_Type", values_to = "Value")
}

# Save results
save(all_results, PCE_summary, Probability_summary, Lambda_summary, estimands_list, file = CONFIG_PATHS$PCE_RESULTS)
message("PCE results saved.")