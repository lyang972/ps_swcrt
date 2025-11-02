# 4_PCE_estimation.R
# PCE estimation with amplitude-scaling k for lambda; identical outputs, cleaner parallel.

load("Analysis_Config.RData")
load(CONFIG_PATHS$PROCESSED_DATA)

posterior_samples = readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)
calibration_results_bayesian = readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)

stopifnot(exists("K_X"), exists("K_C"), exists("D_max"), exists("T_max"),
          exists("sn_sd"), is.finite(sn_sd), sn_sd > 0)

DELTA_LIST = CONFIG_PCE$DELTA_LIST
RHO_SCENARIOS = CONFIG_PCE$RHO_SCENARIOS
LAMBDA_K_SCALING = CONFIG_PCE$LAMBDA_K_SCALING

Posteriors_Aux = calibration_results_bayesian$Posteriors
rho_star_posterior = Posteriors_Aux$rho_star
lambda_z_lower_posterior = as.matrix(Posteriors_Aux$lambda_z_lower)
lambda_minus_z_lower_posterior = as.matrix(Posteriors_Aux$lambda_minus_z_lower)

valid_draws_rho = is.finite(rho_star_posterior) & (rho_star_posterior > -0.999) & (rho_star_posterior < 0.999)
if (nrow(lambda_z_lower_posterior) == length(valid_draws_rho)) {
  valid_draws_lambda_z = rowSums(!is.finite(lambda_z_lower_posterior)) == 0
  valid_draws_lambda_mz = rowSums(!is.finite(lambda_minus_z_lower_posterior)) == 0
  valid_draws = valid_draws_rho & valid_draws_lambda_z & valid_draws_lambda_mz
} else {
  m = min(length(valid_draws_rho), nrow(lambda_z_lower_posterior), nrow(lambda_minus_z_lower_posterior))
  valid_draws = valid_draws_rho[1:m] &
    (rowSums(!is.finite(lambda_z_lower_posterior[1:m, , drop=FALSE])) == 0) &
    (rowSums(!is.finite(lambda_minus_z_lower_posterior[1:m, , drop=FALSE])) == 0)
}
stopifnot(sum(valid_draws) > 0)

rho_star_posterior = rho_star_posterior[valid_draws]
lambda_z_lower_posterior = lambda_z_lower_posterior[valid_draws, , drop = FALSE]
lambda_minus_z_lower_posterior = lambda_minus_z_lower_posterior[valid_draws, , drop = FALSE]
S_aux = length(rho_star_posterior)

# ---------- parallel backend (PSOCK) ----------
n_cores_to_use = if (!is.null(CONFIG_STAN$CORES)) CONFIG_STAN$CORES else max(1, parallel::detectCores())
cl = tryCatch(parallel::makeCluster(n_cores_to_use, type = "PSOCK", outfile = ""), error = function(e) NULL)
if (!is.null(cl)) doSNOW::registerDoSNOW(cl) else foreach::registerDoSEQ()

# ---------- utilities ----------
expit = function(x) { p <- 1/(1+exp(-x)); p[x>709]=1; p[x< -709]=0; p }

N_gh = 20
gh = statmod::gauss.quad(N_gh, kind = "hermite")
gh_nodes = gh$nodes
gh_weights = gh$weights / sqrt(pi)

calculate_covariate_effect = function(covariates, coefficients) {
  if (length(coefficients) > 0 && length(covariates) > 0) {
    v = as.vector(covariates); b = as.vector(coefficients)
    if (length(v) == length(b)) return(sum(v * b))
  }
  0
}

# ---------- parameter mapping ----------
map_params = function(s, posterior, D_max_val, K_X_val, K_C_val) {
  D_max <- as.integer(D_max_val); K_X <- as.integer(K_X_val); K_C <- as.integer(K_C_val)
  ext = function(obj, s) {
    if (is.null(obj)) return(numeric(0))
    d = dim(obj)
    if (is.null(d)) return(obj[s])
    if (length(d) == 2) return(obj[s, ])
    if (length(d) == 3) return(obj[s, , ])
    numeric(0)
  }
  params = list()
  params$eta1  = ext(posterior$eta1, s)
  params$eta2  = ext(posterior$eta2, s)
  params$gamma = if (D_max > 0) ext(posterior$gamma, s) else numeric(0)
  params$beta  = if (D_max > 0) ext(posterior$beta, s)  else numeric(0)
  params$omega1 = if (K_X > 0) ext(posterior$omega1, s) else numeric(0)
  params$psi1   = if (K_X > 0) ext(posterior$psi1, s)   else numeric(0)
  params$omega2 = if (K_C > 0) ext(posterior$omega2, s) else numeric(0)
  params$psi2   = if (K_C > 0) ext(posterior$psi2, s)   else numeric(0)
  params$psi3 = ext(posterior$psi3, s)
  params$psi4 = ext(posterior$psi4, s)
  Sigma_alpha = ext(posterior$Sigma_alpha, s)
  Sigma_phi   = ext(posterior$Sigma_phi, s)
  sigma_eps   = ext(posterior$sigma_epsilon, s)

  if (!is.matrix(Sigma_alpha) || !is.matrix(Sigma_phi) || !is.finite(sigma_eps)) return(NULL)
  VM = Sigma_alpha[1,1] + Sigma_phi[1,1] + sigma_eps^2
  if (!is.finite(VM) || VM <= 0) VM = 1e-9
  params$Var_M = VM
  params$Cov_M_RE2 = Sigma_alpha[1,2] + Sigma_phi[1,2]
  VR2 = Sigma_alpha[2,2] + Sigma_phi[2,2]
  params$Var_RE2_given_M = max(VR2 - (params$Cov_M_RE2^2 / VM), 1e-9)
  params$rho_sensitivity = NA_real_
  params
}

get_history_effects = function(t, d_z, params) {
  t_idx = t - 1
  list(
    E_M_time = if (t_idx >= 1 && t_idx <= length(params$eta1)) params$eta1[t_idx] else 0,
    E_M_duration = if (d_z > 0 && d_z <= length(params$gamma)) params$gamma[d_z] else 0,
    E_Y_time = if (t_idx >= 1 && t_idx <= length(params$eta2)) params$eta2[t_idx] else 0,
    E_Y_duration = if (d_z > 0 && d_z <= length(params$beta))  params$beta[d_z]  else 0
  )
}

# ---------- core functions ----------
calc_E_Y_given_M_Z_X_C = function(t, d_z, m, params, covariates_X, covariates_C) {
  h = get_history_effects(t, d_z, params)
  X_eff = calculate_covariate_effect(covariates_X, params$psi1)
  C_eff = calculate_covariate_effect(covariates_C, params$psi2)
  mu_Y = h$E_Y_time + h$E_Y_duration + X_eff + C_eff + m * params$psi3
  if (d_z > 0) mu_Y = mu_Y + m * params$psi4

  X_eff_M = calculate_covariate_effect(covariates_X, params$omega1)
  C_eff_M = calculate_covariate_effect(covariates_C, params$omega2)
  E_M_z = h$E_M_time + h$E_M_duration + X_eff_M + C_eff_M

  E_RE2_cond = (params$Cov_M_RE2 / params$Var_M) * (m - E_M_z)
  sd_RE2 = sqrt(params$Var_RE2_given_M)
  nodes = E_RE2_cond + sqrt(2) * sd_RE2 * gh_nodes
  sum(gh_weights * expit(mu_Y + nodes))
}

solve_Delta = function(t, d_z, d_z_star, m_z, E_Y_obs, params, X, C, lambda_diff) {
  X_eff_M = calculate_covariate_effect(X, params$omega1)
  C_eff_M = calculate_covariate_effect(C, params$omega2)
  h_z  = get_history_effects(t, d_z, params)
  h_zs = get_history_effects(t, d_z_star, params)
  E_M_z  = h_z$E_M_time  + h_z$E_M_duration  + X_eff_M + C_eff_M
  E_M_zs = h_zs$E_M_time + h_zs$E_M_duration + X_eff_M + C_eff_M

  rho = params$rho_sensitivity; VM = params$Var_M
  if (!is.finite(rho) || abs(rho) > 1) return(NA_real_)
  mu = E_M_zs + rho * (m_z - E_M_z)
  v  = (1 - rho^2) * VM
  if (v <= 1e-9) v = 1e-9
  sd = sqrt(v)
  nodes = mu + sqrt(2) * sd * gh_nodes

  MAXI = 50; TOL = 1e-7
  Ey = pmin(1-1e-9, pmax(1e-9, E_Y_obs))
  Delta = log(Ey / (1 - Ey))
  for (it in 1:MAXI) {
    lp = Delta + lambda_diff * nodes
    sig = expit(lp)
    f  = -E_Y_obs + sum(gh_weights * sig)
    if (abs(f) < TOL) return(Delta)
    fp = sum(gh_weights * sig * (1 - sig))
    if (abs(fp) < 1e-9) break
    Delta = Delta - f / fp
  }
  Delta
}

calculate_PCE_XC = function(t, d_z, d_z_star, I, params, X, C, lambda_z, lambda_mz, N_mc) {
  X_eff_M = calculate_covariate_effect(X, params$omega1)
  C_eff_M = calculate_covariate_effect(C, params$omega2)
  h_z  = get_history_effects(t, d_z, params)
  h_zs = get_history_effects(t, d_z_star, params)
  E_M_z  = h_z$E_M_time  + h_z$E_M_duration  + X_eff_M + C_eff_M
  E_M_zs = h_zs$E_M_time + h_zs$E_M_duration + X_eff_M + C_eff_M

  rho = params$rho_sensitivity; VM = params$Var_M
  if (!is.finite(rho) || abs(rho) > 1) return(NA_real_)

  mu = c(E_M_z, E_M_zs)
  Sigma = matrix(c(VM, rho*VM, rho*VM, VM), 2)
  M = tryCatch(mvtnorm::rmvnorm(N_mc, mean = mu, sigma = Sigma), error = function(e) NULL)
  if (is.null(M)) return(NA_real_)
  m_z  = M[,1]; m_zs = M[,2]
  dM = m_z - m_zs
  sel = (dM >= I[1]) & (dM < I[2])
  if (!any(sel)) return(0)

  m1 = m_z[sel]; m2 = m_zs[sel]
  n <- length(m1)
  Ey1 <- numeric(n); Ey0 <- numeric(n)
  D1  <- numeric(n); D0  <- numeric(n)
  for (k in seq_len(n)) {
    Ey1[k] <- calc_E_Y_given_M_Z_X_C(t, d_z,      m1[k], params, X, C)
    Ey0[k] <- calc_E_Y_given_M_Z_X_C(t, d_z_star, m2[k], params, X, C)
    D1[k]  <- solve_Delta(t, d_z,      d_z_star, m1[k], Ey1[k], params, X, C, lambda_z)
    D0[k]  <- solve_Delta(t, d_z_star, d_z,      m2[k], Ey0[k], params, X, C, lambda_mz)
  }
  if (any(!is.finite(D1)) || any(!is.finite(D0))) return(NA_real_)
  mean(expit(D1 + lambda_z  * m2) - expit(D0 + lambda_mz * m1))
}

calculate_Stratum_Probability = function(I, d_z, d_z_star, params) {
  rho = params$rho_sensitivity; VM = params$Var_M
  if (!is.finite(rho) || abs(rho) > 1) return(NA_real_)
  g1 = if (d_z      > 0 && d_z      <= length(params$gamma)) params$gamma[d_z]      else 0
  g0 = if (d_z_star > 0 && d_z_star <= length(params$gamma)) params$gamma[d_z_star] else 0
  mu = g1 - g0
  var = 2 * (1 - rho) * VM
  if (var <= 1e-9) return(as.numeric(mu >= I[1] && mu < I[2]))
  sd = sqrt(var)
  pnorm(I[2], mean = mu, sd = sd) - pnorm(I[1], mean = mu, sd = sd)
}

# ---------- estimands ----------
estimands_list = list()
for (L in CONFIG_PCE$DURATIONS) {
  if (L > 0 && L <= D_max) {
    min_t = L + 1
    if (min_t <= T_max) {
      for (t in min_t:T_max) {
        estimands_list[[paste0("PCE_d", L, "_t", t)]] =
          list(t=t, d_z=L, d_z_star=0, Duration=L, Type=ifelse(L==1, "stPCE (d=1)", paste0("ltPCE (d=",L,")")))
      }
    }
  }
}
stopifnot(length(estimands_list) > 0)

# subset posterior from main model
N_posterior = length(posterior_samples$psi3)
S_subset = min(CONFIG_PCE$S_SUBSET, N_posterior)
set.seed(123)
idx_S = sample(1:N_posterior, S_subset)

subset_draws = function(x, idx) {
  if (is.null(x)) return(NULL)
  d = dim(x)
  if (is.matrix(x) && nrow(x) >= max(idx)) return(x[idx, , drop=FALSE])
  if (is.array(x)  && length(d) == 3 && d[1] >= max(idx)) return(x[idx, , , drop=FALSE])
  if (is.numeric(x) && length(x) >= max(idx)) return(x[idx])
  NULL
}
posterior_samples_subset = lapply(posterior_samples, subset_draws, idx = idx_S)
posterior_samples_subset = Filter(Negate(is.null), posterior_samples_subset)
rm(posterior_samples); gc()

# province -> individuals
Province_Individual_Map = lapply(as.list(Provinces_List), function(P) {
  Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == P]
})
names(Province_Individual_Map) = as.character(Provinces_List)

N_indiv_target = CONFIG_PCE$N_INDIV_SAMPLE
N_MC = CONFIG_PCE$N_MC_INTEGRATION
LARGE = 1e9

# progress (optional)
progress_fun <- function(n) cat(sprintf("Progress: %d/%d\n", n, S_subset))

# ---------- main loop (parallel over s) ----------
all_results = foreach::foreach(
  s_idx = 1:S_subset, .combine = rbind,
  .packages = c("statmod","mvtnorm","dplyr","triangle"),
  .export = c("map_params","get_history_effects","calc_E_Y_given_M_Z_X_C","solve_Delta",
              "calculate_PCE_XC","calculate_Stratum_Probability","expit","gh_nodes","gh_weights"),
  .options.snow = if (!is.null(cl)) list(progress = progress_fun) else NULL
) %dopar% {
  set.seed(12345 + s_idx)
  s = s_idx
  original_draw = idx_S[s_idx]
  params = map_params(s, posterior_samples_subset, D_max, K_X, K_C)
  if (is.null(params)) return(NULL)

  # draw calibrated rho*, lambda at duration d
  s_aux = sample(1:S_aux, 1)
  rho_star_s = rho_star_posterior[s_aux]

  out = list(); k_res = 1

  for (nm in names(estimands_list)) {
    est = estimands_list[[nm]]
    t = est$t; d1 = est$d_z; d0 = est$d_z_star
    d_idx = d1

    Lz  = if (d_idx > 0 && d_idx <= ncol(lambda_z_lower_posterior))      lambda_z_lower_posterior[s_aux, d_idx]      else 0
    Lmz = if (d_idx > 0 && d_idx <= ncol(lambda_minus_z_lower_posterior)) lambda_minus_z_lower_posterior[s_aux, d_idx] else 0
    if (!is.finite(Lz))  Lz  = 0
    if (!is.finite(Lmz)) Lmz = 0
    if (!is.finite(rho_star_s) || abs(rho_star_s) >= 1) rho_star_s = 0.5

    for (rho_scn in RHO_SCENARIOS) {
      rho_val = if (rho_scn == "Calibrated_Stochastic") rho_star_s else as.numeric(sub("Fixed_", "", rho_scn))
      if (!is.finite(rho_val) || abs(rho_val) > 1) next
      params$rho_sensitivity = rho_val

      for (DELTA in DELTA_LIST) {
        delta_std = DELTA / sn_sd
        I1 = c(-delta_std,  delta_std)
        I2 = c(-LARGE,     -delta_std)
        I3 = c( delta_std,  LARGE)
        strata = list(PCE1_Dissociative=I1, PCE2_Associative_Neg=I2, PCE3_Associative_Pos=I3)

        for (Ksc in LAMBDA_K_SCALING) {
          lam_z  = Ksc * Lz
          lam_mz = Ksc * Lmz

          for (str_nm in names(strata)) {
            I = strata[[str_nm]]
            prob_str = calculate_Stratum_Probability(I, d1, d0, params)

            for (P in names(Province_Individual_Map)) {
              idxP = Province_Individual_Map[[P]]
              if (length(idxP) == 0) next
              if (length(idxP) > N_indiv_target) idxP = sample(idxP, N_indiv_target)

              Pc = numeric(length(idxP))
              for (k in seq_along(idxP)) {
                i = idxP[k]
                X_i = if (K_X > 0 && i <= nrow(X_matrix)) X_matrix[i, ] else numeric(0)
                j = Individual_to_Cluster_Map$j_idx[Individual_to_Cluster_Map$i_idx == i][1]
                C_j = if (K_C > 0 && j <= nrow(C_matrix)) C_matrix[j, ] else numeric(0)

                Pc[k] = calculate_PCE_XC(t, d1, d0, I, params, X_i, C_j, lam_z, lam_mz, N_MC)
              }
              PCE_C = mean(Pc, na.rm = TRUE)

              out[[k_res]] = data.frame(
                Draw = original_draw, Estimand = nm, Stratum = str_nm,
                Rho_Scenario = rho_scn, Province = P, PCE = PCE_C,
                Delta = DELTA, Lambda_K = Ksc,
                Stratum_Probability = prob_str,
                Lambda_z = lam_z, Lambda_minus_z = lam_mz
              )
              k_res = k_res + 1
            }
          }
        }
      }
    }
  }
  do.call(rbind, out)
}

if (!is.null(cl)) { parallel::stopCluster(cl); foreach::registerDoSEQ() }

# ---------- summarize ----------
PCE_summary <- data.frame()
Probability_summary <- data.frame()
Lambda_summary <- data.frame()

if (!is.null(all_results) && nrow(all_results) > 0) {
  defs = dplyr::bind_rows(estimands_list, .id="Estimand") %>% dplyr::select(Estimand, Time = t, Duration, Type)

  PCE_summary = all_results %>%
    dplyr::filter(is.finite(PCE)) %>%
    dplyr::left_join(defs, by = "Estimand") %>%
    dplyr::group_by(Estimand, Stratum, Time, Duration, Type, Province, Delta, Lambda_K, Rho_Scenario) %>%
    dplyr::summarise(
      Mean   = mean(PCE),
      Median = stats::median(PCE),
      Q2.5   = stats::quantile(PCE, 0.025),
      Q97.5  = stats::quantile(PCE, 0.975),
      .groups = "drop"
    )

  Probability_summary = all_results %>%
    dplyr::filter(is.finite(Stratum_Probability), Stratum_Probability >= 0, Stratum_Probability <= 1) %>%
    dplyr::left_join(dplyr::distinct(defs, Estimand, Duration, Type), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Type, Stratum, Delta, Rho_Scenario, Stratum_Probability) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Duration, Type, Stratum, Delta, Rho_Scenario) %>%
    dplyr::summarise(
      Prob_Mean   = mean(Stratum_Probability),
      Prob_Median = stats::median(Stratum_Probability),
      Prob_Q2.5   = stats::quantile(Stratum_Probability, 0.025),
      Prob_Q97.5  = stats::quantile(Stratum_Probability, 0.975),
      .groups = "drop"
    )

  Lambda_summary = all_results %>%
    dplyr::left_join(dplyr::distinct(defs, Estimand, Duration), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Lambda_K, Lambda_z, Lambda_minus_z) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(cols = c(Lambda_z, Lambda_minus_z), names_to = "Lambda_Type", values_to = "Value")
}

save(all_results, PCE_summary, Probability_summary, Lambda_summary, estimands_list,
     file = CONFIG_PATHS$PCE_RESULTS)
