# 4_PCE_estimation.R
# Principal Causal Effect estimation, including lambda scaling (k) and rho scenarios.

stopifnot(exists("CONFIG_PATHS"), exists("CONFIG_PCE"))
stopifnot(file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES))
stopifnot(file.exists(CONFIG_PATHS$CALIBRATION_RESULTS))
posterior_samples <- readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)
calibration_results_bayesian <- readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)

stopifnot(exists("K_X"), exists("K_C"), exists("T_max"), exists("sn_sd"))

DELTA_LIST <- CONFIG_PCE$DELTA_LIST
RHO_SCENARIOS <- CONFIG_PCE$RHO_SCENARIOS
LAMBDA_K_SCALING <- CONFIG_PCE$LAMBDA_K_SCALING

PostAux <- calibration_results_bayesian$Posteriors
rho_star_posterior <- PostAux$rho_star
lambda_z_lower_posterior <- as.matrix(PostAux$lambda_z_lower)
lambda_minus_z_lower_posterior <- as.matrix(PostAux$lambda_minus_z_lower)

valid_rho <- !is.na(rho_star_posterior) & (rho_star_posterior > -0.999) & (rho_star_posterior < 0.999)
S_aux_raw <- length(valid_rho)
stopifnot(nrow(lambda_z_lower_posterior) == S_aux_raw, nrow(lambda_minus_z_lower_posterior) == S_aux_raw)
valid_lz  <- rowSums(is.na(lambda_z_lower_posterior)) == 0
valid_lmz <- rowSums(is.na(lambda_minus_z_lower_posterior)) == 0
valid_draws <- valid_rho & valid_lz & valid_lmz
stopifnot(sum(valid_draws) > 0)

rho_star_posterior        <- rho_star_posterior[valid_draws]
lambda_z_lower_posterior  <- lambda_z_lower_posterior[valid_draws, , drop = FALSE]
lambda_minus_z_lower_posterior <- lambda_minus_z_lower_posterior[valid_draws, , drop = FALSE]
S_aux <- length(rho_star_posterior)

# Parallel backend
n_cores_to_use <- if (exists("CONFIG_STAN") && !is.null(CONFIG_STAN$CORES)) CONFIG_STAN$CORES else max(1, parallel::detectCores())
cl <- tryCatch(parallel::makeCluster(n_cores_to_use, type = "SOCK"), error = function(e) NULL)
if (!is.null(cl)) doSNOW::registerDoSNOW(cl) else foreach::registerDoSEQ()

# Utilities
expit <- function(x) { p <- 1/(1+exp(-x)); p[x>709] <- 1; p[x< -709] <- 0; p }
N_gh <- 20; gh <- statmod::gauss.quad(N_gh, kind = "hermite")
gh_nodes <- gh$nodes; gh_weights <- gh$weights / sqrt(pi)
calc_eff <- function(v, b) if (length(b)>0 && length(v)==length(b)) sum(as.vector(v)*as.vector(b)) else 0

# Map 1 draw of the main posterior to a convenient list of parameters
map_params <- function(s, posterior, D_max_val, K_X_val, K_C_val) {
  D_max <- as.integer(D_max_val); K_X <- as.integer(K_X_val); K_C <- as.integer(K_C_val)
  extr <- function(obj, idx) {
    if (is.null(obj)) return(NULL)
    d <- dim(obj)
    if (is.null(d) || length(d)==1) return(obj[idx])
    if (length(d)==2) return(obj[idx, ])
    if (length(d)==3) return(obj[idx, , ])
    NULL
  }
  p <- list()
  p$eta1 <- extr(posterior$eta1, s); p$eta2 <- extr(posterior$eta2, s)
  p$gamma <- if (D_max>0) extr(posterior$gamma, s) else numeric(0)
  p$beta  <- if (D_max>0) extr(posterior$beta,  s) else numeric(0)
  p$omega1 <- if (K_X>0) extr(posterior$omega1, s) else numeric(0)
  p$psi1   <- if (K_X>0) extr(posterior$psi1,   s) else numeric(0)
  p$omega2 <- if (K_C>0) extr(posterior$omega2, s) else numeric(0)
  p$psi2   <- if (K_C>0) extr(posterior$psi2,   s) else numeric(0)
  p$psi3 <- extr(posterior$psi3, s); p$psi4 <- extr(posterior$psi4, s)
  Sa <- extr(posterior$Sigma_alpha, s); Sp <- extr(posterior$Sigma_phi, s); se <- extr(posterior$sigma_epsilon, s)
  if (is.null(Sa) || is.null(Sp) || is.null(se)) return(NULL)
  p$Var_M <- Sa[1,1] + Sp[1,1] + se^2
  p$Cov_M_RE2 <- Sa[1,2] + Sp[1,2]
  p$Var_RE2 <- Sa[2,2] + Sp[2,2]
  if (is.na(p$Var_M) || p$Var_M < 1e-9) p$Var_M <- 1e-9
  p$Var_RE2_given_M <- p$Var_RE2 - (p$Cov_M_RE2^2 / p$Var_M)
  if (is.na(p$Var_RE2_given_M) || p$Var_RE2_given_M < 1e-9) p$Var_RE2_given_M <- 1e-9
  p$rho_sensitivity <- NA_real_
  p
}

# Time/duration effects helper
hist_eff <- function(t, d_z, p) {
  t_idx <- t-1
  list(
    E_M_time = if (t_idx>=1 && t_idx<=length(p$eta1)) p$eta1[t_idx] else 0,
    E_M_duration = if (d_z>0 && d_z<=length(p$gamma)) p$gamma[d_z] else 0,
    E_Y_time = if (t_idx>=1 && t_idx<=length(p$eta2)) p$eta2[t_idx] else 0,
    E_Y_duration = if (d_z>0 && d_z<=length(p$beta)) p$beta[d_z] else 0
  )
}

# E[Y(z) | M(z), X, C] via Hermite quadrature
E_Y_given_MZXC <- function(t, d_z, m, p, Xv, Cv) {
  h <- hist_eff(t, d_z, p)
  X_effY <- calc_eff(Xv, p$psi1); C_effY <- calc_eff(Cv, p$psi2)
  muY <- h$E_Y_time + h$E_Y_duration + X_effY + C_effY + m*p$psi3
  if (d_z>0) muY <- muY + m*p$psi4
  X_effM <- calc_eff(Xv, p$omega1); C_effM <- calc_eff(Cv, p$omega2)
  E_M_z <- h$E_M_time + h$E_M_duration + X_effM + C_effM
  E_RE2 <- (p$Cov_M_RE2 / p$Var_M) * (m - E_M_z)
  sd_RE2 <- sqrt(p$Var_RE2_given_M)
  nodes <- E_RE2 + sqrt(2)*sd_RE2*gh_nodes
  sum(gh_weights * expit(muY + nodes))
}

# Solve for Delta in the convolution equation (Newton)
solve_Delta <- function(t, d_z, d_z_star, m_z, E_Y_obs, p, Xv, Cv, lambda_diff) {
  X_effM <- calc_eff(Xv, p$omega1); C_effM <- calc_eff(Cv, p$omega2)
  h_z  <- hist_eff(t, d_z, p);      E_M_z  <- h_z$E_M_time  + h_z$E_M_duration  + X_effM + C_effM
  h_zs <- hist_eff(t, d_z_star, p); E_M_zs <- h_zs$E_M_time + h_zs$E_M_duration + X_effM + C_effM
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  E_mstar <- E_M_zs + rho * (m_z - E_M_z)
  Var_mst <- (1 - rho^2) * Var_M; if (Var_mst < 1e-9) Var_mst <- 1e-9
  sd_mst <- sqrt(Var_mst)
  nodes <- E_mstar + sqrt(2)*sd_mst*gh_nodes

  MAX_ITER <- 50; TOL <- 1e-7
  Ey <- min(max(E_Y_obs, 1e-9), 1-1e-9)
  Delta <- log(Ey/(1-Ey))
  for (it in 1:MAX_ITER) {
    lp <- Delta + lambda_diff * nodes; g <- expit(lp)
    f <- -E_Y_obs + sum(gh_weights * g)
    if (abs(f) < TOL) return(Delta)
    fp <- sum(gh_weights * g * (1-g))
    if (abs(fp) < 1e-9) break
    Delta <- Delta - f/fp
  }
  Delta
}

# PCE(X,C) via Monte Carlo over (M(z), M(z*))
calc_PCE_XC <- function(t, d_z, d_z_star, I, p, Xv, Cv, lambda_z, lambda_mz, N_mc) {
  X_effM <- calc_eff(Xv, p$omega1); C_effM <- calc_eff(Cv, p$omega2)
  h_z  <- hist_eff(t, d_z, p);      E_M_z  <- h_z$E_M_time  + h_z$E_M_duration  + X_effM + C_effM
  h_zs <- hist_eff(t, d_z_star, p); E_M_zs <- h_zs$E_M_time + h_zs$E_M_duration + X_effM + C_effM
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  mu <- c(E_M_z, E_M_zs); Sigma <- matrix(c(Var_M, rho*Var_M, rho*Var_M, Var_M), 2)
  M <- tryCatch(mvtnorm::rmvnorm(N_mc, mean = mu, sigma = Sigma), error = function(e) NULL)
  if (is.null(M)) return(NA_real_)
  mz <- M[,1]; mzs <- M[,2]; diff <- mz - mzs
  in_I <- diff >= I[1] & diff < I[2]; if (!any(in_I)) return(0)
  mz <- mz[in_I]; mzs <- mzs[in_I]; nI <- length(mz)
  Dz  <- numeric(nI); Dzs <- numeric(nI)
  for (k in 1:nI) {
    Ey_z  <- E_Y_given_MZXC(t, d_z,      mz[k],  p, Xv, Cv)
    Ey_zs <- E_Y_given_MZXC(t, d_z_star, mzs[k], p, Xv, Cv)
    Dz[k]  <- solve_Delta(t, d_z,      d_z_star, mz[k],  Ey_z,  p, Xv, Cv, lambda_z)
    Dzs[k] <- solve_Delta(t, d_z_star, d_z,      mzs[k], Ey_zs, p, Xv, Cv, lambda_mz)
  }
  Ey1 <- expit(Dz  + lambda_z  * mzs)
  Ey0 <- expit(Dzs + lambda_mz * mz)
  mean(Ey1 - Ey0, na.rm = TRUE)
}

# Probability of belonging to stratum I
prob_stratum <- function(I, d_z, d_z_star, p) {
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  g1 <- if (d_z>0 && d_z<=length(p$gamma)) p$gamma[d_z] else 0
  g0 <- if (d_z_star>0 && d_z_star<=length(p$gamma)) p$gamma[d_z_star] else 0
  mu <- g1 - g0
  v  <- 2*(1-rho)*Var_M
  if (v < 1e-9) return(ifelse(mu >= I[1] && mu < I[2], 1.0, 0.0))
  sd <- sqrt(v)
  pnorm(I[2], mean = mu, sd = sd) - pnorm(I[1], mean = mu, sd = sd)
}

# Estimands: compare duration L vs 0 at times t >= L+1
estimands_list <- list()
for (L in CONFIG_PCE$DURATIONS) {
  if (L > 0 && L <= D_max) {
    min_time <- L + 1
    if (min_time <= T_max) {
      for (t in min_time:T_max) {
        Type <- ifelse(L==1, "stPCE (d=1)", paste0("ltPCE (d=",L,")"))
        estimands_list[[paste0("PCE_d",L,"_t",t)]] <- list(t=t, d_z=L, d_z_star=0, Duration=L, Type=Type)
      }
    }
  }
}

N_posterior <- length(posterior_samples$psi3)
S_subset <- min(CONFIG_PCE$S_SUBSET, N_posterior)
N_ind_smpl <- CONFIG_PCE$N_INDIV_SAMPLE
N_MC <- CONFIG_PCE$N_MC_INTEGRATION
set.seed(123); idx_S <- sample(1:N_posterior, S_subset)

# Province -> individual indices
Province_Individual_Map <- list()
for (P in Provinces_List) {
  idxP <- Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == P]
  if (length(idxP) > 0) Province_Individual_Map[[as.character(P)]] <- idxP
}

# Subset posterior samples to reduce memory/IO in parallel workers
subset_draws <- function(param, idx) {
  if (is.null(param)) return(NULL)
  if (is.matrix(param)) return(param[idx, , drop = FALSE])
  if (is.array(param) && length(dim(param))==3) return(param[idx, , , drop = FALSE])
  if (is.numeric(param) && length(param) >= max(idx)) return(param[idx])
  NULL
}
post_sub <- lapply(posterior_samples, subset_draws, idx = idx_S)
post_sub <- Filter(Negate(is.null), post_sub)
rm(posterior_samples); gc()

LARGE_NUM <- 1e9

# Main loop over posterior draws
all_results <- foreach::foreach(s_idx = 1:S_subset, .combine = rbind,
                                .packages = c("statmod","mvtnorm","dplyr")) %dopar% {
  set.seed(s_idx + 12345)
  tryCatch({
    s <- s_idx; orig <- idx_S[s_idx]
    p <- map_params(s, post_sub, D_max, K_X, K_C); if (is.null(p)) return(NULL)
    results_s <- list(); idx <- 1
    # Sample one calibrated set from aux posteriors
    s_aux_idx <- sample(1:S_aux, 1)
    rho_s <- rho_star_posterior[s_aux_idx]

    for (nm in names(estimands_list)) {
      est <- estimands_list[[nm]]; t <- est$t; d_z <- est$d_z; d0 <- est$d_z_star; d_idx <- d_z
      Lz  <- if (d_idx>0 && d_idx<=ncol(lambda_z_lower_posterior)) lambda_z_lower_posterior[s_aux_idx, d_idx] else 0
      Lmz <- if (d_idx>0 && d_idx<=ncol(lambda_minus_z_lower_posterior)) lambda_minus_z_lower_posterior[s_aux_idx, d_idx] else 0
      p$rho_sensitivity <- rho_s

      for (RHO_SC in RHO_SCENARIOS) {
        rho_val <- if (RHO_SC=="Calibrated_Stochastic") rho_s else if (startsWith(RHO_SC,"Fixed_")) as.numeric(sub("Fixed_","",RHO_SC)) else NA
        if (is.na(rho_val)) next
        p$rho_sensitivity <- rho_val

        for (DELTA in DELTA_LIST) {
          delta_std <- DELTA / sn_sd
          I1 <- c(-delta_std,  delta_std)   # Dissociative
          I2 <- c(-LARGE_NUM, -delta_std)   # Assoc. Negative
          I3 <- c( delta_std,  LARGE_NUM)   # Assoc. Positive
          stratum_list <- list(PCE1_Dissociative=I1, PCE2_Associative_Neg=I2, PCE3_Associative_Pos=I3)

          for (KFAC in LAMBDA_K_SCALING) {
            lam_z  <- KFAC * Lz
            lam_mz <- KFAC * Lmz

            for (str_nm in names(stratum_list)) {
              I <- stratum_list[[str_nm]]
              prob_I <- prob_stratum(I, d_z, d0, p)
              for (P in names(Province_Individual_Map)) {
                idxP <- Province_Individual_Map[[P]]; NP <- length(idxP)
                if (NP == 0) next
                idx_sample <- if (NP > N_ind_smpl) sample(idxP, N_ind_smpl) else idxP
                PCE_i <- numeric(length(idx_sample))
                for (k in seq_along(idx_sample)) {
                  i <- idx_sample[k]
                  Xv <- if (K_X > 0 && i <= nrow(X_matrix)) X_matrix[i, ] else numeric(0)
                  j  <- Individual_to_Cluster_Map$j_idx[Individual_to_Cluster_Map$i_idx == i][1]
                  Cv <- if (K_C > 0 && !is.na(j) && j <= nrow(C_matrix)) C_matrix[j, ] else numeric(0)
                  PCE_i[k] <- calc_PCE_XC(t, d_z, d0, I, p, Xv, Cv, lam_z, lam_mz, N_MC)
                }
                PCE_C <- mean(PCE_i, na.rm = TRUE)
                results_s[[idx]] <- data.frame(
                  Draw=orig, Estimand=nm, Stratum=str_nm,
                  Rho_Scenario=RHO_SC, Province=P, PCE=PCE_C,
                  Delta=DELTA, Lambda_K=KFAC,
                  Stratum_Probability=prob_I,
                  Lambda_z=lam_z, Lambda_minus_z=lam_mz
                )
                idx <- idx + 1
              }
            }
          }
        }
      }
    }
    do.call(rbind, results_s)
  }, error = function(e) NULL)
}

if (!is.null(cl)) { parallel::stopCluster(cl); foreach::registerDoSEQ() }

# Summaries
PCE_summary <- data.frame(); Probability_summary <- data.frame(); Lambda_summary <- data.frame()
if (!is.null(all_results) && nrow(all_results)>0) {
  est_defs <- dplyr::bind_rows(estimands_list, .id="Estimand") %>% dplyr::select(Estimand, Time=t, Duration, Type)
  PCE_summary <- all_results %>%
    dplyr::filter(!is.na(PCE) & !is.nan(PCE)) %>%
    dplyr::left_join(est_defs, by = "Estimand") %>%
    dplyr::group_by(Estimand, Stratum, Time, Duration, Type, Province, Delta, Lambda_K, Rho_Scenario) %>%
    dplyr::summarise(Mean=mean(PCE), Median=stats::median(PCE),
                     Q2.5=stats::quantile(PCE,0.025), Q97.5=stats::quantile(PCE,0.975), .groups='drop')

  Probability_summary <- all_results %>%
    dplyr::filter(!is.na(Stratum_Probability)) %>%
    dplyr::left_join(est_defs %>% dplyr::select(Estimand, Duration, Type) %>% dplyr::distinct(), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Type, Stratum, Delta, Rho_Scenario, Stratum_Probability) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Duration, Type, Stratum, Delta, Rho_Scenario) %>%
    dplyr::summarise(Prob_Mean = mean(Stratum_Probability, na.rm = TRUE),
                     Prob_Median = stats::median(Stratum_Probability, na.rm = TRUE),
                     Prob_Q2.5 = stats::quantile(Stratum_Probability, 0.025, na.rm = TRUE),
                     Prob_Q97.5 = stats::quantile(Stratum_Probability, 0.975, na.rm = TRUE),
                     .groups = 'drop')

  Lambda_summary <- all_results %>%
    dplyr::left_join(est_defs %>% dplyr::select(Estimand, Duration) %>% dplyr::distinct(), by = "Estimand") %>%
    dplyr::select(Draw, Duration, Lambda_K, Lambda_z, Lambda_minus_z) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(cols = c(Lambda_z, Lambda_minus_z), names_to = "Lambda_Type", values_to = "Value")
}

save(all_results, PCE_summary, Probability_summary, Lambda_summary, estimands_list, file = CONFIG_PATHS$PCE_RESULTS)
