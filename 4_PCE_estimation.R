# 4_PCE_estimation.R — amplitude scaling for lambda; PCE posterior summaries

stopifnot(file.exists("Analysis_Config.RData"))
load("Analysis_Config.RData")
stopifnot(file.exists(CONFIG_PATHS$PROCESSED_DATA)); load(CONFIG_PATHS$PROCESSED_DATA)
stopifnot(file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES)); posterior_samples <- readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)
stopifnot(file.exists(CONFIG_PATHS$CALIBRATION_RESULTS)); calib <- readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)

stopifnot(!is.na(K_X), !is.na(K_C), !is.na(D_max), !is.na(T_max), !is.na(sn_sd), sn_sd>0)

DELTA_LIST <- CONFIG_PCE$DELTA_LIST
RHO_SCENARIOS <- CONFIG_PCE$RHO_SCENARIOS
LAMBDA_K_SCALING <- CONFIG_PCE$LAMBDA_K_SCALING

# Aux posteriors
stopifnot(!is.null(calib$Posteriors))
rho_star_posterior <- calib$Posteriors$rho_star
lambda_z_lower_posterior <- calib$Posteriors$lambda_z_lower
lambda_minus_z_lower_posterior <- calib$Posteriors$lambda_minus_z_lower
if (!is.matrix(lambda_z_lower_posterior))        lambda_z_lower_posterior        <- as.matrix(lambda_z_lower_posterior)
if (!is.matrix(lambda_minus_z_lower_posterior))  lambda_minus_z_lower_posterior  <- as.matrix(lambda_minus_z_lower_posterior)

valid_rho <- !is.na(rho_star_posterior) & rho_star_posterior > -0.999 & rho_star_posterior < 0.999
if (nrow(lambda_z_lower_posterior) == length(valid_rho)) {
  valid_lz  <- rowSums(is.na(lambda_z_lower_posterior))==0
  valid_lmz <- rowSums(is.na(lambda_minus_z_lower_posterior))==0
  valid <- valid_rho & valid_lz & valid_lmz
} else {
  m <- min(length(valid_rho), nrow(lambda_z_lower_posterior), nrow(lambda_minus_z_lower_posterior))
  valid <- valid_rho[1:m] & (rowSums(is.na(lambda_z_lower_posterior[1:m,,drop=FALSE]))==0) &
           (rowSums(is.na(lambda_minus_z_lower_posterior[1:m,,drop=FALSE]))==0)
}
stopifnot(sum(valid) > 0)
rho_star_posterior <- rho_star_posterior[valid]
lambda_z_lower_posterior       <- lambda_z_lower_posterior[valid,, drop=FALSE]
lambda_minus_z_lower_posterior <- lambda_minus_z_lower_posterior[valid,, drop=FALSE]
S_aux <- length(rho_star_posterior)

# Parallel setup
n_cores <- if (!is.null(CONFIG_STAN$CORES)) CONFIG_STAN$CORES else max(1, parallel::detectCores())
cl <- tryCatch(parallel::makeCluster(n_cores, outfile = "", type = "SOCK"), error=function(e) NULL)
if (!is.null(cl)) doSNOW::registerDoSNOW(cl) else foreach::registerDoSEQ()

# Utilities
expit <- function(x) { p <- 1/(1+exp(-x)); p[x>709] <- 1; p[x< -709] <- 0; p }
N_gh <- 20
gh <- statmod::gauss.quad(N_gh, kind="hermite")
gh_nodes <- gh$nodes
gh_weights <- gh$weights / sqrt(pi)

calc_eff <- function(v, b) {
  if (length(b)>0 && length(v)>0 && !all(is.na(v))) {
    vv <- as.vector(v); bb <- as.vector(b)
    if (length(vv)==length(bb)) return(sum(vv*bb))
  }
  0
}

map_params <- function(s, posterior, D_max, K_X, K_C) {
  ext <- function(obj, idx) {
    if (is.null(obj)) return(numeric(0))
    dm <- dim(obj)
    if (is.null(dm)) {
      if (idx<=length(obj)) obj[idx] else numeric(0)
    } else if (length(dm)==2) {
      if (idx<=dm[1]) obj[idx,] else numeric(0)
    } else if (length(dm)==3) {
      if (idx<=dm[1]) obj[idx,,] else numeric(0)
    } else numeric(0)
  }
  p <- list()
  p$eta1 <- ext(posterior$eta1, s)
  p$eta2 <- ext(posterior$eta2, s)
  p$gamma<- if (D_max>0) ext(posterior$gamma, s) else numeric(0)
  p$beta <- if (D_max>0) ext(posterior$beta, s) else numeric(0)
  p$omega1 <- if (K_X>0) ext(posterior$omega1, s) else numeric(0)
  p$psi1   <- if (K_X>0) ext(posterior$psi1, s)   else numeric(0)
  p$omega2 <- if (K_C>0) ext(posterior$omega2, s) else numeric(0)
  p$psi2   <- if (K_C>0) ext(posterior$psi2, s)   else numeric(0)
  p$psi3 <- ext(posterior$psi3, s)
  p$psi4 <- ext(posterior$psi4, s)
  Sa <- ext(posterior$Sigma_alpha, s)
  Sp <- ext(posterior$Sigma_phi, s)
  se <- ext(posterior$sigma_epsilon, s)
  if (length(p$eta1)==0 || length(se)==0 || !is.matrix(Sa) || !is.matrix(Sp)) return(NULL)
  p$Var_M <- Sa[1,1] + Sp[1,1] + se^2
  if (is.na(p$Var_M) || p$Var_M < 1e-9) p$Var_M <- 1e-9
  p$Cov_M_RE2 <- Sa[1,2] + Sp[1,2]
  p$Var_RE2   <- Sa[2,2] + Sp[2,2]
  p$Var_RE2_given_M <- p$Var_RE2 - (p$Cov_M_RE2^2 / p$Var_M)
  if (is.na(p$Var_RE2_given_M) || p$Var_RE2_given_M < 1e-9) p$Var_RE2_given_M <- 1e-9
  p$rho_sensitivity <- NA_real_
  p
}

get_hist <- function(t, d, p) {
  t_idx <- t - 1
  list(
    E_M_time = if (t_idx>=1 && t_idx<=length(p$eta1)) p$eta1[t_idx] else 0,
    E_M_duration = if (d>0 && d<=length(p$gamma)) p$gamma[d] else 0,
    E_Y_time = if (t_idx>=1 && t_idx<=length(p$eta2)) p$eta2[t_idx] else 0,
    E_Y_duration = if (d>0 && d<=length(p$beta)) p$beta[d] else 0
  )
}

E_Y_given_M_Z_X_C <- function(t, d, m, p, Xv, Cv) {
  h <- get_hist(t, d, p)
  X_eff <- calc_eff(Xv, p$psi1); C_eff <- calc_eff(Cv, p$psi2)
  mu <- h$E_Y_time + h$E_Y_duration + X_eff + C_eff + m*p$psi3
  if (d>0) mu <- mu + m*p$psi4
  Xm <- calc_eff(Xv, p$omega1); Cm <- calc_eff(Cv, p$omega2)
  E_M <- h$E_M_time + h$E_M_duration + Xm + Cm
  E_RE2 <- (p$Cov_M_RE2 / p$Var_M) * (m - E_M)
  sd_RE2 <- sqrt(p$Var_RE2_given_M)
  nodes <- E_RE2 + sqrt(2) * sd_RE2 * gh_nodes
  sum(gh_weights * expit(mu + nodes))
}

solve_Delta <- function(t, d, d_star, m_z, E_Y_obs, p, Xv, Cv, lambda_diff) {
  Xm <- calc_eff(Xv, p$omega1); Cm <- calc_eff(Cv, p$omega2)
  Ez   <- get_hist(t, d, p);       E_M_z   <- Ez$E_M_time + Ez$E_M_duration + Xm + Cm
  Ezs  <- get_hist(t, d_star, p);  E_M_z_s <- Ezs$E_M_time + Ezs$E_M_duration + Xm + Cm
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  if (is.na(rho) || abs(rho)>1) return(NA_real_)
  E_M_star_cond <- E_M_z_s + rho * (m_z - E_M_z)
  Var_star_cond <- (1 - rho^2) * Var_M
  if (Var_star_cond < 1e-9) Var_star_cond <- 1e-9
  nodes_M_star <- E_M_star_cond + sqrt(2*Var_star_cond) * gh_nodes

  MAX_ITER <- 50; TOL <- 1e-7
  Ey <- pmin(1-1e-9, pmax(1e-9, E_Y_obs))
  Delta <- log(Ey/(1-Ey))
  for (it in 1:MAX_ITER) {
    lp <- Delta + lambda_diff * nodes_M_star
    f  <- -E_Y_obs + sum(gh_weights * expit(lp))
    if (abs(f) < TOL) return(Delta)
    der <- sum(gh_weights * (expit(lp) * (1 - expit(lp))))
    if (abs(der) < 1e-9) break
    Delta <- Delta - f/der
  }
  Delta
}

PCE_XC <- function(t, d, d_star, I, p, Xv, Cv, lambda_z, lambda_mz, N_mc) {
  Xm <- calc_eff(Xv, p$omega1); Cm <- calc_eff(Cv, p$omega2)
  E1 <- get_hist(t, d, p);     Ez <- E1$E_M_time + E1$E_M_duration + Xm + Cm
  E2 <- get_hist(t, d_star, p); Ezs<- E2$E_M_time + E2$E_M_duration + Xm + Cm
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  if (is.na(rho) || abs(rho)>1 || is.na(lambda_z) || is.na(lambda_mz)) return(NA_real_)
  Mu <- c(Ez, Ezs); Sig <- matrix(c(Var_M, rho*Var_M, rho*Var_M, Var_M), 2)
  M <- tryCatch(mvtnorm::rmvnorm(N_mc, mean=Mu, sigma=Sig), error=function(e) NULL)
  if (is.null(M)) return(NA_real_)
  m_z <- M[,1]; m_zs <- M[,2]
  diff <- m_z - m_zs
  in_I <- (diff >= I[1]) & (diff < I[2])
  if (!any(in_I)) return(0)
  mz <- m_z[in_I]; mzs <- m_zs[in_I]
  EYo  <- vapply(seq_along(mz), function(k) E_Y_given_M_Z_X_C(t, d,    mz[k],  p, Xv, Cv), numeric(1))
  EYos <- vapply(seq_along(mzs), function(k) E_Y_given_M_Z_X_C(t, d_star, mzs[k], p, Xv, Cv), numeric(1))
  Dz  <- vapply(seq_along(mz),  function(k) solve_Delta(t, d, d_star, mz[k],  EYo[k],  p, Xv, Cv, lambda_z),  numeric(1))
  Dzs <- vapply(seq_along(mzs), function(k) solve_Delta(t, d_star, d, mzs[k], EYos[k], p, Xv, Cv, lambda_mz), numeric(1))
  if (any(is.na(Dz)) || any(is.na(Dzs))) return(NA_real_)
  mean(expit(Dz  + lambda_z  * mzs) - expit(Dzs + lambda_mz * mz))
}

Stratum_Prob <- function(I, d, d_star, p) {
  rho <- p$rho_sensitivity; Var_M <- p$Var_M
  if (is.na(rho) || abs(rho)>1) return(NA_real_)
  gd <- if (d>0 && d<=length(p$gamma)) p$gamma[d] else 0
  gs <- if (d_star>0 && d_star<=length(p$gamma)) p$gamma[d_star] else 0
  mu <- gd - gs; v <- 2*(1-rho)*Var_M
  if (v < 1e-9) return(if (mu>=I[1] && mu<I[2]) 1 else 0)
  pnorm(I[2], mean=mu, sd=sqrt(v)) - pnorm(I[1], mean=mu, sd=sqrt(v))
}

# Estimands (d vs 0)
estimands_list <- list()
for (L in CONFIG_PCE$DURATIONS) {
  if (L>0 && L<=D_max) {
    min_t <- L + 1
    if (min_t <= T_max) for (t in min_t:T_max) {
      estimands_list[[paste0("PCE_d",L,"_t",t)]] <- list(t=t, d_z=L, d_z_star=0, Duration=L,
                                                         Type=ifelse(L==1,"stPCE (d=1)", paste0("ltPCE (d=",L,")")))
    }
  }
}
stopifnot(length(estimands_list) > 0)

# Index subset
N_post <- length(posterior_samples$psi3)
S_sub  <- min(CONFIG_PCE$S_SUBSET, N_post)
set.seed(123)
idx_S <- sample(1:N_post, S_sub)

Province_Individual_Map <- lapply(Provinces_List, function(P) {
  Individual_to_Cluster_Map$i_idx[Individual_to_Cluster_Map$Province == P]
})
names(Province_Individual_Map) <- as.character(Provinces_List)

subset_draws <- function(x, idx) {
  if (is.null(x)) return(NULL)
  dm <- dim(x)
  if (is.matrix(x) && nrow(x) >= max(idx)) return(x[idx,,drop=FALSE])
  if (is.array(x) && length(dm)==3 && dm[1] >= max(idx)) return(x[idx,,,drop=FALSE])
  if ((is.vector(x) || (is.array(x) && length(dm)==1)) && length(x) >= max(idx)) return(x[idx])
  if (is.null(dm) && is.numeric(x) && length(x) >= max(idx)) return(x[idx])
  NULL
}
posterior_sub <- lapply(posterior_samples, subset_draws, idx = idx_S)
posterior_sub <- Filter(Negate(is.null), posterior_sub)
rm(posterior_samples); gc()

S <- S_sub; N_MC <- CONFIG_PCE$N_MC_INTEGRATION
N_indiv_target <- CONFIG_PCE$N_INDIV_SAMPLE
LARGE <- 1e9

all_results <- foreach::foreach(s_idx = 1:S, .combine = rbind,
                                .packages = c("statmod","mvtnorm","dplyr","triangle")) %dopar% {
  set.seed(s_idx + 12345)
  s <- s_idx
  orig_s <- idx_S[s_idx]
  params <- map_params(s, posterior_sub, D_max, K_X, K_C)
  if (is.null(params)) return(NULL)

  out <- list(); k_out <- 1
  s_aux_idx <- sample(1:S_aux, 1)
  rho_star_s <- rho_star_posterior[s_aux_idx]

  for (nm in names(estimands_list)) {
    est <- estimands_list[[nm]]
    t <- est$t; d <- est$d_z; d0 <- est$d_z_star; d_idx <- d

    Lz_cal <- if (d_idx>0 && d_idx<=ncol(lambda_z_lower_posterior)) lambda_z_lower_posterior[s_aux_idx, d_idx] else 0
    Lmz_cal<- if (d_idx>0 && d_idx<=ncol(lambda_minus_z_lower_posterior)) lambda_minus_z_lower_posterior[s_aux_idx, d_idx] else 0
    if (is.na(Lz_cal))  Lz_cal  <- 0
    if (is.na(Lmz_cal)) Lmz_cal <- 0
    if (is.na(rho_star_s) || abs(rho_star_s)>=1) rho_star_s <- 0.5

    for (RHO_SC in RHO_SCENARIOS) {
      if (RHO_SC=="Calibrated_Stochastic") params$rho_sensitivity <- rho_star_s
      else if (startsWith(RHO_SC,"Fixed_")) {
        rfix <- as.numeric(sub("Fixed_","",RHO_SC))
        if (is.na(rfix) || rfix < -1 || rfix > 1) next
        params$rho_sensitivity <- rfix
      } else next

      for (DELTA in DELTA_LIST) {
        delta_std <- DELTA / sn_sd
        I1 <- c(-delta_std,  delta_std)
        I2 <- c(-LARGE,     -delta_std)
        I3 <- c( delta_std,  LARGE)
        strata <- list(PCE1_Dissociative=I1, PCE2_Associative_Neg=I2, PCE3_Associative_Pos=I3)

        for (K in LAMBDA_K_SCALING) {
          lambda_z <- K * Lz_cal
          lambda_mz<- K * Lmz_cal

          for (stratum_name in names(strata)) {
            I <- strata[[stratum_name]]
            prob_I <- Stratum_Prob(I, d, d0, params)
            if (is.na(prob_I)) next

            for (P in names(Province_Individual_Map)) {
              idxP <- Province_Individual_Map[[P]]
              if (length(idxP)==0) next
              if (length(idxP) > N_indiv_target) idxP <- sample(idxP, N_indiv_target)
              PCE_i <- numeric(length(idxP))
              for (k in seq_along(idxP)) {
                i <- idxP[k]
                Xv <- if (K_X>0 && i <= nrow(X_matrix)) X_matrix[i,] else numeric(0)
                j  <- Individual_to_Cluster_Map$j_idx[Individual_to_Cluster_Map$i_idx==i][1]
                Cv <- if (K_C>0 && j <= nrow(C_matrix)) C_matrix[j,] else numeric(0)
                PCE_i[k] <- PCE_XC(t, d, d0, I, params, Xv, Cv, lambda_z, lambda_mz, N_MC)
              }
              out[[k_out]] <- data.frame(
                Draw=orig_s, Estimand=nm, Stratum=stratum_name,
                Rho_Scenario=RHO_SC, Province=P,
                PCE=mean(PCE_i, na.rm=TRUE),
                Delta=DELTA, Lambda_K=K,
                Stratum_Probability=prob_I,
                Lambda_z=lambda_z, Lambda_minus_z=lambda_mz
              )
              k_out <- k_out + 1
            }
          }
        }
      }
    }
  }
  if (length(out)==0) return(NULL)
  do.call(rbind, out)
}

if (!is.null(cl)) { parallel::stopCluster(cl); foreach::registerDoSEQ() }

PCE_summary <- data.frame()
Probability_summary <- data.frame()
Lambda_summary <- data.frame()

if (!is.null(all_results) && nrow(all_results)>0) {
  defs <- dplyr::bind_rows(estimands_list, .id="Estimand") |>
    dplyr::select(Estimand, Time=t, Duration, Type)

  PCE_summary <- all_results |>
    dplyr::filter(!is.na(PCE) & !is.nan(PCE)) |>
    dplyr::left_join(defs, by="Estimand") |>
    dplyr::group_by(Estimand, Stratum, Time, Duration, Type, Province, Delta, Lambda_K, Rho_Scenario) |>
    dplyr::summarise(
      Mean = mean(PCE), Median = stats::median(PCE),
      Q2.5 = stats::quantile(PCE,0.025), Q97.5 = stats::quantile(PCE,0.975),
      .groups="drop"
    )

  Probability_summary <- all_results |>
    dplyr::filter(!is.na(Stratum_Probability) & !is.nan(Stratum_Probability) &
                    Stratum_Probability>=0 & Stratum_Probability<=1) |>
    dplyr::left_join(defs |> dplyr::select(Estimand, Duration, Type) |> dplyr::distinct(), by="Estimand") |>
    dplyr::select(Draw, Duration, Type, Stratum, Delta, Rho_Scenario, Stratum_Probability) |>
    dplyr::distinct() |>
    dplyr::group_by(Duration, Type, Stratum, Delta, Rho_Scenario) |>
    dplyr::summarise(
      Prob_Mean = mean(Stratum_Probability, na.rm=TRUE),
      Prob_Median = stats::median(Stratum_Probability, na.rm=TRUE),
      Prob_Q2.5 = stats::quantile(Stratum_Probability, 0.025, na.rm=TRUE),
      Prob_Q97.5 = stats::quantile(Stratum_Probability, 0.975, na.rm=TRUE),
      .groups="drop"
    )

  Lambda_summary <- all_results |>
    dplyr::left_join(defs |> dplyr::select(Estimand, Duration) |> dplyr::distinct(), by="Estimand") |>
    dplyr::select(Draw, Duration, Lambda_K, Lambda_z, Lambda_minus_z) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(cols=c(Lambda_z, Lambda_minus_z), names_to="Lambda_Type", values_to="Value")
}

save(all_results, PCE_summary, Probability_summary, Lambda_summary, estimands_list, file = CONFIG_PATHS$PCE_RESULTS)
