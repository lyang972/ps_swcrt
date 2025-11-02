// aux_model_rho.stan
// Bayesian Multivariate Linear Regression for (M_pre, M_post).
// Implements Eq 7 and calculates marginalized rho* (Eq 8) in generated quantities.
// Assumes residual covariance is constant across durations (rho is duration-independent).

data {
  int<lower=1> N; // Total observations
  int<lower=1> D_max; // Max duration (for mean structure adjustment)
  int<lower=0> K_X; // Number of individual covariates
  int<lower=0> K_C; // Number of cluster covariates (dummies)

  matrix[N, 2] M; // [M_pre, M_post]
  matrix[N, K_X] X;
  matrix[N, K_C] C;
  array[N] int<lower=1, upper=D_max> Duration_idx;

  // Pre-calculated Var(X|C) - Required for marginalization (Paper p.18)
  int<lower=1> N_C_levels; // Number of unique C levels (e.g., 2 Provinces)
  // We pass the Var(X|C) for each level of C.
  // If K_X=0, this array will have dimensions [N_C_levels, 0, 0].
  array[N_C_levels] matrix[K_X, K_X] Var_X_given_C;
}

parameters {
  // Duration effects (Matrix: [D_max, 2]). Mean effect for each duration.
  matrix[D_max, 2] Beta_D;

  // Covariate effects (Matrix: [K_X, 2] for alpha_X and beta_X)
  matrix[K_X, 2] Beta_X;
  // Cluster covariate effects (Matrix: [K_C, 2])
  matrix[K_C, 2] Beta_C;

  // Residual covariance (constant across d)
  vector<lower=0>[2] sigma_epsilon;
  cholesky_factor_corr[2] L_Omega_epsilon;
}

transformed parameters {
  matrix[2, 2] L_Sigma_epsilon = diag_pre_multiply(sigma_epsilon, L_Omega_epsilon);
}

model {
  // Priors (Weakly informative)
  to_vector(Beta_D) ~ normal(0, 5);
  if (K_X > 0) to_vector(Beta_X) ~ normal(0, 5);
  if (K_C > 0) to_vector(Beta_C) ~ normal(0, 5);

  sigma_epsilon ~ exponential(1);
  L_Omega_epsilon ~ lkj_corr_cholesky(1);

  // Likelihood (Eq 7)
  // Pre-calculate X/C effects for efficiency
  matrix[N, 2] X_eff = (K_X > 0) ? X * Beta_X : rep_matrix(0.0, N, 2);
  matrix[N, 2] C_eff = (K_C > 0) ? C * Beta_C : rep_matrix(0.0, N, 2);

  for (n in 1:N) {
    vector[2] mu_M;
    int d = Duration_idx[n];

    // Mean structure depends on Duration, X, and C
    mu_M = Beta_D[d, ]' + X_eff[n, ]' + C_eff[n, ]';

    // Multivariate normal likelihood with constant residual covariance
    M[n, ] ~ multi_normal_cholesky(mu_M, L_Sigma_epsilon);
  }
}

generated quantities {
  matrix[2, 2] Sigma_epsilon = multiply_lower_tri_self_transpose(L_Sigma_epsilon);
  // Marginalized correlation conditional on C (rho_c*)
  vector[N_C_levels] rho_star_c;

  // Calculate rho* using Law of Total Covariance (Eq 8)
  if (K_X > 0) {
    // Extract alpha_X (M_pre) and beta_X (M_post)
    vector[K_X] alpha_X = Beta_X[, 1];
    vector[K_X] beta_X_vec = Beta_X[, 2]; // Renamed to avoid conflict

    for (c in 1:N_C_levels) {
      matrix[K_X, K_X] V_X_C = Var_X_given_C[c];

      // Cov(M_pre, M_post | C) = Cov(eps) + alpha_X' * Var(X|C) * beta_X
      real Cov_M_marginal = Sigma_epsilon[1, 2] + alpha_X' * V_X_C * beta_X_vec;
      // Var(M_pre | C) = Var(eps_pre) + alpha_X' * Var(X|C) * alpha_X
      real Var_M_pre_marginal = Sigma_epsilon[1, 1] + alpha_X' * V_X_C * alpha_X;
      // Var(M_post | C) = Var(eps_post) + beta_X' * Var(X|C) * beta_X
      real Var_M_post_marginal = Sigma_epsilon[2, 2] + beta_X_vec' * V_X_C * beta_X_vec;

      if (Var_M_pre_marginal > 1e-9 && Var_M_post_marginal > 1e-9) {
        rho_star_c[c] = Cov_M_marginal / sqrt(Var_M_pre_marginal * Var_M_post_marginal);
        // Clamp correlation to [-1, 1] for numerical stability
        rho_star_c[c] = fmax(-1.0, fmin(1.0, rho_star_c[c]));
      } else {
        rho_star_c[c] = 0.0; // Default if variance is near zero
      }
    }
  } else {
    // If K_X = 0, rho_star is just the residual correlation (independent of C)
    if (Sigma_epsilon[1, 1] > 1e-9 && Sigma_epsilon[2, 2] > 1e-9) {
        real rho_resid = Sigma_epsilon[1, 2] / sqrt(Sigma_epsilon[1, 1] * Sigma_epsilon[2, 2]);
        for (c in 1:N_C_levels) {
            rho_star_c[c] = fmax(-1.0, fmin(1.0, rho_resid));
        }
    } else {
        for (c in 1:N_C_levels) {
            rho_star_c[c] = 0.0;
        }
    }
  }
}