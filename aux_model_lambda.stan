// aux_model_lambda.stan
// Bayesian Logistic Regression with Duration Interactions (Eq 9).
// Used to estimate duration-dependent lower bounds theta_2(d) and zeta_1(d).

data {
  int<lower=1> N;
  int<lower=1> D_max; // Number of durations (1..D_max)
  int<lower=0> K_X;
  int<lower=0> K_C;

  array[N] int<lower=0, upper=1> Y; // Outcome (Y_post or Y_pre)
  vector[N] M_post;
  vector[N] M_pre;
  matrix[N, K_X] X;
  matrix[N, K_C] C;
  // Index for duration (1 to D_max)
  array[N] int<lower=1, upper=D_max> Duration_idx;
}

parameters {
  // Duration-specific coefficients (D_max rows for each coefficient)
  // We parameterize directly using an indexing approach for efficiency.
  vector[D_max] Beta_0; // Intercepts theta_0(d) / zeta_0(d)
  vector[D_max] Beta_M_post; // theta_1(d) or zeta_1(d)
  vector[D_max] Beta_M_pre; // theta_2(d) or zeta_2(d)

  // Covariate effects (constant across d, as per paper description Eq 9)
  vector[K_X] Beta_X;
  vector[K_C] Beta_C;
}

model {
  // Priors (Weakly informative)
  Beta_0 ~ normal(0, 5);
  Beta_M_post ~ normal(0, 5);
  Beta_M_pre ~ normal(0, 5);
  if (K_X > 0) Beta_X ~ normal(0, 5);
  if (K_C > 0) Beta_C ~ normal(0, 5);

  // Likelihood (Eq 9)
  vector[N] eta;

  // Pre-calculate X/C effects
  vector[N] X_eff = (K_X > 0) ? X * Beta_X : rep_vector(0.0, N);
  vector[N] C_eff = (K_C > 0) ? C * Beta_C : rep_vector(0.0, N);

  for (n in 1:N) {
    int d = Duration_idx[n];
    // Linear predictor with duration-specific coefficients
    eta[n] = Beta_0[d] + Beta_M_post[d] * M_post[n] + Beta_M_pre[d] * M_pre[n];

    // Add covariates
    eta[n] += X_eff[n] + C_eff[n];
  }

  Y ~ bernoulli_logit(eta);
}

// The posterior samples for the lower bounds are extracted in R:
// For theta_2(d) (Model Y_post): Extract Beta_M_pre
// For zeta_1(d) (Model Y_pre): Extract Beta_M_post