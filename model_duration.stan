// model_duration.stan
// Implements Equations 4 and 5 with duration-specific effects AND cluster covariates C_j.

data {
  int<lower=1> N;             // Total observations
  int<lower=1> I;             // Individuals
  int<lower=1> J;             // Clusters
  int<lower=1> T_max;         // Max time period
  int<lower=0> D_max;         // Max duration
  int<lower=0> K_X;           // Number of individual covariates
  int<lower=0> K_C;           // Number of cluster covariates

  // Data
  vector[N] M;
  array[N] int<lower=0, upper=1> Y;

  // Indices
  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;
  // Time index starts at 2 (as per preprocessing)
  array[N] int<lower=2, upper=T_max> time_idx;
  array[N] int<lower=0> Duration;

  // Covariates
  matrix[I, K_X] X;
  matrix[J, K_C] C;
}

parameters {
  // Time effects (t=2..T_max, so T_max-1 elements)
  vector[T_max-1] eta1; vector[T_max-1] eta2;

  // Duration effects (D=1..D_max)
  vector[D_max] gamma; vector[D_max] beta;

  // Individual Covariate effects
  vector[K_X] omega1; vector[K_X] psi1;

  // Cluster Covariate effects (Province)
  vector[K_C] omega2; vector[K_C] psi2;

  // M-Y association
  real psi3; // Association when D=0
  real psi4; // Difference in association when D>0

  // Random effects (standardized for Non-Centered Parameterization - NCP)
  matrix[2, J] z_alpha; // Cluster level RE
  matrix[2, I] z_phi;   // Individual level RE

  // Variance components
  real<lower=0> sigma_epsilon;
  // Cluster RE variance and correlation
  vector<lower=0>[2] sigma_alpha;
  cholesky_factor_corr[2] L_Omega_alpha;
  // Individual RE variance and correlation
  vector<lower=0>[2] sigma_phi;
  cholesky_factor_corr[2] L_Omega_phi;
}

transformed parameters {
  matrix[J, 2] alpha;
  matrix[I, 2] phi;
  // Derived Covariance Matrices (for posterior analysis)
  matrix[2, 2] Sigma_alpha;
  matrix[2, 2] Sigma_phi;

  // Non-centered parameterization (NCP) for efficient sampling
  // alpha = L_Sigma_alpha * z_alpha
  // We use diag_pre_multiply to construct the Cholesky factor of the covariance matrix
  matrix[2, 2] L_Sigma_alpha = diag_pre_multiply(sigma_alpha, L_Omega_alpha);
  matrix[2, 2] L_Sigma_phi = diag_pre_multiply(sigma_phi, L_Omega_phi);

  alpha = (L_Sigma_alpha * z_alpha)';
  phi = (L_Sigma_phi * z_phi)';
  
  // Reconstruct full covariance matrices (Sigma = L * L')
  // This makes Sigma_alpha and Sigma_phi available in the posterior samples in R
  Sigma_alpha = multiply_lower_tri_self_transpose(L_Sigma_alpha);
  Sigma_phi = multiply_lower_tri_self_transpose(L_Sigma_phi);
}

model {
  // Priors (Based on Appendix D)
  // Weakly informative Normal(0, 10) on fixed effects. SD = sqrt(10) approx 3.162
  real prior_sd = 3.162;

  target += normal_lpdf(eta1 | 0, prior_sd);
  target += normal_lpdf(eta2 | 0, prior_sd);
  target += normal_lpdf(gamma | 0, prior_sd);
  target += normal_lpdf(beta | 0, prior_sd);

  // Conditional priors for covariates
  if (K_X > 0) {
    target += normal_lpdf(omega1 | 0, prior_sd);
    target += normal_lpdf(psi1 | 0, prior_sd);
  }
  if (K_C > 0) {
    target += normal_lpdf(omega2 | 0, prior_sd);
    target += normal_lpdf(psi2 | 0, prior_sd);
  }

  target += normal_lpdf(psi3 | 0, prior_sd);
  target += normal_lpdf(psi4 | 0, prior_sd);

  // Variance components: Exponential(1)
  target += exponential_lpdf(sigma_epsilon | 1);
  target += exponential_lpdf(sigma_alpha | 1);
  target += exponential_lpdf(sigma_phi | 1);

  // Correlations: LKJ(1) (Uniform over correlation matrices)
  target += lkj_corr_cholesky_lpdf(L_Omega_alpha | 1);
  target += lkj_corr_cholesky_lpdf(L_Omega_phi | 1);

  // NCP priors (Standard Normal)
  target += std_normal_lpdf(to_vector(z_alpha));
  target += std_normal_lpdf(to_vector(z_phi));

  // Likelihood
  {
    vector[N] mu_M;
    vector[N] eta_Y; // Linear predictor for Y

    // Pre-calculate covariate effects for efficiency
    vector[I] X_omega1 = (K_X > 0) ? X * omega1 : rep_vector(0.0, I);
    vector[I] X_psi1 = (K_X > 0) ? X * psi1 : rep_vector(0.0, I);
    vector[J] C_omega2 = (K_C > 0) ? C * omega2 : rep_vector(0.0, J);
    vector[J] C_psi2 = (K_C > 0) ? C * psi2 : rep_vector(0.0, J);

    for (n in 1:N) {
      int i = ii[n];
      int j = jj[n];
      // Adjust time index for accessing eta (time_idx 2 maps to index 1)
      int t_idx = time_idx[n] - 1;
      int d = Duration[n];

      // Model for M (Eq 4)
      mu_M[n] = eta1[t_idx] + X_omega1[i] + C_omega2[j];
      if (d > 0) {
        // Ensure index safety
        if (d <= D_max) {
           mu_M[n] += gamma[d];
        }
      }
      mu_M[n] += alpha[j, 1] + phi[i, 1];

      // Model for Y (Eq 5)
      eta_Y[n] = eta2[t_idx] + X_psi1[i] + C_psi2[j] + M[n] * psi3;
      if (d > 0) {
        if (d <= D_max) {
           eta_Y[n] += beta[d];
        }
        eta_Y[n] += M[n] * psi4;
      }
      eta_Y[n] += alpha[j, 2] + phi[i, 2];
    }

    // M follows a normal distribution
    target += normal_lpdf(M | mu_M, sigma_epsilon);
    // Y follows a Bernoulli distribution with logit link
    target += bernoulli_logit_lpmf(Y | eta_Y);
  }
}