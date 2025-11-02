# 2_fit_model.R

# Load configuration and data (expected to be run from 0_driver.R)
if (!exists("CONFIG_PATHS") || !exists("CONFIG_STAN")) {
  message("Configuration not loaded. Please run via 0_driver.R.")
  if (sys.nframe() > 0) return(invisible(NULL))
}

if (!file.exists(CONFIG_PATHS$PROCESSED_DATA)) {
  message("Processed data not found.")
  if (sys.nframe() > 0) return(invisible(NULL))
}
load(CONFIG_PATHS$PROCESSED_DATA)

# Configure Stan
options(mc.cores = CONFIG_STAN$CORES)
rstan_options(auto_write = TRUE)

stan_file = CONFIG_PATHS$STAN_MODEL

if (!file.exists(stan_file)) {
  message(paste("Stan model file not found:", stan_file))
  if (sys.nframe() > 0) return(invisible(NULL))
}

# Compile the Stan model
message(paste("Compiling Stan model:", stan_file))
# Allow compilation failure to stop execution naturally
model = stan_model(file = stan_file)

# Prepare data for Stan (Ensure matrix dimensions if K=0)
# These checks are necessary for Stan compatibility
if (stan_data$K_X == 0 && (nrow(stan_data$X) != stan_data$I || ncol(stan_data$X) != 0)) {
  stan_data$X = matrix(0, nrow=stan_data$I, ncol=0)
}
if (stan_data$K_C == 0 && (nrow(stan_data$C) != stan_data$J || ncol(stan_data$C) != 0)) {
  stan_data$C = matrix(0, nrow=stan_data$J, ncol=0)
}

# Run the sampling
message(sprintf("Starting Stan sampling (Chains: %d, Iter: %d)... This may take time.",
                CONFIG_STAN$CHAINS, CONFIG_STAN$ITER))

# Allow sampling failure to stop execution naturally
fit = sampling(model, data = stan_data,
               iter = CONFIG_STAN$ITER, warmup = CONFIG_STAN$WARMUP, chains = CONFIG_STAN$CHAINS,
               control = list(adapt_delta = CONFIG_STAN$ADAPT_DELTA, max_treedepth = CONFIG_STAN$MAX_TREEDEPTH))

message("Stan sampling finished.")

# Basic diagnostics
summary_fit = summary(fit)$summary
Rhat_values = summary_fit[, "Rhat"]
if (any(Rhat_values > 1.1, na.rm = TRUE)) {
  message("Warning: Some Rhat values > 1.1. Model convergence may be poor.")
} else {
  message("Model convergence checks (Rhat) look good.")
}

# Extract and save posterior samples
posterior_samples = rstan::extract(fit)
saveRDS(posterior_samples, file = CONFIG_PATHS$POSTERIOR_SAMPLES)
message(paste("Posterior samples saved to", CONFIG_PATHS$POSTERIOR_SAMPLES))