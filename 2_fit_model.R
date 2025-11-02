# 2_fit_model.R — streamlined Stan fit for main model

if (!"rstan" %in% loadedNamespaces()) library(rstan)

if (!exists("CONFIG_PATHS") || !exists("CONFIG_STAN")) {
  if (file.exists("Analysis_Config.RData")) load("Analysis_Config.RData") else stop("Configuration missing.")
}
stopifnot(file.exists(CONFIG_PATHS$PROCESSED_DATA))
load(CONFIG_PATHS$PROCESSED_DATA)

options(mc.cores = CONFIG_STAN$CORES)
rstan_options(auto_write = TRUE)

stan_file <- CONFIG_PATHS$STAN_MODEL
stopifnot(file.exists(stan_file))

# Ensure X/C have correct shapes when K=0
if (stan_data$K_X == 0 && (nrow(stan_data$X) != stan_data$I || ncol(stan_data$X) != 0)) {
  stan_data$X <- matrix(0, nrow=stan_data$I, ncol=0)
}
if (stan_data$K_C == 0 && (nrow(stan_data$C) != stan_data$J || ncol(stan_data$C) != 0)) {
  stan_data$C <- matrix(0, nrow=stan_data$J, ncol=0)
}

model <- stan_model(file = stan_file)
fit <- sampling(
  model, data = stan_data,
  iter = CONFIG_STAN$ITER, warmup = CONFIG_STAN$WARMUP, chains = CONFIG_STAN$CHAINS,
  control = list(adapt_delta = CONFIG_STAN$ADAPT_DELTA, max_treedepth = CONFIG_STAN$MAX_TREEDEPTH)
)

summary_fit <- summary(fit)$summary
if (any(summary_fit[,"Rhat"] > 1.1, na.rm=TRUE)) message("Warning: some Rhat > 1.1")

posterior_samples <- rstan::extract(fit)
saveRDS(posterior_samples, file = CONFIG_PATHS$POSTERIOR_SAMPLES)
