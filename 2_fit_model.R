# 2_fit_model.R
# Fit the main Stan model and save posterior draws.

stopifnot(exists("CONFIG_PATHS"), exists("CONFIG_STAN"))
stopifnot(file.exists(CONFIG_PATHS$PROCESSED_DATA))
load(CONFIG_PATHS$PROCESSED_DATA)

options(mc.cores = CONFIG_STAN$CORES)
rstan_options(auto_write = TRUE)

stan_file <- CONFIG_PATHS$STAN_MODEL
stopifnot(file.exists(stan_file))

model <- stan_model(file = stan_file)

# Ensure X and C have valid dimensions even when K=0
if (stan_data$K_X == 0) stan_data$X <- matrix(0, nrow = stan_data$I, ncol = 0)
if (stan_data$K_C == 0) stan_data$C <- matrix(0, nrow = stan_data$J, ncol = 0)

fit <- sampling(model, data = stan_data,
                iter = CONFIG_STAN$ITER, warmup = CONFIG_STAN$WARMUP, chains = CONFIG_STAN$CHAINS,
                control = list(adapt_delta = CONFIG_STAN$ADAPT_DELTA, max_treedepth = CONFIG_STAN$MAX_TREEDEPTH))

posterior_samples <- rstan::extract(fit)
saveRDS(posterior_samples, file = CONFIG_PATHS$POSTERIOR_SAMPLES)
