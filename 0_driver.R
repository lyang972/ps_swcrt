# 0_driver.R
# Master script for Principal Stratification analysis.

# ==============================================================================
# 0. Setup and Configuration
# ==============================================================================

setwd('~/PROJECT_ps/')

# --- 0.1. Load Packages ---
required_packages <- c("dplyr", "tidyr", "readxl", "tibble", "forcats", "rstan",
                       "statmod", "mvtnorm", "triangle", "ggplot2",
                       "foreach", "parallel", "Matrix", "broom",
                       "RColorBrewer", "gridExtra", "grid", "gtable",
                       "doSNOW", "purrr")

# Quietly load packages
suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

# --- 0.2. File Paths and Configuration ---
CONFIG_PATHS = list(
  # !!! UPDATE THIS PATH TO YOUR DATA FILE !!!
  DATA_PATH = '../data/pmed.1002645.s007.xlsx',
  # Models
  STAN_MODEL = "model_duration.stan",
  STAN_AUX_RHO = "aux_model_rho.stan",
  STAN_AUX_LAMBDA = "aux_model_lambda.stan",
  # Outputs
  PROCESSED_DATA = "Processed_Data.RData",
  POSTERIOR_SAMPLES = "Posterior_Samples.rds",
  CALIBRATION_RESULTS = "Calibration_Results_Bayesian.rds",
  PCE_RESULTS = "PCE_Results.RData",
  ANALYSIS_CONFIG = "Analysis_Config.RData"
)

# --- 0.3. Analysis Flow Control ---
RUN_STEPS = list(
  PREPROCESS = F,
  MODEL_FIT = F,
  CALIBRATION = F,
  PCE_ESTIMATION = F
)

# --- 0.4. Stan Model Fitting Parameters ---
n_cores_available = parallel::detectCores()
if (is.na(n_cores_available) || n_cores_available < 1) n_cores_available = 1

# Main Model Configuration
CONFIG_STAN = list(
  ITER = 4000, WARMUP = 2000, CHAINS = 4,
  CORES = n_cores_available,
  ADAPT_DELTA = 0.95, MAX_TREEDEPTH = 15
)

# Auxiliary Models Configuration
CONFIG_STAN_AUX = list(
  ITER = 4000, WARMUP = 2000, CHAINS = 4,
  CORES = n_cores_available,
  ADAPT_DELTA = 0.90, MAX_TREEDEPTH = 12
)

# --- 0.5. PCE Configuration ---
CONFIG_PCE = list(
  DURATIONS = c(1, 2, 3),
  DELTA_LIST = c(0.5, 1.0, 1.5),
  # "Calibrated_Stochastic" uses the posterior of rho* from the Bayesian aux model
  RHO_SCENARIOS = c("Calibrated_Stochastic", "Fixed_0.8", "Fixed_0.9"),
  # Lambda Sensitivity Analysis (Amplitude Scaling k). k=1.0 is baseline.
  LAMBDA_K_SCALING = c(0.5, 1.0, 1.5, 2.0),
  S_SUBSET = 200,         # Subset of posterior draws
  N_INDIV_SAMPLE = 100,   # Individuals sampled for marginalization
  N_MC_INTEGRATION = 200  # Monte Carlo samples
)

# Primary Analysis Settings
PRIMARY_SETTINGS = list(
  DELTA = 0.5,
  RHO_SCENARIO = "Calibrated_Stochastic",
  LAMBDA_K = 1.0
)

# Load or Save Configuration
if (!any(unlist(RUN_STEPS))) {
  # Load existing config if analysis is not running (for plotting/reporting)
  if (file.exists(CONFIG_PATHS$ANALYSIS_CONFIG)) {
    # Load into a temporary environment to avoid overwriting RUN_STEPS
    config_env <- new.env()
    load(CONFIG_PATHS$ANALYSIS_CONFIG, envir = config_env)
    # Overwrite defaults selectively
    if (exists("CONFIG_PATHS", envir = config_env)) CONFIG_PATHS = config_env$CONFIG_PATHS
    if (exists("CONFIG_STAN", envir = config_env)) CONFIG_STAN = config_env$CONFIG_STAN
    if (exists("CONFIG_STAN_AUX", envir = config_env)) CONFIG_STAN_AUX = config_env$CONFIG_STAN_AUX
    if (exists("CONFIG_PCE", envir = config_env)) CONFIG_PCE = config_env$CONFIG_PCE
    if (exists("PRIMARY_SETTINGS", envir = config_env)) PRIMARY_SETTINGS = config_env$PRIMARY_SETTINGS
    message("Loaded existing analysis configuration.")
  }
} else {
  # If analysis is running, save the current configuration
  save(CONFIG_PATHS, CONFIG_STAN, CONFIG_STAN_AUX, CONFIG_PCE, PRIMARY_SETTINGS, RUN_STEPS, file = CONFIG_PATHS$ANALYSIS_CONFIG)
}

# ==============================================================================
# 1. Run Analysis Pipeline
# ==============================================================================

message(sprintf("\n--- Starting Analysis Pipeline (%s) ---\n", Sys.time()))

# Define the data path for preprocessing
data_path = CONFIG_PATHS$DATA_PATH

# --- Step 1: Preprocessing ---
if (RUN_STEPS$PREPROCESS) {
  message("\n--- Running Step 1: Preprocessing ---")
  source("1_preprocess.R")
}

# --- Step 2: Model Fitting (Main Model) ---
if (RUN_STEPS$MODEL_FIT) {
  message("\n--- Running Step 2: Model Fitting (Stan) ---")
  source("2_fit_model.R")
}

# --- Step 3: Calibration (Bayesian Auxiliary Models) ---
if (RUN_STEPS$CALIBRATION) {
  message("\n--- Running Step 3: Bayesian Calibration ---")
  # Assuming 3_calibration.R exists
  if (file.exists("3_calibration.R")) {
    source("3_calibration.R")
  }
}

# --- Step 4: PCE Estimation ---
if (RUN_STEPS$PCE_ESTIMATION) {
  message("\n--- Running Step 4: PCE Estimation ---")
  # Assuming 4_PCE_estimation.R exists
  if (file.exists("4_PCE_estimation.R")) {
    source("4_PCE_estimation.R")
  }
}

# ==============================================================================
# 2. Generate Reports and Plots
# ==============================================================================
message("\n--- Generating Reports and Plots ---")

# Define a clean theme for the tables
# Check for required packages for PDF reporting
if (requireNamespace("gridExtra", quietly = TRUE) && requireNamespace("grid", quietly = TRUE) && requireNamespace("gtable", quietly = TRUE)) {
  ttheme_header <- gridExtra::ttheme_default(
    core = list(fg_params=list(fontsize=10), bg_params = list(fill = "grey95")),
    colhead = list(fg_params=list(fontsize=11, fontface="bold"), bg_params = list(fill = "grey85"))
  )
  padding <- grid::unit(0.5, "line")
} else {
  ttheme_header = NULL
  message("Note: gridExtra/grid/gtable not available. PDF reports will be skipped.")
}


# Initialize variables
posterior_samples_report = NULL
calibration_results_report = NULL
Calibration_report_data = NULL

# ==============================================================================
# 2.1. Report: Sensitivity Parameter Calibration (Posterior Summary)
# ==============================================================================

# Load Calibration Results
if (file.exists(CONFIG_PATHS$CALIBRATION_RESULTS)) {
  calibration_results_report = readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)
}

if (!is.null(calibration_results_report) && !is.null(calibration_results_report$Posteriors)) {
  
  Posteriors = calibration_results_report$Posteriors
  
  # Helper function to summarize posterior draws (returns a data.frame row)
  summarize_posterior = function(draws, param_name, duration = NA) {
    if (is.null(draws) || length(draws) == 0) return(NULL)
    data.frame(
      Param = param_name,
      Duration = duration,
      Mean = mean(draws, na.rm = TRUE),
      Median = median(draws, na.rm = TRUE),
      SD = sd(draws, na.rm = TRUE),
      Q2.5 = quantile(draws, 0.025, na.rm = TRUE),
      Q97.5 = quantile(draws, 0.975, na.rm = TRUE)
    )
  }
  
  # --- Process Rho* (Duration-independent) ---
  summary_rho = summarize_posterior(Posteriors$rho_star, "Rho* (Calibrated)", "Independent")
  
  # --- Process Lambda (Duration-dependent) ---
  D_max_report_lz = if (is.matrix(Posteriors$lambda_z_lower)) ncol(Posteriors$lambda_z_lower) else 0
  D_max_report_lmz = if (is.matrix(Posteriors$lambda_minus_z_lower)) ncol(Posteriors$lambda_minus_z_lower) else 0
  D_max_report = max(D_max_report_lz, D_max_report_lmz)
  
  summary_lambda_z = data.frame()
  summary_lambda_minus_z = data.frame()
  
  if (D_max_report > 0) {
    for (d in 1:D_max_report) {
      d_str = as.character(d)
      
      # Lambda(z) Calibrated (Posterior of Theta_2(d))
      if (d <= D_max_report_lz) {
        lz_draws = Posteriors$lambda_z_lower[, d]
        summary_lambda_z = rbind(summary_lambda_z,
                                 summarize_posterior(lz_draws, "Lambda_z Calibrated (Theta_2(d))", d_str))
      }
      
      # Lambda(-z) Calibrated (Posterior of Zeta_1(d))
      if (d <= D_max_report_lmz) {
        lmz_draws = Posteriors$lambda_minus_z_lower[, d]
        summary_lambda_minus_z = rbind(summary_lambda_minus_z,
                                       summarize_posterior(lmz_draws, "Lambda_minus_z Calibrated (Zeta_1(d))", d_str))
      }
    }
  }
  
  # Combine results
  Calibration_report_data = bind_rows(summary_rho, summary_lambda_z, summary_lambda_minus_z)
  
  if (nrow(Calibration_report_data) > 0) {
    # Ensure correct sorting
    Calibration_report_data = Calibration_report_data %>%
      mutate(Duration_num = suppressWarnings(as.numeric(as.character(Duration)))) %>%
      mutate(Duration_Sort = ifelse(is.na(Duration_num), 0, Duration_num)) %>%
      arrange(Duration_Sort, Param) %>%
      select(-Duration_Sort, -Duration_num)
    
    # Save CSV
    write.csv(Calibration_report_data, "Report_2_1_Calibration_Summary.csv", row.names = FALSE)
    
    # Generate PDF
    if (!is.null(ttheme_header)) {
      Calib_formatted = Calibration_report_data %>%
        mutate(across(where(is.numeric), ~sprintf("%.4f", .x))) %>%
        mutate(`95% CrI` = paste0("[", Q2.5, ", ", Q97.5, "]")) %>%
        select(Duration, Parameter=Param, Mean, SD, `95% CrI`)
      
      table_grob_calib <- gridExtra::tableGrob(Calib_formatted, rows = NULL, theme = ttheme_header)
      title_text = "Report 2.1: Posterior Summary of Calibrated Sensitivity Parameters"
      title_calib <- grid::textGrob(title_text, gp=grid::gpar(fontsize=14, fontface="bold"))
      
      table_with_title_calib <- gtable::gtable_add_rows(
        table_grob_calib,
        heights = grid::grobHeight(title_calib) + padding * 2,
        pos = 0
      )
      table_with_title_calib <- gtable::gtable_add_grob(
        table_with_title_calib,
        title_calib,
        t = 1, l = 1, r = ncol(table_with_title_calib)
      )
      
      pdf_width = 11
      pdf_height = 3 + 0.3 * nrow(Calib_formatted)
      
      tryCatch({
        pdf("Report_2_1_Calibration_Summary.pdf", width = pdf_width, height = pdf_height)
        grid::grid.draw(table_with_title_calib)
        dev.off()
      }, error = function(e) {
        if (names(dev.cur()) != "null device") dev.off()
      })
    }
  }
}

# ==============================================================================
# 2.2. Report: Key Observed Data Model Parameters (修正后)
# ==============================================================================

# Load D_max and Posterior Samples
D_max = NULL
if (file.exists(CONFIG_PATHS$PROCESSED_DATA)) {
  # Load necessary variables (like D_max) safely
  processed_env <- new.env()
  load(CONFIG_PATHS$PROCESSED_DATA, envir = processed_env)
  if (exists("D_max", envir = processed_env)) D_max = processed_env$D_max
}

if (file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES)) {
  posterior_samples_report = readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)
}

# Simplified extraction assuming known structure of Stan output
if (!is.null(posterior_samples_report)) {
  
  # Helper function for summary (returns a structured data.frame row)
  create_summary_df_row = function(param_name, description, draws) {
    if (is.null(draws) || length(draws) == 0 || all(is.na(draws))) {
      return(data.frame(Parameter = param_name, Description = description, Mean=NA_real_, SD=NA_real_, Q2.5=NA_real_, Q97.5=NA_real_))
    }
    data.frame(
      Parameter = param_name,
      Description = description,
      Mean = mean(draws, na.rm = TRUE),
      SD = sd(draws, na.rm = TRUE),
      Q2.5 = quantile(draws, 0.025, na.rm = TRUE),
      Q97.5 = quantile(draws, 0.975, na.rm = TRUE)
    )
  }
  
  param_summary_list = list()
  
  # Duration effects (Beta and Gamma)
  if (!is.null(D_max) && D_max > 0) {
    for (d in 1:D_max) {
      # Check dimensions before access
      if (!is.null(posterior_samples_report$beta) && ncol(posterior_samples_report$beta) >= d && 
          !is.null(posterior_samples_report$gamma) && ncol(posterior_samples_report$gamma) >= d) {
        beta_draws = posterior_samples_report$beta[, d]
        gamma_draws = posterior_samples_report$gamma[, d]
        param_summary_list[[paste0("beta[", d, "]")]] = create_summary_df_row(paste0("beta[", d, "]"), "Beta_d (Y Duration Effect)", beta_draws)
        param_summary_list[[paste0("gamma[", d, "]")]] = create_summary_df_row(paste0("gamma[", d, "]"), "Gamma_d (M Duration Effect)", gamma_draws)
      }
    }
  }
  
  # M-Y Association
  param_summary_list[["psi3"]] = create_summary_df_row("psi3", "Psi_3 (M-Y Assoc., D=0)", posterior_samples_report$psi3)
  param_summary_list[["psi4"]] = create_summary_df_row("psi4", "Psi_4 (M-Y Assoc. Diff, D>0)", posterior_samples_report$psi4)
  
  # Key Variance Components
  # Check dimensions before access
  if (!is.null(posterior_samples_report$Sigma_alpha) && all(dim(posterior_samples_report$Sigma_alpha)[2:3] == c(2,2))) {
    Sigma_alpha_11_draws = posterior_samples_report$Sigma_alpha[, 1, 1]
    param_summary_list[["Sigma_alpha[1,1]"]] = create_summary_df_row("Sigma_alpha[1,1]", "Sigma^2_alpha1 (Cluster Var M)", Sigma_alpha_11_draws)
  }
  if (!is.null(posterior_samples_report$Sigma_phi) && all(dim(posterior_samples_report$Sigma_phi)[2:3] == c(2,2))) {
    Sigma_phi_11_draws = posterior_samples_report$Sigma_phi[, 1, 1]
    param_summary_list[["Sigma_phi[1,1]"]] = create_summary_df_row("Sigma_phi[1,1]", "Sigma^2_phi1 (Indiv Var M)", Sigma_phi_11_draws)
  }
  
  # Residual Variance
  if ("sigma_epsilon" %in% names(posterior_samples_report)) {
    sigma_epsilon_sq_draws = posterior_samples_report$sigma_epsilon^2
    param_summary_list[["sigma_epsilon^2"]] = create_summary_df_row("sigma_epsilon^2", "Sigma^2_epsilon (Residual Var M)", sigma_epsilon_sq_draws)
  }
  
  if (length(param_summary_list) > 0) {
    
    # Convert list of data frames to a single DataFrame
    Model_Param_Report = dplyr::bind_rows(param_summary_list)
    
    # Calculate Total Variance of M
    var_m_components = c("Sigma^2_alpha1 (Cluster Var M)", "Sigma^2_phi1 (Indiv Var M)", "Sigma^2_epsilon (Residual Var M)")
    
    Var_M_data = Model_Param_Report %>%
      filter(Description %in% var_m_components & !is.na(Mean) & !is.nan(Mean))
    
    # Calculate total variance if all components are present
    if (nrow(Var_M_data) == length(var_m_components)) {
      # Summarize means (Columns Mean, SD, Q2.5, Q97.5 now exist correctly)
      Var_M_total = Var_M_data %>%
        summarise(across(c(Mean, SD, Q2.5, Q97.5), ~sum(., na.rm=TRUE))) %>%
        mutate(Parameter = "Var(M)", Description = "Total Variance of M")
      
      Model_Param_Report = bind_rows(Model_Param_Report, Var_M_total)
    }
    
    # Save CSV
    write.csv(Model_Param_Report, "Report_2_2_Model_Parameters_Summary.csv", row.names = FALSE)
    
    # Generate PDF
    if (!is.null(ttheme_header)) {
      Model_Param_PDF_data = Model_Param_Report %>%
        # Format for display
        dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), NA_character_, sprintf("%.4f", .)))) %>%
        dplyr::mutate(`95% CrI` = paste0("[", Q2.5, ", ", Q97.5, "]")) %>%
        dplyr::select(Description, Mean, SD, `95% CrI`)
      
      table_grob_model <- gridExtra::tableGrob(Model_Param_PDF_data, rows = NULL, theme = ttheme_header)
      title_text_model = "Report 2.2: Posterior Summaries of Key Model Parameters"
      title_model <- grid::textGrob(title_text_model, gp=grid::gpar(fontsize=14, fontface="bold"))
      
      table_with_title_model <- gtable::gtable_add_rows(
        table_grob_model,
        heights = grid::grobHeight(title_model) + padding * 2,
        pos = 0
      )
      table_with_title_model <- gtable::gtable_add_grob(
        table_with_title_model,
        title_model,
        t = 1, l = 1, r = ncol(table_with_title_model)
      )
      
      pdf_width = 10
      pdf_height = 3 + 0.4 * nrow(Model_Param_PDF_data)
      
      tryCatch({
        pdf("Report_2_2_Model_Parameters_Summary.pdf", width = pdf_width, height = pdf_height)
        grid::grid.draw(table_with_title_model)
        dev.off()
      }, error = function(e) {
        if (names(dev.cur()) != "null device") dev.off()
      })
    }
  }
  
  rm(posterior_samples_report)
  gc()
}

# ==============================================================================
# 2.3. Report: Principal Stratum Probabilities
# ==============================================================================

PCE_summary = NULL; Probability_summary = NULL; Lambda_summary = NULL

# Load PCE results
if (file.exists(CONFIG_PATHS$PCE_RESULTS)) {
  # Load results into an environment
  results_env <- new.env()
  load(CONFIG_PATHS$PCE_RESULTS, envir = results_env)
  PCE_summary = results_env$PCE_summary
  Probability_summary = results_env$Probability_summary
  Lambda_summary = results_env$Lambda_summary
  
  if (!is.null(Probability_summary) && nrow(Probability_summary) > 0) {
    
    # Data Preparation for Reporting
    Prob_report_data = Probability_summary %>%
      dplyr::mutate(
        Stratum_Label = case_when(
          Stratum == "PCE1_Dissociative" ~ paste0("PCE1: Dissoc. (|dM| < ", sprintf("%.1f", Delta), ")"),
          Stratum == "PCE2_Associative_Neg" ~ paste0("PCE2: Assoc. Neg (dM < -", sprintf("%.1f", Delta), ")"),
          Stratum == "PCE3_Associative_Pos" ~ paste0("PCE3: Assoc. Pos (dM > ", sprintf("%.1f", Delta), ")"),
          TRUE ~ as.character(Stratum)
        )
      ) %>%
      dplyr::select(Duration, Type, Delta, Rho_Scenario, Stratum_Label, Prob_Mean, Prob_Q2.5, Prob_Q97.5) %>%
      dplyr::arrange(Duration, Delta, Rho_Scenario, Stratum_Label)
    
    if (nrow(Prob_report_data) > 0) {
      
      # Save CSV (Full results)
      write.csv(Prob_report_data, "Report_2_3_Stratum_Probabilities_Full.csv", row.names = FALSE)
      
      # Generate PDF (Primary Settings)
      PRIMARY_DELTA = PRIMARY_SETTINGS$DELTA
      PRIMARY_RHO = PRIMARY_SETTINGS$RHO_SCENARIO
      
      Prob_primary = Prob_report_data %>%
        filter(Delta == PRIMARY_DELTA & Rho_Scenario == PRIMARY_RHO)
      
      if (nrow(Prob_primary) > 0 && !is.null(ttheme_header)) {
        
        Prob_primary_formatted = Prob_primary %>%
          mutate(across(where(is.numeric), ~sprintf("%.4f", .x))) %>%
          mutate(`95% CrI` = paste0("[", Prob_Q2.5, ", ", Prob_Q97.5, "]")) %>%
          select(Duration, Type, Stratum=Stratum_Label, Mean=Prob_Mean, `95% CrI`)
        
        table_grob_prob <- gridExtra::tableGrob(Prob_primary_formatted, rows = NULL, theme = ttheme_header)
        title_text = sprintf("Report 2.3: Principal Stratum Probabilities (Primary: Delta=%.1f, Rho=%s)", PRIMARY_DELTA, PRIMARY_RHO)
        title_prob <- grid::textGrob(title_text, gp=grid::gpar(fontsize=14, fontface="bold"))
        
        table_with_title_prob <- gtable::gtable_add_rows(
          table_grob_prob,
          heights = grid::grobHeight(title_prob) + padding * 2,
          pos = 0
        )
        table_with_title_prob <- gtable::gtable_add_grob(
          table_with_title_prob,
          title_prob,
          t = 1, l = 1, r = ncol(table_with_title_prob)
        )
        
        pdf_width = 11.5
        pdf_height = 3 + 0.3 * nrow(Prob_primary)
        
        tryCatch({
          pdf("Report_2_3_Stratum_Probabilities_Primary.pdf", width = pdf_width, height = pdf_height)
          grid::grid.draw(table_with_title_prob)
          dev.off()
        }, error = function(e) {
          if (names(dev.cur()) != "null device") dev.off()
        })
      }
    }
  }
}

# ==============================================================================
# 2.4. Generate Plots
# ==============================================================================

# Note: The plotting section remains necessarily complex to visualize the multi-dimensional 
# sensitivity analysis results effectively.

if (!is.null(PCE_summary) && nrow(PCE_summary) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  
  # ==========================================================================
  # Define Palettes and Labels
  # ==========================================================================
  
  # Strata Definitions
  strata_levels = c("PCE1 (Dissoc.)", "PCE2 (Assoc. Neg)", "PCE3 (Assoc. Pos)")
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    strata_colors = RColorBrewer::brewer.pal(3, "Dark2")
  } else {
    strata_colors = c("#1B9E77", "#D95F02", "#7570B3") # Fallback
  }
  names(strata_colors) = strata_levels
  
  # Plotmath Labels (Required for subscripts like PCE_D)
  # 1. Mapping for Legends (expression objects)
  strata_labels_legend <- setNames(
    list(expression("PCE"["D"]), expression("PCE"["A-"]), expression("PCE"["A+"])),
    strata_levels
  )
  # 2. Mapping for Facets (strings for label_parsed)
  strata_labels_facet_map <- setNames(
    c("\"PCE\"[\"D\"]", "\"PCE\"[\"A-\"]", "\"PCE\"[\"A+\"]"),
    strata_levels
  )
  strata_labeller_parsed <- ggplot2::as_labeller(
    strata_labels_facet_map,
    default = ggplot2::label_parsed
  )
  
  # Delta Values (Blues)
  delta_values = sort(unique(PCE_summary$Delta))
  n_delta = length(delta_values)
  if (n_delta > 0) {
    # Simplified sequential color selection
    delta_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[5:9])(n_delta)
    names(delta_colors) = as.character(delta_values)
  } else {
    delta_colors = c()
  }
  
  
  # Lambda (k) Scenarios (Purples, highlighting baseline)
  k_values = sort(unique(PCE_summary$Lambda_K))
  n_k = length(k_values)
  if (n_k > 0) {
    k_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples")[5:9])(n_k)
    names(k_colors) = as.character(k_values)
    # Highlight k=1.0 (Baseline) using the green defined for PCE_D
    is_baseline = sapply(k_values, function(x) isTRUE(all.equal(x, 1.0)))
    if (any(is_baseline)) {
      k_colors[as.character(k_values[is_baseline])] = strata_colors[1] 
    }
  } else {
    k_colors = c()
  }
  
  
  # Rho Scenarios
  rho_scenario_levels_all = unique(PCE_summary$Rho_Scenario)
  mean_calibrated_rho_val = NA
  
  # Attempt to get mean calibrated rho from Report 2.1 data if available
  if ("Calibrated_Stochastic" %in% rho_scenario_levels_all && exists("Calibration_report_data")) {
    # Ensure Mean column is numeric before filtering/accessing
    if (is.numeric(Calibration_report_data$Mean)) {
      calib_rho_data = Calibration_report_data %>% filter(Param == "Rho* (Calibrated)")
      if (nrow(calib_rho_data) > 0) mean_calibrated_rho_val = calib_rho_data$Mean[1]
    }
  }
  
  # Create mapping for Rho Scenario labels
  rho_labels_map = list()
  for (scenario in rho_scenario_levels_all) {
    if (scenario == "Calibrated_Stochastic") {
      if (!is.na(mean_calibrated_rho_val)) {
        rho_labels_map[[scenario]] = sprintf("Calibrated (Mean Rho* approx %.3f)", mean_calibrated_rho_val)
      } else {
        rho_labels_map[[scenario]] = "Calibrated (Stochastic)"
      }
    } else if (startsWith(scenario, "Fixed_")) {
      rho_val = sub("Fixed_", "", scenario)
      rho_labels_map[[scenario]] = sprintf("Fixed (Rho = %s)", rho_val)
    } else {
      rho_labels_map[[scenario]] = scenario
    }
  }
  
  # Order and color Rho scenarios
  rho_calibrated = grep("Calibrated", rho_scenario_levels_all, value = TRUE)
  rho_fixed = sort(grep("Fixed_", rho_scenario_levels_all, value = TRUE))
  rho_order = c(rho_calibrated, rho_fixed)
  n_rho = length(rho_order)
  
  if (n_rho > 0) {
    # Use a distinct palette (e.g., Set1 or similar)
    rho_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_rho)
    names(rho_colors) = rho_order
  } else {
    rho_colors = c()
  }
  
  # ==========================================================================
  # Data Preparation
  # ==========================================================================
  
  PCE_plot_data = PCE_summary %>%
    dplyr::mutate(
      # Use consistent short labels defined in strata_levels
      Stratum_Short = factor(dplyr::case_when(
        Stratum == "PCE1_Dissociative" ~ "PCE1 (Dissoc.)",
        Stratum == "PCE2_Associative_Neg" ~ "PCE2 (Assoc. Neg)",
        Stratum == "PCE3_Associative_Pos" ~ "PCE3 (Assoc. Pos)"
      ), levels = strata_levels)
    )
  
  PCE_plot_data = PCE_plot_data %>%
    dplyr::filter(!is.na(Type) & !is.na(Duration)) %>%
    dplyr::mutate(Duration = as.numeric(as.character(Duration)))
  
  # Apply ordering
  if (nrow(PCE_plot_data) > 0) {
    PCE_plot_data = PCE_plot_data %>%
      dplyr::mutate(Type = factor(Type, levels = unique(Type[order(Duration)]))) %>%
      dplyr::mutate(Province = factor(Province))
    
    if (n_rho > 0) PCE_plot_data$Rho_Scenario = factor(PCE_plot_data$Rho_Scenario, levels = rho_order)
    if (n_k > 0) PCE_plot_data$Lambda_K_Factor = factor(PCE_plot_data$Lambda_K, levels = k_values)
  }
  
  n_provinces = length(unique(PCE_plot_data$Province))
  if (n_provinces == 0) n_provinces = 1
  
  PRIMARY_DELTA_VAL = PRIMARY_SETTINGS$DELTA
  PRIMARY_RHO_VAL = PRIMARY_SETTINGS$RHO_SCENARIO
  PRIMARY_LAMBDA_K_VAL = PRIMARY_SETTINGS$LAMBDA_K
  
  # ==========================================================================
  # Base Theme and Device
  # ==========================================================================
  
  theme_science = ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing.x = ggplot2::unit(1.5, "lines"),
      panel.spacing.y = ggplot2::unit(1, "lines"),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      plot.title = ggplot2::element_blank()
    )
  
  # Use Cairo PDF if available for potentially better rendering
  device_to_use <- "pdf"
  if (capabilities("cairo") && exists("cairo_pdf", where = "package:grDevices")) {
    device_to_use <- grDevices::cairo_pdf
  }
  
  # ==========================================================================
  # Plot 1: Comparison of Strata (Primary Analysis Setting) - Figure 1
  # ==========================================================================
  
  PCE_primary = PCE_plot_data %>%
    dplyr::filter(Lambda_K == PRIMARY_LAMBDA_K_VAL &
                    Delta == PRIMARY_DELTA_VAL &
                    Rho_Scenario == PRIMARY_RHO_VAL)
  
  if (nrow(PCE_primary) > 0) {
    p1_strata = ggplot2::ggplot(PCE_primary, ggplot2::aes(x = Time, y = Mean, color = Stratum_Short)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 2) +
      ggplot2::geom_line(ggplot2::aes(group = Stratum_Short), position = ggplot2::position_dodge(width = 0.6), linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.4, position = ggplot2::position_dodge(width = 0.6)) +
      ggplot2::facet_grid(Province ~ Type, scales = "free_x") +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI",
                    color = "Principal Stratum") +
      # Use the plotmath labels for the legend
      ggplot2::scale_color_manual(values = strata_colors, labels = strata_labels_legend) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      theme_science
    
    n_durations = length(unique(PCE_plot_data$Duration))
    plot_width_p1 = 4 + 2.5 * max(1, n_durations)
    plot_height_p1 = 4 * n_provinces + 1.5
    
    ggplot2::ggsave("Plot1_PCE_Strata_Comparison_Primary.pdf", plot = p1_strata, width = plot_width_p1, height = plot_height_p1, device = device_to_use)
  }
  
  # ==========================================================================
  # Plot 2: Sensitivity to Rho (ρ) - Figure 2
  # ==========================================================================
  
  PCE_rho_sensitivity = PCE_plot_data %>%
    dplyr::filter(Type == "stPCE (d=1)" &
                    Delta == PRIMARY_DELTA_VAL &
                    Lambda_K == PRIMARY_LAMBDA_K_VAL)
  
  if (nrow(PCE_rho_sensitivity) > 0 && length(unique(PCE_rho_sensitivity$Rho_Scenario)) > 1) {
    
    # Ensure labels match the levels present in the filtered data
    levels_present <- levels(droplevels(PCE_rho_sensitivity$Rho_Scenario))
    current_rho_labels = sapply(levels_present, function(s) rho_labels_map[[s]])
    
    p2_rho = ggplot2::ggplot(PCE_rho_sensitivity, ggplot2::aes(x = Time, y = Mean, color = Rho_Scenario)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      # Use the plotmath labeller for the facets
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI",
                    color = "Rho Scenario") +
      ggplot2::scale_color_manual(values = rho_colors, labels = current_rho_labels, drop = TRUE) +
      theme_science
    
    n_strata = length(unique(PCE_plot_data$Stratum_Short))
    plot_width_p2 = 2 + 3.5 * max(1, n_strata)
    plot_height_p2 = 4 * n_provinces + 1.8
    
    ggplot2::ggsave("Plot2_PCE_Rho_Sensitivity.pdf", plot = p2_rho, width = plot_width_p2, height = plot_height_p2, device = device_to_use)
  }
  
  # ==========================================================================
  # Plot 3: Sensitivity to Lambda (λ) Amplitude Scaling (k) - Figure 3
  # ==========================================================================
  
  PCE_lambda_k_sensitivity = PCE_plot_data %>%
    dplyr::filter(Type == "stPCE (d=1)" &
                    Delta == PRIMARY_DELTA_VAL &
                    Rho_Scenario == PRIMARY_RHO_VAL)
  
  if (nrow(PCE_lambda_k_sensitivity) > 0 && n_k > 1) {
    
    # Create labels for K, highlighting the baseline
    k_labels = paste0("k=", k_values)
    is_baseline_label = sapply(k_values, function(x) isTRUE(all.equal(x, 1.0)))
    if (any(is_baseline_label)) k_labels[is_baseline_label] = "k=1.0 (Baseline)"
    names(k_labels) = as.character(k_values)
    
    p4_lambda_k = ggplot2::ggplot(PCE_lambda_k_sensitivity, ggplot2::aes(x = Time, y = Mean, color = Lambda_K_Factor)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      # Use the plotmath labeller for the facets
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI",
                    color = "Lambda Scaling (k)") +
      ggplot2::scale_color_manual(values = k_colors, labels = k_labels, drop = TRUE) +
      theme_science
    
    n_strata = length(unique(PCE_plot_data$Stratum_Short))
    plot_width_p4 = 2 + 3.5 * max(1, n_strata)
    plot_height_p4 = 4 * n_provinces + 1.8
    
    # Note: Saving as Plot 4 based on original file structure.
    ggplot2::ggsave("Plot4_PCE_Lambda_K_Sensitivity.pdf", plot = p4_lambda_k, width = plot_width_p4, height = plot_height_p4, device = device_to_use)
  }
  
  # ==========================================================================
  # Plot 4: Sensitivity to Delta (δ) - Figure 4
  # ==========================================================================
  
  PCE_delta_sensitivity = PCE_plot_data %>%
    dplyr::filter(Type == "stPCE (d=1)" &
                    Rho_Scenario == PRIMARY_RHO_VAL &
                    Lambda_K == PRIMARY_LAMBDA_K_VAL) %>%
    dplyr::mutate(Delta_Factor = factor(Delta))
  
  if (nrow(PCE_delta_sensitivity) > 0 && length(unique(PCE_delta_sensitivity$Delta)) > 1) {
    
    p3_delta = ggplot2::ggplot(PCE_delta_sensitivity, ggplot2::aes(x = Time, y = Mean, color = Delta_Factor)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      # Use the plotmath labeller for the facets
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI",
                    color = "Delta Value") +
      ggplot2::scale_color_manual(values = delta_colors, drop = TRUE) +
      theme_science
    
    n_strata = length(unique(PCE_plot_data$Stratum_Short))
    plot_width_p3 = 2 + 3.5 * max(1, n_strata)
    plot_height_p3 = 4 * n_provinces + 1.5
    
    # Note: Saving as Plot 3 based on original file structure.
    ggplot2::ggsave("Plot3_PCE_Delta_Sensitivity.pdf", plot = p3_delta, width = plot_width_p3, height = plot_height_p3, device = device_to_use)
  }
  
  
  # ==========================================================================
  # Plot 5: Visualization of Lambda Distributions (Supplementary)
  # ==========================================================================
  
  if (exists("Lambda_summary") && !is.null(Lambda_summary) && nrow(Lambda_summary) > 0) {
    
    Lambda_plot_data = Lambda_summary %>%
      dplyr::filter(!is.na(Duration)) %>%
      dplyr::mutate(Duration = as.numeric(as.character(Duration)))
    
    if (nrow(Lambda_plot_data) > 0) {
      sorted_durations <- sort(unique(Lambda_plot_data$Duration))
      duration_levels <- paste0("Duration d=", sorted_durations)
      
      Lambda_plot_data = Lambda_plot_data %>%
        # Use plotmath compatible labels: lambda(z)
        dplyr::mutate(Lambda_Type_Label = ifelse(Lambda_Type == "Lambda_z", "lambda(z)", "lambda(-z)"),
                      Duration_Label = factor(paste0("Duration d=", Duration), levels = duration_levels))
      
      # Prepare K factors and labels
      k_values_plot5 = sort(unique(Lambda_plot_data$Lambda_K))
      Lambda_plot_data$Lambda_K_Factor = factor(Lambda_plot_data$Lambda_K, levels = k_values_plot5)
      
      k_labels_plot5 = paste0("k=", k_values_plot5)
      is_baseline_label_p5 = sapply(k_values_plot5, function(x) isTRUE(all.equal(x, 1.0)))
      if (any(is_baseline_label_p5)) k_labels_plot5[is_baseline_label_p5] = "k=1.0 (Baseline)"
      names(k_labels_plot5) = as.character(k_values_plot5)
      
      p5_lambda_dist = ggplot2::ggplot(Lambda_plot_data, ggplot2::aes(x = Value, fill = Lambda_K_Factor, color = Lambda_K_Factor)) +
        ggplot2::geom_density(alpha = 0.4) +
        # Use label_parsed for the lambda(z) facet labels
        ggplot2::facet_grid(Lambda_Type_Label ~ Duration_Label, scales = "free_y", labeller = ggplot2::label_parsed) +
        ggplot2::labs(x = "Lambda Value", y = "Density", fill = "Scaling (k)", color = "Scaling (k)") +
        # Ensure colors match those used in Plot 4
        ggplot2::scale_fill_manual(values = k_colors, labels = k_labels_plot5, drop = TRUE) +
        ggplot2::scale_color_manual(values = k_colors, labels = k_labels_plot5, drop = TRUE) +
        theme_science
      
      n_durations = length(unique(Lambda_plot_data$Duration))
      plot_width_p5 = 3 + 3 * max(1, n_durations)
      plot_height_p5 = 7.0
      
      ggplot2::ggsave("Plot5_Lambda_K_Distributions.pdf", plot = p5_lambda_dist, width = plot_width_p5, height = plot_height_p5, device = device_to_use)
    }
  }
  
  message("\nAnalysis finished. Reports and Plots generated.")
}