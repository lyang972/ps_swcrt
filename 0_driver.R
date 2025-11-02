# 0_driver.R
# Master driver: run the full pipeline or source individual steps as needed.

setwd('~/PROJECT_ps/')

# Dependencies
required_packages <- c(
  "dplyr","tidyr","readxl","tibble","forcats","rstan",
  "statmod","mvtnorm","triangle","ggplot2","foreach","parallel",
  "Matrix","broom","RColorBrewer","gridExtra","grid","gtable",
  "doSNOW","purrr"
)
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

# Paths and configuration
CONFIG_PATHS <- list(
  DATA_PATH = '../data/pmed.1002645.s007.xlsx',
  STAN_MODEL = "model_duration.stan",
  STAN_AUX_RHO = "aux_model_rho.stan",
  STAN_AUX_LAMBDA = "aux_model_lambda.stan",
  PROCESSED_DATA = "Processed_Data.RData",
  POSTERIOR_SAMPLES = "Posterior_Samples.rds",
  CALIBRATION_RESULTS = "Calibration_Results_Bayesian.rds",
  PCE_RESULTS = "PCE_Results.RData",
  ANALYSIS_CONFIG = "Analysis_Config.RData"
)

# Which steps to run from this script (set TRUE to run)
RUN_STEPS <- list(
  PREPROCESS = FALSE,
  MODEL_FIT = FALSE,
  CALIBRATION = FALSE,
  PCE_ESTIMATION = FALSE
)

n_cores_available <- max(1, parallel::detectCores())

# Stan configs
CONFIG_STAN <- list(
  ITER = 4000, WARMUP = 2000, CHAINS = 4, CORES = n_cores_available,
  ADAPT_DELTA = 0.95, MAX_TREEDEPTH = 15
)
CONFIG_STAN_AUX <- list(
  ITER = 4000, WARMUP = 2000, CHAINS = 4, CORES = n_cores_available,
  ADAPT_DELTA = 0.90, MAX_TREEDEPTH = 12
)

# PCE configs
CONFIG_PCE <- list(
  DURATIONS = c(1,2,3),
  DELTA_LIST = c(0.5,1.0,1.5),
  RHO_SCENARIOS = c("Calibrated_Stochastic","Fixed_0.8","Fixed_0.9"),
  LAMBDA_K_SCALING = c(0.5,1.0,1.5,2.0),
  S_SUBSET = 200, N_INDIV_SAMPLE = 100, N_MC_INTEGRATION = 200
)

# Primary display settings used for reports/plots
PRIMARY_SETTINGS <- list(DELTA = 0.5, RHO_SCENARIO = "Calibrated_Stochastic", LAMBDA_K = 1.0)

# Load/save analysis configuration
if (!any(unlist(RUN_STEPS))) {
  if (file.exists(CONFIG_PATHS$ANALYSIS_CONFIG)) {
    cfg <- new.env(); load(CONFIG_PATHS$ANALYSIS_CONFIG, envir = cfg)
    for (nm in c("CONFIG_PATHS","CONFIG_STAN","CONFIG_STAN_AUX","CONFIG_PCE","PRIMARY_SETTINGS"))
      if (exists(nm, envir = cfg)) assign(nm, get(nm, envir = cfg))
  }
} else {
  save(CONFIG_PATHS, CONFIG_STAN, CONFIG_STAN_AUX, CONFIG_PCE, PRIMARY_SETTINGS, RUN_STEPS,
       file = CONFIG_PATHS$ANALYSIS_CONFIG)
}

# Run pipeline
data_path <- CONFIG_PATHS$DATA_PATH
if (RUN_STEPS$PREPROCESS)     source("1_process.R")
if (RUN_STEPS$MODEL_FIT)      source("2_fit_model.R")
if (RUN_STEPS$CALIBRATION)    source("3_calibration.R")
if (RUN_STEPS$PCE_ESTIMATION) source("4_PCE_estimation.R")

# -------------------------- Reports and Plots --------------------------

posterior_samples_report <- NULL
calibration_results_report <- NULL
Calibration_report_data <- NULL

# 2.1 Posterior summary of calibrated sensitivity parameters (rho*, lambda lower bounds)
if (file.exists(CONFIG_PATHS$CALIBRATION_RESULTS)) {
  calibration_results_report <- readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)
}
if (!is.null(calibration_results_report) && !is.null(calibration_results_report$Posteriors)) {
  Posteriors <- calibration_results_report$Posteriors
  summarize_posterior <- function(draws, param_name, duration = NA) {
    if (is.null(draws) || length(draws) == 0) return(NULL)
    data.frame(
      Param = param_name, Duration = duration,
      Mean = mean(draws, na.rm = TRUE),
      Median = median(draws, na.rm = TRUE),
      SD = sd(draws, na.rm = TRUE),
      Q2.5 = quantile(draws, 0.025, na.rm = TRUE),
      Q97.5 = quantile(draws, 0.975, na.rm = TRUE)
    )
  }
  summary_rho <- summarize_posterior(Posteriors$rho_star, "Rho* (Calibrated)", "Independent")
  D_max_report_lz <- if (is.matrix(Posteriors$lambda_z_lower)) ncol(Posteriors$lambda_z_lower) else 0
  D_max_report_lmz <- if (is.matrix(Posteriors$lambda_minus_z_lower)) ncol(Posteriors$lambda_minus_z_lower) else 0
  D_max_report <- max(D_max_report_lz, D_max_report_lmz)
  summary_lambda_z <- data.frame(); summary_lambda_minus_z <- data.frame()
  if (D_max_report > 0) {
    for (d in 1:D_max_report) {
      if (d <= D_max_report_lz) {
        summary_lambda_z <- rbind(summary_lambda_z,
          summarize_posterior(Posteriors$lambda_z_lower[, d], "Lambda_z Calibrated (Theta_2(d))", as.character(d)))
      }
      if (d <= D_max_report_lmz) {
        summary_lambda_minus_z <- rbind(summary_lambda_minus_z,
          summarize_posterior(Posteriors$lambda_minus_z_lower[, d], "Lambda_minus_z Calibrated (Zeta_1(d))", as.character(d)))
      }
    }
  }
  Calibration_report_data <- dplyr::bind_rows(summary_rho, summary_lambda_z, summary_lambda_minus_z)
  if (nrow(Calibration_report_data) > 0) {
    Calibration_report_data <- Calibration_report_data %>%
      dplyr::mutate(Duration_num = suppressWarnings(as.numeric(as.character(Duration))),
                    Duration_Sort = ifelse(is.na(Duration_num), 0, Duration_num)) %>%
      dplyr::arrange(Duration_Sort, Param) %>%
      dplyr::select(-Duration_Sort, -Duration_num)
    write.csv(Calibration_report_data, "Report_2_1_Calibration_Summary.csv", row.names = FALSE)

    if (requireNamespace("gridExtra", quietly = TRUE) &&
        requireNamespace("grid", quietly = TRUE) &&
        requireNamespace("gtable", quietly = TRUE)) {

      ttheme_header <- gridExtra::ttheme_default(
        core = list(fg_params=list(fontsize=10), bg_params = list(fill = "grey95")),
        colhead = list(fg_params=list(fontsize=11, fontface="bold"), bg_params = list(fill = "grey85"))
      )
      padding <- grid::unit(0.5, "line")
      Calib_formatted <- Calibration_report_data %>%
        dplyr::mutate(across(where(is.numeric), ~sprintf("%.4f", .x))) %>%
        dplyr::mutate(`95% CrI` = paste0("[", Q2.5, ", ", Q97.5, "]")) %>%
        dplyr::select(Duration, Parameter=Param, Mean, SD, `95% CrI`)

      g <- gridExtra::tableGrob(Calib_formatted, rows = NULL, theme = ttheme_header)
      title <- grid::textGrob("Report 2.1: Posterior Summary of Calibrated Sensitivity Parameters",
                              gp=grid::gpar(fontsize=14, fontface="bold"))
      gt <- gtable::gtable_add_rows(g, heights = grid::grobHeight(title) + padding * 2, pos = 0)
      gt <- gtable::gtable_add_grob(gt, title, t = 1, l = 1, r = ncol(gt))

      pdf("Report_2_1_Calibration_Summary.pdf", width = 11, height = 3 + 0.3 * nrow(Calib_formatted))
      grid::grid.draw(gt)
      dev.off()
    }
  }
}

# 2.2 Posterior summaries of key model parameters
D_max <- NULL
if (file.exists(CONFIG_PATHS$PROCESSED_DATA)) {
  processed_env <- new.env(); load(CONFIG_PATHS$PROCESSED_DATA, envir = processed_env)
  if (exists("D_max", envir = processed_env)) D_max <- processed_env$D_max
}
if (file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES)) {
  posterior_samples_report <- readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)
}
if (!is.null(posterior_samples_report)) {
  make_row <- function(param, desc, draws) {
    if (is.null(draws) || length(draws)==0 || all(is.na(draws)))
      return(data.frame(Parameter=param, Description=desc, Mean=NA_real_, SD=NA_real_, Q2.5=NA_real_, Q97.5=NA_real_))
    data.frame(Parameter=param, Description=desc,
               Mean=mean(draws, na.rm=TRUE), SD=sd(draws, na.rm=TRUE),
               Q2.5=quantile(draws,0.025,na.rm=TRUE), Q97.5=quantile(draws,0.975,na.rm=TRUE))
  }
  L <- list()
  if (!is.null(D_max) && D_max>0) {
    for (d in 1:D_max) {
      if (!is.null(posterior_samples_report$beta) && ncol(posterior_samples_report$beta)>=d)
        L[[paste0("beta[",d,"]")]]  <- make_row(paste0("beta[",d,"]"),  "Beta_d (Y Duration Effect)",  posterior_samples_report$beta[,d])
      if (!is.null(posterior_samples_report$gamma) && ncol(posterior_samples_report$gamma)>=d)
        L[[paste0("gamma[",d,"]")]] <- make_row(paste0("gamma[",d,"]"), "Gamma_d (M Duration Effect)", posterior_samples_report$gamma[,d])
    }
  }
  L[["psi3"]] <- make_row("psi3","Psi_3 (M-Y Assoc., D=0)", posterior_samples_report$psi3)
  L[["psi4"]] <- make_row("psi4","Psi_4 (M-Y Assoc. Diff, D>0)", posterior_samples_report$psi4)
  if (!is.null(posterior_samples_report$Sigma_alpha) && all(dim(posterior_samples_report$Sigma_alpha)[2:3]==c(2,2)))
    L[["Sigma_alpha[1,1]"]] <- make_row("Sigma_alpha[1,1]","Sigma^2_alpha1 (Cluster Var M)", posterior_samples_report$Sigma_alpha[,1,1])
  if (!is.null(posterior_samples_report$Sigma_phi) && all(dim(posterior_samples_report$Sigma_phi)[2:3]==c(2,2)))
    L[["Sigma_phi[1,1]"]]   <- make_row("Sigma_phi[1,1]","Sigma^2_phi1 (Indiv Var M)", posterior_samples_report$Sigma_phi[,1,1])
  if ("sigma_epsilon" %in% names(posterior_samples_report))
    L[["sigma_epsilon^2"]]  <- make_row("sigma_epsilon^2","Sigma^2_epsilon (Residual Var M)", (posterior_samples_report$sigma_epsilon)^2)

  Model_Param_Report <- dplyr::bind_rows(L)
  var_m_parts <- c("Sigma^2_alpha1 (Cluster Var M)", "Sigma^2_phi1 (Indiv Var M)", "Sigma^2_epsilon (Residual Var M)")
  Var_M_data <- Model_Param_Report %>% dplyr::filter(Description %in% var_m_parts & !is.na(Mean))
  if (nrow(Var_M_data)==length(var_m_parts)) {
    Var_M_total <- Var_M_data %>% dplyr::summarise(across(c(Mean,SD,Q2.5,Q97.5), ~sum(., na.rm=TRUE))) %>%
      dplyr::mutate(Parameter="Var(M)", Description="Total Variance of M")
    Model_Param_Report <- dplyr::bind_rows(Model_Param_Report, Var_M_total)
  }
  write.csv(Model_Param_Report, "Report_2_2_Model_Parameters_Summary.csv", row.names = FALSE)

  if (requireNamespace("gridExtra", quietly = TRUE) &&
      requireNamespace("grid", quietly = TRUE) &&
      requireNamespace("gtable", quietly = TRUE)) {

    ttheme_header <- gridExtra::ttheme_default(
      core = list(fg_params=list(fontsize=10), bg_params = list(fill = "grey95")),
      colhead = list(fg_params=list(fontsize=11, fontface="bold"), bg_params = list(fill = "grey85"))
    )
    padding <- grid::unit(0.5, "line")

    pdf_dat <- Model_Param_Report %>%
      dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.) | is.na(.), NA_character_, sprintf("%.4f", .)))) %>%
      dplyr::mutate(`95% CrI` = paste0("[", Q2.5, ", ", Q97.5, "]")) %>%
      dplyr::select(Description, Mean, SD, `95% CrI`)

    g <- gridExtra::tableGrob(pdf_dat, rows = NULL, theme = ttheme_header)
    title <- grid::textGrob("Report 2.2: Posterior Summaries of Key Model Parameters",
                            gp=grid::gpar(fontsize=14, fontface="bold"))
    gt <- gtable::gtable_add_rows(g, heights = grid::grobHeight(title) + padding * 2, pos = 0)
    gt <- gtable::gtable_add_grob(gt, title, t = 1, l = 1, r = ncol(gt))

    pdf("Report_2_2_Model_Parameters_Summary.pdf", width = 10, height = 3 + 0.4 * nrow(pdf_dat))
    grid::grid.draw(gt); dev.off()
  }
}

# 2.3 Stratum probabilities and 2.4 Figures
PCE_summary <- NULL; Probability_summary <- NULL; Lambda_summary <- NULL
if (file.exists(CONFIG_PATHS$PCE_RESULTS)) {
  res_env <- new.env(); load(CONFIG_PATHS$PCE_RESULTS, envir = res_env)
  PCE_summary <- res_env$PCE_summary
  Probability_summary <- res_env$Probability_summary
  Lambda_summary <- res_env$Lambda_summary

  if (!is.null(Probability_summary) && nrow(Probability_summary)>0) {
    Prob_report_data <- Probability_summary %>%
      dplyr::mutate(Stratum_Label = dplyr::case_when(
        Stratum == "PCE1_Dissociative"    ~ paste0("PCE1: Dissoc. (|dM| < ", sprintf("%.1f", Delta), ")"),
        Stratum == "PCE2_Associative_Neg" ~ paste0("PCE2: Assoc. Neg (dM < -", sprintf("%.1f", Delta), ")"),
        Stratum == "PCE3_Associative_Pos" ~ paste0("PCE3: Assoc. Pos (dM > ", sprintf("%.1f", Delta), ")"),
        TRUE ~ as.character(Stratum))) %>%
      dplyr::select(Duration, Type, Delta, Rho_Scenario, Stratum_Label, Prob_Mean, Prob_Q2.5, Prob_Q97.5) %>%
      dplyr::arrange(Duration, Delta, Rho_Scenario, Stratum_Label)

    write.csv(Prob_report_data, "Report_2_3_Stratum_Probabilities_Full.csv", row.names = FALSE)

    if (requireNamespace("gridExtra", quietly = TRUE) &&
        requireNamespace("grid", quietly = TRUE) &&
        requireNamespace("gtable", quietly = TRUE)) {
      ttheme_header <- gridExtra::ttheme_default(
        core = list(fg_params=list(fontsize=10), bg_params = list(fill = "grey95")),
        colhead = list(fg_params=list(fontsize=11, fontface="bold"), bg_params = list(fill = "grey85"))
      )
      padding <- grid::unit(0.5, "line")
      PRIMARY_DELTA <- PRIMARY_SETTINGS$DELTA
      PRIMARY_RHO   <- PRIMARY_SETTINGS$RHO_SCENARIO
      Prob_primary <- Prob_report_data %>% dplyr::filter(Delta == PRIMARY_DELTA, Rho_Scenario == PRIMARY_RHO)
      if (nrow(Prob_primary) > 0) {
        Prob_primary_fmt <- Prob_primary %>%
          dplyr::mutate(across(where(is.numeric), ~sprintf("%.4f", .x))) %>%
          dplyr::mutate(`95% CrI` = paste0("[", Prob_Q2.5, ", ", Prob_Q97.5, "]")) %>%
          dplyr::select(Duration, Type, Stratum=Stratum_Label, Mean=Prob_Mean, `95% CrI`)
        g <- gridExtra::tableGrob(Prob_primary_fmt, rows = NULL, theme = ttheme_header)
        title <- grid::textGrob(sprintf("Report 2.3: Principal Stratum Probabilities (Primary: Delta=%.1f, Rho=%s)",
                                        PRIMARY_DELTA, PRIMARY_RHO), gp=grid::gpar(fontsize=14, fontface="bold"))
        gt <- gtable::gtable_add_rows(g, heights = grid::grobHeight(title) + padding * 2, pos = 0)
        gt <- gtable::gtable_add_grob(gt, title, t = 1, l = 1, r = ncol(gt))

        pdf("Report_2_3_Stratum_Probabilities_Primary.pdf", width = 11.5, height = 3 + 0.3 * nrow(Prob_primary))
        grid::grid.draw(gt); dev.off()
      }
    }
  }
}

# Figures for PCE summaries
if (!is.null(PCE_summary) && nrow(PCE_summary)>0 && requireNamespace("ggplot2", quietly = TRUE)) {
  strata_levels <- c("PCE1 (Dissoc.)","PCE2 (Assoc. Neg)","PCE3 (Assoc. Pos)")
  strata_colors <- if (requireNamespace("RColorBrewer", quietly = TRUE)) RColorBrewer::brewer.pal(3,"Dark2") else c("#1B9E77","#D95F02","#7570B3")
  names(strata_colors) <- strata_levels
  strata_labels_legend <- setNames(
    list(expression("PCE"["D"]), expression("PCE"["A-"]), expression("PCE"["A+"])),
    strata_levels
  )
  strata_labels_facet_map <- setNames(c("\"PCE\"[\"D\"]", "\"PCE\"[\"A-\"]", "\"PCE\"[\"A+\"]"), strata_levels)
  strata_labeller_parsed <- ggplot2::as_labeller(strata_labels_facet_map, default = ggplot2::label_parsed)

  delta_values <- sort(unique(PCE_summary$Delta))
  k_values <- sort(unique(PCE_summary$Lambda_K))
  rho_scenario_levels_all <- unique(PCE_summary$Rho_Scenario)

  delta_colors <- if (requireNamespace("RColorBrewer", quietly = TRUE))
    stats::setNames(colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")[5:9])(length(delta_values)), as.character(delta_values)) else c()
  k_cols <- if (requireNamespace("RColorBrewer", quietly = TRUE))
    stats::setNames(colorRampPalette(RColorBrewer::brewer.pal(9,"Purples")[5:9])(length(k_values)), as.character(k_values)) else c()
  if (length(k_values)>0) {
    is_base <- sapply(k_values, function(x) isTRUE(all.equal(x, 1.0)))
    if (any(is_base)) k_cols[as.character(k_values[is_base])] <- strata_colors[1]
  }

  mean_calibrated_rho_val <- NA
  if ("Calibrated_Stochastic" %in% rho_scenario_levels_all && exists("Calibration_report_data")) {
    if (is.numeric(Calibration_report_data$Mean)) {
      tmp <- Calibration_report_data %>% dplyr::filter(Param == "Rho* (Calibrated)")
      if (nrow(tmp)>0) mean_calibrated_rho_val <- tmp$Mean[1]
    }
  }
  rho_labels_map <- list()
  for (scenario in rho_scenario_levels_all) {
    rho_labels_map[[scenario]] <-
      if (scenario == "Calibrated_Stochastic") {
        if (!is.na(mean_calibrated_rho_val)) sprintf("Calibrated (Mean Rho* approx %.3f)", mean_calibrated_rho_val) else "Calibrated (Stochastic)"
      } else if (startsWith(scenario,"Fixed_")) sprintf("Fixed (Rho = %s)", sub("Fixed_","",scenario)) else scenario
  }
  rho_order <- c(grep("Calibrated", rho_scenario_levels_all, value = TRUE),
                 sort(grep("Fixed_", rho_scenario_levels_all, value = TRUE)))
  rho_colors <- if (length(rho_order)>0 && requireNamespace("RColorBrewer", quietly = TRUE))
    stats::setNames(colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(length(rho_order)), rho_order) else c()

  PCE_plot_data <- PCE_summary %>%
    dplyr::mutate(Stratum_Short = factor(dplyr::case_when(
      Stratum == "PCE1_Dissociative"    ~ "PCE1 (Dissoc.)",
      Stratum == "PCE2_Associative_Neg" ~ "PCE2 (Assoc. Neg)",
      Stratum == "PCE3_Associative_Pos" ~ "PCE3 (Assoc. Pos)"
    ), levels = strata_levels)) %>%
    dplyr::filter(!is.na(Type) & !is.na(Duration)) %>%
    dplyr::mutate(Duration = as.numeric(as.character(Duration)))

  if (nrow(PCE_plot_data) > 0) {
    PCE_plot_data <- PCE_plot_data %>%
      dplyr::mutate(Type = factor(Type, levels = unique(Type[order(Duration)])),
                    Province = factor(Province))
    if (length(rho_order)>0) PCE_plot_data$Rho_Scenario <- factor(PCE_plot_data$Rho_Scenario, levels = rho_order)
    if (length(k_values)>0) PCE_plot_data$Lambda_K_Factor <- factor(PCE_plot_data$Lambda_K, levels = k_values)
  }
  n_prov <- max(1, length(unique(PCE_plot_data$Province)))

  PRIMARY_DELTA_VAL <- PRIMARY_SETTINGS$DELTA
  PRIMARY_RHO_VAL   <- PRIMARY_SETTINGS$RHO_SCENARIO
  PRIMARY_LAMBDA_K_VAL <- PRIMARY_SETTINGS$LAMBDA_K
  device_to_use <- if (capabilities("cairo") && exists("cairo_pdf", where = "package:grDevices")) grDevices::cairo_pdf else "pdf"

  # Fig 1: Primary analysis across strata
  PCE_primary <- PCE_plot_data %>%
    dplyr::filter(Lambda_K == PRIMARY_LAMBDA_K_VAL, Delta == PRIMARY_DELTA_VAL, Rho_Scenario == PRIMARY_RHO_VAL)
  if (nrow(PCE_primary)>0) {
    p1 <- ggplot2::ggplot(PCE_primary, ggplot2::aes(x = Time, y = Mean, color = Stratum_Short)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 2) +
      ggplot2::geom_line(ggplot2::aes(group = Stratum_Short), position = ggplot2::position_dodge(width = 0.6), linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.4, position = ggplot2::position_dodge(width = 0.6)) +
      ggplot2::facet_grid(Province ~ Type, scales = "free_x") +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI", color = "Principal Stratum") +
      ggplot2::scale_color_manual(values = strata_colors, labels = strata_labels_legend) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position="bottom", strip.text=ggplot2::element_text(face="bold"))
    n_dur <- length(unique(PCE_plot_data$Duration))
    ggplot2::ggsave("Plot1_PCE_Strata_Comparison_Primary.pdf", plot = p1,
                    width = 4 + 2.5 * max(1, n_dur), height = 4 * n_prov + 1.5, device = device_to_use)
  }

  # Fig 2: Rho scenario sensitivity (d = 1)
  PCE_rho <- PCE_plot_data %>%
    dplyr::filter(Type == "stPCE (d=1)", Delta == PRIMARY_DELTA_VAL, Lambda_K == PRIMARY_LAMBDA_K_VAL)
  if (nrow(PCE_rho)>0 && length(unique(PCE_rho$Rho_Scenario))>1) {
    lv <- levels(droplevels(PCE_rho$Rho_Scenario)); labs_now <- sapply(lv, function(s) rho_labels_map[[s]])
    p2 <- ggplot2::ggplot(PCE_rho, ggplot2::aes(x = Time, y = Mean, color = Rho_Scenario)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI", color = "Rho Scenario") +
      ggplot2::scale_color_manual(values = rho_colors, labels = labs_now, drop = TRUE) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom", strip.text=ggplot2::element_text(face="bold"))
    n_strata <- length(unique(PCE_plot_data$Stratum_Short))
    ggplot2::ggsave("Plot2_PCE_Rho_Sensitivity.pdf", plot = p2,
                    width = 2 + 3.5 * max(1, n_strata), height = 4 * n_prov + 1.8, device = device_to_use)
  }

  # Fig 3: Lambda scaling (k) sensitivity (d = 1)
  PCE_k <- PCE_plot_data %>% dplyr::filter(Type=="stPCE (d=1)", Delta==PRIMARY_DELTA_VAL, Rho_Scenario==PRIMARY_RHO_VAL)
  if (nrow(PCE_k)>0 && length(k_values)>1) {
    k_labels <- paste0("k=", k_values); is_base <- sapply(k_values, function(x) isTRUE(all.equal(x,1.0)))
    if (any(is_base)) k_labels[is_base] <- "k=1.0 (Baseline)"; names(k_labels) <- as.character(k_values)
    p4 <- ggplot2::ggplot(PCE_k, ggplot2::aes(x = Time, y = Mean, color = Lambda_K_Factor)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI", color = "Lambda Scaling (k)") +
      ggplot2::scale_color_manual(values = k_cols, labels = k_labels, drop = TRUE) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom", strip.text=ggplot2::element_text(face="bold"))
    n_strata <- length(unique(PCE_plot_data$Stratum_Short))
    ggplot2::ggsave("Plot4_PCE_Lambda_K_Sensitivity.pdf", plot = p4,
                    width = 2 + 3.5 * max(1, n_strata), height = 4 * n_prov + 1.8, device = device_to_use)
  }

  # Fig 4: Delta sensitivity (d = 1)
  PCE_delta <- PCE_plot_data %>% dplyr::filter(Type=="stPCE (d=1)", Rho_Scenario==PRIMARY_RHO_VAL, Lambda_K==PRIMARY_LAMBDA_K_VAL) %>%
    dplyr::mutate(Delta_Factor = factor(Delta))
  if (nrow(PCE_delta)>0 && length(unique(PCE_delta$Delta))>1) {
    p3 <- ggplot2::ggplot(PCE_delta, ggplot2::aes(x = Time, y = Mean, color = Delta_Factor)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Q2.5, ymax = Q97.5), width = 0.5, position = ggplot2::position_dodge(width = 0.8)) +
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller = ggplot2::labeller(Stratum_Short = strata_labeller_parsed)) +
      ggplot2::labs(y = "PCE (Probability Difference) with 95% CrI", color = "Delta Value") +
      ggplot2::scale_color_manual(values = delta_colors, drop = TRUE) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom", strip.text=ggplot2::element_text(face="bold"))
    n_strata <- length(unique(PCE_plot_data$Stratum_Short))
    ggplot2::ggsave("Plot3_PCE_Delta_Sensitivity.pdf", plot = p3,
                    width = 2 + 3.5 * max(1, n_strata), height = 4 * n_prov + 1.5, device = device_to_use)
  }

  # Fig 5: Lambda distribution (supplementary)
  if (!is.null(Lambda_summary) && nrow(Lambda_summary)>0) {
    Lambda_plot <- Lambda_summary %>%
      dplyr::filter(!is.na(Duration)) %>%
      dplyr::mutate(Duration = as.numeric(as.character(Duration)))
    if (nrow(Lambda_plot)>0) {
      dur_lv <- sort(unique(Lambda_plot$Duration))
      Lambda_plot <- Lambda_plot %>%
        dplyr::mutate(Lambda_Type_Label = ifelse(Lambda_Type=="Lambda_z", "lambda(z)", "lambda(-z)"),
                      Duration_Label = factor(paste0("Duration d=", Duration), levels = paste0("Duration d=", dur_lv)))
      kv <- sort(unique(Lambda_plot$Lambda_K))
      Lambda_plot$Lambda_K_Factor <- factor(Lambda_plot$Lambda_K, levels = kv)
      k_labels <- paste0("k=", kv); isb <- sapply(kv, function(x) isTRUE(all.equal(x,1.0)))
      if (any(isb)) k_labels[isb] <- "k=1.0 (Baseline)"; names(k_labels) <- as.character(kv)

      p5 <- ggplot2::ggplot(Lambda_plot, ggplot2::aes(x = Value, fill = Lambda_K_Factor, color = Lambda_K_Factor)) +
        ggplot2::geom_density(alpha = 0.4) +
        ggplot2::facet_grid(Lambda_Type_Label ~ Duration_Label, scales = "free_y", labeller = ggplot2::label_parsed) +
        ggplot2::labs(x="Lambda Value", y="Density", fill="Scaling (k)", color="Scaling (k)") +
        ggplot2::scale_fill_manual(values = k_cols, labels = k_labels, drop = TRUE) +
        ggplot2::scale_color_manual(values = k_cols, labels = k_labels, drop = TRUE) +
        ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom", strip.text=ggplot2::element_text(face="bold"))
      ggplot2::ggsave("Plot5_Lambda_K_Distributions.pdf", plot = p5,
                      width = 3 + 3 * max(1, length(dur_lv)), height = 7, device = device_to_use)
    }
  }
}
