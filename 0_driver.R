# 0_driver.R — streamlined driver
# Run steps, keep outputs identical; lighter logs and comments.

# Optional working directory
# setwd('~/PROJECT_ps/swd_v22/')

# Packages
pkgs <- c(
  "dplyr","tidyr","readxl","tibble","forcats","rstan","statmod","mvtnorm","triangle",
  "ggplot2","foreach","parallel","Matrix","broom","RColorBrewer","gridExtra","grid",
  "gtable","doSNOW","purrr"
)
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))

# Paths and config
CONFIG_PATHS <- list(
  DATA_PATH          = "../data/pmed.1002645.s007.xlsx",
  STAN_MODEL         = "model_duration.stan",
  STAN_AUX_RHO       = "aux_model_rho.stan",
  STAN_AUX_LAMBDA    = "aux_model_lambda.stan",
  PROCESSED_DATA     = "Processed_Data.RData",
  POSTERIOR_SAMPLES  = "Posterior_Samples.rds",
  CALIBRATION_RESULTS= "Calibration_Results_Bayesian.rds",
  PCE_RESULTS        = "PCE_Results.RData",
  ANALYSIS_CONFIG    = "Analysis_Config.RData"
)

RUN_STEPS <- list(
  PREPROCESS    = TRUE,
  MODEL_FIT     = TRUE,
  CALIBRATION   = TRUE,
  PCE_ESTIMATION= TRUE
)

detect_cores_hpc <- function() {
  n <- parallel::detectCores()
  if (is.null(n) || is.na(n) || n < 1) 1L else n
}
n_cores <- detect_cores_hpc()

CONFIG_STAN <- list(ITER=4000, WARMUP=2000, CHAINS=4, CORES=n_cores, ADAPT_DELTA=0.95, MAX_TREEDEPTH=15)
CONFIG_STAN_AUX <- list(ITER=4000, WARMUP=2000, CHAINS=4, CORES=n_cores, ADAPT_DELTA=0.90, MAX_TREEDEPTH=12)

CONFIG_PCE <- list(
  DURATIONS = c(1,2,3),
  DELTA_LIST = c(0.5,1.0,1.5),
  RHO_SCENARIOS = c("Calibrated_Stochastic","Fixed_0.8","Fixed_0.9"),
  LAMBDA_K_SCALING = c(0.5,1.0,1.5,2.0),
  S_SUBSET = 200,
  N_INDIV_SAMPLE = 100,
  N_MC_INTEGRATION = 200
)

PRIMARY_SETTINGS <- list(DELTA=0.5, RHO_SCENARIO="Calibrated_Stochastic", LAMBDA_K=1.0)

# Save/load config
if (!any(unlist(RUN_STEPS))) {
  if (file.exists(CONFIG_PATHS$ANALYSIS_CONFIG)) {
    load(CONFIG_PATHS$ANALYSIS_CONFIG)
  }
} else {
  save(CONFIG_PATHS, CONFIG_STAN, CONFIG_STAN_AUX, CONFIG_PCE, PRIMARY_SETTINGS, file = CONFIG_PATHS$ANALYSIS_CONFIG)
}

# Make data_path visible to preprocess
data_path <- CONFIG_PATHS$DATA_PATH

# Pipeline
if (RUN_STEPS$PREPROCESS)    source("1_preprocess.R")         
if (RUN_STEPS$MODEL_FIT)     source("2_fit_model.R")
if (RUN_STEPS$CALIBRATION)   source("3_calibration.R")
if (RUN_STEPS$PCE_ESTIMATION)source("4_PCE_estimation.R")

# =========================
# Reports and simple plots
# =========================

# Table theme
ttheme_header <- gridExtra::ttheme_default(
  core   = list(fg_params=list(fontsize=10), bg_params=list(fill="grey95")),
  colhead= list(fg_params=list(fontsize=11,fontface="bold"), bg_params=list(fill="grey85"))
)
padding <- unit(0.5,"line")

posterior_samples_report <- NULL
calibration_results_report <- NULL
Calibration_report_data <- NULL

# ----- 2.1 Calibrated sensitivity params (Bayesian summary) -----
if (file.exists(CONFIG_PATHS$CALIBRATION_RESULTS)) {
  calibration_results_report <- readRDS(CONFIG_PATHS$CALIBRATION_RESULTS)
}
if (!is.null(calibration_results_report) && !is.null(calibration_results_report$Posteriors)) {
  Posteriors <- calibration_results_report$Posteriors
  summarize_posterior <- function(draws, name, dur=NA) {
    if (is.null(draws) || length(draws)==0 || all(is.na(draws))) return(NULL)
    data.frame(
      Param = name, Duration = dur,
      Mean = mean(draws,na.rm=TRUE),
      Median = median(draws,na.rm=TRUE),
      SD = sd(draws,na.rm=TRUE),
      Q2.5 = quantile(draws,0.025,na.rm=TRUE),
      Q97.5= quantile(draws,0.975,na.rm=TRUE)
    )
  }
  s_rho <- summarize_posterior(Posteriors$rho_star, "Rho* (Calibrated)", "Independent")
  D_lz  <- if (is.matrix(Posteriors$lambda_z_lower)) ncol(Posteriors$lambda_z_lower) else 0
  D_lmz <- if (is.matrix(Posteriors$lambda_minus_z_lower)) ncol(Posteriors$lambda_minus_z_lower) else 0
  D_max_report <- max(D_lz, D_lmz)

  s_lz <- s_lmz <- data.frame()
  if (D_max_report > 0) for (d in 1:D_max_report) {
    if (d <= D_lz)  s_lz  <- rbind(s_lz,  summarize_posterior(Posteriors$lambda_z_lower[,d],  "Lambda_z Calibrated (Theta_2(d))", d))
    if (d <= D_lmz) s_lmz <- rbind(s_lmz, summarize_posterior(Posteriors$lambda_minus_z_lower[,d], "Lambda_minus_z Calibrated (Zeta_1(d))", d))
  }
  Calibration_report_data <- dplyr::bind_rows(s_rho, s_lz, s_lmz) |>
    dplyr::mutate(Duration_num = suppressWarnings(as.numeric(as.character(Duration))),
                  Duration_Sort = ifelse(is.na(Duration_num),0,Duration_num)) |>
    dplyr::arrange(Duration_Sort, Param) |>
    dplyr::select(-Duration_Sort,-Duration_num)

  if (nrow(Calibration_report_data)>0) {
    write.csv(Calibration_report_data, "Report_2_1_Calibration_Bounds_Bayesian_Summary.csv", row.names=FALSE)

    Calib_fmt <- Calibration_report_data |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~sprintf("%.4f", .x))) |>
      dplyr::mutate(`95% CrI` = paste0("[", Q2.5,", ",Q97.5,"]")) |>
      dplyr::select(Duration, Parameter=Param, Mean, SD, `95% CrI`)

    tab <- gridExtra::tableGrob(Calib_fmt, rows=NULL, theme=ttheme_header)
    title <- grid::textGrob("Report 2.1: Posterior Summary of Calibrated Sensitivity Parameters", gp=grid::gpar(fontsize=14,fontface="bold"))
    tab2 <- gtable::gtable_add_rows(tab, heights = grid::grobHeight(title) + padding*2, pos=0)
    tab2 <- gtable::gtable_add_grob(tab2, title, t=1, l=1, r=ncol(tab2))
    pdf("Report_2_1_Calibration_Bounds_Bayesian_Summary.pdf", width=11, height=3 + 0.3*nrow(Calib_fmt))
    grid::grid.draw(tab2); dev.off()
  }
}

# ----- 2.2 Key model parameters (posterior) -----
D_max <- NULL
if (file.exists(CONFIG_PATHS$PROCESSED_DATA)) {
  env <- new.env(); load(CONFIG_PATHS$PROCESSED_DATA, envir=env)
  if (exists("D_max", envir=env)) D_max <- as.integer(env$D_max)
}
if (file.exists(CONFIG_PATHS$POSTERIOR_SAMPLES)) posterior_samples_report <- readRDS(CONFIG_PATHS$POSTERIOR_SAMPLES)

if (!is.null(posterior_samples_report)) {
  params <- list()
  if (!is.null(D_max) && D_max > 0) {
    for (d in 1:D_max) {
      params[[paste0("beta[",d,"]")]]  <- "Beta_d (Y Duration Effect)"
      params[[paste0("gamma[",d,"]")]] <- "Gamma_d (M Duration Effect)"
    }
  }
  params[["psi3"]] <- "Psi_3 (M-Y Assoc., D=0)"
  params[["psi4"]] <- "Psi_4 (M-Y Assoc. Diff, D>0)"
  params[["Sigma_alpha[1,1]"]] <- "Sigma^2_alpha1 (Cluster Var M)"
  params[["Sigma_phi[1,1]"]]   <- "Sigma^2_phi1 (Indiv Var M)"

  if ("sigma_epsilon" %in% names(posterior_samples_report)) {
    se <- posterior_samples_report$sigma_epsilon
    posterior_samples_report[["sigma_epsilon^2"]] <- if (is.matrix(se)) se[,1]^2 else se^2
    params[["sigma_epsilon^2"]] <- "Sigma^2_epsilon (Residual Var M)"
  }

  parse_idx <- function(nm) {
    s <- regmatches(nm, regexpr("\\[(.*?)\\]", nm))
    if (length(s)==0) return(NULL)
    parts <- strsplit(gsub("^\\[|\\]$","",s), ",")[[1]]
    as.numeric(trimws(parts))
  }
  extract_draws <- function(x, idx) {
    if (is.null(x) || !is.numeric(x)) return(NULL)
    d <- dim(x)
    if (is.null(idx)) return(as.vector(if (!is.null(d) && length(d)==2 && d[2]==1) x[,1] else x))
    if (length(idx)==1) {
      if (!is.null(d) && length(d)==2 && idx <= d[2]) return(as.vector(x[,idx]))
      if (is.null(d) && idx==1) return(as.vector(x))
    }
    if (length(idx)==2 && !is.null(d) && length(d)==3 && idx[1] <= d[2] && idx[2] <= d[3]) return(as.vector(x[,idx[1],idx[2]]))
    NULL
  }
  row_of <- function(nm, desc, draws) {
    if (is.null(draws) || length(draws)==0 || all(is.na(draws))) {
      data.frame(Parameter=nm, Description=desc, Mean=NA_real_, SD=NA_real_, Q2.5=NA_real_, Q97.5=NA_real_)
    } else {
      data.frame(Parameter=nm, Description=desc,
                 Mean=mean(draws,na.rm=TRUE), SD=sd(draws,na.rm=TRUE),
                 Q2.5=quantile(draws,0.025,na.rm=TRUE), Q97.5=quantile(draws,0.975,na.rm=TRUE))
    }
  }

  out <- list()
  for (nm in names(params)) {
    base <- sub("\\[.*$","",nm)
    idx  <- parse_idx(nm)
    dr   <- if (base %in% names(posterior_samples_report)) extract_draws(posterior_samples_report[[base]], idx) else NULL
    out[[nm]] <- row_of(nm, params[[nm]], dr)
  }
  Model_Param_Report <- dplyr::bind_rows(out)

  # Var(M) = cluster + indiv + residual (sum of means)
  parts <- c("Sigma^2_alpha1 (Cluster Var M)","Sigma^2_phi1 (Indiv Var M)","Sigma^2_epsilon (Residual Var M)")
  comp  <- Model_Param_Report |>
    dplyr::filter(Description %in% parts, !is.na(Mean)) |>
    dplyr::summarise(dplyr::across(c(Mean,SD,Q2.5,Q97.5), ~sum(., na.rm=TRUE))) |>
    dplyr::mutate(Parameter="Var(M)", Description="Total Variance of M")
  if (nrow(comp)==1) Model_Param_Report <- dplyr::bind_rows(Model_Param_Report, comp)

  write.csv(Model_Param_Report, "Report_2_2_Model_Parameters_Summary.csv", row.names=FALSE)

  pdf_tbl <- Model_Param_Report |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ifelse(is.nan(.)|is.na(.), NA_character_, sprintf("%.4f",.)))) |>
    dplyr::mutate(`95% CrI` = paste0("[", Q2.5,", ",Q97.5,"]")) |>
    dplyr::select(Description, Mean, SD, `95% CrI`)
  tab <- gridExtra::tableGrob(pdf_tbl, rows=NULL, theme=ttheme_header)
  title <- grid::textGrob("Report 2.2: Posterior Summaries of Key Model Parameters", gp=grid::gpar(fontsize=14,fontface="bold"))
  tab2 <- gtable::gtable_add_rows(tab, heights=grid::grobHeight(title)+padding*2, pos=0)
  tab2 <- gtable::gtable_add_grob(tab2, title, t=1, l=1, r=ncol(tab2))
  pdf("Report_2_2_Model_Parameters_Summary.pdf", width=10, height=3 + 0.4*nrow(pdf_tbl))
  grid::grid.draw(tab2); dev.off()
}

# ----- 2.3 Principal stratum probabilities -----
PCE_summary <- Probability_summary <- Lambda_summary <- NULL
if (file.exists(CONFIG_PATHS$PCE_RESULTS)) {
  env <- new.env(); load(CONFIG_PATHS$PCE_RESULTS, envir=env)
  PCE_summary       <- env$PCE_summary
  Probability_summary <- env$Probability_summary
  Lambda_summary    <- env$Lambda_summary
}
if (!is.null(Probability_summary) && nrow(Probability_summary)>0) {
  Prob_report_data <- Probability_summary |>
    dplyr::mutate(
      Stratum_Label = dplyr::case_when(
        Stratum=="PCE1_Dissociative"     ~ paste0("PCE1: Dissoc. (|dM| < ", sprintf("%.1f", Delta), ")"),
        Stratum=="PCE2_Associative_Neg"  ~ paste0("PCE2: Assoc. Neg (dM < -", sprintf("%.1f", Delta), ")"),
        Stratum=="PCE3_Associative_Pos"  ~ paste0("PCE3: Assoc. Pos (dM > ", sprintf("%.1f", Delta), ")"),
        TRUE ~ as.character(Stratum)
      )
    ) |>
    dplyr::select(Duration, Type, Delta, Rho_Scenario, Stratum_Label, Prob_Mean, Prob_Q2.5, Prob_Q97.5) |>
    dplyr::arrange(Duration, Delta, Rho_Scenario, Stratum_Label)
  write.csv(Prob_report_data, "Report_2_3_Stratum_Probabilities_Full.csv", row.names=FALSE)

  PRIMARY_DELTA <- PRIMARY_SETTINGS$DELTA
  PRIMARY_RHO   <- PRIMARY_SETTINGS$RHO_SCENARIO
  Prob_primary <- dplyr::filter(Prob_report_data, Delta==PRIMARY_DELTA & Rho_Scenario==PRIMARY_RHO)
  if (nrow(Prob_primary)>0) {
    Prob_fmt <- Prob_primary |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~sprintf("%.4f", .x))) |>
      dplyr::mutate(`95% CrI` = paste0("[", Prob_Q2.5,", ",Prob_Q97.5,"]")) |>
      dplyr::select(Duration, Type, Stratum=Stratum_Label, Mean=Prob_Mean, `95% CrI`)
    tab <- gridExtra::tableGrob(Prob_fmt, rows=NULL, theme=ttheme_header)
    ttl <- grid::textGrob(sprintf("Report 2.3: Principal Stratum Probabilities (Delta=%.1f, Rho=%s)", PRIMARY_DELTA, PRIMARY_RHO),
                          gp=grid::gpar(fontsize=14,fontface="bold"))
    tab2 <- gtable::gtable_add_rows(tab, heights=grid::grobHeight(ttl)+padding*2, pos=0)
    tab2 <- gtable::gtable_add_grob(tab2, ttl, t=1, l=1, r=ncol(tab2))
    pdf("Report_2_3_Stratum_Probabilities_Primary.pdf", width=11.5, height=3 + 0.3*nrow(Prob_primary))
    grid::grid.draw(tab2); dev.off()
  }
}

# ----- 2.4 Plots (PCE + Lambda) -----
if (!is.null(PCE_summary) && nrow(PCE_summary)>0 &&
    all(c("Delta","Lambda_K","Rho_Scenario") %in% names(PCE_summary))) {

  okabe_ito_no_yellow <- c("#0072B2","#D55E00","#009E73","#CC79A7","#000000","#7F7F7F","#005AB5","#332288")
  get_no_yellow <- function(n) {
    if (n <= length(okabe_ito_no_yellow)) okabe_ito_no_yellow[1:n]
    else grDevices::hcl(h = seq(10, 350, length.out = n), c = 85, l = 45)
  }
  strata_levels <- c("PCE1 (Dissoc.)","PCE2 (Assoc. Neg)","PCE3 (Assoc. Pos)")
  strata_colors <- if (requireNamespace("RColorBrewer", quietly=TRUE)) RColorBrewer::brewer.pal(3,"Dark2") else c("#1B9E77","#D95F02","#7570B3")
  names(strata_colors) <- strata_levels
  strata_labels_legend <- setNames(list(expression("PCE"["D"]), expression("PCE"["A-"]), expression("PCE"["A+"])), strata_levels)
  strata_labels_facet_map <- setNames(c("\"PCE\"[\"D\"]","\"PCE\"[\"A-\"]","\"PCE\"[\"A+\"]"), strata_levels)
  strata_labeller_parsed <- ggplot2::as_labeller(strata_labels_facet_map, default = ggplot2::label_parsed)

  delta_values <- sort(unique(PCE_summary$Delta))
  delta_colors <- if (length(delta_values)>0) grDevices::hcl(h = seq(200,260,length.out=length(delta_values)), c=85,l=45) else c()
  names(delta_colors) <- as.character(delta_values)

  k_values <- sort(unique(PCE_summary$Lambda_K))
  k_colors <- if (length(k_values)>0) grDevices::hcl(h = seq(280,340,length.out=length(k_values)), c=85,l=45) else c()
  names(k_colors) <- as.character(k_values)
  if ("1" %in% names(k_colors)) k_colors["1"] <- strata_colors[1]  # highlight k=1

  PCE_plot_data <- PCE_summary |>
    dplyr::mutate(
      Stratum_Label = dplyr::case_when(
        Stratum=="PCE1_Dissociative"    ~ paste0("PCE1: Dissociative (|dM| < ", sprintf("%.1f", Delta), ")"),
        Stratum=="PCE2_Associative_Neg" ~ paste0("PCE2: Associative (dM < -", sprintf("%.1f", Delta), ")"),
        Stratum=="PCE3_Associative_Pos" ~ paste0("PCE3: Associative (dM > ", sprintf("%.1f", Delta), ")")
      ),
      Stratum_Short = factor(dplyr::case_when(
        Stratum=="PCE1_Dissociative"    ~ "PCE1 (Dissoc.)",
        Stratum=="PCE2_Associative_Neg" ~ "PCE2 (Assoc. Neg)",
        Stratum=="PCE3_Associative_Pos" ~ "PCE3 (Assoc. Pos)"
      ), levels=strata_levels)
    ) |>
    dplyr::filter(!is.na(Type) & !is.na(Duration)) |>
    dplyr::mutate(Duration = as.numeric(as.character(Duration)))
  if (nrow(PCE_plot_data)>0) {
    PCE_plot_data <- PCE_plot_data |>
      dplyr::mutate(Type = factor(Type, levels = unique(Type[order(Duration)])),
                    Province = factor(Province))
  }

  rho_levels <- unique(PCE_plot_data$Rho_Scenario)
  rho_colors <- setNames(get_no_yellow(length(rho_levels)), rho_levels)

  theme_science <- ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position="bottom",
      panel.spacing.x=ggplot2::unit(1.5,"lines"),
      panel.spacing.y=ggplot2::unit(1,"lines"),
      strip.text.x=ggplot2::element_text(face="bold",size=11),
      strip.text.y=ggplot2::element_text(face="bold",size=11),
      legend.text=ggplot2::element_text(size=10),
      legend.title=ggplot2::element_text(size=11,face="bold"),
      axis.text=ggplot2::element_text(size=9),
      axis.title=ggplot2::element_text(size=10,face="bold"),
      plot.title=ggplot2::element_blank(),
      plot.subtitle=ggplot2::element_blank()
    )

  dev_fun <- if (capabilities("cairo") && exists("cairo_pdf", where="package:grDevices")) grDevices::cairo_pdf else "pdf"

  # Plot 1: primary
  PCE_primary <- dplyr::filter(PCE_plot_data,
                               Lambda_K==PRIMARY_SETTINGS$LAMBDA_K,
                               Delta==PRIMARY_SETTINGS$DELTA,
                               Rho_Scenario==PRIMARY_SETTINGS$RHO_SCENARIO)
  if (nrow(PCE_primary)>0) {
    p1 <- ggplot2::ggplot(PCE_primary, ggplot2::aes(x=Time,y=Mean,color=Stratum_Short))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=0.6), size=2)+
      ggplot2::geom_line(ggplot2::aes(group=Stratum_Short), position=ggplot2::position_dodge(width=0.6), linetype="dashed", alpha=0.6)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Q2.5,ymax=Q97.5), width=0.4, position=ggplot2::position_dodge(width=0.6))+
      ggplot2::facet_grid(Province ~ Type, scales="free_x")+
      ggplot2::labs(y="PCE (Probability Difference) with 95% CrI", color="Principal Stratum")+
      ggplot2::scale_color_manual(values=strata_colors, labels=strata_labels_legend)+
      ggplot2::scale_x_continuous(breaks=function(x){if (length(x)==0) numeric(0) else seq(floor(min(x,na.rm=TRUE)), ceiling(max(x,na.rm=TRUE)), by=1)})+
      theme_science
    n_dur <- length(unique(PCE_plot_data$Duration))
    ggplot2::ggsave("Plot1_PCE_Strata_Comparison_Primary.pdf", plot=p1, width=4+2.5*max(1,n_dur), height=5, device=dev_fun)
  }

  # Plot 2: rho sensitivity (stPCE d=1)
  PCE_rho <- dplyr::filter(PCE_plot_data, Type=="stPCE (d=1)", Delta==PRIMARY_SETTINGS$DELTA, Lambda_K==PRIMARY_SETTINGS$LAMBDA_K)
  if (nrow(PCE_rho)>0 && length(unique(PCE_rho$Rho_Scenario))>1) {
    p2 <- ggplot2::ggplot(PCE_rho, ggplot2::aes(x=Time,y=Mean,color=Rho_Scenario))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=0.8), size=2)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Q2.5,ymax=Q97.5), width=0.5, position=ggplot2::position_dodge(width=0.8))+
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller=ggplot2::labeller(Stratum_Short=strata_labeller_parsed))+
      ggplot2::labs(y="PCE (Probability Difference) with 95% CrI", color="Rho Scenario")+
      ggplot2::scale_color_manual(values=rho_colors, drop=TRUE)+
      theme_science
    ggplot2::ggsave("Plot2_PCE_Rho_Sensitivity.pdf", plot=p2, width=2+3.5*3, height=5, device=dev_fun)
  }

  # Plot 3: delta sensitivity (stPCE d=1)
  PCE_delta <- dplyr::filter(PCE_plot_data, Type=="stPCE (d=1)", Rho_Scenario==PRIMARY_SETTINGS$RHO_SCENARIO, Lambda_K==PRIMARY_SETTINGS$LAMBDA_K) |>
    dplyr::mutate(Delta_Factor = factor(Delta))
  if (nrow(PCE_delta)>0 && length(unique(PCE_delta$Delta))>1) {
    p3 <- ggplot2::ggplot(PCE_delta, ggplot2::aes(x=Time,y=Mean,color=Delta_Factor))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=0.8), size=2)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Q2.5,ymax=Q97.5), width=0.5, position=ggplot2::position_dodge(width=0.8))+
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller=ggplot2::labeller(Stratum_Short=strata_labeller_parsed))+
      ggplot2::labs(y="PCE (Probability Difference) with 95% CrI", color="Delta Value")+
      ggplot2::scale_color_manual(values=delta_colors, drop=TRUE)+
      theme_science
    ggplot2::ggsave("Plot3_PCE_Delta_Sensitivity.pdf", plot=p3, width=2+3.5*3, height=5, device=dev_fun)
  }

  # Plot 4: lambda k sensitivity (stPCE d=1)
  PCE_k <- dplyr::filter(PCE_plot_data, Type=="stPCE (d=1)", Delta==PRIMARY_SETTINGS$DELTA, Rho_Scenario==PRIMARY_SETTINGS$RHO_SCENARIO) |>
    dplyr::mutate(Lambda_K_Factor=factor(Lambda_K))
  if (nrow(PCE_k)>0 && length(unique(PCE_k$Lambda_K))>1) {
    k_labels <- setNames(paste0("k=",k_values), as.character(k_values))
    if ("1" %in% names(k_labels)) k_labels["1"] <- "k=1.0 (Baseline)"
    p4 <- ggplot2::ggplot(PCE_k, ggplot2::aes(x=Time,y=Mean,color=Lambda_K_Factor))+
      ggplot2::geom_point(position=ggplot2::position_dodge(width=0.8), size=2)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Q2.5,ymax=Q97.5), width=0.5, position=ggplot2::position_dodge(width=0.8))+
      ggplot2::facet_grid(Province ~ Stratum_Short, labeller=ggplot2::labeller(Stratum_Short=strata_labeller_parsed))+
      ggplot2::labs(y="PCE (Probability Difference) with 95% CrI", color="Lambda Scaling (k)")+
      ggplot2::scale_color_manual(values=k_colors, labels=k_labels, drop=TRUE)+
      theme_science
    ggplot2::ggsave("Plot4_PCE_Lambda_K_Sensitivity.pdf", plot=p4, width=2+3.5*3, height=5, device=dev_fun)
  }

  # Plot 5: lambda distributions
  if (!is.null(Lambda_summary) && nrow(Lambda_summary)>0) {
    Ldat <- Lambda_summary |>
      dplyr::filter(!is.na(Duration)) |>
      dplyr::mutate(Duration = as.numeric(as.character(Duration)),
                    Lambda_Type_Label = ifelse(Lambda_Type=="Lambda_z","lambda(z)","lambda(-z)"),
                    Duration_Label = factor(paste0("Duration d=", Duration), levels=paste0("Duration d=", sort(unique(Duration)))))
    if (nrow(Ldat)>0 && "Lambda_K" %in% names(Ldat)) {
      Ldat$Lambda_K_Factor <- factor(Ldat$Lambda_K, levels=sort(unique(Ldat$Lambda_K)))
      k_labels <- setNames(paste0("k=",sort(unique(Ldat$Lambda_K))), as.character(sort(unique(Ldat$Lambda_K))))
      if ("1" %in% names(k_labels)) k_labels["1"] <- "k=1.0 (Baseline)"
      p5 <- ggplot2::ggplot(Ldat, ggplot2::aes(x=Value, fill=Lambda_K_Factor, color=Lambda_K_Factor))+
        ggplot2::geom_density(alpha=0.4)+
        ggplot2::facet_grid(Lambda_Type_Label ~ Duration_Label, scales="free_y", labeller=ggplot2::labeller(Lambda_Type_Label=ggplot2::label_parsed))+
        ggplot2::labs(x="Lambda Value", y="Density", fill="Scaling (k)", color="Scaling (k)")+
        ggplot2::scale_fill_manual(values=k_colors, labels=k_labels, drop=TRUE)+
        ggplot2::scale_color_manual(values=k_colors, labels=k_labels, drop=TRUE)+
        theme_science
      ggplot2::ggsave("Plot5_Lambda_K_Distributions.pdf", plot=p5, width=3+3*length(unique(Ldat$Duration)), height=7, device=dev_fun)
    }
  }
}
