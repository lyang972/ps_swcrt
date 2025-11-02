# 1_preprocess.R

# --- 0. Load Configuration ---
# Configuration is expected to be loaded by 0_driver.R
if (!exists("CONFIG_PATHS") || !exists("data_path")) {
  message("Configuration or data path not found. Ensure 0_driver.R is executed.")
  # Stop execution if sourced
  if (sys.nframe() > 0) return(invisible(NULL))
}

# Load Data
if (!is.null(data_path) && file.exists(data_path)) {
  tryCatch({
    tb <- as.data.frame(readxl::read_xlsx(data_path, sheet = 1))
  }, error = function(e) {
    stop(paste("Error reading Excel file:", e$message))
  })
} else {
  message("Data file not found. Skipping preprocessing.")
  if (sys.nframe() > 0) return(invisible(NULL))
}

# Proceed only if data (tb) is loaded
if (exists("tb")) {
  
  # --- 1. Data Transformation ---
  # (Column definitions identical to the original provided code)
  column_question =  c(
    "ID", "Q3", "Q15", "Q17", "Q18", "Q20", "Q1", # Covariates and Location
    "Q144", "Q145", "Q146", "Q147", "Q148", "Q149", "Q150", "Q75", # Baseline M and Y
    "F1Q50", "F1Q51", "F1Q52", "F1Q53", "F1Q54", "F1Q55", "F1Q56", "F1Q7", "F1Q12",
    "F2Q52", "F2Q53", "F2Q54", "F2Q55", "F2Q56", "F2Q57", "F2Q58", "F2Q7", "F2Q12",
    "F3Q67", "F3Q68", "F3Q69", "F3Q70", "F3Q71", "F3Q72", "F3Q73", "F3Q7", "F3Q12",
    "F4Q85", "F4Q86", "F4Q87", "F4Q88", "F4Q89", "F4Q90", "F4Q91", "F4Q7", "F4Q12"
  )
  
  column_name = c(
    "ID", "age", "marital", "education", "income", "orientation", "F0_location",
    "F0_sn1", "F0_sn2", "F0_sn3", "F0_sn4", "F0_sn5", "F0_sn6", "F0_sn7", "F0_HIV_test",
    "F1_sn1", "F1_sn2", "F1_sn3", "F1_sn4", "F1_sn5", "F1_sn6", "F1_sn7", "F1_HIV_test_facility", "F1_HIV_test_self",
    "F2_sn1", "F2_sn2", "F2_sn3", "F2_sn4", "F2_sn5", "F2_sn6", "F2_sn7", "F2_HIV_test_facility", "F2_HIV_test_self",
    "F3_sn1", "F3_sn2", "F3_sn3", "F3_sn4", "F3_sn5", "F3_sn6", "F3_sn7", "F3_HIV_test_facility", "F3_HIV_test_self",
    "F4_sn1", "F4_sn2", "F4_sn3", "F4_sn4", "F4_sn5", "F4_sn6", "F4_sn7", "F4_HIV_test_facility", "F4_HIV_test_self"
  )
  
  if (!all(column_question %in% colnames(tb))) {
    stop("Missing required columns in data.")
  }
  
  data = tb[, column_question]; names(data) = column_name
  data$orientation = ifelse(data$orientation == 4, 3, data$orientation)
  
  # Transform wide to long format
  df = data %>%
    pivot_longer(cols = starts_with("F"),
                 names_to = c('time_idx_R', '.value'),
                 names_pattern = "F([^_]+)_(.*)")
  
  # --- 2. Variable Definition and Cleaning ---
  df = df %>%
    mutate(
      ID = as.character(ID), time_idx_R = as.integer(time_idx_R), time = time_idx_R + 1,
      # Mediator M (Social Norm Calculation)
      sn_raw = rowSums(across(c(sn1, sn2, sn3, sn4, sn5, sn6)), na.rm = FALSE) + (5 - sn7),
      # Set SN to NA if any component is NA
      sn = ifelse(if_any(c(sn1, sn2, sn3, sn4, sn5, sn6, sn7), is.na), NA, sn_raw),
      # Outcome Y (HIV Testing)
      Y_facility = case_when(HIV_test_facility == 1 ~ 1, HIV_test_facility == 2 ~ 0, TRUE ~ NA_real_),
      Y_self = case_when(HIV_test_self == 1 ~ 1, HIV_test_self == 2 ~ 0, TRUE ~ NA_real_),
      Y_baseline = case_when(HIV_test == 1 ~ 1, HIV_test == 2 ~ 0, TRUE ~ NA_real_),
      # Combine testing sources
      Y = ifelse(time > 1, pmax(Y_facility, Y_self, na.rm = FALSE), Y_baseline)
    )
  
  # Covariates (X_ij)
  df = df %>%
    mutate(
      age_raw = as.integer(age),
      Marital = as.factor(marital), Education = as.factor(education),
      Income = as.factor(income), Orientation = as.factor(orientation)
    )
  
  # --- 3. SW-CRT Design and Covariates (C_j) ---
  start_time_map = c("1"=2, "2"=2, "3"=3, "4"=3, "5"=4, "6"=4, "7"=5, "8"=5)
  # Province map: GD (Guangdong) and SD (Shandong)
  province_map = c("1"="GD", "2"="SD", "3"="GD", "4"="SD", "5"="GD", "6"="SD", "7"="GD", "8"="SD")
  
  df = df %>%
    group_by(ID) %>%
    # 'location' was created by pivot_longer from F0_location
    mutate(Cluster_ID = ifelse(time_idx_R == 0, location, NA)) %>%
    tidyr::fill(Cluster_ID, .direction = "downup") %>%
    ungroup() %>%
    mutate(
      Cluster_ID_Str = as.character(Cluster_ID),
      Start_Time = start_time_map[Cluster_ID_Str],
      Province = factor(province_map[Cluster_ID_Str]), # Cj
      Z = ifelse(time >= Start_Time, 1, 0)
    ) %>%
    group_by(ID) %>% arrange(time) %>% mutate(Duration = cumsum(Z)) %>% ungroup()
  
  # --- 4. Data Filtering and Standardization ---
  data_model = df %>%
    filter(!is.na(Y) & !is.na(sn)) %>% filter(time > 1) %>%
    filter(!is.na(age_raw) & !is.na(Marital) & !is.na(Education) & !is.na(Income) & !is.na(Orientation))
  
  if (nrow(data_model) == 0) {
    stop("No data remaining after filtering.")
  }
  
  T_max = max(data_model$time); D_max = max(data_model$Duration)
  
  # Standardize M
  sn_mean = mean(data_model$sn); sn_sd = sd(data_model$sn)
  data_model$M_std = (data_model$sn - sn_mean) / sn_sd
  
  # Center Age
  age_mean = mean(data_model$age_raw, na.rm = TRUE)
  data_model$Age_centered = data_model$age_raw - age_mean
  
  # Create numeric indices for Stan
  data_model = data_model %>%
    mutate(i_idx = as.numeric(fct_inorder(ID)), j_idx = as.numeric(fct_inorder(as.character(Cluster_ID))))
  
  N_individuals = max(data_model$i_idx); N_clusters = max(data_model$j_idx)
  
  Individual_to_Cluster_Map = data_model %>%
    select(i_idx, j_idx, Province) %>% distinct() %>% arrange(i_idx)
  
  Provinces_List = sort(unique(Individual_to_Cluster_Map$Province))
  
  # --- 5. Prepare Covariate Matrices (X_ij and C_j) ---
  X_data = data_model %>% select(i_idx, Age_centered, Marital, Education, Income, Orientation) %>% distinct() %>% arrange(i_idx)
  
  # Simplified matrix creation using model.matrix
  if (nrow(X_data) > 0) {
    X_matrix_full = model.matrix(~ Age_centered + Marital + Education + Income + Orientation, data = X_data)
    if (ncol(X_matrix_full) > 1) {
      X_matrix = X_matrix_full[, -1, drop = FALSE] # Remove intercept
    } else {
      X_matrix = matrix(0, nrow = nrow(X_data), ncol = 0)
    }
  } else {
    X_matrix = matrix(0, nrow = 0, ncol = 0)
  }
  K_X = ncol(X_matrix)
  
  
  C_data = data_model %>% select(j_idx, Province) %>% distinct() %>% arrange(j_idx)
  
  if (nrow(C_data) > 0) {
    C_matrix_full = model.matrix(~ Province, data = C_data)
    if (ncol(C_matrix_full) > 1) {
      C_matrix = C_matrix_full[, -1, drop = FALSE] # Remove intercept
    } else {
      C_matrix = matrix(0, nrow = nrow(C_data), ncol = 0)
    }
  } else {
    C_matrix = matrix(0, nrow = 0, ncol = 0)
  }
  K_C = ncol(C_matrix)
  
  
  # --- 6. Prepare Stan Data List ---
  stan_data = list(
    N = nrow(data_model), I = N_individuals, J = N_clusters, T_max = T_max, D_max = D_max,
    K_X = K_X, K_C = K_C, M = data_model$M_std, Y = as.integer(data_model$Y),
    ii = data_model$i_idx, jj = data_model$j_idx, time_idx = data_model$time, Duration = data_model$Duration,
    X = X_matrix, C = C_matrix
  )
  
  # Save processed data
  save(stan_data, data_model, X_matrix, C_matrix, Individual_to_Cluster_Map, Provinces_List,
       sn_mean, sn_sd, T_max, D_max, K_X, K_C, N_individuals,
       file = CONFIG_PATHS$PROCESSED_DATA)
  
  message("Preprocessing finished.")
  message(sprintf("Summary: N_obs=%d, N_indiv=%d, N_clust=%d, T_max=%d, D_max=%d",
                  nrow(data_model), N_individuals, N_clusters, T_max, D_max))
}