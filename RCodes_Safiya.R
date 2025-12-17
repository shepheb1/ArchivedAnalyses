# n = c(50, 200, 500)
# beta = c(0, 0.5, 2)
# tau = c(0, 0.5, 2)
# Default: n = 50, beta = 0, tau = 0

library(rms)
library(dplyr)
library(purrr)
library(Metrics)

# Data Generating Process (DGP)
generate_data <- function(n = 50,
                          beta = 0,
                          tau = 0,
                          sd_y = 1,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- rnorm(n, mean = 0, sd = 1)
  tlogit <- -0.5 + 0.75 * X
  pA <- plogis(tlogit)
  A <- rbinom(n, size = 1, prob = pA)
  eps <- rnorm(n, mean = 0, sd = sd_y)
  Y0 <- exp(beta * X + eps)
  Y1 <- exp(beta * X + tau + eps)
  Y <- ifelse(A == 1, Y1, Y0)
  data.frame(X = X, A = A, Y = Y)
}

# Function to fit cumulative probability model, predict F and compute Theta
compute_theta <- function(family, formula, data, y_unique = NULL) {
  
  fit <- tryCatch(
    orm(formula, family = family, x = TRUE, y = TRUE, data = data),
    error = function(e) {
      message("Error fitting model (compute_theta): ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(fit)) return(NA_real_)
  
  # Auxiliary function to compute F(y|X)
  get_cdf <- function(fit_obj, newdata_for_predict, y_unique_vec) {
    
    if (!is.numeric(y_unique_vec) || length(y_unique_vec) == 0 || any(is.na(y_unique_vec))) {
      warning("Invalid y_unique_vec provided to get_cdf, returning NA.")
      return(rep(NA_real_, length(y_unique_vec)))
    }
    
    # Predict linear predictors
    lp_result <- tryCatch(
      predict(fit_obj, newdata = newdata_for_predict, se.fit = FALSE, type = "lp"),
      error = function(e) {
        message("Error in predict.orm (type='lp') within get_cdf: ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(lp_result) || any(!is.finite(lp_result))) {
      warning("Invalid or non-finite linear predictors from predict.orm, returning NA for CDF.")
      return(rep(NA_real_, length(y_unique_vec)))
    }
    
      X_for_ExProb <- tryCatch(
      predict(fit_obj, newdata = newdata_for_predict, type = "x"),
      error = function(e) {
        message("Error in predict.orm (type='x') within get_cdf: ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(X_for_ExProb) || any(!is.finite(X_for_ExProb))) {
      warning("Invalid or non-finite design matrix from predict.orm (type='x'), returning NA for CDF.")
      return(rep(NA_real_, length(y_unique_vec)))
    }
    # Ensure X_for_ExProb is a matrix for consistent behavior
    if (!is.matrix(X_for_ExProb)) {
      X_for_ExProb <- as.matrix(X_for_ExProb)
    }
    
    # Get the ExProb function from the fitted object
    d <- ExProb(fit_obj)
    
    # Call the ExProb function with linear predictors, y-values, AND the design matrix X
    fy_result <- tryCatch({
      d(lp_result, y = y_unique_vec, X = X_for_ExProb, conf.int = 0.95)
    }, error = function(e) {
      message("Error in ExProb call (d(lp_result, y=..., X=...)): ", conditionMessage(e))
      return(NULL)
    })
    
    if (is.null(fy_result) || is.null(fy_result$prob)) {
      warning("ExProb did not return valid probabilities, returning NA.")
      return(rep(NA_real_, length(y_unique_vec)))
    }
    
    F_cdf <- 1 - colMeans(fy_result$prob, na.rm = TRUE)
    return(F_cdf)
  }
  
  if (is.null(y_unique)) {
    y_unique <- sort(unique(data$Y))
  }
  
  if (!is.numeric(y_unique) || length(y_unique) < 2 || any(is.na(y_unique))) {
    warning("Not enough unique or valid Y values for compute_theta, returning NA.")
    return(NA_real_)
  }
  
  newdat <- data.frame(X = data$X, A = 1)
  newdat0 <- data.frame(X = data$X, A = 0)
  
  F1 <- get_cdf(fit, newdat, y_unique)
  F0 <- get_cdf(fit, newdat0, y_unique)
  
  if (any(is.na(F1)) || any(is.na(F0))) {
    warning("F1 or F0 contains NAs after get_cdf, cannot compute Theta. Returning NA.")
    return(NA_real_)
  }
  
  p1 <- diff(c(0, F1))
  p0 <- diff(c(0, F0))
  
  if (length(p1) != length(y_unique)) {
    warning("Length mismatch for p1 and y_unique. This suggests F1 might have issues. Returning NA.")
    return(NA_real_)
  }
  if (length(p0) != length(y_unique)) {
    warning("Length mismatch for p0 and y_unique. This suggests F0 might have issues. Returning NA.")
    return(NA_real_)
  }
  
  h_val <- function(y1_val, y0_val) (y1_val > y0_val) + 0.5 * (y1_val == y0_val)
  h_matrix <- outer(y_unique, y_unique, Vectorize(h_val))
  
  p1_col_vec <- matrix(p1, ncol = 1)
  p0_row_vec <- matrix(p0, nrow = 1)
  
  if (length(p1) == nrow(h_matrix) && length(p0) == ncol(h_matrix)) {
    Theta <- sum( (p1_col_vec %*% p0_row_vec) * h_matrix )
  } else {
    warning("Dimensions for Theta calculation are not compatible. Returning NA.")
    Theta <- NA_real_
  }
  
  return(Theta)
}

# Parametric functions 
compute_parametric_log <- function(data, A_col, X_col, Y_col, n, by_interval = 0.01) {
  data$Y_transformed <- log(data[[Y_col]])
  formula_str <- paste0("Y_transformed ~ ", A_col, " + ", X_col)
  
  fit <- tryCatch({ lm(as.formula(formula_str), data = data) }, error = function(e) {
    message("Error in compute_parametric_log: ", e$message); return(NULL) })
  if (is.null(fit)) return(NA_real_)
  
  if (is.null(fit$residuals) || length(fit$residuals) < 2) {
    warning("Not enough residuals for sd() in compute_parametric_log. Returning NA."); return(NA_real_) }
  sigma <- sd(fit$residuals)
  if (is.na(sigma) || sigma == 0) {
    warning("Residual standard deviation is NA or zero in compute_parametric_log. Returning NA."); return(NA_real_) }
  
  ys_min <- min(data$Y_transformed, na.rm = TRUE)
  ys_max <- max(data$Y_transformed, na.rm = TRUE)
  if (!is.finite(ys_min) || !is.finite(ys_max)) {
    warning("Infinite or NA values in Y_transformed range. Returning NA for parametric log."); return(NA_real_) }
  
  ys <- seq(from = ys_min - 0.5 * sd(data$Y_transformed, na.rm=TRUE),
            to = ys_max + 0.5 * sd(data$Y_transformed, na.rm=TRUE), by = by_interval)
  if (length(ys) == 0) { warning("Generated 'ys' sequence is empty. Returning NA for parametric log."); return(NA_real_) }
  
  p1 <- matrix(NA, nrow = n, ncol = length(ys))
  p0 <- matrix(NA, nrow = n, ncol = length(ys))
  newdat_base <- data %>% select(all_of(c(X_col, A_col)))
  
  for (o in 1:n) {
    new1 <- newdat_base[o, , drop = FALSE]; new1[[A_col]] <- 1
    new0 <- newdat_base[o, , drop = FALSE]; new0[[A_col]] <- 0
    pred1 <- tryCatch(predict(fit, newdata = new1), error = function(e) { warning(paste("Predict error for new1:", e$message)); return(NA_real_) })
    pred0 <- tryCatch(predict(fit, newdata = new0), error = function(e) { warning(paste("Predict error for new0:", e$message)); return(NA_real_) })
    if (is.na(pred1) || is.na(pred0)) { p1[o, ] <- NA_real_; p0[o, ] <- NA_real_ } else {
      p1[o, ] <- dnorm(ys, pred1, sigma); p0[o, ] <- dnorm(ys, pred0, sigma) }
  }
  
  pdf1 <- colMeans(p1, na.rm = TRUE); pdf0 <- colMeans(p0, na.rm = TRUE)
  if (any(is.na(pdf1)) || any(is.na(pdf0))) {
    warning("PDFs contain NAs after colMeans. Returning NA for parametric log."); return(NA_real_) }
  
  h.matrix <- matrix(0, nrow = length(ys), ncol = length(ys)); diag(h.matrix) <- 0.5; h.matrix[lower.tri(h.matrix)] <- 1
  if (length(pdf1) == 0 || length(pdf0) == 0 || any(!is.finite(pdf1)) || any(!is.finite(pdf0))) {
    warning("PDFs are invalid for Theta calculation. Returning NA."); return(NA_real_) }
  
  Theta <- sum(pdf1 %*% h.matrix %*% pdf0) * (by_interval^2); return(Theta)
}

compute_parametric_sqrt <- function(data, A_col, X_col, Y_col, n, by_interval = 0.01) {
  data$Y_transformed <- sqrt(data[[Y_col]])
  formula_str <- paste0("Y_transformed ~ ", A_col, " + ", X_col)
  
  fit <- tryCatch({ lm(as.formula(formula_str), data = data) }, error = function(e) {
    message("Error in compute_parametric_sqrt: ", e$message); return(NULL) })
  if (is.null(fit)) return(NA_real_)
  
  if (is.null(fit$residuals) || length(fit$residuals) < 2) {
    warning("Not enough residuals for sd() in compute_parametric_sqrt. Returning NA."); return(NA_real_) }
  sigma <- sd(fit$residuals)
  if (is.na(sigma) || sigma == 0) {
    warning("Residual standard deviation is NA or zero in compute_parametric_sqrt. Returning NA."); return(NA_real_) }
  
  ys_min <- min(data$Y_transformed, na.rm = TRUE)
  ys_max <- max(data$Y_transformed, na.rm = TRUE)
  if (!is.finite(ys_min) || !is.finite(ys_max)) {
    warning("Infinite or NA values in Y_transformed range. Returning NA for parametric sqrt."); return(NA_real_) }
  
  ys <- seq(from = ys_min - 0.5 * sd(data$Y_transformed, na.rm=TRUE),
            to = ys_max + 0.5 * sd(data$Y_transformed, na.rm=TRUE), by = by_interval)
  if (length(ys) == 0) { warning("Generated 'ys' sequence is empty. Returning NA for parametric sqrt."); return(NA_real_) }
  
  p1 <- matrix(NA, nrow = n, ncol = length(ys))
  p0 <- matrix(NA, nrow = n, ncol = length(ys))
  newdat_base <- data %>% select(all_of(c(X_col, A_col)))
  
  for (o in 1:n) {
    new1 <- newdat_base[o, , drop = FALSE]; new1[[A_col]] <- 1
    new0 <- newdat_base[o, , drop = FALSE]; new0[[A_col]] <- 0
    pred1 <- tryCatch(predict(fit, newdata = new1), error = function(e) { warning(paste("Predict error for new1:", e$message)); return(NA_real_) })
    pred0 <- tryCatch(predict(fit, newdata = new0), error = function(e) { warning(paste("Predict error for new0:", e$message)); return(NA_real_) })
    if (is.na(pred1) || is.na(pred0)) { p1[o, ] <- NA_real_; p0[o, ] <- NA_real_ } else {
      p1[o, ] <- dnorm(ys, pred1, sigma); p0[o, ] <- dnorm(ys, pred0, sigma) }
  }
  
  pdf1 <- colMeans(p1, na.rm = TRUE); pdf0 <- colMeans(p0, na.rm = TRUE)
  if (any(is.na(pdf1)) || any(is.na(pdf0))) {
    warning("PDFs contain NAs after colMeans. Returning NA for parametric sqrt."); return(NA_real_) }
  
  h.matrix <- matrix(0, nrow = length(ys), ncol = length(ys)); diag(h.matrix) <- 0.5; h.matrix[lower.tri(h.matrix)] <- 1
  if (length(pdf1) == 0 || length(pdf0) == 0 || any(!is.finite(pdf1)) || any(!is.finite(pdf0))) {
    warning("PDFs are invalid for Theta calculation. Returning NA."); return(NA_real_) }
  
  Theta <- sum(pdf1 %*% h.matrix %*% pdf0) * (by_interval^2); return(Theta)
}


# Main Simulation Function 
run_one_simulation_replication <- function(n, beta, tau, sd_y, B_boot, seed_rep) {
  set.seed(seed_rep)
  
  data <- generate_data(n = n, beta = beta, tau = tau, sd_y = sd_y, seed = seed_rep)
  
  dd_local <- datadist(data)
  old_options_datadist <- options(datadist = "dd_local")
  on.exit({
    options(old_options_datadist)
    if (exists("dd_local", envir = .GlobalEnv)) {
      rm(dd_local, envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign("dd_local", dd_local, envir = .GlobalEnv)
  
  Theta_hat1 <- compute_theta(family = 'probit', formula = Y ~ A + X, data = data)
  Theta_hat2 <- compute_theta(family = 'logistic', formula = Y ~ A + X, data = data)
  Theta_hat3 <- compute_theta(family = 'probit', formula = Y ~ A, data = data)
  
  Theta_hat_parm_log <- compute_parametric_log(data = data, A_col = "A", X_col = "X", Y_col = "Y", n = n)
  Theta_hat_parm_sqrt <- compute_parametric_sqrt(data = data, A_col = "A", X_col = "X", Y_col = "Y", n = n)
  
  boot_theta1_ests <- numeric(B_boot)
  boot_theta2_ests <- numeric(B_boot)
  boot_theta3_ests <- numeric(B_boot)
  boot_parm1_ests <- numeric(B_boot)
  boot_parm2_ests <- numeric(B_boot)
  
  for (b in 1:B_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- data[idx, , drop = FALSE]
    
    boot_dd_local <- datadist(dat_b)
    assign("dd_local", boot_dd_local, envir = .GlobalEnv)
    
    y_unique_boot <- sort(unique(dat_b$Y))
    
    boot_theta1_ests[b] <- compute_theta(family = 'probit', formula = Y ~ A + X, data = dat_b, y_unique = y_unique_boot)
    boot_theta2_ests[b] <- compute_theta(family = 'logistic', formula = Y ~ A + X, data = dat_b, y_unique = y_unique_boot)
    boot_theta3_ests[b] <- compute_theta(family = 'probit', formula = Y ~ A, data = dat_b, y_unique = y_unique_boot)
    boot_parm1_ests[b] <- compute_parametric_log(data = dat_b, A_col = "A", X_col = "X", Y_col = "Y", n = n)
    boot_parm2_ests[b] <- compute_parametric_sqrt(data = dat_b, A_col = "A", X_col = "X", Y_col = "Y", n = n)
  }
  
  calc_ci <- function(point_est, boot_ests) {
    boot_ests_clean <- boot_ests[!is.na(boot_ests)]
    if (length(boot_ests_clean) < 2) { return(c(NA, NA)) }
    se_boot <- sd(boot_ests_clean, na.rm = TRUE)
    if (is.na(se_boot) || se_boot == 0) { return(c(NA, NA)) }
    lower <- point_est - 1.96 * se_boot
    upper <- point_est + 1.96 * se_boot
    c(lower, upper)
  }
  
  ci1 <- calc_ci(Theta_hat1, boot_theta1_ests)
  ci2 <- calc_ci(Theta_hat2, boot_theta2_ests)
  ci3 <- calc_ci(Theta_hat3, boot_theta3_ests)
  ci_parm1 <- calc_ci(Theta_hat_parm_log, boot_parm1_ests)
  ci_parm2 <- calc_ci(Theta_hat_parm_sqrt, boot_parm2_ests)
  
  list(
    n = n, beta = beta, tau = tau, sd_y = sd_y, seed_rep = seed_rep,
    Theta_hat1 = Theta_hat1, Theta_hat2 = Theta_hat2, Theta_hat3 = Theta_hat3,
    Theta_hat_parm_log = Theta_hat_parm_log, Theta_hat_parm_sqrt = Theta_hat_parm_sqrt,
    ci_lower_theta1 = ci1[1], ci_upper_theta1 = ci1[2],
    ci_lower_theta2 = ci2[1], ci_upper_theta2 = ci2[2],
    ci_lower_theta3 = ci3[1], ci_upper_theta3 = ci3[2],
    ci_lower_parm1 = ci_parm1[1], ci_upper_parm1 = ci_parm1[2],
    ci_lower_parm2 = ci_parm2[1], ci_upper_parm2 = ci_parm2[2]
  )
}

# Wrapper for the full Simulation Study 
run_simulation_study <- function(n_sims = 1000,
                                 n = 50,
                                 beta = 0,
                                 tau = 0,
                                 sd_y = 1,
                                 B_boot = 200,
                                 master_seed = 123) {
  set.seed(master_seed)
  seeds_for_reps <- sample(1:1e6, n_sims)
  
  results_list <- purrr::map(1:n_sims, function(i) {
    message("Running simulation replication ", i, " of ", n_sims, "...")
    run_one_simulation_replication(n = n, beta = beta, tau = tau,
                                   sd_y = sd_y, B_boot = B_boot,
                                   seed_rep = seeds_for_reps[i])
  })
  
  results_df <- dplyr::bind_rows(results_list)
  
  return(results_df)
}

### Function to Compute Empirical True Theta

compute_true_theta <- function(beta_val, tau_val, n_truth = 1e6, sd_y = 1, seed_truth = 12345) {
  set.seed(seed_truth)
  x <- rnorm(n_truth, 0, 1) # Confounder X
  
  # Generate exposure variable A
  tlogit <- -0.5 + 0.75 * x
  pA <- plogis(tlogit)
  A <- rbinom(n_truth, 1, pA)
  
  # Outcome under A=1 and A=0
  y1 <- exp(beta_val * x + tau_val + rnorm(n_truth, 0, sd_y))
  y0 <- exp(beta_val * x + rnorm(n_truth, 0, sd_y))
  
  # h function for concordance
  h_func <- function(y1_obs, y0_obs) {
    as.numeric(y1_obs > y0_obs) + 0.5 * as.numeric(y1_obs == y0_obs)
  }
  # For each pair of independent observations, apply h function
  
  true_theta_val <- mean(h_func(y1, y0), na.rm = TRUE)
  
  return(true_theta_val)
}

### Function to Analyze Simulation Results for a Single Scenario

analyze_scenario_results <- function(sim_results_df, true_theta_val) {
  
  # Ensure true_theta_val is replicated for calculations
  true_theta_vec <- rep(true_theta_val, nrow(sim_results_df))
  
  metrics_list <- list()
  
  # List of estimators to analyze
  estimators <- c("Theta_hat1", "Theta_hat2", "Theta_hat3")
  
  for (est_name in estimators) {
    est_values <- sim_results_df[[est_name]]
    ci_lower_name <- paste0("ci_lower_", gsub("Theta_hat", "theta", est_name))
    ci_upper_name <- paste0("ci_upper_", gsub("Theta_hat", "theta", est_name))
    
    ci_lower_values <- sim_results_df[[ci_lower_name]]
    ci_upper_values <- sim_results_df[[ci_upper_name]]
    
    # Remove NAs from calculations
    valid_indices <- which(!is.na(est_values) & !is.na(ci_lower_values) & !is.na(ci_upper_values))
    
    est_values_clean <- est_values[valid_indices]
    ci_lower_values_clean <- ci_lower_values[valid_indices]
    ci_upper_values_clean <- ci_upper_values[valid_indices]
    true_theta_vec_clean <- true_theta_vec[valid_indices]
    
    if (length(est_values_clean) < 2) { 
      metrics_list[[est_name]] <- list(
        Bias = NA,
        SD = NA,
        RMSE = NA,
        CP = NA
      )
      next
    }
    
    bias <- mean(est_values_clean - true_theta_vec_clean, na.rm = TRUE)
    sd_est <- sd(est_values_clean, na.rm = TRUE)
    rmse_est <- Metrics::rmse(actual = true_theta_vec_clean, predicted = est_values_clean)
    
    coverage <- ifelse(ci_lower_values_clean < true_theta_vec_clean & ci_upper_values_clean > true_theta_vec_clean, 1, 0)
    cp <- mean(coverage, na.rm = TRUE)
    
    metrics_list[[est_name]] <- list(
      Bias = bias,
      SD = sd_est,
      RMSE = rmse_est,
      CP = cp
    )
  }
  
  # Flatten the list into a data frame row
  result_row <- list()
  for (est_name in estimators) {
    result_row[[paste0(est_name, "_Bias")]] <- metrics_list[[est_name]]$Bias
    result_row[[paste0(est_name, "_SD")]] <- metrics_list[[est_name]]$SD
    result_row[[paste0(est_name, "_RMSE")]] <- metrics_list[[est_name]]$RMSE
    result_row[[paste0(est_name, "_CP")]] <- metrics_list[[est_name]]$CP
  }
  
  return(as.data.frame(result_row))
}


### Define Scenarios and Loop for Analysis

# Define the scenarios
n_values <- c(50, 200, 500)
beta_values <- c(0, 0.5, 2)
tau_values <- c(0, 0.5, 2)
sd_y_value <- 1 # Assuming sd_y is fixed at 1
n_sims_per_scenario <- 1000 # Number of simulation replications for each scenario
B_boot_value <- 200 # Number of bootstrap samples

# Create a data frame of all combinations of scenarios
scenarios <- expand.grid(
  n = n_values,
  beta = beta_values,
  tau = tau_values
)

all_analysis_results <- list()

message("\n--- Starting Full Simulation Study Across Scenarios ---")

for (i in 1:nrow(scenarios)) {
  current_n <- scenarios$n[i]
  current_beta <- scenarios$beta[i]
  current_tau <- scenarios$tau[i]
  
  message(paste0("\nProcessing Scenario: n=", current_n, ", beta=", current_beta, ", tau=", current_tau))
  
  # Compute the empirical true theta for the current scenario
  true_theta <- compute_true_theta(
    beta_val = current_beta,
    tau_val = current_tau,
    n_truth = 1e6, # Large sample for true theta
    sd_y = sd_y_value,
    seed_truth = 12345 + i # Vary seed for true theta calc slightly
  )
  message("  Empirical True Theta: ", round(true_theta, 4))
  
  # Run the main simulation study for this scenario
  scenario_sim_results <- run_simulation_study(
    n_sims = n_sims_per_scenario,
    n = current_n,
    beta = current_beta,
    tau = current_tau,
    sd_y = sd_y_value,
    B_boot = B_boot_value,
    master_seed = 1000 + i # Vary master seed for replications
  )
  
  # Analyze the results for this scenario
  scenario_analysis <- analyze_scenario_results(scenario_sim_results, true_theta)
  
  # Add scenario parameters to the analysis results
  scenario_analysis$n <- current_n
  scenario_analysis$beta <- current_beta
  scenario_analysis$tau <- current_tau
  scenario_analysis$True_Theta <- true_theta # Add true theta to the table
  
  all_analysis_results[[i]] <- scenario_analysis
}

### Final Table Generation

final_results_table <- dplyr::bind_rows(all_analysis_results)

# Optional: Reorder columns for better readability
final_results_table <- final_results_table %>%
  select(n, beta, tau, True_Theta, everything())

print(final_results_table)

message("\n--- Full Simulation Analysis Complete ---")
