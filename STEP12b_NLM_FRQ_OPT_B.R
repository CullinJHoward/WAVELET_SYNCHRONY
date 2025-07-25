# LIBRARY
library(ggplot2)
library(dplyr)
library(R.utils)
library(circular)

# Set working directory 
POUNDTOWN <- 0

if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\NONLINEAR_PATTERN_MODELING\\PHASE_ANGLE\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\NONLINEAR_PATTERN_MODELING\\PHASE_ANGLE\\'
}
setwd(work_dir)

## IDENTIFY THOSE FAILING QC 
QC_DF <- read.csv("..//..//WAVELET_ANALYSES//FRQ_OPT_CLASSIFICATION_VERSION_B_SUMMARY.csv")
QC_DF$ID <- sprintf("%04d", as.numeric(QC_DF$ID))
FAIL_ID <- QC_DF$ID[QC_DF$OPT_CRIT == "QC_FAIL"]


################################################################################
####################### SET UP THE ITERATIVELY RUN NONLINEAR MODEL 

plot_dir <- "OPT_PHASE_Vb_NLM_PLOTS"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

csv_dir <- "FRQ_OPT_VERSION_B_PHS_LONG"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No CSV files found in the specified directory.")
}

fit_nls <- function(data, start_vals) {
  tryCatch({
    model <- nls(
      PHS ~ a * cos(((2 * pi) / b) * (TIME - c)),
      data = data,
      start = start_vals,
      control = nls.control(maxiter = 50, tol = 1e-05, warnOnly = FALSE),
      na.action = na.omit
    )
    resid <- resid(model)
    rmse <- sqrt(mean(resid^2))
    list(
      model = model,
      rmse = rmse,
      aic = AIC(model),
      bic = BIC(model),
      logLik = as.numeric(logLik(model)),
      residuals = sum(resid^2),
      converged = TRUE,
      params = summary(model)$parameters
    )
  }, error = function(e) {
    cat(sprintf("Model fit error: %s\n", e$message))
    list(
      model = NULL,
      rmse = NA,
      aic = NA,
      bic = NA,
      logLik = NA,
      residuals = NA,
      converged = FALSE,
      params = NULL,
      error = e$message
    )
  })
}

summary_df <- data.frame(
  ID = character(),
  a_start = numeric(),
  b_start = numeric(),
  c_start = numeric(),
  logLik = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  RMSE = numeric(),
  R2 = numeric(),
  a_est = numeric(),
  b_est = numeric(),
  c_est = numeric(),
  a_pval = numeric(),
  b_pval = numeric(),
  c_pval = numeric(),
  approach = character(),
  circ_mean = numeric(),
  stringsAsFactors = FALSE
)

for (csv_file in csv_files) {
  id <- regmatches(basename(csv_file), regexpr("\\d{4}", basename(csv_file)))
  if (length(id) == 0) {
    warning(sprintf("Could not extract 4-digit ID from %s. Skipping.", basename(csv_file)))
    next
  }
  
  # Skip if ID is in FAIL_ID
  if (id %in% FAIL_ID) {
    cat(sprintf("ID %s failed QC. Skipping fitting, marking as QC_FAIL in summary.\n", id))
    qc_row <- as.list(rep("QC_FAIL", ncol(summary_df)))
    names(qc_row) <- names(summary_df)
    qc_row$ID <- id
    summary_df <- rbind(summary_df, qc_row)
    next
  }
  
  cat(sprintf("Processing ID %s...\n", id))
  data <- read.csv(csv_file)
  if (!all(c("PHS", "TIME", "Hz") %in% names(data))) {
    warning(sprintf("CSV %s missing PHS, TIME, or Hz columns. Skipping.", basename(csv_file)))
    next
  }
  
  # Calculate circular mean for this subject
  circ_mean_val <- mean.circular(circular(data$PHS, units="radians", template="none", modulo="2pi"), na.rm=TRUE)
  
  # --- 1. Non-condensed approach ---
  data_noncond <- data
  start_grid_nc <- expand.grid(
    a = c(-2, -1,-.2, -.1, .1, .2, 1, 2),
    b = c(80, 100, 120, 140, 160, 180),
    c = c(50, 75, 100)
  )
  results_nc <- lapply(1:nrow(start_grid_nc), function(i) {
    start_vals <- as.list(start_grid_nc[i, ])
    fit_result <- fit_nls(data_noncond, start_vals)
    if (fit_result$converged) {
      return(data.frame(
        a = start_grid_nc[i, "a"],
        b = start_grid_nc[i, "b"],
        c = start_grid_nc[i, "c"],
        logLik = fit_result$logLik,
        aic = fit_result$aic,
        bic = fit_result$bic,
        rmse = fit_result$rmse,
        residuals = fit_result$residuals,
        a_est = fit_result$params["a", "Estimate"],
        b_est = fit_result$params["b", "Estimate"],
        c_est = fit_result$params["c", "Estimate"],
        a_pval = fit_result$params["a", "Pr(>|t|)"],
        b_pval = fit_result$params["b", "Pr(>|t|)"],
        c_pval = fit_result$params["c", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      ))
    } else {
      NULL
    }
  })
  results_nc <- results_nc[!sapply(results_nc, is.null)]
  if (length(results_nc) == 0) {
    warning(sprintf("No models converged for ID %s (non-condensed). Skipping.", id))
    next
  }
  results_df_nc <- do.call(rbind, results_nc)
  results_df_nc <- results_df_nc %>% arrange(rmse)
  best_nc <- results_df_nc[1, ]
  best_model_nc <- fit_nls(data_noncond, as.list(best_nc[c("a", "b", "c")]))$model
  ss_res_nc <- sum(resid(best_model_nc)^2)
  ss_tot_nc <- sum((data_noncond$PHS - mean(data_noncond$PHS))^2)
  r2_nc <- 1 - ss_res_nc / ss_tot_nc
  
  # --- 2. Condensed approach (3 Hz bands by quantiles) ---
  Hz_vals <- sort(unique(data$Hz))
  if (length(Hz_vals) > 3) {
    Hz_q <- quantile(data$Hz, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    data$Hz_band <- cut(data$Hz,
                        breaks = Hz_q,
                        labels = c("Lower", "Middle", "Upper"),
                        include.lowest = TRUE, right = TRUE)
    data_cond <- data %>%
      group_by(TIME, Hz_band) %>%
      summarise(PHS = mean(PHS, na.rm = TRUE), Hz = mean(Hz, na.rm = TRUE), .groups = 'drop')
  } else {
    data$Hz_band <- as.character(data$Hz)
    data_cond <- data
  }
  start_grid_c <- expand.grid(
    a = c(-2, -1,-.2, -.1, .1, .2, 1, 2),
    b = c(80, 100, 120, 140, 160, 180),
    c = c(50, 75, 100)
  )
  results_c <- lapply(1:nrow(start_grid_c), function(i) {
    start_vals <- as.list(start_grid_c[i, ])
    fit_result <- fit_nls(data_cond, start_vals)
    if (fit_result$converged) {
      return(data.frame(
        a = start_grid_c[i, "a"],
        b = start_grid_c[i, "b"],
        c = start_grid_c[i, "c"],
        logLik = fit_result$logLik,
        aic = fit_result$aic,
        bic = fit_result$bic,
        rmse = fit_result$rmse,
        residuals = fit_result$residuals,
        a_est = fit_result$params["a", "Estimate"],
        b_est = fit_result$params["b", "Estimate"],
        c_est = fit_result$params["c", "Estimate"],
        a_pval = fit_result$params["a", "Pr(>|t|)"],
        b_pval = fit_result$params["b", "Pr(>|t|)"],
        c_pval = fit_result$params["c", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      ))
    } else {
      NULL
    }
  })
  results_c <- results_c[!sapply(results_c, is.null)]
  if (length(results_c) == 0) {
    warning(sprintf("No models converged for ID %s (condensed). Skipping.", id))
    next
  }
  results_df_c <- do.call(rbind, results_c)
  results_df_c <- results_df_c %>% arrange(rmse)
  best_c <- results_df_c[1, ]
  best_model_c <- fit_nls(data_cond, as.list(best_c[c("a", "b", "c")]))$model
  ss_res_c <- sum(resid(best_model_c)^2)
  ss_tot_c <- sum((data_cond$PHS - mean(data_cond$PHS))^2)
  r2_c <- 1 - ss_res_c / ss_tot_c
  
  # Debug print
  cat(sprintf("ID %s: Non-condensed RMSE = %.4f, Condensed RMSE = %.4f\n", id, best_nc$rmse, best_c$rmse))
  cat("Non-condensed best model params:\n")
  print(best_nc)
  cat("Condensed best model params:\n")
  print(best_c)
  
  # --- 3. Compare and select approach ---
  if (best_nc$rmse < best_c$rmse) {
    selected <- "non-condensed"
    best_result <- best_nc
    best_model <- best_model_nc
    r2 <- r2_nc
    plot_data <- data_noncond
    plot_data$Hz_band <- as.character(plot_data$Hz)
  } else {
    selected <- "condensed"
    best_result <- best_c
    best_model <- best_model_c
    r2 <- r2_c
    plot_data <- data_cond
  }
  plot_data$pred <- predict(best_model)
  
  summary_df <- rbind(summary_df, data.frame(
    ID = id,
    a_start = best_result$a,
    b_start = best_result$b,
    c_start = best_result$c,
    logLik = best_result$logLik,
    AIC = best_result$aic,
    BIC = best_result$bic,
    RMSE = best_result$rmse,
    R2 = r2,
    a_est = best_result$a_est,
    b_est = best_result$b_est,
    c_est = best_result$c_est,
    a_pval = best_result$a_pval,
    b_pval = best_result$b_pval,
    c_pval = best_result$c_pval,
    approach = selected,
    circ_mean = circ_mean_val,
    stringsAsFactors = FALSE
  ))
  
  if (!is.null(best_model)) {
    p <- ggplot(plot_data, aes(x = TIME, y = PHS)) +
      geom_point(aes(color = Hz), alpha = 0.5, size = 1.5) +
      geom_line(aes(y = pred, group = Hz_band), color = "black", size = 1) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = paste("NLS Fit for ID", id), 
           x = "TIME", 
           y = "PHASE", 
           color = "Hz") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      annotate("text", x = Inf, y = -Inf, label = sprintf("AIC: %.2f\nBIC: %.2f\nRMSE: %.4f\nR^2: %.3f\nApproach: %s", best_result$aic, best_result$bic, best_result$rmse, r2, selected), hjust = 1.1, vjust = -0.1, size = 4)
    
    pdf_file <- file.path(plot_dir, paste0("OPT_FRQ_PHS_Vb_NLM_PLOT_", id, ".pdf"))
    ggsave(pdf_file, plot = p, width = 8, height = 6)
  }
  
  cat(sprintf("Best fit for ID %s: a=%.3f, b=%.1f, c=%.1f, RMSE=%.4f, Approach: %s\n", 
              id, best_result$a, best_result$b, best_result$c, best_result$rmse, selected))
}

################################################################################
############################# REDUCE THE SUMMARY AND SAVE 

RED_SUMMARY <- summary_df %>%
  select(c(ID, a_est, b_est, c_est, R2, circ_mean)) %>%
  rename(
    FOb_PHSa = a_est, 
    FOb_PHSb = b_est,
    FOb_PHSc = c_est,
    FOb_PHSr2 = R2,
    FOb_PHSm = circ_mean
  )

write.csv(RED_SUMMARY, "..//FRQ_OPT_Vb_NLM_PHS_SUMMARY.csv", row.names = FALSE)
