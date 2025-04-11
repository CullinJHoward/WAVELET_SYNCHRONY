### LIBRARY
library(ggplot2)
library(dplyr)
library(R.utils) 

# Set working directory 

POUNDTOWN <- 0

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\SRCD_25_ANALYSIS\\WAVELET\\'
} else {
  work_dir <- 'F:\\SRCD_25_ANALYSIS\\WAVELET\\'
}

setwd(work_dir)

################################################################################
####################### SET UP THE ITERATIVELY RUN NONLINEAR MODEL 

# Create output directory for plots if it doesn't exist
plot_dir <- "OPT_PHASE_NLME_PLOTS"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Define the folder containing CSV files
csv_dir <- "F:/SRCD_25_ANALYSIS/WAVELET/LONG_OPTIMIZED_PHASE_SERIES"

# List all CSV files
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No CSV files found in the specified directory.")
}

# Define starting value grid
start_grid <- expand.grid(
  a = c(1, 2, 3),                # Amplitude: favor higher
  b = c(80, 100, 120, 140, 160, 180),  # Period: favor lower
  c = c(50, 75, 100)             # Phase shift: favor smaller
)

# Function to fit nls model and return metrics
fit_nls <- function(data, start_vals) {
  tryCatch({
    model <- nls(
      PHASE ~ a * cos(((2 * pi) / b) * (TIME - c)),
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

# Initialize summary data frame
summary_df <- data.frame(
  ID = character(),
  a_start = numeric(),
  b_start = numeric(),
  c_start = numeric(),
  logLik = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  RMSE = numeric(),
  a_est = numeric(),
  b_est = numeric(),
  c_est = numeric(),
  a_pval = numeric(),
  b_pval = numeric(),
  c_pval = numeric(),
  stringsAsFactors = FALSE
)

# Process each CSV file
for (csv_file in csv_files) {
  # Extract ID from filename
  id <- regmatches(basename(csv_file), regexpr("\\d{4}", basename(csv_file)))
  if (length(id) == 0) {
    warning(sprintf("Could not extract 4-digit ID from %s. Skipping.", basename(csv_file)))
    next
  }
  
  cat(sprintf("Processing ID %s...\n", id))
  
  # Load the CSV
  data <- read.csv(csv_file)
  
  # Ensure required columns exist
  if (!all(c("PHASE", "TIME", "PER_CEN") %in% names(data))) {
    warning(sprintf("CSV %s missing PHASE, TIME, or PER_CEN columns. Skipping.", basename(csv_file)))
    next
  }
  
  # Test each combination of starting values
  results <- lapply(1:nrow(start_grid), function(i) {
    start_vals <- as.list(start_grid[i, ])
    fit_result <- fit_nls(data, start_vals)
    if (fit_result$converged) {
      return(data.frame(
        a = start_grid[i, "a"],
        b = start_grid[i, "b"],
        c = start_grid[i, "c"],
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
  
  # Remove NULL entries (failed fits)
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) {
    warning(sprintf("No models converged for ID %s. Skipping.", id))
    next
  }
  
  # Convert to data frame
  results_df <- do.call(rbind, results)
  
  # Rank models: lower RMSE, AIC, BIC, residuals; higher logLik
  results_df <- results_df %>%
    mutate(
      rmse_rank = 1 - (rmse - min(rmse)) / (max(rmse) - min(rmse)),
      aic_rank = 1 - (aic - min(aic)) / (max(aic) - min(aic)),
      bic_rank = 1 - (bic - min(bic)) / (max(bic) - min(bic)),
      logLik_rank = (logLik - min(logLik)) / (max(logLik) - min(logLik)),
      resid_rank = 1 - (residuals - min(residuals)) / (max(residuals) - min(residuals)),
      total_score = 0.3 * rmse_rank + 0.2 * aic_rank + 0.2 * bic_rank + 0.2 * logLik_rank + 0.1 * resid_rank
    ) %>%
    arrange(desc(total_score), desc(a), b, c)
  
  # Select best model
  best_result <- results_df[1, ]
  best_start <- list(a = best_result$a, b = best_result$b, c = best_result$c)
  best_model <- fit_nls(data, best_start)$model
  
  # Add to summary data frame
  summary_df <- rbind(summary_df, data.frame(
    ID = id,
    a_start = best_result$a,
    b_start = best_result$b,
    c_start = best_result$c,
    logLik = best_result$logLik,
    AIC = best_result$aic,
    BIC = best_result$bic,
    RMSE = best_result$rmse,
    a_est = best_result$a_est,
    b_est = best_result$b_est,
    c_est = best_result$c_est,
    a_pval = best_result$a_pval,
    b_pval = best_result$b_pval,
    c_pval = best_result$c_pval
  ))
  
  # Generate and save plot with PER_CEN coloring
  data$pred <- predict(best_model)
  # Create a categorical variable for PER_CEN
  data$PER_CEN_cat <- ifelse(data$PER_CEN == 0, "Mean (0)",
                             ifelse(data$PER_CEN > 0, "Above Mean (>0)", "Below Mean (<0)"))
  data$PER_CEN_cat <- factor(data$PER_CEN_cat, 
                             levels = c("Mean (0)", "Above Mean (>0)", "Below Mean (<0)"))
  
  p <- ggplot(data, aes(x = TIME, y = PHASE)) +
    geom_point(aes(color = PER_CEN_cat), alpha = 0.5, size = 1.5) +
    geom_line(aes(y = pred), color = "black", size = 1) +
    scale_color_manual(values = c("Mean (0)" = "gray50", 
                                  "Above Mean (>0)" = "red", 
                                  "Below Mean (<0)" = "blue")) +
    labs(title = paste("NLS Fit for ID", id), 
         x = "TIME", 
         y = "PHASE", 
         color = "PER_CEN Category") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  pdf_file <- file.path(plot_dir, paste0("OPT_PHS_NLME_PLOT_", id, ".pdf"))
  ggsave(pdf_file, plot = p, width = 8, height = 6)
  
  cat(sprintf("Best fit for ID %s: a=%.1f, b=%.1f, c=%.1f, RMSE=%.4f\n", 
              id, best_result$a, best_result$b, best_result$c, best_result$rmse))
}


# Save summary data frame
write.csv(summary_df, "OPT_NLME_PHASE_SUMMARY.csv", row.names = FALSE)
cat("Summary saved to OPT_NLME_PHASE_SUMMARY.csv\n")
