
### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)
library(slider)
library(tidyr)
library(stringr)
library(circular)
library(mgcv)
library(purrr)

# Set working directory 

POUNDTOWN <- 0

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\ANGLE_FILTERED\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\ANGLE_FILTERED\\'
}

setwd(work_dir)

################################################################################
############################ LOOP THROUGH EACH ANGLE DF AND SUMMARIZE

# LIST ALL FILES IN THE WORKING DIRECTORY
csv_files <- list.files(pattern = "\\.csv$")

# INITIALIZE LISTS FOR TRACKING
summary_list_PHS <- list()
successful_ids_PHS <- character()
failed_ids_PHS <- character()

# LOOP THROUGH ALL CSVs
for (file in csv_files) {
  # EXTRACT PARTICIPANT IDs
  id <- sub(".*ANGLE_([0-9]+)_.*", "\\1", file)
  id <- sprintf("%04d", as.numeric(id)) # Ensure 4-digit ID with leading zeros
  
  # READ IN CSV
  df <- tryCatch(
    {
      read.csv(file)
    },
    error = function(e) {
      failed_ids <<- c(failed_ids, id)
      warning(paste("Failed to read file", file, "for ID", id, ":", e$message))
      return(NULL)
    }
  )
  
  # SKIP IF FILE READ FAILED
  if (is.null(df)) {
    next
  }
  
  # ENSURE FIRST COLUMN IS Period
  if (colnames(df)[1] != "Period") {
    failed_ids_PHS <<- c(failed_ids_PHS, id)
    warning(paste("First column in", file, "is not 'Period'. Skipping file for ID", id))
    next
  }
  
  # CONVERT PERIOD TO HZ
  hz_values <- round(1 / df$Period, 6)
  
  # Extract raw phase matrix (not abs)
  phase_mat <- as.matrix(df[, -1])
  
  # Compute circular mean for each row
  circ_mean <- apply(phase_mat, 1, function(row) {
    row_circ <- circular(row, type = "angles", units = "radians", modulo = "asis")
    mean(row_circ, na.rm = TRUE)
  })
  
  # Compute variance for each row
  row_var <- apply(phase_mat, 1, var, na.rm = TRUE)
  
  # Zero crossings (crossing 0)
  zero_crossings <- apply(phase_mat, 1, function(row) {
    signs <- sign(row)
    sum(diff(signs) != 0, na.rm = TRUE)
  })
  
  # Crossings over circular mean
  crossings_over_circ <- sapply(seq_len(nrow(phase_mat)), function(i) {
    row <- phase_mat[i, ]
    mu <- circ_mean[i]
    signs <- sign(row - mu)
    sum(diff(signs) != 0, na.rm = TRUE)
  })
  
  # Now take ABS for existing phase calculations
  abs_phase_mat <- abs(phase_mat)
  row_sums <- rowSums(abs_phase_mat, na.rm = TRUE)
  non_na_counts <- rowSums(!is.na(abs_phase_mat))
  avg_phase <- ifelse(non_na_counts > 0, row_sums / non_na_counts, NA)
  max_possible <- max(non_na_counts, na.rm = TRUE)
  base_weight <- ifelse(max_possible > 0, (non_na_counts / max_possible)^2, 0)
  
  # Build temp df
  temp_df <- data.frame(
    ID = id,
    Hz = hz_values,
    AVG_PHASE = avg_phase,
    Base_Weight = base_weight,
    Zero_Crossings = zero_crossings,
    Phase_Variance = row_var,
    Circ_Mean = circ_mean,
    Crossings_Over_Circ_Mean = crossings_over_circ
  )
  
  # Inverse row count weighting
  row_counts <- table(temp_df$Hz)
  temp_df$Weight <- temp_df$Base_Weight * (1 / row_counts[as.character(temp_df$Hz)])
  
  # Weighted aggregate for AVG_PHASE
  aggregated_df <- aggregate(
    cbind(AVG_PHASE = AVG_PHASE * Weight, Weight = Weight) ~ Hz + ID,
    data = temp_df,
    FUN = sum,
    na.rm = TRUE
  )
  aggregated_df$AVG_PHASE <- ifelse(
    aggregated_df$Weight > 0,
    aggregated_df$AVG_PHASE / aggregated_df$Weight,
    NA
  )
  
  # Unweighted averages for other metrics
  summary_metrics <- aggregate(
    cbind(Zero_Crossings, Phase_Variance, Circ_Mean, Crossings_Over_Circ_Mean) ~ Hz + ID,
    data = temp_df,
    FUN = mean,
    na.rm = TRUE
  )
  
  # Merge weighted and unweighted summaries
  aggregated_df <- merge(aggregated_df, summary_metrics, by = c("ID", "Hz"))
  
  # Select relevant columns
  aggregated_df <- aggregated_df[, c("ID", "Hz", "AVG_PHASE", "Zero_Crossings", 
                                     "Phase_Variance", "Circ_Mean", 
                                     "Crossings_Over_Circ_Mean")]
  
  # Store results
  summary_list_PHS[[id]] <- aggregated_df
  successful_ids_PHS <<- c(successful_ids_PHS, id)
}

# COMBINE INTO A SINGLE DATAFRAME
ANGLE_SUMMARY <- do.call(rbind, summary_list_PHS)

# CONVERT ID TO CHARACTER TO PRESERVE LEADING ZEROS
ANGLE_SUMMARY$ID <- sprintf("%04d", as.numeric(ANGLE_SUMMARY$ID))

# ENSURE Hz AND AVG_PHASE ARE NUMERIC
ANGLE_SUMMARY$Hz <- as.numeric(ANGLE_SUMMARY$Hz)
ANGLE_SUMMARY$AVG_PHASE <- as.numeric(ANGLE_SUMMARY$AVG_PHASE)

# RESET ROW NAMES
rownames(ANGLE_SUMMARY) <- NULL

## ENSURE ALL DYADS PROCESSED (N= 176 - ID = 1035 cause it's shit)
unique(ANGLE_SUMMARY$ID)

################################################################################
#################### SUMMARIZE CROSS WAVELET POWER MEAN DATA 

# SET WORKING DIRECTORY TO CW POWER FILES
work_dir <- "..\\CW_POWER_FILTERED"
csv_files <- list.files(path = work_dir, pattern = "\\.csv$", full.names = TRUE)

# INITIALIZE LISTS
summary_list <- list()
successful_ids <- character()
failed_ids <- character()

# LOOP THROUGH EACH FILE
for (file in csv_files) {
  # EXTRACT ID FROM FILENAME
  id <- sub(".*([0-9]{4}).*\\.csv$", "\\1", basename(file))
  
  # READ DATA
  df <- tryCatch(
    read.csv(file),
    error = function(e) {
      failed_ids <<- c(failed_ids, id)
      warning(paste("Failed to read file", file, "for ID", id, ":", e$message))
      return(NULL)
    }
  )
  
  if (is.null(df) || colnames(df)[1] != "Period") {
    failed_ids <<- c(failed_ids, id)
    next
  }
  
  # CONVERT PERIOD TO HZ
  hz_values <- round(1 / df$Period, 6)
  
  # COMPUTE MEAN CW POWER ACROSS COLUMNS (excluding Period)
  row_means <- rowMeans(df[, -1], na.rm = TRUE)
  
  # BUILD TEMP DATA FRAME
  temp_df <- data.frame(
    ID = id,
    Hz = hz_values,
    AVG_CW_POWER = row_means
  )
  
  # STORE RESULT
  summary_list[[id]] <- temp_df
  successful_ids <<- c(successful_ids, id)
}

# COMBINE ALL INTO SINGLE DATAFRAME
CW_POWER_SUMMARY <- do.call(rbind, summary_list)

# FORMAT OUTPUT
CW_POWER_SUMMARY$ID <- sprintf("%04d", as.numeric(CW_POWER_SUMMARY$ID))
CW_POWER_SUMMARY$Hz <- as.numeric(CW_POWER_SUMMARY$Hz)
CW_POWER_SUMMARY$AVG_CW_POWER <- as.numeric(CW_POWER_SUMMARY$AVG_CW_POWER)
rownames(CW_POWER_SUMMARY) <- NULL

# OPTIONAL: REPORT OUT
cat("Processed IDs:\n", paste(successful_ids, collapse = ", "), "\n")
cat("Failed IDs:\n", if (length(failed_ids) > 0) paste(failed_ids, collapse = ", ") else "None", "\n")

# ENSURE ALL DYADS PROCESSED 

unique(CW_POWER_SUMMARY$ID)

################################################################################
#################### REGRESS OUT THE INFLUENCE OF FREQUENCY 

# Some of these etimates are highly influences by the sampling frequency. We are handling 
# this by removing the effect of increasing frequency and normalizing both the 
# effect and its influence on the variance using a z-score transformation.
# This should result in corrected estimates for each frequency 


####################### NORMALIZE & CORRECT PHASE ANGLE DATA 

## Zero_Crossings_GAMZ is the normalized score for how many times 0 was crossed 
## in the frequency 

## Cir_Mean_Crossings_GAMZ is the normalized score for how many times the circular
## mean was crossed over in that frequency 

## Phase_Variance_GAMZ is the frequency normalized estimate of the variance observed 
## in phase scores for a given frequency. 

# REDUCE TO YOUR RANGE OF INTEREST TO ENSURE Z-SCORES ARE BASED ON THE PROPER RANGE

ANGLE_SUMMARY_RED <- subset(ANGLE_SUMMARY, Hz >= .04 & Hz <= .50)

### Normalized and Detrended Zero-Crossings 

# 1. Fit GAM to model the nonlinear mean
gam_zc <- gam(Zero_Crossings ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 2. Get residuals from this model
ANGLE_SUMMARY_RED$ZC_GAM_Resid <- resid(gam_zc)

# 3. Fit another GAM to model SD of residuals across Hz
# First, compute absolute residuals (proxy for SD)
ANGLE_SUMMARY_RED$ZC_Abs_Resid <- abs(ANGLE_SUMMARY_RED$ZC_GAM_Resid)
gam_sd <- gam(ZC_Abs_Resid ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 4. Predict mean and SD at each Hz
ANGLE_SUMMARY_RED$ZC_GAM_Pred <- predict(gam_zc, newdata = ANGLE_SUMMARY_RED)
ANGLE_SUMMARY_RED$ZC_GAM_SD <- predict(gam_sd, newdata = ANGLE_SUMMARY_RED)

# 5. Compute final standardized residual
ANGLE_SUMMARY_RED$Zero_Crossings_GAMZ <- (ANGLE_SUMMARY_RED$Zero_Crossings - ANGLE_SUMMARY_RED$ZC_GAM_Pred) / ANGLE_SUMMARY_RED$ZC_GAM_SD

### Normalized and Detrended Phase_Variance  

# 1. Fit GAM to model the nonlinear mean
gam_zv <- gam(Phase_Variance  ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 2. Get residuals from this model
ANGLE_SUMMARY_RED$ZV_GAM_Resid <- resid(gam_zv)

# 3. Fit another GAM to model SD of residuals across Hz
# First, compute absolute residuals (proxy for SD)
ANGLE_SUMMARY_RED$ZV_Abs_Resid <- abs(ANGLE_SUMMARY_RED$ZV_GAM_Resid)
gam_Vsd <- gam(ZV_Abs_Resid ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 4. Predict mean and SD at each Hz
ANGLE_SUMMARY_RED$ZV_GAM_Pred <- predict(gam_zv, newdata = ANGLE_SUMMARY_RED)
ANGLE_SUMMARY_RED$ZV_GAM_SD <- predict(gam_Vsd, newdata = ANGLE_SUMMARY_RED)

# 5. Compute final standardized residual
ANGLE_SUMMARY_RED$Phase_Variance_GAMZ <- (ANGLE_SUMMARY_RED$Phase_Variance - ANGLE_SUMMARY_RED$ZV_GAM_Pred) / ANGLE_SUMMARY_RED$ZV_GAM_SD


### Normalized and Detrended Circular Means Zero-Crossings 

# 1. Fit GAM to model the nonlinear mean
gam_cm <- gam(Crossings_Over_Circ_Mean ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 2. Get residuals from this model
ANGLE_SUMMARY_RED$CM_GAM_Resid <- resid(gam_cm)

# 3. Fit another GAM to model SD of residuals across Hz
# First, compute absolute residuals (proxy for SD)
ANGLE_SUMMARY_RED$CM_Abs_Resid <- abs(ANGLE_SUMMARY_RED$CM_GAM_Resid)
CM_gam_sd <- gam(CM_Abs_Resid ~ s(Hz), data = ANGLE_SUMMARY_RED)

# 4. Predict mean and SD at each Hz
ANGLE_SUMMARY_RED$CM_GAM_Pred <- predict(gam_cm, newdata = ANGLE_SUMMARY_RED)
ANGLE_SUMMARY_RED$CM_GAM_SD <- predict(CM_gam_sd, newdata = ANGLE_SUMMARY_RED)

# 5. Compute final standardized residual
ANGLE_SUMMARY_RED$Cir_Mean_Crossings_GAMZ <- (ANGLE_SUMMARY_RED$Crossings_Over_Circ_Mean - ANGLE_SUMMARY_RED$CM_GAM_Pred) / ANGLE_SUMMARY_RED$CM_GAM_SD


## DROP WORKING VARIABLES 
names(ANGLE_SUMMARY_RED)


ANGLE_SUMMARY_RED <- ANGLE_SUMMARY_RED %>%
  select(c("ID", "Hz", "AVG_PHASE", "Zero_Crossings", "Phase_Variance", "Circ_Mean",
           "Zero_Crossings_GAMZ", "Phase_Variance_GAMZ", "Cir_Mean_Crossings_GAMZ"))


##### VISUALIZE PHASE ANGLE CORRECTIONS  


ANGLE_SUMMARY_RED$SAMPLE <- ifelse(ANGLE_SUMMARY_RED$ID > 1000, "PDM", "DORRY")

ANGLE_SUMMARY_RED <- ANGLE_SUMMARY_RED[ANGLE_SUMMARY_RED$ID != "1035", ]

PLOT_DATA <- subset(ANGLE_SUMMARY_RED, Hz >= .04 & Hz <= .50)

# CREATE THE FACETED PLOT
p <- ggplot(PLOT_DATA, aes(x = Hz, y = Cir_Mean_Crossings_GAMZ, color = ID, group = ID)) +
  geom_line() +               # CONNECT POINTS TO FORM DENSITY-LIKE LINES
  geom_point(size = 0.8) +    # ADD POINTS FOR CLARITY
  theme_minimal() +           # CLEAN THEME
  labs(
    title = "Distribution of AVG_PHASE by Hz for Each Participant",
    x = "Hz",
    y = "AVG_PHASE"
  ) +
  facet_wrap(~SAMPLE, scales = "free_y") +  # FACET BY SAMPLE
  theme(
    legend.position = "none",  # REMOVE LEGEND
    plot.title = element_text(hjust = 0.5)
  )

# DISPLAY PLOT
print(p)


####################### NORMALIZE & CORRECT THE POWER MEAN DATA 

# REDUCE TO YOUR RANGE OF INTEREST TO ENSURE Z-SCORES ARE BASED ON THE PROPER RANGE

CW_POWER_SUMMARY_RED <- subset(CW_POWER_SUMMARY, Hz >= .04 & Hz <= .50)

# CW_POWER_GAMZ is the frequency corrected power mean score

# Fit a GAM to estimate expected mean CW power at each frequency
gam_power <- gam(AVG_CW_POWER ~ s(Hz), data = CW_POWER_SUMMARY_RED)

# Get residuals = deviation from expected power at that frequency
CW_POWER_SUMMARY_RED$CW_POWER_Resid <- resid(gam_power)

# Estimate expected variability (SD) as a function of frequency
CW_POWER_SUMMARY_RED$Abs_Resid <- abs(CW_POWER_SUMMARY_RED$CW_POWER_Resid)
gam_sd <- gam(Abs_Resid ~ s(Hz), data = CW_POWER_SUMMARY_RED)

# Predict expected mean and SD
CW_POWER_SUMMARY_RED$Expected_Mean <- predict(gam_power)
CW_POWER_SUMMARY_RED$Expected_SD <- predict(gam_sd)

# Compute Z-scored, detrended CW power
CW_POWER_SUMMARY_RED$CW_POWER_GAMZ <- (CW_POWER_SUMMARY_RED$AVG_CW_POWER - CW_POWER_SUMMARY_RED$Expected_Mean) / CW_POWER_SUMMARY_RED$Expected_SD

## DROP WORKING VARIABLES 
names(CW_POWER_SUMMARY_RED)


CW_POWER_SUMMARY_RED <- CW_POWER_SUMMARY_RED %>%
  select(c("ID", "Hz", "AVG_CW_POWER", "CW_POWER_GAMZ"))

##### VISUALIZE & COMPARE POWER MEAN CORRECTIONS  


CW_POWER_SUMMARY_RED$SAMPLE <- ifelse(CW_POWER_SUMMARY_RED$ID > 1000, "PDM", "DORRY")

CW_POWER_SUMMARY_RED <- CW_POWER_SUMMARY_RED[CW_POWER_SUMMARY_RED$ID != "1035", ]


# Step 1: Pivot to long format
plot_long <- CW_POWER_SUMMARY_RED %>%
  select(ID, Hz, SAMPLE, CW_POWER_GAMZ, AVG_CW_POWER) %>%
  pivot_longer(cols = c(CW_POWER_GAMZ, AVG_CW_POWER),
               names_to = "Metric", values_to = "Value")

# Step 2: Create a facet grid (rows = Metric, columns = SAMPLE)
p <- ggplot(plot_long, aes(x = Hz, y = Value, color = ID, group = ID)) +
  geom_line() +
  geom_point(size = 0.8) +
  theme_minimal() +
  labs(
    title = "Cross-Wavelet Power by Frequency (Top = GAMZ, Bottom = Raw)",
    x = "Frequency (Hz)",
    y = "Value"
  ) +
  facet_grid(rows = vars(Metric), cols = vars(SAMPLE), scales = "free_y") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10)
  )

print(p)


################################################################################
##################### MERGE COHERENCE AND PHASE DATA 

# Merge on ID and Hz

FULL_DF <- full_join(ANGLE_SUMMARY_RED, CW_POWER_SUMMARY_RED, by = c("ID", "Hz", "SAMPLE"))

hist(FULL_DF$CW_POWER_GAMZ)
hist(FULL_DF$AVG_CW_POWER)


################################################################################
##################### SMOOTH SIGNALS OUT 

# WE ARE GOING TO SMOOTH THE SIGNALS OUT IN AN ATTEMPT TO ELLUCIATE WHERE THE 
# TRUE PEAKS AND VALLIES ARE AND WHICH FLUCTUATIONS CAN BE DISREGARDED. 


# FILTER VARIABLES OF INTEREST

metrics_to_smooth <- c("Zero_Crossings_GAMZ", "Phase_Variance_GAMZ", "Cir_Mean_Crossings_GAMZ" ,"CW_POWER_GAMZ")

## TRIM THE HZ TO THE TIMES OF INTEREST TO ENSURE THAT THE Z-SCORE IS RELATIVE TO 
## APPLICABLE WINDOWS

FULL_DF_RED <- subset(FULL_DF, Hz >= .04 & Hz <= .50)

# Pivot to long format
df_long <- FULL_DF_RED %>%
  select(ID, Hz,
         Zero_Crossings_GAMZ,
         Phase_Variance_GAMZ,
         Cir_Mean_Crossings_GAMZ,
         CW_POWER_GAMZ) %>%
  pivot_longer(
    cols = c(Zero_Crossings_GAMZ, Phase_Variance_GAMZ, Cir_Mean_Crossings_GAMZ, CW_POWER_GAMZ),
    names_to = "Metric",
    values_to = "Value"
  )


df_smooth <- df_long %>%
  filter(Metric %in% metrics_to_smooth) %>%
  group_by(ID, Metric) %>%
  arrange(Hz) %>%
  mutate(Smoothed_Value = stats::loess(Value ~ Hz, span = 0.17)$fitted) %>%
  ungroup()


################################################################################
##################### VISUALIZE AN EXAMPLES OF SMOOTHENING 

# Choose your ID and Metric to plot
selected_id <- "0101"
selected_metric <- "CW_POWER_GAMZ"

# Filter data for that ID and Metric
plot_data <- df_smooth %>%
  filter(ID == selected_id, Metric == selected_metric)

plot_data <- subset(plot_data, Hz >= .04 & Hz <= .50)

# Plot original and smoothed values over Hz
ggplot(plot_data, aes(x = Hz)) +
  geom_line(aes(y = Value), color = "blue", alpha = 0.5, size = 1) +  # raw data points
  geom_line(aes(y = Smoothed_Value), color = "red", size = 2) +       # smoothed line
  theme_minimal() +
  labs(
    title = paste("Original (blue) vs Smoothed (red) for", selected_metric, "ID:", selected_id),
    x = "Frequency (Hz)",
    y = "Value (Z-Score)"
  )

################################################################################
##################### IDENTIFY OPTIMIZED COHERENCE ISALNDS MULTI-MODALLY

##### VISUALIZE AREAS OF OVERLAP BETWEEN CURVES

# Define metrics of interest
metrics_to_plot <- c("Zero_Crossings_GAMZ", "Phase_Variance_GAMZ", "CW_POWER_GAMZ")

# Create output directory
output_dir <- "../OPTIMIZED_METRICS_PLOTS"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each unique ID
for (target_id in unique(df_smooth$ID)) {
  
  # Filter for this ID and metrics, and restrict Hz range
  plot_data <- df_smooth %>%
    filter(ID == target_id, Metric %in% metrics_to_plot) %>%
    filter(Hz > 0.02 & Hz < 0.50)
  
  # Skip if no data left (just in case)
  if (nrow(plot_data) == 0) next
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Hz, y = Smoothed_Value, color = Metric)) +
    geom_line(size = 1.2) +
    labs(
      title = paste("Smoothed Metrics for ID", target_id),
      x = "Frequency (Hz)",
      y = "Z-scored Smoothed Value",
      color = "Metric"
    ) +
    theme_minimal(base_size = 14)
  
  # Save plot
  output_path <- file.path(output_dir, paste0(target_id, "_OPTIMIZED_METRICS_PLOT.pdf"))
  ggsave(output_path, p, width = 8, height = 5)
  
  message("Saved plot for ID: ", target_id)
}

################################################################################
##################### SAVE THIS DATA PRIOR TO LOCATING THE OPTIMAL BAND 

##### Pivot back to wide format 

# Create separate data frames for Value and Smoothed_Value

df_values <- df_smooth %>%
  select(ID, Hz, Metric, Value) %>%
  pivot_wider(names_from = Metric, values_from = Value, names_prefix = "Value_")

df_smoothed <- df_smooth %>%
  select(ID, Hz, Metric, Smoothed_Value) %>%
  pivot_wider(names_from = Metric, values_from = Smoothed_Value, names_prefix = "Smoothed_")

# Join the two wide data frames on ID and Hz
df_wide <- df_values %>%
  full_join(df_smoothed, by = c("ID", "Hz"))

# Reorder by ID
df_wide_ordered <- df_wide %>%
  arrange(ID, Hz)

########## Add back in the unchanged Summary values 

OGs <- FULL_DF_RED %>% 
  select(c("ID", "Hz", "AVG_PHASE", "Zero_Crossings", "Phase_Variance", 
           "AVG_CW_POWER", "SAMPLE"))

FULL_SUMMARY <- df_wide_ordered %>%
  full_join(OGs, by = c("ID","Hz"))

# Rename to my desired conventions

FULL_SUMMARY_NAMED <- FULL_SUMMARY %>%
  rename(PHS_ZERO_CROSS_GAMZ = Value_Zero_Crossings_GAMZ,
         PHS_VAR_GAMZ = Value_Phase_Variance_GAMZ, 
         CW_POWER_GAMZ = Value_CW_POWER_GAMZ,
         PHS_ZERO_CROSS_GAMZ_SM = Smoothed_Zero_Crossings_GAMZ, 
         PHS_VAR_GAMZ_SM = Smoothed_Phase_Variance_GAMZ, 
         CW_POWER_GAMZ_SM = Smoothed_CW_POWER_GAMZ, 
         PHS_ZERO_CROSS = Zero_Crossings, 
         PHS_VAR = Phase_Variance,
         CW_POWER = AVG_CW_POWER)


# Reorder and drop the things I don't need
FULL_SUMMARY_NAMED_RED <- FULL_SUMMARY_NAMED %>%
  select(c("ID", "SAMPLE", "Hz", 
           "CW_POWER", "CW_POWER_GAMZ", "CW_POWER_GAMZ_SM",
           "PHS_VAR", "PHS_VAR_GAMZ", "PHS_VAR_GAMZ_SM",
           "PHS_ZERO_CROSS", "PHS_ZERO_CROSS_GAMZ", "PHS_ZERO_CROSS_GAMZ_SM"))


# SAVE IT! 
write.csv(FULL_SUMMARY_NAMED_RED, "..\\CW_POWER_PHASE_FRQ_SUMMARY.csv", row.names = FALSE)


