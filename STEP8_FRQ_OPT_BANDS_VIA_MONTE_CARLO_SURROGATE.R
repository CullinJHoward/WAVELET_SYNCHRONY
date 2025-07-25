################################################################################
######################### LIBRARY SHTUFFF

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)


################################################################################
######################### SET WORKING DIRECTORY

POUNDTOWN <- 0

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
}

setwd(work_dir)


################################################################################
######################### LOAD SUMMARY DF & SPLIT 

df <- read.csv("CW_POWER_PHASE_FRQ_SUMMARY.csv")
df$ID <- sprintf("%04d", as.numeric(df$ID))  # preserve leading zeros
unique_ids <- unique(df$ID)
dyad_list <- split(df, df$ID)

# SPLIT INTO THEIR OWN DFs
for (id in unique_ids) {
  df_subset <- df %>% filter(ID == id)
  assign(paste0("D_", id), df_subset)
}


################################################################################
######################### QUANTIFY THRESHOLD EACH DOMAIN OF SYCHRONY 


## SET UP FUNCTIONS

# STEP 1: Generate sample-level surrogate with smoothing

generate_sample_surrogate <- function(df_sample, metric) {
  surrogate_sample <- df_sample %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    mutate(Surrogate = sample(.[[metric]], replace = TRUE)) %>%
    arrange(Hz)
  
  # Ensure Hz values are numeric and matched
  surrogate_sample$Hz <- round(surrogate_sample$Hz, 4)
  
  # Apply smoothing (rolling average)
  surrogate_sample <- surrogate_sample %>%
    mutate(Smoothed = zoo::rollapply(Surrogate, width = smooth_window, fill = NA, align = "center", FUN = mean)) %>%
    select(Hz, Smoothed) %>%
    rename(Surrogate = Smoothed)
  
  return(surrogate_sample)
}

# STEP 2: Run surrogates and calculate p-values 

run_surrogates_for_dyad <- function(df_dyad, df_sample, metric, n_surrogates) {
  df_dyad <- df_dyad %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    arrange(Hz)
  
  if (nrow(df_dyad) == 0) {
    message("Skipping ID ", unique(df_dyad$ID), ": no valid data after filtering.")
    return(NULL)
  }
  
  Hz_vals <- round(df_dyad$Hz, 4)
  observed_vals <- df_dyad[[metric]]
  
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df_sample, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # One-tailed test: proportion of surrogate values >= observed
  p_vals <- sapply(1:length(Hz_vals), function(i) {
    mean(surrogate_matrix[i, ] >= observed_vals[i], na.rm = TRUE)
  })
  
  result <- tibble(
    ID = unique(df_dyad$ID),
    Hz = Hz_vals,
    Observed = observed_vals,
    P_value = p_vals
  )
  
  return(result)
}

# STEP 3: Identify contiguous significant bands 

find_all_significant_bands <- function(p_result, alpha = 0.05) {
  p_result <- p_result %>%
    arrange(Hz) %>%
    mutate(Significant = P_value < alpha)
  
  rle_result <- rle(p_result$Significant)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, sig = rle_result$values) %>%
    filter(sig)
  
  if (nrow(run_df) == 0) return(NULL)
  
  bands <- lapply(1:nrow(run_df), function(i) {
    rows <- p_result[run_df$start[i]:run_df$end[i], ]
    tibble(
      ID = rows$ID[1],
      Band_Lower = min(rows$Hz),
      Band_Upper = max(rows$Hz),
      N_Frequencies = nrow(rows),
      Min_P = min(rows$P_value, na.rm = TRUE),
      Peak_Hz = rows$Hz[which.min(rows$P_value)]
    )
  })
  
  return(bind_rows(bands))
}


########################################################
### RUN FUNCTIONS FOR EACH SYNCHRONY METRIC 

### USER INPUT PARAMETERS - Cross-wavelet Power 

metric <- "CW_POWER_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

## RUN 

if (!dir.exists("CW_POWER_SURROGATE_BANDS")) {
  dir.create("CW_POWER_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("CW_POWER Surrogate Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("CW_POWER_SURROGATE_BANDS/", id, "_CW_POWER_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}


### SUMMARY AND SAVE  

SURROGATE_CW_POWER_SUMMARY <- bind_rows(summary_table)

write.csv(SURROGATE_CW_POWER_SUMMARY, "SURROGATE_CW_POWER_SUMMARY.csv", row.names = FALSE)


### USER INPUT PARAMETERS - Phase Zero-Cross  

metric <- "PHS_ZERO_CROSS_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

## RUN 

if (!dir.exists("PHS_ZERO_CROSS_SURROGATE_BANDS")) {
  dir.create("PHS_ZERO_CROSS_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("Phase Zero Cross Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("PHS_ZERO_CROSS_SURROGATE_BANDS/", id, "_PHS_ZERO_CROSS_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}

### SUMMARY AND SAVE

SURROGATE_PHS_ZERO_CROSS_SUMMARY <- bind_rows(summary_table)
write.csv(SURROGATE_PHS_ZERO_CROSS_SUMMARY, "SURROGATE_PHS_ZERO_CROSS_SUMMARY.csv", row.names = FALSE)


### USER INPUT PARAMETERS - PHASE VARIANCE 

metric <- "PHS_VAR_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

## RUN 

if (!dir.exists("PHS_VARIANCE_SURROGATE_BANDS")) {
  dir.create("PHS_VARIANCE_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("Phase Variance Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("PHS_VARIANCE_SURROGATE_BANDS/", id, "_PHS_VARIANCE_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}

### SUMMARY AND SAVE  

SURROGATE_PHS_VARIANCE_SUMMARY <- bind_rows(summary_table)
write.csv(SURROGATE_PHS_VARIANCE_SUMMARY, "SURROGATE_PHS_VARIANCE_SUMMARY.csv", row.names = FALSE)






################################################################################
################################# ARCHIVE ######################################
################################################################################

###########################################################
############ IDENTIFY FREQUENCY OPTIMIZED SYNCHRONY BANDS 
df <- read.csv("CW_POWER_PHASE_FRQ_SUMMARY.csv")
df$ID <- sprintf("%04d", as.numeric(df$ID))  # preserve leading zeros
SURROGATE_CW_POWER_SUMMARY <- read.csv("SURROGATE_CW_POWER_SUMMARY.csv")
SURROGATE_CW_POWER_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_CW_POWER_SUMMARY$ID))  # preserve leading zeros
SURROGATE_PHS_ZERO_CROSS_SUMMARY <- read.csv("SURROGATE_PHS_ZERO_CROSS_SUMMARY.csv")
SURROGATE_PHS_ZERO_CROSS_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_PHS_ZERO_CROSS_SUMMARY$ID))  # preserve leading zeros
SURROGATE_PHS_VARIANCE_SUMMARY <- read.csv("SURROGATE_PHS_VARIANCE_SUMMARY.csv")
SURROGATE_PHS_VARIANCE_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_PHS_VARIANCE_SUMMARY$ID))  # preserve leading zeros

## ADD BACK INTO THE SAMPLE DF 

# Function to mark Hzs inside bands for a given domain
flag_significant_hz <- function(df, summary_df, var_name) {
  df[[var_name]] <- 0  # Default to 0
  
  for (i in seq_len(nrow(summary_df))) {
    id <- summary_df$ID[i]
    lower <- summary_df$Band_Lower[i]
    upper <- summary_df$Band_Upper[i]
    
    df[[var_name]] <- ifelse(
      df$ID == id & df$Hz >= lower & df$Hz <= upper,
      1,
      df[[var_name]]
    )
  }
  
  return(df)
}

# Start with your original long-format dataset
df_annotated <- df

# Apply flagging functions for each surrogate summary
df_annotated <- df_annotated %>%
  flag_significant_hz(SURROGATE_CW_POWER_SUMMARY, "CW_Shz") %>%
  flag_significant_hz(SURROGATE_PHS_ZERO_CROSS_SUMMARY, "PH_ZC_Shz") %>%
  flag_significant_hz(SURROGATE_PHS_VARIANCE_SUMMARY, "PH_VAR_Shz")


############################
## CRITERIA #1

# criteria 1 is the most stringent. It identifies the frequency optimized band as
# overlap between surrogate-defined significant contiguous frequencies. 
# This requires that the CW_POWER band MUST be significant AND overlap with a 
# significant Phase Zero_cross OR Phase variance band (or both!)


## CREATE A VARIABLE THAT NOTES CRITERIA #1 OVERLAP AT A GIVEN Hz

df_annotated$CRIT1_FOS <- ifelse(
  df_annotated$CW_Shz == 1 & 
    (df_annotated$PH_ZC_Shz == 1 | df_annotated$PH_VAR_Shz == 1), 
  1, 0
)

## PLOT EACH DYADS DATA OVERLAYING OVERLAPPING SIGNIFICANT SURROGATE REGIONS 

# Create output folder if not present
if (!dir.exists("FRQ_OPT_CRT1_PLOTS")) {
  dir.create("FRQ_OPT_CRT1_PLOTS")
}


# Initialize summary list
summary_list <- list()

# Helper Function: Overlap Bands (for shading any region with 2+ overlapping signals)
get_overlap_bands <- function(df_dyad) {
  rle_obj <- rle(df_dyad$CRIT1_FOS == 1)
  ends <- cumsum(rle_obj$lengths)
  starts <- c(1, head(ends + 1, -1))
  run_df <- tibble(start = starts, end = ends, value = rle_obj$values, len = rle_obj$lengths) %>%
    filter(value, len >= 3)
  
  if (nrow(run_df) == 0) return(NULL)
  
  bands <- lapply(1:nrow(run_df), function(i) {
    rows <- df_dyad[run_df$start[i]:run_df$end[i], ]
    tibble(xmin = min(rows$Hz), xmax = max(rows$Hz))
  })
  
  bind_rows(bands)
}

# Helper Function: Optimal Band Selection (1 per dyad)
get_optimal_band <- function(df_dyad) {
  rle_obj <- rle(df_dyad$CRIT1_FOS == 1)
  ends <- cumsum(rle_obj$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, value = rle_obj$values, len = rle_obj$lengths) %>%
    filter(value, len >= 3)
  
  if (nrow(run_df) == 0) return(NULL)
  
  band_info <- list()
  
  for (i in seq_len(nrow(run_df))) {
    rows <- df_dyad[run_df$start[i]:run_df$end[i], ]
    
    has_all_three <- any(rows$CW_Shz == 1 & rows$PH_VAR_Shz == 1 & rows$PH_ZC_Shz == 1)
    
    sig_matrix <- rows %>% 
      select(CW_Shz, PH_VAR_Shz, PH_ZC_Shz) %>% 
      summarise_all(~mean(. == 1))
    
    nonsig_domain <- names(sig_matrix)[which.min(sig_matrix)]
    
    nonsig_mean <- mean(rows[[paste0(switch(nonsig_domain,
                                            CW_Shz = "CW_POWER_GAMZ_SM",
                                            PH_VAR_Shz = "PHS_VAR_GAMZ_SM",
                                            PH_ZC_Shz = "PHS_ZERO_CROSS_GAMZ_SM"))]], na.rm = TRUE)
    
    band_info[[i]] <- list(
      xmin = min(rows$Hz),
      xmax = max(rows$Hz),
      has_all_three = has_all_three,
      nonsig_value = nonsig_mean
    )
  }
  
  band_df <- bind_rows(band_info)
  
  if (any(band_df$has_all_three)) {
    best_band <- band_df %>% filter(has_all_three) %>% slice(1)
  } else {
    best_band <- band_df %>% arrange(desc(nonsig_value)) %>% slice(1)
  }
  
  return(tibble(xmin = best_band$xmin, xmax = best_band$xmax))
}


# Loop through each dyad
for (id in unique(df_annotated$ID)) {
  df_dyad <- df_annotated %>% filter(ID == id)
  
  # Shaded bands: all valid 3+ contiguous overlapping points (CRIT1_FOS == 1)
  shade_bands <- get_overlap_bands(df_dyad)
  
  # Optimal band: 1 per dyad based on rules
  optimal_band <- get_optimal_band(df_dyad)
  
  # Red dots for significant values
  sig_cw  <- df_dyad %>% filter(CW_Shz == 1)
  sig_var <- df_dyad %>% filter(PH_VAR_Shz == 1)
  sig_zc  <- df_dyad %>% filter(PH_ZC_Shz == 1)
  
  # Get star placement for optimal band
  star_point <- if (!is.null(optimal_band)) {
    center_hz <- mean(c(optimal_band$xmin, optimal_band$xmax))
    y_max <- max(df_dyad$CW_POWER_GAMZ_SM, df_dyad$PHS_VAR_GAMZ_SM, df_dyad$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
    tibble(Hz = center_hz, y = y_max + 0.2)
  } else {
    NULL
  }
  
  # Build plot
  p <- ggplot(df_dyad, aes(x = Hz)) +
    # Shaded overlap bands
    {if (!is.null(shade_bands)) geom_rect(
      data = shade_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey80", alpha = 0.4
    )} +
    
    # Plot curves
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER")) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VARIANCE")) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "PHASE_ZERO_CROSS")) +
    
    # Red dots for significance
    geom_point(data = sig_cw, aes(x = Hz, y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_var, aes(x = Hz, y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_zc, aes(x = Hz, y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    
    # Star at center of optimal band
    {if (!is.null(star_point)) geom_point(
      data = star_point,
      aes(x = Hz, y = y),
      shape = 8, size = 4, color = "black"
    )} +
    
    scale_color_manual(values = c(
      "CW_POWER" = "blue",
      "PHASE_VARIANCE" = "purple",
      "PHASE_ZERO_CROSS" = "green4"
    )) +
    labs(
      title = paste("Synchrony Domains with Overlap — Dyad", id),
      x = "Hz", y = "Value", color = "Metric"
    ) +
    theme_minimal(base_size = 14)
  
  # Save plot
  ggsave(
    filename = paste0("FRQ_OPT_CRT1_PLOTS/", id, "_CRIT1_OVERLAP_PLOT.pdf"),
    plot = p, width = 8, height = 5
  )

  summary_list[[id]] <- tibble(
    ID = id,
    Any_Overlap = if (!is.null(optimal_band)) "Yes" else "No",
    N_Hz = if (!is.null(optimal_band)) {
      sum(df_dyad$Hz >= optimal_band$xmin & df_dyad$Hz <= optimal_band$xmax)
    } else {
      0
    },
    OPT_LB = if (!is.null(optimal_band)) optimal_band$xmin else NA,
    OPT_UB = if (!is.null(optimal_band)) optimal_band$xmax else NA
  )
}

# Combine and save summary
summary_df <- bind_rows(summary_list)
summary_df$OPT_CRIT <- ifelse(summary_df$Any_Overlap == "Yes", "CRIT.1", "N/A")
  
 
table(summary_df$OPT_CRIT) #46% in Criteria #1

# REDUYCE DF FOR OVERALL SUMMARY
CRIT1_IDS <- summary_df %>%
  filter(OPT_CRIT == "CRIT.1") %>%
  pull(ID)

# REDUCE JUST TO CRIT 1 WINNERS
CRIT1_SUMMARY_DF <- summary_df %>%
  filter(ID %in% CRIT1_IDS)

# REMOVE WORKING VARIABLE
CRIT1_SUMMARY_DF_RED <- CRIT1_SUMMARY_DF %>%
  select(-c(Any_Overlap))

############################
## CRITERIA #2

# In the event that there is no significant overlapping phase band, the optimized
# frequency band will be located within just the CW_POWER significant bands

# Get list of IDs that were not fit by CRIT.1
non_crit1_ids <- summary_df %>%
  filter(OPT_CRIT == "N/A") %>%
  pull(ID)

# Subset df to only those participants
df_CRIT2 <- df_annotated %>%
  filter(ID %in% non_crit1_ids)

## REDUCT TO JUST THOSE IDs WITH AT LEAST 1 CW_POWER Hz is sig

# List the IDs
CRTI2_IDs <- df_CRIT2 %>%
  group_by(ID) %>%
  summarise(has_crit2 = any(CW_Shz == 1, na.rm = TRUE)) %>%
  filter(has_crit2) %>%
  pull(ID)

# REDUCE THE DF AGAIN 
df_CRIT2 <- df_CRIT2 %>%
  filter(ID %in% CRTI2_IDs)

#How many we got 
unique(df_CRIT2$ID) # 63 dyads, 36%

#### ASSESS FIT TO CRITERIA 2 FOR EACH DYAD

# Create output folder if it doesn't exist
if (!dir.exists("FRQ_OPT_CRT2_PLOTS")) {
  dir.create("FRQ_OPT_CRT2_PLOTS")
}

# INITIALIZE LIST 
summary_list_C2 <- list()

# Loop through each ID
unique_ids <- unique(df_CRIT2$ID)

for (id in unique_ids) {
  
  df_sub <- df_CRIT2 %>% filter(ID == id)
  
  # Define CRIT2 region: CW significant + either phase metric >= 0
  df_sub <- df_sub %>%
    mutate(
      CRIT2_REGION = ifelse(
        CW_Shz == 1 & (PHS_VAR_GAMZ_SM >= 0 | PHS_ZERO_CROSS_GAMZ_SM >= 0),
        1, 0
      ),
      REGION_ID = data.table::rleid(CRIT2_REGION) * CRIT2_REGION
    )
  
  # Identify contiguous regions of 3+ Hz
  band_regions <- df_sub %>%
    filter(CRIT2_REGION == 1) %>%
    group_by(REGION_ID) %>%
    filter(n() >= 3) %>%
    summarise(
      xmin = min(Hz),
      xmax = max(Hz),
      n_Hz = n(),
      .groups = "drop"
    )
  
  # Identify optimal band: longest (most Hz)
  if (nrow(band_regions) > 0) {
    optimal_band <- band_regions %>%
      slice_max(n_Hz, with_ties = FALSE)
    
    star_point <- data.frame(
      Hz = mean(c(optimal_band$xmin, optimal_band$xmax)),
      y = max(
        df_sub$CW_POWER_GAMZ_SM,
        df_sub$PHS_VAR_GAMZ_SM,
        df_sub$PHS_ZERO_CROSS_GAMZ_SM,
        na.rm = TRUE
      ) + 0.3
    )
  } else {
    optimal_band <- NULL
    star_point <- NULL
  }
  
  # ---- Plot
  p <- ggplot(df_sub, aes(x = Hz)) +
    # All valid shaded bands
    {if (nrow(band_regions) > 0) {
      geom_rect(data = band_regions,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE,
                fill = "plum1", alpha = 0.3)
    }} +
    # Curves
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER"), size = 1) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VAR"), size = 1) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "ZERO_CROSS"), size = 1) +
    # Significant red dots
    geom_point(data = df_sub %>% filter(CW_Shz == 1),
               aes(y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_VAR_Shz == 1),
               aes(y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_ZC_Shz == 1),
               aes(y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    # Add star for optimal band
    {if (!is.null(star_point)) {
      geom_point(data = star_point, aes(x = Hz, y = y), shape = 8, size = 3, color = "purple")
    }} +
    scale_color_manual(
      values = c("CW_POWER" = "black", "PHASE_VAR" = "blue", "ZERO_CROSS" = "darkgreen"),
      name = "Signal"
    ) +
    ggtitle(paste0("ID: ", id)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("GAMZ (Z-Score)") +
    xlab("Hz")
  
  # Save plot
  ggsave(paste0("FRQ_OPT_CRT2_PLOTS/", id, "_CRIT1_OVERLAP_PLOT.pdf"), p, width = 8, height = 5)
  
  # --- Add to summary
  summary_list_C2[[id]] <- data.frame(
    ID = id,
    Found_Crit2_Band = if (!is.null(optimal_band)) "Yes" else "No",
    OPT_LB = if (!is.null(optimal_band)) optimal_band$xmin else NA,
    OPT_UB = if (!is.null(optimal_band)) optimal_band$xmax else NA,
    N_Hz = if (!is.null(optimal_band)) optimal_band$n_Hz else NA
  )
}

# Combine summary
CRIT2_SUMMARY_DF <- do.call(rbind, summary_list_C2)
table(CRIT2_SUMMARY_DF$Found_Crit2_Band) # 36 fit this criteria, 21%

CRIT2_SUMMARY_DF$OPT_CRIT <- ifelse(CRIT2_SUMMARY_DF$Found_Crit2_Band == "Yes", "CRIT.2", "N/A")

# REDUCE DF FOR OVERALL SUMMARY
CRIT2_IDS <- CRIT2_SUMMARY_DF %>%
  filter(OPT_CRIT == "CRIT.2") %>%
  pull(ID)

# REDUCE JUST TO CRIT 1 WINNERS
CRIT2_SUMMARY_DF <- CRIT2_SUMMARY_DF %>%
  filter(ID %in% CRIT2_IDS)

# REMOVE WORKING VARIABLE
CRIT2_SUMMARY_DF_RED <- CRIT2_SUMMARY_DF %>%
  select(-c(Found_Crit2_Band))

## ADD BACK INTO THE OVERALL SUMMARY DF 

C1C2_OPT_SUM <- rbind(CRIT1_SUMMARY_DF_RED,CRIT2_SUMMARY_DF_RED)

############################
## CRITERIA #3

# In the event that there is no significant overlapping phase band, and 
# there either is no significant CW power band or no sig CW power band with at least
# average levels of phase, then we will examine their significant phase bands 
# and label the regions that of overlap with at least mean levels of CW_POWER 

# Get list of IDs that were not fit by CRIT.1 or CRIT.2
C1C2_OPT_IDs <- unique(C1C2_OPT_SUM$ID) 

# Subset df to only those participants
df_CRIT3 <- df_annotated %>%
  filter(!ID %in% C1C2_OPT_IDs)


#### ASSESS FIT TO CRITERIA 3 FOR EACH DYAD

# Create output folder if it doesn't exist
if (!dir.exists("FRQ_OPT_CRT3_PLOTS")) {
  dir.create("FRQ_OPT_CRT3_PLOTS")
}

# Initialize list for summary data
summary_list_C3 <- list()

# Loop through each ID
unique_ids_C3 <- unique(df_CRIT3$ID)

for (id in unique_ids_C3) {
  
  df_sub <- df_CRIT3 %>% filter(ID == id)
  
  # --- Define CRIT3 region:
  df_sub <- df_sub %>%
    mutate(
      CRIT3_REGION = ifelse(
        (PH_ZC_Shz == 1 | PH_VAR_Shz == 1) & CW_POWER_GAMZ_SM >= 0,
        1, 0
      ),
      REGION_ID = data.table::rleid(CRIT3_REGION) * CRIT3_REGION
    )
  
  # --- Identify 3+ contiguous Hz bands
  band_regions <- df_sub %>%
    filter(CRIT3_REGION == 1) %>%
    group_by(REGION_ID) %>%
    filter(n() >= 3) %>%
    summarise(
      xmin = min(Hz),
      xmax = max(Hz),
      n_Hz = n(),
      .groups = "drop"
    )
  
  # --- Identify optimal band (longest one)
  if (nrow(band_regions) > 0) {
    optimal_band <- band_regions %>%
      slice_max(n_Hz, with_ties = FALSE)
    
    star_point <- data.frame(
      Hz = mean(c(optimal_band$xmin, optimal_band$xmax)),
      y = max(
        df_sub$CW_POWER_GAMZ_SM,
        df_sub$PHS_VAR_GAMZ_SM,
        df_sub$PHS_ZERO_CROSS_GAMZ_SM,
        na.rm = TRUE
      ) + 0.3
    )
  } else {
    optimal_band <- NULL
    star_point <- NULL
  }
  
  # --- Plot
  p <- ggplot(df_sub, aes(x = Hz)) +
    # Shaded purple bands
    {if (nrow(band_regions) > 0) {
      geom_rect(data = band_regions,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE,
                fill = "#88C97E", alpha = 0.3)
    }} +
    # Curves
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER"), size = 1) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VAR"), size = 1) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "ZERO_CROSS"), size = 1) +
    # Red dots for significance
    geom_point(data = df_sub %>% filter(CW_Shz == 1),
               aes(y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_VAR_Shz == 1),
               aes(y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_ZC_Shz == 1),
               aes(y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    # Star for optimal band
    {if (!is.null(star_point)) {
      geom_point(data = star_point, aes(x = Hz, y = y), shape = 8, size = 3, color = "darkgreen")
    }} +
    scale_color_manual(
      values = c("CW_POWER" = "black", "PHASE_VAR" = "blue", "ZERO_CROSS" = "darkgreen"),
      name = "Signal"
    ) +
    ggtitle(paste0("ID: ", id)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("GAMZ (Z-Score)") +
    xlab("Hz")
  
  # --- Save plot
  ggsave(paste0("FRQ_OPT_CRT3_PLOTS/", id, "_CRIT3_PLOT.pdf"), p, width = 8, height = 5)
  
  # --- Add to summary
  summary_list_C3[[id]] <- data.frame(
    ID = id,
    Found_Crit3_Band = if (!is.null(optimal_band)) "Yes" else "No",
    OPT_LB = if (!is.null(optimal_band)) optimal_band$xmin else NA,
    OPT_UB = if (!is.null(optimal_band)) optimal_band$xmax else NA,
    N_Hz = if (!is.null(optimal_band)) optimal_band$n_Hz else NA
  )
}

# --- Final summary
CRIT3_SUMMARY_DF <- do.call(rbind, summary_list_C3)
CRIT3_SUMMARY_DF$OPT_CRIT <- ifelse(CRIT3_SUMMARY_DF$Found_Crit3_Band == "Yes", "CRIT.3", "N/A")

# Optional quick view:
table(CRIT3_SUMMARY_DF$OPT_CRIT) # 28, 16%

# REDUCE DF FOR OVERALL SUMMARY
CRIT3_IDS <- CRIT3_SUMMARY_DF %>%
  filter(OPT_CRIT == "CRIT.3") %>%
  pull(ID)

# REDUCE JUST TO CRIT 1 WINNERS
CRIT3_SUMMARY_DF <- CRIT3_SUMMARY_DF %>%
  filter(ID %in% CRIT3_IDS)

# REMOVE WORKING VARIABLE
CRIT3_SUMMARY_DF_RED <- CRIT3_SUMMARY_DF %>%
  select(-c(Found_Crit3_Band))

## ADD BACK INTO THE OVERALL SUMMARY DF 

C1C2C3_OPT_SUM <- rbind(C1C2_OPT_SUM,CRIT3_SUMMARY_DF_RED)

hist(C1C2C3_OPT_SUM$N_Hz)

################################################################################
##################### PARETO CONTIGUITY THRESHOLDED ZONES 

############################
### USER INPUT PARAMETERS 

metrics <- c("CW_POWER_GAMZ_SM", "PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM")
freq_range <- c(0.06, 0.45)
CW_power_z_threshold <- .5
output_dir <- "PARETO_OPTIMAL_FILTERED"

dir.create(output_dir, showWarnings = FALSE)

# Load data
df <- read.csv("CW_POWER_PHASE_FRQ_SUMMARY.csv")
df$ID <- sprintf("%04d", as.numeric(df$ID))
unique_ids <- unique(df$ID)

# PARETO FRONT CALCULATION

find_filtered_pareto_band <- function(df_dyad, metrics, freq_range, z_thresh = 1) {
  df_dyad <- df_dyad %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    filter(!if_any(all_of(metrics), is.na))
  
  if (nrow(df_dyad) == 0) return(NULL)
  
  # Standardize all metrics
  df_scaled <- df_dyad %>%
    mutate(across(all_of(metrics), ~ scale(.)[, 1])) %>%
    arrange(Hz)
  
  X <- df_scaled[, metrics]
  
  # Identify Pareto-optimal rows
  pareto_idx <- which(!apply(X, 1, function(row_i) {
    any(apply(X, 1, function(row_j) all(row_j >= row_i & row_j != row_i)))
  }))
  
  df_scaled$Pareto <- FALSE
  df_scaled$Pareto[pareto_idx] <- TRUE
  
  # Filter: Must be Pareto AND CW_Power > threshold
  df_scaled$Valid <- df_scaled$Pareto & df_scaled[[metrics[1]]] > z_thresh
  
  # Contiguity: Identify contiguous runs
  rle_result <- rle(df_scaled$Valid)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, valid = rle_result$values) %>%
    filter(valid)
  
  if (nrow(run_df) == 0) return(NULL)
  
  # Find the longest contiguous valid block
  run_df <- run_df %>%
    mutate(length = end - start + 1) %>%
    slice_max(length, n = 1)
  
  final_band <- df_scaled[run_df$start:run_df$end, ]
  
  band_summary <- tibble(
    ID = unique(df_dyad$ID),
    Band_Lower = min(final_band$Hz),
    Band_Upper = max(final_band$Hz),
    N_Hz_in_Band = nrow(final_band)
  )
  
  # Plot
  plot_df <- df_scaled %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value")
  
  p <- ggplot(plot_df, aes(x = Hz, y = Value, color = Metric)) +
    geom_line(size = 1) +
    geom_point(data = df_scaled %>% filter(Valid), aes(x = Hz, y = 0), color = "red", size = 2) +
    geom_vline(xintercept = c(band_summary$Band_Lower, band_summary$Band_Upper),
               linetype = "dashed", color = "blue") +
    labs(title = paste("Dyad", unique(df_dyad$ID), "— Filtered Pareto Band"),
         subtitle = paste0("Band: ", round(band_summary$Band_Lower, 3), " - ", 
                           round(band_summary$Band_Upper, 3), " Hz | N = ", 
                           band_summary$N_Hz_in_Band),
         x = "Frequency (Hz)", y = "Z-Score") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = file.path(output_dir, paste0(unique(df_dyad$ID), "_PARETO_BAND.pdf")),
         plot = p, width = 8, height = 5)
  
  return(band_summary)
}


## RUN IT 

band_summary_list <- list()

for (id in unique_ids) {
  dyad_df <- df %>% filter(ID == id)
  
  result <- find_filtered_pareto_band(dyad_df, metrics, freq_range, CW_power_z_threshold)
  
  if (!is.null(result)) {
    band_summary_list[[id]] <- result
  } else {
    message("Skipped dyad: ", id)
  }
}

PARETO_BAND_SUMMARY <- bind_rows(band_summary_list)
write.csv(PARETO_BAND_SUMMARY, "PARETO_BAND_SUMMARY.csv", row.names = FALSE)
################################################################################
##################### USE THE DATA-DRIVEN WEIGHTS & BAND CRITERIA  

## STEP 1: SPLIT DYADS INTO THEIR OWN DFs (ensuring leading 0's stay in place)

df$ID <- sprintf("%04d", as.numeric(df$ID))

unique_ids <- unique(df$ID)

for (id in unique_ids) {
  df_subset <- df %>% filter(ID == id)
  assign(paste0("D_", id), df_subset)
}
## STEP 2  

find_optimal_band_contiguous <- function(df_dyad, 
                                         metrics = c("CW_POWER_GAMZ_SM", "PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM"),
                                         weights = c(0.35, 0.3, 0.35),
                                         threshold_quantile = 0.66,
                                         plot_result = TRUE,
                                         dyad_id = NULL) {
  stopifnot(all(c("Hz", metrics) %in% names(df_dyad)))
  
  if (is.null(dyad_id)) {
    dyad_id <- if ("ID" %in% names(df_dyad)) unique(df_dyad$ID)[1] else "UNKNOWN"
  }
  
  df_dyad <- df_dyad %>% filter(Hz >= 0.06, Hz <= 0.45)
  df_dyad <- df_dyad %>% filter(!if_any(all_of(metrics), is.na))
  if (nrow(df_dyad) == 0) return(NULL)
  
  # Standardize metrics within dyad and compute composite
  df_scaled <- df_dyad %>%
    arrange(Hz) %>%
    mutate(across(all_of(metrics), ~ scale(.)[,1])) %>%
    mutate(
      Composite = weights[1]*get(metrics[1]) + 
        weights[2]*get(metrics[2]) + 
        weights[3]*get(metrics[3])
    )
  
  threshold <- quantile(df_scaled$Composite, probs = threshold_quantile)
  df_scaled <- df_scaled %>%
    mutate(High = Composite >= threshold)
  
  rle_result <- rle(df_scaled$High)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  run_df <- tibble(start = starts, end = ends, is_high = rle_result$values) %>%
    filter(is_high)
  
  if (nrow(run_df) == 0) return(NULL)
  
  run_df <- run_df %>%
    rowwise() %>%
    mutate(
      avg_composite = mean(df_scaled$Composite[start:end], na.rm = TRUE),
      peak_index = start + which.max(df_scaled$Composite[start:end]) - 1
    ) %>%
    ungroup()
  
  best_run <- run_df %>% slice_max(avg_composite, n = 1)
  best_start <- best_run$start
  best_end <- best_run$end
  peak_idx <- best_run$peak_index
  
  band_range <- range(df_scaled$Hz[best_start:best_end])
  peak_freq <- df_scaled$Hz[peak_idx]
  
  # PIVOT for multi-line plotting
  plot_df <- df_scaled %>%
    select(Hz, all_of(metrics), Composite) %>%
    pivot_longer(cols = -Hz, names_to = "Metric", values_to = "Value")
  
  if (plot_result) {
    p <- ggplot() +
      geom_line(data = plot_df, aes(x = Hz, y = Value, color = Metric), linewidth = 1) +
      geom_vline(xintercept = band_range, linetype = "dashed", color = "red") +
      geom_vline(xintercept = peak_freq, linetype = "solid", color = "blue") +
      labs(title = paste("Dyad", dyad_id, "- Optimized Synchrony Band"),
           subtitle = paste0("Band: ", round(band_range[1], 3), " - ", round(band_range[2], 3), 
                             " Hz | Peak: ", round(peak_freq, 3), " Hz"),
           x = "Frequency (Hz)", y = "Z-scored Value", color = "Metric") +
      scale_color_manual(values = c(
        "Composite" = "black",
        "CW_POWER_GAMZ_SM" = "darkgreen",
        "PHS_VAR_GAMZ_SM" = "steelblue",
        "PHS_ZERO_CROSS_GAMZ_SM" = "purple"
      )) +
      theme_minimal(base_size = 14)
    
    dir.create("CONTIGUOUS_OPT_FRQ_PLOTS", showWarnings = FALSE)
    ggsave(paste0("CONTIGUOUS_OPT_FRQ_PLOTS/", dyad_id, "_CONTIGUOUS_BAND.pdf"),
           plot = p, width = 8, height = 5)
  }
  
  return(tibble(
    ID = dyad_id,
    OPT_PEAK_Hz = peak_freq,
    OPT_LB_Hz = band_range[1],
    OPT_UB_Hz = band_range[2]
  ))
}

## STEP 3: APPLY FUNCTION
band_summary_list <- vector("list", length(unique_ids))

for (i in seq_along(unique_ids)) {
  id <- unique_ids[i]
  dyad_df <- get(paste0("D_", id))
  
  band_summary_list[[i]] <- find_optimal_band_contiguous(
    df_dyad = dyad_df,
    metrics = c("CW_POWER_GAMZ_SM", "PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM"),
    weights = c(0.35, .30, 0.35),
    threshold_quantile = 0.66,
    plot_result = TRUE,
    dyad_id = id
  )
}

## CREATE A SUMMARY TABLE 
OPT_BANDS_TABLE <- bind_rows(band_summary_list)


################################################################################
##################### SAVE DATA 

write.csv(OPT_BANDS_TABLE, "DYNAMIC_OPTIMIZED_FRQ_SUMMARY.csv", row.names = FALSE)


#### ARCHIVE 

################################################################################
##################### ARRIVE AT DATA-DRIVEN COMPOSITE & FREQUENCY BAND CRIERIA


############## QUANTIFY THE BEST MODEL WEIGHTS FOR SYNCHRONY COMPOSITE


# Step 1: Create grid of all possible parameter weights
weight_grid <- expand.grid(
  w1 = seq(0, 1, by = 0.05),
  w2 = seq(0, 1, by = 0.05),
  w3 = seq(0, 1, by = 0.05)
) %>%
  filter(abs(w1 + w2 + w3 - 1) < 1e-6)  # Only keep combos that sum to 1

# Step 2: Define function to evaluate SSE for each weight set
evaluate_weights <- function(df_dyad, id, metrics = c("CW_POWER_GAMZ_SM", "PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM")) {
  # Drop NAs and scale
  df_scaled <- df_dyad %>%
    drop_na(all_of(metrics)) %>%
    mutate(across(all_of(metrics), ~ scale(.)[, 1]))
  
  if (nrow(df_scaled) == 0) return(NULL)
  
  # For each row of the grid, compute SSE
  results <- weight_grid %>%
    mutate(
      SSE = pmap_dbl(list(w1, w2, w3), function(w1, w2, w3) {
        composite <- w1 * df_scaled[[metrics[1]]] +
          w2 * df_scaled[[metrics[2]]] +
          w3 * df_scaled[[metrics[3]]]
        
        sum((composite - df_scaled[[metrics[1]]])^2 +
              (composite - df_scaled[[metrics[2]]])^2 +
              (composite - df_scaled[[metrics[3]]])^2, na.rm = TRUE)
      }),
      ID = id
    )
  
  return(results)
}

# Step 3: Loop through all dyads
all_fit_results <- map_dfr(unique_ids, function(id) {
  dyad_df <- get(paste0("D_", id))
  evaluate_weights(dyad_df, id)
})

# Optional: View best weight set per dyad
best_weights_per_dyad <- all_fit_results %>%
  group_by(ID) %>%
  slice_min(SSE, n = 1) %>%
  arrange(ID)

print(best_weights_per_dyad)

table(best_weights_per_dyad$w1) #.35
table(best_weights_per_dyad$w2) #.30
table(best_weights_per_dyad$w3) #.35


################################################################################
######################### QUANTIFY THRESHOLD EACH DOMAIN OF SYCHRONY 

## Cross-wavelet Power 

#######################################################
### USER INPUT PARAMETERS 

metric <- "CW_POWER_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

#######################################################
### STEP 1: Generate sample-level surrogate with smoothing

generate_sample_surrogate <- function(df_sample, metric) {
  surrogate_sample <- df_sample %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    mutate(Surrogate = sample(.[[metric]], replace = TRUE)) %>%
    arrange(Hz)
  
  # Ensure Hz values are numeric and matched
  surrogate_sample$Hz <- round(surrogate_sample$Hz, 4)
  
  # Apply smoothing (rolling average)
  surrogate_sample <- surrogate_sample %>%
    mutate(Smoothed = zoo::rollapply(Surrogate, width = smooth_window, fill = NA, align = "center", FUN = mean)) %>%
    select(Hz, Smoothed) %>%
    rename(Surrogate = Smoothed)
  
  return(surrogate_sample)
}

#######################################################
### STEP 2: Run surrogates and calculate p-values 

run_surrogates_for_dyad <- function(df_dyad, df_sample, metric, n_surrogates) {
  df_dyad <- df_dyad %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    arrange(Hz)
  
  if (nrow(df_dyad) == 0) {
    message("Skipping ID ", unique(df_dyad$ID), ": no valid data after filtering.")
    return(NULL)
  }
  
  Hz_vals <- round(df_dyad$Hz, 4)
  observed_vals <- df_dyad[[metric]]
  
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df_sample, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # One-tailed test: proportion of surrogate values >= observed
  p_vals <- sapply(1:length(Hz_vals), function(i) {
    mean(surrogate_matrix[i, ] >= observed_vals[i], na.rm = TRUE)
  })
  
  result <- tibble(
    ID = unique(df_dyad$ID),
    Hz = Hz_vals,
    Observed = observed_vals,
    P_value = p_vals
  )
  
  return(result)
}

#######################################################
### STEP 3: Identify contiguous significant bands 

find_all_significant_bands <- function(p_result, alpha = 0.05) {
  p_result <- p_result %>%
    arrange(Hz) %>%
    mutate(Significant = P_value < alpha)
  
  rle_result <- rle(p_result$Significant)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, sig = rle_result$values) %>%
    filter(sig)
  
  if (nrow(run_df) == 0) return(NULL)
  
  bands <- lapply(1:nrow(run_df), function(i) {
    rows <- p_result[run_df$start[i]:run_df$end[i], ]
    tibble(
      ID = rows$ID[1],
      Band_Lower = min(rows$Hz),
      Band_Upper = max(rows$Hz),
      N_Frequencies = nrow(rows),
      Min_P = min(rows$P_value, na.rm = TRUE),
      Peak_Hz = rows$Hz[which.min(rows$P_value)]
    )
  })
  
  return(bind_rows(bands))
}

#######################################################
### STEP 4: Run all and plot 

if (!dir.exists("CW_POWER_SURROGATE_BANDS")) {
  dir.create("CW_POWER_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("CW_POWER Surrogate Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("CW_POWER_SURROGATE_BANDS/", id, "_CW_POWER_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}

#######################################################
### STEP 5: MAKE SUMMARY AND SAVE IT 

SURROGATE_CW_POWER_SUMMARY <- bind_rows(summary_table)

write.csv(SURROGATE_CW_POWER_SUMMARY, "SURROGATE_CW_POWER_SUMMARY.csv", row.names = FALSE)

## Phase Zero-Cross

#######################################################
### USER INPUT PARAMETERS 

metric <- "PHS_ZERO_CROSS_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

#######################################################
### STEP 1: Generate sample-level surrogate with smoothing

generate_sample_surrogate <- function(df_sample, metric) {
  surrogate_sample <- df_sample %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    mutate(Surrogate = sample(.[[metric]], replace = TRUE)) %>%
    arrange(Hz)
  
  # Ensure Hz values are numeric and matched
  surrogate_sample$Hz <- round(surrogate_sample$Hz, 4)
  
  # Apply smoothing (rolling average)
  surrogate_sample <- surrogate_sample %>%
    mutate(Smoothed = zoo::rollapply(Surrogate, width = smooth_window, fill = NA, align = "center", FUN = mean)) %>%
    select(Hz, Smoothed) %>%
    rename(Surrogate = Smoothed)
  
  return(surrogate_sample)
}

#######################################################
### STEP 2: Run surrogates and calculate p-values 

run_surrogates_for_dyad <- function(df_dyad, df_sample, metric, n_surrogates) {
  df_dyad <- df_dyad %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    arrange(Hz)
  
  if (nrow(df_dyad) == 0) {
    message("Skipping ID ", unique(df_dyad$ID), ": no valid data after filtering.")
    return(NULL)
  }
  
  Hz_vals <- round(df_dyad$Hz, 4)
  observed_vals <- df_dyad[[metric]]
  
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df_sample, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # One-tailed test: proportion of surrogate values >= observed
  p_vals <- sapply(1:length(Hz_vals), function(i) {
    mean(surrogate_matrix[i, ] >= observed_vals[i], na.rm = TRUE)
  })
  
  result <- tibble(
    ID = unique(df_dyad$ID),
    Hz = Hz_vals,
    Observed = observed_vals,
    P_value = p_vals
  )
  
  return(result)
}

#######################################################
### STEP 3: Identify contiguous significant bands 

find_all_significant_bands <- function(p_result, alpha = 0.05) {
  p_result <- p_result %>%
    arrange(Hz) %>%
    mutate(Significant = P_value < alpha)
  
  rle_result <- rle(p_result$Significant)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, sig = rle_result$values) %>%
    filter(sig)
  
  if (nrow(run_df) == 0) return(NULL)
  
  bands <- lapply(1:nrow(run_df), function(i) {
    rows <- p_result[run_df$start[i]:run_df$end[i], ]
    tibble(
      ID = rows$ID[1],
      Band_Lower = min(rows$Hz),
      Band_Upper = max(rows$Hz),
      N_Frequencies = nrow(rows),
      Min_P = min(rows$P_value, na.rm = TRUE),
      Peak_Hz = rows$Hz[which.min(rows$P_value)]
    )
  })
  
  return(bind_rows(bands))
}

#######################################################
### STEP 4: Run all and plot 


if (!dir.exists("PHS_ZERO_CROSS_SURROGATE_BANDS")) {
  dir.create("PHS_ZERO_CROSS_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("Phase Zero Cross Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("PHS_ZERO_CROSS_SURROGATE_BANDS/", id, "_PHS_ZERO_CROSS_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}


#######################################################
### STEP 5: MAKE SUMMARY AND SAVE IT 

SURROGATE_PHS_ZERO_CROSS_SUMMARY <- bind_rows(summary_table)
write.csv(SURROGATE_PHS_ZERO_CROSS_SUMMARY, "SURROGATE_PHS_ZERO_CROSS_SUMMARY.csv", row.names = FALSE)


## PHASE VARIANCE 


#######################################################
### USER INPUT PARAMETERS 

metric <- "PHS_VAR_GAMZ_SM"            
n_surrogates <- 1000
alpha_level <- 0.05
freq_range <- c(0.04, 0.50)
smooth_window <- 7   # Must be odd number

#######################################################
### STEP 1: Generate sample-level surrogate with smoothing

generate_sample_surrogate <- function(df_sample, metric) {
  surrogate_sample <- df_sample %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    mutate(Surrogate = sample(.[[metric]], replace = TRUE)) %>%
    arrange(Hz)
  
  # Ensure Hz values are numeric and matched
  surrogate_sample$Hz <- round(surrogate_sample$Hz, 4)
  
  # Apply smoothing (rolling average)
  surrogate_sample <- surrogate_sample %>%
    mutate(Smoothed = zoo::rollapply(Surrogate, width = smooth_window, fill = NA, align = "center", FUN = mean)) %>%
    select(Hz, Smoothed) %>%
    rename(Surrogate = Smoothed)
  
  return(surrogate_sample)
}

#######################################################
### STEP 2: Run surrogates and calculate p-values 

run_surrogates_for_dyad <- function(df_dyad, df_sample, metric, n_surrogates) {
  df_dyad <- df_dyad %>%
    filter(Hz >= freq_range[1], Hz <= freq_range[2]) %>%
    drop_na(all_of(metric)) %>%
    arrange(Hz)
  
  if (nrow(df_dyad) == 0) {
    message("Skipping ID ", unique(df_dyad$ID), ": no valid data after filtering.")
    return(NULL)
  }
  
  Hz_vals <- round(df_dyad$Hz, 4)
  observed_vals <- df_dyad[[metric]]
  
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df_sample, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # One-tailed test: proportion of surrogate values >= observed
  p_vals <- sapply(1:length(Hz_vals), function(i) {
    mean(surrogate_matrix[i, ] >= observed_vals[i], na.rm = TRUE)
  })
  
  result <- tibble(
    ID = unique(df_dyad$ID),
    Hz = Hz_vals,
    Observed = observed_vals,
    P_value = p_vals
  )
  
  return(result)
}

#######################################################
### STEP 3: Identify contiguous significant bands 

find_all_significant_bands <- function(p_result, alpha = 0.05) {
  p_result <- p_result %>%
    arrange(Hz) %>%
    mutate(Significant = P_value < alpha)
  
  rle_result <- rle(p_result$Significant)
  ends <- cumsum(rle_result$lengths)
  starts <- c(1, head(ends + 1, -1))
  
  run_df <- tibble(start = starts, end = ends, sig = rle_result$values) %>%
    filter(sig)
  
  if (nrow(run_df) == 0) return(NULL)
  
  bands <- lapply(1:nrow(run_df), function(i) {
    rows <- p_result[run_df$start[i]:run_df$end[i], ]
    tibble(
      ID = rows$ID[1],
      Band_Lower = min(rows$Hz),
      Band_Upper = max(rows$Hz),
      N_Frequencies = nrow(rows),
      Min_P = min(rows$P_value, na.rm = TRUE),
      Peak_Hz = rows$Hz[which.min(rows$P_value)]
    )
  })
  
  return(bind_rows(bands))
}

#######################################################
### STEP 4: Run all and plot 


if (!dir.exists("PHS_VARIANCE_SURROGATE_BANDS")) {
  dir.create("PHS_VARIANCE_SURROGATE_BANDS")
}

summary_table <- list()

for (id in unique_ids) {
  df_dyad <- dyad_list[[id]]
  
  p_results <- run_surrogates_for_dyad(df_dyad, df_sample = df, metric, n_surrogates)
  if (is.null(p_results)) next
  
  bands <- find_all_significant_bands(p_results, alpha = alpha_level)
  if (is.null(bands)) {
    message("Skipping ID ", id, ": no significant band found.")
    next
  }
  
  summary_table[[id]] <- bands
  
  # Re-generate surrogate matrix for plotting shaded ribbon
  Hz_vals <- p_results$Hz
  surrogate_matrix <- matrix(NA, nrow = length(Hz_vals), ncol = n_surrogates)
  for (i in 1:n_surrogates) {
    surrogate <- generate_sample_surrogate(df, metric)
    surrogate <- surrogate %>% filter(Hz %in% Hz_vals)
    surrogate_matrix[, i] <- surrogate$Surrogate[match(Hz_vals, surrogate$Hz)]
  }
  
  # Calculate 95th percentile (or use mean/SD or max)
  surrogate_df <- tibble(
    Hz = Hz_vals,
    Surrogate_Upper = apply(surrogate_matrix, 1, function(x) quantile(x, 0.95, na.rm = TRUE)),
    Surrogate_Lower = apply(surrogate_matrix, 1, function(x) quantile(x, 0.05, na.rm = TRUE))
  )
  
  plot_df <- p_results %>% left_join(surrogate_df, by = "Hz")
  
  p <- ggplot(plot_df, aes(x = Hz)) +
    geom_ribbon(aes(ymin = Surrogate_Lower, ymax = Surrogate_Upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = Observed), color = "steelblue", size = 1) +
    geom_point(data = plot_df %>% filter(P_value < alpha_level),
               aes(y = Observed), color = "red", size = 2) +
    geom_vline(data = bands, aes(xintercept = Band_Lower), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Band_Upper), linetype = "dashed", color = "darkred") +
    geom_vline(data = bands, aes(xintercept = Peak_Hz), linetype = "solid", color = "blue") +
    labs(title = paste("Phase Variance Significance — Dyad", id),
         subtitle = "Gray shaded area = surrogate 5th–95th percentile",
         x = "Hz", y = "CW_POWER") +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste0("PHS_VARIANCE_SURROGATE_BANDS/", id, "_PHS_VARIANCE_SURROGATE_BAND.pdf"),
         plot = p, width = 8, height = 5)
}


#######################################################
### STEP 5: MAKE SUMMARY AND SAVE IT 

SURROGATE_PHS_VARIANCE_SUMMARY <- bind_rows(summary_table)
write.csv(SURROGATE_PHS_VARIANCE_SUMMARY, "SURROGATE_PHS_VARIANCE_SUMMARY.csv", row.names = FALSE)



