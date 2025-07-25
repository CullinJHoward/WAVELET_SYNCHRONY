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

# Load data

df <- read.csv("CW_POWER_PHASE_FRQ_SUMMARY.csv")
df$ID <- sprintf("%04d", as.numeric(df$ID))
SURROGATE_CW_POWER_SUMMARY <- read.csv("SURROGATE_CW_POWER_SUMMARY.csv")
SURROGATE_CW_POWER_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_CW_POWER_SUMMARY$ID))
SURROGATE_PHS_ZERO_CROSS_SUMMARY <- read.csv("SURROGATE_PHS_ZERO_CROSS_SUMMARY.csv")
SURROGATE_PHS_ZERO_CROSS_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_PHS_ZERO_CROSS_SUMMARY$ID))
SURROGATE_PHS_VARIANCE_SUMMARY <- read.csv("SURROGATE_PHS_VARIANCE_SUMMARY.csv")
SURROGATE_PHS_VARIANCE_SUMMARY$ID <- sprintf("%04d", as.numeric(SURROGATE_PHS_VARIANCE_SUMMARY$ID))


################################################################################
################## DATA WRANGLE

# GET THE SURROGATE BAND DATA INTO LONG FORMAT AND INTO A JOINT DF 

# Function to mark Hzs inside bands for a given domain
flag_significant_hz <- function(df, summary_df, var_name) {
  df[[var_name]] <- 0
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

df_annotated <- df %>%
  flag_significant_hz(SURROGATE_CW_POWER_SUMMARY, "CW_Shz") %>%
  flag_significant_hz(SURROGATE_PHS_ZERO_CROSS_SUMMARY, "PH_ZC_Shz") %>%
  flag_significant_hz(SURROGATE_PHS_VARIANCE_SUMMARY, "PH_VAR_Shz")


################################################################################
## DYNAMIC CRITERIA #1 & #2

# Specify frequency range for candidate bands
freq_min <- 0.06
freq_max <- 0.5

# Create output folder if it doesn't exist
if (!dir.exists("FRQ_OPT_CRT1_CRIT2_PLOTS_A")) {
  dir.create("FRQ_OPT_CRT1_CRIT2_PLOTS_A")
}

# Helper: Find contiguous runs of at least 3 TRUEs
find_contiguous_runs <- function(logical_vec, hz_vec, min_len = 3) {
  rle_obj <- rle(logical_vec)
  ends <- cumsum(rle_obj$lengths)
  starts <- c(1, head(ends + 1, -1))
  bands <- list()
  for (i in seq_along(rle_obj$values)) {
    if (rle_obj$values[i] && rle_obj$lengths[i] >= min_len) {
      idx <- starts[i]:ends[i]
      bands[[length(bands) + 1]] <- range(hz_vec[idx])
    }
  }
  bands
}

# Main loop
summary_list <- list()
unique_ids <- unique(df_annotated$ID)

for (id in unique_ids) {
  df_sub <- df_annotated %>% filter(ID == id) %>% arrange(Hz)
  
  # --- QC check: flag if >45% of CW_POWER_GAMZ_SM <= -2
  frac_low_power <- mean(df_sub$CW_POWER_GAMZ_SM <= -2, na.rm = TRUE)
  qc_fail <- frac_low_power > 0.45
  
  candidate_bands <- list()
  band_info <- list()
  
  # --- Criteria 1: Power significant, phase >= -1 (either phase metric)
  for (phase_col in c("PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM")) {
    crit2_logical <- (df_sub$CW_Shz == 1) & (df_sub[[phase_col]] >= -1)
    bands <- find_contiguous_runs(crit2_logical, df_sub$Hz, min_len = 3)
    for (b in bands) {
      idx <- which(df_sub$Hz >= b[1] & df_sub$Hz <= b[2])
      band_df <- df_sub[idx, ]
      n_Hz <- length(idx)
      # Only keep band if it is entirely within the desired frequency range
      if (min(band_df$Hz) < freq_min || max(band_df$Hz) > freq_max) next
      # Only keep band if all three metrics are >= -1 at every Hz
      if (all(band_df$CW_POWER_GAMZ_SM >= -1, na.rm = TRUE) &&
          all(band_df$PHS_VAR_GAMZ_SM >= -1, na.rm = TRUE) &&
          all(band_df$PHS_ZERO_CROSS_GAMZ_SM >= -1, na.rm = TRUE)) {
        band_info[[length(band_info) + 1]] <- list(
          type = "CRIT.1",
          phase_metric = phase_col,
          xmin = b[1], xmax = b[2],
          n_Hz = n_Hz,
          mean_CW = mean(band_df$CW_POWER_GAMZ_SM, na.rm = TRUE),
          mean_VAR = mean(band_df$PHS_VAR_GAMZ_SM, na.rm = TRUE),
          mean_ZC = mean(band_df$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
        )
      }
    }
  }
  
  # --- Criteria 2: Phase significant, power >= -1 (for each phase metric)
  for (phase_col in c("PH_VAR_Shz", "PH_ZC_Shz")) {
    crit3_logical <- (df_sub[[phase_col]] == 1) & (df_sub$CW_POWER_GAMZ_SM >= -1)
    bands <- find_contiguous_runs(crit3_logical, df_sub$Hz, min_len = 3)
    for (b in bands) {
      idx <- which(df_sub$Hz >= b[1] & df_sub$Hz <= b[2])
      band_df <- df_sub[idx, ]
      n_Hz <- length(idx)
      # Only keep band if it is entirely within the desired frequency range
      if (min(band_df$Hz) < freq_min || max(band_df$Hz) > freq_max) next
      # Only keep band if all three metrics are >= -1 at every Hz
      if (all(band_df$CW_POWER_GAMZ_SM >= -1, na.rm = TRUE) &&
          all(band_df$PHS_VAR_GAMZ_SM >= -1, na.rm = TRUE) &&
          all(band_df$PHS_ZERO_CROSS_GAMZ_SM >= -1, na.rm = TRUE)) {
        band_info[[length(band_info) + 1]] <- list(
          type = "CRIT.2",
          phase_metric = ifelse(phase_col == "PH_VAR_Shz", "PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM"),
          xmin = b[1], xmax = b[2],
          n_Hz = n_Hz,
          mean_CW = mean(band_df$CW_POWER_GAMZ_SM, na.rm = TRUE),
          mean_VAR = mean(band_df$PHS_VAR_GAMZ_SM, na.rm = TRUE),
          mean_ZC = mean(band_df$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
        )
      }
    }
  }
  
  # --- Collate all candidate bands
  band_df <- if (length(band_info) > 0) bind_rows(lapply(band_info, as.data.frame)) else NULL
  CAND_BANDS <- if (!is.null(band_df)) nrow(band_df) else 0
  
  # --- Select optimal band: weighted mean by n_Hz (linear)
  if (!qc_fail && !is.null(band_df) && nrow(band_df) > 0) {
    band_df <- band_df %>%
      mutate(
        weighted_CW = mean_CW * n_Hz,
        reference_phase_mean = dplyr::case_when(
          type == "CRIT.1" & phase_metric == "PHS_VAR_GAMZ_SM" ~ mean_VAR,
          type == "CRIT.1" & phase_metric == "PHS_ZERO_CROSS_GAMZ_SM" ~ mean_ZC,
          type == "CRIT.2" & phase_metric == "PHS_VAR_GAMZ_SM" ~ mean_VAR,
          type == "CRIT.2" & phase_metric == "PHS_ZERO_CROSS_GAMZ_SM" ~ mean_ZC
        ),
        weighted_reference_phase = reference_phase_mean * n_Hz
      )
    band_df <- band_df %>%
      arrange(desc(weighted_CW), desc(weighted_reference_phase))
    optimal <- band_df[1, ]
    OPT_CAT <- optimal$type
    OPT_LB <- optimal$xmin
    OPT_UB <- optimal$xmax
    N_Hz <- optimal$n_Hz
  } else if (qc_fail) {
    OPT_CAT <- "QC_FAIL"
    OPT_LB <- "QC_FAIL"
    OPT_UB <- "QC_FAIL"
    N_Hz <- "QC_FAIL"
    optimal <- NULL
  } else {
    OPT_CAT <- NA
    OPT_LB <- NA
    OPT_UB <- NA
    N_Hz <- NA
    optimal <- NULL
  }
  
  # --- Plotting (unchanged)
  p <- ggplot(df_sub, aes(x = Hz)) +
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER"), size = 1) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VAR"), size = 1) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "ZERO_CROSS"), size = 1) +
    geom_point(data = df_sub %>% filter(CW_Shz == 1), aes(y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_VAR_Shz == 1), aes(y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_ZC_Shz == 1), aes(y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    scale_color_manual(values = c("CW_POWER" = "black", "PHASE_VAR" = "blue", "ZERO_CROSS" = "darkgreen"), name = "Signal") +
    ggtitle(paste0("ID: ", id)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("GAMZ (Z-Score)") +
    xlab("Hz")
  
  # Overlay candidate bands with legend (unchanged)
  if (!is.null(band_df) && nrow(band_df) > 0) {
    band_df$fill_label <- dplyr::case_when(
      band_df$type == "CRIT.1" ~ "Criteria 1: CW_POWER",
      band_df$type == "CRIT.2" & band_df$phase_metric == "PHS_VAR_GAMZ_SM" ~ "Criteria 2a: Phs_Var",
      band_df$type == "CRIT.2" & band_df$phase_metric == "PHS_ZERO_CROSS_GAMZ_SM" ~ "Criteria 2b: Phs_ZeroCross"
    )
    p <- p +
      geom_rect(
        data = band_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_label),
        inherit.aes = FALSE,
        alpha = 0.2
      )
    if (!qc_fail && !is.null(optimal)) {
      star_hz <- mean(c(optimal$xmin, optimal$xmax))
      y_max <- max(df_sub$CW_POWER_GAMZ_SM, df_sub$PHS_VAR_GAMZ_SM, df_sub$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
      p <- p + annotate("point", x = star_hz, y = y_max + 0.3, shape = 8, size = 4, color = "black")
    }
    p <- p + scale_fill_manual(
      name = "Candidate Band",
      values = c(
        "Criteria 1: CW_POWER" = "plum1",
        "Criteria 2a: Phs_Var" = "blue",
        "Criteria 2b: Phs_ZeroCross" = "lightgreen"
      )
    )
  }
  
  # Save plot
  ggsave(paste0("FRQ_OPT_CRT1_CRIT2_PLOTS_A/", id, "_DYN_CRT1&2_OPT_BAND_A.pdf"), p, width = 8, height = 5)
  
  # --- Add to summary
  summary_list[[id]] <- data.frame(
    ID = id,
    OPT_LB = OPT_LB,
    OPT_UB = OPT_UB,
    N_Hz = N_Hz,
    OPT_CRIT = OPT_CAT,
    CAND_BANDS = CAND_BANDS
  )
}

# CREATE SUMMARY
C1C2_OPT_SUMMARY <- do.call(rbind, summary_list)

################################################################################
##  CRITERIA #3

# We will just use their peak CW_POWER, bounded by the average band width

# SEPARATE OUT THE LAST FEW PARTICIPANTS 
NO_C1C2_ID <- C1C2_OPT_SUMMARY %>%
  filter(is.na(OPT_CRIT)) %>%
  pull(ID)

df_CRIT3 <- df_annotated %>%
  filter(ID %in% NO_C1C2_ID)

df_CRIT3 <- subset(df_CRIT3, Hz > .06 & Hz < .5)

## FIX NA ISSUES IN SUMMARY DF 
C1C2_OPT_SUMMARY$N_Hz[C1C2_OPT_SUMMARY$N_Hz == "QC_FAIL"] <- NA
C1C2_OPT_SUMMARY$OPT_LB[C1C2_OPT_SUMMARY$OPT_LB == "QC_FAIL"] <- NA
C1C2_OPT_SUMMARY$OPT_UB[C1C2_OPT_SUMMARY$OPT_UB == "QC_FAIL"] <- NA

C1C2_OPT_SUMMARY$OPT_UB <- as.numeric(C1C2_OPT_SUMMARY$OPT_UB )
C1C2_OPT_SUMMARY$OPT_LB <- as.numeric(C1C2_OPT_SUMMARY$OPT_LB)
C1C2_OPT_SUMMARY$N_Hz <- as.numeric(C1C2_OPT_SUMMARY$N_Hz)

#### THOSE THAT DON'T FIT WILL USE THE FULL SPECTRUM 

C1C2_OPT_SUMMARY$OPT_CRIT[is.na(C1C2_OPT_SUMMARY$OPT_CRIT)] <- "CRIT.3"

#### CREATE OPTIMIZED HZ BANDS AROUND EACH CRITERA 4 DYADS HIGHEST POWER MEAN
## THIS IS JUST THEIR PEAK, WITH THE AVERAGE BAND WIDTH (FROM THE SAMPLE) AROUND IT

# AVERAGE NUMBER OF HZ IN EACH FREQUENCY OPTIMZED BAND 
mean(C1C2_OPT_SUMMARY$N_Hz, na.rm = TRUE) # 7.44

# COMPUTE EACH DYADS CW_POWER PEAK MEAN  
MAX_CW_POWER <- df_CRIT3 %>%
  group_by(ID) %>%
  filter(CW_POWER_GAMZ == max(CW_POWER_GAMZ, na.rm = TRUE)) %>%
  select(ID, Hz, CW_POWER_GAMZ)

# IDENTIFY THE BAND AROUND AROUND EACH PERSONS PEAK 

# Initialize vectors
OPT_LB <- numeric(nrow(MAX_CW_POWER))
OPT_UB <- numeric(nrow(MAX_CW_POWER))

# Loop through each row
for (i in seq_len(nrow(MAX_CW_POWER))) {
  
  id_i <- MAX_CW_POWER$ID[i]
  match_hz <- MAX_CW_POWER$Hz[i]  # Exact Hz to match
  
  # Get the sorted Hz values for this person
  hz_vals <- df_CRIT3 %>%
    filter(ID == id_i) %>%
    arrange(Hz) %>%
    pull(Hz)
  
  # Find the index of the exact match
  match_idx <- which(hz_vals == match_hz)
  
  if (length(match_idx) == 0) {
    # Handle case where match is not found
    OPT_LB[i] <- NA
    OPT_UB[i] <- NA
  } else {
    # Get bounds (3 before and after)
    lb_idx <- max(1, match_idx - 3)
    ub_idx <- min(length(hz_vals), match_idx + 3)
    
    OPT_LB[i] <- hz_vals[lb_idx]
    OPT_UB[i] <- hz_vals[ub_idx]
  }
}

# Assign to data frame
MAX_CW_POWER$OPT_LB <- OPT_LB
MAX_CW_POWER$OPT_UB <- OPT_UB

# Set to 6 rows by default (7.4 was the mean)

MAX_CW_POWER$N_Hz <- 6

##### PLOT THE BANDS

# Create output folder if not present
if (!dir.exists("FRQ_OPT_CRT3_PLOTS_A")) {
  dir.create("FRQ_OPT_CRT3_PLOTS_A")
}

# Merge band bounds with df_CRIT3 for plotting
plot_df <- df_CRIT3
plot_df <- merge(plot_df, MAX_CW_POWER[, c("ID", "OPT_LB", "OPT_UB")], by = "ID", all.x = TRUE)

for (id in unique(plot_df$ID)) {
  df_sub <- plot_df[plot_df$ID == id, ]
  opt_lb <- unique(df_sub$OPT_LB)
  opt_ub <- unique(df_sub$OPT_UB)
  # Significant points
  sig_cw  <- df_sub[df_sub$CW_Shz == 1, ]
  sig_var <- df_sub[df_sub$PH_VAR_Shz == 1, ]
  sig_zc  <- df_sub[df_sub$PH_ZC_Shz == 1, ]
  # Shaded band
  band_rect <- if (!is.na(opt_lb) && !is.na(opt_ub)) {
    data.frame(xmin = opt_lb, xmax = opt_ub)
  } else {
    NULL
  }
  # Star at center of band
  star_point <- if (!is.na(opt_lb) && !is.na(opt_ub)) {
    center_hz <- mean(c(opt_lb, opt_ub))
    y_max <- max(df_sub$CW_POWER_GAMZ_SM, df_sub$PHS_VAR_GAMZ_SM, df_sub$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
    data.frame(Hz = center_hz, y = y_max + 0.2)
  } else {
    NULL
  }
  p <- ggplot(df_sub, aes(x = Hz)) +
    # Shaded optimized band
    {if (!is.null(band_rect)) geom_rect(
      data = band_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gold", alpha = 0.3
    )} +
    # Metric curves
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER")) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VARIANCE")) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "PHASE_ZERO_CROSS")) +
    # Significant points
    geom_point(data = sig_cw, aes(x = Hz, y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_var, aes(x = Hz, y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_zc, aes(x = Hz, y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    # Star at center
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
      title = paste("CRIT.3 Optimized Band — Dyad", id),
      x = "Hz", y = "Value", color = "Metric"
    ) +
    theme_minimal(base_size = 14)
  ggsave(
    filename = paste0("FRQ_OPT_CRT3_PLOTS_A/", id, "_CRIT3_OPT_BAND_A.pdf"),
    plot = p, width = 8, height = 5
  )
} 


##### ADD TO OVERALL SUMMARY 


# SUBSET JUST THE VARIABLES I WANT
CRIT3_BOUNDS <- MAX_CW_POWER %>%
  select(ID, N_Hz, OPT_LB, OPT_UB)

# UPDATE THE EXISTING SUMMARY WITH MY NEW VALUES
C1C2C3_OPT_SUMMARY <- rows_update(C1C2_OPT_SUMMARY, CRIT3_BOUNDS, by = "ID")


#GET COUNTS
table(C1C2C3_OPT_SUMMARY$OPT_CRIT) 

# CRITERIA #1 = 92, 
# CRITERIA 2 = 58, 
# CRITERIA 3 = 10, 
# QC_FAIL = 15

################################################################################
################# SAVE SUMMARY CSV

write.csv(C1C2C3_OPT_SUMMARY, "FRQ_OPT_CLASSIFICATION_VERSION_A_SUMMARY.csv", row.names = FALSE)

