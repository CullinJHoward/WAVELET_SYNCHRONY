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
## DYNAMIC CRITERIA #1 & #2 (VERSION B)

# Specify frequency range for candidate bands
freq_min <- 0.06
freq_max <- 0.5

# Create output folder if it doesn't exist
if (!dir.exists("FRQ_OPT_CRT1_CRIT2_PLOTS_B")) {
  dir.create("FRQ_OPT_CRT1_CRIT2_PLOTS_B")
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

summary_list <- list()
unique_ids <- unique(df_annotated$ID)

for (id in unique_ids) {
  df_sub <- df_annotated %>% filter(ID == id) %>% arrange(Hz)
  
  # --- QC check: flag if >45% of CW_POWER_GAMZ_SM <= -2
  frac_low_power <- mean(df_sub$CW_POWER_GAMZ_SM <= -2, na.rm = TRUE)
  qc_fail <- frac_low_power > 0.45
  
  candidate_bands <- list()
  band_info <- list()
  
  # --- Criteria 1: CW_POWER significant, at least one phase metric always > -1
  crit1_logical <- (df_sub$CW_Shz == 1) & 
    (apply(df_sub[,c("PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM")], 1, function(x) any(x > -1)))
  bands1 <- find_contiguous_runs(crit1_logical, df_sub$Hz, min_len = 3)
  for (b in bands1) {
    idx <- which(df_sub$Hz >= b[1] & df_sub$Hz <= b[2])
    band_df <- df_sub[idx, ]
    n_Hz <- length(idx)
    # Only keep band if it is entirely within the desired frequency range
    if (min(band_df$Hz) < freq_min || max(band_df$Hz) > freq_max) next
    # At least one phase metric always > -1
    if (all(apply(band_df[,c("PHS_VAR_GAMZ_SM", "PHS_ZERO_CROSS_GAMZ_SM")], 1, function(x) any(x > -1)))) {
      band_info[[length(band_info) + 1]] <- list(
        type = "CRIT.1",
        xmin = b[1], xmax = b[2],
        n_Hz = n_Hz,
        mean_CW = mean(band_df$CW_POWER_GAMZ_SM, na.rm = TRUE),
        mean_VAR = mean(band_df$PHS_VAR_GAMZ_SM, na.rm = TRUE),
        mean_ZC = mean(band_df$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
      )
    }
  }
  
  # If no CRIT.1 band, try CRIT.2
  if (length(band_info) == 0) {
    # --- Criteria 2: Phase significant, CW_POWER >= 0
    for (phase_col in c("PH_VAR_Shz", "PH_ZC_Shz")) {
      crit2_logical <- (df_sub[[phase_col]] == 1) & (df_sub$CW_POWER_GAMZ_SM >= 0)
      bands2 <- find_contiguous_runs(crit2_logical, df_sub$Hz, min_len = 3)
      for (b in bands2) {
        idx <- which(df_sub$Hz >= b[1] & df_sub$Hz <= b[2])
        band_df <- df_sub[idx, ]
        n_Hz <- length(idx)
        # Only keep band if it is entirely within the desired frequency range
        if (min(band_df$Hz) < freq_min || max(band_df$Hz) > freq_max) next
        # All CW_POWER >= 0 in band
        if (all(band_df$CW_POWER_GAMZ_SM >= 0, na.rm = TRUE)) {
          band_info[[length(band_info) + 1]] <- list(
            type = "CRIT.2",
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
  }
  
  # --- Collate all candidate bands
  band_df <- if (length(band_info) > 0) bind_rows(lapply(band_info, as.data.frame)) else NULL
  CAND_BANDS <- if (!is.null(band_df)) nrow(band_df) else 0
  
  # --- Select optimal band: longest (most N_Hz rows)
  if (!qc_fail && !is.null(band_df) && nrow(band_df) > 0) {
    band_df <- band_df %>% arrange(desc(n_Hz))
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
  
  # --- Plotting
  p <- ggplot(df_sub, aes(x = Hz)) +
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER"), linewidth = 1) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VAR"), linewidth = 1) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "ZERO_CROSS"), linewidth = 1) +
    geom_point(data = df_sub %>% filter(CW_Shz == 1), aes(y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_VAR_Shz == 1), aes(y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = df_sub %>% filter(PH_ZC_Shz == 1), aes(y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    scale_color_manual(values = c("CW_POWER" = "black", "PHASE_VAR" = "blue", "ZERO_CROSS" = "darkgreen"), name = "Signal") +
    ggtitle(paste0("ID: ", id)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("GAMZ (Z-Score)") +
    xlab("Hz")
  
  # Overlay candidate bands with legend
  if (!is.null(band_df) && nrow(band_df) > 0) {
    # Ensure phase_metric column exists and is character
    if (!"phase_metric" %in% names(band_df)) band_df$phase_metric <- NA_character_
    band_df$phase_metric <- as.character(band_df$phase_metric)
    band_df$fill_label <- dplyr::case_when(
      band_df$type == "CRIT.1" ~ "Criteria 1: CW_POWER",
      band_df$type == "CRIT.2" & band_df$phase_metric == "PH_VAR_Shz" ~ "Criteria 2a: Phs_Var",
      band_df$type == "CRIT.2" & band_df$phase_metric == "PH_ZC_Shz" ~ "Criteria 2b: Phs_ZeroCross",
      TRUE ~ as.character(band_df$type)
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
  ggsave(paste0("FRQ_OPT_CRT1_CRIT2_PLOTS_B/", id, "_DYN_CRT1&2_OPT_BAND_B.pdf"), p, width = 8, height = 5)
  
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
##  CRITERIA #3 (same as Version A)

# SEPARATE OUT THE LAST FEW PARTICIPANTS 
NO_C1C2_ID <- C1C2_OPT_SUMMARY %>%
  filter(is.na(OPT_CRIT)) %>%
  pull(ID)

df_CRIT3 <- df_annotated %>%
  filter(ID %in% NO_C1C2_ID)

df_CRIT3 <- subset(df_CRIT3, Hz > .06 & Hz < .5)

C1C2_OPT_SUMMARY$N_Hz[C1C2_OPT_SUMMARY$N_Hz == "QC_FAIL"] <- NA
C1C2_OPT_SUMMARY$OPT_LB[C1C2_OPT_SUMMARY$OPT_LB == "QC_FAIL"] <- NA
C1C2_OPT_SUMMARY$OPT_UB[C1C2_OPT_SUMMARY$OPT_UB == "QC_FAIL"] <- NA

C1C2_OPT_SUMMARY$OPT_UB <- as.numeric(C1C2_OPT_SUMMARY$OPT_UB )
C1C2_OPT_SUMMARY$OPT_LB <- as.numeric(C1C2_OPT_SUMMARY$OPT_LB)
C1C2_OPT_SUMMARY$N_Hz <- as.numeric(C1C2_OPT_SUMMARY$N_Hz)

C1C2_OPT_SUMMARY$OPT_CRIT[is.na(C1C2_OPT_SUMMARY$OPT_CRIT)] <- "CRIT.3"

mean(C1C2_OPT_SUMMARY$N_Hz, na.rm = TRUE) # 7.44

MAX_CW_POWER <- df_CRIT3 %>%
  group_by(ID) %>%
  filter(CW_POWER_GAMZ == max(CW_POWER_GAMZ, na.rm = TRUE)) %>%
  select(ID, Hz, CW_POWER_GAMZ)

OPT_LB <- numeric(nrow(MAX_CW_POWER))
OPT_UB <- numeric(nrow(MAX_CW_POWER))

for (i in seq_len(nrow(MAX_CW_POWER))) {
  id_i <- MAX_CW_POWER$ID[i]
  match_hz <- MAX_CW_POWER$Hz[i]
  hz_vals <- df_CRIT3 %>%
    filter(ID == id_i) %>%
    arrange(Hz) %>%
    pull(Hz)
  match_idx <- which(hz_vals == match_hz)
  if (length(match_idx) == 0) {
    OPT_LB[i] <- NA
    OPT_UB[i] <- NA
  } else {
    lb_idx <- max(1, match_idx - 3)
    ub_idx <- min(length(hz_vals), match_idx + 3)
    OPT_LB[i] <- hz_vals[lb_idx]
    OPT_UB[i] <- hz_vals[ub_idx]
  }
}

MAX_CW_POWER$OPT_LB <- OPT_LB
MAX_CW_POWER$OPT_UB <- OPT_UB
MAX_CW_POWER$N_Hz <- 6

if (!dir.exists("FRQ_OPT_CRT3_PLOTS_B")) {
  dir.create("FRQ_OPT_CRT3_PLOTS_B")
}

plot_df <- df_CRIT3
plot_df <- merge(plot_df, MAX_CW_POWER[, c("ID", "OPT_LB", "OPT_UB")], by = "ID", all.x = TRUE)

for (id in unique(plot_df$ID)) {
  df_sub <- plot_df[plot_df$ID == id, ]
  opt_lb <- unique(df_sub$OPT_LB)
  opt_ub <- unique(df_sub$OPT_UB)
  sig_cw  <- df_sub[df_sub$CW_Shz == 1, ]
  sig_var <- df_sub[df_sub$PH_VAR_Shz == 1, ]
  sig_zc  <- df_sub[df_sub$PH_ZC_Shz == 1, ]
  band_rect <- if (!is.na(opt_lb) && !is.na(opt_ub)) {
    data.frame(xmin = opt_lb, xmax = opt_ub)
  } else {
    NULL
  }
  star_point <- if (!is.na(opt_lb) && !is.na(opt_ub)) {
    center_hz <- mean(c(opt_lb, opt_ub))
    y_max <- max(df_sub$CW_POWER_GAMZ_SM, df_sub$PHS_VAR_GAMZ_SM, df_sub$PHS_ZERO_CROSS_GAMZ_SM, na.rm = TRUE)
    data.frame(Hz = center_hz, y = y_max + 0.2)
  } else {
    NULL
  }
  p <- ggplot(df_sub, aes(x = Hz)) +
    {if (!is.null(band_rect)) geom_rect(
      data = band_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "gold", alpha = 0.3
    )} +
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER")) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VARIANCE")) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "PHASE_ZERO_CROSS")) +
    geom_point(data = sig_cw, aes(x = Hz, y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_var, aes(x = Hz, y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_zc, aes(x = Hz, y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
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
    filename = paste0("FRQ_OPT_CRT3_PLOTS_B/", id, "_CRIT3_OPT_BAND_B.pdf"),
    plot = p, width = 8, height = 5
  )
}


## COMBINE THE FULL SUMMARY LIST 
CRIT3_BOUNDS <- MAX_CW_POWER %>%
  select(ID, N_Hz, OPT_LB, OPT_UB)


C1C2C3_OPT_SUMMARY <- rows_update(C1C2_OPT_SUMMARY, CRIT3_BOUNDS, by = "ID")

## GET FREQUENCIES 
table(C1C2C3_OPT_SUMMARY$OPT_CRIT)
# CRITERIA 1 = 109
# CRITERIA 2 = 29
# CRITERIA 3 = 22
# QC FAIL = 15

write.csv(C1C2C3_OPT_SUMMARY, "FRQ_OPT_CLASSIFICATION_VERSION_B_SUMMARY.csv", row.names = FALSE)


################################################################################
################### COMPARE VERSION A and VERSION B OPT BANDS 

df_A <- read.csv("FRQ_OPT_CLASSIFICATION_VERSION_A_SUMMARY.csv")
df_A$ID <- sprintf("%04d", as.numeric(df_A$ID))
df_B <- C1C2C3_OPT_SUMMARY

# Assuming df_A and df_B are both loaded and have columns: ID, OPT_LB, OPT_UB

# Rename columns for clarity before merging
df_A_comp <- df_A %>%
  select(ID, OPT_LB, OPT_UB) %>%
  rename(OPT_LB_A = OPT_LB, OPT_UB_A = OPT_UB)

df_B_comp <- df_B %>%
  select(ID, OPT_LB, OPT_UB) %>%
  rename(OPT_LB_B = OPT_LB, OPT_UB_B = OPT_UB)

# Merge by ID
comparison_df <- df_A_comp %>%
  inner_join(df_B_comp, by = "ID")

#### PLOT COMPARISON BETWEEN VERSIONS

# Load data
# (Assume df_annotated, df_A, df_B are already loaded or load them here)
# df_annotated: long data with all Hz rows
# df_A: summary from Version A (columns: ID, OPT_LB, OPT_UB)
# df_B: summary from Version B (columns: ID, OPT_LB, OPT_UB)


# Create output folder if not present
if (!dir.exists("OPT_FRQ_COMPARISON_PLOTS")) {
  dir.create("OPT_FRQ_COMPARISON_PLOTS")
}

for (id in comparison_df$ID) {
  df_sub <- df_annotated %>% filter(ID == id)
  # Significant points
  sig_cw  <- df_sub %>% filter(CW_Shz == 1)
  sig_var <- df_sub %>% filter(PH_VAR_Shz == 1)
  sig_zc  <- df_sub %>% filter(PH_ZC_Shz == 1)
  # Version A band
  band_A <- comparison_df %>% filter(ID == id) %>% select(OPT_LB_A, OPT_UB_A)
  band_rect_A <- if (!is.na(band_A$OPT_LB_A) && !is.na(band_A$OPT_UB_A)) {
    data.frame(xmin = band_A$OPT_LB_A, xmax = band_A$OPT_UB_A, fill_label = "Version A")
  } else {
    NULL
  }
  # Version B band
  band_B <- comparison_df %>% filter(ID == id) %>% select(OPT_LB_B, OPT_UB_B)
  band_rect_B <- if (!is.na(band_B$OPT_LB_B) && !is.na(band_B$OPT_UB_B)) {
    data.frame(xmin = band_B$OPT_LB_B, xmax = band_B$OPT_UB_B, fill_label = "Version B")
  } else {
    NULL
  }
  # Plot
  # Combine for legend
  band_rects <- rbind(band_rect_A, band_rect_B)
  p <- ggplot(df_sub, aes(x = Hz)) +
    # Overlay both bands with fill legend
    {if (!is.null(band_rects)) geom_rect(
      data = band_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_label),
      inherit.aes = FALSE,
      alpha = 0.3
    )} +
    # Metric curves
    geom_line(aes(y = CW_POWER_GAMZ_SM, color = "CW_POWER"), linewidth = 1) +
    geom_line(aes(y = PHS_VAR_GAMZ_SM, color = "PHASE_VARIANCE"), linewidth = 1) +
    geom_line(aes(y = PHS_ZERO_CROSS_GAMZ_SM, color = "PHASE_ZERO_CROSS"), linewidth = 1) +
    # Significant points
    geom_point(data = sig_cw, aes(x = Hz, y = CW_POWER_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_var, aes(x = Hz, y = PHS_VAR_GAMZ_SM), color = "red", size = 1.5) +
    geom_point(data = sig_zc, aes(x = Hz, y = PHS_ZERO_CROSS_GAMZ_SM), color = "red", size = 1.5) +
    scale_color_manual(values = c(
      "CW_POWER" = "blue",
      "PHASE_VARIANCE" = "purple",
      "PHASE_ZERO_CROSS" = "green4"
    )) +
    scale_fill_manual(
      name = "Optimized Band",
      values = c("Version A" = "skyblue", "Version B" = "orange"),
      guide = guide_legend(override.aes = list(alpha = 0.5))
    ) +
    labs(
      title = paste("Optimized Bands Comparison — Dyad", id),
      x = "Hz", y = "Value", color = "Metric"
    ) +
    theme_minimal(base_size = 14)
  ggsave(
    filename = paste0("OPT_FRQ_COMPARISON_PLOTS/", id, "_OPT_BAND_COMPARISON.pdf"),
    plot = p, width = 8, height = 5
  )
} 
