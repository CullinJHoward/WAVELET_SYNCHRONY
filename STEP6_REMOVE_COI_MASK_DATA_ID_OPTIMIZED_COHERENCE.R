### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)
library(WaveletComp)

# Set working directory 


POUNDTOWN <- 0

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
}

setwd(work_dir)

################################################################################
############## REMOVE CONE OF INFLUENCE DATA & EXTRACT KEY INFORMATION

# SET PATH TO RDS FILES 
coherence_dir <- "MATH_RESP_FREE_COHERENCE_DATA"

# LIST ALL RDS FILES TO PROCESS
rds_files <- list.files(coherence_dir, pattern = "\\.RDS$", full.names = TRUE)

## CREATE DIRECTORIES FOR SAVING DATA
if (!dir.exists("COHERENCE_FILTERED")) {
  dir.create("COHERENCE_FILTERED")
}
if (!dir.exists("RIDGE_FILTERED")) {
  dir.create("RIDGE_FILTERED")
}
if (!dir.exists("ANGLE_FILTERED")) {
  dir.create("ANGLE_FILTERED")
}
if (!dir.exists("CW_POWER_FILTERED")) {
  dir.create("CW_POWER_FILTERED")
}

# Define the output directories for saving filtered CSVs
coherence_output_dir <- "COHERENCE_FILTERED"
ridge_output_dir <- "RIDGE_FILTERED"
angle_output_dir <- "ANGLE_FILTERED"
cw_power_output_dir <- "CW_POWER_FILTERED"  # New output directory for CW power

for (rds_file in rds_files) {
  # Load the .RDS file
  data <- readRDS(rds_file)
  
  # Extract the family ID (last 4 digits of the file name)
  family_id <- sub(".*(\\d{4})\\.RDS$", "\\1", basename(rds_file))
  
  # Pull out coi.1 and coi.2
  coi_1 <- data$coi.1
  COI_THRESH <- data$coi.2
  
  # Remove the first two and last two values from both coi.1 and coi.2 
  # (they have 6004 observations and need 6000)
  coi_1 <- coi_1[3:(length(coi_1)-2)]
  COI_THRESH <- COI_THRESH[3:(length(COI_THRESH)-2)]
  
  # Scale coi.2 by 20 to accurately reflect the sampling frequency (20Hz)
  COI_THRESH <- COI_THRESH * 20
  
  # INITIALIZE A MASK OBJECT
  n_rows <- 119
  n_cols <- 6000
  
  # Create a mask for each matrix
  valid_mask <- matrix(NA, nrow = n_rows, ncol = n_cols)
  
  # Loop over each column (observation) and apply the thresholding based on coi.2 values
  for (i in 1:n_cols) {
    # For the ith observation (column), create a mask
    threshold <- COI_THRESH[i]  # The threshold value for this column
    valid_mask[, i] <- ifelse(1:n_rows <= threshold, 1, NA)  # If row index <= threshold, valid (1); otherwise NA
  }
  
  # Filter the matrices based on the mask
  coherence <- data$Coherence
  ridge_xy <- data$Ridge.xy
  angle <- data$Angle
  power_xy <- data$Power.xy  # Pulling in the Power.xy variable
  
  # Apply the mask to filter the data for Coherence, Ridge.xy, Angle, and Power.xy
  coherence_filtered <- coherence * valid_mask
  ridge_filtered <- ridge_xy * valid_mask
  angle_filtered <- angle * valid_mask
  power_filtered <- power_xy * valid_mask  # Apply mask to Power.xy
  
  # Add the 'Period' column to each of the filtered data frames
  period <- data$Period  
  
  coherence_filtered <- cbind(Period = period, coherence_filtered)
  ridge_filtered <- cbind(Period = period, ridge_filtered)
  angle_filtered <- cbind(Period = period, angle_filtered)
  power_filtered <- cbind(Period = period, power_filtered)  # Add Power.xy to filtered data
  
  # Create the column names dynamically using coi.1 values
  period_columns <- paste0("T_", format(coi_1, nsmall = 2))  # Format coi.1 with 2 decimal places (e.g., T_1, T_1.05, ...)
  
  # Assign the new column names (period + the 6000 time points)
  colnames(coherence_filtered)[2:ncol(coherence_filtered)] <- period_columns
  colnames(ridge_filtered)[2:ncol(ridge_filtered)] <- period_columns
  colnames(angle_filtered)[2:ncol(angle_filtered)] <- period_columns
  colnames(power_filtered)[2:ncol(power_filtered)] <- period_columns  # Assign column names for Power.xy
  
  # Save each filtered data frame as a CSV file with appropriate naming and paths
  write.csv(coherence_filtered, file = file.path(coherence_output_dir, paste0("COHERENCE_", family_id, ".csv")), row.names = FALSE)
  write.csv(ridge_filtered, file = file.path(ridge_output_dir, paste0("RIDGE_", family_id, ".csv")), row.names = FALSE)
  write.csv(angle_filtered, file = file.path(angle_output_dir, paste0("ANGLE_", family_id, ".csv")), row.names = FALSE)
  write.csv(power_filtered, file = file.path(cw_power_output_dir, paste0("CW_POWER_", family_id, ".csv")), row.names = FALSE)  # Save Power.xy data
  
  # Optional: You can print a message for each file saved
  message(paste("Saved filtered data for Family ID:", family_id))
}


################################################################################
############# IDENTIFY OPTIMIZED COHERENCE 

# After removing the cone of influence data, we will identify the new regions
# that have the highest power of coherence. I did this by saving the power.xy data 
# this is the "cross-wavelet power (analogous to Fourier cross-frequency power spectrum)"
# for each cell. I then summed them across rows (leaving NAs where the COI was)
# to get updated Power Means. 

# Optimized coherence bands are identified at their peak as the highest power mean 
# frequency, with their bands being the periods that have power means that are greater
# than 75% of all their own power mean scores (visually decided). However, some data is multi-modal, or
# doesn't taper off clearly enough to identify optimal coherence bands.
# So In cases where the difference between the upper and lower
# optimized bands is greater than 20, or less than 5, a new upper and lower bound of +/- 5 was implemented.
## TOO LARGE DIFFERENCE SUGGESTS IT MAY BE MULTI-MODAL AND IT WILL NOT BE A CLEAN ESTIMATE.
## TOO SMALL IS IT IS OVER FOCUSING ON A SMALL RANGE THAT WILL HAVE EXTREMELY (AND LIKELY SPURIOUS) COHERENCE ESTIMATES.

## RANGE SELECTION: I TOOK A SUBSAMPLE OF DYADS NOT NEEDING ADJUSTMENT AND COMPUTED THE AVERAGE SIZE OF THEIR CONFIDENCE BANDS, IT WAS 
## ABOUT 10. SO 20 WOULD BE TWICE THE NORM, AND 5 WOULD BE UNDER HALF THE NORM.THEY ARE COERCED TO THE NORM (10; +/- 5)
## IF THEY VIOLATE THIS RANGE. 

## I also took a second approach of creating optimized bins. Because the VLF starts getting 
## at endocrinology difference (and not ANS), I am taking optimized scores based on if 
## they occur in the lower or higher frequency (both of which under the ANS control)
## this is kind of a happy medium because they are theoretically defined. 


#################### LOAD IN ALL SAVED POWER_MEAN DATA #########################

# LOCATE THE POWER MEAN DATA FROM WITHIN THE COI
folder <- file.path(work_dir, "CW_POWER_FILTERED")

# MAKE A LIST OF THE DFs
csv_files <- list.files(path = folder, pattern = "^CW_POWER_\\d{4}_.*\\.csv$", full.names = TRUE)

# FUNCTION TO LOAD, RENAME, AND CLEAN UP THE DFs
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the family number by extracting the first four digits of the file name
  family_number <- sub("^CW_POWER_(\\d{4})_.*\\.csv$", "\\1", base_name)
  
  # Construct the new variable name with "D_" prefix and no suffix
  new_name <- paste0("D_", family_number)  # This ensures the name is D_####, e.g., D_0002
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Do NOT remove columns with NA values (leave them as time markers)
  # No additional cleaning step needed for NA columns
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
  
  # Return the data for further manipulation if needed (optional)
  return(data)
}

# APPLY THE FUNCTION
loaded_data <- lapply(csv_files, load_and_clean_csv)

################################################################################
###### COMPUTE NEW POWER MEANS AND IDENTIFY OPTIMIZED COHERENCE BANDS 

## LIST ALL DFs that need processing
dyad_names <- ls(pattern = "^D_")

# Create a named list where names are the dyad IDs and values are the actual data frames
ALL_DYADS <- setNames(lapply(dyad_names, get), dyad_names)

################ COMPUTE OPTIMIZED COHERENCE BANDS (FULL RANGE)

# INITIALIZE EMPTY LIST FOR RESULTS SUMMARY
optimized_coherence_summary <- list()

# LOOP THROUGH DYAD DATA
for (dyad_name in names(ALL_DYADS)) {
  
  # GET CURRENT DF NAME 
  df <- ALL_DYADS[[dyad_name]]
  
  # Ensure the first column is named "Period" (for consistency)
  colnames(df)[1] <- "Period"
  
  # FILTER THE PERIOD VALUES TO BE BETWEEN 2 AND 50
  df_filtered <- df[df$Period > 2 & df$Period < 50, ]
  
  # COMPUTE POWER MEANS, SKIPPING THE PERIOD COLUMN
  df_filtered$POWER_MEAN <- rowMeans(df_filtered[, -1], na.rm = TRUE)
  
  # ID THE HIGHEST POWER MEAN AND SAVE IT
  max_power_mean_index <- which.max(df_filtered$POWER_MEAN)
  optimized_coh <- df_filtered$Period[max_power_mean_index]  
  
  # IDENTIFY 80th PERCENTILE OF POWER_MEAN SCORES
  top_80_percentile <- quantile(df_filtered$POWER_MEAN, 0.80, na.rm = TRUE)
  
  # IDENTIFY THE VALUES THAT FIT WITHIN THIS RANGE
  top_80_periods <- df_filtered$Period[df_filtered$POWER_MEAN >= top_80_percentile]
  
  # IDENTIFY THE UPPER AND LOWER BOUNDS OF THE 80 PERCENTILE
  oc_lb <- min(top_80_periods, na.rm = TRUE)
  oc_ub <- max(top_80_periods, na.rm = TRUE)
  
  # IF THE RANGE OF UPPER TO LOWER BOUND IS TOO LARGE (> 20) OR TOO SMALL (< 5), MANUALLY ADJUST BANDS. 
  multi_modal <- 0  # IF NOT REQUIRING ADJUSTMENT, THEN IT IS NOT MULTI MODAL AND IS NOTED AS 0
  if ((oc_ub - oc_lb) > 20 || (oc_ub - oc_lb) < 5) {
    # ADJUST TO +/- 7.5 AROUND THE OPTIMIZED COHERENCE LEVEL
    oc_lb <- optimized_coh - 5
    oc_ub <- optimized_coh + 5
    multi_modal <- 1  # FOR DYADS NEEDING ADJUSTMENT, MARKED AS A 1
  }
  
  # ENFORCE UPPER AND LOWER BOUNDS STAY IN REASONABLE RANGES (NO LOWER THAN 2, NO HIGHER THAN 59)
  oc_lb <- max(oc_lb, 2)
  oc_ub <- min(oc_ub, 59)
  
  # DROP "D_" FROM IDs
  dyad_id <- sub("^D_", "", dyad_name)  
  
  # PUT TOGETHER A SUMMARY DATA FRAME
  summary <- data.frame(
    ID = dyad_id,
    OPT_PERIOD = optimized_coh,
    OPT_LB = oc_lb,
    OPT_UB = oc_ub,
    MULTI_MODAL = multi_modal
  )
  
  # ADD THE SUMMARY TO THE LIST
  optimized_coherence_summary[[dyad_name]] <- summary
}

# COMBINE SUMMARIES TOGETHER
OPTIMIZED_COHERENCE_SUMMARY <- do.call(rbind, optimized_coherence_summary)


################ COMPUTE OPTIMIZED COHERENCE BINS (REDUCED TO HF or LF BAND)

# INITIALIZE EMPTY LIST FOR RESULTS SUMMARY
optimized_BINS_coherence_summary <- list()

# LOOP THROUGH DYAD DATA
for (dyad_name in names(ALL_DYADS)) {
  
  # GET CURRENT DF NAME 
  df <- ALL_DYADS[[dyad_name]]
  
  # Ensure the first column is named "Period" (for consistency)
  colnames(df)[1] <- "Period"
  
  # FILTER THE PERIOD VALUES TO BE BETWEEN 2 AND 25 (this is .5 - .04 Hz [only HF and LF])
  df_filtered <- df[df$Period > 2 & df$Period < 25, ]
  
  # COMPUTE POWER MEANS, SKIPPING THE PERIOD COLUMN
  df_filtered$POWER_MEAN <- rowMeans(df_filtered[, -1], na.rm = TRUE)
  
  # ID THE HIGHEST POWER MEAN AND SAVE IT
  max_power_mean_index <- which.max(df_filtered$POWER_MEAN)
  optimized_coh <- df_filtered$Period[max_power_mean_index]  
  
  # IDENTIFY 80th PERCENTILE OF POWER_MEAN SCORES
  top_80_percentile <- quantile(df_filtered$POWER_MEAN, 0.80, na.rm = TRUE)
  
  # IDENTIFY THE VALUES THAT FIT WITHIN THIS RANGE
  top_80_periods <- df_filtered$Period[df_filtered$POWER_MEAN >= top_80_percentile]
  
  # IDENTIFY THE UPPER AND LOWER BOUNDS OF THE 80 PERCENTILE
  oc_lb <- min(top_80_periods, na.rm = TRUE)
  oc_ub <- max(top_80_periods, na.rm = TRUE)
  
  # IF THE RANGE OF UPPER TO LOWER BOUND IS TOO LARGE (> 20) OR TOO SMALL (< 5), MANUALLY ADJUST BANDS. 
  multi_modal <- 0  # IF NOT REQUIRING ADJUSTMENT, THEN IT IS NOT MULTI MODAL AND IS NOTED AS 0
  if ((oc_ub - oc_lb) > 20 || (oc_ub - oc_lb) < 5) {
    # ADJUST TO +/- 7.5 AROUND THE OPTIMIZED COHERENCE LEVEL
    oc_lb <- optimized_coh - 5
    oc_ub <- optimized_coh + 5
    multi_modal <- 1  # FOR DYADS NEEDING ADJUSTMENT, MARKED AS A 1
  }
  
  # ENFORCE UPPER AND LOWER BOUNDS STAY IN REASONABLE RANGES (NO LOWER THAN 2, NO HIGHER THAN 59)
  oc_lb <- max(oc_lb, 2)
  oc_ub <- min(oc_ub, 59)
  
  # DROP "D_" FROM IDs
  dyad_id <- sub("^D_", "", dyad_name)  
  
  # PUT TOGETHER A SUMMARY DATA FRAME
  summary <- data.frame(
    ID = dyad_id,
    OPT_PERIOD = optimized_coh,
    OPT_LB = oc_lb,
    OPT_UB = oc_ub,
    MULTI_MODAL = multi_modal
  )
  
  # ADD THE SUMMARY TO THE LIST
  optimized_BINS_coherence_summary[[dyad_name]] <- summary
}

# COMBINE SUMMARIES TOGETHER
OPTIMIZED_BINS_COHERENCE_SUMMARY <- do.call(rbind, optimized_BINS_coherence_summary)

################################################################################
########## DATA TRANSFORMATIONS

# PERIOD INTO Hz FOR EASE OF INTERPRETATION

OPTIMIZED_COHERENCE_SUMMARY$OPT_Hz <- 1/OPTIMIZED_COHERENCE_SUMMARY$OPT_PERIOD
OPTIMIZED_BINS_COHERENCE_SUMMARY$BIN_OPT_Hz <- 1/OPTIMIZED_BINS_COHERENCE_SUMMARY$OPT_PERIOD

### CATEGORICAL BINNING DESIGNATIONS FOR OPTIMIZED BINS

OPTIMIZED_BINS_COHERENCE_SUMMARY <- OPTIMIZED_BINS_COHERENCE_SUMMARY %>%
  mutate(OPT_BIN = case_when(
    BIN_OPT_Hz >= 0.04 & BIN_OPT_Hz <= 0.15 ~ "LF",
    BIN_OPT_Hz > 0.15 ~ "HF"
  ))

table(OPTIMIZED_BINS_COHERENCE_SUMMARY$OPT_BIN)

### CREATE NEW LOWER AND UPPER BOUNDARIES FOR OPTIMIZED BINS (LOCKED AT BAND BOUNDS)

OPTIMIZED_BINS_COHERENCE_SUMMARY <- OPTIMIZED_BINS_COHERENCE_SUMMARY %>%
  mutate(
    BIN_UB = case_when(
      OPT_BIN == "LF"  ~ 6.67, # 0.15Hz
      OPT_BIN == "HF"  ~ 2.5 # 0.40Hz 
    ),
    BIN_LB = case_when(
      OPT_BIN == "LF"  ~ 25, # 0.04Hz
      OPT_BIN == "HF"  ~ 6.67 # 0.15Hz
    )
  )

## ADD LEADING ZEROS TO ID SO IT MATCHES DF NAME IDS

OPTIMIZED_COHERENCE_SUMMARY$ID <- sprintf("%04d", as.numeric(OPTIMIZED_COHERENCE_SUMMARY$ID))
OPTIMIZED_BINS_COHERENCE_SUMMARY$ID <- sprintf("%04d", as.numeric(OPTIMIZED_BINS_COHERENCE_SUMMARY$ID))


## REDUCE VARIABLES AND MERGE 

BIN_DATA <- OPTIMIZED_BINS_COHERENCE_SUMMARY %>%
  select(c(ID, BIN_OPT_Hz, OPT_BIN, BIN_UB, BIN_LB))

JOINED <- full_join(OPTIMIZED_COHERENCE_SUMMARY, BIN_DATA, by = "ID")

## REORDER VARIABLES 

names(JOINED)

JOINED <- 
  JOINED[, c("ID", "OPT_Hz", "OPT_PERIOD", 
         "OPT_LB", "OPT_UB", "MULTI_MODAL", 
         "OPT_BIN", "BIN_OPT_Hz", "BIN_UB",  "BIN_LB")]

names(JOINED)

################################################################################
###################### SAVE THE OPTIMIZED VALUE SUMMARY SCORE DF FOR LATER USAGE 

# Define the file path where you want to save the CSV
OUTPUT_SUMMARY <- "F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\OPTIMIZED_COHERENCE_SUMMARY.csv"

# Save the data frame as a CSV file
write.csv(JOINED, file = OUTPUT_SUMMARY, row.names = FALSE)

