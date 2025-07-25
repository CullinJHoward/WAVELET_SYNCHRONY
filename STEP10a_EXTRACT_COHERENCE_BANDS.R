################################################################################
######################### LIBRARY SHTUFFF
library(dplyr)
library(DescTools)
library(ggplot2)
library(tidyr)

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
##################### LOAD IN THE COHERENCE DFs


# LOCATE THE POWER MEAN DATA FROM WITHIN THE COI
folder <- file.path("COHERENCE_FILTERED\\")

# MAKE A LIST OF THE DFs
csv_files <- list.files(path = folder, pattern = "^COHERENCE_\\d{4}_.*\\.csv$", full.names = TRUE)

# FUNCTION TO LOAD, RENAME, AND CLEAN UP THE DFs
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the family number by extracting the first four digits of the file name
  family_number <- sub("^COHERENCE_(\\d{4})_.*\\.csv$", "\\1", base_name)
  
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
#####################  LOAD IN SUMMARY DATA AND CREATE BINNED MODIFICATIONS 

## LOAD IN OPTIMIZED SUMMARY DFs 

#POWER_SUMMARY_DF <- read.csv("OPTIMIZED_COHERENCE_SUMMARY_V2.csv")

OPT_SUM_A <- read.csv("FRQ_OPT_CLASSIFICATION_VERSION_A_SUMMARY.csv")
OPT_SUM_B <- read.csv("FRQ_OPT_CLASSIFICATION_VERSION_B_SUMMARY.csv")
OPT_SUM_A$ID <- sprintf("%04d", as.numeric(OPT_SUM_A$ID))
OPT_SUM_B$ID <- sprintf("%04d", as.numeric(OPT_SUM_B$ID))


################################################################################
############## TRIM DATA TO OPTIMIZED COHERENCE BANDS PREVIOUSLY IDENTIFIED

#VERSION A

# Get all data frames matching "D_####" format
coherence_dfs <- ls(pattern = "^D_\\d{4}$")


## LOOP THROUGH AND REDUCE BASED ON OPTIMIZED COHERENCE WINDOWS 

for (df_name in coherence_dfs) {
  
  dyad_id <- sub("^D_", "", df_name)
  summary_row <- OPT_SUM_A[OPT_SUM_A$ID == dyad_id, ]
  
  if (nrow(summary_row) == 0) {
    cat("Skipping", dyad_id, "- No matching row found.\n")
    next
  }
  
  oc_lb <- summary_row$OPT_LB
  oc_ub <- summary_row$OPT_UB
  
  df <- get(df_name)
  
  # Convert Period to Hz, remove Period, and move Hz to first column
  df$Hz <- round(as.numeric(1/df$Period), 6)
  df$Period <- NULL
  df <- df[, c("Hz", setdiff(names(df), "Hz"))]
  
  # Convert logical columns to numeric while preserving NAs
  df[, sapply(df, is.logical)] <- lapply(df[, sapply(df, is.logical)], as.numeric)
  
  # Convert all remaining non-numeric columns (except Hz) to numeric
  df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))] <- 
    lapply(df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))], as.numeric)
  
  # Debugging: Print Hz range before filtering
  cat("Dyad:", dyad_id, "- Hz range in DF:", min(df$Hz, na.rm = TRUE), "to", max(df$Hz, na.rm = TRUE), "\n")
  cat("Dyad:", dyad_id, "- OC_LB:", oc_lb, "- OC_UB:", oc_ub, "\n")
  
  # Apply filtering
  reduced_df <- df[df$Hz >= oc_lb & df$Hz <= oc_ub, ]
  
  # Debugging: Check if any rows were selected
  cat("Dyad:", dyad_id, "- Rows before filtering:", nrow(df), "- Rows after filtering:", nrow(reduced_df), "\n")
  
  # Assign the reduced data frame to a new object in the environment
  assign(paste0("OPT_A_", dyad_id), reduced_df, envir = .GlobalEnv)
}

## CHECK HOW MANY IT WAS APPLIED TO 

length(ls(pattern = "^OPT_A_"))

#VERSION B


## LOOP THROUGH AND REDUCE BASED ON OPTIMIZED COHERENCE WINDOWS 

for (df_name in coherence_dfs) {
  
  dyad_id <- sub("^D_", "", df_name)
  summary_row <- OPT_SUM_B[OPT_SUM_B$ID == dyad_id, ]
  
  if (nrow(summary_row) == 0) {
    cat("Skipping", dyad_id, "- No matching row found.\n")
    next
  }
  
  oc_lb <- summary_row$OPT_LB
  oc_ub <- summary_row$OPT_UB
  
  df <- get(df_name)
  
  # Convert Period to Hz, remove Period, and move Hz to first column
  df$Hz <- round(as.numeric(1/df$Period), 6)
  df$Period <- NULL
  df <- df[, c("Hz", setdiff(names(df), "Hz"))]
  
  # Convert logical columns to numeric while preserving NAs
  df[, sapply(df, is.logical)] <- lapply(df[, sapply(df, is.logical)], as.numeric)
  
  # Convert all remaining non-numeric columns (except Hz) to numeric
  df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))] <- 
    lapply(df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))], as.numeric)
  
  # Debugging: Print Hz range before filtering
  cat("Dyad:", dyad_id, "- Hz range in DF:", min(df$Hz, na.rm = TRUE), "to", max(df$Hz, na.rm = TRUE), "\n")
  cat("Dyad:", dyad_id, "- OC_LB:", oc_lb, "- OC_UB:", oc_ub, "\n")
  
  # Apply filtering
  reduced_df <- df[df$Hz >= oc_lb & df$Hz <= oc_ub, ]
  
  # Debugging: Check if any rows were selected
  cat("Dyad:", dyad_id, "- Rows before filtering:", nrow(df), "- Rows after filtering:", nrow(reduced_df), "\n")
  
  # Assign the reduced data frame to a new object in the environment
  assign(paste0("OPT_B_", dyad_id), reduced_df, envir = .GlobalEnv)
}

## CHECK HOW MANY IT WAS APPLIED TO 

length(ls(pattern = "^OPT_B_"))

################################################################################
###################### EXTRACT TIME-VARYING OPTIMIZED COHERENCE VALUES 

##################  DOWNSAMPLE TO 1 SECOND INTERVALS 

# CREATE A FUNCTION TO REDUCE THE DF TO 1 SECOND OBSERVATIONS (DOWNSAMPLE FROM 20HZ TO 1 HZ)

# Function to downsample columns by averaging every 20 columns into 1 (20 Hz to 1 Hz), skipping Period column
reduce_columns <- function(df, df_name) {
  # Ensure all columns except Period are numeric
  df[, -1] <- lapply(df[, -1, drop = FALSE], as.numeric)  # Skip column 1
  
  # Define the number of new columns (each 20 columns collapse into 1, excluding Period)
  time_series_cols <- ncol(df) - 1  # Exclude Period column
  num_groups <- floor(time_series_cols / 20)  # Number of 20-column blocks for 1 Hz
  
  # If there are no groups (less than 20 time series columns), return with warning
  if (num_groups == 0) {
    warning(paste("Data frame", df_name, "has fewer than 20 time series columns. No downsampling applied."))
    new_df_name <- paste0("DS_", df_name)
    assign(new_df_name, df, envir = .GlobalEnv)
    return(new_df_name)
  }
  
  # Compute the mean for every 20-column block, starting from column 2
  df_reduced <- as.data.frame(do.call(cbind, lapply(seq_len(num_groups), function(i) {
    start_col <- (i - 1) * 20 + 2  # Start at column 2 (after Period)
    end_col <- start_col + 19      # Fixed 20-column block for 1 Hz
    rowMeans(df[, start_col:end_col, drop = FALSE], na.rm = TRUE)  # Compute mean per row
  })))
  
  # Rename columns to OBS_1, OBS_2, etc.
  colnames(df_reduced) <- paste0("OBS_", seq_len(num_groups))
  
  # Add the Hz column back to the reduced data frame
  df_reduced <- cbind(Hz = df$Hz, df_reduced)
  
  # Rename the data frame with a DS_ prefix
  new_df_name <- paste0("DS_", df_name)
  assign(new_df_name, df_reduced, envir = .GlobalEnv)
  
  return(new_df_name)
}

# Identify all data frames with the OPT_ prefix and apply downsampling
opt_dfs <- ls(pattern = "^OPT_", envir = .GlobalEnv)
if (length(opt_dfs) == 0) {
  stop("No data frames with prefix 'OPT_' found in the environment.")
}

for (df_name in opt_dfs) {
  df <- get(df_name, envir = .GlobalEnv)  # Retrieve data frame
  reduce_columns(df, df_name)  # Apply downsampling function
}

# REMOVE dfs FROM THE ENVIRONMENT TO SAVE SPACE 

# LIST ALL OBJECTS WITH PREFIX OPT_A_ or OPT_B_
objs_to_remove <- ls(pattern = "^OPT_[AB]_")

# REMOVE THEM 
rm(list = objs_to_remove, envir = .GlobalEnv)


################################################################################
################## REMOVE NaNs so that the coherence pattern is complete 

## Here we are dropping NaN rows because we will do a within-dyad nonlinear
## mixed effects model. The NaNs will mess up the estiamtion and having the same
## number of columns is not necessary because at the group level we will just
## use the parameter estimates. 

# Identify all data frames with the DS_OPT_ prefix
ds_opt_dfs <- ls(pattern = "^DS_OPT_", envir = .GlobalEnv)
if (length(ds_opt_dfs) == 0) {
  stop("No data frames with prefix 'DS_OPT_' found in the environment.")
}

# Loop through each DS_OPT_ data frame and remove columns with any NA
for (df_name in ds_opt_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Identify columns with no NA values
  complete_cols <- sapply(df, function(col) !any(is.na(col)))
  
  # Subset the data frame to keep only complete columns
  df_cleaned <- df[, complete_cols, drop = FALSE]
  
  # Check if any columns remain; if not, issue a warning
  if (ncol(df_cleaned) == 0) {
    warning(paste("Data frame", df_name, "has no columns without NA. Keeping original structure."))
    assign(df_name, df, envir = .GlobalEnv)
  } else {
    # Overwrite the original data frame with the cleaned version
    assign(df_name, df_cleaned, envir = .GlobalEnv)
  }
  
  # Optional: Print a message to track progress
  cat(sprintf("Processed %s: %d columns retained out of %d\n", 
              df_name, ncol(df_cleaned), ncol(df)))
}


################################################################################
################## MAKE COLUMNS IN LONG FORMAT ADD EXTRA VARIABLES 


# Loop through each DS_OPT_ data frame and convert to long format
for (df_name in ds_opt_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Extract ID from the suffix (e.g., "1234" from "DS_OPT_1234")
  id <- sub("^DS_OPT_[AB]_", "", df_name)
  
  # Check that df has at least 2 columns (Hz + at least 1 time point)
  if (ncol(df) < 2) {
    warning(paste("Data frame", df_name, "has fewer than 2 columns. Skipping."))
    next
  }
  
  # Convert to long format using tidyr::pivot_longer
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Hz,              # Keep Hz as the clustering column
    names_to = "OBS",            # Temporary name for time point columns
    values_to = "COH"            # Coherence values
  )
  
  # Add the ID column
  df_long$ID <- id
  
  # Group by Hz and assign TIME starting at 0
  df_long <- df_long %>%
    dplyr::group_by(Hz) %>%
    dplyr::mutate(TIME = seq(0, n() - 1)) %>%  # Changed from seq_len(n()) to start at 0
    dplyr::ungroup()
  
  # Remove the temporary OBS column
  df_long$OBS <- NULL
  
  # Reorder and select final columns: ID, Period, TIME, COH
  df_long <- df_long[, c("ID", "Hz", "TIME", "COH")]
  
  # Sort by Hz and TIME for clarity
  df_long <- df_long[order(df_long$Hz, df_long$TIME), ]
  
  # Assign the new data frame with LONG_ prefix
  new_df_name <- paste0("LONG_", df_name)
  assign(new_df_name, df_long, envir = .GlobalEnv)
  
  # Print progress
  cat(sprintf("Converted %s to %s: %d rows, %d columns\n", df_name, new_df_name, nrow(df_long), ncol(df_long)))
}


################################################################################
################# SAVE THE INDIVIDUAL TIME SERIES  

# SAVE OPTIMIZED VERSION A 

# DEFINE TARGET DIRECTORY
output_dir_A <- file.path("..\\NONLINEAR_PATTERN_MODELING\\COHERENCE\\", "FRQ_OPT_VERSION_A_LONG")
output_dir_B <- file.path("..\\NONLINEAR_PATTERN_MODELING\\COHERENCE\\", "FRQ_OPT_VERSION_B_LONG")

# CREATE DIRECTORIES IF THEY DONT EXIST
if (!dir.exists(output_dir_A)) {
  dir.create(output_dir_A, recursive = TRUE)
}
if (!dir.exists(output_dir_B)) {
  dir.create(output_dir_B, recursive = TRUE)
}

# GET ALL OBJECTS 
df_names_A <- ls(pattern = "^LONG_DS_OPT_A_")
df_names_B <- ls(pattern = "^LONG_DS_OPT_B_")

# LOOP AND SAVE ALL VERSION A
for (df_name in df_names_A) {
  
  # Get the actual data frame by name
  df <- get(df_name)
  
  # Extract the last 4 digits from the df name
  id <- sub("^LONG_DS_OPT_A_.*?(\\d{4})$", "\\1", df_name)
  
  # Construct file name
  file_name <- paste0("FRQ_OPT_A_LONG_", id, ".csv")
  
  # Full path
  full_path <- file.path(output_dir_A, file_name)
  
  # Save as CSV
  write.csv(df, full_path, row.names = FALSE)
}

# LOOP AND SAVE ALL VERSION B
for (df_name in df_names_B) {
  
  # Get the actual data frame by name
  df <- get(df_name)
  
  # Extract the last 4 digits from the df name
  id <- sub("^LONG_DS_OPT_B_.*?(\\d{4})$", "\\1", df_name)
  
  # Construct file name
  file_name <- paste0("FRQ_OPT_B_LONG_", id, ".csv")
  
  # Full path
  full_path <- file.path(output_dir_B, file_name)
  
  # Save as CSV
  write.csv(df, full_path, row.names = FALSE)
}


################################################################################
####################### EXTRACT HIGH AND LOW FREQUENCY BANDS

# REMOVE WORKING dfs FROM THE ENVIRONMENT TO SAVE SPACE 

# LIST ALL OBJECTS WITH PREFIX OPT_A_ or OPT_B_
objs_to_remove_2 <- ls(pattern = "^LONG_DS_OPT_[AB]_")
objs_to_remove_3 <- ls(pattern = "^DS_OPT_[AB]_")

# REMOVE THEM 
rm(list = objs_to_remove_2, envir = .GlobalEnv)
rm(list = objs_to_remove_3, envir = .GlobalEnv)

############## HIGH FREQUENCY 

for (df_name in coherence_dfs) {
  dyad_id <- sub("^D_", "", df_name)
  df <- get(df_name)
  
  # Convert Period to Hz, remove Period, and move Hz to first column
  df$Hz <- round(as.numeric(1/df$Period), 6)
  df$Period <- NULL
  df <- df[, c("Hz", setdiff(names(df), "Hz"))]
  
  # Convert logical columns to numeric while preserving NAs
  df[, sapply(df, is.logical)] <- lapply(df[, sapply(df, is.logical)], as.numeric)
  
  # Convert all remaining non-numeric columns (except Hz) to numeric
  df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))] <- 
    lapply(df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))], as.numeric)
  
  # Filter for high frequency band: 0.15 - 0.40 Hz
  hf_df <- df[df$Hz >= 0.15 & df$Hz <= 0.40, ]
  
  # Assign the high frequency data frame to a new object in the environment
  assign(paste0("HF_", dyad_id), hf_df, envir = .GlobalEnv)
}

# Check how many HF_ data frames were created
length(ls(pattern = "^HF_"))

############## LOW FREQUENCY 

for (df_name in coherence_dfs) {
  dyad_id <- sub("^D_", "", df_name)
  df <- get(df_name)
  
  # Convert Period to Hz, remove Period, and move Hz to first column
  df$Hz <- round(as.numeric(1/df$Period), 6)
  df$Period <- NULL
  df <- df[, c("Hz", setdiff(names(df), "Hz"))]
  
  # Convert logical columns to numeric while preserving NAs
  df[, sapply(df, is.logical)] <- lapply(df[, sapply(df, is.logical)], as.numeric)
  
  # Convert all remaining non-numeric columns (except Hz) to numeric
  df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))] <- 
    lapply(df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))], as.numeric)
  
  # Filter for low frequency band: 0.04 - 0.15 Hz
  lf_df <- df[df$Hz >= 0.04 & df$Hz <= 0.15, ]
  
  # Assign the low frequency data frame to a new object in the environment
  assign(paste0("LF_", dyad_id), lf_df, envir = .GlobalEnv)
}

# Check how many LF_ data frames were created
length(ls(pattern = "^LF_"))


################################################################################
###################### EXTRACT TIME-VARYING OPTIMIZED COHERENCE VALUES 

##################  DOWNSAMPLE TO 1 SECOND INTERVALS 


# Identify all data frames with the HF_ prefix and apply downsampling
HF_dfs <- ls(pattern = "^HF_", envir = .GlobalEnv)

for (df_name in HF_dfs) {
  df <- get(df_name, envir = .GlobalEnv)  # Retrieve data frame
  reduce_columns(df, df_name)  # Apply downsampling function
}

# Identify all data frames with the LF_ prefix and apply downsampling
LF_dfs <- ls(pattern = "^LF_", envir = .GlobalEnv)

for (df_name in LF_dfs) {
  df <- get(df_name, envir = .GlobalEnv)  # Retrieve data frame
  reduce_columns(df, df_name)  # Apply downsampling function
}

# REMOVE dfs FROM THE ENVIRONMENT TO SAVE SPACE 

# LIST ALL OBJECTS WITH PREFIX HF_ or LF_
objs_to_remove_4 <- ls(pattern = "^HF_")
objs_to_remove_5 <- ls(pattern = "^LF_")

# REMOVE THEM 
rm(list = objs_to_remove_4, envir = .GlobalEnv)
rm(list = objs_to_remove_5, envir = .GlobalEnv)

################################################################################
################## REMOVE NaNs so that the coherence pattern is complete 

## Here we are dropping NaN rows because we will do a within-dyad nonlinear
## mixed effects model. The NaNs will mess up the estiamtion and having the same
## number of columns is not necessary because at the group level we will just
## use the parameter estimates. 

# Identify all data frames with the DS_ prefix
ds_dfs <- ls(pattern = "^DS_", envir = .GlobalEnv)

# Loop through each DS_OPT_ data frame and remove columns with any NA
for (df_name in ds_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Identify columns with no NA values
  complete_cols <- sapply(df, function(col) !any(is.na(col)))
  
  # Subset the data frame to keep only complete columns
  df_cleaned <- df[, complete_cols, drop = FALSE]
  
  # Check if any columns remain; if not, issue a warning
  if (ncol(df_cleaned) == 0) {
    warning(paste("Data frame", df_name, "has no columns without NA. Keeping original structure."))
    assign(df_name, df, envir = .GlobalEnv)
  } else {
    # Overwrite the original data frame with the cleaned version
    assign(df_name, df_cleaned, envir = .GlobalEnv)
  }
  
  # Optional: Print a message to track progress
  cat(sprintf("Processed %s: %d columns retained out of %d\n", 
              df_name, ncol(df_cleaned), ncol(df)))
}


################################################################################
################## MAKE COLUMNS IN LONG FORMAT ADD EXTRA VARIABLES 


# Loop through each DS_ data frame and convert to long format
for (df_name in ds_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Extract ID from the suffix (e.g., "1234" from "DS_OPT_1234")
  id <- sub("^DS_(LF|HF)_", "", df_name)
  
  # Check that df has at least 2 columns (Hz + at least 1 time point)
  if (ncol(df) < 2) {
    warning(paste("Data frame", df_name, "has fewer than 2 columns. Skipping."))
    next
  }
  
  # Convert to long format using tidyr::pivot_longer
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Hz,              # Keep Hz as the clustering column
    names_to = "OBS",            # Temporary name for time point columns
    values_to = "COH"            # Coherence values
  )
  
  # Add the ID column
  df_long$ID <- id
  
  # Group by Hz and assign TIME starting at 0
  df_long <- df_long %>%
    dplyr::group_by(Hz) %>%
    dplyr::mutate(TIME = seq(0, n() - 1)) %>%  # Changed from seq_len(n()) to start at 0
    dplyr::ungroup()
  
  # Remove the temporary OBS column
  df_long$OBS <- NULL
  
  # Reorder and select final columns: ID, Period, TIME, COH
  df_long <- df_long[, c("ID", "Hz", "TIME", "COH")]
  
  # Sort by Hz and TIME for clarity
  df_long <- df_long[order(df_long$Hz, df_long$TIME), ]
  
  # Assign the new data frame with LONG_ prefix
  new_df_name <- paste0("LONG_", df_name)
  assign(new_df_name, df_long, envir = .GlobalEnv)
  
  # Print progress
  cat(sprintf("Converted %s to %s: %d rows, %d columns\n", df_name, new_df_name, nrow(df_long), ncol(df_long)))
}


################################################################################
################# SAVE THE INDIVIDUAL TIME SERIES  

# SAVE OPTIMIZED VERSION A 

# DEFINE TARGET DIRECTORY
output_dir_HF <- file.path("..\\NONLINEAR_PATTERN_MODELING\\COHERENCE\\", "HI_FRQ_LONG")
output_dir_LF <- file.path("..\\NONLINEAR_PATTERN_MODELING\\COHERENCE\\", "LO_FRQ_LONG")

# CREATE DIRECTORIES IF THEY DONT EXIST
if (!dir.exists(output_dir_HF)) {
  dir.create(output_dir_HF, recursive = TRUE)
}
if (!dir.exists(output_dir_LF)) {
  dir.create(output_dir_LF, recursive = TRUE)
}

# GET ALL OBJECTS 
df_names_HF <- ls(pattern = "^LONG_DS_HF_")
df_names_LF <- ls(pattern = "^LONG_DS_LF_")

# LOOP AND SAVE ALL High frequency
for (df_name in df_names_HF) {
  
  # Get the actual data frame by name
  df <- get(df_name)
  
  # Extract the last 4 digits from the df name
  id <- sub("^LONG_DS_HF_.*?(\\d{4})$", "\\1", df_name)
  
  # Construct file name
  file_name <- paste0("HF_LONG_", id, ".csv")
  
  # Full path
  full_path <- file.path(output_dir_HF, file_name)
  
  # Save as CSV
  write.csv(df, full_path, row.names = FALSE)
}

# LOOP AND SAVE ALL low frequency
for (df_name in df_names_LF) {
  
  # Get the actual data frame by name
  df <- get(df_name)
  
  # Extract the last 4 digits from the df name
  id <- sub("^LONG_DS_LF_.*?(\\d{4})$", "\\1", df_name)
  
  # Construct file name
  file_name <- paste0("LF_LONG_", id, ".csv")
  
  # Full path
  full_path <- file.path(output_dir_LF, file_name)
  
  # Save as CSV
  write.csv(df, full_path, row.names = FALSE)
}
