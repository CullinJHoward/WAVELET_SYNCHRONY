### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)
library(tidyr)
library(circular)

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
####################### EXTRACTING PHASE ESTIMATES 


# LOCATE THE POWER MEAN DATA FROM WITHIN THE COI
folder <- file.path("F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\ANGLE_FILTERED\\")

# MAKE A LIST OF THE DFs
csv_files <- list.files(path = folder, pattern = "^ANGLE_\\d{4}_.*\\.csv$", full.names = TRUE)

# FUNCTION TO LOAD, RENAME, AND CLEAN UP THE DFs
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the family number by extracting the first four digits of the file name
  family_number <- sub("^ANGLE_(\\d{4})_.*\\.csv$", "\\1", base_name)
  
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
############## TRIM DATA TO LOW FREQUENCY BAND

# Create an empty list to store reduced coherence data frames
REDUCED_COHERENCE_DFS <- list()

# Get all data frames matching "D_####" format
coherence_dfs <- ls(pattern = "^D_\\d{4}$")

## LOOP THROUGH AND REDUCE BASED ON FIXED PERIOD WINDOW

for (df_name in coherence_dfs) {
  
  dyad_id <- sub("^D_", "", df_name)
  
  df <- get(df_name)
  
  # Ensure Period is numeric
  df$Period <- as.numeric(df$Period)
  
  # Convert logical columns to numeric while preserving NAs
  df[, sapply(df, is.logical)] <- lapply(df[, sapply(df, is.logical)], as.numeric)
  
  # Convert all remaining non-numeric columns (except Period) to numeric
  df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))] <- 
    lapply(df[, sapply(df, function(col) !is.numeric(col) & !is.logical(col))], as.numeric)
  
  # Debugging: Print Period range before filtering
  cat("Dyad:", dyad_id, "- Period range in DF:", min(df$Period, na.rm = TRUE), "to", max(df$Period, na.rm = TRUE), "\n")
  
  # Apply filtering with fixed bounds
  reduced_df <- df[df$Period >= 6.6 & df$Period <= 25, ]
  
  # Debugging: Check if any rows were selected
  cat("Dyad:", dyad_id, "- Rows before filtering:", nrow(df), "- Rows after filtering:", nrow(reduced_df), "\n")
  
  # Assign the reduced data frame to a new object in the environment
  assign(paste0("LF_", dyad_id), reduced_df, envir = .GlobalEnv)
}


################################################################################
######################  TIME-VARYING LOW-FREQUENCY PHASE CHANGES 


################## DOWNSAMPLE TO 1 SECOND INTERVALS (PHASE DATA)

# Function to downsample phase columns by computing circular mean every 20 columns (20 Hz to 1 Hz), skipping Period column
reduce_columns_phase <- function(df, df_name) {
  # Ensure all columns except Period are numeric
  df[, -1] <- lapply(df[, -1, drop = FALSE], as.numeric)  # Skip column 1 (Period)
  
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
  
  # Compute the circular mean for every 20-column block, starting from column 2
  df_reduced <- as.data.frame(do.call(cbind, lapply(seq_len(num_groups), function(i) {
    start_col <- (i - 1) * 20 + 2  # Start at column 2 (after Period)
    end_col <- start_col + 19      # Fixed 20-column block for 1 Hz
    block <- df[, start_col:end_col, drop = FALSE]  # Extract 20-column block
    
    # Calculate circular mean per row (assuming phase in radians)
    cos_mean <- rowMeans(cos(block), na.rm = TRUE)  # Mean of cosines
    sin_mean <- rowMeans(sin(block), na.rm = TRUE)  # Mean of sines
    atan2(sin_mean, cos_mean)  # Circular mean in radians
  })))
  
  # Rename columns to OBS_1, OBS_2, etc.
  colnames(df_reduced) <- paste0("OBS_", seq_len(num_groups))
  
  # Add the Period column back to the reduced data frame
  df_reduced <- cbind(Period = df$Period, df_reduced)
  
  # Rename the data frame with a DS_ prefix
  new_df_name <- paste0("DS_", df_name)
  assign(new_df_name, df_reduced, envir = .GlobalEnv)
  
  return(new_df_name)
}

# Identify all data frames with the LF_ prefix and apply downsampling for phase
opt_dfs <- ls(pattern = "^LF_", envir = .GlobalEnv)
if (length(opt_dfs) == 0) {
  stop("No data frames with prefix 'LF_' found in the environment.")
}

for (df_name in opt_dfs) {
  df <- get(df_name, envir = .GlobalEnv)  # Retrieve data frame
  reduce_columns_phase(df, df_name)  # Apply downsampling function for phase
}


################################################################################
################## REMOVE NaNs so that the phase pattern is complete 

## Here we are dropping NaN rows because we will do a within-dyad nonlinear
## mixed effects model. The NaNs will mess up the estiamtion and having the same
## number of columns is not necessary because at the group level we will just
## use the parameter estimates. 

# Identify all data frames with the DS_OPT_ prefix
ds_lf_dfs <- ls(pattern = "^DS_LF_", envir = .GlobalEnv)
if (length(ds_lf_dfs) == 0) {
  stop("No data frames with prefix 'DS_LF_' found in the environment.")
}

# Loop through each DS_OPT_ data frame and remove columns with any NA
for (df_name in ds_lf_dfs) {
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
################## REFINE THE PERIOD TO BE TIGHTER 


# Define the target Period values
LF_target_periods <- c(10.5560632861832,  10.9283220540352, 11.3137084989848,
                    11.712685567565, 12.1257325320832,  12.553345566348, 
                    12.9960383416998, 13.4543426440594, 13.928809012738,
                    14.4200074017733, 14.9285278645889, 15.4549812627975)

# Get all data frames with prefix "DS_HF_"
ds_lf_dfs <- ls(pattern = "^DS_LF_")

# Loop through each data frame
for (df_name in ds_lf_dfs) {
  # Get the data frame
  df <- get(df_name)
  
  # Ensure Period is numeric
  df$Period <- as.numeric(df$Period)
  
  # Filter rows where Period is in LF_target_periods (using %in% for exact matching)
  filtered_df <- df[df$Period %in% LF_target_periods, ]
  
  # Debugging: Print number of rows before and after filtering
  cat("Data frame:", df_name, "- Rows before:", nrow(df), "- Rows after:", nrow(filtered_df), "\n")
  
  # Assign the filtered data frame back to the environment with a new prefix (e.g., "FILTERED_DS_HF_")
  assign(paste0("T_", df_name), filtered_df, envir = .GlobalEnv)
}


################################################################################
################## MAKE COLUMNS IN LONG FORMAT ADD EXTRA VARIABLES 


# Identify all data frames with the T_DS_LF_ prefix
ds_lf_dfs <- ls(pattern = "^T_DS_LF_", envir = .GlobalEnv)
if (length(ds_lf_dfs) == 0) {
  stop("No data frames with prefix 'T_DS_LF_' found in the environment.")
}

# Loop through each T_DS_LF_ data frame and convert to long format
for (df_name in ds_lf_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Extract ID from the suffix (e.g., "1234" from "DS_OPT_1234")
  id <- sub("^T_DS_LF_", "", df_name)
  
  # Check that df has at least 2 columns (Period + at least 1 time point)
  if (ncol(df) < 2) {
    warning(paste("Data frame", df_name, "has fewer than 2 columns. Skipping."))
    next
  }
  
  # Convert to long format using tidyr::pivot_longer
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Period,              # Keep Period as the clustering column
    names_to = "OBS",            # Temporary name for time point columns
    values_to = "PHASE"          # PHASE values
  )
  
  # Add the ID column
  df_long$ID <- id
  

  
  # Group by Period and assign TIME starting at 0
  df_long <- df_long %>%
    dplyr::group_by(Period) %>%
    dplyr::mutate(TIME = seq(0, n() - 1)) %>%  # Changed from seq_len(n()) to start at 0
    dplyr::ungroup()
  
  # Remove the temporary OBS column
  df_long$OBS <- NULL
  
  # Reorder and select final columns: ID, Period, TIME, PHASE
  df_long <- df_long[, c("ID", "Period", "TIME", "PHASE")]
  
  # Sort by Period and TIME for clarity
  df_long <- df_long[order(df_long$Period, df_long$TIME), ]
  
  # Assign the new data frame with LONG_ prefix
  new_df_name <- paste0("LONG_", df_name)
  assign(new_df_name, df_long, envir = .GlobalEnv)
  
  # Print progress
  cat(sprintf("Converted %s to %s: %d rows, %d columns", 
              df_name, new_df_name, nrow(df_long), ncol(df_long)))
}


################################################################################
################# SAVE INDIVIDUAL PHASE SERIES 

# Define the folder name
folder_name <- "LONG_LOW_FREQ_PHASE_SERIES"

# Create the folder if it doesn’t exist
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  cat(sprintf("Created folder: %s\n", folder_name))
} else {
  cat(sprintf("Folder %s already exists\n", folder_name))
}

# Identify all data frames with the LONG_DS_OPT_ prefix
long_ds_lf_dfs <- ls(pattern = "^LONG_T_DS_LF_", envir = .GlobalEnv)
if (length(long_ds_lf_dfs) == 0) {
  stop("No data frames with prefix 'LONG_T_DS_LF_' found in the environment.")
}

# Loop through each LONG_T_DS_LF_ data frame and save as CSV
for (df_name in long_ds_lf_dfs) {
  # Retrieve the data frame
  df <- get(df_name, envir = .GlobalEnv)
  
  # Extract the 4-digit ID from the data frame name
  id <- sub(".*LONG_T_DS_LF_(\\d{4}).*", "\\1", df_name)
  
  # Define the CSV file path using the extracted ID
  csv_file <- file.path(folder_name, paste0("LONG_T_DS_LF_", id, ".csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file = csv_file, row.names = FALSE)
  
  # Print progress
  cat(sprintf("Saved %s as %s\n", df_name, csv_file))
}



