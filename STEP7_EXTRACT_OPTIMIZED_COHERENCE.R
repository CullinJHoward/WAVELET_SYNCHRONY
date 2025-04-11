### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)
library(WaveletComp)

# Set working directory 


POUNDTOWN <- 0

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
}

setwd(work_dir)



################################################################################
##################### EXTRACTING COHERENCE ESTIMATES ###########################
################################################################################


### LOAD IN ALL COHERENCE DFS AND THE SUMMARY SCORE DF 

# LOCATE THE POWER MEAN DATA FROM WITHIN THE COI
folder <- file.path(work_dir, "COHERENCE_FILTERED")

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


## LOAD IN OPTIMIZED SUMMARY DF 

POWER_SUMMARY_DF <- read.csv("OPTIMIZED_COHERENCE_SUMMARY.csv")

## ADD LEADING ZEROS TO ID SO IT MATCHES DF NAME IDS

POWER_SUMMARY_DF$ID <- sprintf("%04d", as.numeric(POWER_SUMMARY_DF$ID))



################################################################################
############## TRIM DATA TO OPTIMIZED COHERENCE BANDS PREVIOUSLY IDENTIFIED


# Create an empty list to store reduced coherence data frames
REDUCED_COHERENCE_DFS <- list()

# Get all data frames matching "D_####" format
coherence_dfs <- ls(pattern = "^D_\\d{4}$")


## LOOP THROUGH AND REDUCE BASED ON OPTIMIZED COHERENCE WINDOWS 

for (df_name in coherence_dfs) {
  
  dyad_id <- sub("^D_", "", df_name)
  summary_row <- POWER_SUMMARY_DF[POWER_SUMMARY_DF$ID == dyad_id, ]
  
  if (nrow(summary_row) == 0) {
    cat("Skipping", dyad_id, "- No matching row found.\n")
    next
  }
  
  oc_lb <- summary_row$OPT_LB
  oc_ub <- summary_row$OPT_UB
  
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
  cat("Dyad:", dyad_id, "- OC_LB:", oc_lb, "- OC_UB:", oc_ub, "\n")
  
  # Apply filtering
  reduced_df <- df[df$Period >= oc_lb & df$Period <= oc_ub, ]
  
  # Debugging: Check if any rows were selected
  cat("Dyad:", dyad_id, "- Rows before filtering:", nrow(df), "- Rows after filtering:", nrow(reduced_df), "\n")
  
  # Assign the reduced data frame to a new object in the environment
  assign(paste0("OPT_", dyad_id), reduced_df, envir = .GlobalEnv)
}

## CHECK HOW MANY IT WAS APPLIED TO 

length(ls(pattern = "^OPT_"))


################################################################################
###################### APPROACH 1 - TIME-INVARIANT OPTIMIZED MEAN COHERNECE

# Function to Winsorize values beyond 3 SD
winsorize <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  lower_bound <- max(0, mu - 3 * sigma)  # Ensure coherence values stay within 0-1
  upper_bound <- min(1, mu + 3 * sigma)
  
  # Only apply Winsorization to non-NA values
  x_non_na <- x[!is.na(x)]
  x[!is.na(x)] <- pmin(pmax(x_non_na, lower_bound), upper_bound)
  
  return(x)
}

# Initialize summary dataframe
MEAN_OPT_COH_SUMMARY <- data.frame(ID = character(), W_OPT_MEAN_COH = numeric(), stringsAsFactors = FALSE)

# Loop through all dataframes with prefix "HF_"
for (df_name in ls(pattern = "^OPT_")) {
  df <- get(df_name)  # Retrieve the dataframe
  
  # Extract ID (remove "HF_" prefix)
  ID <- gsub("^OPT_", "", df_name)
  
  # Apply Winsorization to coherence values (excluding the first column, "Period")
  df[ , -1] <- apply(df[ , -1], 2, winsorize)
  
  # Compute row-wise means
  row_means <- rowMeans(df[ , -1], na.rm = TRUE)
  
  # Compute final mean across all rows
  final_mean <- mean(row_means, na.rm = TRUE)
  
  # Store results
  MEAN_OPT_COH_SUMMARY <- rbind(MEAN_OPT_COH_SUMMARY, data.frame(ID = ID, W_OPT_MEAN_COH = final_mean))
}


hist(MEAN_OPT_COH_SUMMARY$W_OPT_MEAN_COH)

### SAVE TIME-INVARIANT OPTIMIZED MEAN SCORES 

# Define the file path where you want to save the CSV
OUTPUT_TINV_HF_COH_MEANS <- file.path(work_dir, "TIME_INV_OPT_COH_SCORES.csv")

# Save the data frame as a CSV file
write.csv(MEAN_OPT_COH_SUMMARY, file = OUTPUT_TINV_HF_COH_MEANS, row.names = FALSE)


################################################################################
###################### APPROACH 2 TIME-VARYING OPTIMIZED COHERENCE VALUES 


################## WINSORIZE AND REDUCE THE DATA TO A SINGLE ROW


# CREATE WINSORIZE (row-wise) FUNCTION TO HANDLE EXTREME OUTLIERS 
winsorize_row <- function(x) {
  # Only apply to numeric values (ignore NAs)
  x_numeric <- x[!is.na(x)]  # Remove NA values for calculation
  
  if (length(x_numeric) > 0) {  # Ensure that there are numeric values to work with
    mu <- mean(x_numeric, na.rm = TRUE)
    sigma <- sd(x_numeric, na.rm = TRUE)
    lower_bound <- max(0, mu - 3 * sigma)
    upper_bound <- min(1, mu + 3 * sigma)
    
    # Winsorize the numeric values, leaving NAs untouched
    x[!is.na(x)] <- pmin(pmax(x_numeric, lower_bound), upper_bound)
  }
  
  return(x)
}

# LOOP THROUGH THE "OPT_" DFS, APPLY WINSORIZING, AND COMPUTE SIMPLE COLUMN MEANS

for (df_name in ls(pattern = "^OPT_")) {
  df <- get(df_name)  # Retrieve the dataframe
  
  # Step 1: Winsorize data row-wise (skip the first column)
  df[ , -1] <- t(apply(df[ , -1], 1, winsorize_row))  # Apply winsorization to rows (ignoring first column)
  
  # Step 2: Compute the simple column mean for each column (ignoring NAs)
  column_means <- colMeans(df[ , -1], na.rm = TRUE)  # Simple mean across rows
  
  # Step 3: Store the result in a new data frame with "RED_" prefix
  mean_df <- data.frame(t(column_means))  # Transpose to make a single row dataframe
  colnames(mean_df) <- colnames(df)[-1]  # Use column names from the original df
  
  # Save the result as a new dataframe with "RED_" prefix
  assign(paste0("RED_", df_name), mean_df)
}


##################  DOWNSAMPLE TO 2 SECOND INTERVALS 


# CREATE A FUNCTION TO REDUCE THE DF TO 2 SECOND OBSERVATIONS (DOWNSAMPLE FROM 20HZ TO .5 HZ)

reduce_columns <- function(df, df_name) {
  # Ensure all columns are numeric
  df[] <- lapply(df, as.numeric)
  
  # Define the number of new columns (each 40 columns collapse into 1)
  num_groups <- ncol(df) / 40
  
  # Compute the mean for every 40-column block while preserving row count
  df_reduced <- as.data.frame(do.call(cbind, lapply(seq_len(num_groups), function(i) {
    start_col <- (i - 1) * 40 + 1
    end_col <- min(start_col + 39, ncol(df))  # Ensure we don't exceed column limits
    rowMeans(df[, start_col:end_col], na.rm = TRUE)  # Compute mean per row
  })))
  
  # RENAME COLUMNS
  colnames(df_reduced) <- paste0("OBS_", seq_len(num_groups))
  
  # RENAME THE DF WITH A TRU_ SUFFIX
  new_df_name <- paste0("TRU_", df_name)
  assign(new_df_name, df_reduced, envir = .GlobalEnv)
  
  return(new_df_name)
}

# APPLY THE DOWNSAMPLING FUNCTION TO ALL DFS WITH THE "RED_" PREFIX
for (df_name in ls(pattern = "^RED_OPT_")) {
  df <- get(df_name)  # Retrieve dataframe
  reduce_columns(df, df_name)  # Apply downsampling function
}


################## MAKE COLUMNS IN LONG FORMAT ADD EXTRA VARIABLES 

# Loop through all dataframes with prefix "TRU_RED_OPT_"
for (df_name in ls(pattern = "^TRU_RED_OPT_")) {
  df <- get(df_name)  # Retrieve the dataframe
  
  # Extract the last four digits after the last underscore in the dataframe name (e.g., TRU_RED_OPT_1234 -> 1234)
  ID <- sub(".*_(\\d{4})$", "\\1", df_name)
  
  # Convert the dataframe to long format (without the TIME column)
  df_long <- data.frame(
    ID = rep(ID, nrow(df)),  # Repeat ID for each row in the dataframe
    OPT_COH_RAW = as.vector(unlist(df))  # Flatten the dataframe to a single column
  )
  
  # Compute the mean and standard deviation of OPT_COH_RAW, ignoring NAs and NaNs
  mean_coh_raw <- mean(df_long$OPT_COH_RAW, na.rm = TRUE)
  sd_coh_raw <- sd(df_long$OPT_COH_RAW, na.rm = TRUE)
  
  # Create OPT_COH_CEN by centering OPT_COH_RAW and replacing NAs/NaNs with 0
  df_long$OPT_COH_CEN <- ifelse(is.na(df_long$OPT_COH_RAW) | is.nan(df_long$OPT_COH_RAW), 
                               0,  # Replace NAs and NaNs with 0 instead of the mean
                               df_long$OPT_COH_RAW - mean_coh_raw)
  
  # Create OPT_COH_STD by standardizing OPT_COH_RAW and replacing NAs/NaNs with 0 (mean)
  df_long$OPT_COH_STD <- ifelse(is.na(df_long$OPT_COH_RAW) | is.nan(df_long$OPT_COH_RAW) | sd_coh_raw == 0, 
                               0,  # Replace NAs and NaNs in OPT_COH_STD with 0 instead of the mean
                               (df_long$OPT_COH_RAW - mean_coh_raw) / sd_coh_raw)
  
  # Assign the resulting long-format dataframe back to the environment with the "LONG_" prefix
  assign(paste0("LONG_", df_name), df_long)
}


################################################################################
################## JOIN ALL LONG DFs INTO ONE LARGE ONE 

# Initialize an empty list to store the data frames
long_dfs <- list()

# Loop through all dataframes with prefix "LONG_"
for (df_name in ls(pattern = "^LONG_TRU_RED_OPT_")) {
  df <- get(df_name)  # Retrieve the dataframe
  
  # Append the dataframe to the list
  long_dfs[[length(long_dfs) + 1]] <- df
}

# Combine all the data frames in the list into one large data frame
OPT_COH_JOINT <- do.call(rbind, long_dfs)


hist(OPT_COH_JOINT$OPT_COH_RAW)

################## ADD SAMPLE VARIABLES  


## RAW MEANS 

# Compute the OPT_MEAN_COH for each ID (mean of COH_RAW) and repeat for all rows within the same ID
OPT_COH_JOINT <- OPT_COH_JOINT %>%
  group_by(ID) %>%
  mutate(OPT_MEAN_COH = mean(OPT_COH_RAW, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after mutation to avoid future issues with grouped data

## GRAND MEAN CENTERED RAW MEANS 

# Compute the grand mean of OPT_MEAN_COH
grand_mean_coh <- mean(OPT_COH_JOINT$OPT_MEAN_COH, na.rm = TRUE)

# Create the MEAN_OPT_COH_CEN variable by centering OPT_MEAN_COH with the grand mean
OPT_COH_JOINT <- OPT_COH_JOINT %>%
  mutate(MEAN_OPT_COH_CEN = OPT_MEAN_COH - grand_mean_coh)


#### KEY

# COH_RAW: The unchanged weighted time series.
# MEAN_COH: The mean value of the COH_RAW variable 
# MEAN_COH_CEN: The grand mean centered MEAN_COH value.
# COH_CEN: The unchanged weighted time series that is centered within person.
#          Nas and NANs are replaced with the mean (0). 
# COH_STD: The unchanged weighted time series that is standardized within person. 
#          Nas and NANs are replaced with the mean (0)

################################################################################
################# SAVE THE SUMMARY TIME SERIES 


### SAVE TIME-INVARIANT OPTIMIZED MEAN SCORES 

# Define the file path where you want to save the CSV
OUTPUT_TVAR_OPT_COH <- file.path(work_dir, "TIME_VAR_OPT_COH_SCORES.csv")

# Save the data frame as a CSV file
write.csv(OPT_COH_JOINT, file = OUTPUT_TVAR_OPT_COH, row.names = FALSE)

