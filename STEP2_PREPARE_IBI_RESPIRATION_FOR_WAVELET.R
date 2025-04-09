
### LIBRARY
library(dplyr)
library(ggplot2)

# Set working directory 


POUNDTOWN <- 1

# Set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\DISSERTATION\\ANALYSIS\\'
}

setwd(work_dir)

################################################################################
####################### LOAD IN IBI DATA #######################################
################################################################################

# Dynamically create the path to the "ALL_IBI_20Hz" subfolder based on work_dir
IBI_folder <- file.path(work_dir, "ALL_IBI_20Hz")

# Get the list of .csv files in the subfolder "ALL_IBI_20Hz"
IBI_csv_files <- list.files(path = IBI_folder, pattern = "\\_IBI_20Hz.csv$", full.names = TRUE)

# Function to create the desired name, load the file, and remove all-NA columns
IBI_load_and_clean_csv <- function(file_name) {
  # Extract the prefix (C_ or P_)
  prefix <- sub("_.*", "", basename(file_name))
  
  # Extract the family number (e.g., "0001") before the "_IBI_20Hz.csv" part
  family_number <- sub("^.*_(\\d+)_IBI_20Hz\\.csv$", "\\1", basename(file_name))
  
  # Create the new name by combining the prefix and family number, then adding "_IBI"
  new_name <- paste0(prefix, "_", family_number, "_IBI")  # Combine prefix and family number with suffix "_IBI"
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Remove columns that are entirely NA
  data <- data[, colSums(!is.na(data)) > 0]
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
}

# Apply the function to all files
lapply(IBI_csv_files, IBI_load_and_clean_csv)


################################################################################
####################### LOAD IN RESPIRATION DATA ###############################
################################################################################

# Dynamically create the path to the "ALL_IBI_20Hz" subfolder based on work_dir
RESP_folder <- file.path(work_dir, "ALL_RESP_20Hz")

# Get the list of .csv files in the subfolder "ALL_IBI_20Hz"
RESP_csv_files <- list.files(path = RESP_folder, pattern = "\\_RESP_20Hz.csv$", full.names = TRUE)

# Function to create the desired name, load the file, and remove all-NA columns
RESP_load_and_clean_csv <- function(file_name) {
  # Extract the prefix (C_ or P_)
  prefix <- sub("_.*", "", basename(file_name))
  
  # Extract the family number (e.g., "0001") before the "_IBI_20Hz.csv" part
  family_number <- sub("^.*_(\\d+)_RESP_20Hz\\.csv$", "\\1", basename(file_name))
  
  # Create the new name by combining the prefix and family number, then adding "_IBI"
  new_name <- paste0(prefix, "_", family_number, "_RESP")  # Combine prefix and family number with suffix "_RESP"
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Remove columns that are entirely NA
  data <- data[, colSums(!is.na(data)) > 0]
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
}

# Apply the function to all files
lapply(RESP_csv_files, RESP_load_and_clean_csv)


################################################################################
############### IDENTIFY EACH PERSONS RESPIRATION AND IBI DATA #################
################################################################################

# Get all object names in the environment
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]


# Extract full prefixes (C_XXXX or P_XXXX) from data frame names
full_prefixes <- unique(gsub("^(C_\\d{4}|P_\\d{4})_.*$", "\\1", all_dfs))

# Identify which families have both IBI and RESP data (both for child or parent)
paired_PHYSIO <- full_prefixes[sapply(full_prefixes, function(prefix) {
  # Check if both IBI and RESP data exist for the same individual (either child or parent)
  ibi_data <- paste0(prefix, "_IBI") %in% all_dfs  # Check for IBI data
  resp_data <- paste0(prefix, "_RESP") %in% all_dfs  # Check for RESP data
  
  # Both IBI and RESP data must exist for the same individual
  ibi_data & resp_data
})]

# Identify non-paired individuals (those that don't have both IBI and RESP data)
non_paired_PHYSIO <- setdiff(full_prefixes, paired_PHYSIO)

## PAIRED DATA
print(paired_PHYSIO)

## UNPAIRED DATA
print(non_paired_PHYSIO)

###################### REMOVE THOSE MISSING IBI OR RESPIRATION

# Create separate lists for IBI and RESP data frames for non-paired individuals
non_paired_IBI <- paste0(non_paired_PHYSIO, "_IBI")  # Data frames with _IBI suffix
non_paired_RESP <- paste0(non_paired_PHYSIO, "_RESP")  # Data frames with _RESP suffix

# Combine both lists (IBI and RESP)
non_paired_PARTS <- c(non_paired_IBI, non_paired_RESP)


# Check which of the non-paired data frames actually exist in the environment
existing_dfs <- intersect(non_paired_PARTS, ls())  # Identify which data frames are actually present

# Remove the existing non-paired data frames from the environment
rm(list = existing_dfs)

################################################################################
################# ENSURE ROW NUMBER MATCHES BETWEEN IBI & RESP #################
################################################################################

# Get all object names in the environment
all_dfs_v2 <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Function to extract the numeric family ID from the data frame name (including "C_" or "P_" prefix)
extract_id <- function(df_name) {
  match <- regmatches(df_name, regexpr("[CP]_\\d{4}", df_name))
  if (length(match) > 0) return(match) else return(NA)
}

# Create a dictionary to store data frames by family ID
df_dict <- list()

# Loop through all data frames and group them by family ID (including prefix)
for (df_name in all_dfs_v2) {
  if (exists(df_name, envir = .GlobalEnv) && is.data.frame(get(df_name))) {
    id <- extract_id(df_name)
    
    if (!is.na(id)) {
      if (!id %in% names(df_dict)) {
        df_dict[[id]] <- list()
      }
      df_dict[[id]][[df_name]] <- get(df_name)
    }
  }
}

# Process and ensure row count matches for IBI and RESP data of each individual
for (id in names(df_dict)) {
  df_names <- names(df_dict[[id]])
  
  # We only process individuals with both IBI and RESP data
  ibi_data <- df_names[grepl("_IBI$", df_names)]
  resp_data <- df_names[grepl("_RESP$", df_names)]
  
  if (length(ibi_data) == 1 && length(resp_data) == 1) {
    # Extract the IBI and RESP data frames for this individual
    df_ibi <- df_dict[[id]][[ibi_data]]
    df_resp <- df_dict[[id]][[resp_data]]
    
    # Ensure both data frames have the same number of rows
    min_rows <- min(nrow(df_ibi), nrow(df_resp))
    df_ibi <- df_ibi[1:min_rows, , drop = FALSE]  # Trim IBI data
    df_resp <- df_resp[1:min_rows, , drop = FALSE]  # Trim RESP data
    
    # Assign the modified data frames back to the environment
    assign(ibi_data, df_ibi, envir = .GlobalEnv)
    assign(resp_data, df_resp, envir = .GlobalEnv)
  }
}

################################################################################
################# MERGE EACH PERSONS IBI AND RESP DATAFRAMES ###################
################################################################################

# Function to merge IBI and RESP data for each individual
merge_ibi_resp <- function() {
  # Loop through all data frames and group them by family ID (including prefix)
  for (df_name in all_dfs) {
    if (exists(df_name, envir = .GlobalEnv) && is.data.frame(get(df_name))) {
      id <- extract_id(df_name)
      
      if (!is.na(id)) {
        # Create a dictionary to store data frames by family ID
        if (!id %in% names(df_dict)) {
          df_dict[[id]] <- list()
        }
        df_dict[[id]][[df_name]] <- get(df_name)
      }
    }
  }
  
  # Process and ensure that the data frames for each individual are paired
  for (id in names(df_dict)) {
    df_names <- names(df_dict[[id]])
    
    # Identify IBI and RESP data frames for the current individual
    ibi_data <- df_names[grepl("_IBI$", df_names)]
    resp_data <- df_names[grepl("_RESP$", df_names)]
    
    if (length(ibi_data) == 1 && length(resp_data) == 1) {
      # Extract the IBI and RESP data frames
      df_ibi <- df_dict[[id]][[ibi_data]]
      df_resp <- df_dict[[id]][[resp_data]]
      
      # Ensure both data frames have the same number of rows
      min_rows <- min(nrow(df_ibi), nrow(df_resp))
      df_ibi <- df_ibi[1:min_rows, , drop = FALSE]  # Trim IBI data
      df_resp <- df_resp[1:min_rows, , drop = FALSE]  # Trim RESP data
      
      # Remove one of the redundant columns (ID or TIME_SEC)
      if ("ID" %in% colnames(df_ibi) && "ID" %in% colnames(df_resp) && all(df_ibi$ID == df_resp$ID)) {
        df_resp <- df_resp[, !(colnames(df_resp) == "ID")]  # Remove ID from RESP data if identical
      }
      
      if ("TIME_SEC" %in% colnames(df_ibi) && "TIME_SEC" %in% colnames(df_resp) && all(df_ibi$TIME_SEC == df_resp$TIME_SEC)) {
        df_resp <- df_resp[, !(colnames(df_resp) == "TIME_SEC")]  # Remove TIME_SEC from RESP data if identical
      }
      
      # Merge the two data frames side by side (column bind)
      merged_df <- cbind(df_ibi, df_resp)
      
      # Create a new name for the merged data frame
      new_name <- paste0(id, "_IBI_RESP")  # Name it based on the ID and indicate IBI_RESP merge
      
      # Assign the merged data frame to the environment
      assign(new_name, merged_df, envir = .GlobalEnv)
    }
  }
}

# Call the function to merge IBI and RESP data for all paired individuals
merge_ibi_resp()


############# REMOVE WORKING DATAFRAMES 

# Get all object names in the environment
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Identify data frames that do not have the "_IBI_RESP" suffix
dfs_to_remove <- all_dfs[!grepl("_IBI_RESP$", all_dfs)]

# Remove the identified data frames from the environment
invisible(lapply(dfs_to_remove, function(df_name) rm(list = df_name, envir = .GlobalEnv)))


rm(df_dict)


################################################################################
############################ SAVE THE DATA FRAMES ##############################
################################################################################

#### Define the output directory

#RIVERS CROSSING
output_dir <- "C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_MERGED_IBI_RESP_20Hz\\"  

#HOME
#output_dir <- "C:\\Users\\0910h\\Dropbox\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_MERGED_IBI_RESP_20Hz\\"


# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Filter for data frames with "_IBI_RESP" suffix
IBI_RESP_DFs <- all_dfs[grepl("_IBI_RESP$", all_dfs)]

# Loop through each data frame with "_IBI_RESP" suffix
for (df_name in IBI_RESP_DFs) {
  # Get the data frame
  df <- get(df_name)
  
  # Construct the full file path for saving, keeping the original name intact
  file_path <- file.path(output_dir, paste0(df_name, ".csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file_path, row.names = FALSE)
  
  # Print a message to confirm saving
  message(paste("Saved:", file_path))
}

