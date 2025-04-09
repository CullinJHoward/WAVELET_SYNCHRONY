
### LIBRARY
library(dplyr)
library(ggplot2)

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
####################### LOAD IN IBI DATA #######################################
################################################################################

# Define the subfolder path within work_dir
folder <- file.path(work_dir, "ALL_CLEAN_IBI_20Hz")

# Get the list of .csv files in the subfolder
csv_files <- list.files(path = folder, pattern = "_CLEAN.csv$", full.names = TRUE)

# Function to load, clean, and rename CSV files
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the prefix (C_ or P_) and retain underscore
  prefix <- sub("(_.*)", "_", base_name)  # Keeps "C_" or "P_"
  
  # Extract the family number (digits before "_IBI_RESP.csv")
  family_number <- sub(".*_(\\d+)_CLEAN.csv$", "\\1", base_name)
  
  # Construct the new variable name with underscore
  new_name <- paste0(prefix, family_number)  # Example: "C_0001" or "P_0002"
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Remove columns that are entirely NA
  data <- data[, colSums(!is.na(data)) > 0]
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
}

# Apply the function to all CSV files
lapply(csv_files, load_and_clean_csv)



################################################################################
################### IDENTIFY DYADS AND INDIVIDUAL DATA  ########################
################################################################################

# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Extract family IDs and role (Parent or Child)
df_info <- data.frame(
  df_name = all_dfs,
  role = sub("_(\\d{4})$", "", all_dfs),  # Extract P_ or C_
  family_id = sub("^[PC]_(\\d{4}).*$", "\\1", all_dfs),  # Extract the 4-digit family ID
  stringsAsFactors = FALSE
)

# Identify unique family IDs
unique_families <- unique(df_info$family_id)

# Create lists for dyads and individuals
dyads <- c()
individuals <- c()

# Loop through each family ID
for (fam_id in unique_families) {
  # Find all data frames related to this family
  family_dfs <- df_info[df_info$family_id == fam_id, ]
  
  # Check if both a parent and child exist
  has_parent <- any(grepl("^P_", family_dfs$df_name))
  has_child <- any(grepl("^C_", family_dfs$df_name))
  
  if (has_parent & has_child) {
    dyads <- c(dyads, fam_id)  # Add to dyads list
  } else {
    individuals <- c(individuals, fam_id)  # Add to individuals list
  }
}

## DYADS
print(dyads)

## INDIVIDUALS 
print(individuals)


###################### REMOVE NON-DYAD DATA

# Remove data frames associated with families in the individuals list
for (fam_id in individuals) {
  # Find all data frames that match this family ID
  dfs_to_remove <- ls()[grepl(paste0("^[PC]_?", fam_id), ls())]
  
  # Remove them from the environment
  rm(list = dfs_to_remove, envir = .GlobalEnv)
}

# ENSURE dyads match the remaning dfs 
ls()
length(ls()[sapply(ls(), function(x) is.data.frame(get(x)))])


################################################################################
############### ENSURE ROW NUMBER MATCHES BETWEEN DYAD MEMBERS #################
################################################################################


# Loop through each family ID in the dyads list
for (family_id in dyads) {
  # Construct the data frame names with the appropriate prefixes
  parent_df_name <- paste0("P_", family_id)
  child_df_name <- paste0("C_", family_id)
  
  # Check if both parent and child data frames exist in the environment
  if (exists(parent_df_name, envir = .GlobalEnv) && exists(child_df_name, envir = .GlobalEnv)) {
    # Retrieve the data frames
    df_parent <- get(parent_df_name, envir = .GlobalEnv)
    df_child <- get(child_df_name, envir = .GlobalEnv)
    
    # Ensure both data frames have the same number of rows
    min_rows <- min(nrow(df_parent), nrow(df_child))
    df_parent <- df_parent[1:min_rows, , drop = FALSE]  # Trim the parent data frame
    df_child <- df_child[1:min_rows, , drop = FALSE]  # Trim the child data frame
    
    # Assign the modified data frames back to the environment
    assign(parent_df_name, df_parent, envir = .GlobalEnv)
    assign(child_df_name, df_child, envir = .GlobalEnv)
  }
}

################################################################################
################# MERGE EACH PERSONS IBI AND RESP DATAFRAMES ###################
################################################################################

# Function to merge IBI and RESP data for dyads
merge_dyad_data <- function() {
  # Loop through all family IDs (dyads) in the list
  for (family_id in dyads) {
    # Construct the data frame names for the parent and child
    parent_df_name <- paste0("P_", family_id)
    child_df_name <- paste0("C_", family_id)
    
    # Check if both parent and child data frames exist in the environment
    if (exists(parent_df_name, envir = .GlobalEnv) && exists(child_df_name, envir = .GlobalEnv)) {
      # Retrieve the data frames for parent and child
      df_parent <- get(parent_df_name, envir = .GlobalEnv)
      df_child <- get(child_df_name, envir = .GlobalEnv)
      
      # Ensure both data frames have the same number of rows (trim if necessary)
      min_rows <- min(nrow(df_parent), nrow(df_child))
      df_parent <- df_parent[1:min_rows, , drop = FALSE]
      df_child <- df_child[1:min_rows, , drop = FALSE]
      
      # Remove redundant 'ID' and 'TIME_SEC' columns from child data frame
      if ("ID" %in% colnames(df_parent) && "ID" %in% colnames(df_child) && all(df_parent$ID == df_child$ID)) {
        df_child <- df_child[, !(colnames(df_child) == "ID")]  # Remove ID from child data
      }
      
      if ("TIME_SEC" %in% colnames(df_parent) && "TIME_SEC" %in% colnames(df_child) && all(df_parent$TIME_SEC == df_child$TIME_SEC)) {
        df_child <- df_child[, !(colnames(df_child) == "TIME_SEC")]  # Remove TIME_SEC from child data
      }
      
      # Merge the parent and child data frames into one family data frame
      merged_df <- cbind(df_parent, df_child)
      
      # Create a new name for the merged family data frame
      new_name <- paste0(family_id, "_DYAD")  # Name it based on the family ID with "_DYAD"
      
      # Assign the merged family data frame back to the environment
      assign(new_name, merged_df, envir = .GlobalEnv)
    }
  }
}

# Call the function to merge IBI and RESP data for all dyads
merge_dyad_data()
View(`0002_DYAD`)

############# REMOVE WORKING DATAFRAMES 

# Get all object names in the environment
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Identify data frames that do not have the "_IBI_RESP" suffix
dfs_to_remove <- all_dfs[!grepl("_DYAD$", all_dfs)]

# Remove the identified data frames from the environment
invisible(lapply(dfs_to_remove, function(df_name) rm(list = df_name, envir = .GlobalEnv)))


################################################################################
############################ SAVE THE DATA FRAMES ##############################
################################################################################

#### Define the output directory

if (POUNDTOWN == 1) {
  output_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\MERGED_DYADIC_IBI_20Hz'
} else {
  output_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\MERGED_DYADIC_IBI_20Hz'
}

# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Loop through each data frame with "_IBI_RESP" suffix
for (df_name in all_dfs) {
  # Get the data frame
  df <- get(df_name)
  
  # Construct the full file path for saving, keeping the original name intact
  file_path <- file.path(output_dir, paste0(df_name, ".csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file_path, row.names = FALSE)
  
  # Print a message to confirm saving
  message(paste("Saved:", file_path))
}

