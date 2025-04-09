################################################################################
##################### LIBRARY
library(forecast)
library(dplyr)
library(ggplot2)


################################################################################
##################### Set working directory

POUNDTOWN <- 0


if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
} else {
  work_dir <- 'D:\\Dropbox\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
}

setwd(work_dir)


################################################################################
##################### LOAD IN THE DATA 

# Define the subfolder path within work_dir
IBI_RESP_folder <- file.path(work_dir, "ALL_MERGED_IBI_RESP_20Hz")

# Get the list of .csv files in the subfolder
csv_files <- list.files(path = IBI_RESP_folder, pattern = "_IBI_RESP.csv$", full.names = TRUE)

# Function to load, clean, and rename CSV files
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the prefix (C_ or P_) and retain underscore
  prefix <- sub("(_.*)", "_", base_name)  # Keeps "C_" or "P_"
  
  # Extract the family number (digits before "_IBI_RESP.csv")
  family_number <- sub(".*_(\\d+)_IBI_RESP.csv$", "\\1", base_name)
  
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
#################### USE Autoregressive Integrated Moving Average (ARIMA)

## I originally tested this using the differencing order and applying it if needed
## not everyone needed it, most were 0. But the ones that did use it completely flattened
## the line. I think it is erroneous considering I have already removed linear 
## and quadratic effects of time on overall drift. So this code does not apply 
## a differencing approach. 


# Get a list of all data frames in the environment
df_list <- ls()

for (df_name in df_list) {
  df <- get(df_name)  # Load the data frame
  
  # Find IBI and RESP column names dynamically
  ibi_col <- grep("IBI_D", names(df), value = TRUE)
  resp_col <- grep("RESP_D", names(df), value = TRUE)
  
  # Ensure we found the required columns
  if (length(ibi_col) == 1 && length(resp_col) == 1) {
    
    # Set differencing order to 0 (no differencing applied)
    d_order <- 0
    
    # Print message for no differencing
    message(paste("No differencing applied to", df_name))
    
    # Fit ARIMAX model with automatic parameter selection, no differencing (d = 0)
    model <- auto.arima(df[[ibi_col]], xreg = df[[resp_col]], d = d_order)
    
    # Extract residuals (IBI with respiration regressed out)
    residuals_name <- paste0(sub("_.*", "", df_name), "_IBI_D_R")  # Create residual column name based on prefix
    df[[residuals_name]] <- residuals(model)
    
    # Save the modified data frame back to the environment
    assign(df_name, df)
    
    cat("Processed:", df_name, "\n")  # Print progress
  } else {
    cat("Skipping", df_name, "- IBI or RESP column not found\n")
  }
}

################################################################################
######################## VISUALIZE SOME EXAMPLES ###############################
################################################################################

#SELECT AN IBI TO INSPECT
TEST_DF <- C_1110

#PLOT THEM
ggplot(reshape2::melt(data.frame(Time = TEST_DF$TIME_SEC, 
                                 Original_IBI = TEST_DF$C_IBI_D, 
                                 Corrected_IBI = TEST_DF$C_IBI_D_R), 
                      id.vars = "Time"), 
       aes(x = Time, y = value, color = variable)) +
  geom_line(aes(linetype = variable), size = 1) +  # Add linetype mapping
  scale_linetype_manual(values = c("dashed", "solid")) +  # Make the Corrected IBI line dashed
  labs(title = "IBI Correction",
       x = "Time (Sec)",
       y = "IBI Value",
       color = "Series") +
  theme_minimal()

################################################################################
############################ SAVE THE DATA FRAMES 

## REMOVE WORKING DFs
rm(df)
rm(model)
rm(TEST_DF)


if (POUNDTOWN == 1) {
  output_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\'
} else {
  output_dir <- 'F:\\DISSERTATION\\ANALYSIS\\WAVELET_ANALYSES\\ALL_CLEAN_IBI_20Hz'
}

# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Loop through each data frame with "_IBI_RESP" suffix
for (df_name in all_dfs) {
  # Get the data frame
  df <- get(df_name)
  
  # Construct the full file path for saving, keeping the original name intact
  file_path <- file.path(output_dir, paste0(df_name, "_CLEAN.csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file_path, row.names = FALSE)
  
  # Print a message to confirm saving
  message(paste("Saved:", file_path))
}
