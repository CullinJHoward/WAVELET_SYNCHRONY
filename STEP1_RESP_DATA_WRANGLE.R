
### LIBRARY
library(dplyr)
library(ggplot2)

# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\WIDE_RESP_DATA\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\DISSERTATION\\ANALYSIS\\WIDE_RESP_DATA\\'
}

setwd(work_dir)

######################################################
########## LOAD IN THE ALL THE RAW RESP DATA #########
######################################################

# Get the list of .csv files in the working directory
csv_files <- list.files(pattern = "\\.csv$")

# Function to create the desired name, load the file, and remove all-NA columns
load_and_clean_csv <- function(file_name) {
  # Extract the prefix (first part before the first underscore)
  prefix <- sub("_.*", "", file_name)
  
  # Extract the numeric suffix (digits before .csv)
  suffix <- sub(".*_(\\d+)\\.csv$", "\\1", file_name)
  
  # Ensure suffix is exactly 4 digits by padding with leading zeros
  suffix <- sprintf("%04d", as.integer(suffix))
  
  new_name <- paste0(prefix, "_", suffix)  # Combine prefix and formatted suffix
  
  # Read the CSV file
  data <- read.csv(file_name, stringsAsFactors = FALSE)
  
  # Remove columns that are entirely NA
  data <- data[, colSums(!is.na(data)) > 0]
  
  # Assign the cleaned data frame to the new name in the environment
  assign(new_name, data, envir = .GlobalEnv)
}

# Apply the function to all files
lapply(csv_files, load_and_clean_csv)


# Ensure numbers are displayed in standard notation
options(scipen = 999)

################################################################################
################## COMBINE THE RESP DATA INTO 1 VARIABLE #######################
################################################################################

#IDENFITY ALL dfs to work through 

data_frames <- ls()[grepl("_\\d{4}$", ls())]

combine_RESP_columns <- function(df_name) {
  # Get the data frame from the environment
  df <- get(df_name)
  
  # Extract the prefix (C or P) from the data frame name
  prefix <- substr(df_name, 1, 1)
  
  # Define the new variable name based on the prefix
  new_var_name <- paste0(prefix, "_RESP")
  
  # Select the ID column
  ID_col <- df$ID
  
  # Select only RESP_* columns (excluding ID)
  RESP_cols <- df[, grepl("^RESP", names(df)), drop = FALSE]
  
  # Check if there are valid columns to process
  if (ncol(RESP_cols) == 0) {
    message(paste("No RESP_* columns found in", df_name, "- Skipping."))
    return(NULL)  # Skip processing if no valid columns
  }
  
  # Remove trailing NAs from each RESP_* column and stack them
  stacked_RESP <- do.call(rbind, lapply(RESP_cols, function(col) {
    valid_index <- max(which(!is.na(col)), na.rm = TRUE) # Find last non-NA
    if (is.finite(valid_index)) {
      col <- col[seq_len(valid_index)]  # Keep meaningful NAs
    }
    data.frame(ID = ID_col[seq_along(col)], RESP_value = col)  # Retain ID
  }))
  
  # Ensure the stacked dataset is not empty
  if (nrow(stacked_RESP) == 0) {
    message(paste("All RESP_* columns in", df_name, "are empty after processing - Skipping."))
    return(NULL)
  }
  
  # Rename RESP_value column dynamically
  names(stacked_RESP)[2] <- new_var_name
  
  # Save the new data frame back to the environment
  assign(paste0(df_name, "_combined"), stacked_RESP, envir = .GlobalEnv)
}

# Apply the function to all data frames in the environment
data_frames <- Filter(function(x) is.data.frame(get(x)), ls()) # Get all data frame names
lapply(data_frames, combine_RESP_columns)


################################################################################
################# REMOVE WIDE RESP DFS and RETAIN THE LONG RESP ################
################################################################################

# Get all objects in the environment
all_objects <- ls()

# Identify data frames with the suffix "_combined"
combined_dfs <- grep("_combined$", all_objects, value = TRUE)

# Remove all other objects
rm(list = setdiff(all_objects, combined_dfs), envir = .GlobalEnv)


## CLEAN UP THE ID VARIABLE 

# Loop through all objects in the environment
for (df_name in ls()) {
  # Check if the object is a data frame
  if (is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Check if the 'ID' column exists in the data frame
    if ("ID" %in% names(df)) {
      # Fill the entire 'ID' column with the first value
      df$ID <- df$ID[1]
      
      # Assign the modified data frame back to the global environment
      assign(df_name, df, envir = .GlobalEnv)
    }
  }
}


################################################################################
#### TRANSFORM THE VALUES TO A COMMON SCALE Z- normalize them ##################
################################################################################


# Z-transform function to standardize numeric values
z_transform_values <- function(df) {
  excluded_dfs <- list()  # Initialize an empty list to track excluded data frames
  
  # Identify the relevant RESP column (C_RESP or P_RESP)
  resp_col <- grep("_RESP$", names(df), value = TRUE)  # This will match both C_RESP or P_RESP
  
  if (length(resp_col) == 1) {  # If we find exactly one RESP column (either C_RESP or P_RESP)
    col_name <- resp_col  # The RESP column name
    col_sd <- sd(df[[col_name]], na.rm = TRUE)  # Calculate the standard deviation
    
    if (col_sd == 0 || sum(!is.na(df[[col_name]])) < 2) {
      # If standard deviation is zero or not enough valid data, skip the column and mark for exclusion
      excluded_dfs <<- c(excluded_dfs, df_name)  # Add the df name to the excluded list
      return(df)  # Return the unmodified data frame (no transformation)
    } else {
      # Z-transform (subtract mean and divide by standard deviation)
      df[[col_name]] <- (df[[col_name]] - mean(df[[col_name]], na.rm = TRUE)) / col_sd
    }
  } else {
    # If no RESP column is found or more than one, leave the df unchanged
    excluded_dfs <<- c(excluded_dfs, df_name)
    return(df)
  }
  
  return(df)
}

# Get names of all objects in the environment that match the "_combined" suffix
df_names <- ls()[grepl("_combined$", ls())]

# Initialize a list to track which data frames have no variability
excluded_dfs_list <- list()

# Apply Z-transformation to the data frames
for (df_name in df_names) {
  df <- get(df_name)
  if (is.data.frame(df)) {
    result_df <- z_transform_values(df)
    
    # If there were excluded columns (no variability), add to the exclusion list
    if (length(result_df$excluded_dfs) > 0) {
      excluded_dfs_list[[df_name]] <- result_df$excluded_dfs
    } else {
      # Update the data frame in the environment if transformed
      assign(df_name, result_df, envir = .GlobalEnv)
    }
  }
}

# Print out the excluded data frames that had no variability
print("Excluded data frames (no variability in RESP column):")
print(excluded_dfs_list)


################################################################################
############################# VISUALIZE SOME DATA ##############################
################################################################################

## RESET ROW NUMBERS


# Reset row numbers for each data frame
for (df_name in combined_dfs) {
  df <- get(df_name)
  rownames(df) <- 1:nrow(df)  # Reset the row names to sequential numbers
  assign(df_name, df, envir = .GlobalEnv)  # Save the modified df back to the environment
}


## VISUIALIZE BREATHING 

P_1118_combined$TIME <- as.numeric(rownames(P_1118_combined))


TEST <- P_1118_combined[1:500, ]

ggplot(TEST, aes(x = TIME, y = P_RESP)) + 
  geom_point() +  # For a line plot
  # geom_point() +  # Uncomment if you prefer a scatter plot
  labs(x = "TIME", y = "RESPIRATION", title = "Time vs. P_RESP") +
  theme_minimal()  # Optional for a clean theme

## remove my TIME variables

P_1118_combined$TIME <- NULL


################################################################################
################### ADD TIME, THEN TEST AND REMOVE TIME TRENDS #################
################################################################################

#REMOVE WORKING DFs
rm(df)
rm(TEST)
rm(result_df)
rm(excluded_dfs_list)

## ADD A TIME_SEC columns that counts the running time in seconds. 

# Loop through all data frames in the list all_dfs
for (df_name in combined_dfs) {
  # Check if the object exists and is a data frame
  if (exists(df_name, envir = .GlobalEnv) && is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Add the TIME_SEC column
    df$TIME_SEC <- seq(0, by = 0.05, length.out = nrow(df))
    
    # Assign the updated data frame back to the environment
    assign(df_name, df, envir = .GlobalEnv)
  }
}


## TEST FOR LINEAR AND QUADRATIC TIME TRENDS ##

# Initialize an empty data frame for storing results (if not already initialized)
if (!exists("trend_results")) {
  trend_results <- data.frame(
    DataFrame = character(),
    ID_Name = character(),
    Linear_Slope = numeric(),
    Quadratic_Slope = numeric(),
    Linear_P_Value = numeric(),
    Quadratic_P_Value = numeric(),
    Better_Model = character(),
    stringsAsFactors = FALSE
  )
}

# Loop through all data frames in the environment whose names are in combined_dfs
for (name in combined_dfs) {
  # Check if the object exists and is a data frame
  if (exists(name, envir = .GlobalEnv) && is.data.frame(get(name))) {
    df <- get(name)  # Retrieve the data frame
    
    # Print the column names to ensure the correct ones are being used
    print(paste("Checking columns for selected data frame:", name))
    print(names(df))  # Print the column names of the data frame
    
    # Check for required columns
    if (!any(grepl("_RESP$", names(df)))) {
      warning(paste("Data frame", name, "does not contain an RESP column. Skipping."))
      next
    }
    
    # Dynamically identify the IBI column (either C_IBI or P_IBI)
    resp_col <- grep("_RESP$", names(df), value = TRUE)[1]  # Get the first match for IBI column
    print(paste("Using RESP column:", resp_col))  # Print which IBI column is being used
    
    # Ensure the data frame has the 'TIME_SEC' column
    if (!"TIME_SEC" %in% names(df)) {
      warning(paste("Data frame", name, "is missing the 'TIME_SEC' column. Skipping."))
      next
    }
    
    # Fit linear model
    lm_linear <- lm(as.formula(paste(resp_col, "~ TIME_SEC")), data = df)
    
    # Fit quadratic model
    lm_quadratic <- lm(as.formula(paste(resp_col, "~ TIME_SEC + I(TIME_SEC^2)")), data = df)
    
    # Extract coefficients and p-values
    linear_slope <- coef(lm_linear)["TIME_SEC"]
    linear_p <- summary(lm_linear)$coefficients["TIME_SEC", "Pr(>|t|)"]
    
    quadratic_slope <- coef(lm_quadratic)["I(TIME_SEC^2)"]
    quadratic_p <- summary(lm_quadratic)$coefficients["I(TIME_SEC^2)", "Pr(>|t|)"]
    
    # Compare models (AIC or F-test)
    aic_linear <- AIC(lm_linear)
    aic_quadratic <- AIC(lm_quadratic)
    better_model <- ifelse(aic_quadratic < aic_linear, "Quadratic", "Linear")
    
    # Append results with the ID (which is the data frame name) for reference
    trend_results <- rbind(trend_results, data.frame(
      DataFrame = name,
      ID_Name = name,
      Linear_Slope = linear_slope,
      Quadratic_Slope = quadratic_slope,
      Linear_P_Value = linear_p,
      Quadratic_P_Value = quadratic_p,
      Better_Model = better_model,
      stringsAsFactors = FALSE
    ))
  }
}


# Print the final results
print(trend_results)

## REMOVE LINEAR/QUADRATIC TIME TRENDS BASED ON MODEL FIT ##

# Loop through each row in the trend_results data frame
for (i in 1:nrow(trend_results)) {
  # Extract the data frame name (ID) and best model choice
  df_name <- trend_results$ID_Name[i]
  better_model <- trend_results$Better_Model[i]
  
  # Load the data frame by name
  df <- get(df_name)
  
  # Dynamically identify the RESP column (either C_RESP or P_RESP)
  resp_col <- grep("_RESP$", names(df), value = TRUE)[1]
  
  # Ensure the data frame has the 'TIME_SEC' column
  if (!"TIME_SEC" %in% names(df)) {
    warning(paste("Data frame", df_name, "is missing the 'TIME_SEC' column. Skipping."))
    next
  }
  
  # Remove rows with missing values in TIME_SEC or RESP column
  df_clean <- df[!is.na(df[[resp_col]]) & !is.na(df$TIME_SEC), ]
  
  # Perform regression based on the chosen model (linear or quadratic)
  if (better_model == "Linear") {
    # Fit the linear model
    lm_linear <- lm(as.formula(paste(resp_col, "~ TIME_SEC")), data = df_clean)
    
    # Predict for all rows in df, handling missing values
    df$DETREND_RESP <- df[[resp_col]] - predict(lm_linear, newdata = df)
    
  } else if (better_model == "Quadratic") {
    # Fit the quadratic model
    lm_quadratic <- lm(as.formula(paste(resp_col, "~ TIME_SEC + I(TIME_SEC^2)")), data = df_clean)
    
    # Predict for all rows in df, handling missing values
    df$DETREND_RESP <- df[[resp_col]] - predict(lm_quadratic, newdata = df)
  }
  
  # Modify the DETREND_RESP column name based on the prefix (P_ or C_)
  if (grepl("^C_", df_name)) {
    new_col_name <- paste0("C_", sub("^C_", "", sub("_RESP$", "_RESP_D", resp_col)))  # C_RESP_D
  } else if (grepl("^P_", df_name)) {
    new_col_name <- paste0("P_", sub("^P_", "", sub("_RESP$", "_RESP_D", resp_col)))  # P_RESP_D
  } else {
    # If no matching prefix, use the default name
    new_col_name <- sub("_RESP$", "_RESP_D", resp_col)  # _RESP_D
  }
  
  # Rename the detrended RESP column
  names(df)[names(df) == "DETREND_RESP"] <- new_col_name
  
  # Assign the modified data frame back to the global environment
  assign(df_name, df, envir = .GlobalEnv)
}


## VISUALIZE THE RESULTS ##

## TAKE A SMALLER TIME SET
TEST <- C_1110_combined[1:500, ]

ggplot(TEST, aes(x = TIME_SEC)) +
  geom_line(aes(y = C_RESP, color = "Original"), size = 1) +  # Original IBI series
  geom_line(aes(y = C_RESP_D, color = "Detrended"), size = 1, linetype = "dashed") +  # Detrended IBI series
  labs(title = "Original vs. Detrended IBI Series",
       x = "Time (seconds)",
       y = "RESPIRATION",
       color = "Series") +
  scale_color_manual(values = c("Original" = "blue", "Detrended" = "red")) +
  theme_minimal()


## REMOVE WORKING DFs
rm(TEST)
rm(trend_results)
rm(df)
rm(df_clean)
rm(lm_linear)
rm(lm_quadratic)




################################################################################
################################ SAVE DFs ######################################
################################################################################

#### Define the output directory

#RIVERS CROSSING
output_dir <- "C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_RESP_20Hz\\"  

#HOME
#output_dir <- "C:\\Users\\0910h\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_RESP_20Hz\\"



# Filter for data frames with the "_combined" suffix
RESP_df_names <- ls()[grepl("_combined$", ls())]

# Loop through each data frame with "_combined" suffix
for (df_name in RESP_df_names) {
  # Get the data frame
  df <- get(df_name)
  
  # Remove the "_combined" suffix from the data frame name
  file_name <- sub("_combined$", "", df_name)
  
  # Add the "_INT_20Hz" suffix
  file_name <- paste0(file_name, "_RESP_20Hz")
  
  # Construct the full file path for saving
  file_path <- file.path(output_dir, paste0(file_name, ".csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file_path, row.names = FALSE)
  
  # Print a message to confirm saving
  message(paste("Saved:", file_path))
}

