

### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)

# Set working directory 


POUNDTOWN <- 1

# set working directory
if (POUNDTOWN == 1) {
  work_dir <- 'C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\WIDE_IBI_DATA\\'
} else {
  work_dir <- 'C:\\Users\\0910h\\Dropbox\\DISSERTATION\\ANALYSIS\\WIDE_IBI_DATA\\'
}

setwd(work_dir)

######################################################
########## LOAD IN THE ALL THE RAW IBI DATA ##########
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


######################################################
######## COMBINE THE IBI DATA INTO 1 VARIABLE ########
######################################################

data_frames <- Filter(function(x) is.data.frame(get(x)), ls())

# Function to process each data frame
combine_ibi_columns <- function(df_name) {
  # Get the data frame from the environment
  df <- get(df_name)
  
  # Extract the prefix (C or P) from the data frame name
  prefix <- substr(df_name, 1, 1)
  
  # Define the new variable name based on the prefix
  new_var_name <- paste0(prefix, "_IBI")
  
  # Select only IBI_* columns
  ibi_cols <- df[, grepl("^IBI_", names(df)), drop = FALSE]
  
  # Initialize the combined IBI vector
  combined_ibi <- numeric()
  
  if (ncol(ibi_cols) > 0) { # Ensure there are IBI columns
    for (i in 1:ncol(ibi_cols)) {
      # Remove NAs from the current column
      current_column <- ibi_cols[[i]][!is.na(ibi_cols[[i]])]
      
      if (length(current_column) > 0) { # Process only non-empty columns
        if (length(combined_ibi) > 0) { # Combine with the previous column
          last_value <- tail(combined_ibi, 1) # Get the last value of combined_ibi
          first_value <- current_column[1] # Get the first value of the current column
          
          if (last_value + first_value < 1350) {
            # Combine the last value with the first value of the current column
            combined_ibi[length(combined_ibi)] <- last_value + first_value
            # Add the rest of the current column
            combined_ibi <- c(combined_ibi, current_column[-1])
          } else {
            # Directly append the current column without combining
            combined_ibi <- c(combined_ibi, current_column)
          }
        } else {
          # First column, no need to combine
          combined_ibi <- c(combined_ibi, current_column)
        }
      }
    }
  }
  
  # Extract the ID column from the original data frame and repeat it
  id_column <- rep(df$ID, length.out = length(combined_ibi))
  
  # Create a new data frame with ID and the combined IBI column
  combined_df <- data.frame(ID = id_column, combined_ibi)
  names(combined_df)[2] <- new_var_name # Rename the IBI column dynamically
  
  # Save the new data frame back to the environment with a new name
  assign(paste0(df_name, "_combined"), combined_df, envir = .GlobalEnv)
}

# Apply the function to all data frames in the environment
data_frames <- Filter(function(x) is.data.frame(get(x)), ls()) # Get all data frame names
lapply(data_frames, combine_ibi_columns)




################################################################################
################# REMOVE WIDE IBI DFS and RETAIN THE LONG IBIs #################
################################################################################

# Get all objects in the environment
all_objects <- ls()

# Identify data frames with the suffix "_combined"
combined_dfs <- grep("_combined$", all_objects, value = TRUE)

# Remove all other objects
rm(list = setdiff(all_objects, combined_dfs), envir = .GlobalEnv)


################################################################################
########### CLEAN UP THE ID VARIABLE AND ADD A TIME VARIABLE ###################
################################################################################

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




## ADD A TIME VARIABLE BY SUMMING ALL THE IBIS

# Loop through each data frame with _combined suffix
for (df_name in combined_dfs) {
  # Get the data frame
  df <- get(df_name)
  
  # Identify the IBI column (the one that ends with "_IBI")
  ibi_column <- grep("_IBI$", names(df), value = TRUE)
  
  # Check if there is exactly one IBI column
  if (length(ibi_column) == 1) {
    # Calculate the cumulative time using the identified IBI column
    df$CUM_TIME <- cumsum(df[[ibi_column]])
    
    # Save the updated data frame back to the environment
    assign(df_name, df, envir = .GlobalEnv)
  } else {
    message(paste0("No '_IBI' column found or multiple '_IBI' columns found in ", df_name, ". Skipping."))
  }
}

### CONVERT MILLISECOND TIME TO TIME IN SECONDS 

# Loop through all data frames in the environment
for (df_name in ls()) {
  # Check if the object is a data frame
  if (is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Check if it contains the CUM_TIME column
    if ("CUM_TIME" %in% colnames(df)) {
      # Convert CUM_TIME to seconds and add it as a new column
      df$CUM_TIME_SEC <- df$CUM_TIME / 1000
      
      # Update the data frame in the environment
      assign(df_name, df)
    }
  }
}


## REMOVE WORKING DF
rm(df)

################################################################################
######################### REMOVE FIRST IBI #####################################
################################################################################

## REMOVE THE FIRST OBSERVATION BECAUSE IT DOESN'T HAVE A GOOD RELATIVE TIME DISTANCE 
## IT LEAdS TO WEIRD STARTING VALUES THAT ARE CERTAINLY SPECIOUS 


# Loop through all data frames in the environment
for (df_name in ls()) {
  # Check if the object is a data frame
  if (is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Identify the IBI column by suffix
    IBI_column <- names(df)[grepl("_IBI$", names(df))]
    
    # Ensure that exactly one IBI column is found
    if (length(IBI_column) == 1) {
      # Set the first IBI value to NA
      df[[IBI_column]][1] <- NA
      
      # Update the data frame in the environment
      assign(df_name, df)
    }
  }
}

################################################################################
############# IDENTIFY AND REMOVE EXCESSIVE/OUTLIER IBIs #######################
################################################################################

# We need to take out IBIs that are clearly incorrect. I used the threshold of 
# 3 SD above or below that persons mean. Visually, it seemed to be the best. They
# are marked as NA and estimated in the interpolation step.


BAD_IBIs_TO_NA <- function(df) {
  # Identify the IBI column by suffix
  IBI_column <- names(df)[grepl("_IBI$", names(df))]
  
  # Ensure that exactly one IBI column is found
  if (length(IBI_column) != 1) {
    stop("Error: Could not uniquely identify the IBI column in the data frame.")
  }
  
  # Extract the IBI values
  IBI_values <- df[[IBI_column]]
  
  # Compute mean and standard deviation for the individual's IBI series
  IBI_mean <- mean(IBI_values, na.rm = TRUE)
  IBI_sd <- sd(IBI_values, na.rm = TRUE)
  
  # Define dynamic thresholds (3 SD above and below the mean)
  lower_bound <- IBI_mean - 3 * IBI_sd
  upper_bound <- IBI_mean + 3 * IBI_sd
  
  # Set values outside these bounds to NA
  df[[IBI_column]][IBI_values < lower_bound | IBI_values > upper_bound] <- NA
  
  return(df)
}

# Get a list of all data frames in the environment
data_frames <- ls()  # List all objects in the environment
data_frames <- data_frames[sapply(data_frames, function(x) is.data.frame(get(x)))]

# Apply the function to each data frame to remove excessive IBIs
cleaned_dfs <- lapply(data_frames, function(df_name) {
  df <- get(df_name)
  BAD_IBIs_TO_NA(df)
})

# Optionally, save the cleaned data frames back to the environment
names(cleaned_dfs) <- paste0("CLEAN_", data_frames)
list2env(cleaned_dfs, .GlobalEnv)



######################################################
############# INTERPOLATE DATA TO 20hz ###############
######################################################

## USING A CUBIC SPLINE INTERPOLATION 

# Define a function to interpolate IBI data to 4 Hz using cubic spline
CUBIC_SPLINE_INT_IBI <- function(df) {
  # Identify the IBI column by suffix
  IBI_column <- names(df)[grepl("_IBI$", names(df))]
  
  # Ensure that exactly one IBI column is found
  if (length(IBI_column) != 1) {
    stop("Error: Could not uniquely identify the IBI column in the data frame.")
  }
  
  # Extract the original time and IBI series
  time_original <- df$CUM_TIME_SEC
  IBI_original <- df[[IBI_column]]
  
  # Create a new sequence of time points at 20 Hz (every 0.05 seconds)
  time_new <- seq(min(time_original), max(time_original), by = 0.05)
  
  # Perform cubic spline interpolation
  spline_interpolation <- spline(x = time_original, y = IBI_original, xout = time_new, method = "natural")
  
  # Create a new data frame with the interpolated data
  interpolated_df <- data.frame(
    ID = unique(df$ID),       # Keep the participant ID
    TIME_SEC = time_new,          # New regular time points
    IBI = spline_interpolation$y  # Interpolated IBI values
  )
  
  return(interpolated_df)
}

# Get a list of all the CLEAN data frames in the environment
data_frames <- ls()  # List all objects in the environment
data_frames <- data_frames[grepl("^CLEAN_", data_frames)]  

# Apply interpolation to each data frame and save the results
interpolated_dfs <- lapply(data_frames, function(df_name) {
  df <- get(df_name)
  CUBIC_SPLINE_INT_IBI(df)
})

# Optionally assign new names to the interpolated data frames
names(interpolated_dfs) <- paste0(data_frames, "_INTERPOLATED")

# Save interpolated data frames back to the environment
list2env(interpolated_dfs, .GlobalEnv)

###  REMOVE WORKING FILES 

# Get all objects in the environment
all_objects <- ls()

# Identify data frames with the suffix "_combined"
combined_dfs <- grep("_INTERPOLATED$", all_objects, value = TRUE)

# Remove all other objects
rm(list = setdiff(all_objects, combined_dfs), envir = .GlobalEnv)


################################################################################
############################ CENTER THE IBI SERIES  ############################
################################################################################


## CENTER THE TIME VARIABLE AT THE START 

# Loop through all data frames in the environment
for (df_name in ls()) {
  # Check if the object is a data frame
  if (is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Check if it contains the TIME_SEC column
    if ("TIME_SEC" %in% colnames(df)) {
      # Subtract the first value of TIME_SEC from all values in TIME_SEC
      df$TIME_SEC <- df$TIME_SEC - df$TIME_SEC[1]
      
      # Update the data frame in the environment
      assign(df_name, df)
    }
    
  }
}

## REMOVE WORKING DFs

rm(df)


################################################################################
########################### PLOT SOME EXAMPLES #################################
################################################################################


# Assuming your data frame is called df
ggplot(CLEAN_P_0011_combined_INTERPOLATED, aes(x = TIME_SEC, y = IBI)) +
  geom_line() +   # To create a line plot; use geom_point() for scatter plot
  labs(x = "Time (Seconds)", y = "C_IBI (Milliseconds)", title = "Time vs. C_IBI") +
  theme_minimal()  # For a clean plot theme



################################################################################
######################## ADD THE PREFIX TO THE IBI SERIES ######################
################################################################################

## ADD THE CORRECT PREFIX TO THE IBI SERIES 

# Loop through all data frames in the environment
for (df_name in ls()) {
  # Check if the object is a data frame
  if (is.data.frame(get(df_name))) {
    # Get the data frame
    df <- get(df_name)
    
    # Extract the identifier after the "CLEAN_" prefix
    identifier <- sub("^CLEAN_(.)_.*", "\\1", df_name)
    
    # Check if the data frame contains the IBI column
    if ("IBI" %in% colnames(df)) {
      # Rename the IBI column to include the identifier
      colnames(df)[colnames(df) == "IBI"] <- paste0(identifier, "_IBI")
      
      # Update the data frame in the environment
      assign(df_name, df)
    }
  }
}


## REMOVE WORKING dfs

rm(df)


################################################################################
######################## TEST AND REMOVE TIME TRENDS ###########################
################################################################################


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

# Loop through all data frames in the environment that end with _INT
for (name in ls()) {
  # Check if the object is a data frame and ends with "_INT"
  if (is.data.frame(get(name)) && grepl("_INTERPOLATED$", name)) {
    df <- get(name)  # Get the data frame by name
    
    # Print the column names to ensure the correct ones are being used
    print(paste("Checking columns for selected data frame:", name))
    print(names(df))  # Print the column names of the data frame
    
    # Check for required columns
    if (!any(grepl("_IBI$", names(df)))) {
      warning(paste("Data frame", name, "does not contain an IBI column. Skipping."))
      next
    }
    
    # Dynamically identify the IBI column (either C_IBI or P_IBI)
    ibi_col <- grep("_IBI$", names(df), value = TRUE)[1]  # Get the first match for IBI column
    print(paste("Using IBI column:", ibi_col))  # Print which IBI column is being used
    
    # Ensure the data frame has the 'TIME_SEC' column
    if (!"TIME_SEC" %in% names(df)) {
      warning(paste("Data frame", name, "is missing the 'TIME_SEC' column. Skipping."))
      next
    }
    
    # Fit linear model
    lm_linear <- lm(as.formula(paste(ibi_col, "~ TIME_SEC")), data = df)
    
    # Fit quadratic model
    lm_quadratic <- lm(as.formula(paste(ibi_col, "~ TIME_SEC + I(TIME_SEC^2)")), data = df)
    
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

## REMOVE LINEAR/QUADRATIC TIME TRENDS ##

# Loop through each row in the trend_results data frame
for (i in 1:nrow(trend_results)) {
  # Extract the data frame name (ID) and best model choice
  df_name <- trend_results$ID_Name[i]
  better_model <- trend_results$Better_Model[i]
  
  # Load the data frame by name
  df <- get(df_name)
  
  # Dynamically identify the IBI column (either C_IBI or P_IBI)
  ibi_col <- grep("_IBI$", names(df), value = TRUE)[1]
  
  # Ensure the data frame has the 'TIME_SEC' column
  if (!"TIME_SEC" %in% names(df)) {
    warning(paste("Data frame", df_name, "is missing the 'TIME_SEC' column. Skipping."))
    next
  }
  
  # Perform regression based on the chosen model (linear or quadratic)
  if (better_model == "Linear") {
    # Fit the linear model
    lm_linear <- lm(as.formula(paste(ibi_col, "~ TIME_SEC")), data = df)
    
    # Subtract the predicted values from the original IBI values to detrend
    df$detrended_IBI <- df[[ibi_col]] - predict(lm_linear)
    
  } else if (better_model == "Quadratic") {
    # Fit the quadratic model
    lm_quadratic <- lm(as.formula(paste(ibi_col, "~ TIME_SEC + I(TIME_SEC^2)")), data = df)
    
    # Subtract the predicted values from the original IBI values to detrend
    df$detrended_IBI <- df[[ibi_col]] - predict(lm_quadratic)
  }
  
  # Modify the detrended IBI column name based on the prefix (P_ or C_)
  if (grepl("^C_", df_name)) {
    new_col_name <- paste0("C_", sub("^C_", "", sub("_IBI$", "_IBI_D", ibi_col)))  # C_IBI_D
  } else if (grepl("^P_", df_name)) {
    new_col_name <- paste0("P_", sub("^P_", "", sub("_IBI$", "_IBI_D", ibi_col)))  # P_IBI_D
  } else {
    # If no matching prefix, use the default name (can be customized if needed)
    new_col_name <- sub("_IBI$", "_IBI_D", ibi_col)  # _IBI_D
  }
  
  # Rename the detrended IBI column
  names(df)[names(df) == "detrended_IBI"] <- new_col_name
  
  # Assign the modified data frame back to the global environment
  assign(df_name, df)  # This updates the data frame in the environment
}


## VISUALIZE THE RESULTS ##


# Assuming 'C_1002_INT' is your data frame and it already has 'detrended_IBI' column
# Create a plot with original and detrended IBI series

ggplot(CLEAN_C_1036_combined_INTERPOLATED, aes(x = TIME_SEC)) +
  geom_line(aes(y = C_IBI, color = "Original"), size = 1) +  # Original IBI series
  geom_line(aes(y = C_IBI_D, color = "Detrended"), size = 1, linetype = "dashed") +  # Detrended IBI series
  labs(title = "Original vs. Detrended IBI Series",
       x = "Time (seconds)",
       y = "IBI (ms)",
       color = "Series") +
  scale_color_manual(values = c("Original" = "blue", "Detrended" = "red")) +
  theme_minimal()


## REMOVE WORKIN dfS

rm(df)
rm(lm_linear)
rm(lm_quadratic)
rm(trend_results)
rm(CLEAN_df_INTERPOLATED)


################################################################################
############################ SAVE THE DATA FRAMES ##############################
################################################################################

#### Define the output directory

#RIVERS CROSSING
output_dir <- "C:\\Users\\cjh37695\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_IBI_20Hz\\"  

#HOME
#output_dir <- "C:\\Users\\0910h\\Dropbox\\Dropbox\\DISSERTATION\\ANALYSIS\\ALL_IBI_20Hz\\"


# Get all object names in the environment that are data frames
all_dfs <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Filter for data frames with the "CLEAN_" prefix and "_INTERPOLATED" suffix
CLEAN_INTERPOLATED_DFs <- all_dfs[grepl("^CLEAN_.*_INTERPOLATED$", all_dfs)]


# Loop through each data frame with "_INTERPOLATED" suffix
for (df_name in CLEAN_INTERPOLATED_DFs) {
  # Get the data frame
  df <- get(df_name)
  
  # Remove the "CLEAN_" prefix and the "_combined_INTERPOLATED" suffix from the data frame name
  file_name <- sub("^CLEAN_", "", df_name)  # Remove "CLEAN_" prefix
  file_name <- sub("_combined_INTERPOLATED$", "", file_name)  # Remove "_combined_INTERPOLATED" suffix
  
  # Add the "_INT_20Hz" suffix
  file_name <- paste0(file_name, "_IBI_20Hz")
  
  # Construct the full file path for saving
  file_path <- file.path(output_dir, paste0(file_name, ".csv"))
  
  # Save the data frame as a CSV
  write.csv(df, file_path, row.names = FALSE)
  
  # Print a message to confirm saving
  message(paste("Saved:", file_path))
}
