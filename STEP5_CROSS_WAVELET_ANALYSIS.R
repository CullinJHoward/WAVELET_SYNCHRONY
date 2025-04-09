
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
############## LOAD IN ALL CLEANED DYADIC DATA

# Define the subfolder path within work_dir
folder <- file.path(work_dir, "MERGED_DYADIC_IBI_20Hz")

# Get the list of .csv files in the subfolder
csv_files <- list.files(path = folder, pattern = "_DYAD.csv$", full.names = TRUE)

# Function to load, clean, and rename CSV files
load_and_clean_csv <- function(file_name) {
  # Extract the base file name (without the full path)
  base_name <- basename(file_name)
  
  # Extract the family number by removing "_DYAD.csv" and keeping only the family ID
  family_number <- sub("_DYAD.csv$", "", base_name)
  
  # Construct the new variable name with "D_" prefix and no suffix
  new_name <- paste0("D_", family_number)  # This ensures the name is D_####, e.g., D_1100
  
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
############ REDUCE TO JUST THE TASK WANTED

# THE MATH TASK IS 5 MINUTES (i.e., 6,000 observations [rows] @ 20 Hz)
## HOWEVER, I WILL REMOVE ANYONE WITH LESS THAN 7500 ROWS BECAUSE THEIR LENGTHS 
## ARE ANOMALOUS AND I DON'T TRUST THEY DID THE FULL TASK. 

# Get all objects in the environment that start with "D_"
df_names <- ls(pattern = "^D_")

# Identify DYADS with 7500 rows or less
ODD_DFs <- df_names[sapply(df_names, function(df_name) {
  obj <- get(df_name)  # Retrieve the object
  is.data.frame(obj) && nrow(obj) <= 7500  # Check if it's a data frame and has ≤ 7500 rows
})]

# Remove ODD DFs from the environment
rm(list = ODD_DFs)

# UPDATE THE LIST OF DFs 
df_names <- ls(pattern = "^D_")

# Loop through each data frame and trim to first 6,000 rows
for (df_name in df_names) {
  df <- get(df_name)  # Retrieve the data frame
  if (is.data.frame(df)) {  # Ensure it's a data frame
    assign(df_name, df[1:min(6000, nrow(df)), ], envir = .GlobalEnv)  # Trim and save back
  }
}


################################################################################
############ CROSS WAVELET COHERENCE ANALYSIS 


# List all data frames in the environment with the prefix "DYAD_"
dyad_dfs <- ls(pattern = "^D_")


################## DATA WITH RESPIRATION REGRESSED OUT #########################

# Create a directory for saving coherence data if it doesn't exist

if (!dir.exists("MATH_RESP_FREE_COHERENCE_DATA")) {
  dir.create("MATH_RESP_FREE_COHERENCE_DATA")
}
if (!dir.exists("MATH_RESP_FREE_COHERENCE_PLOTS")) {
  dir.create("MATH_RESP_FREE_COHERENCE_PLOTS")
}

# Loop through each dyad data frame
for (dyad in dyad_dfs) {
  # Get the dyad ID (remove the "DYAD_" prefix)
  dyad_id <- gsub("^D_", "", dyad)
  
  # Perform cross-wavelet coherence analysis
  RESP_FREE_WC <- analyze.coherency(
    get(dyad),
    my.pair = c("C_IBI_D_R", "P_IBI_D_R"), #D = TIME DETRENDED, R = RESPIRATION REGRESSED OUT 
    dt = 1/20,
    dj = .05, # CHANGE TO .01 FOR MORE PRECISION (MUCH LONGER COMPUTATIONAL TIME)
    lowerPeriod = 1,
    upperPeriod = 60,
    make.pval = TRUE,
    n.sim = 50
  )
  
  # Save the coherence analysis result as an .RDS file
  saveRDS(RESP_FREE_WC, file = paste0("MATH_RESP_FREE_COHERENCE_DATA/", dyad_id, "_MATH_RESP_FREE_COH_DATA.RDS"))
  
  # Save the cross-wavelet power spectrum plot (with COI respected) to a PDF
  pdf(file = paste0("MATH_RESP_FREE_COHERENCE_PLOTS/", dyad_id, "_MATH_RESP_FREE_CrossWaveletPower.pdf"))
  wc.image(
    RESP_FREE_WC,
    main = paste("Dyad ID:", dyad_id, "Cross-Wavelet Power Spectrum"),
    legend.params = list(lab = paste0(dyad_id, " cross-wavelet power levels")),
    timelab = "Seconds",
    periodlab = "Time Period",
    plot.coi = TRUE  # Include the cone of influence
  )
  dev.off() # Close the PDF device
  
  # Save the average power spectrum plot to a PDF
  pdf(file = paste0("MATH_RESP_FREE_COHERENCE_PLOTS/", dyad_id, "_MATH_RESP_FREE_AveragePowerSpectrum.pdf"))
  wc.avg(
    RESP_FREE_WC,
    siglvl = 0.01, 
    sigcol = "red", 
    sigpch = 20,
    periodlab = "Time Period",
    main = paste("Dyad ID:", dyad_id, "Average Power Spectrum")
  )
  dev.off() # Close the PDF device
}

################## DATA WITH RESPIRATION STILL PRESENT #########################

# Create a directory for saving coherence data if it doesn't exist

if (!dir.exists("MATH_RESP_INC_COHERENCE_DATA")) {
  dir.create("MATH_RESP_INC_COHERENCE_DATA")
}
if (!dir.exists("MATH_RESP_INC_COHERENCE_PLOTS")) {
  dir.create("MATH_RESP_INC_COHERENCE_PLOTS")
}

# Loop through each dyad data frame
for (dyad in dyad_dfs) {
  # Get the dyad ID (remove the "DYAD_" prefix)
  dyad_id <- gsub("^D_", "", dyad)
  
  # Perform cross-wavelet coherence analysis
  RESP_FREE_WC <- analyze.coherency(
    get(dyad),
    my.pair = c("C_IBI_D", "P_IBI_D"), #D = TIME DETRENDED
    dt = 1/20,
    dj = 0.05, # CHANGE TO .01 FOR MORE PRECISION (MUCH LONGER COMPUTATIONAL TIME)
    lowerPeriod = 1,
    upperPeriod = 60,
    make.pval = TRUE,
    n.sim = 50
  )
  
  # Save the coherence analysis result as an .RDS file
  saveRDS(RESP_FREE_WC, file = paste0("MATH_RESP_INC_COHERENCE_DATA/", dyad_id, "_MATH_RESP_INC_COH_DATA.RDS"))
  
  # Save the cross-wavelet power spectrum plot (with COI respected) to a PDF
  pdf(file = paste0("MATH_RESP_INC_COHERENCE_PLOTS/", dyad_id, "_MATH_RESP_INC_CrossWaveletPower.pdf"))
  wc.image(
    RESP_FREE_WC,
    main = paste("Dyad ID:", dyad_id, "Cross-Wavelet Power Spectrum"),
    legend.params = list(lab = paste0(dyad_id, " cross-wavelet power levels")),
    timelab = "Seconds",
    periodlab = "Time Period",
    plot.coi = TRUE  # Include the cone of influence
  )
  dev.off() # Close the PDF device
  
  # Save the average power spectrum plot to a PDF
  pdf(file = paste0("MATH_RESP_INC_COHERENCE_PLOTS/", dyad_id, "_MATH_RESP_INC_AveragePowerSpectrum.pdf"))
  wc.avg(
    RESP_FREE_WC,
    siglvl = 0.01, 
    sigcol = "red", 
    sigpch = 20,
    periodlab = "Time Period",
    main = paste("Dyad ID:", dyad_id, "Average Power Spectrum")
  )
  dev.off() # Close the PDF device
}





