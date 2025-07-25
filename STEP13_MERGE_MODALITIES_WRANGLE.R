### LIBRARY
library(dplyr)
library(DescTools)
library(ggplot2)
library(MplusAutomation)
library(tidyr)
library(psych)
library(naniar)

# Set working directory 
POUNDTOWN <- 0

if (POUNDTOWN == 1) {
  work_dir <- 'D:\\DISSERTATION\\ANALYSIS\\'
} else {
  work_dir <- 'F:\\DISSERTATION\\ANALYSIS\\'
}
setwd(work_dir)

################################################################################
####################### LOAD IN SYNCHRONY DATA 


FOa_COH <- read.csv("NONLINEAR_PATTERN_MODELING\\FRQ_OPT_Va_NLM_COH_SUMMARY.csv")
FOb_COH <- read.csv("NONLINEAR_PATTERN_MODELING\\FRQ_OPT_Vb_NLM_COH_SUMMARY.csv")
HF_COH <- read.csv("NONLINEAR_PATTERN_MODELING\\HF_NLM_COH_SUMMARY.csv")
LF_COH <- read.csv("NONLINEAR_PATTERN_MODELING\\LF_NLM_COH_SUMMARY.csv")
FOa_PHS <- read.csv("NONLINEAR_PATTERN_MODELING\\FRQ_OPT_Va_NLM_PHS_SUMMARY.csv")
FOb_PHS <- read.csv("NONLINEAR_PATTERN_MODELING\\FRQ_OPT_Vb_NLM_PHS_SUMMARY.csv")
HF_PHS <- read.csv("NONLINEAR_PATTERN_MODELING\\HF_NLM_PHS_SUMMARY.csv")
LF_PHS <- read.csv("NONLINEAR_PATTERN_MODELING\\LF_NLM_PHS_SUMMARY.csv")

# ENSURE 4 DIGIT ID 
FOa_COH$ID <- sprintf("%04d", as.numeric(FOa_COH$ID))
FOb_COH$ID <- sprintf("%04d", as.numeric(FOb_COH$ID))
HF_COH$ID <- sprintf("%04d", as.numeric(HF_COH$ID))
LF_COH$ID <- sprintf("%04d", as.numeric(LF_COH$ID))
FOa_PHS$ID <- sprintf("%04d", as.numeric(FOa_PHS$ID))
FOb_PHS$ID <- sprintf("%04d", as.numeric(FOb_PHS$ID))
HF_PHS$ID <- sprintf("%04d", as.numeric(HF_PHS$ID))
LF_PHS$ID <- sprintf("%04d", as.numeric(LF_PHS$ID))


## MERGE ALL PHYSIOLOGY TOGETHER
FULL_PHY <- FOa_COH %>%
  full_join(FOb_COH, by = "ID") %>%
  full_join(HF_COH, by = "ID") %>%
  full_join(LF_COH, by = "ID") %>%
  full_join(FOa_PHS, by = "ID") %>%
  full_join(FOb_PHS, by = "ID") %>%
  full_join(HF_PHS, by = "ID") %>%
  full_join(LF_PHS, by = "ID") 
  
# REMOVE PARTICIPANTS WHO FAILED QC 

## QUANTITATIVE FAIL 
QC_DF <- read.csv("WAVELET_ANALYSES//FRQ_OPT_CLASSIFICATION_VERSION_A_SUMMARY.csv")
QC_DF$ID <- sprintf("%04d", as.numeric(QC_DF$ID))
QUANT_FAIL_ID <- QC_DF$ID[QC_DF$OPT_CRIT == "QC_FAIL"]

## QUALITATIVE FAIL
QUAL_FAIL_ID <- c("0031", "0046", "0098", "1011", "1018", "1019", "1035", "1083")

# COMBIN QUAL/QUANT QC FAILS
QC_FAIL_ID <- union(QUANT_FAIL_ID, QUAL_FAIL_ID)

## REMOVE THE BAD PHYSIO DATA FROM THE FULL SET 

FULL_PHY_CLEAN <- FULL_PHY  %>%
  filter(!ID %in% QC_FAIL_ID)

# FIX THE STRUCTURE SO THAT NUMBERS ARE NUMERIC 
FULL_PHY_CLEAN[ , 2:41] <- FULL_PHY_CLEAN[ , 2:41] |> 
  lapply(function(x) round(as.numeric(x), 3))


################################################################################
####################### WRANGLE PHYSIO DATA 

## VISUALIZE SOME DATA 


# LIST VARIABLES OF INTEREST 
my_vars <- c("A_FOa_A", "A_FOa_B", "A_FOa_C", "A_FOb_A", "A_FOb_B",
             "A_FOb_C", "A_HF_A", "A_HF_B", "A_HF_C",   "A_LF_A",     
             "A_LF_B", "A_LF_C")  
df <- FULL_PHY_CLEAN_ABS  # Replace with your actual data frame

# Step 1: Convert to long format for just the selected variables
df_long <- df %>%
  select(all_of(my_vars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Step 2: Compute min and max for each variable
min_max_labels <- df_long %>%
  group_by(Variable) %>%
  summarize(
    min_val = min(Value, na.rm = TRUE),
    max_val = max(Value, na.rm = TRUE),
    label = paste0("Min: ", round(min_val, 2), "\nMax: ", round(max_val, 2))
  )

# Step 3: Merge labels into long data
df_long <- left_join(df_long, min_max_labels, by = "Variable")

# Step 4: Plot
ggplot(df_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ Variable, scales = "free") +
  geom_text(
    data = min_max_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal() +
  labs(title = "Histograms with Min/Max")

## ABSOLUTE VALUES

# THIS IS BECAUSE ALL THE ESTIMATES CAN BE POSITIVE OR NEGATIVE. TO TEST IF IT IS 
# THE DEGREE OF THESE PATTERNS WE NEED THEM ALL THE SAME DIRECTION. WE'LL TAKE 
# THEIR ABOSLUTE VALUE. WE CAN ALSO TEST IF POSITIVE/NEGATIVE MATTERS WHICH WOULD 
# BE WHEN THEY MOVE (TEMPORALLY) BUT THAT IS NOT MY CURRENT AIM

# LIST VARIABLES TO ABSOLUTE VALUE
my_vars <- c("FOa_A",     "FOa_B",     "FOa_C",
             "FOb_A",     "FOb_B",     "FOb_C",
             "HF_A",      "HF_B",      "HF_C",
             "LF_A",      "LF_B",      "LF_C",
             "FOa_PHSa",  "FOa_PHSb",  "FOa_PHSc",
             "FOb_PHSa", "FOb_PHSb",  "FOb_PHSc",
             "HF_PHSa",   "HF_PHSb",   "HF_PHSc",
             "LF_PHSa",   "LF_PHSb",   "LF_PHSc")


# TAKE THEIR ABS VALUES 
FULL_PHY_CLEAN_ABS <- FULL_PHY_CLEAN %>%
  mutate(across(all_of(my_vars), ~ abs(.), .names = "A_{.col}"))

# MULITPLE THE AMPLITUDE VALUES BY 100 TO KEEP EVERYTHING ON A SIMILAR SCALE
names(FULL_PHY_CLEAN_ABS)
my_vars_T <- c("A_FOa_A", 
             "A_FOb_A", 
             "A_HF_A", 
             "A_LF_A",
             "A_FOa_PHSa", 
             "A_FOb_PHSa",
             "A_HF_PHSa", 
             "A_LF_PHSa",
             "FOa_MU",
             "FOb_MU",
             "HF_MU",
             "LF_MU",
             "FOa_PHSm",
             "FOb_PHSm",
             "HF_PHSm",
             "LF_PHSm")


# ULTIPLY BY 100
FULL_PHY_CLEAN_ABS_T <- FULL_PHY_CLEAN_ABS %>%
  mutate(across(all_of(my_vars_T), ~ . * 100, .names = "T_{.col}"))

## WINSORIZE VARIABLES

# LOSE MASKING 
detach("package:DescTools", unload = TRUE, character.only = TRUE)
detach("package:ggplot2", unload = TRUE, character.only = TRUE)
detach("package:MplusAutomation", unload = TRUE, character.only = TRUE)
detach("package:tidyr", unload = TRUE, character.only = TRUE)
detach("package:psych", unload = TRUE, character.only = TRUE)
#BRING BACK THE PACKAGE YOU NEED 
library(psych)


FULL_PHY_CLEAN_ABS_T_WINS <- FULL_PHY_CLEAN_ABS_T %>% 
  mutate(WT_A_FOa_A = psych::winsor(T_A_FOa_A, trim = 0.05)) %>%
  mutate(WA_FOa_B = psych::winsor(A_FOa_B, trim = 0.05)) %>%
  mutate(WA_FOa_C = psych::winsor(A_FOa_C, trim = 0.05)) %>%
  mutate(WT_A_FOb_A = psych::winsor(T_A_FOb_A, trim = 0.05)) %>%
  mutate(WA_FOb_B = psych::winsor(A_FOb_B, trim = 0.05)) %>%
  mutate(WA_FOb_C = psych::winsor(A_FOb_C, trim = 0.05)) %>%
  mutate(WT_A_HF_A = psych::winsor(T_A_HF_A, trim = 0.05)) %>%
  mutate(WA_HF_B = psych::winsor(A_HF_B, trim = 0.05)) %>%
  mutate(WA_HF_C = psych::winsor(A_HF_C, trim = 0.05)) %>%
  mutate(WT_A_LF_A = psych::winsor(T_A_LF_A, trim = 0.05)) %>%
  mutate(WA_LF_B = psych::winsor(A_LF_B, trim = 0.05)) %>%
  mutate(WA_LF_C = psych::winsor(A_LF_C, trim = 0.05)) %>%
  mutate(WT_A_FOa_PHSa = psych::winsor(T_A_FOa_PHSa, trim = 0.05)) %>%
  mutate(WA_FOa_PHSb = psych::winsor(A_FOa_PHSb, trim = 0.05)) %>%
  mutate(WA_FOa_PHSc = psych::winsor(A_FOa_PHSc, trim = 0.05)) %>%
  mutate(WT_A_FOb_PHSa = psych::winsor(T_A_FOb_PHSa, trim = 0.05)) %>%
  mutate(WA_FOb_PHSb = psych::winsor(A_FOb_PHSb, trim = 0.05)) %>%
  mutate(WA_FOb_PHSc = psych::winsor(A_FOb_PHSc, trim = 0.05)) %>%
  mutate(WT_A_HF_PHSa = psych::winsor(T_A_HF_PHSa, trim = 0.05)) %>%
  mutate(WA_HF_PHSb = psych::winsor(A_HF_PHSb, trim = 0.05)) %>%
  mutate(WA_HF_PHSc = psych::winsor(A_HF_PHSc, trim = 0.05)) %>%
  mutate(WT_A_LF_PHSa = psych::winsor(T_A_LF_PHSa, trim = 0.05)) %>%
  mutate(WA_LF_PHSb = psych::winsor(A_LF_PHSb, trim = 0.05)) %>%
  mutate(WA_LF_PHSc = psych::winsor(A_LF_PHSc, trim = 0.05))

## SEE HOW IT DID 

#BRING THEM BACK (LOLZ)
library(dplyr)
library(DescTools)
library(ggplot2)
library(MplusAutomation)
library(tidyr)

# LIST VARIABLES OF INTEREST 
WINS_VARS <- c("WT_A_FOa_A",     "WA_FOa_B",     "WA_FOa_C",
             "WT_A_FOb_A",     "WA_FOb_B",     "WA_FOb_C",
             "WT_A_HF_A",      "WA_HF_B",      "WA_HF_C",
             "WT_A_LF_A",      "WA_LF_B",      "WA_LF_C",
             "WT_A_FOa_PHSa",  "WA_FOa_PHSb",  "WA_FOa_PHSc",
             "WT_A_FOb_PHSa", "WA_FOb_PHSb",  "WA_FOb_PHSc",
             "WT_A_HF_PHSa",   "WA_HF_PHSb",   "WA_HF_PHSc",
             "WT_A_LF_PHSa",   "WA_LF_PHSb",   "WA_LF_PHSc")  


# Step 1: Convert to long format for just the selected variables
df_long <- FULL_PHY_CLEAN_ABS_T_WINS %>%
  select(all_of(WINS_VARS)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Step 2: Compute min and max for each variable
min_max_labels <- df_long %>%
  group_by(Variable) %>%
  summarize(
    min_val = min(Value, na.rm = TRUE),
    max_val = max(Value, na.rm = TRUE),
    label = paste0("Min: ", round(min_val, 2), "\nMax: ", round(max_val, 2))
  )

# Step 3: Merge labels into long data
df_long <- left_join(df_long, min_max_labels, by = "Variable")

# Step 4: Plot
ggplot(df_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ Variable, scales = "free") +
  geom_text(
    data = min_max_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal() +
  labs(title = "Histograms with Min/Max")


## CENTER VARIABLES FOR MODERATION 
names(FULL_PHY_CLEAN_ABS_T_WINS)

# LIST OF VARIABLES TO CENTER 
CENTER_VARS <- c("WT_A_FOa_A", "WA_FOa_B", "WA_FOa_C",     
                    "WT_A_FOb_A", "WA_FOb_B", "WA_FOb_C",
                    "WT_A_HF_A",  "WA_HF_B",  "WA_HF_C",      
                    "WT_A_LF_A",  "WA_LF_B",  "WA_LF_C",
                    "WT_A_FOa_PHSa", "WA_FOa_PHSb",   "WA_FOa_PHSc",  
                    "WT_A_FOb_PHSa", "WA_FOb_PHSb",   "WA_FOb_PHSc",
                    "WT_A_HF_PHSa",  "WA_HF_PHSb",    "WA_HF_PHSc",   
                    "WT_A_LF_PHSa",  "WA_LF_PHSb",    "WA_LF_PHSc",
                    "T_FOa_MU", "T_FOb_MU", "T_HF_MU", "T_LF_MU",
                    "T_FOa_PHSm", "T_FOb_PHSm", "T_HF_PHSm", "T_LF_PHSm")

# MEAN CENTER VARIABLES 
FULL_PHY_CLEAN_ABS_T_WINS_CENT <- FULL_PHY_CLEAN_ABS_T_WINS %>%
  mutate(across(all_of(CENTER_VARS), 
                ~ .x - mean(.x, na.rm = TRUE), 
                .names = "C_{.col}"))

## REDUCE TO THE VARIABLES FOR THE ANALYSIS

PHYIO_RED <- FULL_PHY_CLEAN_ABS_T_WINS_CENT %>%
  select(c("ID", "C_WT_A_FOa_A", "C_WA_FOa_B", "C_WA_FOa_C",     
           "C_WT_A_FOb_A", "C_WA_FOb_B", "C_WA_FOb_C",
           "C_WT_A_HF_A",  "C_WA_HF_B",  "C_WA_HF_C",      
           "C_WT_A_LF_A",  "C_WA_LF_B",  "C_WA_LF_C",
           "C_WT_A_FOa_PHSa", "C_WA_FOa_PHSb",   "C_WA_FOa_PHSc",  
           "C_WT_A_FOb_PHSa", "C_WA_FOb_PHSb",   "C_WA_FOb_PHSc",
           "C_WT_A_HF_PHSa",  "C_WA_HF_PHSb",    "C_WA_HF_PHSc",   
           "C_WT_A_LF_PHSa",  "C_WA_LF_PHSb",    "C_WA_LF_PHSc",
           "C_T_FOa_MU", "C_T_FOb_MU", "C_T_HF_MU", "C_T_LF_MU",
           "C_T_FOa_PHSm", "C_T_FOb_PHSm", "C_T_HF_PHSm", "C_T_LF_PHSm"))

## REVERSE NAMES FOR MPLUS EASE

PHYIO_RED_RENAMED <- PHYIO_RED %>%
  rename(FOa_A_atwC = C_WT_A_FOa_A,
         FOa_B_awC = C_WA_FOa_B, 
         FOa_C_awC = C_WA_FOa_C, 
         FOb_A_atwC = C_WT_A_FOb_A, 
         FOb_B_awC = C_WA_FOb_B, 
         FOb_C_awC = C_WA_FOb_C, 
         HF_A_atwC = C_WT_A_HF_A,
         HF_B_awC = C_WA_HF_B, 
         HF_C_awC = C_WA_HF_C,
         LF_A_atwC = C_WT_A_LF_A, 
         LF_B_awC = C_WA_LF_B,
         LF_C_awC = C_WA_LF_C,
         FOa_PHSa_atwC = C_WT_A_FOa_PHSa,
         FOa_PHSb_awC = C_WA_FOa_PHSb,
         FOa_PHSc_awC = C_WA_FOa_PHSc,
         FOb_PHSa_atwC = C_WT_A_FOb_PHSa,
         FOb_PHSb_awC = C_WA_FOb_PHSb,
         FOb_PHSc_atwC = C_WA_FOb_PHSc,
         HF_PHSa_atwC = C_WT_A_HF_PHSa,
         HF_PHSb_awC = C_WA_HF_PHSb,
         HF_PHSc_awC = C_WA_HF_PHSc,  
         LF_PHSa_atwC = C_WA_LF_PHSb,
         LF_PHSc_atwC = C_WA_LF_PHSc,
         FOa_MU_tC = C_T_FOa_MU, 
         FOb_MU_tC = C_T_FOb_MU, 
         HF_MU_tC  = C_T_HF_MU, 
         LF_MU_tC = C_T_LF_MU,
         FOa_PHSm_tC = C_T_FOa_PHSm, 
         FOb_PHSm_tC = C_T_FOb_PHSm,
         HF_PHSm_tC = C_T_HF_PHSm, 
         LF_PHSm_tC = C_T_LF_PHSm
         )

names(PHYIO_RED_RENAMED)


################################################################################
####################### LOAD IN SURVEY DATA 

# DORRY ONLY 
DORRY_SURVEY <- read.csv("SURVEY_MEASUREMENT//DISSERTATION_ALL_SURVEY_DATA_DORRY_7.6.25.csv")

## DROP UNNECESSARY DORRY DATA (FROM PRIOR FACTOR FORMATION)
names(DORRY_SURVEY)

#GET W2 AGE
AGE_W2 <- read.csv("SURVEY_MEASUREMENT\\DORRY_ONLY\\SURVEY_DATA\\DORRY WAVE 2_CHILD_AGE.csv")

AGE_W2 <- AGE_W2 %>%
  rename(C_AGE2 = DEM3PCW2)

DORRY_SURVEY_RED <- DORRY_SURVEY %>%
  select(c("ID", "P_SEX", "P_EDU", "C_SEX", "C_AGE", "C_RACE", "MAR_STAT", "HOUSE_NUMW1",
           "INCOME_W1", "SUBSEScm", "SUBSESus",
           "C_pTHR_Y", "C_pUNP_Y", "C_pDEP_Y", "C_pCUM_Y", "C_pTHR_P", "C_pDEP_P",         
           "C_pCUM_P", "Q_C_pTHR_Y", "Q_C_pUNP_Y", "Q_C_pDEP_Y", "Q_C_pCUM_Y",       
           "Q_C_pTHR_P", "Q_C_pDEP_P", "Q_C_pCUM_P", "Pfac_P1", "Pfac_P2",          
           "Pfac_Y1", "Pfac_Y2", "Pfac_PLI", "Pfac_PLC", "Pfac_YLI",         
           "Pfac_YLC", "SOACC_1", "SOACC_2", "SOACC_LI", "SOACC_LC",         
           "PERTK_1", "EMPCN_1", "PERTK_2", "EMPCN_2", "PERTK_LI",         
           "PERTK_LC", "EMPCN_LI", "EMPCN_LC")) %>%
  full_join(AGE_W2, by = "ID")

DORRY_SURVEY_RED <- DORRY_SURVEY_RED %>%
  rename(C_AGE1 = C_AGE)

names(DORRY_SURVEY_RED)

################################################################################
###################### MERGING DATA 

#### REDUCE TO ONLY THE DORRY PHYSIO DATA 

# MAKE ID NUMERIC
PHYIO_RED_RENAMED$ID <- as.numeric(PHYIO_RED_RENAMED$ID)

# REDUCE TO JUST DORRY (ID 1000+)
DORRY_PHY <- subset(PHYIO_RED_RENAMED, ID > 999)

## MERGE SURVEY AND PHYSIO

ALL_DATA_ONLY_PHYSIO <- left_join(DORRY_PHY, DORRY_SURVEY_RED, by = "ID")


## TEST IF I SHOULD USE THE FULL SAMPLE BY TESTING MCAR 
# We can use the full sample if we would like. It is missing MCAR 


ALL_DATA_FULL_SAMPLE <-  left_join(DORRY_SURVEY_RED, DORRY_PHY, by = "ID")

#VARIABLES TO TEST 
mcar_vars <- ALL_DATA_FULL_SAMPLE[, c("C_pUNP_Y", "C_pDEP_Y",  "C_pTHR_Y",
                           "Pfac_Y1", "Pfac_Y2",
                           "FOa_MU_tC", "FOb_MU_tC", "HF_MU_tC", "LF_MU_tC",      
                           "FOa_PHSm_tC", "FOb_PHSm_tC", "HF_PHSm_tC",  "LF_PHSm_tC",
                           "C_SEX", "C_RACE" , "C_AGE", "SUBSESus", "INCOME_W1")]


# Run Little's MCAR test

mcar_test(mcar_vars) #STAT = 173 (147), p = .07

################################################################################
############################# SAVE DATA 

# DEFINE FILE PATH
OUTPUT_PATH <- file.path(work_dir, "/STRUCTURAL_EQUATION_MODELS//")

# SAVE DATA 

#ONLY COMPLETE 
write.csv(ALL_DATA_ONLY_PHYSIO, file = paste0(OUTPUT_PATH,"DORRY_onlyCOMPLETE_DATA_7.23.25.csv"), row.names = FALSE, na="")

prepareMplusData(ALL_DATA_ONLY_PHYSIO,paste0(OUTPUT_PATH,"DORRY_onlyCOMPLETE_DATA_7.23.25.dat"))

#FULL SAMPLE 
write.csv(ALL_DATA_FULL_SAMPLE, file = paste0(OUTPUT_PATH,"DORRY_FULL_DATA_7.23.25.csv"), row.names = FALSE, na="")

prepareMplusData(ALL_DATA_FULL_SAMPLE,paste0(OUTPUT_PATH,"DORRY_FULL_DATA_7.23.25.dat"))

################################################################################
####################### GET SOME KEY DESCRIPTIVES 

## AGE 
mean(ALL_DATA_FULL_SAMPLE$C_AGE1, na.rm = TRUE) # FULL 12.89
mean(ALL_DATA_ONLY_PHYSIO$C_AGE1, na.rm = TRUE) # FULL 13.03

mean(ALL_DATA_FULL_SAMPLE$C_AGE2, na.rm = TRUE) # FULL 15.01
sd(ALL_DATA_FULL_SAMPLE$C_AGE2, na.rm = TRUE) # FULL 1.07
mean(ALL_DATA_ONLY_PHYSIO$C_AGE2, na.rm = TRUE) # FULL 15.28
sd(ALL_DATA_ONLY_PHYSIO$C_AGE2, na.rm = TRUE) # FULL .93

## PHYSIO 

#COHERENCE 
mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WT_A_FOa_A, na.rm = TRUE) #6.67
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WT_A_FOa_A, na.rm = TRUE) #1.82

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_B, na.rm = TRUE) #66.57
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_B, na.rm = TRUE) # 51.42

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_C, na.rm = TRUE) #110.79
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_C, na.rm = TRUE) #63.21

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$T_FOa_MU, na.rm = TRUE) #83.77 
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$T_FOa_MU, na.rm = TRUE) #8.71

#PHASE
mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WT_A_FOa_PHSa, na.rm = TRUE) #94.87
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WT_A_FOa_PHSa, na.rm = TRUE) #29.15

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_PHSb, na.rm = TRUE) #105.83
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_PHSb, na.rm = TRUE) # 93.99

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_PHSc, na.rm = TRUE) #120.97
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$WA_FOa_PHSc, na.rm = TRUE) #80.82

mean(FULL_PHY_CLEAN_ABS_T_WINS_CENT$T_FOa_PHSm, na.rm = TRUE) #318.08
sd(FULL_PHY_CLEAN_ABS_T_WINS_CENT$T_FOa_PHSm, na.rm = TRUE) #181.36

## DIMENSIONAL ADVERSITY 

# FULL
mean(DORRY_SURVEY$pTHR_Y, na.rm = TRUE) #13.34
sd(DORRY_SURVEY$pTHR_Y, na.rm = TRUE) #5.89

mean(DORRY_SURVEY$pDEP_Y, na.rm = TRUE) #15.29
sd(DORRY_SURVEY$pDEP_Y, na.rm = TRUE) #10.52

mean(DORRY_SURVEY$pUNP_Y, na.rm = TRUE) #9.27
sd(DORRY_SURVEY$pUNP_Y, na.rm = TRUE) #6.27

mean(DORRY_SURVEY$pCUM_Y, na.rm = TRUE) #12.91
sd(DORRY_SURVEY$pCUM_Y, na.rm = TRUE) #5.90

#SUBSAMPLE
DESC_DORRY <- left_join(ALL_DATA_ONLY_PHYSIO, DORRY_SURVEY, by = "ID")

mean(DESC_DORRY$pTHR_Y, na.rm = TRUE) #14.02
sd(DESC_DORRY$pTHR_Y, na.rm = TRUE) #4.82

mean(DESC_DORRY$pDEP_Y, na.rm = TRUE) #15.43
sd(DESC_DORRY$pDEP_Y, na.rm = TRUE) #11.24

mean(DESC_DORRY$pUNP_Y, na.rm = TRUE) #9.64
sd(DESC_DORRY$pUNP_Y, na.rm = TRUE) #6.14

mean(DESC_DORRY$pCUM_Y, na.rm = TRUE) #13.51
sd(DESC_DORRY$pCUM_Y, na.rm = TRUE) #5.61
