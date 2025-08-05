# --------------------------------------------------------------------------
#### Tasks #### - Emily Dorado

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".
# Done!

# 2. Create Project Folder
# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.
# Done! 

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# ---------------------------------------------------------------------------
# 3. Download "patient_info.csv" dataset from GitHub repository

# load the dataset into your R environment
# patient_info <- read.csv(file.choose()) # patient_info
patient_info <- read.csv("clean_data/patient_info_clean.csv")

# Inspect the structure of the dataset using appropriate R functions
str(patient_info) # structure of the dataset
summary(patient_info) # summary of the columns 
head(patient_info) # displays rows

# Identify variables with incorrect or inconsistent data types.
  # gender should be a factor not chr
  # diagnosis should be a factor not chr
  # smoker should be a factor not chr 

# Convert variables to appropriate data types where needed
patient_info$gender <- as.factor(patient_info$gender)
patient_info$diagnosis <- as.factor(patient_info$diagnosis)
patient_info$smoker <- as.factor(patient_info$smoker)

str(patient_info) #run this to make sure there's corrections

# Create a new variable for smoking status as a binary factor:
patient_info$smoker_binary <- ifelse(patient_info$smoker == "Yes", 1, 0)
patient_info$smoker_binary <- as.factor(patient_info$smoker_binary)

str(patient_info) # run this to make sure its correct) 
# looking at the table, i also see the new variable and 0/1s 

# 1 for "Yes", 0 for "No" #### reminder of what they mean 

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
write.csv(patient_info, file = "clean_data/patient_info_clean.csv")

#Save your R script in your script folder with name "class_Ib"
  # click save, rename the file, save it in scripts
# Upload "class_Ib" R script into your GitHub repository