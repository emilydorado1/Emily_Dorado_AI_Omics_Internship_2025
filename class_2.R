# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
# ===================================================================

# --------------------
# Topics in this module
# --------------------
#   1. Operators in R
#   2. Data Structures in R
#   3. User-defined Functions in R
#   4. Automating Workflows with for-Loops
# --------------------------------------------------------------------------------------------------
#### 1. Operators in R ####

# Operators are special symbols in R.
# They instruct R to perform actions on values or variables.
# You will use them for assignment, calculations, comparisons, and logical tests.

# ----------------------
# Assignment Operators
# ----------------------
# Used to store values inside variables

# <- (Most common, rightward assignment operator)
height <- c(1.75, 1.76, 1.82, 1.67)

# -> (Same as above, but leftward assignment operator)
c(68, 78, 85, 75) -> weight

# = ( also assigns, used in function arguments)
smoking_status = c("Yes", "No", "No", "Yes")

# ----------------------
# Arithmetic Operators
# ----------------------
# Perform basic math: +, -, *, /, ^
#   +  addition
#   -  subtraction
#   *  multiplication
#   /  division
#   ^  exponent (to the power of)

# let's BMI using weight and height

BMI <- weight/(height^2)
BMI

# Note: R applies operations element-wise when variables are vectors.
# This is called "vectorization", every weight is divided by every height squared.


# ---------------------
# Comparison Operators
# ---------------------
# Comparison operators ask TRUE/FALSE questions about values.
# They do not calculate a number, they return a logical output.(TRUE or FALSE)
# They compare values

# > greater than
BMI > 25 
# < less than
BMI < 18.5
# >= greater than or equal to
height >= 1.75
# <= less than or equal to
weight <= 65
# == equal to
smoking_status == "No"
# != not equal to 
smoking_status != "No"
# In R 
# yes = TRUE 
# no = FALSE
# -------------------
# Logical Operators
# -------------------
# Combine multiple conditions using:
#   &   AND
#   |   OR
#   !   NOT

# & AND (both must be TRUE)
# Is the patient overweight AND a smoker
# BMI cutoff = 25
(BMI > 25) & (smoking_status == "Yes")
(BMI < 25) & (smoking_status == "Yes")

# | OR (at least one must be TRUE)
# Is the patient overrweight or a smoker 
(BMI > 25) | (smoking_status == "Yes")
BMI 
smoking_status

# | NOT (reverse the condition)
# Is the pateint Not a smoker
!smoking_status == "No"
# condiiton = yes
# output = FALSE

# ------------------------
# 2. Data Structures in R
# ------------------------
# Data structures are how R organizes information.
# We commonly use:
#   1. Vectors
#   2. Lists
#   3. Matrices
#   4. Data Frames

# ---------
# Vectors
# ---------
# Simplest structure in R. Stores a sequence of values of the SAME type.

# - Numeric vector
num_vec <- c(1, 2, 3, 4)
class(num_vec)
# numeric vectors used to perform mathematical calculation
# - charecter vector
chrc_vector <- c("gene1", "gene2", "gene3")

# - logical vector 
logical_vector <- c(TRUE, FALSE, TRUE)

mix_vector <- c("gene1", 1, "gene2", 2)
mean(mix_vector)

# we can extract values from vectors using indexing with []
num_vec[2]
num_vec[2:4] # : indicates sequence

# you can only combine vectors of equal sequence
# you can treat vectors as columns or rows
# by column 
vec_col <- cbind(num_vec, chrc_vector)

# ---------
# Lists
# ---------
# Unlike vectors, lists can store different types together.

all_vectors <- list(num_vec, chrc_vector, logical_vector)

# save your raw_data
# save processed data
# results

# we access elements with [[ ]]
all_vectors[[2]]

# ---------
# Matrices
# ---------
# A 2D table where all values must be the same type.
# exampleL gene ecpression matrix where rows are genes and
# columns are samples

my_matrix <- matrix(1:9, nrow = 3, ncol = 3)

my_matrix
# By default R fills the matrix column wise
# we change using byrow = TRUE
my_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE)
my_matrix
# we access elements with [row, column]
my_matrix[2, 3]
my_matrix[2,]
my_matrix[#rows, #columns]

# ---------
# Data Frames
# ---------
# Most important structure for biological datasets.
# Columns can have different types (numeric, character, factor).

data <- data.frame(
  patient_id = c("P1", "P2", "P3"),
  age = c(65, 78, NA), 
  diagnosis = c("cancer", "diabetes", "cancer")
)

print(data)

# --------------------
# Dataset Assessment
# --------------------
# Before analyzing, inspect your dataset to understand its structure.
str(data) # structure of the dataset
head(data) # first 6 rows
head(data, n=2)
tail(data) # last 6 rows
tail(data, n =1)
dim(data) # rows & column
names(data)

# data frame are indexed like matrix with mroe flexibilty
# acess a column directly
data$patient_id
data[c(1,3), c(2,3)]

# ----------------
# Missing Values
# ----------------
# Real data often contains missing values (NA).
# You must check and handle them before analysis.

# qw can: 
# - detect them: is.na()
is.na(data)
# - count them: sum(is.na())
sum(is.na(data))
# - missing values by column
colSums(is.na(data))
# - missing values by rows
rowSums(is.na(data))
# - remove them: na.omit()
# remove rows with missing values
clean_data_1 <- na.omit(data)
clean_data_1
# - remove colum with missing value
clean_data_2 <- data[, colSums(is.na(data))==0]
clean_data_2
# - replace them: input with 0 or column mean
clean_data_3 <- data
clean_data_3[is.na(clean_data_3)] <- 0
clean_data_3
data

clean_data_4 <- data
clean_data_4[is.na(clean_data_4)] <- mean(data$age, na.rm = TRUE)
clean_data_4

# ------------------------------
# Summary of Data Structures:
# ------------------------------

#   - Vectors: simple sequences of same data type
#   - Lists: mix of different data types
#   - Matrices: numeric tables
#   - Data Frames: mixed-type tables 

#--------------------
# 3. Functions in R
#--------------------
# Functions let us wrap code into reusable blocks.

# function is  reusable block of code 
# Why use functions?
#   - Avoid repetition
#   - Organize and simplify code
#   - Reuse across projects (save it for later use)
#   - Share with others

# A function in R has 4 key parts:
#   1. Name         -> the name you give to the function
#   2. Arguments    -> the inputs you provide to the function
#   3. Body         -> the set of operations the function performs
#   4. Return Value -> the output the function gives back

# Example: A function to calculate Body Mass Index (BMI)
# 1. Function Name: calculate_BMI
# 2. Arguments: weight (in kg), height (in meters)
# 3. Body: performs BMI calculation e.g   # Formula: BMI = weight / (height^2)
# 4. Return Value: the BMI value

# create a function to calculate BMI
# first argument = weight
# second argument = height
calculate_BMI <- function(weight, height, age){
  # operation we want to perform
  bmi <- weight/(height^2)
  return(bmi)
}

# call the function
calculate_BMI(weight = 60, height = 1.75)
calculate_BMI(weight = weight, height = height)
calculate_BMI(60)  

# ----------------------------
# Lazy evaluation in R
# ----------------------------
# If your function has three arguments, but the body only uses two,
# R does not force you to supply the third argument for the calculation.
# Example: 'age' is defined as an argument, but not used in the formula.

# if your function is expecting two arguments than one two would show up
# Even though 'age' exists as an argument, it is ignored because it is not used
calculate_BMI(60, 1.65)

#---------
# Summary:
#---------
# Functions help us package logic once and apply it to different inputs.

# ----------------------------------
# 4. Automating Workflows with for-Loop
# ----------------------------------
# Suppose you have multiple datasets and you want to:
#   - import them,
#   - check missing values,
#   - clean columns,
#   - compute BMI,
#   - and save results.
#
# Instead of repeating steps for each file, we use loops.

# -----------------------
# Typical loop workflow:
# -----------------------
# 1. Define input and output folders 
input_dir <- "raw_data"
output_dir <- "results"
# create output folder if not already exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
# the results folder will be created if already not exist

# 2. list which files to process
files_to_process <- c("BMI_data_1.csv", "BMI_data_2.csv") # name of files

# 3. prepare empty list to store results in R
result_list <- list()

# 4. For each file: 
#   - import data 
#   - handle NA values
#   - calculate BMI using calculate-BMI function
#   - save results (both as CSV and inside R list)
# 

#for (file_names in files_to_process) {
#  cat("\nProcessing:", file_names, "\n")
#  input_file_path <- file.path(input_dir, file_names)
  # import dataset
#  data <- read.csv(input_file_path, header = TRUE)
#  cat("file imported. Checking for missing values.. \n")
  # handling missing values
#  if("height" %in% names(data)){
#    missing_count <- sum(is.na(data$height))
 #   
 #   cat("Missing values in 'height':", missing_count, "\n")
#    data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
  }
 # if("weight" %in% names(data)){
  #  missing_count <- sum(is.na(data$weight))
    
   # cat("Missing values in 'weight':", missing_count, "\n")
   # data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
  }
  # calculate BMI
  #data$bmi <- calculate_BMI(data$weight, data$height)
 # cat("BMI has been calculated successfully.\n")
  
  # save results in R
  #result_list[[file_names]] <- data
  # save results in results folder
  #output_file_path <- file.path(output_dir, paste0("BMI_results", file_names))
  #write.csv(data, output_file_path, row.names = FALSE)
 # cat("Results saved to:", output_file_path, "\n")
#}

for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  
  # Build the full path
  input_file_path <- file.path(input_dir, file_name)
  
  # Check if the file exists
  if (!file.exists(input_file_path)) {
    cat("WARNING: File not found ->", input_file_path, "\n")
    next  # Skip this file and move on
  }
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported successfully.\n")
  
  # Handle missing values for height
  if ("height" %in% names(data)) {
    missing_height <- sum(is.na(data$height))
    cat("Missing 'height' values:", missing_height, "\n")
    data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
  }
  
  # Handle missing values for weight
  if ("weight" %in% names(data)) {
    missing_weight <- sum(is.na(data$weight))
    cat("Missing 'weight' values:", missing_weight, "\n")
    data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
  }
  
  # Use your pre-defined calculate_BMI function
  data$bmi <- calculate_BMI(data$weight, data$height)
  cat("BMI calculated successfully.\n")
  
  # Save results inside R
  result_list[[file_name]] <- data
  
  # Save results to the results folder
  output_file_path <- file.path(output_dir, paste0("BMI_results_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
}

# the loop repeats until all files are processed. 

results_1 <- result_list[[1]]
results_2 <- result_list[[2]]
View(results_1)
View(results_2)

# --------
# Summary:
# --------
# Loops automate repetitive work 
# making your workflow faster 
# consistent, and reproducible

###########################################################################
# --------------------------
# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise

# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv

# Each file contains three columns: 
# Gene_Id	
# padj	
# logFC

####################################################################
# 1. Define Input and Output Folders
input_dir  <- "raw_data"
output_dir <- "Results"  # Match assignment capitalization

# Check input folder
if (!dir.exists(input_dir)) {
  stop("ERROR: The 'raw_data' folder was not found.\nPlease make sure it exists next to your script.")
}

# Create Results folder if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created Results folder at:", output_dir, "\n")
}

# 2. List the Files to Process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

available_files <- list.files(input_dir)
cat("Files found in raw_data:\n")
print(available_files)

missing_files <- setdiff(files_to_process, available_files)
if (length(missing_files) > 0) {
  stop("ERROR: The following files are missing from raw_data:\n",
       paste(missing_files, collapse = "\n"))
}

# 3. Define classify_gene() function
# Returns exact labels requested by the assignment.
classify_gene <- function(logFC, padj) {
  # Replace missing padj with 1 per instructions
  if (is.na(padj)) padj <- 1
  
  if (logFC > 1 && padj < 0.05) {
    "Upregulated"
  } else if (logFC < -1 && padj < 0.05) {
    "Downregulated"
  } else {
    "Not_Significant"
  }
}

# Quick tests
stopifnot(classify_gene( 2, 0.01) == "Upregulated")
stopifnot(classify_gene(-2, 0.01) == "Downregulated")
stopifnot(classify_gene( 0.5, 0.10) == "Not_Significant")

# 4. Prepare an Empty List to Store Results
result_list <- list()

# 5. Process the Datasets
for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  input_file_path <- file.path(input_dir, file_name)
  
  if (!file.exists(input_file_path)) {
    cat("WARNING: File not found", input_file_path, "\n")
    next
  }
  
  # Import data
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported successfully!\n")
  
  # Basic column check (as described in the assignment)
  required_cols <- c("Gene_Id", "padj", "logFC")
  if (!all(required_cols %in% names(data))) {
    stop("ERROR: One or more required columns are missing. Expected: ",
         paste(required_cols, collapse = ", "))
  }
  
  # Classify genes (NA handling happens inside classify_gene)
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  cat("Gene classification completed.\n")
  
  # Save processed dataset
  result_list[[file_name]] <- data
  output_file_path <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  # Print summary counts with table()
  cat("\nSummary for", file_name, ":\n")
  status_tab <- table(data$status)
  print(status_tab)
  
  # Optional: combined Significant count (Up + Down)
  significant_count <- sum(data$status %in% c("Upregulated", "Downregulated"))
  cat("Significant (Up + Down):", significant_count, "\n")
}

# 6. Simple summary table across datasets
summary_table <- data.frame(
  Dataset = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Not_Significant = integer(),
  Significant = integer(),
  stringsAsFactors = FALSE
)

for (file_name in files_to_process) {
  data <- result_list[[file_name]]
  counts <- table(data$status)
  
  up      <- ifelse("Upregulated"     %in% names(counts), counts["Upregulated"], 0)
  down    <- ifelse("Downregulated"   %in% names(counts), counts["Downregulated"], 0)
  not_sig <- ifelse("Not_Significant" %in% names(counts), counts["Not_Significant"], 0)
  
  summary_table <- rbind(
    summary_table,
    data.frame(
      Dataset = file_name,
      Upregulated = as.integer(up),
      Downregulated = as.integer(down),
      Not_Significant = as.integer(not_sig),
      Significant = as.integer(up + down)
    )
  )
}

print(summary_table)

# Check how many have padj < 0.05 and |logFC| > 1 separately
# for (file_name in files_to_process) {
  # d <- result_list[[file_name]]
  # cat("\n", file_name, "\n")
  # cat("padj < 0.05: ", sum(d$padj < 0.05, na.rm = TRUE), "\n")
  # cat("|logFC| > 1 : ", sum(abs(d$logFC) > 1, na.rm = TRUE), "\n")
  #cat("Both (significant): ",
     # sum(d$padj < 0.05 & abs(d$logFC) > 1, na.rm = TRUE), "\n")
#}

###### DONE!

cat("\n--- Assignment 2 Completed Successfully ---\n")
cat("Processed files saved in:", output_dir, "\n")
cat("Summary table printed above.\n")
