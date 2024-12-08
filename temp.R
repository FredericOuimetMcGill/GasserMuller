# Set the working directory to the current directory
setwd(".")

# Specify the file names
file_list <- c("raw_ISE_MC_results_456_1_20.csv",
               "raw_ISE_MC_results_456_21_40.csv",
               "raw_ISE_MC_results_456_41_60.csv",
               "raw_ISE_MC_results_456_61_80.csv",
               "raw_ISE_MC_results_456_81_100.csv")

# Read the first file to get the header and initialize the combined dataframe
combined_data <- read.csv(file_list[1], header = TRUE, stringsAsFactors = FALSE)

# For each subsequent file, read the data and append it to combined_data
for (i in 2:length(file_list)) {
  temp_data <- read.csv(file_list[i], header = TRUE, stringsAsFactors = FALSE)
  combined_data <- rbind(combined_data, temp_data)
}

# Write out the combined file
write.csv(combined_data, "raw_ISE_MC_results_456.csv", row.names = FALSE)
