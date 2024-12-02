# install.packages("jsonlite")
library(jsonlite)
# Specify the path to your JSON file
json_file <- "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/all_dart_reports/Report-DLan23-8067/metadata.json"
json_file <- "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/all_dart_reports/Report-DLan22-7500/Report_DLan22_7500_7160_6495_6370/metadata.json"

# Read the JSON data into an R variable
data <- fromJSON(json_file)

# Now, 'data' contains your JSON data as an R object

targetids <- data[["targetslist"]]

table(targetids %in% colnames(gt))
"2049605" %in% colnames(gt)

z <- data.frame(data[["datafiles"]][["Report_DLan23-8067_4_moreOrders_SNP_2.csv"]][["sampleheaders"]][["headers"]]) %>% t()

targetid_nsw <- z[,c(5,7)] 
rownames(targetid_nsw) <- NULL
colnames(targetid_nsw) <- c("sample","targetid")
targetid_nsw <- as.data.frame(targetid_nsw)
targetid_nsw$sample

table(targetid_nsw$targetid %in% colnames(gta))



targetids <- data[["targetslist"]]
table(targetids %in% colnames(gt))


#########

# Directory path containing CSV files
directory_path <- "/Users/eilishmcmaster/Documents/LantCama/LantCama/dart_raw/lantana_camara_meta"

# List all CSV files in the directory
csv_files <- list.files(directory_path, pattern = ".csv", full.names = TRUE)

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through CSV files and read them into data frames
for (file in csv_files) {
  data_frames[[file]] <- read.csv(file, header = TRUE)
}

# Combine all data frames into a single data frame using do.call and rbind
combined_data <- do.call(rbind, data_frames)
rownames(combined_data) <- NULL
# Print the combined data frame
print(head(combined_data))

table(combined_data$targetid %in%  colnames(gta))

c2 <- merge(m2, combined_data, by.x="sample", by.y="genotype", all.x=TRUE)

write.csv(c2, "/Users/eilishmcmaster/Documents/LantCama/LantCama/meta_targetid.csv")
