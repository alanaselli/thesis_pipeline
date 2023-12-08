library(AlphaSimR)
library(dplyr)

# ---- new_subdir ----
# Create subfolders

new_subdir = function(subdir_names){  # Provide sub directory names in a vector
    for (subdir in subdir_names) {      # For each name provided
        if (!dir.exists(subdir)) {        # Check if folder already exists
            dir.create(subdir)              # If not, create folder
            print(paste("Directory", subdir, "successfully created"))
        } else{
            print("Directory already exists")
        }
    }
}