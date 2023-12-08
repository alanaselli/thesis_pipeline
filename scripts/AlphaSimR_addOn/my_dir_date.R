library(AlphaSimR)
library(dplyr)

# ---- new_dir_date ----
# Create a new directory with today's date

new_dir_date = function(prefix = "sim"){
    # Create directory's name
    newdir <- paste(prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
    
    # Create new directory
    if (!dir.exists(newdir)) {        # Check if folder already exists
        dir.create(newdir)              # If not, create folder
        setwd(newdir)                   # Set to new directory
        print(paste("Directory", newdir, "successfully created"))
    } else{
        print("Directory already exists")
    }
}