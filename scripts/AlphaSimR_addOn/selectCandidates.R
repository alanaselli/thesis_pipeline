library(AlphaSimR)
library(dplyr)

# ---- selectCandidates ----
selectCandidates = function(pop, 
                            file_name, 
                            append=TRUE, 
                            method=1, # 1=EBV; 2=GEBV+F; 3=GEBV+Fg
                            top_ebv,  # c(nMales, nFemales)
                            # int or fraction
                            max_F = NULL # Max inbreeding (num)
){
    
    # Fit RR-BLUP model for genomic predictions
    BLUP = GBLUP_AlphaSimR(pop)
    
    # Read EBVs
    BLUPF90_EBVs = read.table("05_BLUPF90/solutions.orig", header = T,
                              colClasses = c("integer","integer",
                                             "integer","character",
                                             "numeric")
    )
    
    # Read Fped
    Fped = read.table("05_BLUPF90/renf90.inb",
                      col.names = c("ID","Fped","delete"),
                      colClasses = c("character","numeric","numeric")
    )
    Fped = Fped[,c(1,2)]
    
    # Read Fg
    Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                    col.names = c("ID","Fg"),
                    colClasses = c("character","numeric")
    )
    
    # Merge dataframes
    merged_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                        by.x = "ID", by.y = "original_id",
                        how = "left")
    merged_data = merge(merged_data, Fped,
                        by = "ID", how = "left")
    merged_data = merge(merged_data, Fg,
                        by = "ID", how = "left")
    
    # Save metrics to scenario file
    if (isTRUE(append)) {
        write.table(merged_data, file_name, sep = ",", 
                    append = T, col.names = F,
                    quote = F, row.names = F)
    } else {
        write.table(merged_data, file_name, sep = ",", 
                    append = F, col.names = T,
                    quote = F, row.names = F)
    }
    
    # Check correlations
    cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and GV: ", 
                          round(cor(merged_data$EBV, merged_data$GV),2)))
    cli_alert_info(paste0("\nCorrelation BLUPF90 EBV and GV: ", 
                          round(cor(merged_data$solution, merged_data$GV),2)))
    cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and BLUPF90 EBV: ", 
                          round(cor(merged_data$EBV, merged_data$solution),2)))
    cli_alert_info(paste0("\nCorrelation Fped and Fg: ", 
                          round(cor(merged_data$Fped, merged_data$Fg),2)))
    
    # If a fraction was passed for EBV, calculate the number of
    # animals to select
    # males
    if (top_ebv[1]<1) {
        nMales = nrow(merged_data[merged_data$sex == "M",])
        selectMales = ceiling(nMales * top_ebv[1]) # round up
    } else {selectMales = top_ebv[1]}
    
    # females
    if (top_ebv[2]<1) {
        nFemales = nrow(merged_data[merged_data$sex == "F",])
        selectFemales = ceiling(nFemales * top_ebv[2]) # round up
    } else {selectFemales = top_ebv[2]}
    
    # Method 1
    # Select animals based on BLUPF90 GEBV
    if (method == 1) {
        males = merged_data %>%
            select(ID, sex, solution) %>% 
            filter(sex == "M") %>%
            arrange(desc(solution)) %>% 
            slice_head(n=selectMales)
        
        females = merged_data %>%
            select(ID, sex, solution) %>%
            filter(sex == "F") %>%
            arrange(desc(solution)) %>% 
            slice_head(n=selectFemales)
        
    } else if (method == 2) {
        # Select animals based on BLUPF90 EBV and Fped
        males = merged_data %>% 
            select(ID, sex, solution, Fped) %>%
            filter(sex == "M") %>% 
            filter(Fped <= max_F) %>% 
            arrange(desc(solution)) %>% 
            slice_head(n=selectMales)
        
        females = merged_data %>% 
            select(ID, sex, solution, Fped) %>%
            filter(sex == "F") %>% 
            filter(Fped <= max_F) %>% 
            arrange(desc(solution)) %>% 
            slice_head(n=selectFemales)
        
    } else if (method == 3) {
        # Select animals based on BLUPF90 EBV and Fg
        males = merged_data %>% 
            select(ID, sex, solution, Fg) %>%
            filter(sex == "M") %>% 
            filter(Fg <= max_F) %>% 
            arrange(desc(solution)) %>% 
            slice_head(n=selectMales)
        
        females = merged_data %>% 
            select(ID, sex, solution, Fg) %>%
            filter(sex == "F") %>% 
            filter(Fg <= max_F) %>% 
            arrange(desc(solution)) %>% 
            slice_head(n=selectFemales)
        
    } else {cli_alert_danger("\nSelection method not identified!\n")}
    
    male_parents = pop[pop@id %in% as.character(males$ID)]
    female_parents = pop[pop@id %in% as.character(females$ID)]
    
    # Check the number of selected animals
    if (male_parents@nInd != selectMales) {
        cli_alert_warning("The number of selected males is different than the required.")
    }
    
    if (female_parents@nInd != selectFemales){
        cli_alert_warning("The number of selected females is different than the required.")
    }
    
    cli_alert_info(paste0("Mean EBV of the population: ",
                          round(mean(merged_data$solution), 2)))
    cli_alert_info(paste0("Mean EBV of selected males: ",
                          round(mean(males$solution), 2)))
    cli_alert_info(paste0("Mean EBV of selected females: ",
                          round(mean(females$solution), 2)))
    
    parents = c(male_parents, female_parents)
    
    return(parents)
}