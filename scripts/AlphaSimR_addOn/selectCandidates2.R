library(AlphaSimR)
library(dplyr)

# ---- selectCandidates2 ----
selectCandidates2 = function(pop, 
                             genotypes, # runBLUPF90
                             min_gen,   # runBLUPF90
                             file_name, 
                             append=TRUE, 
                             method=1, # 1=EBV+F; 2=GEBV+Fg; 3=GEBV+Froh
                             top_ebv,  # c(nMales, nFemales)
                             # int or fraction
                             max_F = NULL, # Max inbreeding (num)
                             gen,
                             OriginalRecords,
                             pathScenario
){
    # Select candidates based on fake progeny (pedigree or genomic)
    
    # Append or not
    if (isTRUE(append)) {
        colnames_par = FALSE
    } else {colnames_par = TRUE}
    
    # Create fake data
    if (method == 1){
        fake_pedigree = makeFakePed(pop)
    } else {
        fake_pedigree = makeFakeHaplos(pathScenario, pop)
        
        # Append fake haplos to candidates' genotypes
        system(paste0("./append_snp_file.sh ",
                      pathScenario,genotypes," ",
                      pathScenario,"fakeHaplos.txt"))
    }
    
    # Add data to scenario/pedigree.txt
    write.table(fake_pedigree, paste0(pathScenario,"pedigree.txt"), 
                append = T, col.names = F, 
                quote = F, row.names = F)
    
    # Add generation number
    fake_pedigree$generation = gen
    
    # Run BLUPF90
    if (method==1){ 
        # Run without genomic information
        runBLUPF90(ped_name,
                   pedigree = paste0(pathScenario,"pedigree.txt"),
                   param_card="renum.txt",
                   min_gen)
    } else {
        # Run with genomic information
        runBLUPF90(ped_name,
                   pedigree = paste0(pathScenario,"pedigree.txt"),
                   param_card="renum_genomic.txt",
                   min_gen)
    }
    
    # Read EBVs
    BLUPF90_EBVs = read.table("05_BLUPF90/solutions.orig", header = T,
                              colClasses = c("integer","integer",
                                             "integer","character",
                                             "numeric"))
    
    # Read Fped
    Fped = read.table("05_BLUPF90/renf90.inb",
                      col.names = c("ID","Fped","delete"),
                      colClasses = c("character","numeric","numeric"))
    Fped = Fped[,c(1,2)]
    
    if (method == 1) {
        # Run BLUP again to obtain Fg
        runBLUPF90(ped_name,
                   param_card="renum_genomic.txt",
                   min_gen)
    }
    
    # Read Fg
    Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                    col.names = c("ID","Fg"),
                    colClasses = c("character","numeric")
    )  
    
    # Fit RR-BLUP model for genomic predictions (candidates)
    BLUP = GBLUP_AlphaSimR(pop)
    
    # Merge dataframes (candidates)
    candidates_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                            by.x = "ID", by.y = "original_id",
                            all.x = TRUE)
    candidates_data = merge(candidates_data, Fped,
                            by = "ID", all.x = TRUE)
    candidates_data = merge(candidates_data, Fg,
                            by = "ID", all.x = TRUE)
    
    # Save metrics to scenario file
    write.table(candidates_data, paste0(pathScenario,"candidates_metrics.txt"),
                append = append, col.names = colnames_par,
                quote = F, row.names = F)
    
    # Merge dataframes (progeny)
    progeny_data = merge(fake_pedigree, BLUPF90_EBVs[,c("original_id","solution")],
                         by.x = "ID", by.y = "original_id",
                         all.x = TRUE)
    progeny_data = merge(progeny_data, Fped,
                         by = "ID", all.x = TRUE)
    progeny_data = merge(progeny_data, Fg,
                         by = "ID", all.x = TRUE)
    
    all_data = rabind(candidates_data,progeny_data)
    
    # Add data to scenario/all_data.txt
    write.table(all_data, paste0(pathScenario,"all_data.txt"), 
                append = append, col.names = colnames_par,
                quote = F, row.names = F)
    
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