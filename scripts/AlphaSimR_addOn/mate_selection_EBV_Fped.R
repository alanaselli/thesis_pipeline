mate_selection_EBV_Fped = function(pop,
                                   year,
                                   scenario_folder,
                                   pre_selection_males_porc = 0.8, # percentage to remove
                                   Fped_percentage = 0.2, # percentage to keep
                                   male_groups = c(40,10),
                                   female_groups = c(300,250,200,150,100),
                                   nMatings = 1000,
                                   nMatings_per_sire=20
                                   ){
    fakePed = makeFakePed(pop)
    
    fakePed[,c(4:7)] = 0
    fakePed$gen = year
    
    system(command=paste0("cp ",scenario_folder,"pedigree.txt ",
                          scenario_folder,"pedigree_with_fakes.txt"))
    
    write.table(fakePed,paste0(scenario_folder,"pedigree_with_fakes.txt"),
                col.names = F, append = T,
                row.names = F, quote = F)
    
    system(command=paste0("./scripts/extract_ped.sh ",
                           scenario_folder,
                           "pedigree_with_fakes.txt 05_BLUPF90/dat1.txt 05_BLUPF90/ped1.txt 10"))
    
    # Run BLUPF90 without genomic data
    runBLUPF90(param_card="renum.txt")
    
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
    
    # Run BLUPF90 with genomic data
    runBLUPF90(param_card="renum_genomic.txt")
    
    # Fit RR-BLUP model for genomic predictions (candidates)
    BLUP = GBLUP_AlphaSimR(pop)
    
    # Read Fg
    Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                    col.names = c("ID","Fg"),
                    colClasses = c("character","numeric"))
    
    # Merge dataframes (candidates)
    candidates_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                            by.x = "ID", by.y = "original_id",
                            all.x = TRUE)
    candidates_data = merge(candidates_data, Fped,
                            by = "ID", all.x = TRUE)
    candidates_data = merge(candidates_data, Fg,
                            by = "ID", all.x = TRUE)
    
    rm(BLUP,Fg)
    
    cli_alert_info(paste0("\nTotal number of candidates: ",nrow(candidates_data)))
    
    # Save metrics to scenario file
    write.table(candidates_data, paste0(scenario_folder,"candidates_metrics.txt"),
                append = F, col.names = T,
                quote = F, row.names = F)
    
    # Merge dataframes (progeny)
    df = merge(fakePed[,c(1:3,8)], BLUPF90_EBVs[,c("original_id","solution")],
                         by.x = "ID", by.y = "original_id",
                         all.x = TRUE)
    df = merge(df, Fped, by = "ID", all.x = TRUE)
    
    rm(fakePed, BLUPF90_EBVs, Fped)
    
    write.table(df,paste0(scenario_folder,"fakeProgeny.txt"),
                append = F, col.names = T,
                quote = F, row.names = F)
    
    # Pre-selection of male candidates
    best_EBVs = candidates_data %>% 
        filter(sex == "M") %>% 
        slice_max(solution, prop = pre_selection_males_porc) %>%
        select(ID)
    
    df = df %>% 
        filter(sire %in% best_EBVs$ID)
    
    rm(best_EBVs)
    
    # Assign groups
    candidates_year = data.frame(id=pop@id, 
                                 year=unname(unlist(pop@misc))) %>% 
        mutate(group = dense_rank(desc(year))) %>% 
        select(!year)
    
    df = merge(df, candidates_year, 
                by.x = "dam", by.y = "id") %>% 
                rename(dam_group = group)
    df = merge(df, candidates_year,
                by.x = "sire", by.y = "id") %>% 
                rename(sire_group = group)
    rm(candidates_year)
    
    # 1st remove matings with high Fped
    # 2nd rank matings by EBV and Fped
    df = df %>% 
        slice_min(Fped, prop = Fped_percentage) %>%  # take only X% of matings with lowest Fped
        arrange(desc(round(solution,2)),Fped)
    
    cli_alert_info(paste0("Number of female candidates after Fped filter: ",
                          length(unique(df$dam))))
    cli_alert_info(paste0("Number of male candidates after Fped filter: ",
                          length(unique(df$sire))))
    
    # Select matings
    matings = df[0,]
    
    group_counts = data.frame(group=c(1:length(male_groups),
                                      1:length(female_groups)), 
                              sex=c(rep("M",length(male_groups)),
                                    rep("F",length(female_groups))), 
                              count=0,
                              max=c(male_groups,female_groups))
    
    while (nrow(matings)<nMatings) {
        # Select the first sire
        sel_sire = df[1,'sire']
        
        # Select the first N matings of that sire
        sire_matings = df %>% 
            filter(sire == sel_sire) %>% 
            slice_head(n = nMatings_per_sire)
        
        # Remove dams from selection pool
        df = df[!df$dam %in% sire_matings$dam,]
        
        # Remove sire from selection pool
        df = df[!df$sire == sel_sire,]
        
        # Add count to male group
        male_group = group_counts$sex == "M" & group_counts$group == sire_matings[1,'sire_group']
        group_counts$count[male_group] = group_counts$count[male_group]+1
        
        # Add counts to female groups
        dam_counts = table(sire_matings$dam_group)
        female_group = group_counts$group %in% names(dam_counts) & group_counts$sex == "F"
        group_counts$count[female_group] = group_counts$count[female_group] + dam_counts
        
        matings = rbind(matings, sire_matings)
        
        # Check constraints
        # Max group counts
        for (s in c("M","F")) {
            if (s == "M") {
                groups = 1:length(male_groups)
                } else {groups = 1:length(female_groups)}
            for (g in groups) {
                g_count = group_counts$count[group_counts$sex == s & group_counts$group == g]
                g_max = group_counts$max[group_counts$sex == s & group_counts$group == g]
                if (g_count > g_max) {
                    if (s == "M") {
                        df = df[!df$sire_group == g,]
                    } else {df = df[!df$dam_group == g,]}
                }
            }
        }
    }
    
    # Check for repeated females
    if (length(unique(matings$dam)) < 1000) {
        cli_alert_danger("\nRepeated females!\n")
    }
    
    # Breed new generation
    new_gen = makeCross(pop, as.matrix(matings[,c(2,1)]))
    
    new_gen = setMisc(x = new_gen,
                      node = "yearOfBirth",
                      value = year)
    
    # Merge new generation and older candidates
    older_candidates = c(unique(matings$sire[!matings$sire_group == male_groups[length(male_groups)]]),
                         unique(matings$dam[!matings$dam_group == female_groups[length(female_groups)]]))
    
    candidates_data = candidates_data %>% 
        filter(ID %in% older_candidates)
    
    candidates = mergePops(list(new_gen,pop[pop@id %in% older_candidates]))
    
    return(candidates)
}