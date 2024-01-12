mate_selection_GEBV_Fg = function(pop = recentPop,
                                   year,
                                   scenario_folder="scenario_02/",
                                   ped_ROH="recent.ped",
                                   pre_selection_males_porc = 0.8, # percentage to remove
                                   Fg_percentage = 0.2, # percentage to keep
                                   male_groups = c(40,10),
                                   female_groups = c(300,250,200,150,100),
                                   nMatings = 1000,
                                   nMatings_per_sire=20,
                                   append = TRUE
){
    if (isTRUE(append)) {
        col.names = FALSE
    } else {col.names = TRUE}
    
    fakePed = makeFakeHaplos(dirToSave = scenario_folder,
                             pop=pop)
    fakePed$sire = sub("\\_.*", "", fakePed$sire)
    fakePed$dam = sub("\\_.*", "", fakePed$dam)
    fakePed[,c(4:7)] = 0
    fakePed$gen = year
    
    system(command=paste0("cp ",scenario_folder,"pedigree.txt ",
                          scenario_folder,"pedigree_with_fakes.txt"))
    
    # Append fake progenies in the pedigree
    write.table(fakePed,paste0(scenario_folder,"pedigree_with_fakes.txt"),
                col.names = FALSE, append = TRUE,
                row.names = F, quote = F)
    
    # Convert raw progeny file to plink format - only in 3rd scenario
    system(command=paste0("scripts/raw_to_plink.sh ",
                          scenario_folder, "fakeHaplos_raw.txt ",
                          scenario_folder, "fakeHaplos.ped"))
    
    
    # Generate files for BLUPF90
    system(command=paste0("./scripts/extract_ped.sh ",
                          scenario_folder,
                          "pedigree_with_fakes.txt 05_BLUPF90/dat1.txt 05_BLUPF90/ped1.txt 10"))
    
    # Fit RR-BLUP model for genomic predictions (candidates)
    BLUP = GBLUP_AlphaSimR(pop)
    
    # For BLUPF90 I need all generations (recent)
    # Merge genomic data from last generations with fake progeny (for BLUPF90)
    system(command=paste0("cp ",scenario_folder,"recent.ped ",
                          scenario_folder,"genotype_with_fakes.ped"))
    system(command=paste0('plink --cow --ped ',
                          scenario_folder,'genotype_with_fakes.ped --map ',
                          scenario_folder,'new_map.map --merge ',
                          scenario_folder,'fakeHaplos.ped ',
                          scenario_folder,'new_map.map --make-bed --recode --allow-no-sex --out ',
                          scenario_folder,'genotype_with_fakes'))
    
    # Convert to BLUPF90 format
    system(command = paste0("scripts/prepare_snp_file.sh ",
                            scenario_folder,"genotype_with_fakes.ped ",
                            scenario_folder,"new_map.map 05_BLUPF90/snp_file.txt"))

    # Run BLUPF90 with genomic data
    runBLUPF90(param_card="renum_genomic.txt")
    
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
    
    # Read Fg
    Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                    col.names = c("ID","Fg"),
                    colClasses = c("character","numeric"))
    
    # ROH (in this scenario, only candidates - for evaluation purposes)
    FROH = ROH_analyses(pop = pop,
                        generation = year,
                        ped = paste0(scenario_folder,ped_ROH),
                        scenario = "sc_02",
                        map="01_genotypes/new_map.map",
                        save_to="03_ROH/")
    
    # Merge dataframes (candidates)
    candidates_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                            by.x = "ID", by.y = "original_id",
                            all.x = TRUE)
    candidates_data = merge(candidates_data, Fped,
                            by = "ID", all.x = TRUE)
    candidates_data = merge(candidates_data, Fg,
                            by = "ID", all.x = TRUE)
    candidates_data = merge(candidates_data, FROH[,c(1,4)],
                            by.x = "ID", by.y = "id", all.x = TRUE)
    
    rm(BLUP,FROH)
    
    cli_alert_info(paste0("\nTotal number of candidates: ",nrow(candidates_data)))
    
    # Save metrics to scenario file
    write.table(candidates_data, paste0(scenario_folder,"candidates_metrics.txt"),
                append = append, col.names = col.names,
                quote = F, row.names = F)
    
    # Define Fg_threshold
    Fg_threshold = mean(candidates_data$Fg) + mean(candidates_data$Fg)*Fg_percentage
    
    # ROH progeny
    # FROH = ROH_analyses(pop = (scenario_folder),
    #                     generation = year,
    #                     ped = paste0(scenario_folder,ped_ROH),
    #                     map="01_genotypes/new_map.map",
    #                     save_to="03_ROH/")
    
    # Merge dataframes (progeny)
    df = merge(fakePed[,c(1:3,8)], BLUPF90_EBVs[,c("original_id","solution")],
               by.x = "ID", by.y = "original_id",
               all.x = TRUE)
    df = merge(df, Fped, by = "ID", all.x = TRUE)
    df = merge(df, Fg, by = "ID", all.x = TRUE)
    # df = merge(df, FROH[,c(1,4)], by = "ID", all.x = TRUE)
    
    rm(fakePed, BLUPF90_EBVs, Fped, Fg)
    
    write.table(df,paste0(scenario_folder,"fakeProgeny.txt"),
                append = append, col.names = col.names,
                quote = F, row.names = F)
    df2 = df
    df = df2
    
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
    
    # Summarise max Fg and mean GEBV for each mating
    # 1st select X% of matings with the lowest Fg
    # 2nd rank matings by EBV and Fg
    
    # df = df %>%
    #     slice_min(Fg, prop = Fg_percentage) %>%  # take only X% of matings with lowest Fg
    #     arrange(desc(round(solution,2)),Fg)
    
    df = df %>% 
        group_by(sire, dam) %>% 
        summarise(max_Fg = max(Fg), 
                  mean_GEBV = mean(solution),
                  dam_group = mean(dam_group),
                  sire_group = mean(sire_group)) %>% 
        slice_min(max_Fg, prop = Fg_percentage) %>% 
        arrange(desc(round(mean_GEBV,2)),max_Fg)
    
    df = as.data.frame(df)
    
    cli_alert_info(paste0("Number of female candidates after Fg filter: ",
                          length(unique(df$dam))))
    cli_alert_info(paste0("Number of male candidates after Fg filter: ",
                          length(unique(df$sire))))
    
    df3 = df
    df = df3
    # Select matings
    matings = df[0,]
    
    group_counts = data.frame(group=c(1:length(male_groups),
                                      1:length(female_groups)), 
                              sex=c(rep("M",length(male_groups)),
                                    rep("F",length(female_groups))), 
                              count=0,
                              max=c(male_groups,female_groups))
    
    max_matings_per_sire = 20
    min_matings_per_sire = 5
    sire_use = data.frame(sire=character(0),
                          group=integer(0),
                          use=integer(0))
    
    # ----
    
    best_dams = df[!duplicated(df$dam),]
    best_dams = best_dams %>% 
        mutate(rank = dense_rank(desc(mean_GEBV))) %>% 
        select(dam, dam_group, rank)
    
    best_dams = merge(best_dams, n_matings_each_dam)
    
    # ----
    
    # Rank matings
    df$rank = 1:nrow(df)
    df = df[order(df$rank),]
    
    df_backup = df
    
    while (nrow(matings)<nMatings) {
        # N matings per dam
        n_matings_each_dam = df %>% 
            group_by(dam) %>% 
            summarise(n=n()) %>% 
            select(dam, n)
        
        df = df %>% 
            select_if(!names(.) %in% c('n'))
        
        df = merge(df, n_matings_each_dam)
        
        # Select matings with n == 1 and within the best matings (n_matings_left)
        n_matings_left = nMatings - nrow(matings)
        selected_matings = df %>% 
            filter(n==min(df$n[df$rank <= n_matings_left]) & rank <= n_matings_left)
        
        matings = rbind(matings, selected_matings)
        
        # Remove other matings for the selected dams
        df = df[!df$dam %in% selected_matings$dam,]
        
        # Add counts to male groups
        selected_sire = data.frame(sire=selected_matings$sire,
                                   group=selected_matings$sire_group)
        selected_sire = selected_sire %>% 
            group_by(sire,group) %>% 
            summarise(n=n())
        sire_use = rbind(sire_use,selected_sire)
        
        male_counts = table(selected_sire$group)
        male_group = group_counts$group %in% names(male_counts) & group_counts$sex == "M"
        group_counts$count[male_group] = group_counts$count[male_group] + male_counts
        
        # Add counts to female groups
        dam_counts = table(selected_matings$dam_group)
        female_group = group_counts$group %in% names(dam_counts) & group_counts$sex == "F"
        group_counts$count[female_group] = group_counts$count[female_group] + dam_counts
        
        # Check constraints
        # Max group counts
        for (sex in c("M","F")) {
            if (sex == "M") {
                groups = 1:length(male_groups)
            } else {groups = 1:length(female_groups)}
            for (group in groups) {
                g_count = group_counts$count[group_counts$sex == sex & group_counts$group == group]
                g_max = group_counts$max[group_counts$sex == sex & group_counts$group == group]
                if (g_count >= g_max) {
                    if (sex == "M") {
                        df = df[!df$sire_group == group,]
                    } else {df = df[!df$dam_group == group,]}
                }
            }
        }
        
        # Max sire use
    }
    
    # ----
    
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
        for (sex in c("M","F")) {
            if (sex == "M") {
                groups = 1:length(male_groups)
            } else {groups = 1:length(female_groups)}
            for (group in groups) {
                g_count = group_counts$count[group_counts$sex == sex & group_counts$group == group]
                g_max = group_counts$max[group_counts$sex == sex & group_counts$group == group]
                if (g_count >= g_max) {
                    if (sex == "M") {
                        df = df[!df$sire_group == group,]
                    } else {df = df[!df$dam_group == group,]}
                }
            }
        }
    }
    
    # Check for repeated females
    if (length(unique(matings$dam)) < sum(female_groups)) {
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
