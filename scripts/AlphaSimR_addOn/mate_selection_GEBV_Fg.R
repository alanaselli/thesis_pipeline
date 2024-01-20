mate_selection_GEBV_Fg = function(pop = recentPop,
                                   year,
                                   scenario_folder="scenario_02/",
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
    
    # Append fake progeny to the pedigree
    system(command=paste0("cp ",scenario_folder,"pedigree.txt ",
                          scenario_folder,"pedigree_with_fakes.txt"))
    
    write.table(fakePed,paste0(scenario_folder,"pedigree_with_fakes.txt"),
                col.names = FALSE, append = TRUE,
                row.names = F, quote = F)
    
    # Convert raw progeny file to plink format
    system(command=paste0("scripts/raw_to_plink.sh ",
                          scenario_folder, "fakeHaplos_raw.txt ",
                          scenario_folder, "fakeHaplos.ped"))
    
    # Generate files for BLUPF90
    system(command=paste0("./scripts/extract_ped.sh ",
                          scenario_folder,
                          "pedigree_with_fakes.txt 05_BLUPF90/dat1.txt 05_BLUPF90/ped1.txt 10"))
    
    cli_alert_info("\nStart evaluation of candidates.\n")
    
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
    FROH = ROH_analyses(pop_name = "scenario_02",
                        generation = year,
                        ped = paste0(scenario_folder,"candidates.ped"),
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
    candidates_data$gen = year
    
    rm(BLUP,FROH)
    
    cli_alert_info(paste0("\nTotal number of candidates: ",nrow(candidates_data)))
    
    # Save metrics to scenario file
    write.table(candidates_data, paste0(scenario_folder,"candidates_metrics.txt"),
                append = append, col.names = col.names,
                quote = F, row.names = F)
    
    # Merge dataframes (progeny)
    df = merge(fakePed[,c(1:3,8)], BLUPF90_EBVs[,c("original_id","solution")],
               by.x = "ID", by.y = "original_id",
               all.x = TRUE)
    df = merge(df, Fped, by = "ID", all.x = TRUE)
    df = merge(df, Fg, by = "ID", all.x = TRUE)
    
    rm(fakePed, BLUPF90_EBVs, Fped, Fg)
    
    write.table(df,paste0(scenario_folder,"fakeProgeny.txt"),
                append = append, col.names = col.names,
                quote = F, row.names = F)
    
    # Pre-selection of male candidates
    best_EBVs = candidates_data %>% 
        filter(sex == "M") %>% 
        slice_max(solution, prop = pre_selection_males_porc) %>%
        select(ID)
    
    df = df %>% 
        filter(sire %in% best_EBVs$ID)
    
    cli_alert_info(paste0("\nTotal number of male candidates after pre-selection: ",
                          length(unique(df$sire)),".\n"))
    
    rm(best_EBVs)
    
    # Assign groups
    candidates_year = data.frame(id=pop@id, 
                                 year=unname(unlist(pop@misc))) %>% 
        dplyr::mutate(group = dense_rank(desc(year))) %>% 
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
    # 2nd rank matings by GEBV and Fg
    
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
    
    # Select matings
    matings = df[0,]
    
    group_counts = data.frame(group=c(1:length(male_groups),
                                      1:length(female_groups)), 
                              sex=c(rep("M",length(male_groups)),
                                    rep("F",length(female_groups))), 
                              count=0,
                              max=c(male_groups,female_groups))
    
    # ---- Define priorities ----
    df$priority = 0
    best_dams = df[!duplicated(df$dam),] %>% 
        slice_head(n=sum(female_groups))
    df$priority[!df$dam %in% best_dams$dam] = df$priority[!df$dam %in% best_dams$dam] + 1
    
    n_matings_each_dam = df %>% 
        group_by(dam) %>% 
        summarise(n = n())
    
    more_than_2 = df$dam %in% n_matings_each_dam$dam[n_matings_each_dam$n > 2]
    only_2 = df$dam %in% n_matings_each_dam$dam[n_matings_each_dam$n == 2]
    df[more_than_2,'priority'] = df[more_than_2,'priority']+2
    df[only_2,'priority'] = df[only_2,'priority']+1
    
    while (nrow(matings)<nMatings) {
        # Order by priority
        df = df[order(df$priority),]
        
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
                        #df = df[!df$sire_group == group,]
                        df[df$sire_group == group, 'priority'] = df[df$sire_group == group, 'priority'] + 1
                    } else {
                        #df = df[!df$dam_group == group,]
                        df[df$dam_group == group,'priority'] = df[df$dam_group == group,'priority'] + 1
                        }
                }
            }
        }
    }
    
    # Check for repeated females
    if (length(unique(matings$dam)) < sum(female_groups)) {
        cli_alert_danger("\nRepeated females!\n")
    }
    
    # Report selected groups
    final_dam_groups = as.data.frame(table(matings$dam_group))
    invisible(apply(final_dam_groups, 1, 
          function(row) cli_alert_info(paste0("\n",row["Freq"], 
                                              " Females selected from group ", 
                                              row["Var1"]), "\n")))
    
    final_sire_groups = as.data.frame(matings %>% 
        group_by(sire) %>% 
        summarise(n = n(), sire_group=mean(sire_group)))
    
    invisible(apply(as.data.frame(table(final_sire_groups$sire_group)), 1, 
                    function(row) cli_alert_info(paste0("\n",row["Freq"], 
                                                        " Males selected from group ", 
                                                        row["Var1"]), "\n")))
    cli_alert_info(paste0("\nMean of ", round(mean(final_sire_groups$n)),
                          " (+/- ", round(sd(final_sire_groups$n)),
                          ") progeny per sire.\n"))
    
    # Breed new generation
    cli_alert_info("\nBreeding new generation.\n")
    new_gen = makeCross(pop, as.matrix(matings[,c(2,1)]))
    
    new_gen = setMisc(x = new_gen,
                      node = "yearOfBirth",
                      value = year)
    
    # Merge new generation and older candidates
    older_candidates = c(unique(matings$sire[!matings$sire_group == length(male_groups)]),
                         unique(matings$dam[!matings$dam_group == length(female_groups)]))
    
    candidates_data = candidates_data %>% 
        filter(ID %in% older_candidates)
    
    candidates = mergePops(list(new_gen,pop[pop@id %in% older_candidates]))
    
    return(candidates)
}
