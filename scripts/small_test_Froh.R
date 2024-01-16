library(sys)
library(AlphaSimR)
library(dplyr)
library(cli)

AlphaSimR_addOn = list.files("scripts/AlphaSimR_addOn", full.names = T)
for (script in AlphaSimR_addOn) {
    source(script)
}
rm(AlphaSimR_addOn, script)

cli_h1("\nInitiating data simulation\n")

geno_path = "01_genotypes/"

# ---- Generate founder genomes ----
start_time = Sys.time()

times = data.frame(activity = c("start_simulation"),
                   timestamp=c(format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cli_h2("\nGenerating founder genomes.\n")

founderGenomes = runMacs2(
    nInd = 100,
    nChr = 5,
    Ne = 200,
    segSites = c(1000,800,700,600,500)
)

SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
#SP$setTrackPed(TRUE)
SP$addSnpChip(nSnpPerChr = c(950,780,700,600,500))

# Add simple additive trait
SP$addTraitA(nQtlPerChr = c(50,20,0,0,0), mean = 0, var = 1)
SP$setVarE(h2 = 0.3)

# ---- Generate initial populations ----
founderPop = newPop(founderGenomes)
year = 0
founderPop = setMisc(x = founderPop,
                     node = "yearOfBirth",
                     value = year)

rec_data(paste0(geno_path,"pedigree.txt"), founderPop, 
         "Founder", year, append = FALSE)

# ---- Expand population ----

cli_h2("\nExpanding founder population.\n")

nCrosses = 100
expandedPop = randCross(
    pop = founderPop,
    nCrosses = nCrosses,
    nProgeny = 1
)
rm(founderPop, founderGenomes)
year = year + 1
expandedPop = setMisc(x = expandedPop,
                     node = "yearOfBirth",
                     value = year)

rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
         "Expanded", year, append = TRUE)

cli_progress_bar("Expanding population", total = 5)
for (gen in 1:5) {
    nCrosses = round(nCrosses*2)
    expandedPop = randCross(
        pop = expandedPop,
        nCrosses = nCrosses,
        nProgeny = 1
    )
    year = year + 1
    # expandedPop = setMisc(x = expandedPop,
    #                       node = "yearOfBirth",
    #                       value = year)
    
    rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
             "Expanded", year, append = TRUE)
    cli_progress_update()
}
cli_alert_info(paste0("\n",expandedPop@nInd, 
                      " individuals in the last expansion generation.\n"))

# ---- Recent Population ----

cli_h2("\nGenerating recent population.\n")

males_each_year = c(3,2)
females_each_year = c(20,15,10,5)
nCrosses = 50

recent_pop_list = makeRecentPop(previous_pop = expandedPop, 
                          males_each_year = males_each_year, 
                          females_each_year = females_each_year, 
                          nCrosses = nCrosses, 
                          year = year,
                          years_of_breeding = 10,
                          return_breeding_groups = TRUE,
                          rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                                "Recent",
                                                TRUE))
year = year + 10

# Progeny and parents except for oldest groups
recentPop = mergePops(recent_pop_list[c(1,2,4:6)])

rm(expandedPop, recent_pop_list)

# Save genotypes
writePlink(recentPop, paste0(geno_path,"recent"))

times = rbind(times, c("end_simulation", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Change positions in map file
adjust_pos("01_genotypes/recent.map",
           c(950,780,700,600,500), c(1.2e8,1e8,9e7,7e7,5e7))

# ---- Start new selection process ----

# GEBV and Froh

scenario_folder = "scenario_03/"

# Copy pedigree file to scenario folder
system(command=paste0("cp 01_genotypes/pedigree.txt ",scenario_folder))

# Copy genotypes to scenario folder
system(command=paste0("cp 01_genotypes/recent.ped ",scenario_folder))
system(command=paste0("cp 01_genotypes/new_map.map ",scenario_folder))

scenario_03 = recentPop

for (gen in 1:10) {
    
    times = rbind(times, c(paste0("scenario3_","gen",gen), 
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    cli_alert_info(paste0("\nInitiating generation ",gen,
                          " of scenario 3 (EBV + F).\n"))
    year = year + 1
    
    if (gen == 1) {append = FALSE} else {append = TRUE}
    
    scenario_03 = mate_selection_GEBV_Froh(pop = scenario_03,
                                          year = year,
                                          scenario_folder = scenario_folder,
                                          pre_selection_males_porc = 0.8, # percentage to remove
                                          Froh_percentage = 0.2, # percentage to keep
                                          male_groups = males_each_year,
                                          female_groups = females_each_year,
                                          nMatings = nCrosses,
                                          nMatings_per_sire=10,
                                          append = append)
    
    # Take only the last generation to add records
    last_gen = data.frame(id=scenario_03@id, 
                            year=unname(unlist(scenario_03@misc)))
    last_gen = last_gen[last_gen$year == year,]
    record_data = scenario_03[scenario_03@id %in% last_gen$id]
    
    # Record pedigree
    rec_data(paste0(scenario_folder,"pedigree.txt"), 
             record_data, scenario_folder, 
             year, append = TRUE)
    
    # Save genotypes
    writePlink(record_data, paste0(scenario_folder,"scenario_03"))
    
    # Merge genotypes
    system(command=paste0('plink --cow --ped ',
                          scenario_folder,'recent.ped --map ',
                          scenario_folder,'new_map.map --merge ',
                          scenario_folder,'scenario_03.ped ',
                          scenario_folder,'new_map.map --make-bed --recode --out ',
                          scenario_folder,'recent'))
    
    # Prepare SNP file
    system(command=paste0("scripts/prepare_snp_file.sh ",
                          scenario_folder,"recent.ped ",
                          "01_genotypes/new_map.map 05_BLUPF90/snp_file.txt"))
    
    cli_alert_success(paste0("\nGeneration ",gen,
                            " of scenario 1 (EBV + F) is complete.\n"))
}

times = rbind(times, c("end_script", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
write.csv(times,"simulation_timestamps.csv", quote = F, row.names = F)

end_time = Sys.time()
total_time = end_time-start_time

cli_alert_info(paste0("\nRunning time of the simulation: ", total_time,"\n"))
cli_alert_success("All processes of the simulation are completed.")
