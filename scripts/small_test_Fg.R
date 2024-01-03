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
cli_h2("\nGenerating founder genomes.\n")

founderGenomes = runMacs2(
    nInd = 100,
    nChr = 1,
    Ne = 200,
    segSites = 1000
)

SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
#SP$setTrackPed(TRUE)
SP$addSnpChip(nSnpPerChr = 950)

# Add simple additive trait
SP$addTraitA(nQtlPerChr = c(50), mean = 0, var = 1)
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

cli_progress_bar("Expanding population", total = 10)
for (gen in 1:10) {
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

recent_pop_list = makeRecentPop(previous_pop = expandedPop, 
                          males_each_year = c(40,10), 
                          females_each_year = c(300,250,200,150,100), 
                          nCrosses = 1000, 
                          year = year,
                          years_of_breeding = 10,
                          return_breeding_groups = TRUE,
                          rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                                "Recent",
                                                TRUE))
year = year + 10

# Progeny and parents except for oldest groups
recentPop = mergePops(recent_pop_list[c(1,2,4:7)])

rm(expandedPop, recent_pop_list)

# Save genotypes
writePlink(recentPop, paste0(geno_path,"recent"))

# Change positions in map file
adjust_pos("01_genotypes/recent.map",
           c(950), c(1e8))

# ---- Start new selection process ----

# GEBV and Fg

scenario_folder = "scenario_02/"

# Copy pedigree file to scenario folder
system(command=paste0("cp 01_genotypes/pedigree.txt ",scenario_folder))

# Copy genotypes to scenario folder
system(command=paste0("cp 01_genotypes/recent.ped ",scenario_folder))
system(command=paste0("cp 01_genotypes/recent.ped ",scenario_folder,"last_gen.ped"))
system(command=paste0("cp 01_genotypes/new_map.map ",scenario_folder))

# Prepare SNP file of recent
#system(command="./scripts/prepare_snp_file.sh 01_genotypes/recent.ped 01_genotypes/new_map.map 05_BLUPF90/snp_file.txt")

scenario_02 = recentPop

for (gen in 1:10) {
    cli_alert_info(paste0("\nInitiating generation ",gen,
                          " of scenario 2 (EBV + F).\n"))
    year = year + 1
    
    if (gen == 1) {append = FALSE} else {append = TRUE}
    
    scenario_01 = mate_selection_GEBV_Fg(pop = scenario_02,
                                          year = year,
                                          scenario_folder = scenario_folder,
                                          ped_ROH="recent.ped",
                                          pre_selection_males_porc = 0.8, # percentage to remove
                                          Fped_percentage = 0.2, # percentage to keep
                                          male_groups = c(40,10),
                                          female_groups = c(300,250,200,150,100),
                                          nMatings = 1000,
                                          nMatings_per_sire=20,
                                          append = append)
    
    # Take only the last generation to add records
    last_gen = data.frame(id=scenario_01@id, 
                            year=unname(unlist(scenario_01@misc)))
    last_gen = last_gen[last_gen$year == year,]
    record_data = scenario_01[scenario_01@id %in% last_gen$id]
    
    # Record pedigree
    rec_data("scenario_01/pedigree.txt", record_data, 
             "scenario_01", year, append = TRUE)
    
    # Save genotypes
    writePlink(record_data, "scenario_01/scenario_01")
    
    # Merge genotypes
    system(command='plink --cow --ped scenario_01/recent.ped --map scenario_01/new_map.map --merge scenario_01/scenario_01.ped scenario_01/new_map.map --make-bed --recode --out scenario_01/recent')
    
    # Prepare SNP file
    system(command="scripts/prepare_snp_file.sh scenario_01/recent.ped 01_genotypes/new_map.map 05_BLUPF90/snp_file.txt")
    
    cli_alert_success(paste0("\nGeneration ",gen,
                            " of scenario 1 (EBV + F) is complete.\n"))
}

end_time = Sys.time()
total_time = end_time-start_time

cli_alert_info(paste0("\nRunning time of the simulation: ", total_time,"\n"))
cli_alert_success("All processes of the simulation are completed.")