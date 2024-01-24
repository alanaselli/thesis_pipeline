#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
    stop("The following arguments must be provided: number_of_generations, males_each_year, females_each_year, nCrosses, folder_name.", 
         call.=FALSE)
}

number_of_generations = args[1]
males_each_year = args[2]
females_each_year = args[3]
nCrosses = args[4]
folder_name = args[5]

dir.create(file.path(folder_name), showWarnings = TRUE)

library(sys)
library(AlphaSimR)
library(dplyr)
library(cli)

start_time = proc.time()

AlphaSimR_addOn = list.files("scripts/AlphaSimR_addOn", full.names = T)
for (script in AlphaSimR_addOn) {
    source(script)
}
rm(AlphaSimR_addOn, script)

cli_h1("\nInitiating data simulation\n")

geno_path = "01_genotypes/"

# ---- Generate founder genomes ----

times = data.frame(activity = c("start_simulation"),
                   timestamp=c(format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cli_h2("\nGenerating founder genomes.\n")

segSites = c(2100,1975,1850,1710,1700)
nSnpPerChr = c(2000,1900,1800,1700,1700)
nQtlPerChr = c(100,75,50,10,0)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2374996/table/T4/?report=objectonly
chr_size = c(1.46e8, 1.25e8, 1.16e8, 1.1e8, 1.18e8)

founderGenomes = runMacs2(
    nInd = 100,
    nChr = 5,
    Ne = 200,
    segSites = segSites
)

SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$addSnpChip(nSnpPerChr = nSnpPerChr)

# Add simple additive trait
SP$addTraitA(nQtlPerChr = nQtlPerChr, mean = 0, var = 1)
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
    
    rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
             "Expanded", year, append = TRUE)
    cli_progress_update()
}
cli_alert_info(paste0("\n",expandedPop@nInd, 
                      " individuals in the last expansion generation.\n"))

# ---- Recent Population ----

cli_h2("\nGenerating recent population.\n")

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

times = rbind(times, c("end_base_simulation", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Change positions in map file
adjust_pos("01_genotypes/recent.map",
           nSnpPerChr, chr_size)

# ---- Start new selection process ----

# Split time
split_time = year

# ---- Traditional EBV and Fped ----

scenario_folder = file.path(folder_name,"scenario_01")
scenario_name = "EBV_Fped"

dir.create(file.path(scenario_folder), showWarnings = TRUE)

times = rbind(times, c(paste0("start_scenario_",scenario_name), 
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cli_h2(paste0("Starting scenario ", scenario_name,".\n"))

# Copy files
for (FILE in c("01_genotypes/pedigree.txt",
               "01_genotypes/recent.ped",
               "01_genotypes/new_map.map")) {
    system(command=paste0("cp ",FILE," ",scenario_folder))
}
system(command=paste0("cp 01_genotypes/recent.ped ",
                      scenario_folder,"candidates.ped"))

scenario_pop = recentPop

for (gen in 1:number_of_generations) {
    
    times = rbind(times, c(paste0("scenario_",scenario_name,"_gen_",gen), 
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    cli_alert_info(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1 
    
    if (gen == 1) {append = FALSE} else {append = TRUE}
    
    scenario_pop = mate_selection_EBV_Fped(pop = scenario_pop,
                                          year = year,
                                          scenario_folder = scenario_folder,
                                          pre_selection_males_porc = 0.8, # percentage to remove
                                          Fped_percentage = 0.2, # percentage to keep
                                          male_groups = males_each_year,
                                          female_groups = females_each_year,
                                          nMatings = nCrosses,
                                          nMatings_per_sire=10,
                                          append = append)
    
    # Take only the last generation to add records
    last_gen = data.frame(id=scenario_pop@id, 
                          year=unname(unlist(scenario_pop@misc)))
    last_gen = last_gen[last_gen$year == year,]
    record_data = scenario_pop[scenario_pop@id %in% last_gen$id]
    
    # Record pedigree
    rec_data(paste0(scenario_folder,"pedigree.txt"), 
             record_data, scenario_folder, 
             year, append = TRUE)
    
    # Save genotypes
    writePlink(record_data, paste0(scenario_folder,"last_gen"))
    writePlink(scenario_pop, paste0(scenario_folder,"candidates"))
    
    # Merge genotypes
    system(command=paste0('plink --cow --ped ',
                          scenario_folder,'recent.ped --map ',
                          scenario_folder,'new_map.map --merge ',
                          scenario_folder,'last_gen.ped ',
                          scenario_folder,'new_map.map --make-bed --recode --out ',
                          scenario_folder,'recent'))
    
    cli_alert_success(paste0("\nGeneration ",gen,
                             " of scenario ",scenario_name,
                             " is complete.\n"))
}

# ---- GEBV and Fg ----

# Rewind time
year = split_time

scenario_folder = file.path(folder_name,"scenario_02")
scenario_name = "GEBV_Fg"

dir.create(file.path(scenario_folder), showWarnings = TRUE)

times = rbind(times, c(paste0("start_scenario_",scenario_name), 
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cli_h2(paste0("Starting scenario ", scenario_name,".\n"))

# Copy files
for (FILE in c("01_genotypes/pedigree.txt",
               "01_genotypes/recent.ped",
               "01_genotypes/new_map.map")) {
    system(command=paste0("cp ",FILE," ",scenario_folder))
}
system(command=paste0("cp 01_genotypes/recent.ped ",
                      scenario_folder,"candidates.ped"))

scenario_pop = recentPop

for (gen in 1:number_of_generations) {
    
    times = rbind(times, c(paste0("scenario_",scenario_name,"_gen_",gen), 
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    cli_alert_info(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1
    
    if (gen == 1) {append = FALSE} else {append = TRUE}
    
    scenario_pop = mate_selection_GEBV_Fg(pop = scenario_pop,
                                            year = year,
                                            scenario_folder = scenario_folder,
                                            pre_selection_males_porc = 0.8, # percentage to remove
                                            Fg_percentage = 0.2, # percentage to keep
                                            male_groups = males_each_year,
                                            female_groups = females_each_year,
                                            nMatings = nCrosses,
                                            nMatings_per_sire=10,
                                            append = append)
    
    # Take only the last generation to add records
    last_gen = data.frame(id=scenario_pop@id, 
                          year=unname(unlist(scenario_pop@misc)))
    last_gen = last_gen[last_gen$year == year,]
    record_data = scenario_pop[scenario_pop@id %in% last_gen$id]
    
    # Record pedigree
    rec_data(paste0(scenario_folder,"pedigree.txt"), 
             record_data, scenario_folder, 
             year, append = TRUE)
    
    # Save genotypes
    writePlink(record_data, paste0(scenario_folder,"last_gen"))
    writePlink(scenario_pop, paste0(scenario_folder,"candidates"))
    
    # Merge genotypes
    system(command=paste0('plink --cow --ped ',
                          scenario_folder,'recent.ped --map ',
                          scenario_folder,'new_map.map --merge ',
                          scenario_folder,'last_gen.ped ',
                          scenario_folder,'new_map.map --make-bed --recode --out ',
                          scenario_folder,'recent'))
    
    cli_alert_success(paste0("\nGeneration ",gen,
                             " of scenario ",scenario_name,
                             " is complete.\n"))
}

# ---- GEBV and Froh ----

# Rewind time
year = split_time

scenario_folder = file.path(folder_name,"scenario_03")
scenario_name = "GEBV_Froh"

dir.create(file.path(scenario_folder), showWarnings = TRUE)

times = rbind(times, c(paste0("start_scenario_",scenario_name), 
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cli_h2(paste0("Starting scenario ", scenario_name,".\n"))

# Copy files
for (FILE in c("01_genotypes/pedigree.txt",
               "01_genotypes/recent.ped",
               "01_genotypes/new_map.map")) {
    system(command=paste0("cp ",FILE," ",scenario_folder))
}
system(command=paste0("cp 01_genotypes/recent.ped ",
                      scenario_folder,"candidates.ped"))

scenario_pop = recentPop

for (gen in 1:number_of_generations) {
    
    times = rbind(times, c(paste0("scenario_",scenario_name,"_gen_",gen), 
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    cli_alert_info(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1
    
    if (gen == 1) {append = FALSE} else {append = TRUE}
    
    scenario_pop = mate_selection_GEBV_Froh(pop = scenario_pop,
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
    last_gen = data.frame(id=scenario_pop@id, 
                            year=unname(unlist(scenario_pop@misc)))
    last_gen = last_gen[last_gen$year == year,]
    record_data = scenario_pop[scenario_pop@id %in% last_gen$id]
    
    # Record pedigree
    rec_data(paste0(scenario_folder,"pedigree.txt"), 
             record_data, scenario_folder, 
             year, append = TRUE)
    
    # Save genotypes
    writePlink(record_data, paste0(scenario_folder,"last_gen"))
    writePlink(scenario_pop, paste0(scenario_folder,"candidates"))
    
    # Merge genotypes
    system(command=paste0('plink --cow --ped ',
                          scenario_folder,'recent.ped --map ',
                          scenario_folder,'new_map.map --merge ',
                          scenario_folder,'last_gen.ped ',
                          scenario_folder,'new_map.map --make-bed --recode --out ',
                          scenario_folder,'recent'))
    
    cli_alert_success(paste0("\nGeneration ",gen,
                            " of scenario ",scenario_name,
                            " is complete.\n"))
}

times = rbind(times, c("end_script", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
write.csv(times,"simulation_timestamps.csv", quote = F, row.names = F)

end_time = proc.time()
total_time = end_time-start_time

cli_alert_info(paste0("\nRunning time of the simulation: ", total_time," seconds.\n"))
cli_alert_success("All processes of the simulation have been completed.")
