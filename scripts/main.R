#!/usr/bin/env Rscript
library("optparse")
library(sys)
library(AlphaSimR)
library(dplyr)
library(cli)

option_list = list(
    make_option(c("-g", "--nGenerations"), type="integer", default=NULL, 
                help="number of generations in scenarios", metavar="character"),
    make_option(c("-m", "--males_each_year"), type="character", default=NULL, 
                help="number of males each year separated by spaces"),
    make_option(c("-f", "--females_each_year"), type="character", default=NULL, 
                help="number of females each year separated by spaces"),
    make_option(c("-c", "--nCrosses"), type="integer", default=NULL, 
                help="number of crosses", metavar="integer"),
    make_option(c("-n", "--folder_name"), type="character", default=NULL, 
                help="folder name to store scenarios", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

number_of_generations = opt$nGenerations
males_each_year = as.double(unlist(strsplit(opt$males_each_year, " ")))
females_each_year = as.double(unlist(strsplit(opt$females_each_year, " ")))
nMatings = opt$nCrosses
folder_name = opt$folder_name

cli_alert(paste0("\nArguments passed:\n"))
cli_alert(paste0("\nnumber_of_generations = ", number_of_generations))
cli_alert(paste0("\nmales_each_year = "))
cli_alert(paste0("\nfemales_each_year = "))
cli_alert(paste0("\nnCrosses = ", nMatings))
cli_alert(paste0("\nfolder_name = ", folder_name))

dir.create(file.path(folder_name), showWarnings = FALSE)

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

#segSites = c(4000,3500,3200,2500,2500)
bp = c(1.58e8,1.37e8,1.21e8)
histNe = c(1000, 2000, 4000, 7500, 10000)
histGen = c(50, 100, 500, 1000, 1500)
nSnpPerChr = c(3000,2500,2200,2000,1900)
nQtlPerChr = c(100,75,50,10,0)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2374996/table/T4/?report=objectonly
chr_size = c(1.46e8, 1.25e8, 1.16e8, 1.1e8, 1.18e8)

founderGenomes = runMacs2(
    nInd = 200,
    nChr = 5,
    Ne = 200,
    bp = bp,
    histNe = histNe,
    histGen = histGen
    #segSites = segSites
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

nCrosses_e = 100
expandedPop = randCross(
    pop = founderPop,
    nCrosses = nCrosses_e,
    nProgeny = 1
)
rm(founderPop, founderGenomes)
year = year + 1
expandedPop = setMisc(x = expandedPop,
                     node = "yearOfBirth",
                     value = year)

rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
         "Expanded", year, append = TRUE)

cli_progress_bar("Expanding population", total = 30)
for (gen in 1:30) {
    nCrosses_e = round(nCrosses_e*1.2)
    expandedPop = randCross(
        pop = expandedPop,
        nCrosses = nCrosses_e,
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
                          nCrosses = nMatings, 
                          year = year,
                          years_of_breeding = 10,
                          return_breeding_groups = TRUE,
                          rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                                "Recent",
                                                TRUE))
year = year + 10

# Progeny and parents except for oldest groups
recentPop = mergePops(recent_pop_list[c(1,2,4:6)])
cli_alert(paste0("\nN recent: ", recentPop@nInd))

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

scenario_folder = file.path(folder_name,"scenario_01/")
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
    
    cli_h2(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1 
    
    if (gen == 1) {
        append = FALSE
    } else {
            append = TRUE
            }
    
    if (gen == number_of_generations) {
        last_gen = TRUE
    } else {
        last_gen = FALSE
        }
    
    scenario_pop = mate_selection_EBV_Fped(pop = scenario_pop,
                                          year = year,
                                          scenario_folder = scenario_folder,
                                          pre_selection_males_porc = 0.8, # percentage to remove
                                          Fped_percentage = 0.5, # percentage to keep
                                          male_groups = males_each_year,
                                          female_groups = females_each_year,
                                          nMatings = nMatings,
                                          nMatings_per_sire=10,
                                          append = append,
                                          last_gen = last_gen)
    
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

scenario_folder = file.path(folder_name,"scenario_02/")
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
    
    cli_h2(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1
    
    if (gen == 1) {
        append = FALSE
    } else {
        append = TRUE
    }
    
    if (gen == number_of_generations) {
        last_gen = TRUE
    } else {
        last_gen = FALSE
    }
    
    scenario_pop = mate_selection_GEBV_Fg(pop = scenario_pop,
                                            year = year,
                                            scenario_folder = scenario_folder,
                                            pre_selection_males_porc = 0.8, # percentage to remove
                                            Fg_percentage = 0.5, # percentage to keep
                                            male_groups = males_each_year,
                                            female_groups = females_each_year,
                                            nMatings = nMatings,
                                            nMatings_per_sire=10,
                                            append = append,
                                          last_gen = last_gen)
    
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

scenario_folder = file.path(folder_name,"scenario_03/")
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
    
    cli_h2(paste0("\nInitiating generation ",gen,
                          " of scenario ", scenario_name, ".\n"))
    year = year + 1
    
    if (gen == 1) {
        append = FALSE
    } else {
        append = TRUE
    }
    
    if (gen == number_of_generations) {
        last_gen = TRUE
    } else {
        last_gen = FALSE
    }
    
    scenario_pop = mate_selection_GEBV_Froh(pop = scenario_pop,
                                          year = year,
                                          scenario_folder = scenario_folder,
                                          pre_selection_males_porc = 0.8, # percentage to remove
                                          Froh_percentage = 0.5, # percentage to keep
                                          male_groups = males_each_year,
                                          female_groups = females_each_year,
                                          nMatings = nMatings,
                                          nMatings_per_sire=10,
                                          append = append,
                                          last_gen = last_gen)
    
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
write.csv(times,file.path(folder_name,"simulation_timestamps.csv"), quote = F, row.names = F)

end_time = proc.time()
total_time = end_time-start_time

cli_alert_info(paste0("\nRunning time of the simulation: ", total_time," seconds.\n"))
cli_alert_success("All processes of the simulation have been completed.")
