library(sys)
library(AlphaSimR)
library(dplyr)
library(cli)
source("scripts/AlphaSimR_addOn.R")

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
SP$addSnpChip(nSnpPerChr = 950)

# Add simple additive trait
SP$addTraitA(nQtlPerChr = c(50), mean = 0, var = 1)
SP$setVarE(h2 = 0.3)

# ---- Generate initial populations ----
founderPop = newPop(founderGenomes)
founderPop@id = add_prefix(founderPop, "A")
year = 0

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
expandedPop@id = add_prefix(expandedPop, "B")
year = year + 1

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
    expandedPop@id = add_prefix(expandedPop, "B")
    year = year + 1
    rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
             "Expanded", year, append = TRUE)
    cli_progress_update()
}
cli_alert_info(paste0("\n",expandedPop@nInd, 
                      " individuals in the last expansion generation.\n"))

# ---- Recent Population ----

cli_h2("\nGenerating recent population.\n")
year = year + 1
recentPop = makeRecentPop(previous_pop = expandedPop, 
                          males_each_year = c(400,100), 
                          females_each_year = c(250,150,75,25), 
                          nCrosses = 1000, 
                          year = year,
                          years_of_breeding = 10,
                          addPrefix = "C",
                          rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                                "Recent",
                                                TRUE))
recentPop = recentPop[[1]]

rm(expandedPop)

# Save genotypes
writePlink(recentPop, paste0(geno_path,"recent"))

# Change positions in map file
source("scripts/convert_SNP_pos.R")

# ---- Start new selection process ----
cli_alert_info("\nRunning evaluation for the last recent generation.\n")
# Run BLUPF90
runBLUPF90("recent.ped",
           min_gen=year-10) # This function is very specific to my data

year_1 = year+1
year_2 = year+1

# ---- Scenario 1 ----
cli_h3("\nInitiating scenario 1.\n")
# Select candidates
parents_1 = selectCandidates(recentPop, 
                           "scenario_1/only_EBV.csv", 
                           FALSE, method=1,
                           top_ebv=c(50,250))

# Perform matings
# 1000 crosses to generate 500 females and select 250 of them
candidates_1 = randCross(parents_1, nCrosses = 1000)

# Save pedigree data
cli_alert_info("\nRecording data for gen 1 of scenario 1.\n")
rec_data(paste0(geno_path,"pedigree.txt"), 
         candidates_1, "scenario_1", year_1, append = TRUE)

# Save genotypes
writePlink(candidates_1, paste0(geno_path,"sc_1_gen_1"))

# Continue this process for 10 generations
for (i in 1:10) {
    cli_alert_info(paste0("\nInitiating gen ",i," of scenario 1.\n"))
    year_1=year_1+1
    ped_name = paste0("sc_1_gen_",i)
    
    # Run BLUPF90
    runBLUPF90(paste0(ped_name,".ped"),
               min_gen=year_1-10)
    
    # Select candidates
    parents_1 = selectCandidates(candidates_1, 
                                 paste0("scenario_1/",
                                        ped_name,
                                        ".csv"), 
                                 FALSE, method=1,
                                 top_ebv=c(50,250))
    
    # Perform matings
    candidates_1 = randCross(parents_1, nCrosses = 1000)
    
    # Save pedigree data
    rec_data(paste0(geno_path,"pedigree.txt"), 
             candidates_1, "scenario_1", year_1, append = TRUE)
    
    # Save genotypes
    ped_name = paste0("sc_1_gen_",i+1)
    writePlink(candidates_1, paste0(geno_path,ped_name))
}

rm(parents_1, candidates_1)
cli_alert_success("\nScenario 1 completed.\n")

# ---- Scenario 2 ----
cli_h3("\nInitiating scenario 2.\n")
# Select candidates
parents_2 = selectCandidates(recentPop, 
                             "scenario_2/EBV_Fg.csv", 
                             FALSE, method=2,
                             top_ebv=c(50,250),
                             Fg_threshold = 1.05)

# Perform matings
candidates_2 = randCross(parents_2, nCrosses = 1000)

# Save pedigree data
cli_alert_info("\nRecording data for gen 1 of scenario 2.\n")
rec_data(paste0(geno_path,"pedigree.txt"), 
         candidates_2, "scenario_2", year_2, append = TRUE)

# Save genotypes
writePlink(candidates_2, paste0(geno_path,"sc_2_gen_2"))

# Continue this process for 10 generations
for (i in 1:10) {
    cli_alert_info(paste0("\nInitiating gen ",i," of scenario 2.\n"))
    year_2=year_2+1
    ped_name = paste0("sc_2_gen_",i)
    
    # Run BLUPF90
    runBLUPF90(paste0(ped_name,".ped"),
               min_gen=year_2-10)
    
    # Select candidates
    parents_2 = selectCandidates(candidates_2, 
                                 "scenario_2/only_EBV.csv", 
                                 TRUE, method=2,
                                 top_ebv=c(50,250),
                                 Fg_threshold = 1.05)
    
    # Perform matings
    candidates_2 = randCross(parents_2, nCrosses = 1000)
    
    # Save pedigree data
    rec_data(paste0(geno_path,"pedigree.txt"), 
             candidates_2, "scenario_2", year_2, append = TRUE)
    
    # Save genotypes
    ped_name = paste0("sc_2_gen_",i+1)
    writePlink(candidates_2, paste0(geno_path,ped_name))
}
cli_alert_success("\nScenario 2 completed.\n")
cli_alert_success("All processes of the simulation are completed.")

