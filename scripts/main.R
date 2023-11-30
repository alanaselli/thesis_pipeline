library(sys)
library(AlphaSimR)
library(cli)
source("scripts/AlphaSimR_addOn.R")

geno_path = "01_genotypes/"

cli_h1("\nInitiating data simulation\n")

# ---- Generate founder genomes ----
start_time = Sys.time()
cli_h2("\nGenerating founder genomes.\n")

founderGenomes = runMacs2(
    nInd = 2000,
    nChr = 5,
    Ne = 200,
    histNe = c(1000, 2000, 4000, 7500, 10000),
    histGen = c(50, 100, 500, 1000, 1500)
)

SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$addSnpChip(nSnpPerChr = 3000)

# Add simple additive trait
SP$addTraitA(nQtlPerChr = c(50, 30, 15, 5, 0), mean = 0, var = 1)
SP$setVarE(h2 = 0.3)

# Save QTL effects
write.table(paste(SP$traits), paste0(geno_path,"QTL_effects.txt"))

# Save QTL map
QTL_map = getQtlMap(trait = 1, sex = "A")
write.csv(QTL_map, paste0(geno_path,"QTL_map.csv"), row.names = F, quote = F)
rm(QTL_map)

# ---- Generate initial populations ----
founderPop = newPop(founderGenomes)
year = 0
rec_data(paste0(geno_path,"pedigree.txt"), founderPop, 
         "Founder", year, append = FALSE)

cli_alert_info("\nWriting PLINK file for founder_sample_100\n")
founder_sample_100 = selectInd(founderPop, nInd = 100, use = "rand")
writePlink(founder_sample_100, paste0(geno_path,"founder_sample_100"))

cli_alert_info("\nWriting founder QTLs\n")
writePlink(founder_sample_100,paste0(geno_path,"founder_QTL"), useQtl = TRUE)

rm(founderGenomes,founder_sample_100)

# Change positions in map file
source("scripts/convert_SNP_pos.R")

# ---- Expand population ----

cli_h2("\nExpanding founder population.\n")

nCrosses = 2000
expandedPop = randCross(
    pop = founderPop,
    nCrosses = nCrosses,
    nProgeny = 1
)
rm(founderPop)

year = year + 1
rec_data(paste0(geno_path,"pedigree.txt"), expandedPop, 
         "Expanded", year, append = TRUE)

cli_progress_bar("Expanding population", total = 70)
for (gen in 1:70) {
    nCrosses = round(nCrosses*1.05)
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

cli_alert_info("\nWriting PLINK file for expanded_sample_100\n")
expanded_sample_100 = selectInd(expandedPop, nInd = 100, use = "rand")
writePlink(expanded_sample_100, paste0(geno_path,"expanded_sample_100"))

cli_alert_info("\nWriting expanded QTLs\n")
writePlink(expanded_sample_100, paste0(geno_path,"expanded_QTL"), useQtl = TRUE)

rm(expanded_sample_100)

# ---- Recent Population ----

cli_h2("\nGenerating recent population.\n")
year = year + 1
recentPop = makeRecentPop(previous_pop = expandedPop, 
                     males_each_year = c(200,50), 
                     females_each_year = c(1500,1250,1000,750,500), 
                     nCrosses = 5000, 
                     year = year,
                     years_of_breeding = 30,
                     rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                           "Recent",
                                           TRUE))
recentPop = recentPop[[1]]

cli_alert_info("\nWriting PLINK file for recent\n")
writePlink(recentPop, paste0(geno_path,"recent"))

cli_alert_info("\nWriting PLINK file for recent_sample_100\n")
recentPop_sample_100 = selectInd(recentPop, nInd = 100, use = "rand")
writePlink(recentPop_sample_100, paste0(geno_path,"recent_sample_100"))

cli_alert_info("\nWriting recent Pop QTLs\n")
writePlink(recentPop_sample_100, paste0(geno_path,"recent_QTL"), useQtl = TRUE)

rm(recentPop_sample_100)

cli_alert_success("\nMain simulation completed.\n")

# ---- Start new selection process ----
cli_alert_info("\nRunning evaluation for the last recent generation.\n")

# Run BLUPF90
runBLUPF90("recent.ped",
           min_gen=year-10) # This function is very specific to my data

year_1 = year+1

# Select candidates
parents_1 = selectCandidates(pop=recentPop, 
                             file_name=paste0("scenario_1/select_by_EBV.csv"), 
                             append=FALSE, 
                             method=1,
                             top_ebv=c(250,1250))

# Perform matings
# 10000 crosses to generate 5000 females and select 2500 of them
candidates_1 = randCross(parents_1, nCrosses = 5000)

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
    parents_1 = selectCandidates(pop=candidates_1, 
                                 file_name=paste0("scenario_1/select_by_EBV.csv"), 
                                 append=TRUE, 
                                 method=1,             # 1=EBV; 2=EBV+Fg
                                 top_ebv=c(250,1250))  # c(nMales, nFemales)
    
    # Perform matings
    candidates_1 = randCross(parents_1, nCrosses = 5000)
    
    # Save pedigree data
    rec_data(paste0(geno_path,"pedigree.txt"), 
             candidates_1, "scenario_1", year_1, append = TRUE)
    
    # Save genotypes
    ped_name = paste0("sc_1_gen_",i+1)
    writePlink(candidates_1, paste0(geno_path,ped_name))
}

end_time = Sys.time()
total_time = end_time-start_time

cli_alert_info(paste0("\nRunning time of the simulation: ", total_time,"\n"))
cli_alert_success("\nAll processes of the simulation are completed.\n")