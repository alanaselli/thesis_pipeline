library(sys)
library(AlphaSimR)
library(cli)
source("AlphaSimR_addOn.R")

cli_h1("\nInitiating data simulation\n")

# ---- Generate founder genomes ----
start_time = Sys.time()
cli_h2("\nGenerating founder genomes.\n")

founderGenomes = runMacs2(
    nInd = 2000,
    nChr = 4,
    Ne = 200,
    histNe = c(1000, 2000, 4000, 7500, 10000),
    histGen = c(50, 100, 500, 1000, 1500)
)

SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$addSnpChip(nSnpPerChr = 2000)

# Add simple additive trait
SP$addTraitA(nQtlPerChr = c(30, 15, 5, 0), mean = 0, var = 1)
SP$setVarE(h2 = 0.3)

# Save QTL effects
write.table(paste(SP$traits), "01_genotypes/QTL_effects.txt")

# Save QTL map
QTL_map = getQtlMap(trait = 1, sex = "A")
write.csv(QTL_map, "01_genotypes/QTL_map.csv", row.names = F, quote = F)
rm(QTL_map)

# ---- Generate initial populations ----
founderPop = newPop(founderGenomes)
founderPop@id = add_prefix(founderPop, "A")
year = 0
rec_data("01_genotypes/pedigree.txt", founderPop, 
         "Founder", year, append = FALSE)

cli_alert_info("\nWriting PLINK file for founder_sample_100\n")
founder_sample_100 = selectInd(founderPop, nInd = 100, use = "rand")
writePlink(founder_sample_100, "01_genotypes/founder_sample_100")

cli_alert_info("\nWriting founder_phased\n")
founder_phased = getPhasedHaplo(founder_sample_100)
write.table(founder_phased, "01_genotypes/founder_phased.txt", 
            quote = F, col.names = F)

cli_alert_info("\nWriting founder QTLs\n")
writePlink(founder_sample_100, "01_genotypes/founder_QTL", useQtl = TRUE)

rm(founderGenomes,founder_sample_100,founder_phased)

# ---- Expand for 100 generations ----

cli_h2("\nExpanding founder population.\n")

nCrosses = 2000
expandedPop = randCross(
    pop = founderPop,
    nCrosses = nCrosses,
    nProgeny = 1
)
rm(founderPop)
expandedPop@id = add_prefix(expandedPop, "B")
year = year + 1
rec_data("01_genotypes/pedigree.txt", expandedPop, 
         "Expanded", year, append = TRUE)

cli_progress_bar("Expanding population", total = 50)
for (gen in 1:50) {
    nCrosses = round(nCrosses*1.05)
    expandedPop = randCross(
        pop = expandedPop,
        nCrosses = nCrosses,
        nProgeny = 1
    )
    expandedPop@id = add_prefix(expandedPop, "B")
    year = year + 1
    rec_data("01_genotypes/pedigree.txt", expandedPop, 
             "Expanded", year, append = TRUE)
    cli_progress_update()
} # 115335 individuals in the final generation (11)
cli_alert_info(paste0("\n",expandedPop@nInd, 
                      " individuals in the last expansion generation.\n"))

cli_alert_info("\nWriting PLINK file for expanded_sample_100\n")
expanded_sample_100 = selectInd(expandedPop, nInd = 100, use = "rand")
writePlink(expanded_sample_100, "01_genotypes/expanded_sample_100")

cli_alert_info("\nWriting expanded_phased\n")
expanded_phased = getPhasedHaplo(expanded_sample_100)
write.table(expanded_phased, "01_genotypes/expanded_phased.txt", 
            quote = F, col.names = F)

cli_alert_info("\nWriting expanded QTLs\n")
writePlink(expanded_sample_100, "01_genotypes/expanded_QTL", useQtl = TRUE)

rm(expanded_sample_100, expanded_phased)

# ---- Recent Population ----

cli_h2("\nGenerating recent population A.\n")
year = year + 1
recentPop = makeRecentPop(previous_pop = expandedPop, 
                     males_each_year = c(400,100), 
                     females_each_year = c(5000,2500,1500,750,250), 
                     nCrosses = 10000, 
                     year = year,
                     years_of_breeding = 19,
                     addPrefix = "C",
                     rec_data_param = list("01_genotypes/pedigree.txt",
                                           "RecentA",
                                           TRUE))
recentPop = recentPop[[1]]

# Fit RR-BLUP model for genomic predictions
ans = RRBLUP(recentPop, simParam=SP)
recentPop = setEBV(recentPop, ans, simParam=SP)

cli_alert_info("\nCalculating and writing EBVs\n")
BLUP = data.frame(ID = pop@id, 
                  EBV = pop@ebv, 
                  GV = pop@gv)
names(BLUP) = c('ID','EBV','GV')

write.table(BLUP,"BLUP_AlphaSimR.txt", row.names = F, quote = F)

cli_alert_info("\nWriting PLINK file for recent\n")
writePlink(recentPop, "01_genotypes/recent")

cli_alert_info("\nWriting PLINK file for recent_sample_100\n")
recentPop_sample_100 = selectInd(recentPop, nInd = 100, use = "rand")
writePlink(recentPop_sample_100, "01_genotypes/recent_sample_100")

cli_alert_info("\nWriting recent_phased\n")
recentPop_phased = getPhasedHaplo(recentPop_sample_100)
write.table(recentPop_phased, "01_genotypes/recent_phased.txt", 
            quote = F, col.names = F)

cli_alert_info("\nWriting Pop A QTLs\n")
writePlink(recentPop_sample_100, "01_genotypes/recent_QTL", useQtl = TRUE)

rm(recentPop_sample_100, recentPop_phased)

cli_alert_success("\nSimulation completed.\n")