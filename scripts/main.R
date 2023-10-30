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
rec_data("01_simulations/pedigree.txt", founderPop, 
         "Founder", year, append = FALSE)

cat("\nWriting PLINK file for founder_sample_100\n")
founder_sample_100 = selectInd(founderPop, nInd = 100, use = "rand")
writePlink(founder_sample_100, "01_genotypes/founder_sample_100")

cat("\nWriting founder_phased\n")
founder_phased = getPhasedHaplo(founder_sample_100)
write.table(founder_phased, "01_simulations/founder_phased.txt", 
            quote = F, col.names = F)

cat("\nWriting founder QTLs\n")
writePlink(founder_sample_100, "01_genotypes/founder_QTL", useQtl = TRUE)

rm(founderGenomes,founder_sample_100,founder_sample_10,founder_phased)

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
rec_data("01_simulations/pedigree.txt", expandedPop, 
         "Expanded", year, append = TRUE)

cli_progress_bar("Expanding population", total = 10)
for (gen in 1:100) {
    nCrosses = round(nCrosses*1.5)
    expandedPop = randCross(
        pop = expandedPop,
        nCrosses = nCrosses,
        nProgeny = 1
    )
    expandedPop@id = add_prefix(expandedPop, "B")
    year = year + 1
    rec_data("01_simulations/pedigree.txt", expandedPop, 
             "Expanded", year, append = TRUE)
    cli_progress_update()
} # 115335 individuals in the final generation (11)
cli_alert_info(paste0("\n",expandedPop@nInd, 
                      " individuals in the last expansion generation.\n"))

cat("\nWriting PLINK file for expanded_sample_100\n")
expanded_sample_100 = selectInd(expandedPop, nInd = 100, use = "rand")
writePlink(expanded_sample_100, "01_genotypes/expanded_sample_100")

cat("\nWriting expanded_phased\n")
expanded_phased = getPhasedHaplo(expanded_sample_100)
write.table(expanded_phased, "01_simulations/expanded_phased.txt", 
            quote = F, col.names = F)

cat("\nWriting expanded QTLs\n")
writePlink(expanded_sample_100, "01_genotypes/expanded_QTL", useQtl = TRUE)

rm(expanded_sample_100, expanded_sample_10, expanded_phased)

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
                     rec_data_param = list("01_simulations/pedigree.txt",
                                           "RecentA",
                                           TRUE))
recentPop = recentPop[[1]]

cat("\nWriting PLINK file for recent\n")
writePlink(recentPop, "01_genotypes/recent")

cat("\nWriting PLINK file for recent_sample_100\n")
recentPop_sample_100 = selectInd(recentPop, nInd = 100, use = "rand")
writePlink(recentPop_sample_100, "01_genotypes/recent_sample_100")

cat("\nWriting recent_phased\n")
recentPop_sample_10 = selectInd(recentPop, nInd = 10, use = "rand")
recentPop_phased = getPhasedHaplo(recentPop_sample_100)
write.table(recentPop_phased, "01_simulations/recent_phased.txt", 
            quote = F, col.names = F)

cat("\nWriting Pop A QTLs\n")
writePlink(recentPop_sample_100, "01_genotypes/recent_QTL", useQtl = TRUE)

rm(recentPop_sample_100, recentPop_sample_10, recentPop_phased)