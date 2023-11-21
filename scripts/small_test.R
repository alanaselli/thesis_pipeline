library(AlphaSimR)
library(dplyr)
source("AlphaSimR_addOn.R")

setwd("../simulations/2023_11_21_small_test/")

geno_path = "01_genotypes/"

# ---- Generate founder genomes ----

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

# ---- Expand population ----

nCrosses = 100
expandedPop = randCross(
    pop = founderPop,
    nCrosses = nCrosses,
    nProgeny = 1
)
rm(founderPop, founderGenomes)
expandedPop@id = add_prefix(expandedPop, "B")
year = year + 1

for (gen in 1:10) {
    nCrosses = round(nCrosses*2)
    expandedPop = randCross(
        pop = expandedPop,
        nCrosses = nCrosses,
        nProgeny = 1
    )
    expandedPop@id = add_prefix(expandedPop, "B")
    year = year + 1
}

# ---- Recent Population ----

year = year + 1
recentPop = makeRecentPop(previous_pop = expandedPop, 
                          males_each_year = c(400,100), 
                          females_each_year = c(250,150,75,25), 
                          nCrosses = 1000, 
                          year = year,
                          years_of_breeding = 10,
                          addPrefix = "C",
                          rec_data_param = list(paste0(geno_path,"pedigree.txt"),
                                                "RecentA",
                                                TRUE))
recentPop = recentPop[[1]]

rm(expandedPop)

# Save genotypes
writePlink(recentPop, paste0(geno_path,"recent"))

# Change positions in map file
source("../../scripts/convert_SNP_pos_small_test.R")

# ---- Start new selection process ----

# Run BLUPF90
runBLUPF90() # This function is very specific to my data

# Select candidates
parents = selectCandidates(BLUP, recentPop)

# Perform matings
# 1000 crosses to generate 500 females and select 250 of them
candidates = randCross2(female_parents, male_parents, nCrosses = 1000)


