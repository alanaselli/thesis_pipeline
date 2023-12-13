library(sys)
library(AlphaSimR)
library(dplyr)
library(cli)

AlphaSimR_addOn = list.files("scripts/AlphaSimR_addOn/", full.names = T)
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
#founderPop@id = add_prefix(founderPop, "A")
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
#expandedPop@id = add_prefix(expandedPop, "B")
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
    #expandedPop@id = add_prefix(expandedPop, "B")
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
# Include sires and dams from older generations
writePlink(recentPop, paste0(geno_path,"recent"))

# Change positions in map file
adjust_pos("01_genotypes/recent.map",
           c(950), c(1e8))

# Prepare SNP file of recent
system(command="./scripts/prepare_snp_file.sh 01_genotypes/recent.ped 01_genotypes/new_map.map 05_BLUPF90/snp_file.txt")

# ---- Start new selection process ----
year = year + 1

# Traditional EBV and F
fakePed = makeFakePed(recentPop)

fakePed[,c(4:7)] = 0
fakePed$gen = year

system(command="cp 01_genotypes/pedigree.txt scenario_01/")

write.table(fakePed,"scenario_01/pedigree.txt",
            col.names = F, append = T,
            row.names = F, quote = F)

system(command="./scripts/extract_ped.sh scenario_01/pedigree.txt 05_BLUPF90/dat1.txt 05_BLUPF90/ped1.txt 10")

# Run BLUPF90 without genomic data
runBLUPF90(param_card="renum.txt")

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

# Run BLUPF90 with genomic data
runBLUPF90(param_card="renum_genomic.txt")

# Fit RR-BLUP model for genomic predictions (candidates)
BLUP = GBLUP_AlphaSimR(recentPop)

# Read Fg
Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                col.names = c("ID","Fg"),
                colClasses = c("character","numeric"))

# Merge dataframes (candidates)
candidates_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                        by.x = "ID", by.y = "original_id",
                        all.x = TRUE)
candidates_data = merge(candidates_data, Fped,
                        by = "ID", all.x = TRUE)

# Save metrics to scenario file
write.table(candidates_data, "scenario_01/candidates_metrics.txt",
            append = F, col.names = T,
            quote = F, row.names = F)

# Merge dataframes (progeny)
progeny_data = merge(fakePed[,c(1:3,8)], BLUPF90_EBVs[,c("original_id","solution")],
                     by.x = "ID", by.y = "original_id",
                     all.x = TRUE)
progeny_data = merge(progeny_data, Fped,
                     by = "ID", all.x = TRUE)

write.table(progeny_data,"scenario_01/fakeProgeny.txt",
            append = F, col.names = T,
            quote = F, row.names = F)

# Pre-selection of male candidates
worse_EBVs = candidates_data %>% 
    filter(sex == "M") %>% 
    slice_min(solution, prop = 0.2) %>% # worst 20% by BLUPF90 EBV
    select(ID)
    
progeny_data = progeny_data %>% 
    filter(!sire %in% worse_EBVs$ID)

# Assign groups
candidates_year = data.frame(id=recentPop@id, 
                             year=unname(unlist(recentPop@misc))) %>% 
    mutate(group = dense_rank(desc(year))) %>% 
    select(!year)

progeny_data = merge(progeny_data, candidates_year, 
                     by.x = "dam", by.y = "id") %>% 
    rename(dam_group = group)
progeny_data = merge(progeny_data, candidates_year,
                     by.x = "sire", by.y = "id") %>% 
    rename(sire_group = group)


# 1st remove matings with high Fped
# 2nd rank matings by EBV and Fped
Fped_percentage = 0.2
progeny_data = progeny_data %>% 
    slice_min(Fped, prop = Fped_percentage) %>%  # take only 20% of matings with lowest Fped
    arrange(desc(round(solution,2)),Fped)

length(unique(progeny_data$dam)) # there are 1400 females left (good)

# Split each female group in a new df
#list_of_dfs = split(progeny_data, progeny_data$dam_group)
#matings_per_group = c(300,250,200,150,100)

# Select matings
number_of_matings = 0
matings = progeny_data[0,]

#sire_counts = data.frame(sire=character(0), count=integer(0))

group_counts = data.frame(group=c(1,2,1:5), 
                          sex=c("M","M",rep("F",5)), 
                          count=0,
                          max=c(40,10,c(300,250,200,150,100)))
#df = list_of_dfs[[1]]
df = progeny_data
while (nrow(matings)<1000) {
    # Select the first sire
    sel_sire = df[1,'sire']
    
    # Select the first 20 matings of that sire
    sire_matings = df %>% 
        filter(sire == sel_sire) %>% 
        slice_head(n=20)
    
    # Remove dams from selection pool
    df = df[!df$dam %in% sire_matings$dam,]
    
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
    for (s in c("M","F")) {
        if (s == "M") {groups = 1:2} else {groups = 1:5}
        for (g in groups) {
            g_count = group_counts$count[group_counts$sex == s & group_counts$group == g]
            g_max = group_counts$max[group_counts$sex == s & group_counts$group == g]
            if (g_count > g_max) {
                if (s == "M") {
                    df = df[!df$sire_group == g,]
                } else {df = df[!df$dam_group == g,]}
            }
        }
    }
}

# Check for repeated females
which(duplicated(matings$dam))

cli_alert_success("All processes of the simulation are completed.")