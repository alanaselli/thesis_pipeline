library(detectRUNS)

folder_name = "50_50/"
scenario = "03"

ped = paste0("simulations/",folder_name,'scenario_',scenario,'/last_gen.ped')
map = paste0("simulations/",folder_name,'scenario_',scenario,'/new_map.map')
out = paste0("simulations/",folder_name,'scenario_',scenario,"/last_gen")

# ---- Consecutive runs approach ----

consecutiveRuns <- consecutiveRUNS.run(
    genotypeFile = ped,
    mapFile = map,
    minSNP = 20,
    maxGap = 500000,
    minLengthBps = 10^6, # 1.000.000
    maxOppRun = 1,
    maxMissRun = 1
)

write.table(consecutiveRuns, paste0(out,"_consecutiveRuns.txt"), 
            quote = F, row.names = F)

# ---- Summary statistics on detected runs ----

summaryList_CR <- summaryRuns(
    runs = consecutiveRuns, mapFile = map, genotypeFile = ped, 
    Class = 2, snpInRuns = TRUE)

# Save Froh
FROH = summaryList_CR$result_Froh_genome_wide

write.table(FROH, paste0(out,"_FROH.txt"), 
            quote = F, row.names = F)

# Save SNP in run
SNPRUN = summaryList_CR$SNPinRun

write.table(SNPRUN, paste0(out,"_SNPRUN.txt"), 
            quote = F, row.names = F)

# topRuns <- tableRuns(
#     runs =  consecutiveRuns, genotypeFile = ped, mapFile = map, 
#     threshold = 0.5)

# ---- Plots ----
plot_Runs(runs = consecutiveRuns, 
          savePlots = T, outputName = out)
# Sliding windows approach captures more small sized ROHs

plot_StackedRuns(runs = consecutiveRuns, 
                 savePlots = T, outputName = out)

plot_SnpsInRuns(
    runs = consecutiveRuns[consecutiveRuns$chrom==1,], 
    genotypeFile = ped, 
    mapFile = map, savePlots = T, 
    outputName = out)

plot_manhattanRuns(
    runs = consecutiveRuns[consecutiveRuns$group=='1',], 
    genotypeFile = ped, mapFile = map,
    savePlots = T, outputName = paste0(out,"_manhattanRuns"))
