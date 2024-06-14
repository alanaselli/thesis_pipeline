library(detectRUNS)

folder_name = "01_genotypes/"

ped = paste0(folder_name,'recent.ped')
map = paste0(folder_name,'new_map.map')
out = paste0("03_ROH/","recent")

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

consecutiveRuns %>% 
    group_by(chrom) %>% 
    mutate(lengthMps = lengthBps/1000000) %>% 
    summarise(max(lengthMps))

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

round(summaryList_CR$summary_ROH_count_chr/98, 2)
summaryList_CR$summary_ROH_percentage_chr
mean(summaryList_CR$result_Froh_chromosome_wide$Chr_1,3)
mean(summaryList_CR$result_Froh_chromosome_wide$Chr_2,3)
mean(summaryList_CR$result_Froh_chromosome_wide$Chr_3,3)
mean(summaryList_CR$result_Froh_chromosome_wide$Chr_4,3)
mean(summaryList_CR$result_Froh_chromosome_wide$Chr_5,3)
summaryList_CR$summary_ROH_count
summaryList_CR$summary_ROH_percentage
mean(summaryList_CR$result_Froh_class$Froh_Class_0)
mean(summaryList_CR$result_Froh_class$Froh_Class_2)
mean(summaryList_CR$result_Froh_class$Froh_Class_4)
sum(summaryList_CR$result_Froh_class$Froh_Class_8, na.rm = T)/length(which(!is.na(summaryList_CR$result_Froh_class$Froh_Class_8)))
sum(summaryList_CR$result_Froh_class$Froh_Class_16, na.rm = T)/length(which(!is.na(summaryList_CR$result_Froh_class$Froh_Class_16)))

summaryList_CR$result_Froh_class$Froh_Class_8[is.na(summaryList_CR$result_Froh_class$Froh_Class_8)] = 0
summaryList_CR$result_Froh_class$Froh_Class_16[is.na(summaryList_CR$result_Froh_class$Froh_Class_16)] = 0
mean(summaryList_CR$result_Froh_class$Froh_Class_8)
mean(summaryList_CR$result_Froh_class$Froh_Class_16)

# ---- Plots ----
plot_Runs(runs = consecutiveRuns, 
          #savePlots = T, outputName = out
          )
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
    savePlots = F, outputName = paste0(out,"_manhattanRuns"))
