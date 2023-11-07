library(detectRUNS)
library(cli)

pop_names = tools::file_path_sans_ext(list.files(path = "01_genotypes/", pattern = "*_100.ped"))

pop_names = c(pop_names, "recent")
#pop_names = c("recent_MAF")

for (FILE in pop_names) {
    
    cli_alert_info(paste0("\nStarting ROH analysis for ", FILE,"\n"))
    
    ped = paste0('01_genotypes/',FILE,'.ped')
    map = "01_genotypes/new_map.map"
    #map = paste0('01_genotypes/',FILE,'.map')
    out = paste0('03_ROH/',FILE)
    
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
}

cli_alert_success("\nROH analyses completed.\n")
