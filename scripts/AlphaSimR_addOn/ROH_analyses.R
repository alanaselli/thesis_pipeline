library(detectRUNS)

ROH_analyses = function(pop_name,
                        generation,
                        scenario,
                        ped,
                        map="01_genotypes/new_map.map",
                        save_to="03_ROH/",
                        save_metrics=FALSE){
    
    cli_alert_info(paste0("\nStarting ROH analysis for ",pop_name,".\n"))
    
    out = paste0(save_to,scenario,"_",generation)
    
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
    if (isTRUE(save_metrics)) {
        write.table(consecutiveRuns, paste0(out,"_consecutiveRuns.txt"), 
                    quote = F, row.names = F)
    }
    
    # ---- Summary statistics on detected runs ----
    
    summaryList_CR <- summaryRuns(
        runs = consecutiveRuns, mapFile = map, genotypeFile = ped, 
        Class = 2, snpInRuns = TRUE)
    
    # Save Froh
    FROH = summaryList_CR$result_Froh_genome_wide
    
    if (isTRUE(save_metrics)) {
        write.table(FROH, paste0(out,"_FROH.txt"), 
                    quote = F, row.names = F)
    }
    
    # Save SNP in run
    SNPRUN = summaryList_CR$SNPinRun
    
    if (isTRUE(save_metrics)) {
        write.table(SNPRUN, paste0(out,"_SNPinRun.txt"), 
                    quote = F, row.names = F)
    }
    
    cli_alert_success("\nROH analyses completed.\n")
    
    return(FROH)
}
