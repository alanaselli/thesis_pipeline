library(AlphaSimR)
library(dplyr)

# ---- SelectionWrapper2 ----
SelectionWrapper2 = function(InitialPop,
                             InitialGenotypes="01_genotypes/recent.ped",
                             nInd, # c(nMales,nFemales)
                             nCrosses,
                             nGenerations,
                             SelectionMethod,
                             year,
                             path_Genotypes="01_genotypes/", 
                             pathScenario="scenario_01/",
                             OriginalRecords="01_genotypes/pedigree.txt",
                             scenarioName="01"
){
    # The wrapper to be used with selectCandidates2
    
    # Create a new pedigree file for the new scenario
    system(paste0("awk '{ print $1, $2, $3 }' ", OriginalRecords, " > ", 
                  pathScenario, "pedigree.txt"))
    
    # 
    
    # Select candidates from the initial population
    initialParents = selectCandidates2(pop=InitialPop, 
                                       append=FALSE, 
                                       method=SelectionMethod,
                                       top_ebv=nInd)
    
    # Perform matings
    candidates = randCross(initialParents, nCrosses = nCrosses)
    
    # Save pedigree data
    cli_alert_info(paste0("\nRecording data for gen 1 of ",
                          scenarioName,".\n"))
    rec_data(pathAllRecords, candidates, scenarioName, year, append = TRUE)
    
    # Save genotypes
    writePlink(candidates,
               paste0(path_Genotypes,"sc_",scenarioName,"_gen_1"))
    
    # Continue this process for N generations
    for (i in 2:nGenerations+1) {
        cli_alert_info(paste0("\nInitiating gen ",i,
                              " of scenario ", scenarioName,".\n"))
        year=year+1
        cli_alert(paste0("Year: ",year))
        pedParents = paste0("sc_",scenarioName,"_gen_",i-1)
        
        # Run BLUPF90
        runBLUPF90(paste0(pedParents,".ped"),
                   min_gen=year-10)
        
        # Select candidates
        parents = selectCandidates(pop=candidates, 
                                   file_name=pathNewRecords, 
                                   append=TRUE, 
                                   method=SelectionMethod,
                                   top_ebv=nInd)
        
        # Perform matings
        candidates = randCross(parents, nCrosses = nCrosses)
        
        # Save pedigree data
        rec_data(pathAllRecords, candidates, scenarioName, year, append = TRUE)
        
        # Save genotypes
        ped_name = paste0("sc_",scenarioName,"_gen_",i)
        writePlink(candidates, paste0(geno_path,ped_name))
    }
    cli_alert_success(paste0("/nSelection process for scenario ",
                             scenarioName, " is complete./n"))
}