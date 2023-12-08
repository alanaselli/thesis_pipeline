library(AlphaSimR)
library(dplyr)

# ---- makeRecentPop ----
makeRecentPop = function(
        previous_pop,    
        males_each_year,
        females_each_year,
        nCrosses,
        nProgeny = 1,
        years_of_breeding,
        return_breeding_groups = F,
        return_last_generations = 1,
        year=1,
        rec_data_param = list("file_to_write", 
                              "pop_name", 
                              "append_to_existing_file"),
        addPrefix = NULL
){
    #' @title Generate recent generations.
    #' 
    #' @description Make recent generations with overlapping generations and selection.
    #' 
    #' @param previous_pop A population to sample the first parents from.
    #' @param males_each_year A vector of the number o males to use from generations ago,
    #' starting from the current generation.
    #' @param females_each_year A vector of the number o females to use from generations ago,
    #' starting from the current generation.
    #' @param nCrosses Number of crosses to make in each generation.
    #' @param nProgeny Number of progeny per cross. Default is 1.
    #' @param years_of_breeding The number of generations to be bred, assuming that 1 
    #' generation = 1 year. The minimum number must be 2.
    #' @param return_breeding_groups True returns male and female breeding groups. Default is
    #'  False.
    #' @param return_last_generations Number of generations to be returned. Default is 1.
    #' @param year Starting year of the breeding process. Default = 1.
    #' @param rec_data_param List with parameters to write record file. If Default, no file 
    #' will be written.
    #' @param addPrefix Provide a string to be added at the beginning of the population's
    #' IDs. Default is NULL.
    
    # Check if total number of males and females is less or equal to the number of crosses.
    try(if(females_each_year[1] > nCrosses/2) 
        stop("Not enough crosses to produce the number of female candidates required."))
    
    try(if(males_each_year[1] > nCrosses/2) 
        stop("Not enough crosses to produce the number of male candidates required."))
    
    # Check if the number of males and females is decreasing for older generations
    if (length(males_each_year) > 1) {
        try(if(males_each_year[1]<males_each_year[2]) 
            stop("The number of males in older generations is greater than in younger generations."))
    }
    
    if (length(females_each_year) > 1) {
        try(if(females_each_year[1]<females_each_year[2]) 
            stop("The number of females in older generations is greater than in younger generations."))
    }
    
    # Check if years_of_breeding is greater or equal to 2
    try(if (years_of_breeding < 2)
        stop("Years of breeding must be equal or greater than 2."))
    
    cat(paste0("\nWorking on generation number ",1,"/",years_of_breeding))
    
    # Select first generation of males
    males = selectInd(pop = previous_pop, nInd = sum(males_each_year), use = "pheno", sex = "M")
    
    pos_start = 1
    sires = c()
    list_of_sires = c()
    for (n in 1:length(males_each_year)) {
        pos_end=pos_start+males_each_year[n]-1
        
        siresN = males[pos_start:pos_end]
        
        siresN = setMisc(x = siresN,
                         node = "yearOfBirth",
                         value = year-n)
        
        sires = c(siresN, sires)
        list_of_sires = c(list_of_sires, siresN)
        assign(paste0("sires",n), siresN)
        
        pos_start = pos_end+1
    }
    
    # checking total number of sires
    try(if(sires@nInd!=sum(males_each_year)) 
        stop("Different number of sires from the specified."))
    
    # Checking for overlapped sires
    try(if(isTRUE(duplicated(sires@id)))
        stop("Repeated sires selected."))
    
    # Select first generation of females
    females = selectInd(pop = previous_pop, 
                        nInd = sum(females_each_year), 
                        use = "pheno", sex = "F")
    
    pos_start = 1
    dams = c()
    list_of_dams = c()
    for (n in 1:length(females_each_year)) {
        pos_end=pos_start+females_each_year[n]-1
        
        damsN = females[pos_start:pos_end]
        
        damsN = setMisc(x = damsN,
                        node = "yearOfBirth",
                        value = year-n)
        
        dams = c(damsN, dams)
        list_of_dams = c(list_of_dams,damsN)
        assign(paste0("dams",n), damsN)
        
        pos_start = pos_end+1
    }
    
    # checking total number of dams
    try(if(dams@nInd!=sum(females_each_year)) 
        stop("Different number of dams from the specified."))
    
    # Checking for overlapped dams
    try(if(isTRUE(duplicated(dams@id)))
        stop("Repeated dams selected."))
    
    # Cross first generation of recent generations
    year = year + 1
    recent_pop = randCross2(dams, sires, nCrosses = nCrosses, nProgeny = nProgeny)
    recent_pop = setMisc(x = recent_pop,
                         node = "yearOfBirth",
                         value = year)
    if(!is.null(addPrefix)){
        recent_pop@id = add_prefix(recent_pop,addPrefix)
    }
    if(rec_data_param[1] != "file_to_write"){
        rec_data(rec_data_param[[1]], 
                 recent_pop, 
                 rec_data_param[[2]], 
                 year, 
                 append = rec_data_param[[3]])
    }
    
    # Simulate generations of breeding
    # Define which generation to return
    return_gen = years_of_breeding - return_last_generations + 1
    return_pop_list = list()
    for (gen in (2):(years_of_breeding)) {
        
        year = year+1
        
        cat(paste0("\nWorking on generation number ",gen,"/",years_of_breeding))
        
        sires = c()
        dams = c()
        
        # Select sires
        for (n in length(list_of_sires):1) {
            if (n > 1) {
                siresN = selectInd(list_of_sires[[n-1]], 
                                   nInd = males_each_year[n])
                assign(paste0("sires",n), siresN)
                sires = c(siresN, sires)
                list_of_sires[[n]] = siresN
            } else{
                siresN = selectInd(recent_pop,
                                   nInd = males_each_year[n], 
                                   sex = "M")
                sires = c(siresN, sires)
                list_of_sires[[n]] = siresN
            }
        }
        
        # Select dams    
        for (n in length(list_of_dams):1) {
            if (n > 1) {
                damsN = selectInd(list_of_dams[[n-1]], 
                                  nInd = females_each_year[n])
                assign(paste0("dams",n), damsN)
                dams = c(damsN, dams)
                list_of_dams[[n]] = damsN
            } else{
                damsN = selectInd(recent_pop,
                                  nInd = females_each_year[n], 
                                  sex = "F")
                dams = c(damsN, dams)
                list_of_dams[[n]] = damsN
            }
        }
        
        # Breed
        recent_pop = randCross2(dams, sires, nCrosses = nCrosses, nProgeny = nProgeny)
        recent_pop = setMisc(x = recent_pop,
                             node = "yearOfBirth",
                             value = year)
        if(!is.null(addPrefix)){
            recent_pop@id = add_prefix(recent_pop,addPrefix)
        }
        
        if(rec_data_param[1] != "file_to_write"){
            rec_data(rec_data_param[[1]], 
                     recent_pop, 
                     rec_data_param[[2]], 
                     year, 
                     append = TRUE)
        }
        
        if (gen >= return_gen){
            return_pop_list = c(return_pop_list, assign(as.character(gen), recent_pop), use.names=T)
        }
    }
    cat("\nSimulation of recent generations is complete.")
    
    return(return_pop_list)
}