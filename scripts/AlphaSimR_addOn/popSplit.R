library(AlphaSimR)
library(dplyr)

# ---- popSplit ----
# Splits a population into 2
popSplit = function(pop_to_split, nMales_each_pop, nFemales_each_pop){
    
    #' @title Split a population.
    #' 
    #' @description Splits a population into 2, ensuring that each population
    #' derives from different parents.
    #' 
    #' @param pop_to_split is an AlphaSimR population.
    #' @param nMales_each_pop is a scalar of the number of males in each new population.
    #' @param nFemales_each_pop is a scalar of the number of females in each new population. 
    
    cat("\nStarting population split")
    
    cat("\nSampling males")
    randSample_males = selectInd(pop = pop_to_split,
                                 nInd = nMales_each_pop*2,
                                 use = "rand",
                                 sex = "M")
    
    cat("\nSampling females")
    randSample_females = selectInd(pop = pop_to_split,
                                   nInd = nFemales_each_pop*2,
                                   use = "rand",
                                   sex = "F")
    
    cat("\nCreating male populations")
    males_pop1 = randSample_males[1:nMales_each_pop]
    males_pop2 = randSample_males[(nMales_each_pop+1):(nMales_each_pop*2)]
    
    cat("\nCreating female populations")
    females_pop1 = randSample_females[1:nFemales_each_pop]
    females_pop2 = randSample_females[(nFemales_each_pop+1):(nFemales_each_pop*2)]
    
    cat("\nMerging population 1")
    Pop1 = mergePops(list(males_pop1,females_pop1))
    cat("\nMerging population 2")
    Pop2 = mergePops(list(males_pop2,females_pop2))
    
    # Testing
    if (length(intersect(Pop1@id, Pop2@id)) != 0) {
        cat("\nWARNING: There are overlaping ids!")
    }
    
    return(list("Pop1"=Pop1,"Pop2"=Pop2))
}