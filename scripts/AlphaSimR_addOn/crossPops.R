library(AlphaSimR)
library(dplyr)

# ---- crossPops ----
crossPops = function(
        pop1,
        pop2,
        nMales,
        nFemales,
        nCrosses,
        nProgeny = 1
){
    
    #' @title Cross populations to form 1/2 descendants.
    #'
    #' @description Cross populations making sure that each parent is being bred
    #' with an individual from the other population.
    #' 
    #' @param pop1 The first population
    #' @param pop2 The second population
    #' @param nMales A vector with the number of males from each population
    #' @param nFemales A vector with the number of females from each population
    #' @param nCrosses A vector with the number of crosses F1xM2 and F2xM1 to be performed
    #' @param nProgeny Number of progeny per cross. Default is 1.
    
    # Select males from pop1
    malesA = selectInd(pop = pop1, nInd = nMales[1], sex = "M", returnPop = FALSE)
    
    # Select females from pop1
    femalesA = selectInd(pop = pop1, nInd = nFemales[1], sex = "F", returnPop = FALSE)
    
    # Select males from pop2
    malesB = selectInd(pop = pop2, nInd = nMales[2], sex = "M", returnPop = FALSE)
    
    # Select females from pop2
    femalesB = selectInd(pop = pop2, nInd = nFemales[2], sex = "F", returnPop = FALSE)
    
    # Make crosses
    cross1 = randCross2(females = pop1, males = pop2, nCrosses = nCrosses[1],
                        femaleParents = femalesA, maleParents = malesB,
                        nProgeny = nProgeny)
    cross2 = randCross2(females = pop2, males = pop1, nCrosses = nCrosses[2],
                        femaleParents = femalesB, maleParents = malesA,
                        nProgeny = nProgeny)
    
    # Merge crosses
    crossed_pop = c(cross1, cross2)
    
    return(crossed_pop)
}