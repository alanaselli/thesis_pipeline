library(AlphaSimR)
library(dplyr)

# ---- add_prefix ----

add_prefix = function(pop, prefix){
    new_IDs = paste0(rep(prefix,nInd(pop)),"_",pop@id)
    return(new_IDs)
}