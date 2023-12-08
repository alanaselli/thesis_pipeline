library(AlphaSimR)
library(dplyr)

# ---- makeFakePed ----
makeFakePed = function(pop) {
    cli_alert_info("\nStart fake pedigree\n")
    males = pop[pop@sex == "M"]
    females = pop[pop@sex == "F"]
    
    fakePed = expand.grid(males@id, females@id,
                          KEEP.OUT.ATTRS=FALSE,
                          stringsAsFactors=FALSE)
    names(fakePed) = c("sire", "dam")
    
    fakePed = fakePed %>% 
        dplyr::mutate(ID = paste0(sire,"_",dam)) %>% 
        select(ID, sire, dam)
    
    cli_alert_success("\nFake pedigree completed.\n")
    
    return(fakePed)
}