library(AlphaSimR)
library(dplyr)

# ---- rec_data ----
# Write general data from simulation

rec_data = function(file, pop, pop_name, gen, append = TRUE) {
    
    ped_rec = cbind(pop@id, pop@father, pop@mother, pop_name,
                    pop@sex, pop@gv, pop@pheno, gen)
    colnames(ped_rec) = c("ID", "sire", "dam", "population",
                          "sex", "gv", "pheno", "generation")
    if (isTRUE(append)){
        write.table(ped_rec, file, append = TRUE, col.names = F,
                    quote = F, row.names = F)
    }
    else {
        write.table(ped_rec, file, append = FALSE, quote = F, row.names = F)
    }
}