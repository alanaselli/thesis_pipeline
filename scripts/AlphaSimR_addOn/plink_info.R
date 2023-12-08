library(AlphaSimR)
library(dplyr)

# ---- plink_info ----
# Create first 6 columns of plink.ped

plink_info = function(file, pop, pop_name, append = TRUE) {
    ped_rec = cbind(pop_name, pop@id, pop@father, pop@mother,
                    pop@sex, pop@pheno)
    colnames(ped_rec) = c("population", "ID", "sire", "dam", "sex", "pheno")
    if (isTRUE(append)){
        write.table(ped_rec, file, append = TRUE, col.names = F,
                    quote = F, row.names = F)
    }
    else {
        write.table(ped_rec, file, append = FALSE, col.names = F, quote = F, row.names = F)
    }
}