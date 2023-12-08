library(AlphaSimR)
library(dplyr)

# ---- getPhasedHaplo ----
# Convert haplotypes to phased genotype

getPhasedHaplo = function(pop){
    
    haplotypes_F = pullSnpHaplo(pop, haplo = 1)
    haplotypes_M = pullSnpHaplo(pop, haplo = 2)
    
    phased = matrix(data = NA, 
                    nrow = nrow(haplotypes_F),
                    ncol = ncol(haplotypes_F))
    
    for (n in 1:nrow(haplotypes_F)) {   # n is the individual (row)
        for (m in 1:ncol(haplotypes_F)){  # m is the snp (col)
            phased[n,m] = 
                ifelse(haplotypes_F[n,m] == 0 & haplotypes_M[n,m] == 0, 0,
                       ifelse(haplotypes_F[n,m] == 0 & haplotypes_M[n,m] == 1, 3,
                              ifelse(haplotypes_F[n,m] == 1 & haplotypes_M[n,m] == 0, 4,
                                     ifelse(haplotypes_F[n,m] == 1 & haplotypes_M[n,m] == 1, 2, 1)
                              )
                       )
                )
        }
    }
    rownames(phased) = sub('\\_1$', '', rownames(haplotypes_F))
    
    phased_collapsed = apply(phased, 1, paste, collapse = "")
    
    return(phased_collapsed)
}