library(AlphaSimR)
library(dplyr)

# ---- GBLUP_AlphaSimR ----
GBLUP_AlphaSimR = function(pop){
    # Fit RR-BLUP model for genomic predictions
    ans = RRBLUP(pop, simParam=SP)
    pop = setEBV(pop, ans, simParam=SP)
    
    df = data.frame(ID = pop@id, 
                    sex = pop@sex,
                    pheno = pop@pheno,
                    EBV = pop@ebv, 
                    GV = pop@gv)
    names(df) = c('ID','sex','pheno','EBV','GV')
    return(df)
}