library(AlphaSimR)
library(dplyr)

# ---- makeFakeHaplos ----
makeFakeHaplos = function(dirToSave,
                          pop=NULL, 
                          malePop=NULL,
                          femalePop=NULL){
    
    cli_alert_info("\nStart fake haplotypes\n")
    
    if (!is.null(pop)) {
        haplo_males = pullSnpHaplo(pop[pop@sex == "M"])
        haplo_females = pullSnpHaplo(pop[pop@sex == "F"])
    }else{
        haplo_males = pullSnpHaplo(malePop)
        haplo_females = pullSnpHaplo(femalePop)
    }
    
    # Create all combinations of males and females
    combinations <- expand.grid(m = rownames(haplo_males), 
                                f = rownames(haplo_females),
                                stringsAsFactors = FALSE,
                                KEEP.OUT.ATTRS = FALSE)
    
    # Create progeny IDs
    combinations$ID <- paste0(sub("_.*", "", combinations$m), "_",
                              sub("_.*", "", combinations$f), "_",
                              sub(".*_", "", combinations$m),
                              sub(".*_", "", combinations$f))
    
    # Add corresponding rows and convert to matrix
    # 0 0 -> 0
    # 3 0 -> 3
    # 0 1 -> 1
    # 3 1 -> 4
    
    haplo_males[haplo_males==1] = 3
    
    haplo_matrix = haplo_males[combinations$m, ] + haplo_females[combinations$f, ]
    
    cli_alert_info("\nHaplotype combinations were created.\n")
    
    rownames(haplo_matrix) = combinations$ID
    
    geno = apply(haplo_matrix, 1, function(row) paste(row, collapse = ""))
    
    # Save the genotypes
    write.table(geno,
                paste0(dirToSave,"fakeHaplos_raw.txt"),
                quote=F, row.names = T, 
                col.names = F)
    
    cli_alert_success("\nFake haplotypes created.\n")
    
    # Return pedigree
    return(combinations %>% 
               select(ID,m,f) %>% 
               rename(sire=m,dam=f))
}
