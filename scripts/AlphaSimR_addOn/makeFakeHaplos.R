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
    
    # Add corresponding rows and convert to matrix (0125 genotype)
    # 0 0 -> 0
    # 1 0 -> 1
    # 1 1 -> 2
    haplo_matrix <- haplo_males[combinations$m, ] + haplo_females[combinations$f, ]
    
    # Create a data frame with ID and the collapsed genotypes
    fakeHaplos = data.frame(ID = combinations$ID,
                            geno = apply(haplo_matrix, 1, 
                                         function(row) paste(row, collapse = "")))
    # Order by ID
    fakeHaplos = fakeHaplos[order(fakeHaplos$ID),]
    
    # Save the genotypes in the appropriate format for BLUPF90
    write.table(fakeHaplos,
                paste0(dirToSave,"fakeHaplos.txt"),
                quote=F, row.names = F, 
                col.names = F, sep = "\t")
    
    cli_alert_success("\nFake haplotypes created.\n")
    
    # Return pedigree
    return(combinations %>% 
               select(ID,m,f) %>% 
               rename(sire=m,dam=f))
}