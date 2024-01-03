library(AlphaSimR)
library(dplyr)

# ---- makeFakeHaplos ----
makeFakeHaplos = function(dirToSave,
                          pop=NULL, 
                          malePop=NULL,
                          femalePop=NULL,
                          BLUPF90_format=FALSE){
    
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
    
    if (isFALSE(BLUPF90_format)) {
        cli_alert_info("\nPreparing plink file.\n")
        # Create plink format genotypes
        # 0 0 -> 1 1
        # 0 1 -> 1 2
        # 1 0 -> 2 1
        # 1 1 -> 2 2
        # genotypes = matrix(ncol=ncol(haplo_males),
        #                    nrow=nrow(combinations))
        # 
        # for (row_index in 1:nrow(combinations)) {
        #     genotypes[row_index,]=paste0(haplo_males[combinations$m[row_index], ], 
        #                                  haplo_females[combinations$f[row_index], ])
        # }
        # 
        # genotypes = data.frame(ID = combinations$ID,
        #                        geno = apply(genotypes, 1, 
        #                                     function(row) paste(row, collapse = "")))
        
        geno = apply(combinations, 1, function(row) {
            paste0(haplo_males[row[["m"]], ], haplo_females[row[["f"]], ])
        })
        
        geno = t(geno)
        
        genotypes = data.frame(ID = combinations$ID,
                               geno = apply(geno, 1,
                                            function(row) paste(row, collapse = "")))
        
        write.table(genotypes,
                    paste0(dirToSave,"fakeGenotypes.txt"),
                    col.names = FALSE, quote = FALSE, row.names = FALSE)
        
        cli_alert_info("\nGenotypes successfuly created.\n")
        
        system(command = paste0("scripts/genotypes_to_plink.sh ",
                                dirToSave,"fakeGenotypes.txt ",
                                dirToSave,"fakeGenotypes.ped"))
        cli_alert_success("\nGenotypes successfuly written in plink format.\n")
    } else {
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
    }
    
    # Return pedigree
    return(combinations %>% 
               select(ID,m,f) %>% 
               rename(sire=m,dam=f))
}
