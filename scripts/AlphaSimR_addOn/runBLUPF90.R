library(cli)

# ---- runBLUPF90 ----
runBLUPF90 = function(param_card, # "renum.txt"
                    genotypes="01_genotypes/recent.ped",
                    plink_format=TRUE, # prepare_snp_file.sh
                    pedigree="01_genotypes/pedigree.txt",
                    min_gen){

    # Prepare files for RENUMF90
    setwd("05_BLUPF90/")

    # cli_alert_info("\nInitiating data preparation for BLUPF90.\n")
    # extract_ped = paste0("./extract_ped.sh ",
    #                         "../",pedigree,
    #                         " dat1.txt ",
    #                         "ped1.txt ",
    #                         min_gen)
    # system(command = extract_ped)

    # if (isTRUE(plink_format)){
    #     prepare_snp_file = paste0("./prepare_snp_file.sh ",
    #                                 paste0("../",genotypes," "),
    #                                 "../01_genotypes/new_map.map ",
    #                                 "snp_file.txt")
    #     system(command = prepare_snp_file)
    # } else {
    #     cli_alert_info("No genotypic file was provided.")
    # }

    # Run RENUMF90
    cli_alert_info("\nRunning RENUMF90.\n")

    system(command = "renumf90", input = param_card)

    # Run BLUPF90
    cli_alert_info("\nRunning BLUPF90.\n")
    system(command = "blupf90", input = "renf90.par")

    setwd("../")
    cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}