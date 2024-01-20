library(cli)

# ---- runBLUPF90 ----
runBLUPF90 = function(param_card, # "renum.txt"
                    genotypes="01_genotypes/recent.ped",
                    plink_format=TRUE, # prepare_snp_file.sh
                    pedigree="01_genotypes/pedigree.txt",
                    min_gen){

    # Prepare files for RENUMF90
    setwd("05_BLUPF90/")

    # Run RENUMF90
    cli_alert_info("\nRunning RENUMF90.\n")

    system(command = "renumf90", input = param_card)

    # Run BLUPF90
    cli_alert_info("\nRunning BLUPF90.\n")
    system(command = "blupf90 --dense", input = "renf90.par")

    setwd("../")
    cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}