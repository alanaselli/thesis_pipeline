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
    invisible(processx::run(command = "renumf90", c(param_card),
                  wd = "05_BLUPF90/", echo_cmd = T, spinner = T,
                  stdout = ""))

    # Run BLUPF90
    invisible(processx::run(command = "blupf90", c("renf90.par", "--dense"),
                  wd = "05_BLUPF90/", echo_cmd = T, spinner = T,
                  stdout = ""))
    
    cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}