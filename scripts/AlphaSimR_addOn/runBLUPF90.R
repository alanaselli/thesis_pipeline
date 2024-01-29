library(cli)
library(processx)

# ---- runBLUPF90 ----
runBLUPF90 = function(param_card, # "renum.txt"
                    genotypes="01_genotypes/recent.ped",
                    plink_format=TRUE, # prepare_snp_file.sh
                    pedigree="01_genotypes/pedigree.txt",
                    min_gen){
    
    # http://nce.ads.uga.edu/wiki/doku.php?id=faq.segfault
    system(command="ulimit -s unlimited")
    system(command="export OMP_STACKSIZE=64M")
    
    setwd("05_BLUPF90/")
    
    # Run RENUMF90
    # invisible(processx::run(command = "renumf90", c(param_card),
    #               wd = "05_BLUPF90/", echo_cmd = T, spinner = T,
    #               stdout = ""))
    cli_alert_info("\nStarting renumf90.\n")
    system(command = "renumf90", input = param_card)

    # Run BLUPF90
    # invisible(processx::run(command = "blupf90", c("renf90.par", "--dense"),
    #               wd = "05_BLUPF90/", echo_cmd = T, spinner = T,
    #               stdout = ""))
    cli_alert_info("\nStarting renumf90.\n")
    system(command = "blupf90 --dense", input = "renf90.par")
    
    setwd("../")
    
    cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}