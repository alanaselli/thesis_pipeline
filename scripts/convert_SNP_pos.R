library(dplyr)
library(cli)

cli_h1("\nInitiating map correction.\n")

n_snps = c(950)
n_snps = c(3000, 3000, 3000, 3000, 3000)
#chr_lengths = c(1e8)
chr_lengths = c(1.58e8, 1.58e8, 1.58e8, 1.58e8, 1.58e8)

# maps_list = tools::file_path_sans_ext(list.files(path = "01_genotypes", 
#                                                  pattern = "*.map",
#                                                  full.names = T))

options(scipen = 999)

plink_map1 = read.table("01_genotypes/founder_sample_100.map", header = F, 
                       col.names = c("chr", "snp", "dist_cM", "dist_bp"))
plink_map1$dist_cM = plink_map1$dist_cM/10

dist_bp = c()
for (n in 1:length(chr_lengths)) {
  new_dist = plink_map1[plink_map1$chr==n,3]*chr_lengths[n]
  dist_bp = c(dist_bp,new_dist)
}

plink_map2=plink_map1
plink_map2['dist_bp'] = floor(dist_bp)
plink_map2['dist_cM'] = round(plink_map2['dist_cM'],4)

write.table(plink_map2, "01_genotypes/new_map.map", quote = F, row.names = F, col.names = F)

cli_alert_success("\nMap correction completed.\n")
