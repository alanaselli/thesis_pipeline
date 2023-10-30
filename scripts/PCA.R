library(flashpcaR)
library(dplyr)
library(ggplot2)
library(cli)

pop_names = tools::file_path_sans_ext(list.files(path = "01_genotypes", 
                                                 pattern = "*_100.ped",
                                                 full.names = T))

# Merge populations

merge_file = c()
for (FILE in pop_names) {
  merge_file = c(merge_file, paste0(FILE, ".ped", " 01_genotypes/new_map.map"))
}

write.table(merge_file, "merge_file.txt", row.names = F, col.names = F, quote = F)

system(command = 'plink --cow --merge-list merge_file.txt --make-bed --out 01_genotypes/merged')

plink_path = "01_genotypes/merged"

f3 <- flashpca(X = plink_path)

fam = read.table('01_genotypes/merged.fam', 
                 col.names = c('group','id','V3','V4','sex','pheno'))

# PCA is a method that:
#  Measures how each variable is associated with one another using a Covariance matrix
#  Understands the directions of the spread of our data using Eigenvectors
#  Brings out the relative importance of these directions using Eigenvalues

# a numeric vector. The eigenvalues of X X' / m.
f3$values
write.table(f3$values, '04_PCA/eigenvalues.csv', 
            row.names = F, col.names = F, quote = F, sep = ',')

# a numeric matrix. The eigenvectors of X X' / m.

vectors_ids = cbind(fam, f3$vectors) %>% 
  select(!c(V3,V4)) %>% 
  dplyr::mutate(group=rep(pop_names,each=100))

write.csv(vectors_ids, '04_PCA/eigenvectors.csv', row.names = F, quote = F)

# Plots

ggplot(vectors_ids, aes(vectors_ids[,4], vectors_ids[,5], col = group)) +
  geom_point(size=2, alpha=0.5) + theme_light() +
  xlab(paste0("PC1 (", signif(f3$values[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(f3$values[2], 3), "%)"))
ggsave(paste0('04_PCA/PCA.png'), 
       width = 2000, height = 1000, dpi = 320, units = "px")

cli_alert_success("\nPCA analyses completed.\n")