# ---- import packages ----
library(dplyr)
library(stringr)
library(ggplot2)
library(cli)

cli_h1("\nInitiating LD analyses.\n")

pop_names = tools::file_path_sans_ext(list.files(path = "01_genotypes/", 
                                                 pattern = "*_sample_100.ped"))
cli_alert_info("\nPopulations found:")
for (pop in pop_names) {
  cat("\n",pop)
}

cli_h2("\nCalculating LD for each population.\n")
# ---- Run PLINK ----
for (FILE in pop_names){
  out = paste0('02_LD/',FILE)
  
  command = paste0('plink --ped 01_genotypes/', FILE, '.ped',
                   ' --map 01_genotypes/new_map.map',
                   ' --make-bed --cow --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 1000000 --out ', 
                   out)
  system(command = command)
}
LD_files = tools::file_path_sans_ext(list.files(path = "02_LD/", pattern = "*.ld"))
cli_alert_info("\nLD files found:\n")
for (f in LD_files) {
  cat(f,"\n")
}

for (FILE in pop_names) {
  cli_h2(paste0("\nInitiating summarization for population ", FILE, ".\n"))
  df = read.table(paste0('02_LD/',FILE,".ld"), header = T)
  
  # ---- adjacent LD metrics ----
  
  adjacent_r2 = df %>% 
    filter(CHR_A==CHR_B) %>%
    group_by(CHR_A) %>% 
    arrange(BP_A, BP_B) %>% 
    filter(!duplicated(SNP_A)) %>% 
    group_by(CHR_A) %>% 
    summarise(size_MB = round(max(BP_B)/10^6,1), 
              SNP_n = n()+1, 
              mean_distance_kb = round(mean(BP_B-BP_A)/10^3,1),
              SD_distance_kb = round(sd(BP_B-BP_A)/10^3,2),
              mean_r2 = round(mean(R2),2),
              SD_r2 = round(sd(R2),2),
              median_r2 = round(median(R2),2))
  
  write.csv(adjacent_r2, 
            paste0('02_LD/',FILE,"_table_01.csv"),
            row.names = F, quote = F)
  
  # ---- autossomic r2 ----
  for (chromosome in unique(df$CHR_A)) {
    autossomic_r2 = df %>% 
      filter(CHR_A==chromosome & CHR_B==chromosome) %>% 
      dplyr::mutate(dist = (BP_B-BP_A)/1000) %>% 
      select(dist, R2) %>% 
      dplyr::mutate(distance_kb = case_when(
        dist <= 50 ~ '0-50',
        50 < dist & dist <= 100 ~ '050-100',
        100 < dist & dist <= 200 ~ '100-200',
        200 < dist & dist <= 300 ~ '200-300',
        300 < dist & dist <= 400 ~ '300-400',
        400 < dist & dist <= 500 ~ '400-500',
        500 < dist & dist <= 1000 ~ '500-1000'
      )) %>% 
      filter(distance_kb!='NA') %>% 
      select(distance_kb, R2) %>% 
      group_by(distance_kb) %>% 
      summarise(N = n(),
                mean_r2 = round(mean(R2),2),
                SD_r2 = round(sd(R2),2),
                median_r2 = round(median(R2),2),
                perc_r2_gt_01 = round(length(R2[R2>0.1])/n(),2),
                perc_r2_gt_03 = round(length(R2[R2>0.3])/n(),2),
                perc_r2_gt_05 = round(length(R2[R2>0.5])/n(),2)) %>% 
      arrange(distance_kb)
    
    write.csv(autossomic_r2, 
              paste0('02_LD/',FILE,'_chr',chromosome,".csv"),
              row.names = F, quote = F)
  }
  
  # all chromosomes
  autossomic_r2 = df %>% 
    filter(CHR_A==CHR_B) %>% 
    dplyr::mutate(dist = (BP_B-BP_A)/1000) %>% 
    select(dist, R2) %>% 
    dplyr::mutate(distance_kb = case_when(
      dist <= 50 ~ '0-50',
      50 < dist & dist <= 100 ~ '050-100',
      100 < dist & dist <= 200 ~ '100-200',
      200 < dist & dist <= 300 ~ '200-300',
      300 < dist & dist <= 400 ~ '300-400',
      400 < dist & dist <= 500 ~ '400-500',
      500 < dist & dist <= 1000 ~ '500-1000'
    )) %>% 
    filter(distance_kb!='NA') %>% 
    select(distance_kb, R2) %>% 
    group_by(distance_kb) %>% 
    summarise(N = n(),
              mean_r2 = round(mean(R2),2),
              SD_r2 = round(sd(R2),2),
              median_r2 = round(median(R2),2),
              perc_r2_gt_01 = round(length(R2[R2>0.1])/n(),2),
              perc_r2_gt_03 = round(length(R2[R2>0.3])/n(),2),
              perc_r2_gt_05 = round(length(R2[R2>0.5])/n(),2)) %>% 
    arrange(distance_kb)
  
  write.csv(autossomic_r2, 
            paste0('02_LD/',FILE,'_all_chr.csv'),
            row.names = F, quote = F)
  
  
  # ---- Plot LD 100kb ----
  
  # read LD summary FILE
  for (CHR in unique(df$CHR_A)) {
    dfr = df %>% 
      filter(CHR_A==CHR) %>% 
      dplyr::mutate(dist = (BP_B-BP_A)/1000) %>% 
      select(dist, R2)
    
    # group distances into 10 kb intervals
    dfr$distc_a <- cut(dfr$dist,
                       breaks = seq(from=0,to=100,by=1),
                       labels = seq(from=1,to=100,by=1), 
                       include.lowest = TRUE)
    
    dfr$distc_b <- cut(dfr$dist,
                       breaks = seq(from=100,to=1000,by=10),
                       labels = seq(from=105,to=995,by=10), 
                       include.lowest = TRUE)
    
    dfr$distc_c <- cut(dfr$dist,
                       breaks = seq(from=1000,to=10000,by=100),
                       labels = seq(from=1050,to=9950,by=100), 
                       include.lowest = TRUE)
    
    dfr$distc_d <- cut(dfr$dist,
                       breaks = seq(from=0,to=10000,by=100),
                       labels = seq(from=50,to=9950,by=100), 
                       include.lowest = TRUE)
    
    # compute mean r2 within the blocks
    dfr1 <- dfr %>% 
      group_by(distc_a) %>% 
      summarise(mean=mean(R2))
    
    dfr2 <- dfr %>% 
      group_by(distc_b) %>% 
      summarise(mean=mean(R2))
    
    dfr3 <- dfr %>% 
      group_by(distc_c) %>% 
      summarise(mean=mean(R2))
    
    dfr4 <- dfr %>% 
      group_by(distc_d) %>% 
      summarise(mean=mean(R2))
    
    # plot a
    ggplot()+
      geom_point(data=dfr1,
                 aes(x=as.integer(as.character(distc_a)),y=mean),
                 size=0.4,colour="grey20")+
      geom_line(data=dfr1,
                aes(x=as.integer(as.character(distc_a)),y=mean),
                size=0.3,alpha=0.5,colour="grey40")+
      labs(x="Distance (kb)",y=expression(LD~(r^{2})),
           title = paste0('Linkage Disequilibrium Decay on chromosome ', CHR),
           subtitle = paste0("Population ", FILE)
      )+
      ylim(0,1)+
      scale_x_continuous(breaks = seq(0,100,by=10),
                         limits = c(0,100)) +
      theme_bw()
    
    ggsave(paste0('02_LD/',FILE, "_chr",CHR,"_a.png"), 
           width = 2400, height = 1000, dpi = 320, units = "px")
    
    # plot b
    ggplot()+
      geom_point(data=dfr2,
                 aes(x=as.integer(as.character(distc_b)),y=mean),
                 size=0.4,colour="grey20")+
      geom_line(data=dfr2,
                aes(x=as.integer(as.character(distc_b)),y=mean),
                size=0.3,alpha=0.5,colour="grey40")+
      labs(x="Distance (kb)",y=expression(LD~(r^{2})),
           title = paste0('Linkage Disequilibrium Decay on chromosome ', CHR),
           subtitle = paste0("Population ", FILE)
      )+
      ylim(0,0.2)+
      scale_x_continuous(breaks = seq(100,1000,by=100),
                         limits = c(100,1000)) +
      theme_bw()
    
    ggsave(paste0('02_LD/',FILE, "_chr",CHR,"_b.png"), 
           width = 2400, height = 1000, dpi = 320, units = "px")
    
    # plot c
    ggplot()+
      geom_point(data=dfr3,
                 aes(x=as.integer(as.character(distc_c)),y=mean),
                 size=0.4,colour="grey20")+
      geom_line(data=dfr3,
                aes(x=as.integer(as.character(distc_c)),y=mean),
                size=0.3,alpha=0.5,colour="grey40")+
      labs(x="Distance (kb)",y=expression(LD~(r^{2})),
           title = paste0('Linkage Disequilibrium Decay on chromosome ', CHR),
           subtitle = paste0("Population ", FILE)
      )+
      scale_x_continuous(breaks = seq(1000,10000,by=1000), 
                         limits = c(1000,10000)) +
      scale_y_continuous(breaks = seq(0,0.12,by=0.02),
                         limits = c(0,0.12)) +
      theme_bw()
    
    ggsave(paste0('02_LD/',FILE, "_chr",CHR,"_c.png"), 
           width = 2400, height = 1000, dpi = 320, units = "px")
    
    # plot d
    ggplot()+
      geom_point(data=dfr4,
                 aes(x=as.integer(as.character(distc_d)),y=mean),
                 size=0.4,colour="grey20")+
      geom_line(data=dfr4,
                aes(x=as.integer(as.character(distc_d)),y=mean),
                size=0.3,alpha=0.5,colour="grey40")+
      labs(x="Distância (kb)",y=expression(LD~(r^{2})),
           title = paste0('Linkage Disequilibrium Decay on chromosome ', CHR),
           subtitle = paste0("Population ", FILE))+
      ylim(0,0.4)+
      scale_x_continuous(breaks = c(100, seq(1000, 10000, by=1000)), 
                         limits = c(0,10000)) +
      theme_bw()
    
    ggsave(paste0('02_LD/',FILE, "_chr",CHR,"_d.png"), 
           width = 2400, height = 1000, dpi = 320, units = "px")
    
  }
  
}

system(command = 'rm 02_LD/*.bed 02_LD/*.bim 02_LD/*.fam 02_LD/*.log')

cli_alert_success("\nLD analyses completed.\n")
