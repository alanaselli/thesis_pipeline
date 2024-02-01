library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

folder_name = "50_50/"

analyze_df = function(folder_name="",
                      scenario){
    print("--------")
    
    # Import data from all scenarios
    df = read.table(paste0("simulations/",folder_name,
                           "scenario_",scenario,"/candidates_metrics.txt"),
                    header = T)
    
    print(summary(df))
    
    # Number of unique IDs
    print(paste0("Number of unique IDs: ", nrow(df[!duplicated(df$ID),])))
    # Number of times IDs repeat at 1, 2, 3, or 4 times
    print(paste0("Number of times IDs repeat: ",table(table(df$ID))))
    
    # Save repeated IDs
    
    repeated_IDs = df %>% 
        select(ID, gen) %>% 
        group_by(ID) %>% 
        summarise(n = n()) %>% 
        filter(n > 1)
    
    # Variation in metrics for the same animal
    assign(paste0('summ_ind_',scenario),df %>% 
        filter(ID %in% repeated_IDs$ID) %>% 
        group_by(ID) %>% 
        summarise(n = n(), across(pheno:Froh_genome,  
                                  list(mean = mean, 
                                       sd = sd))),
        envir = .GlobalEnv)
    
    # Variation in metrics excluding repeated records
    assign(paste0('summ_',scenario), df %>% 
        filter(!duplicated(ID)) %>% 
        group_by(gen) %>% 
        summarise(n = n(), 
                  across(pheno:Froh_genome,
                         list(mean = mean, 
                              sd = sd))),
        envir = .GlobalEnv)
}

analyze_df(folder_name,"01")
analyze_df(folder_name,"02")
analyze_df(folder_name,"03")

# Plots
make_lineplot = function(folder_name, trait){
    df1 = summ_01[,c('gen', trait)] %>% 
        filter(gen > 15)
    names(df1) = c('gen', 'y')
    df2 = summ_02[,c('gen', trait)] %>% 
        filter(gen > 15)
    names(df2) = c('gen', 'y')
    df3 = summ_03[,c('gen', trait)] %>% 
        filter(gen > 15)
    names(df3) = c('gen', 'y')
    
    g = ggplot(df1, aes(gen)) + 
        geom_point(aes(y=y, color = "EBV + Fped")) +
        geom_smooth(aes(y=y, color = "EBV + Fped"),
                    formula = y ~ x, method = "loess") +
        geom_point(aes(y=df2$y, color = "GEBV + Fg")) +
        geom_smooth(aes(y=df2$y, color = "GEBV + Fg"),
                    formula = y ~ x, method = "loess") +
        geom_point(aes(y=df3$y, color = "GEBV + Froh")) +
        geom_smooth(aes(y=df3$y, color = "GEBV + Froh"),
                    formula = y ~ x, method = "loess") +
        scale_colour_manual("", 
                            breaks = c("EBV + Fped", "GEBV + Fg", "GEBV + Froh"),
                            values = c("red", "green", "blue")) +
        labs(title = paste0("Evolution of ",trait)) +
        ylab(trait) +
        theme_bw()
    
    ggsave(paste0("plots/",folder_name,trait,".png"),
           width = 20, height = 10, units = "cm", device = "png")
    
    return(g)
}

dir.create(file.path("plots", folder_name), showWarnings = FALSE)

for (trait in c(names(summ_01)[c(3,5,7,9,11,13,15)])){
    make_lineplot(folder_name,trait)
}

# ---- Expected progeny x reality ----
folder_name = "50_50/"
scenario = "02"
df_real = read.table(paste0("simulations/",folder_name,
                       "scenario_",scenario,"/candidates_metrics.txt"),
                header = T)
df_fake = read.table(paste0("simulations/",folder_name,
                            "scenario_",scenario,"/fakeProgeny.txt"),
                     header = T)
pedigree = read.table(paste0("simulations/",folder_name,
                                       "scenario_",scenario,"/pedigree.txt"),
                                header = T)

# Choose one generation
generation = 65

parents = df_real[df_real$gen == generation,]

progeny_pedigree = pedigree[pedigree$generation == generation,]
progeny_pedigree$sire_dam = paste0(progeny_pedigree$sire, "_", progeny_pedigree$dam)

real_progeny_IDs = progeny_pedigree[progeny_pedigree$sire %in% parents$ID & progeny_pedigree$dam %in% parents$ID, 'ID']
real_progeny = df_real[df_real$gen == generation+1 & df_real$ID %in% real_progeny_IDs,]
real_progeny = merge(real_progeny, pedigree[,c('ID','sire','dam')], 
                     by = "ID", all.x = T)
real_progeny$sire_dam = paste0(real_progeny$sire,"_",real_progeny$dam)

real_progeny = real_progeny %>% 
    mutate(sire_dam_code = dense_rank(sire_dam))

real_progeny = real_progeny %>% 
    mutate(sire_code = dense_rank(pick(sire, Froh_genome))) # REPLACE TRAIT

parents = parents[parents$ID %in% real_progeny$sire | parents$ID %in% real_progeny$dam,]

df_fake$sire_dam = paste0(df_fake$sire,"_",df_fake$dam)
fake_progeny = df_fake[df_fake$sire_dam %in% real_progeny$sire_dam & df_fake$gen == generation,]

table(fake_progeny$sire)/4
table(real_progeny$sire)

real_progeny$progeny = "Real"
fake_progeny$progeny = "Fake"

all_progeny = rbind(real_progeny[,names(real_progeny) %in% names(fake_progeny)], 
                    fake_progeny[,names(fake_progeny) %in% names(real_progeny)])
all_progeny = (all_progeny) %>% 
    mutate(sire_dam_code = dense_rank(sire_dam))

all_progeny = merge(all_progeny, real_progeny[,c('sire_dam_code','sire_code')],
                    by='sire_dam_code')

sire_code = all_progeny %>% 
    group_by(sire) %>% 
    summarise(code = median(sire_code))

all_progeny = merge(all_progeny, sire_code, by = "sire")

all_progeny$sire = as.factor(all_progeny$sire)

g = ggplot(all_progeny, aes(x = sire_code, color = sire)) +
    geom_point(aes(y = Froh_genome, shape = progeny), size=3) + # REPLACE TRAIT
    scale_shape_manual(values=c(1, 19))+
    scale_x_continuous(name = "Sire", breaks = sire_code$code, labels = sire_code$sire) +
    labs(title = paste0("Real vs Fake progeny scenario ",scenario)) +
    theme_bw()
    
ggsave(paste0("plots/",folder_name,"expectation_vs_reality_",scenario,"_Froh",".png"),
       width = 20, height = 10, units = "cm", device = "png")

# ---- Violin plots of generations ----
folder_name = "50_50/"
scenario = "03"

df_real = read.table(paste0("simulations/",folder_name,
                            "scenario_",scenario,"/candidates_metrics.txt"),
                     header = T)

df_real$gen = as.factor(df_real$gen)

g = ggplot(df_real, aes(y = Fg, x = gen)) +
    geom_violin() +
    labs(x="Generation", 
         title = paste0("Violin plot for scenario ",scenario)) +
    theme_bw()

ggsave(paste0("plots/",folder_name,"violin_",scenario,"_Froh.png"),
       width = 20, height = 10, units = "cm", device = "png")

