library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

folder_name = "75_50_3/rep_01/"

# ---- N animals per generation ----
n_animals = read.table(paste0("simulations/",folder_name,"scenario_03/n_animals.txt"),
                       row.names=NULL)[,c(2:5)]
n_animals = n_animals[,c(1,3,4)]

n_animals = n_animals %>% 
    pivot_longer(cols = c("n_fakes","n_train"),
                 names_to = "class",
                 values_to = "count") %>% 
    mutate(gen = gen-min(gen)) %>% 
    mutate(class = case_when(
        class == "n_fakes" ~ "Hypothetical progenies",
        class == "n_train" ~ "Training population"
    ))

ggplot(data=n_animals, aes(x=gen, y=count, group = class, color = class)) +
    geom_point()+
    geom_smooth(aes(color = class), method = "loess") +
    scale_colour_manual("", 
                        breaks = c("Hypothetical progenies", 
                                   "Training population"),
                        values = c("orange", "purple")) +
    ylab("Number of individuals") +
    xlab("Generation") +
    theme_bw() +
    theme(legend.position = c(0.84, 0.125),
          legend.background = element_rect(fill='transparent'))

ggsave(paste0("plots/",folder_name,"n_animals",".png"),
       width = 15, height = 10, units = "cm", device = "png")

# ---- analyze_df ----
analyze_df = function(folder_name="",
                      scenario){
    print("--------")
    
    # Import data from all scenarios
    df = read.table(paste0("simulations/",folder_name,
                           "scenario_",scenario,"/candidates_metrics.txt"),
                    header = T)
    
    assign(paste0('df_',scenario),df,envir = .GlobalEnv)
    
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

# ---- Line Plots ----
make_lineplot = function(folder_name, trait){
    
    df1 = summ_01[,c('gen', trait)] %>%
        mutate(gen=gen-min(gen))
    names(df1) = c('gen', 'y')
    df2 = summ_02[,c('gen', trait)] %>%
        mutate(gen=gen-min(gen))
    names(df2) = c('gen', 'y')
    df3 = summ_03[,c('gen', trait)] %>%
        mutate(gen=gen-min(gen))
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
        scale_x_continuous(breaks = seq(0,75,5)) +
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


# ---- Scatter plots ----
make_scatter_plot = function(folder_name,
                      trait1,
                      trait2,
                      generation,
                      scenario){
    
    # Import data from all scenarios
    DF = read.table(paste0("simulations/",folder_name,
                           "scenario_",scenario,"/candidates_metrics.txt"),
                    header = T)
    
    DF = DF[,c('ID','gen', trait1, trait2)] %>% 
        filter(!duplicated(ID)) %>% 
        filter(gen == generation) %>% 
        select(trait1, trait2)
    names(DF) = c("trait1","trait2")
    
    g = ggplot(DF, aes(trait1, trait2)) + 
        geom_point() +
        labs(x = trait1, y = trait2,
             title = paste0(trait1," and ",trait2," for generation ",generation),
             subtitle = paste0("Scenario ",scenario))+
        theme_bw()
    ggsave(paste0("plots/",folder_name,"scatter_",
                  trait1,"_",trait2,
                  "_sc_",scenario,"_gen_",generation,".png"),
           width = 20, height = 10, units = "cm", device = "png")
}

# ---- Expected progeny x reality ----
folder_name = "75_50_3/rep_01/"
scenario = "02"

analyze_df(folder_name,"01")
analyze_df(folder_name,"02")
analyze_df(folder_name,"03")

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
generation = max(summ_01$gen)-1

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
    mutate(sire_code = dense_rank(pick(sire, Fg))) # REPLACE TRAIT

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
    geom_point(aes(y = Fg, shape = progeny), size=3) + # REPLACE TRAIT
    scale_shape_manual(values=c(1, 19))+
    scale_x_continuous(name = "Sire", breaks = sire_code$code, labels = sire_code$sire) +
    labs(title = paste0("Real vs Fake progeny scenario ",scenario)) +
    theme_bw()
    
ggsave(paste0("plots/",folder_name,"expectation_vs_reality_",scenario,"_Fg",".png"),
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

