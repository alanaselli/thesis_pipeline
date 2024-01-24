library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

analyze_df = function(folder_name="",
                      scenario){
    print("--------")
    
    # Import data from all scenarios
    df = read.table(paste0(folder_name,"scenario_",scenario,"/candidates_metrics.txt"),
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
                                       sd = sd), na.rm = TRUE)),
        envir = .GlobalEnv)
    
    # Variation in metrics excluding repeated records
    assign(paste0('summ_',scenario), df %>% 
        filter(!duplicated(ID)) %>% 
        group_by(gen) %>% 
        summarise(n = n(), 
                  across(pheno:Froh_genome,
                         list(mean = mean, 
                              sd = sd), na.rm = TRUE)),
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

dir.create(file.path("plots/", folder_name), showWarnings = FALSE)

for (trait in c(names(summ_01)[c(3,5,7,9,11,13,15)])){
    make_lineplot(folder_name,trait)
}
