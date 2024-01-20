library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

analyze_df = function(scenario){
    print("--------")
    
    # Import data from all scenarios
    df = read.table(paste0("scenario_",scenario,"/candidates_metrics.txt"),
                    header = T)
    
    pedigree = read.table(paste0("scenario_",scenario,"/pedigree.txt"),
                          header = T)
    
    # Include generation
    df = merge(df, pedigree[,c(1,8)], all.x = T)
    
    # Number of unique IDs
    print(paste0("Number of unique IDs: ", nrow(df[!duplicated(df$ID),])))
    # Number of times IDs repeat at 1, 2, 3, or 4 times
    print(paste0("Number of times IDs repeat: ",table(table(df$ID))))
    
    # Save repeated IDs
    
    repeated_IDs = df %>% 
        select(ID, generation) %>% 
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
        group_by(generation) %>% 
        summarise(n = n(), 
                  across(pheno:Froh_genome,
                         list(mean = mean, 
                              sd = sd), na.rm = TRUE)),
        envir = .GlobalEnv)
}

analyze_df("01")
analyze_df("02")
analyze_df("03")

# Plots
make_lineplot = function(trait){
    df1 = summ_01[,c('generation', trait)] %>% 
        filter(generation > 15)
    names(df1) = c('generation', 'y')
    df2 = summ_02[,c('generation', trait)] %>% 
        filter(generation > 15)
    names(df2) = c('generation', 'y')
    df3 = summ_03[,c('generation', trait)] %>% 
        filter(generation > 15)
    names(df3) = c('generation', 'y')
    
    # Scenario 1
    fig <- plot_ly(data=df1,
                   x = ~generation, 
                   y = ~y,
                   type = 'scatter', mode = 'lines',
                   line = list(color='rgb(0,100,80)'),
                   name = 'EBV + Fped')
    
    # Scenario 2
    fig <- fig %>% add_trace(data=df2,
                             x = ~generation, 
                             y = ~y,
                             type = 'scatter', mode = 'lines',
                             line = list(color='rgb(19, 16, 234)'),
                             name = 'GEBV + Fg')
    
    # Scenario 3
    fig <- fig %>% add_trace(data=df3,
                             x = ~generation, 
                             y = ~y, 
                             type = 'scatter', mode = 'lines',
                             line = list(color='red'),
                             name = 'GEBV + Froh')
    
    fig <- fig %>% layout(title = paste0("Evolution of ",trait))
    
    htmlwidgets::saveWidget(
        widget = fig, #the plotly object
        file = paste0("plots/line_",trait,".html"), #the path & file name
        selfcontained = TRUE #creates a single html file
    )
    
    return(fig)
}

make_lineplot('pheno_mean')
make_lineplot('EBV_mean')
make_lineplot('GV_mean')
make_lineplot('solution_mean')
make_lineplot('Fped_mean')
make_lineplot('Fg_mean')
make_lineplot('Froh_genome_mean')
