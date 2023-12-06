library(AlphaSimR)
library(dplyr)

# ---- new_dir_date ----
# Create a new directory with today's date

new_dir_date = function(prefix = "sim"){
  # Create directory's name
  newdir <- paste(prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
  
  # Create new directory
  if (!dir.exists(newdir)) {        # Check if folder already exists
    dir.create(newdir)              # If not, create folder
    setwd(newdir)                   # Set to new directory
    print(paste("Directory", newdir, "successfully created"))
  } else{
    print("Directory already exists")
  }
}

# ---- new_subdir ----
# Create subfolders

new_subdir = function(subdir_names){  # Provide sub directory names in a vector
  for (subdir in subdir_names) {      # For each name provided
    if (!dir.exists(subdir)) {        # Check if folder already exists
      dir.create(subdir)              # If not, create folder
      print(paste("Directory", subdir, "successfully created"))
    } else{
      print("Directory already exists")
    }
  }
}

# ---- rec_data ----
# Write general data from simulation

rec_data = function(file, pop, pop_name, gen, append = TRUE) {
  
  ped_rec = cbind(pop@id, pop@father, pop@mother, pop_name,
                  pop@sex, pop@gv, pop@pheno, gen)
  colnames(ped_rec) = c("ID", "sire", "dam", "population",
                        "sex", "gv", "pheno", "generation")
  if (isTRUE(append)){
    write.table(ped_rec, file, append = TRUE, col.names = F,
                quote = F, row.names = F)
  }
  else {
    write.table(ped_rec, file, append = FALSE, quote = F, row.names = F)
  }
}

# ---- plink_info ----
# Create first 6 columns of plink.ped

plink_info = function(file, pop, pop_name, append = TRUE) {
  ped_rec = cbind(pop_name, pop@id, pop@father, pop@mother,
                  pop@sex, pop@pheno)
  colnames(ped_rec) = c("population", "ID", "sire", "dam", "sex", "pheno")
  if (isTRUE(append)){
    write.table(ped_rec, file, append = TRUE, col.names = F,
                quote = F, row.names = F)
  }
  else {
    write.table(ped_rec, file, append = FALSE, col.names = F, quote = F, row.names = F)
  }
}

# ---- adjust_pos ----
# Adjust positions on SNPmap for different chromosome lengths

adjust_pos = function(map_name, n_snps, chr_lengths){
  cli_h1("\nInitiating map correction.\n")
  
  options(scipen = 999)
  
  plink_map1 = read.table(map_name, header = F, 
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
}

# ---- getPhasedHaplo ----
# Convert haplotypes to phased genotype

getPhasedHaplo = function(pop){
  
  haplotypes_F = pullSnpHaplo(pop, haplo = 1)
  haplotypes_M = pullSnpHaplo(pop, haplo = 2)
  
  phased = matrix(data = NA, 
                  nrow = nrow(haplotypes_F),
                  ncol = ncol(haplotypes_F))
  
  for (n in 1:nrow(haplotypes_F)) {   # n is the individual (row)
    for (m in 1:ncol(haplotypes_F)){  # m is the snp (col)
      phased[n,m] = 
        ifelse(haplotypes_F[n,m] == 0 & haplotypes_M[n,m] == 0, 0,
               ifelse(haplotypes_F[n,m] == 0 & haplotypes_M[n,m] == 1, 3,
                      ifelse(haplotypes_F[n,m] == 1 & haplotypes_M[n,m] == 0, 4,
                             ifelse(haplotypes_F[n,m] == 1 & haplotypes_M[n,m] == 1, 2, 1)
                      )
               )
        )
    }
  }
  rownames(phased) = sub('\\_1$', '', rownames(haplotypes_F))
  
  phased_collapsed = apply(phased, 1, paste, collapse = "")
  
  return(phased_collapsed)
}

# ---- popSplit ----
# Splits a population into 2
popSplit = function(pop_to_split, nMales_each_pop, nFemales_each_pop){
  
  #' @title Split a population.
  #' 
  #' @description Splits a population into 2, ensuring that each population
  #' derives from different parents.
  #' 
  #' @param pop_to_split is an AlphaSimR population.
  #' @param nMales_each_pop is a scalar of the number of males in each new population.
  #' @param nFemales_each_pop is a scalar of the number of females in each new population. 
  
  cat("\nStarting population split")
  
  cat("\nSampling males")
  randSample_males = selectInd(pop = pop_to_split,
                               nInd = nMales_each_pop*2,
                               use = "rand",
                               sex = "M")
  
  cat("\nSampling females")
  randSample_females = selectInd(pop = pop_to_split,
                                 nInd = nFemales_each_pop*2,
                                 use = "rand",
                                 sex = "F")
  
  cat("\nCreating male populations")
  males_pop1 = randSample_males[1:nMales_each_pop]
  males_pop2 = randSample_males[(nMales_each_pop+1):(nMales_each_pop*2)]
  
  cat("\nCreating female populations")
  females_pop1 = randSample_females[1:nFemales_each_pop]
  females_pop2 = randSample_females[(nFemales_each_pop+1):(nFemales_each_pop*2)]
  
  cat("\nMerging population 1")
  Pop1 = mergePops(list(males_pop1,females_pop1))
  cat("\nMerging population 2")
  Pop2 = mergePops(list(males_pop2,females_pop2))
  
  # Testing
  if (length(intersect(Pop1@id, Pop2@id)) != 0) {
    cat("\nWARNING: There are overlaping ids!")
  }
  
  return(list("Pop1"=Pop1,"Pop2"=Pop2))
}


# ---- Add prefix ----

add_prefix = function(pop, prefix){
  new_IDs = paste0(rep(prefix,nInd(pop)),"_",pop@id)
  return(new_IDs)
}

# ---- makeRecentPop ----
makeRecentPop = function(
    previous_pop,    
    males_each_year,
    females_each_year,
    nCrosses,
    nProgeny = 1,
    years_of_breeding,
    return_breeding_groups = F,
    return_last_generations = 1,
    year=1,
    rec_data_param = list("file_to_write", 
                          "pop_name", 
                          "append_to_existing_file"),
    addPrefix = NULL
){
  #' @title Generate recent generations.
  #' 
  #' @description Make recent generations with overlapping generations and selection.
  #' 
  #' @param previous_pop A population to sample the first parents from.
  #' @param males_each_year A vector of the number o males to use from generations ago,
  #' starting from the current generation.
  #' @param females_each_year A vector of the number o females to use from generations ago,
  #' starting from the current generation.
  #' @param nCrosses Number of crosses to make in each generation.
  #' @param nProgeny Number of progeny per cross. Default is 1.
  #' @param years_of_breeding The number of generations to be bred, assuming that 1 
  #' generation = 1 year. The minimum number must be 2.
  #' @param return_breeding_groups True returns male and female breeding groups. Default is
  #'  False.
  #' @param return_last_generations Number of generations to be returned. Default is 1.
  #' @param year Starting year of the breeding process. Default = 1.
  #' @param rec_data_param List with parameters to write record file. If Default, no file 
  #' will be written.
  #' @param addPrefix Provide a string to be added at the beginning of the population's
  #' IDs. Default is NULL.
  
  # Check if total number of males and females is less or equal to the number of crosses.
  try(if(females_each_year[1] > nCrosses/2) 
    stop("Not enough crosses to produce the number of female candidates required."))
  
  try(if(males_each_year[1] > nCrosses/2) 
    stop("Not enough crosses to produce the number of male candidates required."))
  
  # Check if the number of males and females is decreasing for older generations
  if (length(males_each_year) > 1) {
    try(if(males_each_year[1]<males_each_year[2]) 
      stop("The number of males in older generations is greater than in younger generations."))
  }
  
  if (length(females_each_year) > 1) {
    try(if(females_each_year[1]<females_each_year[2]) 
      stop("The number of females in older generations is greater than in younger generations."))
  }
  
  # Check if years_of_breeding is greater or equal to 2
  try(if (years_of_breeding < 2)
    stop("Years of breeding must be equal or greater than 2."))
  
  cat(paste0("\nWorking on generation number ",1,"/",years_of_breeding))
  
  # Select first generation of males
  males = selectInd(pop = previous_pop, nInd = sum(males_each_year), use = "pheno", sex = "M")
  
  pos_start = 1
  sires = c()
  list_of_sires = c()
  for (n in 1:length(males_each_year)) {
    pos_end=pos_start+males_each_year[n]-1
    
    siresN = males[pos_start:pos_end]
    
    siresN = setMisc(x = siresN,
                     node = "yearOfBirth",
                     value = year-n)
    
    sires = c(siresN, sires)
    list_of_sires = c(list_of_sires, siresN)
    assign(paste0("sires",n), siresN)
    
    pos_start = pos_end+1
  }
  
  # checking total number of sires
  try(if(sires@nInd!=sum(males_each_year)) 
    stop("Different number of sires from the specified."))
  
  # Checking for overlapped sires
  try(if(isTRUE(duplicated(sires@id)))
    stop("Repeated sires selected."))
  
  # Select first generation of females
  females = selectInd(pop = previous_pop, 
                      nInd = sum(females_each_year), 
                      use = "pheno", sex = "F")
  
  pos_start = 1
  dams = c()
  list_of_dams = c()
  for (n in 1:length(females_each_year)) {
    pos_end=pos_start+females_each_year[n]-1
    
    damsN = females[pos_start:pos_end]
    
    damsN = setMisc(x = damsN,
                    node = "yearOfBirth",
                    value = year-n)
    
    dams = c(damsN, dams)
    list_of_dams = c(list_of_dams,damsN)
    assign(paste0("dams",n), damsN)
    
    pos_start = pos_end+1
  }
  
  # checking total number of dams
  try(if(dams@nInd!=sum(females_each_year)) 
    stop("Different number of dams from the specified."))
  
  # Checking for overlapped dams
  try(if(isTRUE(duplicated(dams@id)))
    stop("Repeated dams selected."))
  
  # Cross first generation of recent generations
  year = year + 1
  recent_pop = randCross2(dams, sires, nCrosses = nCrosses, nProgeny = nProgeny)
  recent_pop = setMisc(x = recent_pop,
                       node = "yearOfBirth",
                       value = year)
  if(!is.null(addPrefix)){
    recent_pop@id = add_prefix(recent_pop,addPrefix)
  }
  if(rec_data_param[1] != "file_to_write"){
    rec_data(rec_data_param[[1]], 
             recent_pop, 
             rec_data_param[[2]], 
             year, 
             append = rec_data_param[[3]])
  }
  
  # Simulate generations of breeding
  # Define which generation to return
  return_gen = years_of_breeding - return_last_generations + 1
  return_pop_list = list()
  for (gen in (2):(years_of_breeding)) {
    
    year = year+1
    
    cat(paste0("\nWorking on generation number ",gen,"/",years_of_breeding))
    
    sires = c()
    dams = c()
    
    # Select sires
    for (n in length(list_of_sires):1) {
      if (n > 1) {
        siresN = selectInd(list_of_sires[[n-1]], 
                           nInd = males_each_year[n])
        assign(paste0("sires",n), siresN)
        sires = c(siresN, sires)
        list_of_sires[[n]] = siresN
      } else{
        siresN = selectInd(recent_pop,
                           nInd = males_each_year[n], 
                           sex = "M")
        sires = c(siresN, sires)
        list_of_sires[[n]] = siresN
      }
    }
    
    # Select dams    
    for (n in length(list_of_dams):1) {
      if (n > 1) {
        damsN = selectInd(list_of_dams[[n-1]], 
                          nInd = females_each_year[n])
        assign(paste0("dams",n), damsN)
        dams = c(damsN, dams)
        list_of_dams[[n]] = damsN
      } else{
        damsN = selectInd(recent_pop,
                          nInd = females_each_year[n], 
                          sex = "F")
        dams = c(damsN, dams)
        list_of_dams[[n]] = damsN
      }
    }
    
    # Breed
    recent_pop = randCross2(dams, sires, nCrosses = nCrosses, nProgeny = nProgeny)
    recent_pop = setMisc(x = recent_pop,
                         node = "yearOfBirth",
                         value = year)
    if(!is.null(addPrefix)){
      recent_pop@id = add_prefix(recent_pop,addPrefix)
    }
    
    if(rec_data_param[1] != "file_to_write"){
      rec_data(rec_data_param[[1]], 
               recent_pop, 
               rec_data_param[[2]], 
               year, 
               append = TRUE)
    }
    
    if (gen >= return_gen){
      return_pop_list = c(return_pop_list, assign(as.character(gen), recent_pop), use.names=T)
    }
  }
  cat("\nSimulation of recent generations is complete.")
  
  return(return_pop_list)
}

# ---- crossPops ----
crossPops = function(
    pop1,
    pop2,
    nMales,
    nFemales,
    nCrosses,
    nProgeny = 1
    ){
  
  #' @title Cross populations to form 1/2 descendants.
  #'
  #' @description Cross populations making sure that each parent is being bred
  #' with an individual from the other population.
  #' 
  #' @param pop1 The first population
  #' @param pop2 The second population
  #' @param nMales A vector with the number of males from each population
  #' @param nFemales A vector with the number of females from each population
  #' @param nCrosses A vector with the number of crosses F1xM2 and F2xM1 to be performed
  #' @param nProgeny Number of progeny per cross. Default is 1.
  
  # Select males from pop1
  malesA = selectInd(pop = pop1, nInd = nMales[1], sex = "M", returnPop = FALSE)
  
  # Select females from pop1
  femalesA = selectInd(pop = pop1, nInd = nFemales[1], sex = "F", returnPop = FALSE)
  
  # Select males from pop2
  malesB = selectInd(pop = pop2, nInd = nMales[2], sex = "M", returnPop = FALSE)
  
  # Select females from pop2
  femalesB = selectInd(pop = pop2, nInd = nFemales[2], sex = "F", returnPop = FALSE)
  
  # Make crosses
  cross1 = randCross2(females = pop1, males = pop2, nCrosses = nCrosses[1],
                      femaleParents = femalesA, maleParents = malesB,
                      nProgeny = nProgeny)
  cross2 = randCross2(females = pop2, males = pop1, nCrosses = nCrosses[2],
                      femaleParents = femalesB, maleParents = malesA,
                      nProgeny = nProgeny)
  
  # Merge crosses
  crossed_pop = c(cross1, cross2)
  
  return(crossed_pop)
}

# ---- runBLUPF90 ----
runBLUPF90 = function(param_card, # "renum.txt"
                      ped_file=NULL,
                      min_gen){
  
  # Prepare files for RENUMF90
  setwd("05_BLUPF90/")
  
  cli_alert_info("\nInitiating data preparation for BLUPF90.\n")
  extract_ped = paste0("./extract_ped.sh ",
                       "../01_genotypes/pedigree.txt ",
                       "dat1.txt ",
                       "ped1.txt ",
                       min_gen)
  system(command = extract_ped)
  
  if (!is.null(ped_file)){
    prepare_snp_file = paste0("./prepare_snp_file.sh ",
                              paste0("../01_genotypes/",ped_file," "),
                              "../01_genotypes/new_map.map ",
                              "snp_file.txt")
    system(command = prepare_snp_file)
  } else {
    cli_alert_info("No genotypic file was provided.")
  }
  
  # Run RENUMF90
  cli_alert_info("\nRunning RENUMF90.\n")
  
  system(command = "renumf90", input = param_card)
  
  # Run BLUPF90
  cli_alert_info("\nRunning BLUPF90.\n")
  system(command = "blupf90", input = "renf90.par")
  
  setwd("../")
  cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}

# ---- GBLUP_AlphaSimR ----
GBLUP_AlphaSimR = function(pop){
  # Fit RR-BLUP model for genomic predictions
  ans = RRBLUP(pop, simParam=SP)
  pop = setEBV(pop, ans, simParam=SP)
  
  df = data.frame(ID = pop@id, 
                    sex = pop@sex,
                    pheno = pop@pheno,
                    EBV = pop@ebv, 
                    GV = pop@gv)
  names(df) = c('ID','sex','pheno','EBV','GV')
  return(df)
}

# ---- selectCandidates ----
selectCandidates = function(pop, 
                            file_name, 
                            append=TRUE, 
                            method=1, # 1=EBV; 2=GEBV+F; 3=GEBV+Fg
                            top_ebv,  # c(nMales, nFemales)
                                      # int or fraction
                            max_F = NULL # Max inbreeding (num)
                            ){
  
  # Fit RR-BLUP model for genomic predictions
  BLUP = GBLUP_AlphaSimR(pop)
  
  # Read EBVs
  BLUPF90_EBVs = read.table("05_BLUPF90/solutions.orig", header = T,
                            colClasses = c("integer","integer",
                                           "integer","character",
                                           "numeric")
                            )
  
  # Read Fped
  Fped = read.table("05_BLUPF90/renf90.inb",
                    col.names = c("ID","Fped","delete"),
                    colClasses = c("character","numeric","numeric")
  )
  Fped = Fped[,c(1,2)]
  
  # Read Fg
  Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                          col.names = c("ID","Fg"),
                          colClasses = c("character","numeric")
                          )
  
  # Merge dataframes
  merged_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                      by.x = "ID", by.y = "original_id",
                      how = "left")
  merged_data = merge(merged_data, Fped,
                      by = "ID", how = "left")
  merged_data = merge(merged_data, Fg,
                      by = "ID", how = "left")
  
  # Save metrics to scenario file
  if (isTRUE(append)) {
    write.table(merged_data, file_name, sep = ",", 
                append = T, col.names = F,
                quote = F, row.names = F)
  } else {
    write.table(merged_data, file_name, sep = ",", 
                append = F, col.names = T,
                quote = F, row.names = F)
  }
  
  # Check correlations
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and GV: ", 
                        round(cor(merged_data$EBV, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation BLUPF90 EBV and GV: ", 
                        round(cor(merged_data$solution, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and BLUPF90 EBV: ", 
                        round(cor(merged_data$EBV, merged_data$solution),2)))
  cli_alert_info(paste0("\nCorrelation Fped and Fg: ", 
                        round(cor(merged_data$Fped, merged_data$Fg),2)))
  
  # If a fraction was passed for EBV, calculate the number of
  # animals to select
  # males
  if (top_ebv[1]<1) {
    nMales = nrow(merged_data[merged_data$sex == "M",])
    selectMales = ceiling(nMales * top_ebv[1]) # round up
  } else {selectMales = top_ebv[1]}
  
  # females
  if (top_ebv[2]<1) {
    nFemales = nrow(merged_data[merged_data$sex == "F",])
    selectFemales = ceiling(nFemales * top_ebv[2]) # round up
  } else {selectFemales = top_ebv[2]}
  
  # Method 1
  # Select animals based on BLUPF90 GEBV
  if (method == 1) {
    males = merged_data %>%
      select(ID, sex, solution) %>% 
      filter(sex == "M") %>%
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>%
      select(ID, sex, solution) %>%
      filter(sex == "F") %>%
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
  } else if (method == 2) {
    # Select animals based on BLUPF90 EBV and Fped
    males = merged_data %>% 
      select(ID, sex, solution, Fped) %>%
      filter(sex == "M") %>% 
      filter(Fped <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>% 
      select(ID, sex, solution, Fped) %>%
      filter(sex == "F") %>% 
      filter(Fped <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
  } else if (method == 3) {
    # Select animals based on BLUPF90 EBV and Fg
    males = merged_data %>% 
      select(ID, sex, solution, Fg) %>%
      filter(sex == "M") %>% 
      filter(Fg <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>% 
      select(ID, sex, solution, Fg) %>%
      filter(sex == "F") %>% 
      filter(Fg <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
    } else {cli_alert_danger("\nSelection method not identified!\n")}
  
  male_parents = pop[pop@id %in% as.character(males$ID)]
  female_parents = pop[pop@id %in% as.character(females$ID)]
  
  # Check the number of selected animals
  if (male_parents@nInd != selectMales) {
    cli_alert_warning("The number of selected males is different than the required.")
  }
  
  if (female_parents@nInd != selectFemales){
    cli_alert_warning("The number of selected females is different than the required.")
  }
  
  cli_alert_info(paste0("Mean EBV of the population: ",
                        round(mean(merged_data$solution), 2)))
  cli_alert_info(paste0("Mean EBV of selected males: ",
                        round(mean(males$solution), 2)))
  cli_alert_info(paste0("Mean EBV of selected females: ",
                        round(mean(females$solution), 2)))
  
  parents = c(male_parents, female_parents)
  
  return(parents)
}

# ---- makeFakePed ----
makeFakePed = function(pop) {
  cli_alert_info("\nStart fake pedigree\n")
  males = pop[pop@sex == "M"]
  females = pop[pop@sex == "F"]
  
  fakePed = expand.grid(males@id, females@id,
                        KEEP.OUT.ATTRS=FALSE,
                        stringsAsFactors=FALSE)
  names(fakePed) = c("sire", "dam")
  
  fakePed = fakePed %>% 
    dplyr::mutate(ID = paste0(sire,"_",dam)) %>% 
    select(ID, sire, dam)
  
  cli_alert_success("\nFake pedigree completed.\n")

  return(fakePed)
}

# ---- makeFakeHaplos ----
makeFakeHaplos = function(dirToSave,
                          pop=NULL, 
                          malePop=NULL,
                          femalePop=NULL){
  
  cli_alert_info("\nStart fake haplotypes\n")
  
  if (!is.null(pop)) {
    haplo_males = pullSnpHaplo(pop[pop@sex == "M"])
    haplo_females = pullSnpHaplo(pop[pop@sex == "F"])
  }else{
    haplo_males = pullSnpHaplo(malePop)
    haplo_females = pullSnpHaplo(femalePop)
  }
  
  # Create all combinations of males and females
  combinations <- expand.grid(m = rownames(haplo_males), 
                              f = rownames(haplo_females),
                              stringsAsFactors = FALSE)
  
  # Create progeny IDs
  combinations$ID <- paste0(sub("_.*", "", combinations$m), "_",
                            sub("_.*", "", combinations$f), "_",
                            sub(".*_", "", combinations$m),
                            sub(".*_", "", combinations$f))
  
  
  # Add corresponding rows and convert to matrix (0125 genotype)
  # 0 0 -> 0
  # 1 0 -> 1
  # 1 1 -> 2
  haplo_matrix <- haplo_males[combinations$m, ] + haplo_females[combinations$f, ]
  
  # Create a data frame with ID and the collapsed genotypes
  fakeHaplos = data.frame(ID = combinations$ID,
                          geno = apply(haplo_matrix, 1, 
                                       function(row) paste(row, collapse = "")))
  # Order by ID
  fakeHaplos = fakeHaplos[order(fakeHaplos$ID),]
  
  # Save the genotypes in the appropriate format for BLUPF90
  write.table(fakeHaplos,
              paste0(dirToSave,"fakeHaplos.txt"),
              quote=F, row.names = F, 
              col.names = F, sep = "\t")
  
  if (!is.null(pop)) {
    cli_alert_success(paste0("Fake haplotypes for population ",
                             pop," is completed."))
  }else{
    cli_alert_success(paste0("Fake haplotypes for populations ",
                             malePop, " and ", femalePop,
                             " is completed."))
  }
}

# ---- SelectionWrapper ----
SelectionWrapper = function(InitialPop,
                            nInd, # c(nMales,nFemales)
                            nCrosses,
                            nGenerations,
                            SelectionMethod,
                            year,
                            path_Genotypes,
                            pathNewRecords,
                            pathAllRecords, # path_to/pedigree.txt
                            scenarioName
                            ){
  # Run BLUPF90 for the initial population
  runBLUPF90("recent.ped",
             min_gen=year-10)
  
  # Select candidates from the initial population
  initialParents = selectCandidates(pop=InitialPop, 
                                   file_name=pathNewRecords, 
                                   append=FALSE, 
                                   method=SelectionMethod,
                                   top_ebv=nInd)
  
  # Perform matings
  candidates = randCross(initialParents, nCrosses = nCrosses)
  
  # Save pedigree data
  cli_alert_info(paste0("\nRecording data for gen 1 of ",
                        scenarioName,".\n"))
  rec_data(pathAllRecords, candidates, scenarioName, year, append = TRUE)
  
  # Save genotypes
  writePlink(candidates,
             paste0(path_Genotypes,"sc_",scenarioName,"_gen_1"))
  
  # Continue this process for N generations
  for (i in 2:nGenerations+1) {
    cli_alert_info(paste0("\nInitiating gen ",i,
                          " of scenario ", scenarioName,".\n"))
    year=year+1
    cli_alert(paste0("Year: ",year))
    pedParents = paste0("sc_",scenarioName,"_gen_",i-1)
    
    # Run BLUPF90
    runBLUPF90(paste0(pedParents,".ped"),
               min_gen=year-10)
    
    # Select candidates
    parents = selectCandidates(pop=candidates, 
                               file_name=pathNewRecords, 
                               append=TRUE, 
                               method=SelectionMethod,
                               top_ebv=nInd)
    
    # Perform matings
    candidates = randCross(parents, nCrosses = nCrosses)
    
    # Save pedigree data
    rec_data(pathAllRecords, candidates, scenarioName, year, append = TRUE)
    
    # Save genotypes
    ped_name = paste0("sc_",scenarioName,"_gen_",i)
    writePlink(candidates, paste0(geno_path,ped_name))
  }
  cli_alert_success(paste0("/nSelection process for scenario ",
                           scenarioName, " is complete./n"))
}

# ---- selectCandidates2 ----
selectCandidates2 = function(pop, 
                             ped_name, # runBLUPF90
                             min_gen,  # runBLUPF90
                             file_name, 
                             append=TRUE, 
                             method=1, # 1=EBV+F; 2=GEBV+Fg; 3=GEBV+Froh
                             top_ebv,  # c(nMales, nFemales)
                                       # int or fraction
                             max_F = NULL # Max inbreeding (num)
){
  # Select candidates based of fake progenies (pedigree or genomic)
  
  if (method == 1){
    fake_ped = makeFakePed(pop)
  }
  
  # Fit RR-BLUP model for genomic predictions
  BLUP = GBLUP_AlphaSimR(pop)
  
  # Run BLUPF90
  if (method==1){ 
    # Run without genomic information
    runBLUPF90(ped_name,
               param_card="renum.txt",
               min_gen)
  } else {
    # Run with genomic information
    runBLUPF90(ped_name,
               param_card="renum_genomic.txt",
               min_gen)
  }
  
  # Read EBVs
  BLUPF90_EBVs = read.table("05_BLUPF90/solutions.orig", header = T,
                            colClasses = c("integer","integer",
                                           "integer","character",
                                           "numeric"))
  
  # Read Fped
  Fped = read.table("05_BLUPF90/renf90.inb",
                    col.names = c("ID","Fped","delete"),
                    colClasses = c("character","numeric","numeric"))
  Fped = Fped[,c(1,2)]
  
  if (method == 1) {
    # Run BLUP again to obtain Fg
    runBLUPF90(ped_name,
               param_card="renum_genomic.txt",
               min_gen)
  }
  
  # Read Fg
  Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                  col.names = c("ID","Fg"),
                  colClasses = c("character","numeric")
  )
  
  # Merge dataframes
  merged_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                      by.x = "ID", by.y = "original_id",
                      how = "left")
  merged_data = merge(merged_data, Fped,
                      by = "ID", how = "left")
  merged_data = merge(merged_data, Fg,
                      by = "ID", how = "left")
  
  # Save metrics to scenario file
  if (isTRUE(append)) {
    write.table(merged_data, file_name, sep = ",", 
                append = T, col.names = F,
                quote = F, row.names = F)
  } else {
    write.table(merged_data, file_name, sep = ",", 
                append = F, col.names = T,
                quote = F, row.names = F)
  }
  
  # Check correlations
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and GV: ", 
                        round(cor(merged_data$EBV, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation BLUPF90 EBV and GV: ", 
                        round(cor(merged_data$solution, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and BLUPF90 EBV: ", 
                        round(cor(merged_data$EBV, merged_data$solution),2)))
  cli_alert_info(paste0("\nCorrelation Fped and Fg: ", 
                        round(cor(merged_data$Fped, merged_data$Fg),2)))
  
  # If a fraction was passed for EBV, calculate the number of
  # animals to select
  # males
  if (top_ebv[1]<1) {
    nMales = nrow(merged_data[merged_data$sex == "M",])
    selectMales = ceiling(nMales * top_ebv[1]) # round up
  } else {selectMales = top_ebv[1]}
  
  # females
  if (top_ebv[2]<1) {
    nFemales = nrow(merged_data[merged_data$sex == "F",])
    selectFemales = ceiling(nFemales * top_ebv[2]) # round up
  } else {selectFemales = top_ebv[2]}
  
  # Method 1
  # Select animals based on BLUPF90 GEBV
  if (method == 1) {
    males = merged_data %>%
      select(ID, sex, solution) %>% 
      filter(sex == "M") %>%
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>%
      select(ID, sex, solution) %>%
      filter(sex == "F") %>%
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
  } else if (method == 2) {
    # Select animals based on BLUPF90 EBV and Fped
    males = merged_data %>% 
      select(ID, sex, solution, Fped) %>%
      filter(sex == "M") %>% 
      filter(Fped <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>% 
      select(ID, sex, solution, Fped) %>%
      filter(sex == "F") %>% 
      filter(Fped <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
  } else if (method == 3) {
    # Select animals based on BLUPF90 EBV and Fg
    males = merged_data %>% 
      select(ID, sex, solution, Fg) %>%
      filter(sex == "M") %>% 
      filter(Fg <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>% 
      select(ID, sex, solution, Fg) %>%
      filter(sex == "F") %>% 
      filter(Fg <= max_F) %>% 
      arrange(desc(solution)) %>% 
      slice_head(n=selectFemales)
    
  } else {cli_alert_danger("\nSelection method not identified!\n")}
  
  male_parents = pop[pop@id %in% as.character(males$ID)]
  female_parents = pop[pop@id %in% as.character(females$ID)]
  
  # Check the number of selected animals
  if (male_parents@nInd != selectMales) {
    cli_alert_warning("The number of selected males is different than the required.")
  }
  
  if (female_parents@nInd != selectFemales){
    cli_alert_warning("The number of selected females is different than the required.")
  }
  
  cli_alert_info(paste0("Mean EBV of the population: ",
                        round(mean(merged_data$solution), 2)))
  cli_alert_info(paste0("Mean EBV of selected males: ",
                        round(mean(males$solution), 2)))
  cli_alert_info(paste0("Mean EBV of selected females: ",
                        round(mean(females$solution), 2)))
  
  parents = c(male_parents, female_parents)
  
  return(parents)
}