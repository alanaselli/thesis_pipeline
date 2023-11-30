library(AlphaSimR)

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

adjust_pos = function(SNPmap, chr_len_vector){
  for (x in 1:nrow(SNPmap)) {
    SNPmap$chr_len_vector[x] = chr_len_vector[SNPmap$chr[x]]
  }
  
  SNPmap = SNPmap %>% 
    mutate(pos_2 = round(pos*chr_len_vector)) %>% 
    select(chr,id,pos,pos_2)
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
runBLUPF90 = function(ped_file, min_gen){
  
  # Prepare files for RENUMF90
  setwd("05_BLUPF90/")
  
  cli_alert_info("\nInitiating data preparation for BLUPF90.\n")
  extract_ped = paste0("./extract_ped.sh ",
                       "../01_genotypes/pedigree.txt ",
                       "dat1.txt ",
                       "ped1.txt ",
                       min_gen)
  
  prepare_snp_file = paste0("./prepare_snp_file.sh ",
                            paste0("../01_genotypes/",ped_file," "),
                            "../01_genotypes/new_map.map ",
                            "snp_file.txt")
  
  system(command = extract_ped)
  
  system(command = prepare_snp_file)
  
  # Run RENUMF90
  cli_alert_info("\nRunning RENUMF90.\n")
  
  system(command = "./renumf90", input = "renum.txt")
  
  # Run BLUPF90
  cli_alert_info("\nRunning BLUPF90.\n")
  system(command = "./blupf90", input = "renf90.par")
  
  setwd("../")
  cli_alert_success("\nAnalyses with BLUPF90 completed.\n")
}

# ---- selectCandidates ----
selectCandidates = function(pop, 
                            file_name, 
                            append=TRUE, 
                            method=1, # 1=EBV; 2=EBV+Fg
                            top_ebv,  # c(nMales, nFemales)
                                      # int or fraction
                            Fg_threshold = NULL # Max inbreeding
                            ){
  
  # Fit RR-BLUP model for genomic predictions
  ans = RRBLUP(pop, simParam=SP)
  pop = setEBV(pop, ans, simParam=SP)
  
  BLUP = data.frame(ID = pop@id, 
                    sex = pop@sex,
                    pheno = pop@pheno,
                    EBV = pop@ebv, 
                    GV = pop@gv)
  names(BLUP) = c('ID','sex','pheno','EBV','GV')
  
  # Read EBVs
  BLUPF90_EBVs = read.table("05_BLUPF90/solutions.orig", header = T
                            # colClasses = c("integer","integer",
                            #                "integer","character",
                            #                "numeric")
                            )
  
  # Read Fg
  Fg = read.table("05_BLUPF90/DiagGOrig.txt",
                  col.names = c("ID","Fg")
                  #colClasses = c("character","numeric")
                  )
  
  # Merge dataframes
  merged_data = merge(BLUP, BLUPF90_EBVs[,c("original_id","solution")],
                      by.x = "ID", by.y = "original_id",
                      how = "left")
  merged_data = merge(merged_data, Fg,
                      by = "ID", how = "left")
  
  write.csv(merged_data, file_name, 
            append = append,
            quote = F, row.names = F)
  
  # Check correlations
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and GV: ", 
                        round(cor(merged_data$EBV, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation BLUPF90 EBV and GV: ", 
                        round(cor(merged_data$solution, merged_data$GV),2)))
  cli_alert_info(paste0("\nCorrelation AlphaSimR EBV and BLUPF90 EBV: ", 
                        round(cor(merged_data$EBV, merged_data$solution),2)))
  
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
  
  
  # Select animals based on BLUPF90 EBV
  if (method == 1) {
    males = merged_data %>%
      arrange(desc(solution)) %>% 
      filter(sex == "M") %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>%
      arrange(desc(solution)) %>% 
      filter(sex == "F") %>% 
      slice_head(n=selectFemales)
  }
  
  # Select animals based on BLUPF90 EBV and Fg
  if (method == 2) {
    
    # Check how many males and females pass the Fg filter
    check = merged_data %>% 
      filter(Fg <= Fg_threshold) %>% 
      group_by(sex) %>% 
      summarise(n=n())
    
    cli_alert_info(paste0(check$n[2]," males pass the Fg filter."))
    cli_alert_info(paste0(check$n[1]," females pass the Fg filter."))
    
    males = merged_data %>%
      arrange(desc(solution)) %>% 
      filter(sex == "M" & Fg <= Fg_threshold) %>% 
      slice_head(n=selectMales)
    
    females = merged_data %>%
      arrange(desc(solution)) %>% 
      filter(sex == "F" & Fg <= Fg_threshold) %>% 
      slice_head(n=selectFemales)
  }
  
  # Check the number of selected animals
  male_parents = pop[pop@id %in% as.character(males$ID)]
  try(if(male_parents@nInd != selectMales) 
    stop("The number of selected males is different than the required."))
  
  female_parents = pop[pop@id %in% as.character(females$ID)]
  try(if(female_parents@nInd != selectFemales) 
    stop("The number of selected females is different than the required."))
  
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
  males = pop[pop@sex == "M"]
  females = pop[pop@sex == "F"]
  
  fakePed = expand.grid(males@id, females@id,
                        KEEP.OUT.ATTRS=FALSE,
                        stringsAsFactors=FALSE)
  names(fakePed) = c("sire", "dam")
  
  fakePed = fakePed %>% 
    dplyr::mutate(ID = paste0(sire,"_",dam)) %>% 
    select(ID, sire, dam)

  return(fakePed)
}

# ---- makeFakeHaplos ----
makeFakeHaplos = function(dirToSave,
                          pop=NULL, 
                          malePop=NULL,
                          femalePop=NULL){
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