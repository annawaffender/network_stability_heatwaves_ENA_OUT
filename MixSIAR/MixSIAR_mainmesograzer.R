
#### Stable isotope analysis using the MixSIAR package


rm(list=ls(all=TRUE)) # clear workspace


## load needed packages (eventually install.packages)
library(MixSIAR)
library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(rjags)

## Set paths 
your_path <- "/Users/anna/Anna/Karriere/PhD/Scotti Lab/1_WP" # replace with your path to Github folder 
output_path <- file.path(your_path,"/Github/network_stability_heatwaves/MixSIAR/SIA_dataframes")
save_path <- file.path(your_path,"/Github/network_stability_heatwaves/MixSIAR/output")
setwd(output_path)


##########################################################################################
##											                   ##
## Staple isotope data: main mesograzer (MM) & primary producer (PP)       ##
## reference: Ito etal 2024.               ##
##											                   ##
##########################################################################################

SI_Ito_path <- paste0(output_path, "/", "SI_HW_Ito_etal2024.csv")
SI_Ito <- read.csv(SI_Ito_path)


## change names of columns
colnames(SI_Ito)[colnames(SI_Ito) == "d.15N.14N"] <- "d15N" # change column name
colnames(SI_Ito)[colnames(SI_Ito) == "Area.All"] <- "Area_All_N"
colnames(SI_Ito)[colnames(SI_Ito) == "d.13C.12C"] <- "d13C"
colnames(SI_Ito)[colnames(SI_Ito) == "Area.All.1"] <- "Area_All_C"
colnames(SI_Ito)[colnames(SI_Ito) == "Area.All.2"] <- "Area_All_S"
colnames(SI_Ito)[colnames(SI_Ito) == "d.34S.32S"] <- "d34S"
colnames(SI_Ito)[colnames(SI_Ito) == "corr.d.34S.32S"] <- "corr_d34S"
colnames(SI_Ito)[colnames(SI_Ito) == "taxa"] <- "Species"


## adapt species' names 
SI_Ito$Species <- gsub("Epiphytes", "filamentous_algae", SI_Ito$Species) 
SI_Ito <- SI_Ito %>%mutate(Species = gsub(" ", "_", Species))



##########################################################################################
##											                     ##
## Staple isotope data: POC                  ##
## reference: Mittermayr etal, 2014 (doi.pangaea.de/10.1594/PANGAEA.848529) ##
##											                     ##
##########################################################################################

## revised staple isotope data from Mittermayr etal. 2014 (doi.pangaea.de/10.1594/PANGAEA.848529) 

SI_Mittermayr_path <- paste0(output_path, "/", "SI_Mittermayr2014_corrected_AW.csv")
SI_Mittermayr <- read.csv(SI_Mittermayr_path, head=TRUE, sep=";")

## change names of cleaned df columns 
colnames(SI_Mittermayr)[colnames(SI_Mittermayr) == "δ15N....air._new"] <- "d15N" # change column name
colnames(SI_Mittermayr)[colnames(SI_Mittermayr) == "δ13C....PDB._new"] <- "d13C"
colnames(SI_Mittermayr)[colnames(SI_Mittermayr) == "δ34S....CDT."] <- "d34S"
colnames(SI_Mittermayr)[colnames(SI_Mittermayr) == "Temp...C."] <- "Temp"
colnames(SI_Mittermayr)[colnames(SI_Mittermayr) == "Date.Time"] <- "Date"

## select only the need species, overlapping Mittermayr 2024 SI data and species of species_compartments
SI_Mittermayr_filtered <- SI_Mittermayr %>% filter(Species %in% 
                                                     c("Arenicola marina", "Cerastoderma edule", "Corophium volutator","Eteone longa",
                                                       "Harmothoe imbricata","Macoma balthica","Mya truncata","Mytilus edulis","Nereis pelagica",
                                                       "Acartia sp.","Pseudocalanus sp.","Temora sp.","Oithona sp.","Seston"))

## rename the species to match the species names in the file species_compartments
SI_Mittermayr_rename <- SI_Mittermayr_filtered %>%
  dplyr::mutate(Species = dplyr::recode(
    Species,
      "Arenicola marina"      = "Arenicola_marina",
      "Cerastoderma edule"    = "Cerastoderma_edule",
      "Corophium volutator"   = "Corophium",
      "Eteone longa"          = "Eteone_longa",
      "Harmothoe imbricata"   = "Harmothoe_imbricata",
      "Macoma balthica"       = "Macoma_balthica",
      "Mya truncata"          = "Mya_truncata",
      "Mytilus edulis"        = "Mytilus_edulis",
      "Nereis pelagica"       = "Nereis",
      "Acartia sp."           = "Acartia_sp",
      "Pseudocalanus sp."     = "Pseudocalanus_sp",
      "Temora sp."            = "Temora_sp",
      "Oithona sp."           = "Oithona_sp"))


## chosen dates of SI_Mittermayr samples 
Arenicola_Polychaeta_d <- c("08.03.11") # only those two replicates
Eteone_Polychaeata_o <- c("22.03.11","07.04.11","05.05.11") # only 3  replicates
Cerastoderma_Bivalve <- c("21.07.11","18.08.11","30.08.11")
Corophium_Amphipod_d <- c("21.04.11","23.06.11","30.08.11") # only 4 replicates
Harmothoe_Polychaeata_o <- c("21.04.11") # 1 rep
Nereis_Polychaeata_o <- c("04.08.11")  # or this:  c("2011-07-07")
Macoma_Bivalve <- c("22.03.11","07.04.11","21.04.11","05.05.11") 
Mya_Bivalve <- c("21.04.11","29.09.11")  # 2 rep
Mytilus_Bivalve <- c("07.07.11","21.07.11","04.08.11")
Acartia_Copepoda <-c("07.07.11","21.07.11")
Oithona_Copepoda <- c("19.05.11","09.06.11","23.06.11")
Pseudocalanus_Copoepoda <- c("08.03.11","22.03.11","19.05.11","30.08.11") # 6 rep
Tempora_Copepoda <- c("19.05.11","09.06.11","18.08.11","30.08.11") # 5 rep
Seston_dates <- c("07.07.11", "21.07.11","04.08.11")# 6 rep (with only July only 4 rep)

## chosen species sub-samples 
SI_Mittermayr_dates <- subset(SI_Mittermayr_rename, (Species == "Arenicola_marina" & Date %in% Arenicola_Polychaeta_d ) 
                              | (Species == "Eteone_longa" & Date %in% Eteone_Polychaeata_o)
                              | (Species == "Cerastoderma_edule" & Date %in% Cerastoderma_Bivalve)
                              | (Species == "Corophium" & Date %in% Corophium_Amphipod_d)
                              | (Species == "Harmothoe_imbricata" & Date %in% Harmothoe_Polychaeata_o)
                              | (Species == "Nereis" & Date %in% Nereis_Polychaeata_o)
                              | (Species == "Macoma_balthica" & Date %in% Macoma_Bivalve)
                              | (Species == "Mya_truncata" & Date %in% Mya_Bivalve)
                              | (Species == "Mytilus_edulis" & Date %in% Mytilus_Bivalve)
                              | (Species == "Acartia_sp" & Date %in% Acartia_Copepoda)
                              | (Species == "Oithona_sp" & Date %in% Oithona_Copepoda)
                              | (Species == "Pseudocalanus_sp" & Date %in% Pseudocalanus_Copoepoda)
                              | (Species == "Temora_sp" & Date %in% Tempora_Copepoda)
                              | (Species == "Seston" & Date %in% Seston_dates))


## create one SI_POC_month df, one for each heatwave (triplicate the dfs)
SI_Mittermayr_0HW <- SI_Mittermayr_dates
SI_Mittermayr_0HW$treatment <- "control"
SI_Mittermayr_1HW <- SI_Mittermayr_dates
SI_Mittermayr_1HW$treatment <- "1HW"
SI_Mittermayr_3HW <- SI_Mittermayr_dates
SI_Mittermayr_3HW$treatment <- "3HW"

## Combine all treatment data frames into one
SI_Mittermayr_HWs <- rbind(SI_Mittermayr_0HW, SI_Mittermayr_1HW, SI_Mittermayr_3HW)

## subset from the mesograzer and POC dataframes the columns: Species, d13C, d15N, d34S, treatment
SI_Mittermayr_subset <- SI_Mittermayr_HWs %>% dplyr::select(Species, d13C, d15N, d34S, treatment)
SI_Ito_subset <- SI_Ito %>% dplyr::select(Species, d13C, d15N, d34S, treatment)

## merge all SI dataframes
## nodes = compartments
SI_all_nodes <- bind_rows(SI_Ito_subset,SI_Mittermayr_subset)


################################################
##											  ##
##        Prepare data for MixSIAR           ##
##											  ##
################################################

## add data on species compartment identity 
groups_id_path <- paste0(output_path, "/", "species_compartments.csv")
groups_id <- read.csv(groups_id_path, head=TRUE, sep=";")

## merge the two datasets by species 
SI_all_nodes_id <- plyr::join(SI_all_nodes,groups_id, by = "Species") 

## check, if there are NAs, if yes, clean new df, deleting missing values in Node column  
if (any(is.na(SI_all_nodes_id$Node))) {
  SI_all_nodes_id <- SI_all_nodes_id[-which(is.na(SI_all_nodes_id$Node)), ] 
}


## Create subsets of big dataframe; with only the HW groups; only nodes and SI
HW_treatments <- c("control", "1HW", "3HW")
SI_HW_dataset <- subset(SI_all_nodes_id, SI_all_nodes_id$treatment %in% HW_treatments) # only rows where: treatment = HW_treatments 

SI_HW_13nodes_df <- data.frame(treatment = SI_HW_dataset$treatment,
                               node = SI_HW_dataset$Node,
                               d13C = SI_HW_dataset$d13C,
                               d15N = SI_HW_dataset$d15N,
                               d34S = SI_HW_dataset$d34S)

## Remove the nodes that are not used in the analysis (use only nodes that have at least 3 food sources)
SI_HW_df <- SI_HW_13nodes_df[!(SI_HW_13nodes_df$node %in% c("Amphipod_deposit","Bivalve_filter","Polychaeta_deposit","Polychaeta_omnivore","Copepoda")), ]



#################################################
##											  ##
## Import binary feeding network (MM, PP, POC) ##
##											  ##
#################################################

## MM: Main mesograzer
## PP: Primary producer
## POC: particulate organic carbon (POC)

binary_foodweb_path <- paste0(output_path, "/", "binary_network_MM_PP_POC.csv")
binary_foodweb <- read.csv(binary_foodweb_path, header = T,sep = ";")
rownames(binary_foodweb) <- binary_foodweb$X
binary_foodweb$X <-NULL


## Organizing the data for MixSIAR ###
HW_prey_list <- HW_prey_list_comp <- HW_TEF <- as.list(rep(NA, ncol(binary_foodweb)))

## produce a list of the prey names for each compartment: include in the list the staple isotope values
for (i in 1:ncol(binary_foodweb)) {
  Prey_list <- rownames(binary_foodweb)[binary_foodweb[,i] == 1]
  HW_prey_list[[(i)]] <- Prey_list
}

HW_prey_list <- setNames(HW_prey_list, rownames(binary_foodweb))


## prepare a list containing data frames with the prey data
for (j in 1:ncol(binary_foodweb)) {
  Prey_data <- subset(SI_HW_df, SI_HW_df$node %in% HW_prey_list[[j]])
  HW_prey_list_comp[[(j)]] <- Prey_data #add the list of prey data in the main list
}

HW_prey_list_comp <- setNames(HW_prey_list_comp, rownames(binary_foodweb))



#################################################
##											  ##
##                   Run MixSIAR               ##
##											  ##
#################################################

## save data in save_path
prey_list_file <- paste(save_path, "HW_prey_list_comp_experiment.RData", sep = "/")
save(HW_prey_list_comp, file = prey_list_file) 


## MixSIAR 
options(device = function(...) pdf(file = NULL))

for (m in 4:(ncol(binary_foodweb) - 1)) {
  
  setwd(save_path)
  
  group_id <- colnames(binary_foodweb)[m]
  file3 <- paste0("dataset_experiment_", group_id)
  dir.create(paste(file3, sep = "/"))
  
  # Consumer data
  target_cons <- subset(SI_HW_df, SI_HW_df$node %in% group_id)
  Consumer_SI <- data.frame(
    d13C = target_cons$d13C,
    d15N = target_cons$d15N,
    d34S = target_cons$d34S,
    treatment = target_cons$treatment
  )
  file_consumer <- paste(save_path, file3, "Consumer_SI.csv", sep = "/")
  #write.csv(Consumer_SI, file = file_consumer, row.names = FALSE)
  
  # Source data
  Source_SI <- data.frame(
    source = HW_prey_list_comp[[m]]$node,
    d13C = HW_prey_list_comp[[m]]$d13C,
    d15N = HW_prey_list_comp[[m]]$d15N,
    d34S = HW_prey_list_comp[[m]]$d34S,
    treatment = HW_prey_list_comp[[m]]$treatment
  )
  file_source <- paste(save_path, file3, "Source_SI.csv", sep = "/")
  #write.csv(Source_SI, file = file_source, row.names = FALSE)
  
  # Aggregate sources
  mean_source <- aggregate(. ~ source + treatment, Source_SI, mean)
  colnames(mean_source) <- c("source", "treatment", "Meand13C", "Meand15N", "Meand34S")
  sd_source <- aggregate(. ~ source + treatment, Source_SI, sd)
  colnames(sd_source) <- c("source", "treatment", "SDd13C", "SDd15N", "SDd34S")
  meand_sd_source <- plyr::join(mean_source, sd_source, by = c("source", "treatment"))
  Source_SI$n <- rep(1, nrow(Source_SI))
  n_source <- aggregate(n ~ source + treatment, Source_SI, sum)
  source_agg <- plyr::join(meand_sd_source, n_source, by = c("source", "treatment"))
  file_source_agg <- paste(save_path, file3, "source_agg.csv", sep = "/")
  #write.csv(source_agg, file = file_source_agg, row.names = FALSE)
  
  # fixed trophic enrichment factors (TEF) 
  source_label <- unique(source_agg$source)
  TEF.C <- 0.8 # reference: Zanden and Rasmussen 2001
  TEF.N <- 2.6 # reference: Brauns et al. 2018
  TEF.S <- 0.3 # reference: Mittermayr et al. 2014
  TEF.sd <- 0
  discrmination_data <- data.frame(
    source = source_label,
    Meand13C = rep(TEF.C, length(source_label)),
    SDd13C = rep(TEF.sd, length(source_label)),
    Meand15N = rep(TEF.N, length(source_label)),
    SDd15N = rep(TEF.sd, length(source_label)),
    Meand34S = rep(TEF.S, length(source_label)),
    SDd34S = rep(TEF.sd, length(source_label))
  )
  
  if (nrow(discrmination_data) > 1) {
    file_TEF <- paste(save_path, file3, "discrmination_data.csv", sep = "/")
    #write.csv(discrmination_data, file = file_TEF, row.names = FALSE)
    
    setwd(paste(save_path, file3, sep = "/"))
    
    mix <- load_mix_data(
      filename = file_consumer,
      iso_names = c("d13C", "d15N", "d34S"),
      factors = "treatment",
      fac_random = TRUE,
      fac_nested = FALSE,
      cont_effects = NULL
    )
    
    source <- load_source_data(
      filename = file_source_agg,
      source_factors = "treatment",
      conc_dep = FALSE,
      data_type = "means", 
      mix
    )
    
    discr <- load_discr_data(filename = file_TEF, mix)
    
    # Suppress plot window and generate isospace plot to file only if needed
    pdf(NULL)  # Suppress Quartz
    plot_data(
      filename = "isospace_plot", 
      plot_save_pdf = TRUE,  # Set to TRUE if you want to save the file
      plot_save_png = FALSE, 
      mix, source, discr
    )
    dev.off()
    
    if (mix$n.iso == 2) {
      calc_area(source = source, mix = mix, discr = discr)
    }
    
    model_filename <- "MixSIAR_model.txt"
    resid_err <- TRUE
    process_err <- FALSE
    #write_JAGS_model(model_filename, resid_err, process_err, mix, source)
    
    jags.1 <- run_model(run = "extreme", mix, source, discr, model_filename, alpha.prior = 1) # run length can be adapted
    
    output_options <- list(
      summary_save = TRUE,
      summary_name = "summary_statistics",
      sup_post = TRUE,
      plot_post_save_pdf = FALSE,
      plot_post_save_png = FALSE,
      sup_pairs = TRUE,
      plot_pairs_save_pdf = FALSE,
      plot_pairs_save_png = FALSE,
      sup_xy = TRUE,
      plot_xy_save_pdf = FALSE,
      plot_xy_save_png = FALSE,
      gelman = FALSE,
      heidel = FALSE,
      geweke = FALSE,
      diag_save = FALSE,
      diag_name = "diagnostics",
      indiv_effect = FALSE,
      diag_save_ggmcmc = FALSE
    )
    
    output_JAGS(jags.1, mix, source, output_options)
  }
}


## analyse output: 
## check: 'isospace_plots' for model fit 
## use: 'summary_statistics.txt' for mean staple isotope values

