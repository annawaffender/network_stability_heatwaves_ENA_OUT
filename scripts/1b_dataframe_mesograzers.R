

#### Create dataframe of all mesograzers: (1) I_o, (2) A_o, (3) G_h 

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.package)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(multcomp)


## Set paths 
your_path <- "your_path/Github" # replace with your path to Github folder  
correct_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/dataframes")
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/outputs")


############################################
##										 ##
##   load data frames and set parameters  ##
##										 ##
############################################

tank <- c("A1", "A2", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2")
treat <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

tank_area <- 1.53 

## file with various sheets for:
## (1) LWR equations (intercept and slope)
## (2) metabolic parameters (P/B and R/B ratios, assimilation efficiency, and intrinsic mortality)
## (3) respiration rates
file_path_source_param <- paste0(correct_path, "/", "1_source_data.xlsx")
source_data_1 <- as.data.frame(read_xlsx(file_path_source_param, sheet = "a_b_estimates"))
source_data_2 <- as.data.frame(read_xlsx(file_path_source_param, sheet = "B_PB_alpha"))
source_data_R <- as.data.frame(read_xlsx(file_path_source_param, sheet = "respiration"))


## a and b estimates (a, b) for length-dry weight (LWR) regressions
##
## (1) Idotea balthica
IS_a <- source_data_1[source_data_1$taxa == "I_o","a"]
IS_b <- source_data_1[source_data_1$taxa == "I_o","b"]

## (2) Gammarus sp.
GS_a <- source_data_1[source_data_1$taxa == "A_o","a"]
GS_b <- source_data_1[source_data_1$taxa == "A_o","b"]

## (3) Littorina littorea
LL_a <- source_data_1[source_data_1$taxa == "G_h","a"]
LL_b <- source_data_1[source_data_1$taxa == "G_h","b"]


## carbon fractions in DW (c_fract)
##
## (1) Idotea balthica
IS_B_0HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "0HW","c_fract"]
IS_B_1HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "1HW","c_fract"]
IS_B_3HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "3HW","c_fract"]

## (2) Gammarus sp.
GS_B_0HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "0HW","c_fract"]
GS_B_1HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "1HW","c_fract"]
GS_B_3HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "3HW","c_fract"]

## (3) Littorina littorea
LL_B_0HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "0HW","c_fract"]
LL_B_1HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "1HW","c_fract"]
LL_B_3HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "3HW","c_fract"]


## assimilation efficiency (alpha) ~ HWs
##
## (1) Idotea balthica
IS_alpha_0HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "0HW","alpha"]
IS_alpha_1HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "1HW","alpha"]
IS_alpha_3HW <- source_data_2[source_data_2$taxa == "I_o" & source_data_2$treatment == "3HW","alpha"]

## (2) Gammarus sp.
GS_alpha_0HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "0HW","alpha"]
GS_alpha_1HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "1HW","alpha"]
GS_alpha_3HW <- source_data_2[source_data_2$taxa == "A_o" & source_data_2$treatment == "3HW","alpha"]
##
## (3) Littorina littorea
LL_alpha_0HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "0HW","alpha"]
LL_alpha_1HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "1HW","alpha"]
LL_alpha_3HW <- source_data_2[source_data_2$taxa == "G_h" & source_data_2$treatment == "3HW","alpha"]


############################################
##										 ##
##   Littorina littorea  biomass ##
##										 ##
############################################

## import data on DW biomass of Littorina littorea estimated for the whole tank from what sieved in 5 beakers
##
## control = A1, A2, D1, D2
## 3HW = B1 (missing data), B2, E1, E2
## 1HW = C1, C2, F1, F2

## species: sediment community
## response variable: biomass
## method: five beakers were sieved with 1000 µm and fixed in Formalin on 17.08.2015
## unit: mg of dry weight (DW)
## beaker volume: 1.4 L
## total sediment volume: 67 L
## tank area: 1.53 m2

## calculated the SFDW (shell-free dry weight) biomass of Littorina littorea using the LWR (length-weight regression)
## in each tank, summed it with the SFDW biomass estimated for the whole tank from the beakers

### SFDW vector of Littorina in the whole tank, estimated from the content in the 5 beakers
Littorina_infauna <- paste0(correct_path, "/" ,"infauna_taxa_biomass_Ito_etal2024.xlsx")
Littorina_infauna_s0 <- as.data.frame(read_xlsx(Littorina_infauna, sheet = "DW_tank"))
Littorina_infauna_v <- unname(unlist(Littorina_infauna_s0[which(as.character(Littorina_infauna_s0$species) == "Littorina littorea"),-1])) 

## convert SFDW into carbon content
Littorina_inf_df <- data.frame(tank,treat,c(Littorina_infauna_v))
colnames(Littorina_inf_df) <- c("tank", "treat", "DW")
carbon <- rep(NA,nrow(Littorina_inf_df))

for(i in 1:nrow(Littorina_inf_df)){
	tr_ll <- Littorina_inf_df[i,"treat"]
	nm_ll <- "Littorina littorea"
	rw_ll <- which(source_data_2$species == nm_ll & source_data_2$treatment == tr_ll)
	conv_f <- source_data_2[rw_ll,"c_fract"]
	## print(conv_f)
	carbon[i] <- Littorina_inf_df$DW[i] * conv_f
}

Littorina_inf_df$carbon <- carbon

#####################################################
##											  ##
## Carbon flow calculations (see: Ito et al., 2024) ##
##											  ##
#####################################################

## path to Excel file
file_path_meso <- paste0(correct_path, "/", "free-living_taxa_sizes_Ito_etal2024.xlsx")
sheets <- c("A1", "A2", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2")

## initialize an empty list to store the data from each sheet
list_of_data <- list()

## loop through each sheet, read data and add tank label
for(sheet in sheets){
  data <- read_excel(file_path_meso, sheet = sheet)
  ## add a column for tank label
  data$tank <- sheet
  ## add the dataframe to the list
  list_of_data[[sheet]] <- data
}

## combine all sheets into one dataframe
combined_data <- bind_rows(list_of_data)

## reshape the data frame from "wide" to "long" to create one column for all length measurements 
all.taxa_data <- combined_data %>%
  pivot_longer(
    cols = -tank, 			## exclude the 'tank' column from pivoting
    names_to = "species",	## name the new species column
    values_to = "length_µm"	## name the new length column
  )

## change species names by removing "(length [µm])"
all.taxa_data <- all.taxa_data %>%
  mutate(species = str_replace(species, "\\(length \\[µm\\]\\)", ""))

## remove all NAs 
all.taxa_data <- all.taxa_data %>%
  drop_na()

## add column length_mm
all.taxa_data$length_mm <- all.taxa_data$length_µm/1000

## add a new 'treatment' column based on the 'tank' column 
all.taxa_data$treatment <- ifelse(all.taxa_data$tank %in% c("A1", "A2", "D1", "D2"), "0HW",
                                  ifelse(all.taxa_data$tank %in% c("C1", "C2", "F1", "F2"), "1HW",
                                         ifelse(all.taxa_data$tank %in% c("B2", "E1", "E2"), "3HW", NA)))

## add a new 'temperature' column based on the 'treatment' column 
all.taxa_data$temperature <- ifelse(all.taxa_data$treatment %in% c("0HW"), "20",
                                    ifelse(all.taxa_data$treatment %in% c("1HW", "3HW"), "25", NA))

## subtract main mesograzers
main.mesograzer_data <- all.taxa_data[all.taxa_data$species %in% c("Gammarus sp ", "Idotea sp ", "Littorina littorea "), ]
main.mesograzer_data$species <- gsub(" $", "", main.mesograzer_data$species)	## remove the spaces from the species names 

## replace species names using ifelse()
main.mesograzer_data$species <- ifelse(main.mesograzer_data$species == "Gammarus sp", "A_o",
                                       ifelse(main.mesograzer_data$species == "Idotea sp", "I_o",
                                              ifelse(main.mesograzer_data$species == "Littorina littorea", "G_h", 
                                                     main.mesograzer_data$species)))


## add columns with estimates of intercept (a) and slope (b) to calculate the dry weight biomass using LWR
## source: Supp. SI.1 in Ito et al. (2024); equation format: log10(mass) = a + b·log10(length) --> DW(mg)/length(mm)

main.mesograzer_data$estimate_a <- ifelse(main.mesograzer_data$species == "I_o", IS_a,
                                          ifelse(main.mesograzer_data$species == "A_o", GS_a,
                                                 ifelse(main.mesograzer_data$species == "G_h", LL_a, NA)))

main.mesograzer_data$estimate_b <- ifelse(main.mesograzer_data$species == "I_o", IS_b,
                                          ifelse(main.mesograzer_data$species == "A_o", GS_b,
                                                 ifelse(main.mesograzer_data$species == "G_h", LL_b, NA)))

## add column with length-dry weight regression calculation: DW_mg (see: Ito et al. 2024, Supp. SI.1)
main.mesograzer_data$DW_mg <- 10^(main.mesograzer_data$estimate_a + main.mesograzer_data$estimate_b * log10(main.mesograzer_data$length_mm))

## add column with dry weight in µm
main.mesograzer_data$DW_µg <- main.mesograzer_data$DW_mg * 1000

## add column with biomass in carbon
## mean carbon fractions (species average values per treatment at the end of experiment; see the file: SI_HW 2015):
## carbon content values match those used to normalize respiration rates and calculate R/B from incubations

## conversion factor for Littorina littorea (DW --> SFDW; see Rumohr 1987, p.18 Hydrobia sp.)
cf_LL_DW_SFDW <- 0.156

main.mesograzer_data$B_µg <- ifelse(main.mesograzer_data$species == "I_o" & main.mesograzer_data$treatment == "0HW",
                                     main.mesograzer_data$DW_µg * IS_B_0HW, 
                                     ifelse(main.mesograzer_data$species == "I_o" & main.mesograzer_data$treatment == "1HW",
                                            main.mesograzer_data$DW_µg * IS_B_1HW,
                                            ifelse(main.mesograzer_data$species == "I_o" & main.mesograzer_data$treatment == "3HW",
                                                   main.mesograzer_data$DW_µg * IS_B_3HW, 
                                                   ##
                                                   ifelse(main.mesograzer_data$species == "A_o" & main.mesograzer_data$treatment == "0HW",
                                                          main.mesograzer_data$DW_µg * GS_B_0HW, 
                                                          ifelse(main.mesograzer_data$species == "A_o" & main.mesograzer_data$treatment == "1HW",
                                                                 main.mesograzer_data$DW_µg * GS_B_1HW,
                                                                 ifelse(main.mesograzer_data$species == "A_o" & main.mesograzer_data$treatment == "3HW",
                                                                        main.mesograzer_data$DW_µg * GS_B_3HW,  
                                                                        ##
                                                                        ifelse(main.mesograzer_data$species == "G_h" & main.mesograzer_data$treatment == "0HW",
                                                                               main.mesograzer_data$DW_µg * LL_B_0HW * cf_LL_DW_SFDW,				 
                                                                               ifelse(main.mesograzer_data$species == "G_h" & main.mesograzer_data$treatment == "1HW",
                                                                                      main.mesograzer_data$DW_µg * LL_B_1HW * cf_LL_DW_SFDW,				
                                                                                      ifelse(main.mesograzer_data$species == "G_h" & main.mesograzer_data$treatment == "3HW",
                                                                                             main.mesograzer_data$DW_µg * LL_B_3HW * cf_LL_DW_SFDW, NA)))))))))	

main.mesograzer_data$B_mg.m2 <- (main.mesograzer_data$B_µg/1000)/tank_area 
mesograzers <- data.frame(main.mesograzer_data)

#write.csv(mesograzers, file =file.path(output_path,"df_mesograzers_individ.csv"), row.names = FALSE)


##########################################################
##													  ##
##  database carbon biomasses per tank  ##
##													  ##
#########################################################

## create a dataframe from the tibble object and summarize the carbon biomasses of the three groups
tank_n <- unique(mesograzers$tank)	## tank names
tank_c <- length(tank_n)			## tank count

meso_n <- unique(mesograzers$species)	## main mesograzer names
meso_c <- length(meso_n)				## mesograzer count

## create a new data frame with tank_n as first column, unique treatment (each treatment only once) and three empty place holder columns
mesograzers_summary_s1 <- data.frame(tank_n, unlist(lapply(tank_n,function(x)unique(mesograzers[
  which(mesograzers$tank==x),"treatment"]))), rep(0,tank_c), rep(0,tank_c), rep(0,tank_c))

colnames(mesograzers_summary_s1) <- c("tank", "treatment", meso_n)

## sum of the carbon biomass for each mesograzer in the 11 tanks
for(i in 1:tank_c){
  for(j in 1:meso_c){
    mesograzers_summary_s1[i, (j+2)] <- sum(mesograzers[which(as.character(mesograzers$tank) ==
                                                             tank_n[i] & as.character(mesograzers$species) == meso_n[j]),"B_mg.m2"])
  }
}



## adding the carbon biomass of Littorina estimated from the 5 beakers to the one calculated from LWR
mesograzers_summary <- mesograzers_summary_s1

for(i in 1:nrow(mesograzers_summary_s1)){
	row_ll_sel <- which(as.character(Littorina_inf_df$tank) == as.character(mesograzers_summary$tank[i]))
	mesograzers_summary[i,"G_h"] <- mesograzers_summary_s1[i,"G_h"] + Littorina_inf_df[row_ll_sel,"carbon"]
}



##########################################################
##													  ##
## final dataframe: B, alpha, PB, RB, P, R, E,C,G  ##
##													  ##
##########################################################

n_tanks <- length(sheets)
treat <- as.character(mesograzers_summary$treatment)

mesograzers_full_dataset <- data.frame(rbind(
  cbind(rep("Gammarus sp.",n_tanks), rep("A_o",n_tanks), sheets, treat),
  cbind(rep("Idotea balthica",n_tanks), rep("I_o",n_tanks), sheets, treat),
  cbind(rep("Littorina littorea",n_tanks), rep("G_h",n_tanks), sheets, treat)))

mesograzers_B <- c(as.numeric(mesograzers_summary$A_o),
                   as.numeric(mesograzers_summary$I_o),
                   as.numeric(mesograzers_summary$G_h))

mesograzers_full_dataset$B <- mesograzers_B

colnames(mesograzers_full_dataset) <- c("species", "taxa", "tank", "treat", "B")

n_row <- nrow(mesograzers_full_dataset)

alpha_v <- PB_v <- RB_v <- rep(NA,n_row)	## create vectors of length n_row
count <- 1

for(i in 1:n_row){
  taxa_id <- mesograzers_full_dataset[i,"taxa"]
  treat_id <- mesograzers_full_dataset[i,"treat"]
  tank_id <- mesograzers_full_dataset[i,"tank"]
  ##
  meta_row <- which(source_data_2$taxa == taxa_id & source_data_2$treatment == treat_id)
  resp_row <- which(source_data_R$taxa == taxa_id & source_data_R$tank == tank_id)
  ##
  alpha_v[count] <- source_data_2[meta_row,"alpha"]
  PB_v[count] <- source_data_2[meta_row,"PB"]
  RB_v[count] <- mean(source_data_R[resp_row,"RB"])
  count <- count + 1
}

mesograzers_full_dataset$alpha <- alpha_v
mesograzers_full_dataset$PB <- PB_v
mesograzers_full_dataset$RB <- RB_v

## calculation of the production, starting from P/B ratios and taxa biomass
mesograzers_full_dataset$P <- mesograzers_full_dataset$B * mesograzers_full_dataset$PB

## calculation of the respiration, starting from R/B ratios and taxa biomass
mesograzers_full_dataset$R <- mesograzers_full_dataset$B * mesograzers_full_dataset$RB

## calculation of the gross production: G = P + R
mesograzers_full_dataset$G <- mesograzers_full_dataset$P + mesograzers_full_dataset$R

## calculation of the consumption, starting from assimilation efficiency and gross production
mesograzers_full_dataset$C <- mesograzers_full_dataset$G/mesograzers_full_dataset$alpha

## calculation of the egestion as a difference between consumption and gross production: E = C - G
mesograzers_full_dataset$E <- mesograzers_full_dataset$C - mesograzers_full_dataset$G

## save the database on biomass, metabolic rates, and full quantification of the flows
mesograzers_full_dataset_names <- mesograzers_full_dataset
mesograzers_full_dataset_names$taxa <- c(rep("A_o",n_tanks), rep("I_o",n_tanks), rep("G_h",n_tanks))

#write.csv(mesograzers_full_dataset_names, file =file.path(output_path,"df_mesograzers_tank.csv"), row.names = FALSE)



