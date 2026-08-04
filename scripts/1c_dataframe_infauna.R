

#### Create dataframes for mesograzers and infauna: lengths and biomasses 

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
your_path <- "your_path/Github"  # replace with your path to Github folder 
correct_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/dataframes")
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/outputs")
setwd(correct_path)

############################################
##										 ##
##   load data frames and set parameters  ##
##										 ##
############################################

## define the carbon content of the four species and load parameters for LWR
##  Microdeutopus sp.
##  Hydrobia sp.
##  Jaera albifrons
##  Mytilus edulis

tank_area <- 1.53 
tanks <- c("A1", "A2", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2")
treat <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

## file with various sheets for:
## (1) LWR equations (intercept and slope)
## (2) metabolic parameters (P/B and R/B ratios, assimilation efficiency, and intrinsic mortality)
## (3) respiration rates
file_path_source_param <- paste0(correct_path, "/", "1_source_data.xlsx")
source_data_1 <- as.data.frame(read_xlsx(file_path_source_param, sheet = "a_b_estimates"))
source_data_2 <- as.data.frame(read_xlsx(file_path_source_param, sheet = "B_PB_alpha"))
source_data_R <- as.data.frame(read_xlsx(file_path_source_param, sheet = "respiration"))

file_path_infauna <- paste0(correct_path, "/", "free-living_taxa_sizes_Ito_etal2024.xlsx")
sheets <- c("A1", "A2", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2")

infauna <- paste0(correct_path, "/" ,"infauna_taxa_biomass_Ito_etal2024.xlsx")
infauna_df1 <- as.data.frame(read_xlsx(infauna, sheet = "DW_tank"))

## adding biomasses calculated using LWR 
infauna_df2 <- read.csv(file = file.path(correct_path, "df_infauna_tank_biomass.csv"))
colnames(infauna_df2) <- c("tank", "treatment", "Microdeutopus sp.", "Hydrobia ulvae", "Jaera albifrons", "Mytilus edulis")
names2i <- c("Microdeutopus sp.", "Hydrobia ulvae", "Jaera albifrons", "Mytilus edulis")


############################################
##										 ##
##   Mesograzer & Infauna length              ##
##										 ##
############################################

## initialize an empty list to store the data from each sheet
list_of_data <- list()

## loop through each sheet, read data and add tank label
for(sheet in sheets){
  data <- read_excel(file_path_infauna, sheet = sheet)
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

#write.csv(all.taxa_data, file=file.path(output_path,"df_mesograzer_infauna_indiv_size.csv"), row.names =FALSE)




############################################
##										 ##
##   Infauna length (A_d; G_d; I_d; B_f)  ##
##										 ##
############################################

## subtract infauna
main.infauna_data <- all.taxa_data[all.taxa_data$species %in% c("Microdeutopus sp ", "Hydrobia sp ", "Jaera albifrons ", "Mytilus edulis "), ]
main.infauna_data$species <- gsub(" $", "", main.infauna_data$species)	## remove the spaces from the species names 

## replace species names using ifelse()
main.infauna_data$species <- ifelse(main.infauna_data$species == "Microdeutopus sp", "A_d",
                                    ifelse(main.infauna_data$species == "Hydrobia sp", "G_d",
                                           ifelse(main.infauna_data$species == "Jaera albifrons", "I_d",
                                                  ifelse(main.infauna_data$species == "Mytilus edulis", "B_f", NA))))


#write.csv(main.infauna_data, file =file.path(output_path,"df_infauna_indiv_size.csv"), row.names = FALSE)




############################################
##										 ##
##    Biomass: Length-weight regression   ##
##										 ##
############################################

## add columns with estimates of intercept (a) and slope (b) to calculate the dry weight biomass using LWR
## source: Rumohr et al. (1987) and references therein
## equation format: log10(mass) = a + b·log10(length) --> DW(mg)/length(mm)
## except for Microdeutopus sp.: µg/mm and natural logarithm (log_e) is used

main.infauna_data$estimate_a <- ifelse(main.infauna_data$species == "A_d", source_data_1[which(source_data_1$taxa == "A_d"),"a"],
                                       ifelse(main.infauna_data$species == "G_d", source_data_1[which(source_data_1$taxa == "G_d"),"a"],
                                              ifelse(main.infauna_data$species == "I_d", source_data_1[which(source_data_1$taxa == "I_d"),"a"],
                                                     ifelse(main.infauna_data$species == "B_f", source_data_1[which(source_data_1$taxa == "B_f"),"a"], NA))))

main.infauna_data$estimate_b <- ifelse(main.infauna_data$species == "A_d", source_data_1[which(source_data_1$taxa == "A_d"),"b"],
                                       ifelse(main.infauna_data$species == "G_d", source_data_1[which(source_data_1$taxa == "G_d"),"b"],
                                              ifelse(main.infauna_data$species == "I_d", source_data_1[which(source_data_1$taxa == "I_d"),"b"],
                                                     ifelse(main.infauna_data$species == "B_f", source_data_1[which(source_data_1$taxa == "B_f"),"b"], NA))))

## add column with length-dry weight regression calculation: DW_mg (see: Rumohr et al. 1987)
main.infauna_data$DW_mg <- rep(NA,nrow(main.infauna_data))
for(i in 1:nrow(main.infauna_data)){
  if(main.infauna_data$species[i] == "A_d"){
    main.infauna_data$DW_mg[i] <- (1/1000) * exp(1)^(main.infauna_data$estimate_a[i] + main.infauna_data$estimate_b[i] *log(main.infauna_data$length_mm[i]))}
  else{
    main.infauna_data$DW_mg[i] <- 10^(main.infauna_data$estimate_a[i] + main.infauna_data$estimate_b[i] * log10(main.infauna_data$length_mm[i]))
  }
}

## add column with dry weight in µm
main.infauna_data$DW_µg <- main.infauna_data$DW_mg * 1000

## add a column with the full species name
main.infauna_data$species_full <- rep(NA, nrow(main.infauna_data))
for(i in 1:nrow(main.infauna_data)){
  if(main.infauna_data$species[i] == "A_d")main.infauna_data$species_full[i] <- "Microdeutopus sp."
  else{
    if(main.infauna_data$species[i] == "G_d")main.infauna_data$species_full[i] <- "Hydrobia ulvae"
    else{
      if(main.infauna_data$species[i] == "I_d")main.infauna_data$species_full[i] <- "Jaera albifrons"
      else{
        if(main.infauna_data$species[i] == "B_f")main.infauna_data$species_full[i] <- "Mytilus edulis"
        else main.infauna_data$species_full[i] <- NA
      }
    }
  }
}

## add column with biomass in carbon; mean carbon fractions are stored in file "1_source_data.xlsx"; sheet "B_PB_alpha" 
main.infauna_data$BC_µg <- rep(NA,nrow(main.infauna_data))
for(i in 1:nrow(main.infauna_data)){
  if(is.na(main.infauna_data$species_full[i])==FALSE){
    row_select <- which(source_data_2$species == main.infauna_data$species_full[i] &
                          source_data_2$treatment == main.infauna_data$treatment[i])
  }
  C_conv_factor <- source_data_2[row_select,"c_fract"]
  main.infauna_data$BC_µg[i] <- main.infauna_data$DW_µg[i] * C_conv_factor
}

main.infauna_data$B_mg.m2 <- (main.infauna_data$BC_µg/1000)/tank_area 


## create a dataframe from the tibble object and summarize the carbon biomasses of the three groups

infauna_sediment <- data.frame(main.infauna_data)
tank_n <- sheets			## tank names
tank_c <- length(tank_n)	## tank count
infauna_n <- unique(main.infauna_data$species_full)	## main infauna names
infauna_c <- length(infauna_n)						## infauna count

## create a new data frame with tank_n, unique treatment, four empty place holder columns
treatment <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")
infauna_summary_s1 <- data.frame(tank_n, treatment, rep(0,tank_c), rep(0,tank_c), rep(0,tank_c), rep(0,tank_c))
colnames(infauna_summary_s1) <- c("tank", "treatment", infauna_n)

## sum of the carbon biomass for each mesograzer in the 11 tanks
for(i in 1:tank_c){
  for(j in 1:infauna_c){
    infauna_summary_s1[i, (j+2)] <- sum(infauna_sediment[which(as.character(infauna_sediment$tank) ==
                                                                 tank_n[i] & as.character(infauna_sediment$species_full) == infauna_n[j]),"B_mg.m2"])
  }
}


###################################################################
##										 ##
##   create infauna dataframe: B, alpha, PB,RB, P, R, G, C, E  ##
##										 ##
##################################################################

infauna_df3 <- infauna_df1

row_vect <- which(is.na(match(infauna_df3$species,names2i))==FALSE)
for(i in 1:length(row_vect)){
  for(j in 2:ncol(infauna_df3)){
    infauna_df3[row_vect[i],j] <- infauna_df3[row_vect[i],j] + infauna_df2[(j-1),infauna_df3$species[row_vect[i]]]
  }
}

infauna_df <- infauna_df3[-which(infauna_df1$species=="Littorina littorea"),] # exclude Littorina littorea (not part of infauna)

name_species <- infauna_df[,1]
n_species <- length(name_species)
taxa_species <- c("A_d", "G_d", "I_d", "B_f", "P_d", "A_d", "B_f", "B_f", "P_o", "P_d", "P_d")

df_long_infauna <- data.frame(
  species = rep(infauna_df[,1], each = ncol(infauna_df[,-1])),	## repeat species names
  taxa = rep(taxa_species, each = ncol(infauna_df[,-1])),		## repeat compartment names
  tank = rep(tanks, n_species),									## tank names
  treat = rep(treat, n_species),								## treatments
  DW = as.vector(t(infauna_df[,-1]))							## flatten row-wise
)


## carbon content from empirical staple isotope data (see Ito etal. 2024)
C_SIA <- paste0(correct_path, "/" ,"SIA_HW2015_Ito_etal2024.xlsx")
C_SIA_df <- as.data.frame(read_xlsx(C_SIA, sheet = "consumers_clean"))
C_taxa <- unique(C_SIA_df$taxa)
C_trea <- unique(C_SIA_df$treatment)
CM <- matrix(rep(NA, length(C_taxa)*length(C_trea)), nrow = length(C_taxa)) # create empty matrix for carbon content of all species and treatments
rownames(CM) <- C_taxa
colnames(CM) <- c("0HW", "1HW", "3HW", "field")


## fill empty CM matrix with carbon values 
for(i in 1:length(C_taxa)){
  for(j in 1:length(C_trea)){
    mtc <- which(C_SIA_df$taxa == rownames(CM)[i] & C_SIA_df$trea == colnames(CM)[j])
    if(length(mtc)!=0){
      CM[i,j] <- mean(C_SIA_df[mtc,"C_fraction"])
    }
  }
}
round(CM,4)


## create vectors with appropriate carbon content, R/B and P/B ratios, and assimilation efficiency (alpha)
all_ele <- nrow(df_long_infauna)
CC_v <- PB_v <- RB_v <- alpha_v <- rep(NA,all_ele)

for(i in 1:all_ele){
  sp_sel <- as.character(df_long_infauna$species[i]) # species selected (i)
  tr_sel <- as.character(df_long_infauna$treat[i]) # treatment selected (i)
  ta_sel <- as.character(df_long_infauna$tank[i]) # tank selected
  
  f2r <- which(source_data_2$species == sp_sel & source_data_2$treatment == tr_sel) # find rows of source_data_2 that match conditions 
  frr <- which(source_data_R$species == sp_sel & source_data_R$tank == ta_sel) # find rows of source_data_r that match conditions 
  
  if(length(f2r)==1){
    CC_v[i] <- source_data_2[f2r,"c_fract"]
    PB_v[i] <- source_data_2[f2r,"PB"]
    alpha_v[i] <- source_data_2[f2r,"alpha"]
  }
  
  if(length(frr)==1){
    RB_v[i] <- source_data_R[frr,"RB"]
  }
}


## completing the database to be consistent with the one of mesograzers_full_dataset
df_long_infauna_s1 <- data.frame(df_long_infauna, CC_v, PB_v, RB_v, alpha_v)

## calculation of the biomass, expressed as mgC m-2 d-1 from the product of DW * C/DW
df_long_infauna_s1$B <- df_long_infauna_s1$DW * df_long_infauna_s1$CC_v

df_long_infauna_s1$alpha <- df_long_infauna_s1$alpha_v
df_long_infauna_s1$PB <- df_long_infauna_s1$PB_v
df_long_infauna_s1$RB <- df_long_infauna_s1$RB_v

## calculation of the production, starting from P/B ratios and taxa biomass
df_long_infauna_s1$P <- df_long_infauna_s1$B * df_long_infauna_s1$PB

## calculation of the respiration, starting from R/B ratios and taxa biomass
df_long_infauna_s1$R <- df_long_infauna_s1$B * df_long_infauna_s1$RB

## calculation of the gross production: G = P + R
df_long_infauna_s1$G <- df_long_infauna_s1$P + df_long_infauna_s1$R

## calculation of the consumption, starting from assimilation efficiency and gross production
df_long_infauna_s1$C <- df_long_infauna_s1$G/df_long_infauna_s1$alpha

## calculation of the egestion as a difference between consumption and gross production: E = C - G
df_long_infauna_s1$E <- df_long_infauna_s1$C - df_long_infauna_s1$G
df_long_infauna_final <- df_long_infauna_s1[,-c(5:9)]

#write.csv(df_long_infauna_final, file =file.path(output_path,"df_infauna_tank.csv"), row.names = FALSE)


