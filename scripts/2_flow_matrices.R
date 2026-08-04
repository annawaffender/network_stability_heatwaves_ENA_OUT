

#### Create flow matrices

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.packages)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

## Set paths 
your_path <- "your_path/Github"  # replace with your path to Github folder 
source_file_path <-  file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/dataframes")
flow_matrix_save <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/outputs/flow_matrices")
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/outputs")
setwd(output_path)


############################################
##										 ##
##   load data frames and set parameters  ##
##										 ##
############################################

compartment_names <- c("ZM", "FV", "FA", "A_d", "A_o", "I_d", "I_o", "G_d", "G_h", "B_f", "P_d", "P_o", "POC")
tank <- c("A1", "A2", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2")
treat <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

macrophytes_nonliving_df <- data.frame(read.csv(file = "df_macrophytes_POC.csv"))
mesograzers_df <- data.frame(read.csv(file = "df_mesograzers_tank.csv"))
infauna_df <- data.frame(read.csv(file = "df_infauna_tank.csv"))
consumers_df <- rbind(mesograzers_df, infauna_df)

macrophytes_exports <- read.csv(file.path(source_file_path, "macrophytes_exports.csv")) # export file for the three primary producers (mortalities)




###################################
##						##
##   compartment attributes      ##
##						##
###################################

## build database where, for each compartment, the following properties are listed for the different tanks:
## biomass (B)
## consumption (C)
## production (P)
## respiration (R)
## egestion (E)
## assimilation efficiency (alpha)
## mortality (M)

## list the names of all compartments
all_compartments_old_order <- c(unique(as.character(macrophytes_nonliving_df$taxa)),
unique(as.character(mesograzers_df$taxa)),
unique(as.character(infauna_df$taxa)))

all_compartments <- c(all_compartments_old_order[1:3],all_compartments_old_order[5:13],all_compartments_old_order[4]) # set POC compartment as the last of the food web
n_compartments <- length(all_compartments) # total number of compartments in the food web

## prepare the df for all compartments 
compartments_df <- data.frame(unlist(lapply(all_compartments,function(x)rep(x,length(tank)))),
rep(tank,n_compartments),rep(treat,n_compartments))
colnames(compartments_df) <- c("compartments", "tank", "treat")
compartments_df$M <- compartments_df$alpha <- compartments_df$E <- compartments_df$R <- compartments_df$P <- compartments_df$C <- compartments_df$B <- rep(NA,nrow(compartments_df))

select_PP_POC<- unique(as.character(macrophytes_nonliving_df$taxa))


# loop through each row of compartments_df, update columns based on whether compartments is PP (primary producer) or POC (particulate organic carbon)
for(i in 1:nrow(compartments_df)){
	{
	  # checks if current compartment is in select_PP_POC vector
	if(is.na(match(as.character(compartments_df$compartments[i]),select_PP_POC))==FALSE){
		row_selected <- which(is.na(match(as.character(macrophytes_nonliving_df$taxa),
		as.character(compartments_df$compartments[i])))==FALSE & is.na(match(as.character(macrophytes_nonliving_df$tank),
		as.character(compartments_df$tank[i])))==FALSE)
		# if in select_PP_POC, calculate and assign these values: 
		compartments_df$B[i] <- sum(macrophytes_nonliving_df[row_selected,"B"])
		compartments_df$C[i] <- sum(macrophytes_nonliving_df[row_selected,"GPP"])
		compartments_df$P[i] <- sum(macrophytes_nonliving_df[row_selected,"NPP"])
		compartments_df$R[i] <- sum(macrophytes_nonliving_df[row_selected,"R"])
		## compartments_df$E[i] <- sum(macrophytes_nonliving_df[row_selected,"E"])
	  }
	  # if not PP or POC, find matching row in consumers_df, where taxa and tank match 
	else{
		row_selected <- which(is.na(match(as.character(consumers_df$taxa),
		as.character(compartments_df$compartments[i])))==FALSE & is.na(match(as.character(consumers_df$tank),
		as.character(compartments_df$tank[i])))==FALSE)
		##
		# calculate and assign these values:  
		compartments_df$B[i] <- sum(consumers_df[row_selected,"B"])
		compartments_df$C[i] <- sum(consumers_df[row_selected,"C"])
		compartments_df$P[i] <- sum(consumers_df[row_selected,"P"])
		compartments_df$R[i] <- sum(consumers_df[row_selected,"R"])
		compartments_df$E[i] <- sum(consumers_df[row_selected,"E"])
		##
		# for non PP or POC, calculate alpha
		alpha_v <- rep(0,length(row_selected))          # alpha_v as zero vector with length of row_selected
		alpha_sel <- consumers_df[row_selected,"alpha"] # select alpha from consumer_df
		B_sel <- consumers_df[row_selected,"B"]         # select B from consumer_df
		
		# fill alpha_v with weighted alpha values (alpha_sel *B_selected/B_total compartment)
		for(k in 1:length(row_selected))if(compartments_df$B[i]!=0)alpha_v[k] <- alpha_sel[k] * (B_sel[k]/compartments_df$B[i])
		compartments_df$alpha[i] <- sum(alpha_v)
		}
	}
}

#write.csv(compartments_df, file = "df_network_compartments.csv", row.names = FALSE)



###################################
##					 ##
##   flow matrices   ##
##					 ##
###################################

## build a binary matrix from the edgelist (including feeding preferences)
file_path_source_EL <- paste0(source_file_path, "/", "1_source_data.xlsx")
source_data_EL <- as.data.frame(read_xlsx(file_path_source_EL, sheet = "edgelist"))

binary <- matrix(rep(0,length(compartment_names)^2), nrow = length(compartment_names))
rownames(binary) <- colnames(binary) <- compartment_names

## fill binary matrix: for extracted resource and consumer values from row i --> assign 1  
for(i in 1:nrow(source_data_EL)){
	binary[source_data_EL[i,"resource"],source_data_EL[i,"consumer"]] <- 1
}

## prepare a list where the flow matrices of each tank will be stored 
flows_list <- as.list(rep(NA,length(tank)))
names(flows_list) <- tank

## import feeding preferences of mesograzers
file_path_source_SIA <- paste0(source_file_path, "/", "1_source_data.xlsx")
source_data_SIA <- as.data.frame(read_xlsx(file_path_source_SIA, sheet = "SIA"))

## start adding the carbon flows to mesograzers (MG) as their feeding preferences were determined with staple isotope analysis (SIA)
## consumption (C)* feed_pref (SIA data)
MG <- c("A_o", "I_o", "G_h") 
MG_n <- length(MG)
for(k in 1:length(tank)){
	flows_list[[k]] <- binary # assign binary matrix to each tank
	tk_ID <- tank[k]
	##
	for(i in 1:n_compartments){ # iterate over each compartments
		for(j in 1:MG_n){  # iterate over each species
			selected_col <- which(compartment_names == MG[j]) # selected col, which corresponds to MG[j]
			
			if(flows_list[[k]][i,selected_col]!=0){ # check if feeding relationship exists 
				total_consump <- compartments_df[which(as.character(compartments_df$compartments) == MG[j] & as.character(compartments_df$tank) == tk_ID),"C"] # total consumption (C) for predator MG[j] in tank
				feed_pref <- source_data_SIA[which(source_data_SIA$resource == compartment_names[i] & source_data_SIA$predator == MG[j] & # feeding pref. of MG[j] consuming compartment_name[i] under given treat,(from 1_source_data_PoPo, SIA) 
				source_data_SIA$treat == treat[k]),"SIA"]
				flows_list[[k]][i,selected_col] <- total_consump * feed_pref # replace binary value (1), with weighted interaction strength 
			}
		}
	}
}


## add the carbon flows to the other consumers (OC); feeding preferences determined according to donor-controlled criteria
## consumption (C)* feed_pref (compartment biomass/total biomass)
OC <- c("A_d", "I_d", "G_d", "B_f", "P_d", "P_o")
OC_n <- length(OC)
flows_list_s1 <- flows_list

for(k in 1:length(tank)){
	tk_ID <- tank[k]
	##
	for(i in 1:n_compartments){ # loop through all compartments
		for(j in 1:OC_n){  # loop through "other consumer"
			selected_col <- which(compartment_names == OC[j])
			all_res <- compartment_names[which(flows_list[[k]][,selected_col]!=0)]
			tot_bio <- sum(compartments_df[which(is.na(match(as.character(compartments_df$compartments),all_res)) == FALSE &
			as.character(compartments_df$tank) == tk_ID),"B"])
			
			# adjust consumption flow based on resource proportions 
			if(flows_list[[k]][i,selected_col]!=0){ # check for none zero interaction 
				total_consump <- compartments_df[which(as.character(compartments_df$compartments) == OC[j] & as.character(compartments_df$tank) == tk_ID),"C"] # total consumption of OC[j] from compartments_df
				feed_pref <- (compartments_df[which(as.character(compartments_df$compartments) == compartment_names[i] & # feeding preference as proportion of total biomass contributed by compartment_name[i]
				as.character(compartments_df$tank) == tk_ID),"B"])/tot_bio
				flows_list[[k]][i,selected_col] <- total_consump * feed_pref # update flow: Consumption*Preference 
			}
		}
	}
}

## check the column sum equals the total consumption of consumers
for(i in 1:length(flows_list)){
	if(i==1)sum_col <- apply(flows_list[[i]],2,sum) # summarize colums sums (total consumption per consumer per tank), for the first tank start with 1, for others add tank numbers
	else sum_col <- rbind(sum_col, apply(flows_list[[i]],2,sum)) # 2: apply column wise
}

flows_list

#saveRDS(flows_list, file=file.path(output_path,"list_flows_matrices.rds"))

## summary of all flows
colnames(sum_col) <- compartment_names
rownames(sum_col) <- tank
sum_col



###################################################################
##										##
##   ADD columns Z (GPP), R (Respiration), E (Export/ Egestion)	##
 
##										##
##################################################################


flows_list_Z_R_E <- flows_list # create a new list to store flow matrices 
hete <- c("A_d", "A_o", "I_d", "I_o", "G_d", "G_h", "B_f", "P_d", "P_o") # heterotroph compartments in the flow matrices
l_hete <- length(hete) #  total number of heterotroph compartments in each tank

## add the egestion flows of heterotrophs to each flow matrix
for(k in 1:length(flows_list_Z_R_E)){
  for(i in 1:l_hete){
    egest <- compartments_df[which(as.character(compartments_df$tank) == tank[k] & as.character(compartments_df$compartments) == hete[i]),"E"]
    flows_list_Z_R_E[[k]][hete[i],"POC"] <- egest
  }
}

auto <- c("ZM", "FV", "FA") # primary producer compartments in the flow matrices (autotrophs)
l_auto <- length(auto) # total number of autotroph compartments in each tank
living <- c(auto, hete) # total living compartments in the flow matrices
l_living <- length(living) # total number of living compartments in each tank


## completing the database with three vectors summarizing import (Z = GPP), respiration (R) and export (E) flows
Z <- R <- E <- rep(0,length(all_compartments_old_order))
flows_list_Z_R_E_1 <- flows_list_Z_R_E

for(k in 1:length(flows_list_Z_R_E)){
  flows_list_Z_R_E[[k]] <- cbind(flows_list_Z_R_E_1[[k]], Z, R, E)
  ##
  ## add GPP and export for the three autotroph compartments
  for(i in 1:l_auto){
    GPP <- compartments_df[which(as.character(compartments_df$tank) == tank[k] & as.character(compartments_df$compartments) == auto[i]),"C"]
    flows_list_Z_R_E[[k]][auto[i],"Z"] <- GPP
  
    export <- macrophytes_exports[which(as.character(macrophytes_exports$tank) == tank[k] &
                                          as.character(macrophytes_exports$group) == auto[i]),"mgC_m2_day"]
    flows_list_Z_R_E[[k]][auto[i],"E"] <- export
  }
  ## add R for all living compartments
  for(j in 1:l_living){
    resp <- compartments_df[which(as.character(compartments_df$tank) == tank[k] & as.character(compartments_df$compartments) == living[j]),"R"]
    flows_list_Z_R_E[[k]][living[j],"R"] <- resp
  }
}

flows_list_Z_R_E

#saveRDS(flows_list_Z_R_E, file=file.path(output_path,"list_flows_Z_R_E_matrices.rds"))



