

#### Create dataframe of all macrophytes (B, PB, RB, NPP, R, GPP) 



rm(list=ls(all=TRUE))	

## load needed packages 
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

## Set paths 
your_path <- "your_path/Github" # replace with your path to Github folder 


###########################################
##										 ##
##   load data frames  ##
##										 ##
###########################################

## B (biomass): macrophytes and the nonliving compartment (particulate organic carbon: POC)
file_path_macrophytes_B <- paste0(correct_path, "/", "macrophytes_biomass_Ito_etal2024.xlsx")
macrophytes_B_df <- as.data.frame(read_xlsx(file_path_macrophytes_B, sheet = "biomass_complete"))

## GPP (general primary production ), R (respiration) 
file_path_macrophytes_M <- paste0(correct_path, "/", "macrophytes_NPP_RESP_Ito_etal2024.xlsx")
macrophytes_M_df <- as.data.frame(read_xlsx(file_path_macrophytes_M, sheet = "dataset_macrophytes"))

## combined database
macrophtyes_nonliving_df <- data.frame(macrophytes_B_df$species, macrophytes_B_df$taxa,
macrophytes_B_df$tank, macrophytes_B_df$treatment, macrophytes_B_df$tank_BC_mg_m2)
colnames(macrophtyes_nonliving_df) <- c("species", "taxa", "tank", "treat", "B")

tot_el_PPNL <- nrow(macrophtyes_nonliving_df)
PB_v <- RB_v <- rep(NA,tot_el_PPNL)

macro_rows <- 36 # total number of rows in the database corresponding to the three macrophytes in the 12 tanks
POC_rows <- 11 # total number of rows in the database corresponding to POC in the 11 tanks (tank B1 data unavailable)

for(i in 1:tot_el_PPNL){
	taxa_i <- macrophtyes_nonliving_df[i,"taxa"]
	tank_i <- macrophtyes_nonliving_df[i,"tank"]

	row_mnl <- which(macrophytes_M_df$group == taxa_i & macrophytes_M_df$tank == tank_i)

	if(length(row_mnl)!=0){
		PB_v[i] <- macrophytes_M_df[row_mnl,"NPP_mgC_mgC_day"]
		RB_v[i] <- macrophytes_M_df[row_mnl,"RESP_mgC_mgC_day"]
	}
}
macrophtyes_nonliving_df1 <- data.frame(macrophtyes_nonliving_df, PB_v, RB_v)
colnames(macrophtyes_nonliving_df1) <- c("species", "taxa", "tank", "treat", "B", "PB", "RB")

## add average treatment values in the case of missing NPP (net primary production) and R (respiration) rates for FV (Fucus versicolor) and FA (filamentous algae)
PP_POC_df <- macrophtyes_nonliving_df1
NPP_NA <- R_NA <- matrix(rep(NA,6), nrow = 3)
colnames(NPP_NA) <- colnames(R_NA) <- c("FV", "FA")
rownames(NPP_NA) <- rownames(R_NA) <- c("0HW", "1HW", "3HW")


## Fucus vesciculosus (FV)
NPP_NA["0HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "0HW"),"PB"],na.rm=TRUE)
NPP_NA["1HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "1HW"),"PB"],na.rm=TRUE)
NPP_NA["3HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "3HW"),"PB"],na.rm=TRUE)

R_NA["0HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "0HW"),"RB"],na.rm=TRUE)
R_NA["1HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "1HW"),"RB"],na.rm=TRUE)
R_NA["3HW","FV"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FV" & PP_POC_df$treat == "3HW"),"RB"],na.rm=TRUE)

## filamentous algae (FA)
NPP_NA["0HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "0HW"),"PB"],na.rm=TRUE)
NPP_NA["1HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "1HW"),"PB"],na.rm=TRUE)
NPP_NA["3HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "3HW"),"PB"],na.rm=TRUE)
##
R_NA["0HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "0HW"),"RB"],na.rm=TRUE)
R_NA["1HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "1HW"),"RB"],na.rm=TRUE)
R_NA["3HW","FA"] <- mean(PP_POC_df[which(PP_POC_df$taxa == "FA" & PP_POC_df$treat == "3HW"),"RB"],na.rm=TRUE)

for(i in 1:macro_rows){
	if(is.na(macrophtyes_nonliving_df1[i,"PB"])==TRUE){
		ta_n <- as.character(macrophtyes_nonliving_df1[i,"taxa"])
		tr_n <- as.character(macrophtyes_nonliving_df1[i,"treat"])
		macrophtyes_nonliving_df1[i,"PB"] <- NPP_NA[tr_n,ta_n]
	}
	##
	if(is.na(macrophtyes_nonliving_df1[i,"RB"])==TRUE){
		ta_n <- as.character(macrophtyes_nonliving_df1[i,"taxa"])
		tr_n <- as.character(macrophtyes_nonliving_df1[i,"treat"])
		macrophtyes_nonliving_df1[i,"RB"] <- R_NA[tr_n,ta_n]
	}	
}


###########################################
##										 ##
##   merge data and add all parameters ##
##										 ##
###########################################

## calculation of the NPP, starting from P/B (PB) ratios and taxa biomass
# NPP= PB*B; R= RB*B --> flows 
# NPP = P 
NPP_v0 <- macrophtyes_nonliving_df1$PB[1:macro_rows] * macrophtyes_nonliving_df1$B[1:macro_rows]
macrophtyes_nonliving_df1$NPP <- c(NPP_v0, rep(NA,POC_rows))

## calculation of the R, starting from R/B ratios and taxa biomass
R_v0 <- macrophtyes_nonliving_df1$RB[1:macro_rows] * macrophtyes_nonliving_df1$B[1:macro_rows]
macrophtyes_nonliving_df1$R <- c(R_v0, rep(NA,POC_rows))

## calculation of the GPP = NPP + R
GPP_v0 <- NPP_v0 + R_v0
macrophtyes_nonliving_df1$GPP <- c(GPP_v0, rep(NA,POC_rows))

#write.csv(macrophtyes_nonliving_df1, file =file.path(output_path,"df_macrophytes_POC.csv"), row.names = FALSE)
