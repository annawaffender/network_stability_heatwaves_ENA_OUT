
## clean the workspace 

rm(list = ls())

## Set paths 
your_path <- "your_path/Github"  # replace with your path to Github folder 
source_file_path <-  file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_out/dataframes")
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_out/outputs")
setwd(output_path)


## load the ggplot2 library and the file with ENA functions
library("ggplot2")
source("ENAwithR.R")

## load the flux matrices and the three vectors of import, export & respiration
MF_raw <- readRDS("list_flows_Z_R_E_matrices.rds")
el <- length(MF_raw)

## flux matrices only (not needed for the balancing)
# MF_raw_onlyM <- readRDS("list_flows_matrices.rds")

## vector where treatments are stored
tr <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

## vector where tank names are stored
tn <- names(MF_raw)

## total number of compartments
co <- nrow(MF_raw[[1]])

################################
##							  ##
##   exports from consumers   ##
##							  ##
################################

## create the list where fluxes will be changed by changes the column of fluxes to detritus from consumers, which will become exports
MF_raw_NE <- as.list(rep(NA,el))

for(k in 1:el){ 
	MF_raw_NE[[k]] <- MF_raw[[k]]
	MF_raw_NE[[k]][4:co,"E"] <- MF_raw_NE[[k]][4:co,"POC"]
	MF_raw_NE[[k]][4:co,"POC"] <- rep(0,length(4:co))
}
names(MF_raw_NE) <- tn

## calculation of the four potential version to account for weighted fluxes as a proxy of connectance
HF_v <- rep(NA,el)
names(HF_v) <- tn
IC_v <- FW_v <- EC_v <- HF_v

MF_raw_NE_T <- as.list(rep(NA,el))

for(k in 1:el){
	MM <- MF_raw_NE[[k]]
	zero <- rep(0,nrow(MM))
	MAT <- MM[,1:nrow(MM)]
	colnames(MAT) <- rownames(MAT) <- rownames(MF_raw_NE[[k]])
	MF_raw_NE_T[[k]] <- MAT
	HF_v[k] <- DC(zero,MAT,zero,zero)/TST(zero,MAT,zero,zero)
	IC_v[k] <- I.C(MAT)
	FW_v[k] <- FW.C(MAT,1)
	EC_v[k] <- effconn.m(zero,MAT,zero,zero)
}

NS <- data.frame(tn,tr,HF_v,IC_v,FW_v,EC_v)
colnames(NS) <- c("tank","hw","HF","IC","FWC","EC")
NS$hw <- factor(NS$hw, levels = c("0HW", "1HW", "3HW"))


## save the RDS object
names(MF_raw_NE) <- names(MF_raw_NE_T) <- tn
#saveRDS(MF_raw_NE, file = "MF_raw_NE_ZTER.rds")
#saveRDS(MF_raw_NE_T, file = "MF_raw_NE_T.rds")



################################
##					   ##
## ENA: OUT.balance function ##
##					   ##
###############################

M_out <- M_out_T <- as.list(rep(NA,el))

## vectors where to store the indices obtained using the matrices balanced with the output algorithm
HF_ou <- rep(NA,el)
names(HF_ou) <- tn
IC_ou <- FW_ou <- HF_ou

for(k in 1:el){
	MM <- MF_raw_NE[[k]]
	ZZZ <- MM[,"Z"]
	MAT <- MM[,1:nrow(MM)]
	EEE <- MM[,"E"]
	RRR <- MM[,"R"]
	
	zero <- rep(0,nrow(MAT))
	
	## output-based balancing algorithm
	OUBM <- OUT.balance(ZZZ,MAT,EEE,RRR)
	ZZZ_OU <- OUBM[nrow(OUBM),1:co]
	MAT_OU <- OUBM[1:co,1:co]
	EEE_OU <- OUBM[1:co,(co+1)]
	RRR_OU <- OUBM[1:co,(co+2)]
	colnames(MAT_OU) <- rownames(MAT_OU) <- rownames(MF_raw_NE[[k]])
	names(ZZZ_OU) <- names(EEE_OU) <- names(RRR_OU) <- rownames(MF_raw_NE[[k]])
	HF_ou[k] <- DC(zero,MAT_OU,zero,zero)/TST(zero,MAT_OU,zero,zero)
	IC_ou[k] <- I.C(MAT_OU)
	FW_ou[k] <- FW.C(MAT_OU,1)
	M_out[[k]] <- cbind(MAT_OU,ZZZ_OU,RRR_OU,EEE_OU)
	M_out_T[[k]] <- MAT_OU
	colnames(M_out[[k]]) <- colnames(MF_raw_NE[[k]])
	rownames(M_out[[k]]) <- rownames(MF_raw_NE[[k]])
	
}

## OUT.balance: best balancing algorithm to use, in the absence of other changes to the flux matrix (i.e., without fixing GPP and POC export). 
## It preserves all fluxes and maintains the relative differences in the matrix of intercompartmental exchanges
NS_OU <- data.frame(tn,tr,HF_ou,IC_ou,FW_ou)
colnames(NS_OU) <- c("tank","hw","HF","IC","FWC")

## save the RDS object
names(M_out) <- names(M_out_T) <- tn
#saveRDS(M_out, file = "M_out_ZTER.rds")
#saveRDS(M_out_T, file = "M_out_T.rds")

