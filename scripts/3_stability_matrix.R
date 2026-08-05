

#### Calculate local stability

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.package)
library(readxl) 
library(dplyr)
library(crayon)

## Set paths
your_path <- "your_path/Github/network_stability_heatwaves_ENA_M_OUT/"  # replace with your path to Github folder 
correct_path <-  file.path(your_path,"dataframes")
output_path <- file.path(your_path,"outputs")
setwd(correct_path)



############################################
##										 ##
##   load data frames and set parameters  ##
##										 ##
############################################
## Load data frames
compartments_B <- read.csv(file.path(correct_path,"compartments_biomass.csv")) # all biomasses per tank and taxa
flow_tank_matrices  <- readRDS(file.path(output_path,"M_out_T.rds")) # flow matrices
network_compartments <- read.csv(file.path(output_path,"df_network_compartments.csv")) # all data 

file_path_source_param <- paste0(correct_path, "/", "1_source_data.xlsx")
source_data_m <- as.data.frame(read_xlsx(file_path_source_param, sheet = "mortality"))
source_data_a <- as.data.frame(read_xlsx(file_path_source_param, sheet = "B_PB_alpha"))


## Set parameters
num_compartments <- length(unique(network_compartments$compartments))
num_treat <- length(unique(network_compartments$treat))
num_tank <- length(unique(network_compartments$tank))

# taxa <- unique(network_compartments$compartments)
taxa <- c("ZM", "FV", "FA", "A_d", "A_o", "I_d", "I_o", "G_d", "G_h", "B_f", "P_d", "P_o", "POC") # reorder taxa to make it all match 
tanks <- unique(network_compartments$tank)
treat <- unique(network_compartments$treat)




#############################
##
## BIOMASS MATRIX ##
##
#############################

## Sum biomass data for different species of the same taxa and tank 
compartments_sumB <- aggregate(B~tank+taxa,data=compartments_B,sum) # 

# Create empty biomass matrix (matrix_B_empty)
matrix_B_empty <- matrix(NA, nrow = num_tank, ncol = num_compartments)
rownames(matrix_B_empty) <- tanks
colnames(matrix_B_empty) <- taxa
matrix_B_empty

matrix_B <- matrix_B_empty

## Fill empty biomass matrix with the corresponding biomass sums 
for (i in seq_len(nrow(compartments_sumB))) {
  # Get the current row's values
  tank_id <- compartments_sumB$tank[i]
  taxa_id <- compartments_sumB$taxa[i]
  B_sum <- compartments_sumB$B[i]
  
  # Check if tank and taxa exist in matrix
  if (tank_id %in% rownames(matrix_B) && taxa_id %in% colnames(matrix_B)) {
    matrix_B[tank_id, taxa_id] <- B_sum
  }
}
matrix_B





#############################
##
##  PE - MATRIX  ##
##
#############################

# Calculate production efficiency (PE)
network_compartments$PE <- with(network_compartments, P / (P + R))

# Create empty PE matrix
matrix_PE <- matrix(NA,nrow = num_tank,ncol = num_compartments)
rownames(matrix_PE) <- tanks
colnames(matrix_PE) <- taxa


## Fill empty PE matrix with the corresponding PE values
for (i in seq_len(nrow(network_compartments))) {
  # Get the current row's values
  tank_id <- network_compartments$tank[i]
  comp_id <- network_compartments$compartments[i]
  pe_val <- network_compartments$PE[i]
  
  # Check if tank and taxa exist in matrix
  if (tank_id %in% rownames(matrix_PE) &&
      comp_id %in% colnames(matrix_PE)) {
    matrix_PE[tank_id, comp_id] <- pe_val
  }
}
matrix_PE



#############################
##
##  alpha - MATRIX ##
##
#############################

source_data_alpha <- network_compartments[, c("compartments","tank","treat","alpha")]
source_data_alpha

# Create empty alpha matrix (matrix_alpha_empty)
matrix_alpha_empty <- matrix(NA, nrow = num_tank, ncol = num_compartments)
rownames(matrix_alpha_empty) <- tanks
colnames(matrix_alpha_empty) <- taxa
matrix_alpha_empty

matrix_alpha <- matrix_alpha_empty



## Fill empty alpha matrix with the corresponding alpha values
for (i in seq_len(nrow(source_data_alpha))) {
  # Get the current row's values
  tank_id <- source_data_alpha$tank[i]
  compartments_id <- source_data_alpha$compartments[i]
  alpha_values <- source_data_alpha$alpha[i]
  
  # Check if tank and taxa exist in matrix
  if (tank_id %in% rownames(matrix_alpha) && compartments_id %in% colnames(matrix_alpha)) {
    matrix_alpha[tank_id, compartments_id] <- alpha_values
  }
}
matrix_alpha




#############################
##
##  MORTALITIES ##
##
#############################

## source data file: mortality 

## Zostera marina
ZM_m_0HW <- source_data_m[source_data_m$taxa == "ZM" & source_data_m$treatment == "0HW","m"]
ZM_m_1HW <- source_data_m[source_data_m$taxa == "ZM" & source_data_m$treatment == "1HW","m"]
ZM_m_3HW <- source_data_m[source_data_m$taxa == "ZM" & source_data_m$treatment == "3HW","m"]

## Fucus vesiculosus 
FV_m_0HW <- source_data_m[source_data_m$taxa == "FV" & source_data_m$treatment == "0HW","m"]
FV_m_1HW <- source_data_m[source_data_m$taxa == "FV" & source_data_m$treatment == "1HW","m"]
FV_m_3HW <- source_data_m[source_data_m$taxa == "FV" & source_data_m$treatment == "3HW","m"]

## filamentous algae ´ 
FA_m_0HW <- source_data_m[source_data_m$taxa == "FA" & source_data_m$treatment == "0HW","m"]
FA_m_1HW <- source_data_m[source_data_m$taxa == "FA" & source_data_m$treatment == "1HW","m"]
FA_m_3HW <- source_data_m[source_data_m$taxa == "FA" & source_data_m$treatment == "3HW","m"]

## Amphipod o. - Gammarus
A_o_m_0HW <- source_data_m[source_data_m$taxa == "A_o" & source_data_m$treatment == "0HW","m"]
A_o_m_1HW <- source_data_m[source_data_m$taxa == "A_o" & source_data_m$treatment == "1HW","m"]
A_o_m_3HW <- source_data_m[source_data_m$taxa == "A_o" & source_data_m$treatment == "3HW","m"]

## Isopod o.- Idotea balthica
I_o_m_0HW <- source_data_m[source_data_m$taxa == "I_o" & source_data_m$treatment == "0HW","m"]
I_o_m_1HW <- source_data_m[source_data_m$taxa == "I_o" & source_data_m$treatment == "1HW","m"]
I_o_m_3HW <- source_data_m[source_data_m$taxa == "I_o" & source_data_m$treatment == "3HW","m"]

## Gastropod h. - Littorina littorea
G_h_m_0HW <- source_data_m[source_data_m$taxa == "G_h" & source_data_m$treatment == "0HW","m"]
G_h_m_1HW <- source_data_m[source_data_m$taxa == "G_h" & source_data_m$treatment == "1HW","m"]
G_h_m_3HW <- source_data_m[source_data_m$taxa == "G_h" & source_data_m$treatment == "3HW","m"]

## Amphipod d. - Corophium volutator
A_d_m_0HW <- source_data_m[source_data_m$taxa == "A_d" & source_data_m$treatment == "0HW","m"]
A_d_m_1HW <- source_data_m[source_data_m$taxa == "A_d" & source_data_m$treatment == "1HW","m"]
A_d_m_3HW <- source_data_m[source_data_m$taxa == "A_d" & source_data_m$treatment == "3HW","m"]

## Gastropod d. - Hydrobia ulvae
G_d_m_0HW <- source_data_m[source_data_m$taxa == "G_d" & source_data_m$treatment == "0HW","m"]
G_d_m_1HW <- source_data_m[source_data_m$taxa == "G_d" & source_data_m$treatment == "1HW","m"]
G_d_m_3HW <- source_data_m[source_data_m$taxa == "G_d" & source_data_m$treatment == "3HW","m"]

## Isopod d. - Jaera albifrons
I_d_m_0HW <- source_data_m[source_data_m$taxa == "I_d" & source_data_m$treatment == "0HW","m"]
I_d_m_1HW <- source_data_m[source_data_m$taxa == "I_d" & source_data_m$treatment == "1HW","m"]
I_d_m_3HW <- source_data_m[source_data_m$taxa == "I_d" & source_data_m$treatment == "3HW","m"]

## Bivale f. - 	Mytilus edulis
B_f_m_0HW <- source_data_m[source_data_m$taxa == "B_f" & source_data_m$treatment == "0HW","m"]
B_f_m_1HW <- source_data_m[source_data_m$taxa == "B_f" & source_data_m$treatment == "1HW","m"]
B_f_m_3HW <- source_data_m[source_data_m$taxa == "B_f" & source_data_m$treatment == "3HW","m"]

## Polychaeta d. - 	Marenzelleria sp.
P_d_m_0HW <- source_data_m[source_data_m$taxa == "P_d" & source_data_m$treatment == "0HW","m"]
P_d_m_1HW <- source_data_m[source_data_m$taxa == "P_d" & source_data_m$treatment == "1HW","m"]
P_d_m_3HW <- source_data_m[source_data_m$taxa == "P_d" & source_data_m$treatment == "3HW","m"]

## Polychaeta o. - 	 Hediste diversicolor
P_o_m_0HW <- source_data_m[source_data_m$taxa == "P_o" & source_data_m$treatment == "0HW","m"]
P_o_m_1HW <- source_data_m[source_data_m$taxa == "P_o" & source_data_m$treatment == "1HW","m"]
P_o_m_3HW <- source_data_m[source_data_m$taxa == "P_o" & source_data_m$treatment == "3HW","m"]

## POC
POC_m_0HW <- source_data_m[source_data_m$taxa == "POC" & source_data_m$treatment == "0HW","m"]
POC_m_1HW <- source_data_m[source_data_m$taxa == "POC" & source_data_m$treatment == "1HW","m"]
POC_m_3HW <- source_data_m[source_data_m$taxa == "POC" & source_data_m$treatment == "3HW","m"]


## Summarize the mortalities of each treatment
mortality_vector_0HW <- c(ZM_m_0HW,FV_m_0HW,FA_m_0HW,A_d_m_0HW,A_o_m_0HW,I_d_m_0HW,I_o_m_0HW,G_d_m_0HW,G_h_m_0HW,B_f_m_0HW,P_d_m_0HW,P_o_m_0HW,POC_m_0HW)        
mortality_vector_1HW <- c(ZM_m_1HW,FV_m_1HW,FA_m_1HW,A_d_m_1HW,A_o_m_1HW,I_d_m_1HW,I_o_m_1HW,G_d_m_1HW,G_h_m_1HW,B_f_m_1HW,P_d_m_1HW,P_o_m_1HW,POC_m_1HW)
mortality_vector_3HW <- c(ZM_m_3HW,FV_m_3HW,FA_m_3HW,A_d_m_3HW,A_o_m_3HW,I_d_m_3HW,I_o_m_3HW,G_d_m_3HW,G_h_m_3HW,B_f_m_3HW,P_d_m_3HW,P_o_m_3HW,POC_m_3HW)


## Create a mortality matrix
mortality_matrix <- rbind(mortality_vector_0HW,mortality_vector_1HW ,mortality_vector_3HW)
colnames(mortality_matrix) <- taxa
rownames(mortality_matrix) <- treat
mortality_matrix





#############################
##
## INTERACTION STRENGTHS  ##
##
#############################

## Create an empty matrix: number compartments x number compartments 
compartments_matrix <- matrix(rep(0,num_compartments^2),nrow = num_compartments)
colnames(compartments_matrix) <- rownames(compartments_matrix) <- taxa
compartments_matrix

## Create an list, using the compartments_matrix: interaction matrices - a list of 11 empty matrices, one for each tank to be filled with the interactions 
orig_interaction_tank_matrices <- as.list(rep(NA,num_tank))
for(i in 1:length(orig_interaction_tank_matrices)) orig_interaction_tank_matrices[[i]] <- compartments_matrix
names(orig_interaction_tank_matrices) <- tanks
orig_interaction_tank_matrices


## Fill list of interaction matrices  
for (i in 1:length(flow_tank_matrices)) {      # loop through each tank matrix
  current_matrix <- flow_tank_matrices[[i]]
  tank <- names(flow_tank_matrices)[i]         # Retrieve tank names at index i 
  
  for (row in 1:nrow(current_matrix)) {        # Loop through each row: prey species 
    for(col in 1:ncol(current_matrix)){        # Loop through each column: predator species
      if (current_matrix[row, col]!=0){        # If hit non-zero-value, continue with calculations 
        B_prey <- matrix_B[tank, row]          # Retrieve biomass of prey from each tank and species (row) 
        B_pred <- matrix_B[tank, col]          # Retrieve biomass of predator from each tank and species (col) 
        PE_pred <- matrix_PE[tank, col]        # Retrieve PE values
        alpha_pred <- matrix_alpha[tank, col]  # retrieve alpha values
        Flow <- current_matrix[row, col]       # Extract flow between prey and pred. from current_matrix 
        
        orig_interaction_tank_matrices[[i]][row, col] <- (Flow/B_prey) *alpha_pred *PE_pred   # Update previously created orig_interaction_tank_matrices with +impact of prey on pred  * multiply by assimilation efficiency (see: Gaedke etal.; equ.2) 
        orig_interaction_tank_matrices[[i]][col, row] <- -Flow/B_pred}            # Update previously created orig_interaction_tank_matrices with -impact of pred on prey 
    }}}

orig_interaction_tank_matrices

interaction_tank_matrices<- orig_interaction_tank_matrices


#############################
##
## INTERACTION MATRICES  ##
##
#############################

## Assign mortalities to the diagonals of the interaction matrices
ref_tanks <- tanks
ref_treatments <- network_compartments$treat[1:11]


for (i in 1:length(interaction_tank_matrices)) {            # Loop through the list of 11 matrices 
  tank <- names(interaction_tank_matrices)[i]               # Retrieve tank names at index i 
  row_i <- which(ref_tanks==tank)                           # Find index in ref_tanks, where tank name matches the current processed tank
  row_i2 <- ref_treatments[row_i]                           # Retrieve corresponding treatment type from ref_treatment using the index row_i 
  diag(interaction_tank_matrices[[i]])<- -mortality_matrix[which(rownames(mortality_matrix) ==row_i2), ] # set -mortality values as diagonal of each i matrix 
}

interaction_tank_matrices


## override the intrinsic mortality values of P_o with neg. interactions from cannibalism 
for (i in seq_along(interaction_tank_matrices)) {
  mat <- interaction_tank_matrices[[i]]
  orig_mat <- orig_interaction_tank_matrices[[i]]
  
  p_o_index <- which(colnames(mat) == "P_o")   # Find the index of the "P_o" column

  if (length(p_o_index) == 1) {
    diag(mat)[p_o_index] <- diag(orig_mat)[p_o_index]     # Replace the diagonal element for P_o with original value
    interaction_tank_matrices[[i]] <- mat     # Update the matrix back in the list
  } else {
    warning(paste("P_o column not found in matrix", i))
  }
}

interaction_tank_matrices
#saveRDS(interaction_tank_matrices, file=file.path(output_path,"list_interactions_matrices.rds"))



##############################
##
## STABILITY CALUCLATIONS    ##
## WITH CHANGING MORTALITIES ##
##
##############################


multipliers <- seq(from = 0.015, to = 0.001, by = -0.0001)
num_mulipliers <- length(multipliers)

## Create stability matrix: one empty matrix, to be filled with the stability results  
stability_matrix <- matrix(NA, nrow=num_mulipliers, ncol=num_tank)
colnames(stability_matrix) <- tanks
rownames(stability_matrix) <- multipliers
stability_matrix 


for (i in 1:length(multipliers )){                     
  for(j in 1:length(interaction_tank_matrices)) {       # Loop through each interaction matrix, one for each tank 
    interaction_matrix_j <- interaction_tank_matrices[[j]]
    diag(interaction_matrix_j) <- diag(interaction_matrix_j)*multipliers[i]  # multiply (mortalities) in the matrix diagonal with the multipliers
    
    eigenvalues <- eigen(interaction_matrix_j)$values  # Compute eigenvalues ($values) of interaction_matrix_j
    {         
      if(length(which(Re(eigenvalues)>0)>=1)) {        
        stability_matrix[i,j] <- 0                  # If at least one eigenvalue is positive: assign 0 -> system is unstable
      }
      else {
        stability_matrix[i,j] <- 1                  # If all eigenvalues are non-positive:assign 1 -> system is stable
      }}}}

stability_matrix

#write.csv(stability_matrix, file=file.path(output_path,"stability_matrix.csv"))

