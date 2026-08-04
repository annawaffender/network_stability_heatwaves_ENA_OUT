
#### Sensitivity analysis: interaction simulations & stability calculations 

rm(list = ls(all = TRUE))  # clear workspace 

## load needed packages (eventually install.package)
library(dplyr)

## Set paths
your_path <- "your_path/Github"  # replace with your path to Github folder 
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_out/outputs")
setwd(output_path)


###################################
##						##
##   Load data & set parameters ##
##						##
###################################

flow_tank_matrices <- readRDS(file.path(output_path,"M_out_T.rds")) # flow matrices
interaction_tank_matrices <- readRDS(file.path(output_path,"list_interactions_matrices.rds"))
stability_matrix <- read.csv(file.path(output_path,"stability_matrix.csv"),row.names = 1)

set.seed(1000)
runs <- 1000

treat <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

## reorder to match compartment order in both matrices
compartments_reordered <- c("ZM", "FV", "FA", "A_d", "A_o", "I_d", "I_o", "G_d", "G_h", "B_f", "P_d", "P_o", "POC") # see script 3_flow_matrices
orig_interaction_matrices <- lapply(interaction_tank_matrices, function(mat) {mat[, compartments_reordered]})




###################################
##						##
##   Mortality multiplier ##
##						##
###################################

## set mortality_multipliers: the last stable multiplier for each tank  
mortality_multipliers <- apply(stability_matrix, 2, function(col) {rownames(stability_matrix)[which(diff(col) == -1)[1]]}) # Find smallest row name where each column transitions from 1 to 0: apply function to each column (2= col wise); diff(col) == -1): logical vector where transitions occurs 
mortality_multipliers <- setNames(as.numeric(mortality_multipliers), names(mortality_multipliers)) # Convert to vector

all_stability_results <- data.frame(tank = character(0), run = integer(0),row = integer(0), column = integer(0), stability = integer(0)) # Initialize data frame to store all results for all tanks
interaction_matrices <- orig_interaction_matrices 

## Apply the unique multipliers to the diagonal of each single_matrix (mortalities)
for(tank in names(interaction_matrices)) {
  single_matrix <- orig_interaction_matrices[[tank]]
  x <- mortality_multipliers[[tank]]
  diag(interaction_matrices[[tank]]) <- diag(single_matrix) * x ## * 1.05 --> Adjust diagonal; increase multiplier value by 1% (1.01)/ 5% (1.05)/ 10% (1.10)/ ...
}



###################################
##						##
##   Simulation runs   ##
##						##
###################################

## same constant change at each run for positive and negative effects of the interaction strengths 
## constant values to multiply intervals of the interactions strengths

upper_m <- 1.5 
lower_m <- 0.5 

## simulations with modified pairwise interactions 
for (tank in names(interaction_matrices)) {                # Loop through each matrix in interaction_matrices
  orig_interaction_matrix <- interaction_matrices[[tank]]  # keep original matrix to reset after each run
  orig_flow_matrix <- flow_tank_matrices[[tank]]
  modified_matrix <- interaction_matrices[[tank]]
  
  # modified_interactions <- list()  # Initialize a list to store modified versions of this matrix for each run
  
  nonzero_flows <- length(which(orig_flow_matrix != 0))           # Nonzero flows (interactions)
  # stability_sensitivity <- matrix(rep(0, nonzero_flows * runs), nrow = nonzero_flows)  # Initialize matrix for stability results
  
  count <- 1  # Keeps track of rows in summary matrix
  
  for (i in 1:nrow(orig_flow_matrix)) {  # Loop through rows (i)
    for (j in 1:ncol(orig_flow_matrix)) {  # Loop through columns (j)
      if (orig_flow_matrix[i, j] != 0) {  # only include non-zero elements of flow matrix # if avoid testing the impact of the self-loop P_o-P_o add:& i != j
        for(run in 1:runs){ 
          
          # Modify matrix values
          multiplier <- runif(1, min = lower_m, max = upper_m)
          
		      modified_matrix[i, j] <- multiplier * abs(orig_interaction_matrix[i, j])
          modified_matrix[j, i] <- -1 * multiplier * abs(orig_interaction_matrix[j, i])

          # Check eigenvalues 
          eigen_values <- eigen(modified_matrix)$values
          
          # Check stability: if at least one eigenvalue has a positive real part, system is unstable
          stability <- ifelse(length(which(Re(eigen_values) > 0)) >= 1, 1, 0)  # 1=unstable, 0=stable
          
          all_stability_results <- rbind(all_stability_results, data.frame(tank = tank,run = run, row = i, column = j, stability = stability, multiplier = multiplier))  # Store all results 
          
        }
        
        count <- count + 1  # Increment count after each interaction
        modified_matrix <- interaction_matrices[[tank]]  # Reset matrix to original for next iteration
      }
    }
  }
}

## Print final stability results and interaction table
all_stability_results  # if column stability: 0 = stable; 1 =unstable 




###################################
##						##
##   Summarize outputs     ##
##						##
###################################

## interaction - stability outputs 
sum_stability_results <- as.data.frame(all_stability_results %>% group_by(tank, row, column) %>%
                                         summarise(stable = sum(stability == 0),unstable = sum(stability == 1),.groups = "drop"))


# rename row and column elements to replace numbers by species; based on flow_tank_matrices or orig_interaction_matrices
final_interaction_stability <- sum_stability_results

final_interaction_stability$row <- ifelse(final_interaction_stability$row == 1, "ZM",
                                          ifelse(final_interaction_stability$row == 2, "FV",
                                                 ifelse(final_interaction_stability$row == 3, "FA",
                                                        ifelse(final_interaction_stability$row == 4, "A_d", 
                                                               ifelse(final_interaction_stability$row == 5, "A_o", 
                                                                      ifelse(final_interaction_stability$row == 6, "I_d",
                                                                             ifelse(final_interaction_stability$row == 7, "I_o",
                                                                                    ifelse(final_interaction_stability$row == 8, "G_d",       
                                                                                           ifelse(final_interaction_stability$row == 9, "G_h",
                                                                                                  ifelse(final_interaction_stability$row == 10, "B_f",       
                                                                                                         ifelse(final_interaction_stability$row == 11, "P_d", 
                                                                                                                ifelse(final_interaction_stability$row == 12, "P_o",
                                                                                                                       ifelse(final_interaction_stability$row ==13,"POC",
                                                                                                                              final_interaction_stability$row)))))))))))))

final_interaction_stability$column <-  ifelse(final_interaction_stability$column == 1, "ZM",
                                              ifelse(final_interaction_stability$column == 2, "FV",
                                                     ifelse(final_interaction_stability$column == 3, "FA",
                                                            ifelse(final_interaction_stability$column == 4, "A_d", 
                                                                   ifelse(final_interaction_stability$column == 5, "A_o", 
                                                                          ifelse(final_interaction_stability$column == 6, "I_d",
                                                                                 ifelse(final_interaction_stability$column == 7, "I_o",
                                                                                        ifelse(final_interaction_stability$column == 8, "G_d",       
                                                                                               ifelse(final_interaction_stability$column == 9, "G_h",
                                                                                                      ifelse(final_interaction_stability$column == 10, "B_f",       
                                                                                                             ifelse(final_interaction_stability$column == 11, "P_d", 
                                                                                                                    ifelse(final_interaction_stability$column == 12, "P_o",
                                                                                                                           ifelse(final_interaction_stability$column ==13,"POC",
                                                                                                                                  final_interaction_stability$column)))))))))))))

final_interaction_stability
#write.csv(final_interaction_stability,"stability_sensitivity_interactions.csv",row.names = FALSE)


