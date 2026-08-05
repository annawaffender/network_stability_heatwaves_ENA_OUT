
#### Figure 3: Body size - frequency distribution  
#### Plots and Wasserstein distance and permutation test analysis 

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.packages)
library(readxl)
library(readr)
library(ggplot2)
library(dplyr)
library(transport)


## Set paths
your_path <- "your_path/Github/network_stability_heatwaves_ENA_M_OUT/"  # replace with your path to Github folder 
output_path <- file.path(your_path,"outputs")


####################################################
##						##
##  STATISTICS: MAIN MESOGRAZER: Sizes frequency  ##
##						##
####################################################

## Load and prepare data 
MM_individual_path <- paste0(output_path, "/", "df_mesograzers_individ.csv")
MM_individual_df <- as.data.frame(read.csv(MM_individual_path))

MM_average_path <- paste0(output_path, "/", "df_mesograzers_tank.csv")
MM_average_df <- as.data.frame(read.csv(MM_average_path))

MM_individual_df$log_length_mm <- log(MM_individual_df$length_mm) # add log transformation of length_mm
MM_average_df$log_C <- log(MM_average_df$C) # add log transformation of Consumption 
MM_average_df$log_R <- log(MM_average_df$R) # add log transformation of Respiration  
MM_average_df$log_E <- log(MM_average_df$E) # add log transformation of Egestion  

## Filter data for each species for species based statistics 
Ao_individual_data <- MM_individual_df %>% filter(species == "A_o")
Gh_individual_data <- MM_individual_df %>% filter(species == "G_h")
Io_individual_data <- MM_individual_df %>% filter(species == "I_o")


## Filter data for species, treatment specific body sizes 
Ao_size_0HW <- Ao_individual_data$length_mm[Ao_individual_data$treatment=="0HW"]
Ao_size_1HW <- Ao_individual_data$length_mm[Ao_individual_data$treatment=="1HW"]
Ao_size_3HW <- Ao_individual_data$length_mm[Ao_individual_data$treatment=="3HW"]

Gh_size_0HW <- Gh_individual_data$length_mm[Gh_individual_data$treatment=="0HW"]
Gh_size_1HW <- Gh_individual_data$length_mm[Gh_individual_data$treatment=="1HW"]
Gh_size_3HW <- Gh_individual_data$length_mm[Gh_individual_data$treatment=="3HW"]

Io_size_0HW <- Io_individual_data$length_mm[Io_individual_data$treatment=="0HW"]
Io_size_1HW <- Io_individual_data$length_mm[Io_individual_data$treatment=="1HW"]
Io_size_3HW <- Io_individual_data$length_mm[Io_individual_data$treatment=="3HW"]




#############################################
##						##
##   Wasserstein distance & permutation test ##
##						##
#############################################

## 1-Wasserstein distance (Earth Mover's Distance): distance 1 --> one dimensional distance (true for most trait value distributions)

## MAIN MESOGRAZER 
## Initialize a results table

results_table <- data.frame(Comparison = character(),Wasserstein_Distance = numeric(),p_Value = numeric(),stringsAsFactors = FALSE)

#######################
##
## Amphipod omnivore ##
##
#######################

## Ao_individual data: 0HW ~1HW ------------------------------------------------------
## compute the 1-Wasserstein distance
wasserstein_Ao_0HW_1HW <- wasserstein1d(Ao_size_0HW,Ao_size_1HW)
print(wasserstein_Ao_0HW_1HW)

## set num_permutations to generate the null distribution (= what the distance would look like if there's no group difference)
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## shuffle and split the combined data (the 2 HW distributions) into two new groups with the same sizes as originals
## calculate wasserstein distances for shuffled groups 
## store wasserstein distances of shuffles distributions in null_distribution 
for(i in 1:num_permutations){
  full_distribution <- c(Ao_size_0HW,Ao_size_1HW)
  indices <- sample(1:length(full_distribution), length(Ao_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-value: proportion of permutated distances greater or equal to observed distances
## if < 0.05 = observed wasserstein distance significantly different 
p_value_Ao_0HW_1HW <- sum(wasserstein_Ao_0HW_1HW <= null_distribution) / num_permutations 
p_value_Ao_0HW_1HW

## store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Ao_0HW_vs_1HW", Wasserstein_Distance = wasserstein_Ao_0HW_1HW, p_Value = p_value_Ao_0HW_1HW))


## Ao_individual data: 0HW ~3HW ------------------------------------------------------
## Compute the 1-Wasserstein distance 
wasserstein_Ao_0HW_3HW <- wasserstein1d(Ao_size_0HW,Ao_size_3HW)
print(wasserstein_Ao_0HW_3HW)

## set up num_permutations 
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## permutation test of wasserstein distances  
for(i in 1:num_permutations){
  full_distribution <- c(Ao_size_0HW,Ao_size_3HW)
  indices <- sample(1:length(full_distribution), length(Ao_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-values
p_value_Ao_0HW_3HW <- sum(wasserstein_Ao_0HW_3HW <= null_distribution) / num_permutations 
p_value_Ao_0HW_3HW

## store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Ao_0HW_vs_3HW", Wasserstein_Distance = wasserstein_Ao_0HW_3HW, p_Value = p_value_Ao_0HW_3HW))


#######################
##
## Gastropod herbivore ##
##
#######################


## Gh_individual data: 0HW ~1HW ------------------------------------------------------
## compute the 1-Wasserstein distance 
wasserstein_Gh_0HW_1HW <- wasserstein1d(Gh_size_0HW,Gh_size_1HW)
print(wasserstein_Gh_0HW_1HW)

## set up num_permutations 
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## permutation test of wasserstein distances  
for(i in 1:num_permutations){
  full_distribution <- c(Gh_size_0HW,Gh_size_1HW)
  indices <- sample(1:length(full_distribution), length(Gh_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-value
p_value_Gh_0HW_1HW <- sum(wasserstein_Gh_0HW_1HW <= null_distribution) / num_permutations 
p_value_Gh_0HW_1HW

## store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Gh_0HW_vs_1HW", Wasserstein_Distance = wasserstein_Gh_0HW_1HW, p_Value = p_value_Gh_0HW_1HW))



## Gh_individual data: 0HW ~3HW ------------------------------------------------------
## compute the 1-Wasserstein distance 
wasserstein_Gh_0HW_3HW <- wasserstein1d(Gh_size_0HW,Gh_size_3HW)
print(wasserstein_Gh_0HW_3HW)

## set up num_permutations 
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## permutation test of wasserstein distances  
for(i in 1:num_permutations){
  full_distribution <- c(Gh_size_0HW,Gh_size_3HW)
  indices <- sample(1:length(full_distribution), length(Gh_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-value
p_value_Gh_0HW_3HW <- sum(wasserstein_Gh_0HW_3HW <= null_distribution) / num_permutations 
p_value_Gh_0HW_3HW

## store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Gh_0HW_vs_3HW", Wasserstein_Distance = wasserstein_Gh_0HW_3HW, p_Value = p_value_Gh_0HW_3HW))
results_table



#######################
##
## Isopod omnivore ##
##
#######################

## Io_individual data: 0HW ~1HW ------------------------------------------------------
## compute the 1-Wasserstein distance
wasserstein_Io_0HW_1HW <- wasserstein1d(Io_size_0HW,Io_size_1HW)
print(wasserstein_Io_0HW_1HW)

## set up num_permutations 
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## permutation test of wasserstein distances  
for(i in 1:num_permutations){
  full_distribution <- c(Io_size_0HW,Io_size_1HW)
  indices <- sample(1:length(full_distribution), length(Io_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-value
p_value_Io_0HW_1HW <- sum(wasserstein_Io_0HW_1HW <= null_distribution) / num_permutations 
p_value_Io_0HW_1HW

## store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Io_0HW_vs_1HW", Wasserstein_Distance = wasserstein_Io_0HW_1HW, p_Value = p_value_Io_0HW_1HW))
results_table


## Io_individual data: 0HW ~3HW ------------------------------------------------------
## compute the 1-Wasserstein distance 
wasserstein_Io_0HW_3HW <- wasserstein1d(Io_size_0HW,Io_size_3HW)
print(wasserstein_Io_0HW_3HW)

## set up num_permutations 
num_permutations <- 10000
null_distribution <- rep(NA, num_permutations)

## permutation test of wasserstein distances  
for(i in 1:num_permutations){
  full_distribution <- c(Io_size_0HW,Io_size_3HW)
  indices <- sample(1:length(full_distribution), length(Io_size_0HW), replace = F)
  x_HW <- full_distribution[indices]
  y_HW <- full_distribution[-indices]
  null_distribution[i] <- wasserstein1d(x_HW, y_HW) 
}

## calculate p-value
p_value_Io_0HW_3HW <- sum(wasserstein_Io_0HW_3HW <= null_distribution) / num_permutations 
p_value_Io_0HW_3HW

# store result in the table
results_table <- rbind(results_table, data.frame(Comparison = "Io_0HW_vs_3HW", Wasserstein_Distance = wasserstein_Io_0HW_3HW, p_Value = p_value_Io_0HW_3HW))
results_table


# count decimal places
#results_table$p_Value <- round(results_table$p_Value, 6)



##############################################
##						##
## PLOT: MAIN MESOGRAZER: Sizes frequency  ##
##						##
#############################################

## Load and prepare data 
MM_individual_path <- paste0(output_path, "/", "df_mesograzers_individ.csv")
MM_individual_df <- as.data.frame(read.csv(MM_individual_path))

MM_average_path <- paste0(output_path, "/", "df_mesograzers_tank.csv")
MM_average_df <- as.data.frame(read.csv(MM_average_path))

MM_individual_df$log_length_mm <- log(MM_individual_df$length_mm) # add log transformation of length_mm
MM_average_df$log_C <- log(MM_average_df$C) # add log transformation of Consumption 
MM_average_df$log_R <- log(MM_average_df$R) # add log transformation of Respiration  
MM_average_df$log_E <- log(MM_average_df$E) # add log transformation of Egestion  


## Filter data for each species for species based statistics 
Ao_individual_data <- MM_individual_df %>% filter(species == "A_o")
Gh_individual_data <- MM_individual_df %>% filter(species == "G_h")
Io_individual_data <- MM_individual_df %>% filter(species == "I_o")



#######################################
##						##
##   Plot: taxa size categories ~ HWs ##
##						##
#######################################

## Plot: each treatment 
## Ao_individual_data
histo_Ao <-ggplot(data = Ao_individual_data, aes(x = length_mm, fill = treatment)) +
  geom_histogram(binwidth = 2.5, boundary = 0, color = "black", alpha=1) +  #alpha=0.6
  facet_wrap(~ treatment) +
  scale_fill_manual(values = c("0HW" = "#1F5673",  "1HW" = "#E3B505", "3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10), limits =c(0,35)) +
  guides(fill="none")+
  labs(x = "\nBody length (mm)", y = "Counts\n") +
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(),strip.text = element_blank()) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        #axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25),
        axis.text.x  = element_blank(),   # remove x-axis text
        axis.ticks.x = element_blank())   
histo_Ao


## Io_individual_data
histo_Io <- ggplot(data = Io_individual_data, aes(x = length_mm, fill = treatment)) +
  geom_histogram(binwidth = 2.5, boundary = 0, color = "black",alpha=1) + # alpha=0.6
  facet_wrap(~ treatment) +
  scale_fill_manual(values = c(
    "0HW" = "#1F5673",  
    "1HW" = "#E3B505",  
    "3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10), limits =c(0,35)) +
  guides(fill="none")+
  labs(x = "\nBody length (mm)", y = "Counts\n") +
  theme_light(base_size = 24) +
  #theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.grid = element_blank(),strip.text = element_blank(), legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        #axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25),
        axis.text.x  = element_blank(),   # remove x-axis text
        axis.ticks.x = element_blank())   
histo_Io


## G_h_individual_data
histo_Gh <-ggplot(data = Gh_individual_data, aes(x = length_mm, fill = treatment)) +
  geom_histogram(binwidth = 2.5, boundary = 0, color = "black", alpha=1) + #alpha=0.6
  facet_wrap(~ treatment) +
  scale_fill_manual(values = c(
    "0HW" = "#1F5673",  
    "1HW" = "#E3B505",  
    "3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10), limits =c(0,35)) +
  guides(fill="none")+
  labs(x = "\nBody length (mm)", y = "Counts\n") +
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(),strip.text = element_blank(), legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25))
    
histo_Gh



#################################### MEAN COUNT NUMBERS + TANK INDIVIDUAL POINTS ########################

###########################################################
##						##
##   Plot: taxa size categories ~ HWs + INDIVIDUAL TANKS ##
##						##
###########################################################


binwidth <- 2.5


##
## Ao individuals
##
# tank-level counts per bin
Ao_tank_counts <- Ao_individual_data %>%
  mutate(bin = floor(length_mm / binwidth) *           # assign each individual to size bin 
           binwidth, bin_mid = bin + binwidth / 2) %>% # binmid: position where bar should be plotted on x-axis (bin middle)
  group_by(treatment, tank, bin_mid) %>%               # count individuals within each tank and bin 
  summarise(count = n(), .groups = "drop")             # each row: one tank x one bin size (values plotted as points)

# mean count across tanks (for bars)                  
Ao_mean_counts <- Ao_tank_counts %>%                  # average count across tanks for each treatment (bar height: mean count across treatment; points: indiv. tanks counts)
  group_by(treatment, bin_mid) %>%
  summarise(mean_count = mean(count),.groups = "drop")

## plot
histo_Ao <- ggplot() +
  
  # mean bars
  geom_col(data = Ao_mean_counts,
           aes(x = bin_mid, y = mean_count, fill = treatment),
           color = "black",
           width = binwidth,
           alpha = 1) +
  
  # individual tank points
  geom_point(
    data = Ao_tank_counts,
    aes(x = bin_mid, y = count),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    stroke = 0.4,
    color = "black",
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(width = 0.15, height = 0)) +
  
  facet_wrap(~treatment) +
  
  scale_fill_manual(values = c("0HW" = "#1F5673","1HW" = "#E3B505","3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10),limits = c(0, 35)) +
  
  labs(x = "\nBody length (mm)", y = "Mean counts\n") +

theme_light(base_size = 24) +
  theme(panel.grid = element_blank(),strip.text = element_blank(), legend.position = "none") +

theme(axis.title.x = element_text(size = 28),   
      axis.title.y = element_text(size = 28),
      #axis.text.x  = element_text(size = 25),   
      axis.text.y  = element_text(size = 25),
      axis.text.x  = element_blank(),   # remove x-axis text
      axis.ticks.x = element_blank()) 

histo_Ao



##
## Io individuals
##
# tank-level counts per bin
Io_tank_counts <- Io_individual_data %>%
  mutate(bin = floor(length_mm / binwidth) *           # assign each individual to size bin 
           binwidth, bin_mid = bin + binwidth / 2) %>% # binmid: position where bar should be plotted on x-axis (bin middle)
  group_by(treatment, tank, bin_mid) %>%               # count individuals within each tank and bin 
  summarise(count = n(), .groups = "drop")             # each row: one tank x one bin size (values plotted as points)

# mean count across tanks (for bars)                  
Io_mean_counts <- Io_tank_counts %>%                  # average count across tanks for each treatment (bar height: mean count across treatment; points: indiv. tanks counts)
  group_by(treatment, bin_mid) %>%
  summarise(mean_count = mean(count),.groups = "drop")

## plot
histo_Io <- ggplot() +
  
  # mean bars
  geom_col(data = Io_mean_counts,
           aes(x = bin_mid, y = mean_count, fill = treatment),
           color = "black",
           width = binwidth,
           alpha = 1) +
  
  # individual tank points
  geom_point(
    data = Io_tank_counts,
    aes(x = bin_mid, y = count),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    stroke = 0.4,
    color = "black",
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(width = 0.15, height = 0)) +
  
  facet_wrap(~treatment) +
  
  scale_fill_manual(values = c("0HW" = "#1F5673","1HW" = "#E3B505","3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10),limits = c(0, 35)) +
  
  labs(x = "\nBody length (mm)", y = "Mean counts\n") +
  
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(),strip.text = element_blank(), legend.position = "none") +
  
  theme(axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        #axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25),
        axis.text.x  = element_blank(),   # remove x-axis text
        axis.ticks.x = element_blank()) 

histo_Io




##
## Gh individuals
##
# tank-level counts per bin
Gh_tank_counts <- Gh_individual_data %>%
  mutate(bin = floor(length_mm / binwidth) *           # assign each individual to size bin 
           binwidth, bin_mid = bin + binwidth / 2) %>% # binmid: position where bar should be plotted on x-axis (bin middle)
  group_by(treatment, tank, bin_mid) %>%               # count individuals within each tank and bin 
  summarise(count = n(), .groups = "drop")             # each row: one tank x one bin size (values plotted as points)

# mean count across tanks (for bars)                  
Gh_mean_counts <- Gh_tank_counts %>%                  # average count across tanks for each treatment (bar height: mean count across treatment; points: indiv. tanks counts)
  group_by(treatment, bin_mid) %>%
  summarise(mean_count = mean(count),.groups = "drop")

## plot
histo_Gh <- ggplot() +
  
  # mean bars
  geom_col(data = Gh_mean_counts,
           aes(x = bin_mid, y = mean_count, fill = treatment),
           color = "black",
           width = binwidth,
           alpha = 1) +
  
  # individual tank points
  geom_point(
    data = Gh_tank_counts,
    aes(x = bin_mid, y = count),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    stroke = 0.4,
    color = "black",
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(width = 0.15, height = 0)) +
  
  facet_wrap(~treatment) +
  
  scale_fill_manual(values = c("0HW" = "#1F5673","1HW" = "#E3B505","3HW" = "#D95D39")) +
  scale_x_continuous(breaks = seq(0, 30, by = 10),limits = c(0, 35)) +
  labs(x = "\nBody length (mm)", y = "Mean counts\n") +
  
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(),strip.text = element_blank(), legend.position = "none") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25))

histo_Gh



