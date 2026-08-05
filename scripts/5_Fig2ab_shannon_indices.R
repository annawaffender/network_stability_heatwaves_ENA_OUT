
#### Figure 2a,b: Shannon index biomasses and flows
#### Statistical analysis and plots

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.packages)
library(ggplot2)
library(vegan)
library(dplyr)
library(MASS)

## Set paths
your_path <- "your_path/Github/network_stability_heatwaves_ENA_M_OUT/"  # replace with your path to Github folder 
output_path <- file.path(your_path,"outputs")


#######################################
##						##
## LOAD AND PREPARE DATA   ##
##						##
#######################################

treat <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")
HW_col = c("#1F5673", "#E3B505", "#D95D39")
HW_col_light <-  c("#94BAD7", "#FFE0A0", "#FFBBB1")

## Biomasses
compartments_biomasses_orig <- read.csv(file.path(output_path, "df_network_compartments.csv")) # export file for the three primary producers (mortalities)
compartments_biomasses <- compartments_biomasses_orig %>% dplyr::select(compartments, tank, treat, B)
abundances_df <- read.csv(file.path(output_path,"df_mesograzer_infauna_indiv_size.csv"))

## Flows
flow_tank_matrices  <- readRDS(file.path(output_path,"M_out_ZTER.rds")) # flow matrices with columns Z,T,E,R



#######################################
##						##
## Figure 2a:  SHANNON: BIOMASSES    ##
##						##
######################################

shannon_B <- compartments_biomasses %>%
  group_by(tank, treat) %>%
  summarise(shannon_values_B = diversity(B, index = "shannon"),.groups = "drop")
shannon_B

#write.csv(shannon_B, file=file.path(output_path,"shannon_Hb.csv"), row.names=FALSE)

#### STATS ####
## lm: Gaussian
lm_shannon_B <- lm(shannon_values_B ~ treat, data = shannon_B)  
summary(lm_shannon_B)

AIC(lm_shannon_B)
par(mfrow=c(2,2))
plot(lm_shannon_B)
par(mfrow=c(1,1))


#### PLOT ####
shannon_B_mean <- shannon_B %>% group_by(treat) %>%
  summarise(mean_shannon = mean(shannon_values_B, na.rm =TRUE),
            sd_shannon = sd(shannon_values_B, na.rm =TRUE),
            .groups = "drop")


plot_shannon_biomass <- ggplot() + 
  geom_point(data = shannon_B_mean, aes(x = treat, y = mean_shannon, color= treat), 
             size = 6, alpha=1) +
  geom_errorbar(data = shannon_B_mean, aes(x = treat, ymin = mean_shannon - sd_shannon, ymax = mean_shannon + sd_shannon, color = treat), 
                width = 0.2, size=1, alpha=1) +
  geom_jitter(data = shannon_B, 
              aes(x = treat, y = shannon_values_B, fill = treat), 
              shape = 21, color = "black", width = 0.15, size = 4.25, alpha = 1) +
  scale_fill_manual(values = HW_col_light) +  # fill for jitter points 
  scale_color_manual(values = HW_col) +
  scale_x_discrete(labels = c("0HW" = "0","1HW" = "1","3HW" = "3")) +
  labs(x = "\nHeatwaves",y = expression(paste("Biodiversity (", italic(H[B]), ")"), "\n"))+
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  #theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25))  
plot_shannon_biomass


plot_shannon_biomass_narrow <- plot_shannon_biomass +
  theme(aspect.ratio = 8/6)   # height / desired panel width
plot_shannon_biomass_narrow




######################################
##						##
## Figure 2b:  SHANNON: FLOWS       ##
##						##
######################################

## Apply shannon index calculation to all matrices
M <- flow_tank_matrices
calc_shannon <- function(M) {
  M <- M[, 1:13]         # ensure it's a matrix
  p <- M / sum(M)           # normalize to relative flows
  p <- p[p > 0]             # remove zeros
  shannon <- -sum(p * log(p))    # Shannon entropy in bits
  return(shannon)
}

shannon_vector <- sapply(flow_tank_matrices, calc_shannon)
shannon_vector


## Convert vector to dataframe 
shannon_flow_df <- data.frame(tank = names(shannon_vector),   # names of the vector
                         shannon_values_flows = as.numeric(shannon_vector))  # the numeric values

shannon_flow_df$treat <- treat
shannon_flow_df

#write.csv(shannon_flow_df, file=file.path(output_path,"shannon_Hf.csv"), row.names=FALSE)



#### STATS ####
## lm: Gaussian
lm_shannon_flow <- glm(shannon_values_flows ~ treat, family= gaussian, data = shannon_flow_df)  
summary(lm_shannon_flow)

AIC(lm_shannon_flow)
par(mfrow=c(2,2))
plot(lm_shannon_flow)
par(mfrow=c(1,1))



#### PLOTS ####
shannon_flow_mean <- shannon_flow_df %>% group_by(treat) %>%
  summarise(mean_shannon = mean(shannon_values_flows, na.rm =TRUE),
            sd_shannon = sd(shannon_values_flows, na.rm =TRUE),
            .groups = "drop")

plot_shannon_flow <- ggplot() + 
  geom_point(data = shannon_flow_mean, aes(x = treat, y = mean_shannon, color= treat), size = 5) +
  geom_errorbar(data = shannon_flow_mean, 
                aes(x = treat, ymin = mean_shannon - sd_shannon, ymax = mean_shannon + sd_shannon, color = treat), 
                width = 0.2, size=1) +
  geom_jitter(data = shannon_flow_df, 
              aes(x = treat, y = shannon_values_flows, fill = treat), 
              shape = 21, color = "black", width = 0.15, size = 4, alpha = 1) +
  scale_fill_manual(values = HW_col_light) +  # fill for jitter points
  scale_color_manual(values = HW_col) +
  scale_x_discrete(labels = c("0HW" = "0","1HW" = "1","3HW" = "3")) +
  #labs(x = "\nHeatwaves", y = expression(paste("Flow diversity (", italic(F), ")"), "\n"))+
  labs(x = "\nHeatwaves",y = expression(paste("Flow diversity (", italic(H[F]), ")"), "\n"))+
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25))

plot_shannon_flow


plot_shannon_flow_narrow <- plot_shannon_flow +
  theme(aspect.ratio = 8/6)   # height / desired panel width
plot_shannon_flow_narrow

