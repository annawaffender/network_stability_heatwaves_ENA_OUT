
#### Figure 4: Summary of the (a) carbon flows, (b) interaction strength, (c) impact on stability for each resource - consumer feeding interaction
#### Analyse (a,b,c) separately for the summarized feeding categories: detritivory (POC); herbivory (ZM, FV, FA); Polychaeta omnivorous (Po; see supporting information)

rm(list = ls(all = TRUE))  # clear workspace

## load needed packages (eventually install.packages)
library(lmerTest)
library(lme4)
library(DHARMa)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

## Set paths
your_path <- "your_path/Github"  # replace with your path to Github folder 
output_path <- file.path(your_path,"/Github/network_stability_heatwaves_ENA_M_OUT/outputs")



#######################################
##						##
## LOAD AND PREPARE DATA   ##
##						##
#######################################

treatment <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")
tank_treat <- c(A1 = "0HW", A2 = "0HW", D1 = "0HW", D2 = "0HW",B2 = "3HW", E1 = "3HW", E2 = "3HW",C1 = "1HW", C2 = "1HW", F1 = "1HW", F2 = "1HW")

flow_matrices <- readRDS(file.path(output_path,"M_out_T.rds"))
interaction_matrices <- readRDS(file.path(output_path,"list_interactions_matrices.rds"))
pair_stability <- read.csv(file.path(output_path,"stability_sensitivity_interactions.csv"))

## convert interaction_matrices to long format
df_interactions <- imap_dfr(interaction_matrices, function(mat, tank_name) {
  as.data.frame(mat) %>%
    rownames_to_column("resource") %>%
    pivot_longer(cols = -resource,
                 names_to = "consumer",
                 values_to = "interactions") %>%
    filter(interactions != 0, resource != consumer) %>%
    mutate(tank = tank_name,
           treat = tank_treat[tank_name])
})

## add treatment to pair_stability 
pair_stability <- pair_stability %>% mutate(treat = tank_treat[tank])




#######################################
##						##
## Figure 4a:  CARBON FLOWS        ##
##						##
######################################

## divide carbon flows into the flow categories of the following feeding categories: detritivory (POC); herbivory (ZM, FV, FA); Polychaeta omnivorous (Po)
## analyse the flows separately for separate feeding categories  
df_POC_flows <- data.frame(tank = names(flow_matrices), POC_flows_sum = sapply(flow_matrices, function(m) {sum(m[nrow(m), 1:12])}))
df_herbivory_flows <- data.frame(tank = names(flow_matrices), herbiv_flows_sum = sapply(flow_matrices, function(m) {sum(m[1:3, 1:12])}))
df_Po_flows <- data.frame(tank = names(flow_matrices), Po_flows_sum = sapply(flow_matrices, function(m) {sum(m[c(4, 5, 8:12),12 ])}))

## combine the data frames by tank
df_combined_1 <- merge(df_POC_flows,df_herbivory_flows, by = "tank")
df_combined <- merge(df_combined_1,df_Po_flows, by="tank")
df_combined$ratio_detriv_herbiv <- df_combined$POC_flows_sum / df_combined$herbiv_flows_sum
df_combined$treat <- treatment 
df_combined


## DETRITIVORY ------------------------------------------------------
POC_lm <- lm(POC_flows_sum~treat, data=df_combined)
summary(POC_lm)

## check model fit
par(mfrow=c(2,2))
plot(POC_lm)
par(mfrow=c(1,1))
AIC(POC_lm)


## HERBIVORY ------------------------------------------------------
herbiv_lm <- lm(herbiv_flows_sum~treat, data=df_combined)
summary(herbiv_lm)

## check model fit
par(mfrow=c(2,2))
plot(herbiv_lm)
par(mfrow=c(1,1))
AIC(herbiv_lm)


## DETRITIVORY / HERBIVORY ------------------------------------------------------
detriv_herbiv_lm <- glm(ratio_detriv_herbiv~treat,family= Gamma,data=df_combined)
summary(detriv_herbiv_lm)

## check model fit
par(mfrow=c(2,2))
plot(detriv_herbiv_lm)
par(mfrow=c(1,1))
AIC(detriv_herbiv_lm)


## POLYCHAETA ONMIVOROUS (see supporting information Fig. S2) ------------------
Po_glm <- glm(Po_flows_sum~treat, family=Gamma, data=df_combined)
summary(Po_glm)

## check model fit
par(mfrow=c(2,2))
plot(Po_glm)
par(mfrow=c(1,1))
AIC(Po_glm)



#######################################
##						##
## Figure 4b:  INTERACTION STRENGTH ##
##						##
######################################

## positive interactions: BOTTOM-UP control 
## negative interaction: TOP-DOWN control


## DETRITIVORY ------------------------------------------------------
POC_pos_interact <- df_interactions[df_interactions$resource == "POC",] 
POC_neg_interact <- df_interactions[df_interactions$consumer == "POC",] 

## POSITVE INTERACTIONS: BOTTOM-UP control of POC
## GLM GAMMA: with inverse link reverses the intuitive sign of coefficients
#glm_pos_POC <- glmmTMB(interactions~treat + (1|tank), family=Gamma, data=POC_pos_interact)
#summary(glm_pos_POC)

# exclude random effect, as variance near zero 
glm_pos_POC_no_re <- glmmTMB(interactions ~ treat,family = Gamma(link = "log"),data = POC_pos_interact) 
summary(glm_pos_POC_no_re)


## check model fit
sim_res_nb <- simulateResiduals(glm_pos_POC_no_re)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glm_pos_POC_no_re)
fitted_vals <- fitted(glm_pos_POC_no_re)
plot(fitted_vals, resid_vals,xlab = "Fitted values",ylab = "Residuals",main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)


## NEGATIVE INTERACTIONS: TOP-DOWN control on POC  
# GLM GAMMA: with inverse link reverses the intuitive sign of coefficients.
glm_neg_POC <- glmmTMB((interactions*-1)~treat + (1|tank), family=Gamma, data=POC_neg_interact)
summary(glm_neg_POC)

## check model fit
sim_res_nb <- simulateResiduals(glm_neg_POC)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glm_neg_POC)
fitted_vals <- fitted(glm_neg_POC)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)



## HERBIVORY ------------------------------------------------------
Herbiv_pos_interact <- df_interactions[df_interactions$resource %in% c("ZM", "FA", "FV"), ]
Herbiv_neg_interact <- df_interactions[df_interactions$consumer %in% c("ZM", "FA", "FV"), ]

## POSITIVE INTERACTIONS: BOTTOM-UP control of ZM,FV,FA
## glmmTMB with log link: matches skewed data
glmm_pos_Herbiv <- glmmTMB(log(interactions) ~ treat + (1|tank), data = Herbiv_pos_interact)
summary(glmm_pos_Herbiv)

## check model fit
sim_res_nb <- simulateResiduals(glmm_pos_Herbiv)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glmm_pos_Herbiv)
fitted_vals <- fitted(glmm_pos_Herbiv)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)


## NEGATIVE INTERACTIONS: TOP-DOWN control on ZM,FV,FA
glm_neg_Herbiv <- glmmTMB((interactions*-1) ~treat + (1|tank), family=Gamma(link="log"), data=Herbiv_neg_interact)
summary(glm_neg_Herbiv)

## check model fit
sim_res_nb <- simulateResiduals(glm_neg_Herbiv)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glm_neg_Herbiv)
fitted_vals <- fitted(glm_neg_Herbiv)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)




## POLYCHAETA ONMIVOROUS (see supporting information Fig. S2) ------------------
Po_pos_interact <- df_interactions[df_interactions$consumer == "P_o" & df_interactions$resource != "POC" & df_interactions$resource != "P_o", ]
Po_neg_interact <- df_interactions[df_interactions$resource == "P_o"& df_interactions$consumer != "POC", ]


## POSITIVE INTERACTIONS: BOTTOM-UP control of Polychaeta omnivorous (Po)
## use mean values per tank, since within tank values are all the same 
Po_pos_interact_means <- Po_pos_interact %>% group_by(tank, treat) %>% dplyr::summarise(interactions_mean = mean(interactions),.groups = "drop")
lm_Po_pos_mean <- lm(interactions_mean~treat, data= Po_pos_interact_means)
summary(lm_Po_pos_mean)

## check model fit
sim_res_nb <- simulateResiduals(lm_Po_pos_mean)
plotQQunif(sim_res_nb)
resid_vals <- residuals(lm_Po_pos_mean)
fitted_vals <- fitted(lm_Po_pos_mean)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)


## NEGATIVE INTERACTIONS: TOP-DOWN control on Polychaeta omnivorous (Po)
## use mean values per tank, since within tank values are all the same 
Po_neg_interact_means <- Po_neg_interact %>%group_by(tank, treat) %>% dplyr::summarise(interactions_mean = mean(interactions), .groups = "drop")
lm_Po_neg_mean <- lm(interactions_mean~treat, data= Po_neg_interact_means)
summary(lm_Po_neg_mean)

## check model fit 
sim_res_nb <- simulateResiduals(lm_Po_neg_mean)
plotQQunif(sim_res_nb)
resid_vals <- residuals(lm_Po_neg_mean)
fitted_vals <- fitted(lm_Po_neg_mean)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)





###########################################
##						##
## Figure 4c:  INTERACTIONS ON STABILITY ##
##						##
###########################################

## prepare data
POC_pair_stability <- pair_stability[pair_stability$row=="POC",]
Herbiv_pair_stability <- pair_stability[pair_stability$row %in% c("ZM","FA","FV"), ]
Po_pair_stability <- pair_stability[pair_stability$column == "P_o" & pair_stability$row != "POC", ]


## DETRITIVORY ------------------------------------------------------
## glmmTMB: negative binomial distribution (no poisson since overdispersion) with random effect 
## truncated_binom1: since there are no zeros in the data set (otherwise model will expect zeros and misfit model)
glmm_POC_stab <- glmmTMB(stable ~ treat + (1|tank), family = truncated_nbinom1, data = POC_pair_stability)
summary(glmm_POC_stab)


## check model fit
sim_res_nb <- simulateResiduals(glmm_POC_stab)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glmm_POC_stab)
fitted_vals <- fitted(glmm_POC_stab)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)





## HERBIVORY ------------------------------------------------------
## glmmTMB: negative binomial distribution (no poisson since overdispersion) with random effect 
glmm_Herbiv_stab <- glmmTMB(stable ~ treat + (1|tank),family = truncated_nbinom1, data = Herbiv_pair_stability)
summary(glmm_Herbiv_stab)

## check model fit
sim_res_nb <- simulateResiduals(glmm_Herbiv_stab)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glmm_Herbiv_stab)
fitted_vals <- fitted(glmm_Herbiv_stab)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)




## POLYCHAETA ONMIVOROUS (see supporting information Fig. S2) ------------------
## glmmTMB: negative binomial distribution (no poisson since overdispersion) with random effect 
glmm_Po_stab <- glmmTMB(stable ~ treat + (1|tank),family = nbinom1,data = Po_pair_stability)
summary(glmm_Po_stab)

## check model fit 
sim_res_nb <- simulateResiduals(glmm_Po_stab)
plotQQunif(sim_res_nb)
resid_vals <- residuals(glmm_Po_stab)
fitted_vals <- fitted(glmm_Po_stab)
plot(fitted_vals, resid_vals,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

