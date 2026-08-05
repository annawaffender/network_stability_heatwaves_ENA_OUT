
#### Figure 2c: Local stability 
#### Statistical analysis and plots

rm(list=ls(all=TRUE))	## clear workspace

## load needed packages (eventually install.packages)
library(tidyverse)


## Set paths
your_path <- "your_path/Github/network_stability_heatwaves_ENA_M_OUT/"  # replace with your path to Github folder 
output_path <- file.path(your_path,"outputs")


#######################################
##						##
## LOAD AND PREPARE DATA   ##
##						##
#######################################

treatments <- c("0HW", "0HW", "3HW", "1HW", "1HW", "0HW", "0HW", "3HW", "3HW", "1HW", "1HW")

matrix_PoPo_interact <- read.csv(file.path(output_path,"stability_matrix.csv"), header= TRUE,row.names=1)
shift_points <- as.numeric(row.names(matrix_PoPo_interact))  # Convert row names (shift points) to numeric
treatments_matrix <- matrix_PoPo_interact  # Extract the treatment matrix (rest of the data excluding the row names)
column_names <- colnames(treatments_matrix)  # Define the treatment labels for each column
names(treatments) <- column_names  # Associate treatments with column names

get_shift_index <- function(column) {
  which(column == 0)[1]  #Find the row index where each treatment shifts from 1 to 0
}

shift_indexes <- sapply(treatments_matrix, get_shift_index) # Apply the function to each treatment column

# Create a dataframe for treatments and shift points
shift_df_PoPo_interact <- data.frame(Treatment = treatments[names(shift_indexes)],
                                     ShiftCounts = shift_indexes,
                                     ShiftPoint = shift_points[shift_indexes])


## save file for further analysis 
shift_df_PoPo_interact$Tank <- rownames(shift_df_PoPo_interact)
rownames(shift_df_PoPo_interact) <- NULL
#write.csv(shift_df_PoPo_interact, file=file.path(output_path,"stability_shift_counts.csv"), row.names=FALSE)



#######################################
##						##
## STATISTICS: POISSON DISTRIBUTION  ##
##						##
#######################################

shift_df_PoPo_interact$Treatment <- as.factor(shift_df_PoPo_interact$Treatment)

## POISSON: use glm with poisson distribution since ShiftCounts is count data
glm_model_PoPo_interact <- glm(ShiftCounts ~ Treatment,family = poisson,data = shift_df_PoPo_interact)
summary(glm_model_PoPo_interact)

par(mfrow = c(2, 2))
plot(glm_model_PoPo_interact)
par(mfrow = c(1, 1))




#######################################
##						##
## PLOT: LOCAL STABILITY ##
##						##
#######################################

data <- matrix_PoPo_interact 

## Convert row names (mortality values) into a column
data_long <- data %>% rownames_to_column(var = "Mortality") %>% mutate(Mortality = as.numeric(Mortality)) %>%
  pivot_longer(cols = -Mortality, names_to = "tank", values_to = "Stability")

data_long <- data_long %>%
  mutate(Treatment = case_when(tank %in% c("A1","A2","D1","D2") ~ "0HW",
                               tank %in% c("C1","C2","F1","F2") ~ "1HW",
                               tank %in% c("B2","E1","E2") ~ "3HW"))

df_tank <- data_long %>% group_by(Treatment, tank) %>%
  summarise(stable = mean(Stability == 1)*100,   # from fraction to percentage
            unstable = mean(Stability == 0)*100, # from fraction to percentage
            .groups = "drop")

df_treat <- df_tank %>% group_by(Treatment) %>%
  summarise(mean_stable_treatment = mean(stable),
            sd_stable_treatment = sd(stable),
            mean_unstable_treatment = mean(unstable),
            sd_unstable_treatment = sd(unstable),
            .groups = "drop")

df_summary <- df_tank %>% left_join(df_treat, by = "Treatment")


## Colour palette
treatment_colors <- data.frame(Treatment = c("0HW", "1HW", "3HW"),
                               stable_col = c("#1F5673", "#E3B505", "#D95D39"),   # darker for stable
                               unstable_col = c("#94BAD7", "#FFE0A0", "#FFBBB1"))  # lighter for unstable

df_summary <- df_summary %>% left_join(treatment_colors, by = "Treatment") # Merge colors into summary data frame

HW_col = c("#1F5673", "#E3B505", "#D95D39")
HW_col_light <-  c("#94BAD7", "#FFE0A0", "#FFBBB1")

## ggplot
df_tank_mean_stab <- df_tank %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    mean_stab   = mean(stable, na.rm = TRUE),
    mean_unstab = mean(unstable, na.rm = TRUE),
    sd_stab     = sd(stable, na.rm = TRUE),
    sd_unstab   = sd(unstable, na.rm = TRUE),
    .groups = "drop")

plot_stability <- ggplot() +
  geom_point(data = df_tank_mean_stab, aes(x = Treatment, y = mean_stab, color= Treatment), size = 5) +
  geom_errorbar(data = df_tank_mean_stab, 
                aes(x = Treatment, ymin = mean_stab - sd_stab, ymax = mean_stab + sd_stab, color = Treatment), 
                width = 0.2, size=1) +
  geom_jitter(data = df_tank, 
              aes(x = Treatment, y = stable, fill = Treatment), 
              shape = 21, color = "black", width = 0.15, size = 4.25, alpha = 1) +
  
  scale_fill_manual(values = HW_col_light) +  
  scale_color_manual(values = HW_col) +
  scale_x_discrete(labels = c("0HW" = "0","1HW" = "1","3HW" = "3")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = '\nHeatwaves', y = 'Ecosystem stability (%)\n') +
  theme_light(base_size = 24) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  #theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(axis.title.x = element_text(size = 28),   
        axis.title.y = element_text(size = 28),
        axis.text.x  = element_text(size = 25),   
        axis.text.y  = element_text(size = 25)) 
plot_stability  

plot_stability_narrow <- plot_stability +
  theme(aspect.ratio = 8/6)   # height / desired panel width
plot_stability_narrow


