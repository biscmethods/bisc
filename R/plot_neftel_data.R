R_path <- here::here("R")
path_temp <- here::here('temp')
path_data <- here::here('data')
library(tidyverse)

group2 <- readRDS(  file.path(path_data,"neftel_mn_group2.rds"))

group2_cells <- as_tibble(read.table(file.path(path_data,"Neftel2019", "Group2","cells.txt"),
                           sep = " ", header = TRUE))
# View(group2_cells)

cell_names1 <- rownames(group2)

cell_names2 <- group2_cells$cell_name

str(cell_names1)
str(cell_names2)

colnames(group2)
names(group2_cells)




# Yes, the relative meta-module score mentioned in Figure 3F can be derived
# from the variables in your dataset. The equation (\log_2(|SC1 - SC2| + 1))
# suggests that you need to calculate the logarithm of the absolute difference
# between two meta-module scores plus one. In your dataset, the meta-module
# scores could correspond to the variables like
# “MESlike2”, “MESlike1”, “AClike”, “OPClike”, “NPClike1”, and “NPClike2”

rmms <- function(a,b){
  log2(abs(a - b) + 1)
}
group2_cells %>% mutate(
rmms1_x = rmms(NPClike2, OPClike),
rmms1_y = rmms(NPClike2, MESlike2),

rmms2_x = rmms1_x,
rmms2_y = rmms(OPClike, AClike),

rmms3_x = rmms(AClike,MESlike2),
rmms3_y = rmms2_y ,

rmms4_x = rmms3_x,
rmms4_y = rmms1_y,
) -> dat

dat_list <- list(
  dat2 = dat %>% select( rmms_x = rmms1_x, rmms_y = rmms1_y),
  dat1 = dat %>% select( rmms_x = rmms2_x, rmms_y = rmms2_y),
  dat3 = dat %>% select( rmms_x = rmms3_x, rmms_y = rmms3_y),
  dat4 = dat %>% select( rmms_x = rmms4_x, rmms_y = rmms4_y)
)

# Bind all data frames together
dat_long <- bind_rows(dat_list, .id = "id")

# Create the plot
ggplot(dat_long, aes(x = rmms_x, y = rmms_y)) +
  geom_point() +
  facet_wrap(~id, ncol = 2) +
  theme_bw() +
  labs(x = "rmms x", y = "rmms y", title = "2x2 Grid of Plots")


dat %>%
  ggplot( aes(x = rmms1, y = rmms2, color = cycling_score)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_point() + # This adds the scatterplot points
  theme_minimal()



# first, try to recreate plot 3F

group2_cells %>%
  ggplot( aes(x = tSNE1, y = tSNE2)) +
  geom_point() + # This adds the scatterplot points
  theme_minimal() + # Optional: This gives a minimalistic theme to the plot
  labs(title = "tSNE Scatterplot", x = "tSNE1", y = "tSNE2")

group2_cells %>%
  ggplot( aes(x =MESlike1, y = MESlike2)) +
  geom_point() + # This adds the scatterplot points
  theme_minimal() # Optional: This gives a minimalistic theme to the plot


group2_cells %>%
  filter(malignant == "yes") %>%
  ggplot( aes(x = tSNE1, y = tSNE2, color = cell_name)) +
  geom_point() + # This adds the scatterplot points
  theme_minimal() + # Optional: This gives a minimalistic theme to the plot
  labs(title = "tSNE Scatterplot", x = "tSNE1", y = "tSNE2")


names(group2_cells)


