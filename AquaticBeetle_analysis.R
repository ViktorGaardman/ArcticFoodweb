library(readxl)
library(readr)
library(tidyverse)
library(bipartite)
library(patchwork)
library(forcats)
library(igraph)
library(DHARMa)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

#1.Find the distribution of individuals across time

df_24 <- read.csv("Foodweb_samples_2024.csv", sep = ";") |> 
  filter(Group == "Aquatic_predators") |> 
  filter(!Species == "Chironomid_larvae") |> 
  mutate(Stage = case_when(
    Species %in% "C._dolabratus_larvae" ~ "Larvae",
    Species %in% c("C._dolabratus", "C._Dolabratus",
                   "Colymbetes_dolabratus") ~ "Adult")
  )


df_25 <- read.csv("FoodwebSamples_2025.csv", sep = ";")|> 
  filter(Group == "Aquatic_predators") |> 
  filter(!Species == "Gyrinus_opacus") |> 
  mutate(Stage = case_when(
    Size.Gender %in% c("XL", "Large",
                       "Medium", "Small") ~ "Larvae",
    Size.Gender %in% c("Male", "male", "Female",
                   "female", "Adult_female", "Adult_male") ~ "Adult")
  )

df_tot <- rbind(df_24, df_25)
#Make sure all dates have six digits
df_tot$Date <- str_pad(as.character(df_tot$Date),
                       width = 6,
                       side = "left",
                       pad = "0")

# Convert to Date format (DDMMYY)
df_tot$Date <- as.Date(df_tot$Date, format = "%d%m%y")

# Extract day of year (0-365) and year
df_tot$doy <- as.integer(format(df_tot$Date, "%j")) - 1
df_tot$year <- as.integer(format(df_tot$Date, "%Y"))

df_tot <- df_tot |> 
  filter(!Location %in% c("Near_pond_12", "Pond_16",
                          "Pond_20", "Pond_25", "Pond_7"))

daily_data <- df_tot %>%
  filter(!is.na(Stage)) %>%
  group_by(Location, year, doy) %>%
  summarise(
    adults = sum(Stage == "Adult"),
    larvae = sum(Stage == "Larvae"),
    total = adults + larvae,
    .groups = "drop"
  )

mod <- glmmTMB::glmmTMB(
  cbind(adults, larvae) ~ doy + I(doy^2),
  data = daily_data,
  family = glmmTMB::betabinomial()
)

mod2 <- glmmTMB::glmmTMB(
  cbind(adults, larvae) ~
    doy,
  data = daily_data,
  family = glmmTMB::betabinomial()
)

AIC(mod, mod2)

res <- simulateResiduals(mod)
plot(res)
testDispersion(res)

# Create prediction dataset
pred_data <- expand.grid(
  doy = seq(150, 225, by = 1)
)


# Predictions
pred <- predict(
  mod,
  newdata = pred_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

# Add predictions and confidence intervals
pred_data$adult <- plogis(pred$fit)

pred_data$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
pred_data$upper <- plogis(pred$fit + 1.96 * pred$se.fit)

# Larval predictions
pred_data$larvae <- 1 - pred_data$adult

plot_predictions <- pred_data %>%
  mutate(
    larvae_lower = 1 - upper,
    larvae_upper = 1 - lower
  ) %>%
  pivot_longer(
    cols = c(adult, larvae),
    names_to = "life_stage",
    values_to = "proportion"
  ) %>%
  mutate(
    lower_stage = ifelse(
      life_stage == "adult", lower, larvae_lower
    ),
    upper_stage = ifelse(
      life_stage == "adult", upper, larvae_upper
    )
  )

#Calculate size proportions
size_data <- df_tot %>%
  filter(
    Stage == "Larvae",
    !is.na(Size.Gender)
  ) %>%
  group_by(doy) %>%
  summarise(
    small = sum(Size.Gender %in% c("Small", "Medium", "medium")),
    large = sum(Size.Gender %in% c("Large", "XL")),
    total = n(),
    .groups = "drop"
  ) %>%
  mutate(
    small = small / total,
    large = large / total
  ) |> 
  mutate(
    sum = small + large
  )

size_long <- size_data %>%
  pivot_longer(
    cols = c(small, large),
    names_to = "size_class",
    values_to = "proportion"
  )

#Make a line for day 182 and day 212 (start and end of July)
larvaldev_plot <- ggplot() +
  # Observed size-class proportions
  geom_col(
    data = size_long,
    aes(
      x = doy,
      y = proportion,
      fill = size_class
    ),
    width = 1,
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = 153,
    linetype = "dashed",
    colour = "#D55E00",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 183,
    linetype = "dashed",
    colour = "#D55E00",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 213,
    linetype = "dashed",
    colour = "#D55E00",
    linewidth = 0.8
  ) +
  # Model confidence interval
  geom_ribbon(
    data = pred_data,
    aes(
      x = doy,
      ymin = 1 - upper,
      ymax = 1 - lower
    ),
    fill = "black",
    alpha = 0.15
  ) +
  # Modelled larval proportion
  geom_line(
    data = pred_data,
    aes(
      x = doy,
      y = larvae
    ),
    colour = "black",
    linewidth = 1.2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    name = "Proportion of larvae"
  ) +
  scale_fill_manual(
    values = c(
      small = "#E69F00",
      large = "#009E73"
    ),
    breaks = c("small", "large"),
    name = "Larval size class"
  ) +
  labs(
    x = "Day of year"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(larvaldev_plot, filename = "aquaticpred_devplot.png",
       dpi = 300, width = 8, height = 5)

#2. 
#Check if diet differs between larval stages and between larvae and adults
#Get diet data
#Will require re-cleaning of specifically the C.dolabratus data
#We can work directly with the edgelists to figure this out
#some type of chi-sqr test for the distribution of FOO or meanRRA
