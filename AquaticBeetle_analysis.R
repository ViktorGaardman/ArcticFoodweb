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
                   "female", "Adult_female", "Adult_male",
                   "Unknown") ~ "Adult")
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
  doy = seq(150, 250, by = 1)
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
  group_by(doy, year) %>%
  summarise(
    small = sum(Size.Gender %in% c("Small", "Medium", "medium")),
    large = sum(Size.Gender %in% c("Large", "XL")),
    adult = sum(Size.Gender %in% c("female", "male", "Female",
                                   "Male", "Adult_female", 
                                   "Adult_male", "Unknown")),
    total = n(),
    .groups = "drop"
  ) %>%
  mutate(
    small = small / total,
    large = large / total,
    adult = adult / total
  ) |> 
  mutate(
    sum = small + large + adult
  )

size_long <- size_data %>%
  pivot_longer(
    cols = c(small, large, adult),
    names_to = "size_class",
    values_to = "proportion"
  ) 

size_long$year <- as.factor(size_long$year)

#Offset doy 183 2024 to not overlap with 2025
size_long$doy[size_long$year == 2024 & size_long$doy == 183] <- 182

#Make a line for day 182 and day 212 (start and end of July)
larvaldev_plot_year <- ggplot() +
  # Observed size-class proportions
  geom_col(
    data = size_long,
    aes(
      x = doy,
      y = proportion,
      fill = size_class,
      alpha = year
    ),
    width = 1
  ) +
  geom_vline(
    xintercept = 153,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 183,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 213,
    linetype = "dashed",
    colour = "black",
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
    breaks = seq(0, 1, 0.25),
    expand = c(0,0),
    name = "Proportion of larvae"
  ) +
  scale_fill_manual(
    values = c(
      small = "#E69F00",
      large = "#009E73",
      adult = "#CC79A7"
    ),
    breaks = c("small", "large", "adult"),
    name = "Development stage"
  ) +
  scale_alpha_manual(
    name = "Year",
    values = c(
      "2024" = 1,
      "2025" = 0.5
    )
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

ggsave(larvaldev_plot_year, filename = "aquaticpred_devplot_year.png",
       dpi = 300, width = 8, height = 5)

#2. 
#Check if diet differs between larval stages and between larvae and adults
#Get diet data
#Will require re-cleaning of specifically the C.dolabratus data
#We can work directly with the edgelists to figure this out
#some type of chi-sqr test for the distribution of FOO or meanRRA

matrix24 <- read_xlsx("mergedData_wide_ReadCounts_2024.xlsx")
matrix25 <- read_xlsx("mergedData_wide_ReadCounts_2025.xlsx")

matrix24 <- matrix24 |> 
  select(-Kingdom)

#combine
long24 <- pivot_longer(matrix24,
                       -c(Phylum:Assigned),
                       names_to = "Sample",
                       values_to = "Reads")

long25 <- pivot_longer(matrix25,
                       -c(Phylum:Assigned),
                       names_to = "Sample",
                       values_to = "Reads")

combined_df <- bind_rows(long24, long25)


#Remove the plants and worms for now
combined_df <- combined_df |> 
  filter(!Phylum %in% c("Streptophyta", "Annelida"))

#Remove prey species that don't exist in Greenland
combined_df <- combined_df |> 
  filter(!Species %in% c("Hypera_diversipunctata", "Otiorhynchus_singularis",
                         "Phyllobius_pyri", "Strophosoma_capitatum"))

#Remove prey families that don't exist in Greenland
combined_df <- combined_df %>%
  filter(!Family %in% c("Blattidae", "Formicidae", "Elateridae", "Cantharidae",
                        "Lasiocampidae", "Elmidae", "Xylophagidae", "Pentatomidae",
                        "Rhagionidae", "Clubionidae", "Anyphaenidae",
                        "Chaoboridae", "Clusiidae", "Dolomedidae", 
                        "Ephemerellidae", "Gryllidae", "Leuctridae",
                        "NA", "Nemouridae", "Ptychopteridae",
                        "Salticidae", "Sarcophagidae"))

#Keep only the beetle samples
#1. Find beetle samples
beetle_foodweb <- combined_df %>%
  filter(Sample %in% df_tot$ID)

beetle_foodweb <- beetle_foodweb |> 
  left_join(
    df_tot %>%
      select(ID, Stage, Size.Gender),
    by = c("Sample" = "ID")
  )

beetle_foodweb <- beetle_foodweb |> 
  mutate(Size = case_when(
    Size.Gender %in% c("Medium", "medium", "Small") ~ "Small",
  Size.Gender %in% c("XL", "Large") ~ "Large",
  Size.Gender %in% c("Male", "male", "Female",
                     "female", "Adult_female", "Adult_male",
                     "Unknown") ~ "Adult")
)

family_pa <- beetle_foodweb %>%
  filter(!Family %in% "Dytiscidae") |> 
  group_by(Sample, Size, Family) |> 
  summarise(
    Present = as.integer(sum(Reads, na.rm = TRUE) > 0),
    .groups = "drop"
  ) %>%
  group_by(Family) |> 
  filter(sum(Present, na.rm = TRUE) > 0) %>%
  ungroup() |> 
  group_by(Sample) |> 
  filter(sum(Present, na.rm = TRUE) > 0) |> 
  ungroup()

#Compare size classes, also compare adult males-females if there is
#enough data
#First group prey by family level and samples by beetle Size
#Then remove all interactions that are always 0 and self-reads
dyt_samples <- beetle_foodweb %>%
  filter(Family == "Dytiscidae", Reads > 0) %>%
  pull(Sample) %>%
  unique()

beetle_foodweb <- beetle_foodweb %>%
  filter(Sample %in% dyt_samples)

family_foodweb <- beetle_foodweb %>%
  group_by(Size, Family) %>%
  summarise(
    N_present = n_distinct(Sample[Reads > 0]),
    .groups = "drop"
  ) %>%
  left_join(
    beetle_foodweb %>%
      group_by(Size) %>%
      summarise(
        N_samples = n_distinct(Sample),
        .groups = "drop"
      ),
    by = "Size"
    ) %>%
  mutate(
    FOO = N_present / N_samples
  )

family_foodweb <- family_foodweb %>%
  group_by(Family) |> 
  filter(sum(FOO, na.rm = TRUE) > 0) %>%
  ungroup() |> 
  filter(!Family %in% "Dytiscidae")

beetle_matrix_FOO <- family_foodweb %>%
  select(Size, Family, FOO) %>%
  pivot_wider(
    id_cols = Size,
    names_from = Family,
    values_from = FOO,
    values_fill = 0
  )%>%
  column_to_rownames("Size")
  

#Convert to matrix

png(
  "beetle_visweb.png",
  width = 1200,
  height = 500,
  res = 300
)

visweb(
  as.matrix(beetle_matrix_FOO),
  box.col = "cornflowerblue",
  square = "b",
  labsize = 0.6
)

dev.off()
