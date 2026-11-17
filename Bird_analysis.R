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

##1. Calculate relative abundance and relative biomass of flying prey

prey24 <- read.csv("MalaiseSamples_2024.csv", sep = ";")
prey25 <- read.csv("MalaiseSamples_2025.csv", sep = ";")

prey_raw <- rbind(prey24, prey25)

prey_raw <- prey_raw |> 
  mutate(
    Year = substr(as.character(Date_Start), 1, 4))

prey_raw <- prey_raw |> 
  mutate(Month = case_when(
    Week %in% 23:26 ~ "June",
    Week %in% 27:30 ~ "July",
    Week %in% 31:35 ~ "August"))

prey_sum <- prey_raw |>     
  group_by(Year, Month, Order, Taxon) |> 
  summarise(
    TotalCount = sum(Count)
  ) |> 
  ungroup()

#Remove ground-dwelling groups and microlepidoptera
prey_sum <- prey_sum |> 
  filter(!Order %in% c("Arenea", "Opiliones", "Coleoptera")) |> 
  filter(!Taxon %in% c("Microlepidoptera", "Lygaeidae"))

#Add weight data
weight_df <- read.csv("Insect_dryweight.csv", sep = ";")

mean_biomass <- weight_df |> 
  group_by(Taxon) |> 
  summarize(
    mean_weight = mean(Weight.mg.)
  ) |> 
  select(Taxon, mean_weight) |> 
  ungroup()


#2. Create full viswebs
#ADD KNOWN SPECIES TO THE MIX FROM BOTH MALAISE AND PITFALL
#TO SHOW MISSING LINKS
#IT SEEMS SOME PITFALL 2025 DATA IS NOW AVAILABLE TOO

files_RRA <- list.files(
  pattern = "^edgelist_.*_20.*_Fam_RRA\\.csv$"
)

read_edgelist <- function(file) {
  read.csv(file, sep = ",")
}

edgelists_RRA <- files_RRA |>
  set_names(
    basename(files_RRA) |>
      gsub("\\.csv$", "", x = _) |>
      gsub("^edgelist_", "", x = _) |>
      tolower()
  ) |>
  map(\(f) read_edgelist(f)[, 2:6])

visdata_RRA_monthly <- function(edgelist,
                                year,
                                month) {
  
  #Convert to wide format datasets
  matrixRRA <- edgelist |> 
    select(Predator, Prey, MeanRRA) |> 
    pivot_wider(
      names_from = "Predator",
      values_from = "MeanRRA",
      id_cols = "Prey",
      values_fill = 0
    )
  
  matrixRRA <- matrixRRA |> 
    column_to_rownames("Prey") |> 
    as.matrix()
  
}

all_visdata_RRA <- purrr::imap(edgelists_RRA, function(edgelist, nm) {
  
  parts <- strsplit(nm, "_")[[1]]
  
  month <- tools::toTitleCase(parts[1])
  year <- parts[2]
  
  visdata_RRA_monthly(
    edgelist,
    year,
    month
  )
  
})

visweb(t(all_visdata_RRA[[5]]))

#RRA


#Read edgelist data

files <- list.files(
  pattern = "^edgelist_.*_20.*_Fam_RRA\\.csv$"
)

read_edgelist <- function(file) {
  read.csv(file, sep = ",")
}

edgelists <- files |>
  set_names(
    basename(files) |>
      gsub("\\.csv$", "", x = _) |>
      gsub("^edgelist_", "", x = _) |>
      tolower()
  ) |>
  map(\(f) read_edgelist(f)[, 2:6])

#Calculate biomass and count per predator-prey interaction

process_edgelist <- function(edgelist,
                             prey_sum,
                             mean_biomass,
                             year,
                             month) {
  
  # Prey matrix for the correct month/year
  prey_matrix <- prey_sum |>
    select(-Order) |>
    filter(
      Year == year,
      Month == month
    ) |>
    select(-c(Year, Month)) |>
    pivot_wider(
      names_from = Taxon,
      values_from = TotalCount,
      values_fill = 0
    ) |>
    as.matrix()
  
  # Food web matrix
  web_matrix <- edgelist |>
    select(Predator, Prey, MeanRRA) |> 
    pivot_wider(
      names_from = Predator,
      values_from = MeanRRA,
      values_fill = 0
    ) |>
    column_to_rownames("Prey") |>
    as.matrix()
  
  web_matrix <- t(web_matrix)
  
  # Keep only shared taxa
  shared_taxa <- intersect(
    colnames(prey_matrix),
    colnames(web_matrix)
  )
  
  prey_matrix <- prey_matrix[, shared_taxa, drop = FALSE]
  web_matrix  <- web_matrix[, shared_taxa, drop = FALSE]
  
  # Relative abundance
  prey_clean <- as.data.frame(prey_matrix) |>
    pivot_longer(
      everything(),
      names_to = "Taxon",
      values_to = "Count"
    ) |>
    left_join(
      mean_biomass,
      by = "Taxon"
    ) |>
    mutate(
      Biomass = Count * mean_weight
    )
  
  # Interaction strengths
  interactions_long <- as.data.frame(web_matrix) |>
    rownames_to_column("Predator") |>
    pivot_longer(
      -Predator,
      names_to = "Taxon",
      values_to = "MeanRRA"
    )
  
  # Final output
  interactions_long |>
    left_join(
      prey_clean |>
        select(
          Taxon,
          Count,
          Biomass
        ),
      by = "Taxon"
    ) |>
    mutate(
      Count = replace_na(Count, 0),
      Biomass = replace_na(Biomass, 0)
    )
}

combined_data <- purrr::imap(edgelists, function(edgelist, nm) {
  
  parts <- strsplit(nm, "_")[[1]]
  
  month <- tools::toTitleCase(parts[1])
  year <- parts[2]
  
  process_edgelist(
    edgelist,
    prey_sum,
    mean_biomass,
    year,
    month
  )
  
})

combined_dates <- bind_rows(
  combined_data,
  .id = "Dataset"
)

#FOO

#Read edgelist data

files <- list.files(
  pattern = "^edgelist_.*_20.*_Fam\\.csv$"
)

read_edgelist <- function(file) {
  read.csv(file, sep = ",")
}

edgelists <- files |>
  set_names(
    basename(files) |>
      gsub("\\.csv$", "", x = _) |>
      gsub("^edgelist_", "", x = _) |>
      tolower()
  ) |>
  map(\(f) read_edgelist(f)[, 2:6])

#Calculate biomass and count per predator-prey interaction

process_edgelist <- function(edgelist,
                             prey_sum,
                             mean_biomass,
                             year,
                             month) {
  
  # Prey matrix for the correct month/year
  prey_matrix <- prey_sum |>
    select(-Order) |>
    filter(
      Year == year,
      Month == month
    ) |>
    select(-c(Year, Month)) |>
    pivot_wider(
      names_from = Taxon,
      values_from = TotalCount,
      values_fill = 0
    ) |>
    as.matrix()
  
  
  shared_taxa <- intersect(
    colnames(prey_matrix),
    unique(edgelist$Prey)
  )
  
  interactions_long <- edgelist |>
    filter(Prey %in% shared_taxa) |>
    select(
      Predator,
      Prey,
      FOO,
      TotalPresence,
      N
    ) |>
    rename(Taxon = Prey)
  
  # Relative abundance
  prey_clean <- as.data.frame(prey_matrix) |>
    pivot_longer(
      everything(),
      names_to = "Taxon",
      values_to = "Count"
    ) |>
    left_join(
      mean_biomass,
      by = "Taxon"
    ) |>
    mutate(
      Biomass = Count * mean_weight
    )
  
  # Final output
  interactions_long |>
    left_join(
      prey_clean |>
        select(
          Taxon,
          Count,
          Biomass
        ),
      by = "Taxon"
    ) |>
    mutate(
      Count = replace_na(Count, 0),
      Biomass = replace_na(Biomass, 0)
    )
}

combined_data_FOO <- purrr::imap(edgelists, function(edgelist, nm) {
  
  parts <- strsplit(nm, "_")[[1]]
  
  month <- tools::toTitleCase(parts[1])
  year <- parts[2]
  
  process_edgelist(
    edgelist,
    prey_sum,
    mean_biomass,
    year,
    month
  )
  
})


combined_dates_FOO <- bind_rows(
  combined_data_FOO,
  .id = "Dataset"
)


#COMBINE DATASETS

#1. Exchange datasets for month and year columns

FOO_clean <- combined_dates_FOO |>
  mutate(
    Month = str_extract(Dataset, "^[^_]+"),
    Year  = str_extract(Dataset, "\\d{4}")
  ) |> 
  select(-Dataset)

RRA_clean <- combined_dates |>
  mutate(
    Month = str_extract(Dataset, "^[^_]+"),
    Year  = str_extract(Dataset, "\\d{4}")
  ) |> 
  select(-Dataset)

#2. left join datasets

df <- RRA_clean |>
  left_join(
    FOO_clean,
    by = c("Predator", "Taxon", "Count", "Biomass", "Year", "Month")
  )

#3 group birds and spiders

df <- df |> 
  mutate(Group = case_when(
    Predator %in% c("Lapland_bunting", "Snow_bunting",
                    "Wheatear", "Redpoll", "Ptarmigan") ~ "Birds",
    Predator %in% c("Lycosidae", "Araneidae",
                    "Dictynidae", "Linyphiidae", "Phalangiidae",
                    "Philodromidae", "Thomisidae") ~ "Spiders",
    Predator %in% c("Dytiscidae", "Gyrinidae") ~ "Aquatic_predators",
    Predator %in% "Chironomidae" ~ "Predatory_midges",
    Predator %in% "Scathophagidae" ~ "Predatory_flies"
  ))

#Drop all zero biomass
df_clean <- df |> 
  filter(Biomass > 0)

df_clean <- df_clean |> 
  mutate(
    Biomass_log = log(Biomass),
    Count_log = log(Count))

#Make year and month into factors
df_clean <- df_clean |> 
  mutate(
    across(
      c("Month", "Year", "Predator", "Taxon"), as.factor)
  )


#4 fit the models RRA
#We don't add Prey as a factor due to lack of data
#We use ziformula to account for zeros - 
#Zeros may be true but also a consequence of read thresholds.
#Interactions and spline did not improve the models. If they did,
#Those models are used instead.
#compare AIC with normal model to see if log model is better

ggplot(data = subset(df_clean, Predator == "Redpoll"), aes(x = Biomass_log,
                                                           y = MeanRRA)) +
  geom_point()

bird_df <- subset(df_clean, Group == "Birds") |> 
  droplevels()

RRA_bio_mod <- glmmTMB::glmmTMB(MeanRRA ~ poly(Biomass_log, 2) +
                                  Month +
                                  Year +
                                  Predator,
                                ziformula = ~1,
                                dispformula = ~Month,
                                data = bird_df,
                                family = glmmTMB::beta_family())


RRA_bio_spline <- glmmTMB::glmmTMB(
  MeanRRA ~ splines::ns(Biomass_log, df = 3) +
    Month +
    Year +
    Predator,
  ziformula = ~1,
  dispformula = ~Month,
  data = subset(df_clean, Group == "Birds"),
  family = glmmTMB::beta_family()
)

AIC(RRA_bio_mod, RRA_bio_mod_pred)
#No great difference- stick with the simpler quadratic model

car::Anova(RRA_bio_mod, type = "2")
#Differences between months and a non-linear relationship between 
#biomass and meanRRA!

res_RB <- simulateResiduals(RRA_bio_mod)

plot(res_RB)
testDispersion(res_RB)
testUniformity(res_RB)
plotResiduals(res_RB, subset(df_clean, Group == "Birds")$Predator)
plotResiduals(res_RB, subset(df_clean, Group == "Birds")$Month)
plotResiduals(res_RB, subset(df_clean, Group == "Birds")$Year)

#Predict

library(marginaleffects)

bio_values <- seq(
  min(bird_df$Biomass_log, na.rm = TRUE),
  max(bird_df$Biomass_log, na.rm = TRUE),
  length.out = 100
)

pred_RRA_bio <- avg_predictions(
  RRA_bio_mod,
  newdata = datagrid(
    Biomass_log = bio_values,
    grid_type = "counterfactual"
  ),
  by = "Biomass_log",
  re.form = NA
)

RRA_bio_plot <- ggplot(pred_RRA_bio, aes(x = Biomass_log, y = estimate)) +
  #  geom_point(
  #    data = subset(df_clean, Group == "Birds"),
  #    aes(x = Biomass_log,
  #        y = MeanRRA)
  #  ) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_classic() +
  labs(
    x = "Log biomass",
    y = "Predicted MeanRRA"
  )



#Count mod
#Basic mod best again
RRA_count_mod <- glmmTMB::glmmTMB(MeanRRA ~ poly(Count_log, 2) + 
                                    Month +
                                    Year +
                                    Predator,
                                  ziformula = ~1,
                                  dispformula = ~Month,
                                  data = bird_df,
                                  family = glmmTMB::beta_family())

car::Anova(RRA_count_mod, type = "2")

#Same results as for biomass!

res_RC <- simulateResiduals(RRA_count_mod)

plot(res_RC)
testDispersion(res_RC)
testUniformity(res_RC)
plotResiduals(res_RC, bird_df$Predator)
plotResiduals(res_RC, bird_df$Month)
plotResiduals(res_RC, bird_df$Year)


count_values <- seq(
  min(bird_df$Count_log, na.rm = TRUE),
  max(bird_df$Count_log, na.rm = TRUE),
  length.out = 100
)

pred_RRA_count <- avg_predictions(
  RRA_count_mod,
  newdata = datagrid(
    Count_log = count_values,
    grid_type = "counterfactual"
  ),
  by = "Count_log",
  re.form = NA
)

RRA_count_plot <- ggplot(pred_RRA_count, aes(x = Count_log, y = estimate)) +
  #  geom_point(
  #    data = subset(df_clean, Group == "Birds"),
  #    aes(x = Biomass_log,
  #        y = MeanRRA)
  #  ) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_classic() +
  labs(
    x = "Log abundance",
    y = "Predicted MeanRRA"
  )


#5 fit the models FOO
#Betabinomial regression to handle additional variability
#around the underlying probability
#Poly still the best model

ggplot(bird_df, aes(x = Biomass_log,
                    y = FOO)) +
  geom_point()

FOO_bio_mod <- glmmTMB::glmmTMB(cbind(TotalPresence, N - TotalPresence) ~ 
                                  poly(Biomass_log, 2) +
                                  Month +
                                  Year +
                                  Predator,
                                dispformula = ~ Month,
                                data = bird_df,
                                family = glmmTMB::betabinomial())

AIC(FOO_bio_mod,FOO_bio_poly)

car::Anova(FOO_bio_mod, type = "2")


#Model still is not perfect...

res_FB <- simulateResiduals(FOO_bio_mod)

plot(res_FB)
testDispersion(res_FB)
testUniformity(res_FB)
plotResiduals(res_FB, bird_df$Predator)
plotResiduals(res_FB, bird_df$Month)
plotResiduals(res_FB, bird_df$Year)

pred_FOO_bio <- avg_predictions(
  FOO_bio_mod,
  newdata = datagrid(
    Biomass_log = bio_values,
    grid_type = "counterfactual"
  ),
  by = "Biomass_log",
  re.form = NA
)

#Decide which prey taxa to color
# Find the 10 most common taxa
top10 <- bird_df |>
  count(Taxon, sort = TRUE) |>
  slice_head(n = 8) |>
  pull(Taxon) |> 
  droplevels()

# Create a grouping variable for plotting
plot_bird <- bird_df |>
  mutate(
    Highlight = case_when(
      Taxon %in% top10 ~ as.character(Taxon),
      TRUE ~ "Other"
    )
  )

#remove Gemoetridae since they are hard to represent
plot_bird$Highlight[plot_bird$Highlight == "Geometridae"] <- "Other"


FOO_bio_plot <- ggplot(pred_FOO_bio, aes(x = Biomass_log, y = estimate)) +
  geom_point(
    data = plot_bird,
    aes(x = Biomass_log,
        y = FOO,
        color = Highlight,
        alpha = Highlight),
    size = 3
  ) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c(
      Anthomyiidae = "#E69F00", 
        Chironomidae = "#56B4E9",
        Culicidae = "#009E73",
        Empididae = "#F0E442",
      Ichneumonidae = "#0072B2",
      Scatophagidae = "#D55E00",
      Sciaridae = "#CC79A7",
      Other = "black"
    )
  ) +
  scale_alpha_manual(
    values = c(
      setNames(rep(1, length(top10)), top10),
      Other = 0.5
    ),
    guide = "none"
  ) +
  theme_classic() +
  labs(
    x = "Log biomass",
    y = "Predicted FOO",
    color = "Prey taxon"
  ) +
  theme(
    legend.position = "none"
  )

FOO_bio_plot

#Count mod
#Here, a linear relationship was the best fit
#No dispersion formula needed either

ggplot(bird_df, aes(x = Count_log,
                    y = FOO)) +
  geom_point()

FOO_count_mod <- glmmTMB::glmmTMB(cbind(TotalPresence, N - TotalPresence) ~ 
                                    Count_log +
                                    Month + 
                                    Year +
                                    Predator,
                                  data = bird_df,
                                  family = glmmTMB::betabinomial())

AIC(FOO_count_mod, FOO_count_org)
car::Anova(FOO_count_log, type = "2")

res_FC <- simulateResiduals(FOO_count_mod)

plot(res_FC)
testDispersion(res_FC)
testUniformity(res_FC)
plotResiduals(res_FC, bird_df$Predator)
plotResiduals(res_FC, bird_df$Month)
plotResiduals(res_FC, bird_df$Year)

pred_FOO_count <- avg_predictions(
  FOO_count_mod,
  newdata = datagrid(
    Count_log = count_values,
    grid_type = "counterfactual"
  ),
  by = "Count_log",
  re.form = NA
)

FOO_count_plot <- ggplot(pred_FOO_count, aes(x = Count_log, y = estimate)) +
  geom_point(
    data = plot_bird,
    aes(x = Count_log,
        y = FOO,
        color = Highlight,
        alpha = Highlight),
    size = 3
  ) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c(
      Anthomyiidae = "#E69F00", 
      Chironomidae = "#56B4E9",
      Culicidae = "#009E73",
      Empididae = "#F0E442",
      Ichneumonidae = "#0072B2",
      Scatophagidae = "#D55E00",
      Sciaridae = "#CC79A7",
      Other = "black"
    )
  ) +
  scale_alpha_manual(
    values = c(
      setNames(rep(1, length(top10)), top10),
      Other = 0.5
    ),
    guide = "none"
  ) +
  theme_classic() +
  labs(
    x = "Log abundance",
    y = "Predicted FOO",
    color = "Prey taxon"
  )

FOO_count_plot

#Combine plots

combined_preds <- (RRA_bio_plot | RRA_count_plot) /
  (FOO_bio_plot | FOO_count_plot)

ggsave(combined_preds, filename = "birdpredation_plots_clr.png",
       dpi = 300, height = 8, width = 12)

