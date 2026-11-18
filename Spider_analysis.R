library(readxl)
library(readr)
library(tidyverse)
library(bipartite)
library(patchwork)
library(forcats)
library(vegan)
library(DHARMa)


rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

#1. Read metadata for spiders
df_24 <- read.csv("Foodweb_samples_2024.csv", sep = ";") |> 
  filter(Group == "Spiders") |> 
  select(ID, Size.Gender, Location, Date) |> 
  mutate(Location = fct_recode(
    factor(Location),
                     "Plot_1" = "Spider_1",
                    "Plot_2" = "Spider_2" 
  ))


df_25 <- read.csv("FoodwebSamples_2025.csv", sep = ";")|> 
  filter(Group == "Spiders") |> 
  select(ID, Size.Gender, Location, Date)

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

#Read diet matrix

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
                        "Salticidae", "Sarcophagidae", "Lycaenidae"))

#Keep only the beetle samples
#1. Find beetle samples
spider_foodweb <- combined_df %>%
  filter(Sample %in% df_tot$ID)

spider_foodweb <- spider_foodweb |> 
  left_join(
    df_tot %>%
      select(ID, Location, doy, year),
    by = c("Sample" = "ID")
  )

#Add predator names to the samples
#Sum reads across remaining primers
summed_reads <- spider_foodweb %>%
  group_by(Sample, Family) %>%
  summarise(total_reads = sum(Reads), .groups = "drop")

levels(as.factor(summed_reads$Family))

#Identify families that cannot be spider predators (all non-spiders)
non_predator_species <- c(
  "Apidae",
  "Anatidae",
  "Anthomyiidae",
  "Aphididae",
  "Baetidae",
  "Braconidae",
  "Byrrhidae",
  "Carabidae",
  "Calliphoridae" ,
  "Canidae",
  "Carabidae", 
  "Cecidomyiidae",
  "Ceratopogonidae",
  "Cervidae",
  "Chironomidae",
  "Chloropidae",
  "Cicadellidae",
  "Coccinellidae",
  "Coleophoridae",
  "Corvidae",
  "Crambidae",
  "Culicidae",
  "Curculionidae",
  "Dolichopodidae",
  "Dytiscidae",
  "Ectopsocidae",
  "Empididae",
  "Entomobryidae" ,
  "Ephydridae",
  "Erebidae" ,
  "Eriocraniidae",
  "Fringillidae",
  "Gelechiidae" ,
  "Geometridae",
  "Gyrinidae",
  "Heleomyzidae",
  "Hemerobiidae",
  "Hominidae",
  "Ichneumonidae",
  "Isotomidae" ,
  "Laridae",
  "Lebertiidae" ,
  "Limnephilidae",
  "Limoniidae",
  "Lygaeidae",
  "Muscicapidae",
  "Muscidae",
  "Mycetophilidae",
  "Noctuidae",
  "Orchesellidae",
  "Oribatulidae",
  "Pediciidae",
  "Phocidae",
  "Phoridae",
  "Phryganeidae",
  "Piophilidae",
  "Procellariidae",
  "Psychodidae" ,
  "Psyllidae",
  "Pyralidae",
  "Scathophagidae",
  "Sciaridae",
  "Scythrididae",
  "Simuliidae" ,
  "Sphaeroceridae",
  "Staphylinidae",
  "Syrphidae",
  "Tachinidae" ,
  "Tenthredinidae",
  "Tetragnathidae",
  "Tipulidae",
  "Tortricidae",
  "Vespidae" 
)

dominant_spider_sp <- summed_reads %>%
  filter(!Family %in% non_predator_species) %>% 
  group_by(Sample) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample, Family)
dominant_spider_sp$Predator <- dominant_spider_sp$Family

#Check count per family
summary <- dominant_spider_sp |> 
  group_by(Predator) |> 
  summarize(n = n_distinct(Sample))

#Gnaphosidae has only one sample - consider removing

spider_df <- spider_foodweb |> 
  left_join(dominant_spider_sp |> 
              select(Sample, Predator),
            by = "Sample")

#Group data by prey family
spider_df <- spider_df |> 
  select(Family, Sample, Reads, Location ,doy, year, Predator) |> 
  group_by(Family, Sample) |> 
  mutate(SumReads = sum(Reads, na.rm = TRUE)) |> 
  select(-Reads)

#Filter out all self-reads
spider_df <- spider_df |> 
  filter(Predator != Family)

#Filter out all families with no reads across any samples
spider_df <- spider_df |> 
  group_by(Family) |> 
  filter(sum(SumReads, na.rm = TRUE) > 0) |> 
  ungroup()

#look at prey families
levels(as.factor(spider_df$Family))

#Check if some samples had no prey reads
no_reads <- spider_df |>
  group_by(Sample) |>
  summarise(TotalReads = sum(SumReads, na.rm = TRUE)) |>
  filter(TotalReads == 0)

#Filter those samples
spider_df_clean <- spider_df |> 
  filter(!Sample %in% no_reads$Sample)

spider_df_clean <- spider_df_clean |> 
  group_by(Sample) |> 
  mutate(
    RRA = SumReads / sum(SumReads)
  )

test <- spider_df_clean |> 
  group_by(Sample) |> 
  summarize(TotRRA = sum(RRA))

#VISWEB
family_grouping <- spider_df_clean %>%
  group_by(Predator, Family) %>%
  summarise(
    N_present = n_distinct(Sample[SumReads > 0]),
    MeanRRA = mean(RRA),
    .groups = "drop"
  ) %>%
  left_join(
    spider_df %>%
      group_by(Predator) %>%
      summarise(
        N_samples = n_distinct(Sample),
        .groups = "drop"
      ),
    by = "Predator"
  ) %>%
  mutate(
    FOO = N_present / N_samples
  ) |> 
  filter(!Predator %in% c("Gnaphosidae", "Dictynidae"))

#Turn to simple matrix for visweb
visweb_mat_FOO <- family_grouping |> 
  pivot_wider(
    values_from = FOO,
    names_from = Family,
    values_fill = 0,
    id_cols = Predator
  )%>%
  column_to_rownames("Predator")

visweb_mat_RRA <- family_grouping |> 
  pivot_wider(
    values_from = MeanRRA,
    names_from = Family,
    values_fill = 0,
    id_cols = Predator
  )%>%
  column_to_rownames("Predator")

visweb(as.matrix(visweb_mat_RRA),
       labsize = 0.6)

visweb(as.matrix(visweb_mat_FOO),
       labsize = 0.6)

#3. Visualize diet composition across sites for wolfspiders

spider_df_clean

#4. Run permanova to investigate if differences are significant
