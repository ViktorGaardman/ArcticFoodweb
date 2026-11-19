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
  mutate(
    Location = fct_collapse(
      factor(Location),
      Plot_1 = c("Spider_1", "Spider1"),
      Plot_2 = c("Spider_2", "Spider2")
    )
  )


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
  select(Family, Sample, Reads, Location, doy, year, Predator) |>
  group_by(Family, Sample, Location, doy, year, Predator) |>
  summarise(
    Reads = sum(Reads, na.rm = TRUE),
    .groups = "drop"
  )

#Filter out all self-reads
spider_df <- spider_df |> 
  filter(Predator != Family)

#Filter out all families with no reads across any samples
spider_df <- spider_df |> 
  group_by(Family) |> 
  filter(sum(Reads, na.rm = TRUE) > 0) |> 
  ungroup()

#look at prey families
levels(as.factor(spider_df$Family))

#Check if some samples had no prey reads
no_reads <- spider_df |>
  group_by(Sample) |>
  summarise(TotalReads = sum(Reads, na.rm = TRUE)) |>
  filter(TotalReads == 0)

#Filter those samples
spider_df_clean <- spider_df |> 
  filter(!Sample %in% no_reads$Sample)

spider_df_clean <- spider_df_clean |> 
  group_by(Sample) |> 
  mutate(
    RRA = Reads / sum(Reads)
  )

test <- spider_df_clean |> 
  group_by(Sample) |> 
  summarize(TotRRA = sum(RRA))

#VISWEB
family_grouping <- spider_df_clean %>%
  group_by(Predator, Family) %>%
  summarise(
    N_present = n_distinct(Sample[Reads > 0]),
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

#Make sure the order of predators are the same
pred_order <- c("Araneidae", 
                "Lycosidae",
                "Phalangiidae", "Philodromidae",
                "Thomisidae")
prey_order <- c("Araneidae", "Dictynidae", "Linyphiidae",
                "Lycosidae",
                "Phalangiidae", "Philodromidae",
                "Thomisidae", "Chironomidae", 
                "Culicidae", "Muscidae", "Scathophagidae",
                 "Sciaridae", "Simuliidae", "Coccinellidae",
                "Dytiscidae",
                "Gyrinidae", "Geometridae")

visweb_mat_FOO <- visweb_mat_FOO[pred_order, prey_order]
visweb_mat_RRA <- visweb_mat_RRA[pred_order, prey_order]

for(taxon in intersect(rownames(visweb_mat_FOO),
                       colnames(visweb_mat_FOO))) {
  visweb_mat_FOO[taxon, taxon] <- NA
}

visweb(as.matrix(visweb_mat_FOO),
       labsize = 0.6,
       type = "none",
       NA.col = "firebrick")

for(taxon in intersect(rownames(visweb_mat_RRA),
                       colnames(visweb_mat_RRA))) {
  visweb_mat_RRA[taxon, taxon] <- NA
}

visweb(as.matrix(visweb_mat_RRA),
       labsize = 0.6,
       type = "none",
       NA.col = "firebrick")

png(
  "spider_visweb_RRA.png",
  width = 1200,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_mat_RRA),
       labsize = 0.6,
       type = "none",
       NA.col = "firebrick")

dev.off()

png(
  "spider_visweb_FOO.png",
  width = 1200,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_mat_FOO),
       labsize = 0.6,
       type = "none",
       NA.col = "firebrick")

dev.off()

#Remember to try and color self-boxes


#3. Visualize diet composition across sites for Lycosidae and Aranidae

with(
  spider_df_clean %>% distinct(Sample, Predator, Location),
  table(Predator, Location)
)

spider_df_clean <-spider_df_clean |> 
  mutate(
    Habitat = fct_collapse(
      factor(Location),
      "Dry" = c("Plot_1", "Plot_5"),
      "River" = "Plot_3",
      "Wet" = c("Plot_2", "Plot_4")
    )
  )

#Filter only Lycosidae and Araneidae
perm_df <- spider_df_clean |> 
  filter(Predator %in% c("Lycosidae"))

#Remove all zero interactions from this dataset
#Check if some samples had no prey reads
no_reads <- perm_df |>
  group_by(Family) |> 
  summarise(TotalReads = sum(Reads, na.rm = TRUE)) |>
  filter(TotalReads == 0)

#Also remove families with only one/two appearances
small_fam <- perm_df |>
  filter(Reads > 0) |>
  group_by(Family) |>
  summarise(TotalPres = n_distinct(Sample)) |> 
  filter(TotalPres < 4)

#Filter those samples
perm_df_clean <- perm_df |> 
  filter(!Family %in% no_reads$Family) |> 
  filter(!Family %in% small_fam$Family)

#Remove any sampls with no more reads
zero_samp <- perm_df_clean |> 
  group_by(Sample) |> 
  summarize(TotReads = sum(Reads)) |> 
  filter(TotReads == 0)

perm_df_clean <- perm_df_clean |> 
  filter(!Sample %in% zero_samp$Sample)

comm_matrix <- perm_df_clean %>%
  select(Sample, Family, RRA) %>%
  pivot_wider(
    id_cols = Sample,
    names_from = Family,
    values_from = RRA,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample")

meta <- perm_df_clean %>%
  distinct(Sample, Habitat, Predator, Location, doy, year)

meta <- meta[match(rownames(comm_matrix), meta$Sample), ]

#4. Run NMDS + permanova to investigate if differences are significant

nmds <- metaMDS(comm_matrix,
                distance = "jaccard",
                k = 2, trymax = 100)

nmds_bray <- metaMDS(comm_matrix,
                distance = "bray",
                k = 2, trymax = 100)

plot_df <- data.frame(
  NMDS1 = nmds_bray$points[,1],
  NMDS2 = nmds_bray$points[,2],
  Location = meta$Habitat
)

ggplot(plot_df,
       aes(NMDS1, NMDS2, colour = Location)) +
  geom_point(size = 3)

##Permanova with year as a random factor
#permutations <- with(meta, how(nperm=999, blocks = Location))

#Calculate distance matrix
dist_matrix <- vegdist(comm_matrix, method = "bray")

###Fit permanova model
adonis2(dist_matrix ~ Habitat  + doy + year,
        by = "margin",
        data=meta)

pairwise.adonis2(dist_matrix ~ Habitat, data=meta,
                 permutations=999,  p.adjust.m="none")

#Check assumption of homogeneity of multivariate dispersion
anova(betadisper(dist_matrix, meta$Location))

#Visualize as viswebs again

visweb_Lyc <- spider_df_clean %>%
  filter(Predator == "Lycosidae") |> 
  group_by(Family, Habitat) %>%
  summarise(
    N_present = n_distinct(Sample[Reads > 0]),
    MeanRRA = mean(RRA),
    .groups = "drop"
  ) %>%
  left_join(
    spider_df_clean %>%
      group_by(Habitat) %>%
      summarise(
        N_samples = n_distinct(Sample),
        .groups = "drop"
      ),
    by = "Habitat"
  ) %>%
  mutate(
    FOO = N_present / N_samples
  ) 

visweb_FOO_Lyc <- visweb_Lyc |> 
  pivot_wider(
    values_from = FOO,
    names_from = Family,
    values_fill = 0,
    id_cols = Habitat
  )%>%
  column_to_rownames("Habitat")

visweb_RRA_Lyc <- visweb_Lyc |> 
  pivot_wider(
    values_from = MeanRRA,
    names_from = Family,
    values_fill = 0,
    id_cols = Habitat
  )%>%
  column_to_rownames("Habitat")

png(
  "lyc_visweb_FOO.png",
  width = 800,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_FOO_Lyc),
       labsize = 0.4,
       type = "none")

dev.off()

png(
  "lyc_visweb_RRA.png",
  width = 800,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_RRA_Lyc),
       labsize = 0.4,
       type = "none")

dev.off()

visweb_ara <- spider_df_clean %>%
  filter(Predator == "Araneidae") |> 
  group_by(Family, Habitat) %>%
  summarise(
    N_present = n_distinct(Sample[Reads > 0]),
    MeanRRA = mean(RRA),
    .groups = "drop"
  ) %>%
  left_join(
    spider_df_clean %>%
      group_by(Habitat) %>%
      summarise(
        N_samples = n_distinct(Sample),
        .groups = "drop"
      ),
    by = "Habitat"
  ) %>%
  mutate(
    FOO = N_present / N_samples
  ) 

visweb_FOO_ara <- visweb_ara |> 
  pivot_wider(
    values_from = FOO,
    names_from = Family,
    values_fill = 0,
    id_cols = Habitat
  )%>%
  column_to_rownames("Habitat")

visweb_RRA_ara <- visweb_ara |> 
  pivot_wider(
    values_from = MeanRRA,
    names_from = Family,
    values_fill = 0,
    id_cols = Habitat
  )%>%
  column_to_rownames("Habitat")

png(
  "ara_visweb_FOO.png",
  width = 800,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_FOO_ara),
       labsize = 0.4,
       type = "none")

dev.off()

png(
  "ara_visweb_RRA.png",
  width = 800,
  height = 500,
  res = 300
)

visweb(as.matrix(visweb_RRA_ara),
       labsize = 0.4,
       type = "none")

dev.off()

#Try dbRDA

cap_mod <- capscale(
  comm_matrix ~ Habitat + Predator + year,
  data = meta,
  distance = "jaccard"
)

anova(cap_mod)
anova(cap_mod, by = "terms")

plot(cap_mod)
