library(tidyverse)
library(ggjoy)

setwd("~/Box Sync/GHdata/BgeVars/freq/")

vars <- read.table("BGE.VB.AORO.tsv", header = TRUE) %>%
  filter(CHROM != "mitochondrion,")
colnames(vars)[5] <- "RO"
colnames(vars)[1] <- "Scaffold"
# filter out ratios that wouldn't occur in octoploid SNPs
df <- vars %>%
  mutate(AP = AO / (AO + RO), RP = RO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

ggplot(df, aes(AP)) +
  geom_histogram()

map <- read.table("scaffold_mapping.txt", sep = " ", header = TRUE)

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
name <- strsplit(as.character(df$Scaffold), split = "_")
name <- matrix(unlist(name), ncol = 3, byrow = TRUE)
name <- data.frame(name)
name <- rbind(name, name[753093,])
df$Scaffold <- name$X3

df <- left_join(df, map) %>%
  filter(LGs != "NA")

df2 <-  filter(df, AO + RO > 100, AO + RO < 10000)
all_vars <- df2

save(all_vars, file = "all_vars.rda")
load("allele.rda")
load("all.df.rda")

df3 <- inner_join(all_vars, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000) %>%
  sample_frac(0.05)
  
dot.p <- ggplot(df3, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme(legend.position = "none")
dot.p

joy.p <- ggplot(df3, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme(legend.position = "none")
joy.p
