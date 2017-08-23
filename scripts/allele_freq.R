library(tidyverse)
library(ggjoy)

setwd("~/Box Sync/GHdata/BgeVars/")

vars <- read.table("snp_filter.AO-RO.tsv", header = TRUE) %>%
  filter(CHROM != "mitochondrion,")
colnames(vars)[6] <- "RO"
colnames(vars)[1] <- "Scaffold"
df <- vars %>%
  filter(AO > 0, RO > 0) %>%
  mutate(AP = AO / (AO + RO), RP = RO / (AO + RO))

map <- read.table("scaffold_mapping.txt", sep = " ", header = TRUE)

df <- left_join(df, map) %>%
  filter(LGs != "NA")

df <-  filter(df, AO + RO > 400, AO > 100)

save(df, file = "allele.rda")
load("allele.rda")
load("all.df.rda")

df2 <- inner_join(df, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000) %>%
  sample_frac(0.25)
  
dot.p <- ggplot(df2, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_discrete(name = "Alternate Allele Fraction", limits = c(0, 0.25, 0.5, 0.75, 1)) +
  ylim(0,1) +
  theme(legend.position = "none")
dot.p

joy.p <- ggplot(df2, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_discrete(name = "Density of Alternative Allele Fraction", limits = c(0, 0.25, 0.5, 0.75, 1, 1.25)) + 
  xlim(0,1.25) + 
  theme(legend.position = "none")
joy.p
