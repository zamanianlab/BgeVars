library(tidyverse)
library(ggjoy)

setwd("~/Box Sync/GHdata/BgeVars/freq/")

map <- read.table("scaffold_mapping.txt", sep = " ", header = TRUE)

##############################
######   ALL VARIANTS   ######
##############################

vars <- read.table("../bcftools.dp4.vcf", col.names = c("Scaffold", "Position", "Reference", "Alternative", "DP4"), sep = " ") %>%
  filter(Scaffold != "mitochondrion,") %>%
  separate(DP4, c("F_REF", "R_REF", "F_ALT", "R_ALT"), sep = ",")

# filter out ratios that wouldn't occur in octoploid SNPs
df <- vars %>% transform(F_REF = as.numeric(F_REF), R_REF = as.numeric(R_REF), F_ALT = as.numeric(F_ALT), R_ALT = as.numeric(R_ALT)) %>%
  mutate(RO = F_REF + R_REF, AO = F_ALT + R_ALT) %>%
  mutate(RP = RO / (AO + RO), AP = AO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

ggplot(df, aes(AP)) +
  geom_histogram()

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
name <- strsplit(as.character(df$Scaffold), split = "_")
name <- matrix(unlist(name), ncol = 3, byrow = TRUE)
name <- data.frame(name)
df$Scaffold <- name$X3

df2 <- left_join(df, map) %>%
  filter(LGs != "NA")

df3 <-  filter(df2, AO + RO > 100, AO + RO < 10000)
all_vars <- df3

save(all_vars, file = "all_vars.rda")
load("all_vars.rda")
load("all.df.rda")

df4 <- inner_join(all_vars, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000)

dot.p <- ggplot(df4, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme(legend.position = "none")
dot.p

joy.p <- ggplot(df4, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme(legend.position = "none")
joy.p

hist.p <- ggplot(df4, aes(AP)) +
  geom_histogram() +
  facet_wrap(~LGs)
hist.p

##############################
######       SNPS       ######
##############################

snps <- read.table("../bcftools.snp.dp4.vcf", col.names = c("Scaffold", "Position", "Reference", "Alternative", "DP4"), sep = " ") %>%
  filter(Scaffold != "mitochondrion,") %>%
  separate(DP4, c("F_REF", "R_REF", "F_ALT", "R_ALT"), sep = ",")

# filter out ratios that wouldn't occur in octoploid SNPs
snp.df <- snps %>% transform(F_REF = as.numeric(F_REF), R_REF = as.numeric(R_REF), F_ALT = as.numeric(F_ALT), R_ALT = as.numeric(R_ALT)) %>%
  mutate(RO = F_REF + R_REF, AO = F_ALT + R_ALT) %>%
  mutate(RP = RO / (AO + RO), AP = AO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

ggplot(snp.df, aes(AP)) +
  geom_histogram()

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
snp_name <- strsplit(as.character(snp.df$Scaffold), split = "_")
snp_name <- matrix(unlist(snp_name), ncol = 3, byrow = TRUE)
snp_name <- data.frame(snp_name)
snp.df$Scaffold <- snp_name$X3

snp.df2 <- left_join(snp.df, map) %>%
  filter(LGs != "NA")

snp.df3 <-  filter(snp.df2, AO + RO > 100, AO + RO < 10000)

save(snp.df3, file = "snps.rda")
load("snps.rda")
load("all.df.rda")

snp.df4 <- inner_join(snp.df3, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000)

dot.p <- ggplot(snp.df4, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme(legend.position = "none")
dot.p

joy.p <- ggplot(snp.df4, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme(legend.position = "none")
joy.p

hist.p <- ggplot(snp.df4, aes(AP)) +
  geom_histogram() +
  facet_wrap(~LGs)
hist.p

##############################
######     FILTER 1     ######
##############################

fil1 <- read.table("../bcftools.snp.fil1.dp4.vcf", col.names = c("Scaffold", "Position", "Reference", "Alternative", "DP4"), sep = " ") %>%
  filter(Scaffold != "mitochondrion,") %>%
  separate(DP4, c("F_REF", "R_REF", "F_ALT", "R_ALT"), sep = ",")

# filter out ratios that wouldn't occur in octoploid SNPs
fil1.df <- fil1 %>% transform(F_REF = as.numeric(F_REF), R_REF = as.numeric(R_REF), F_ALT = as.numeric(F_ALT), R_ALT = as.numeric(R_ALT)) %>%
  mutate(RO = F_REF + R_REF, AO = F_ALT + R_ALT) %>%
  mutate(RP = RO / (AO + RO), AP = AO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

ggplot(fil1.df, aes(AP)) +
  geom_histogram()

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
fil1_name <- strsplit(as.character(fil1.df$Scaffold), split = "_")
fil1_name <- matrix(unlist(fil1_name), ncol = 3, byrow = TRUE)
fil1_name <- data.frame(fil1_name)
fil1.df$Scaffold <- fil1_name$X3

fil1.df2 <- left_join(fil1.df, map) %>%
  filter(LGs != "NA")

fil1.df3 <-  filter(fil1.df2, AO + RO > 100, AO + RO < 10000)

save(fil1.df3, file = "fil1.rda")
load("fil1.rda")
load("all.df.rda")

fil1.df4 <- inner_join(fil1.df3, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000)

dot.p <- ggplot(fil1.df4, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme(legend.position = "none")
dot.p

ggsave("allele_point.pdf", dot.p, width = 15, height = 10)

joy.p <- ggplot(fil1.df4, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme(legend.position = "none")
joy.p

hist.p <- ggplot(fil1.df4, aes(AP)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  facet_wrap(~LGs, scales = "free_y")
hist.p


##############################
######     FILTER 2     ######
##############################

fil2 <- read.table("../bcftools.snp.fil3.dp4.vcf", col.names = c("Scaffold", "Position", "Reference", "Alternative", "DP4"), sep = " ") %>%
  filter(Scaffold != "mitochondrion,") %>%
  separate(DP4, c("F_REF", "R_REF", "F_ALT", "R_ALT"), sep = ",")

# filter out ratios that wouldn't occur in octoploid SNPs
fil2.df <- fil2 %>% transform(F_REF = as.numeric(F_REF), R_REF = as.numeric(R_REF), F_ALT = as.numeric(F_ALT), R_ALT = as.numeric(R_ALT)) %>%
  mutate(RO = F_REF + R_REF, AO = F_ALT + R_ALT) %>%
  mutate(RP = RO / (AO + RO), AP = AO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

ggplot(fil2.df, aes(AP)) +
  geom_histogram()

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
fil2_name <- strsplit(as.character(fil2.df$Scaffold), split = "_")
fil2_name <- matrix(unlist(fil2_name), ncol = 3, byrow = TRUE)
fil2_name <- data.frame(fil2_name)
fil2.df$Scaffold <- fil2_name$X3

fil2.df2 <- left_join(fil2.df, map) %>%
  filter(LGs != "NA")

fil2.df3 <-  filter(fil2.df2, AO + RO > 200, AO + RO < 10000)

save(fil2.df3, file = "fil2.rda")
load("fil2.rda")
load("all.df.rda")

fil2.df4 <- inner_join(fil2.df3, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000)

dot.p <- ggplot(fil2.df4, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme(legend.position = "none")
dot.p

ggsave("allele_point.pdf", dot.p, width = 15, height = 10)

joy.p <- ggplot(fil2.df4, aes(x = AP, y = LGs)) +
  geom_joy(aes(colour = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme(legend.position = "none")
joy.p

hist.p <- ggplot(fil2.df4, aes(AP)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  facet_wrap(~LGs, scales = "free_y")
hist.p

ggsave("allele_hist.pdf", hist.p, width = 15, height = 10)

