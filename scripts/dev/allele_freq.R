library(tidyverse)
library(ggridges)
library(viridis)
library(cowplot)

setwd("~/Box Sync/GHdata/BgeVars/freq/")

map <- read.table("../scaffold_mapping.txt", sep = " ", header = TRUE)
load("all.df.rda")

##################################
######   VARIANT ANALYSIS   ######
##################################

fil2 <- read.table("../bcftools.snp.fil3.dp4.vcf", col.names = c("Scaffold", "Position", "Reference", "Alternative", "DP4"), sep = " ") %>%
  filter(Scaffold != "mitochondrion,") %>%
  separate(DP4, c("F_REF", "R_REF", "F_ALT", "R_ALT"), sep = ",")

# remove ratios that wouldn't occur in octoploid SNPs
fil2.df <- fil2 %>% transform(F_REF = as.numeric(F_REF), R_REF = as.numeric(R_REF), F_ALT = as.numeric(F_ALT), R_ALT = as.numeric(R_ALT)) %>%
  mutate(RO = F_REF + R_REF, AO = F_ALT + R_ALT) %>%
  mutate(RP = RO / (AO + RO), AP = AO / (AO + RO)) %>%
  filter(AP > 0.125, AP < 0.875)

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
fil2_name <- strsplit(as.character(fil2.df$Scaffold), split = "_")
fil2_name <- matrix(unlist(fil2_name), ncol = 3, byrow = TRUE)
fil2_name <- data.frame(fil2_name)
fil2.df$Scaffold <- fil2_name$X3

fil2.df2 <- left_join(fil2.df, map) %>%
  filter(LGs != "NA")

# remove sites with lower than 50 and greater than 10000 reads
fil2.df3 <-  filter(fil2.df2, AO + RO > 50, AO + RO < 10000)

save(fil2.df3, file = "fil2.rda")
load("fil2.rda")

# remove variants on scaffolds less than 10000 in length
fil2.df4 <- inner_join(fil2.df3, all.df) %>%
  select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000) %>%
  filter(!Scaffold %in% c("Scaffold33", "Scaffold24", "Scaffold4970")) #remove scaffolds with greatly inflated SNV counts

lg_lengths <- all.df %>%
  group_by(LGs) %>%
  summarise(Length = sum(Scaffold_Length))

lg_summary <- fil2.df4 %>%
  group_by(LGs) %>%
  summarise(SNVs = n()) %>%
  left_join(., lg_lengths) %>%
  mutate(SNV_per_bp = SNVs/Length) %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other"))

lg_plot <- ggplot(lg_summary) +
  geom_point(aes(x = LGs, y = SNV_per_bp, group = LGs), size = 5) +
  labs(x = "LGs, arranged from longest to shortest", y = "SNV per base pair", x = "") + 
  coord_flip() +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none")
lg_plot

scaffold_summary <- fil2.df4 %>%
  group_by(Scaffold) %>%
  summarise(SNVs = n()) %>%
  left_join(., all.df) %>%
  select(Scaffold, SNVs, Scaffold_Length) %>%
  mutate(SNV_per_bp = SNVs/Scaffold_Length)

scaffold_plot <- ggplot(scaffold_summary) +
  geom_point(aes(x = Scaffold, y = SNV_per_bp), size = 5) +
  labs(x = "LGs, arranged from longest to shortest", y = "SNV per base pair", x = "") + 
  # coord_flip() +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none")
scaffold_plot

scaffold_plot <- ggplot(scaffold_summary) +
  geom_point(aes(x = Scaffold_Length, y = SNVs), size = 5) +
  labs(x = "LGs, arranged from longest to shortest", y = "SNV per base pair", x = "") + 
  # coord_flip() +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none")
scaffold_plot

ggsave("SNV_frequency.pdf", summary_plot, width = 8, height = 6, units = "in")

df.p <- group_by(fil2.df4, LGs) %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other")) %>%
  sample_n(300)

dot.p <- ggplot(df.p, aes(x = LGs, y = AP, group = LGs)) +
  geom_jitter(aes(colour = LGs)) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Fraction", 
                     breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
dot.p

ggsave("allele_point.pdf", dot.p, width = 11, height = 8.5, units = "in")
ggsave("allele_point.png", dot.p, width = 11, height = 8.5, units = "in")

joy.p <- ggplot(df.p, aes(x = AP, y = LGs)) +
  geom_density_ridges(aes(fill = LGs)) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", 
                     breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none")
joy.p

hist.p <- ggplot(fil2.df4, aes(AP)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(name = "Density of Alternate Allele Fraction", 
                     breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) + 
  theme_bw(base_size = 12, base_family = "Helvetica") +
  facet_wrap(~LGs, scales = "free_y")
hist.p

ggsave("allele_hist.pdf", hist.p, width = 15, height = 10)

final_plot <- plot_grid(dot.p, joy.p, labels = c("A", "B"), rel_widths = c(1,2), ncol = 1)








