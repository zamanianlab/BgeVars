library(ggridges)
library(viridis)
library(cowplot)
library(data.table)
library(stats)
library(ggrepel)
library(tidyverse)

# Script for generating figures for Bge genome manuscript

#########################################################################################################
######################                                                             ######################
######################                        Figure 1                             ######################
######################                                                             ######################
#########################################################################################################

setwd("~/GitHub/BgeVars/auxillary/coverage")

# read in bedgraph files (Multi GB files, not included here. Instead, load .rda file below.)
snail.file = "snail.bga.bedgraph"
bge.file = "CA301ANXX.bga.bedgraph"
snail.data <- fread(snail.file)
colnames(snail.data) <- c("Scaffold", "Start", "Stop", "Coverage")
bge.data <- fread(bge.file)
colnames(bge.data) <- c("Scaffold", "Start", "Stop", "Coverage")

## summarize data
# Locus_Coverage is the coverage weighted by length 
# the average coverage for the entire contig is the sum of Locus_Coverage divided by the Contig_Length
snail.summ.data <- snail.data %>%
  mutate(Locus_Coverage = Coverage*(Stop-Start+1)) %>%
  group_by(Scaffold) %>%
  summarise(Scaffold_Coverage = sum(Locus_Coverage)/(max(Stop) + 1), Scaffold_Length = (max(Stop) + 1)) %>%
  arrange(desc(Scaffold_Length))

bge.summ.data <- bge.data %>%
  mutate(Locus_Coverage = Coverage*(Stop-Start+1)) %>%
  group_by(Scaffold) %>%
  summarise(Scaffold_Coverage = sum(Locus_Coverage)/max(Stop+1), Scaffold_Length = (max(Stop + 1))) %>%
  arrange(desc(Scaffold_Length))

load("snail.summ.data.rda")
load("bge.summ.data.rda")

snail.df <- snail.summ.data
bge.df <- bge.summ.data

# rename columns to prepare for merging
colnames(snail.df)[2] <- "Snail_Scaffold_Coverage"
colnames(bge.df)[2] <- "Bge_Scaffold_Coverage"

# merge data into a single data frame
merge.df <- left_join(bge.df, snail.df)

# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
name <- strsplit(as.character(merge.df$Scaffold), split = "_")
name <- matrix(unlist(name), ncol = 3, byrow = TRUE)
name <- data.frame(name)
merge.df$Scaffold <- name$X3
LG_mapping <- read.table("scaffold_mapping.txt", sep = " ", header = TRUE)
all.df <- left_join(merge.df, LG_mapping)

# Change NAs to "Other"
all.df$LGs <- as.character(all.df$LGs)
all.df$LGs[is.na(all.df$LGs)] <- "Other"

# Add new group, "Ambiguous", as those contigs that map to multiple LGs
all.df$LGs <- gsub(".*,.*", "Ambiguous", all.df$LGs)

# normalize coverage measurements
all.df.norm <- dplyr::select(all.df, Scaffold, Scaffold_Length, Snail_Scaffold_Coverage, Bge_Scaffold_Coverage, LGs) %>%
  mutate(Snail_Scaffold_Coverage_Norm = Snail_Scaffold_Coverage * Scaffold_Length / sum(Snail_Scaffold_Coverage*Scaffold_Length)) %>%
  mutate(Bge_Scaffold_Coverage_Norm = Bge_Scaffold_Coverage * Scaffold_Length / sum(Bge_Scaffold_Coverage*Scaffold_Length))
all.df.norm.m <- gather(all.df.norm, key = Coverage_Type, value = Value, Snail_Scaffold_Coverage, Bge_Scaffold_Coverage, Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm) %>%
  separate(Coverage_Type, into = c("ID", "Coverage_Type"), sep = "_Scaffold_Coverage")
dot.df <- dplyr::select(all.df.norm, -Bge_Scaffold_Coverage_Norm, -Snail_Scaffold_Coverage_Norm)

# remove small scaffolds and LGs that aren't the largest 18
dot.m <- dot.df %>%
  filter(Scaffold_Length > 10000) %>%
  arrange(LGs) %>% group_by(LGs) %>%
  mutate(Index = row_number()) %>%
  melt(id = c("Scaffold", "Scaffold_Length", "LGs", "Index")) %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other"))

# calculate summary statistics
dot.summ <- dot.m %>% 
  group_by(LGs, variable) %>%
  summarise(med = median(value), ave = sum(Scaffold_Length * value) / sum(Scaffold_Length)) %>%
  mutate(lm = log2(med))

# prepare for plotting
bge.dot.summ <- subset(dot.summ, variable == "Bge_Scaffold_Coverage")
linkage_order <- dot.summ %>%
  arrange(desc(ave)) %>%
  filter(variable == "Bge_Scaffold_Coverage")
linkage_order <- linkage_order$LGs
dot.m$LGs <- factor(dot.m$LGs, levels = c(linkage_order))

# plot
fig1 <- ggplot(data = subset(dot.m, variable == "Bge_Scaffold_Coverage"), aes(x = LGs, y = value, group = LGs)) +
  geom_jitter(aes(alpha = 0.6, color = "black")) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  geom_hline(data = subset(dot.summ, variable == "Bge_Scaffold_Coverage"),
             color = "black",
             aes(yintercept = ave)) +
  geom_label(data = bge.dot.summ, 
             size = 2.5, 
             aes(x = LGs, y = -5, label = plyr::round_any(ave, 0.001, ceiling))) +
  ylim(-10, 250) +
  labs(y = "Coverage") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# fig1

save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig1.pdf", fig1, base_width = 18, base_aspect_ratio = .5, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig1.eps", fig1, base_width = 18, base_aspect_ratio = .5, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig1.png", fig1, base_width = 18, base_aspect_ratio = .5, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig1.tiff", fig1, base_width = 18, base_aspect_ratio = .5, units = "in", dpi = 400)

#########################################################################################################
######################                                                             ######################
######################                        Figure 2                             ######################
######################                                                             ######################
#########################################################################################################

source("http://bioconductor.org/biocLite.R")
library("topGO")
library("xlsx")

setwd("~/Box Sync/GHdata/BgeVars/snpeff/")

# function for capitilizing text
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

##### Figure 2A #####

# read in output from snpeff
df <- read.csv("snpEff_genes.txt", sep = "\t", header = TRUE) %>%
  dplyr::select(-BioType, -GeneId)

# split into impact and effects
impact <- dplyr::select(df, TranscriptId, HIGH, LOW, MODERATE, MODIFIER)
effects <- dplyr::select(df, TranscriptId, X3_prime_UTR_variant:upstream_gene_variant)
effects.m <- reshape2::melt(effects)

# classify effects
# remove "frameshift_variant", "exon_loss_variant",  because trimmed variants to only SNVs in filtering strategy
high <- c("chromosome_number_variation", 
          "rare_amino_acid_variant", 
          "splice_acceptor_variant", 
          "splice_donor_variant",
          "start_lost",
          "stop_gained",
          "stop_lost",
          "transcript_ablation")

# remove "disruptive_inframe_deletion", "disruptive_inframe_insertion", "conservative_inframe_deletion", "conservative_inframe_insertion",
# because trimmed variants to only SNVs in filtering strategy
moderate <- c("coding_sequence_variant",
              "missense_variant",
              "regulatory_region_ablation",
              "splice_region_variant",
              "TFBS_ablation")

low <- c("X5_prime_UTR_premature_start_codon_gain_variant",
         "initiator_codon_variant",
         "start_retained",
         "stop_retained_variant",
         "synonymous_variant",
         "non_canonical_start_codon")

modifier <- c("X3_prime_UTR_variant",
              "X5_prime_UTR_variant",
              "coding_sequence_variant",
              "conserved_intergenic_variant",
              "conserved_intron_variant",
              "downstream_gene_variant",
              "exon_variant",
              "feature_elongation",
              "feature_truncation",
              "gene_variant",
              "intergenic_region",
              "intragenic_variant",
              "intron_variant",
              "mature_miRNA_variant",
              "miRNA",
              "NMD_transcript_variant",
              "non_coding_transcript_exon_variant",
              "non_coding_transcript_variant",
              "regulatory_region_amplification",
              "regulatory_region_variant",
              "TF_binding_site_variant",
              "TFBS_amplification",
              "transcript_amplification",
              "transcript_variant",
              "upstream_gene_variant")

high.df <- data.frame(x = high, y = rep("HIGH", length(high)))
moderate.df <- data.frame(x = moderate, y = rep("MODERATE", length(moderate)))
low.df <- data.frame(x = low, y = rep("LOW", length(low)))
modifier.df <- data.frame(x = modifier, y = rep("MODIFIER", length(modifier)))

impact_factors <- rbind(high.df, moderate.df, low.df, modifier.df) %>%
  dplyr::rename(variable = x)

effects.p <- left_join(effects.m, impact_factors) %>%
  dplyr::rename(Type = variable, Impact = y)

impact_factors <- dplyr::rename(impact_factors, Type = variable, Impact = y)

# count the number of genes that have eact SNV effect and classify by impact; cutoff number of genes at 1000 to help visualization
effects.temp <- effects.p
effects.temp$value[effects.temp$value >= 1] <- 1
effects.sum <-  group_by(effects.temp, Type) %>%
  summarise(Sum = sum(value)) %>%
  left_join(., impact_factors)
effects.sum$Sum[effects.sum$Sum >= 2000] <- 2000

# plot
fig2a <- ggplot(effects.sum) + geom_point(aes(x = Type, y = Sum, color = Impact), size = 5) +
  scale_y_continuous(limits = c(0, 2000), breaks = c(0, 1000, 2000), labels = c("0", "1000", ">2000")) +
  scale_x_discrete(limits=c("start_lost", "stop_gained", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "missense_variant",  "splice_region_variant", "X5_prime_UTR_premature_start_codon_gain_variant", "synonymous_variant", "stop_retained_variant", "non_canonical_start_codon", "initiator_codon_variant", "X5_prime_UTR_variant", "X3_prime_UTR_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant", "intron_variant", "downstream_gene_variant"),
                   labels = c("Start Lost", "Stop Gained", "Stop Lost", "Splice Acceptor Variant", "Splice Donor Variant", "Missense", "Splice Region Variant", "5' UTR Premature Start Codon", "Synonymous Mutation", "Stop Retained", "Non-Canonical Start Codon", "Initiator Codon Variant", "5' UTR Variant", "3' UTR Variant", "Upstream Gene Variant", "Non-Coding Transcript Variant", "Intron Variant", "Downstream Gene Variant")) +
  ylab("Number of Genes") +
  xlab("Type of Variant") +
  coord_flip() +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5))
fig2a

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig2A.pdf", fig2a, height = 10, width = 8)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig2A.png", fig2a, height = 10, width = 8)

##### Figure 2B #####

# import B. glabrata genes mapped to GO terms and named factor of all B. glabrata genes with HIGH impact genes marked 1
geneID2GO <- readMappings("geneID2GO")
geneNames <- names(geneID2GO)
allGenes <- read.csv("HIGH_genes.txt", sep = "\t", header = FALSE)
highGenes <- as.factor(allGenes$V2)
names(highGenes) <- allGenes$V1

# create topGOdata class for enrichment analyses
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = highGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = highGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = highGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)

# run Fisher enrichment
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
Fisher_CC <- getSigGroups(GOdata_CC, test.stat)
Fisher_BP <- getSigGroups(GOdata_BP, test.stat)
Fisher_MF <- getSigGroups(GOdata_MF, test.stat)

# organize results in table
Results_CC <- GenTable(GOdata_CC, classic = Fisher_CC, ranksOf = "classic", topNodes = 10) %>%
  mutate(Class = "Cellular Component")
Results_BP <- GenTable(GOdata_BP, classic = Fisher_BP, ranksOf = "classic", topNodes = 10) %>%
  mutate(Class = "Biological Process")
Results_MF <- GenTable(GOdata_MF, classic = Fisher_MF, ranksOf = "classic", topNodes = 10) %>%
  mutate(Class = "Molecular Function")
df2 <- rbind(Results_CC, Results_BP, Results_MF)
df2$Term <- capwords(df2$Term)

fig2b <- ggplot(df2) + 
  geom_point(aes(x = Term, y = -log10(as.numeric(classic))), color = "black", size = 5) +
  facet_grid(Class ~ ., scales = "free") +
  labs(y = "-log10(p-value)", x = "GO Term") + 
  coord_flip() +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
fig2b

ggplot2::ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig2B.pdf", fig2b, width = 11, height = 8.5, units = "in")
ggplot2::ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig2B.png", fig2b, width = 11, height = 8.5, units = "in")

##### Combine Figure 2 #####

fig2 <- plot_grid(fig2a, fig2b, labels = c("a", "b"), nrow = 1, rel_widths = c(1, 1.2))
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig2.pdf", fig2, base_width = 15, base_height = 6.5)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig2.eps", fig2, base_width = 15, base_height = 6.5)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig2.png", fig2, base_width = 15, base_height = 6.5)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig2.tiff", fig2, base_width = 15, base_height = 6.5, dpi = 400)

#########################################################################################################
######################                                                             ######################
######################                        Figure 3                             ######################
######################                                                             ######################
#########################################################################################################

setwd("~/Box Sync/GHdata/BgeVars/coverage")

# get columns for overlaying distributions of LG coverage
hist.df <- filter(all.df.norm, Scaffold_Length > 10000) %>%
  dplyr::select(Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm, LGs)
hist.m <- melt(hist.df, id ="LGs") %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other"))

##### Figure 3A #####

# plot total coverage
all.m <- mutate(hist.m, N = rep("Group"))

fig3a <- ggplot(data = all.m, aes(x = log2(value), fill = variable)) +
  geom_density() +
  scale_x_continuous(limits = c(-20, -10)) +
  labs(x = "log2(Normalized Coverage)", y = "Density") +
  scale_fill_viridis(discrete = TRUE, 
                     option = "viridis",
                     direction = -1,
                     name = "Genome",
                     labels = c(expression(italic("B. glabrata")), "Bge3 Cell Line")) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")
fig3a

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig3A.pdf", fig3a, width = 6, height = 6)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig3A.png", fig3a, width = 6, height = 6)

##### Figure 3B #####

# perform KS test
test.df <- filter(all.df.norm, Scaffold_Length > 10000) %>%
  dplyr::select(Scaffold, LGs, Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm) %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other")) %>%
  gather(key = variable, value = value, -LGs, -Scaffold) %>%
  group_by(LGs, variable) %>%
  summarise(value = list(value)) %>%
  spread(variable, value) %>%
  group_by(LGs) %>%
  mutate(p_value = ks.test(unlist(Snail_Scaffold_Coverage_Norm), unlist(Bge_Scaffold_Coverage_Norm))$p.value,
         t_value = ks.test(unlist(Snail_Scaffold_Coverage_Norm), unlist(Bge_Scaffold_Coverage_Norm))$statistic) 

# test.df$p_value <- as.factor(test.df$_p)

# plot coverage grouped by LG
fig3b <- ggplot(data = hist.m, aes(x = log2(value), y = as.factor(LGs), fill = variable)) +
  geom_density_ridges(scale = 0.9) +
  annotate("text", size = 3, x = -10.5, y = 1.5:18.5, label = paste("p =", plyr::round_any(test.df$p_value, 0.01, ceiling)), color = ifelse(test.df$p_value > 0.05, "black", "red")) +
  scale_x_continuous(limits = c(-20, -10)) +
  labs(x = "log2(Normalized Coverage)", y = "") +
  scale_fill_viridis(discrete = TRUE, 
                     option = "viridis", 
                     direction = -1,
                     name = "Genome",
                     labels = c(expression(italic("B. glabrata")), "Bge3 Cell Line")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "bottom",
        legend.title=element_blank())
fig3b

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig3B.pdf", fig3b, width = 8.5, height = 11)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig3B.png", fig3b, width = 8.5, height = 11)

##### Figure 3C #####

setwd("~/Box Sync/GHdata/BgeVars/freq/")

map <- read.table("../scaffold_mapping.txt", sep = " ", header = TRUE)
load("all.df.rda")

# variant filtration (large vcf files are not included, please load .rda file below)
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

load("fil2.rda")

# remove variants on scaffolds less than 10000 in length
fil2.df4 <- inner_join(fil2.df3, all.df) %>%
  dplyr::select(Scaffold, Scaffold_Length, LGs, AP) %>%
  filter(Scaffold_Length > 10000) %>%
  filter(!Scaffold %in% c("Scaffold33", "Scaffold24", "Scaffold4970")) #remove scaffolds with greatly inflated SNV counts

df.p <- group_by(fil2.df4, LGs) %>%
  filter(!LGs %in% c("LGs", "LGt", "LGu", "LGv", "LGw", "LGx", "Ambiguous", "Other")) #%>%
  # sample_n(300)

fig3c <- ggplot(df.p, aes(x = LGs, y = AP, group = LGs)) +
  geom_violin(colour = "black", fill = "grey40") +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  scale_y_continuous(name = "Alternate Allele Frequency", 
                     breaks = c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
fig3c

##### Combine Figure 3 #####

top_row <- plot_grid(fig3a, fig3b, labels = c("a", ""), rel_widths = c(1.4,3))
fig3 <- plot_grid(top_row, NULL, fig3c, labels = c("", "", "b"), ncol = 1, rel_heights = c(1.4, 0.05, 1))

save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig3.pdf", fig3, base_width = 10.5, base_height = 9, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig3.eps", fig3, base_width = 10.5, base_height = 9, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig3.png", fig3, base_width = 10.5, base_height = 9, units = "in")
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig3.tiff", fig3, base_width = 10.5, base_height = 9, dpi = 400)

#########################################################################################################
######################                                                             ######################
######################                        Figure 4                             ######################
######################                                                             ######################
#########################################################################################################

library(magick)

setwd("~/Box Sync/GHdata/BgeVars/karyotype/")

# manually populate karyotype data
Cell <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25)
Count <- c(62, 60, 65, 62, 63, 62, 62, 62, 59, 58, 62, 57, 63, 64, 62, 59, 61, 67, 60)
A <- c(7, 4, 10, 7, 11, 2, 8, 12, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
B <- c(14, 16, 12, 11, 18, 8, 13, 11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
C <- c(9, 7, 10, 10, 5, 8, 11, 2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
D <- c(11, 8, 7, 13, 8, 15, 11, 10, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
E <- c(8, 6, 6, 6, 5, 9, 5, 9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
F <- c(7, 15, 15, 9, 12, 15, 10, 9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
Unassigned <- c(6, 4, 5, 5, 4, 5, 5, 9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# function for calculating mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##### Figure 4A #####

# create dataframe
df3 <- data_frame(Cell, Count, A, B, C, D, E, F, Unassigned)

# plot
fig4a <- ggplot(df3, aes(Cell, Count)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6, alpha = 0.8) +
  # coord_flip() +
  geom_label(size = 3, 
             x = 13, 
             y = 69.75, 
             label = paste("Bge3 Mode = ", Mode(Count))) +
  geom_label(size = 3, 
             x = 13, 
             y = 68.75, 
             label = paste("Bge1 Mode = 63")) +
  geom_label(size = 3, 
             x = 13, 
             y = 67.75, 
             label = paste("Bge2 Mode = 67")) +
  scale_y_continuous(limits = c(55, 70)) +
  labs(x = "", y = "Chromosome Count") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
fig4a

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig4A.pdf", fig4a, height = 15, width = 5)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig4A.png", fig4a, height = 15, width = 5)

##### Figure 4B #####

# get group data
df4 <- df3 %>% 
  dplyr::select(-Cell, -Count) %>%
  gather(Group, Count, na.rm = TRUE) %>%
  group_by(Group)

df4$Group <- as.factor(df4$Group)

# calculate mean and import mean from Odoemelam 2009
df5 <- summarise(df4, mean(Count)) %>%
  dplyr::rename("Bge3" = "mean(Count)") %>%
  mutate(Bge1 = c(6.18, 10.36, 6.24, 19.42, 13.76, 2.8, 4.04), Bge2 = c(6.4, 11.04, 6.8, 20.58, 14.54, 2.14, 3.7)) %>%
  gather(Strain, Mean, -Group)

Bge1 <- subset(df5, Strain %in% "Bge1")
Bge2 <- subset(df5, Strain %in% "Bge2")
Bge3 <- subset(df5, Strain %in% "Bge3")

fig4b <- ggplot(df4, aes(Group, Count), group = Group, colour = group) +
  geom_boxplot(alpha = 0.7) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6, alpha = 0.8) +
  geom_label(size = 3,
    data = Bge3, 
    fill = "white",
    colour = "black",
    aes(x = Group, 
        y = 23.75, 
        label = paste(Strain, "Mean =", plyr::round_any(Mean, 0.01, ceiling)))) +
  geom_label(size = 3, 
    data = Bge1, 
    fill = "white",
    colour = "black",
    aes(x = Group, 
        y = 22, 
        label = paste(Strain, "Mean =", plyr::round_any(Mean, 0.01, ceiling)))) +
  geom_label(size = 3, 
    data = Bge2, 
    fill = "white",
    colour = "black",
    aes(x = Group, 
        y = 20.25, 
        label = paste(Strain, "Mean =", plyr::round_any(Mean, 0.01, ceiling)))) +
  ylim(0, 25) +
  labs(y = "") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none")
fig4b

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig4B.pdf", fig4b, height = 10, width = 14)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/Fig4B.png", fig4b, height = 10, width = 14)

##### Combine Figure 4 #####

fig4c <- "~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/Figures/Fig4C.pdf"
fig4c <- image_read_pdf(fig4c)
fig4c <- ggdraw() + draw_image(fig4c)

fig4d <- "~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/Figures/Fig4D.pdf"
fig4d <- image_read_pdf(fig4d)
fig4d <- ggdraw() + draw_image(fig4d)

top_row <- plot_grid(fig4a, fig4b, labels = c("a", "b"), rel_widths = c(1,4))
bottom_row <- plot_grid(NULL, NULL, labels = c("c", "d"), rel_widths = c(2,2), label_size = 15)
fig4 <- plot_grid(top_row, NULL, bottom_row, labels = c("", "", ""), ncol = 1, rel_heights = c(1.5, 0.1, 1), align = 'v', axis ='l')

save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig4_edit.pdf", fig4, base_height = 10, base_width = 12)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig4_edit.eps", fig4, base_height = 10, base_width = 12)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig4_edit.png", fig4, base_height = 10, base_width = 12)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Figures/Fig4_edit.tiff", fig4, base_height = 10, base_width = 12, units = "in", dpi = 400)

#########################################################################################################
######################                                                             ######################
######################                  Supplementary Figure 2                     ######################
######################                                                             ######################
#########################################################################################################

setwd("~/Box Sync/GHdata/BgeVars/coverage")

snail.dot.summ <- subset(dot.summ, variable == "Snail_Scaffold_Coverage")

# plot
sfig2 <- ggplot(data = subset(dot.m, variable == "Snail_Scaffold_Coverage"), aes(x = LGs, y = value, group = LGs)) +
  geom_jitter(aes(alpha = 0.6, color = "black")) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  geom_hline(data = subset(dot.summ, variable == "Snail_Scaffold_Coverage"), 
             color = "black",
             aes(yintercept = ave)) + 
  geom_label(data = snail.dot.summ, 
             size = 2.5, 
             aes(x = LGs, y = -5, label = paste(plyr::round_any(ave, 0.001, ceiling)))) +
  ylim(c(-10, 40)) +
  labs(y = "Coverage") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
sfig2

save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Additional_Files/Additional_file_2.pdf", sfig2, base_width = 18, base_aspect_ratio = .5, units = "in")


#########################################################################################################
######################                                                             ######################
######################                  Supplementary Figure 3                     ######################
######################                                                             ######################
#########################################################################################################

setwd("~/Box Sync/GHdata/BgeVars/coverage")

sfig3a <- ggplot(data = subset(dot.m, variable == "Bge_Scaffold_Coverage"), aes(x = LGs, y = log2(value), group = LGs)) +
  geom_jitter(aes(alpha = 0.6, colour = "black")) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  geom_hline(data = subset(dot.summ, variable == "Bge_Scaffold_Coverage"),
             color = "black",
             aes(yintercept = log2(ave))) +
  geom_label(data = bge.dot.summ, 
             size = 2.5, 
             aes(x = LGs, y = -1, label = paste(plyr::round_any(log2(ave), 0.001, ceiling)))) +
  ylim(-2, 10) +
  labs(x= "Bge3", y = "log2(Coverage)") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
sfig3a

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/SFig3A.pdf", sfig3a, width = 18, height = 6)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/SFig3A.png", sfig3a, width = 18, height = 6)

sfig3b <- ggplot(data = subset(dot.m, variable == "Snail_Scaffold_Coverage"), aes(x = LGs, y = log2(value), group = LGs)) +
  geom_jitter(aes(alpha = 0.6, colour = "black")) +
  facet_wrap(~ LGs, scale = "free_x", nrow = 1) +
  geom_hline(data = subset(dot.summ, variable == "Snail_Scaffold_Coverage"), 
             color = "black",
             aes(yintercept = log2(ave))) + 
  geom_label(data = snail.dot.summ, 
             size = 2.5, 
             aes(x = LGs, y = -1, label = paste(plyr::round_any(log2(ave), 0.001, ceiling)))) +
  ylim(c(-2, 10)) +
  labs(x = expression(italic("B. glabrata")), y = "log2(Coverage)") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
sfig3b

ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/SFig3B.pdf", sfig3b, width = 18, height = 6)
ggsave("~/Box Sync/ZamanianLab/Manuscripts/2018_Collaborations/2018-Bge_cells/Figures/SFig3B.png", sfig3b, width = 18, height = 6)

sfig3 <- plot_grid(sfig3a, sfig3b, labels = c("a", "b"), nrow = 2)
save_plot("~/Box Sync/ZamanianLab/Manuscripts/2018-Collaborations/2018-Bge_cells/BMC/Additional_Files/Additional_file_3.pdf", sfig3, base_width = 18, base_aspect_ratio = 0.5)











