library(tidyverse)
library(ggjoy)
library(data.table)
library(stats)

setwd("~/Box Sync/GHdata/BgeVars/")

# read in bedgraph files
snail.file = "snail.bga.bedgraph"
bge.file = "CA301ANXX.bga.bedgraph"
snail.data <- fread(snail.file)
colnames(snail.data) <- c("Scaffold", "Start", "Stop", "Coverage")
bge.data <- fread(bge.file)
colnames(bge.data) <- c("Scaffold", "Start", "Stop", "Coverage")


# Locus_Coverage is the coverage weighted by length 
# The average coverage for the entire contig is the sum Locus_Coverage divided by the Contig_Length
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

save(snail.summ.data, file = "snail.summ.data.rda")
save(bge.summ.data, file = "bge.summ.data.rda")

################################################################################

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
# Add new group, "Ambiguous", as those contigs that are in multiple LGs
all.df$LGs <- gsub(".*,.*", "Ambiguous", all.df$LGs)

save(all.df, file = "all.df.rda")

load("all.df.rda")

# normalize coverage measurements
scale.factor <- (sum(all.df$Snail_Scaffold_Coverage) + sum(all.df$Bge_Scaffold_Coverage) / 2)
all.df.norm <- select(all.df, Scaffold, Scaffold_Length, Snail_Scaffold_Coverage, Bge_Scaffold_Coverage, LGs) %>%
  mutate(Snail_Scaffold_Coverage_Norm = Snail_Scaffold_Coverage * Scaffold_Length / sum(Snail_Scaffold_Coverage*Scaffold_Length)) %>%
  mutate(Bge_Scaffold_Coverage_Norm = Bge_Scaffold_Coverage * Scaffold_Length / sum(Bge_Scaffold_Coverage*Scaffold_Length))

s.test <- sum(all.df.rescale$Snail_Scaffold_Coverage_Rescale)
b.test <- sum(all.df.rescale$Bge_Scaffold_Coverage_Rescale)

# get percentage of contigs and length lost when filtering by C_length
filter.df.norm <- filter(all.df.norm, Scaffold_Length > 10000)
percent_filtered <- nrow(filter.df.norm) / nrow(all.df.norm)

# add index row and remove unwanted LGs
unwanted <- c("Ambiguous", "Other", "LGx", "LGw", "LGv", "LGu", "LGt", "LGs")
plot.df <- filter.df.norm %>%
  filter(!LGs %in% unwanted) %>%
  arrange(LGs, desc(Bge_Scaffold_Coverage_Norm))
plot.dt <- data.table(plot.df)
plot.dt[, agg.Scaffold_Length := cumsum(Scaffold_Length)]
plot.df <- data.frame(plot.dt)

# melt for plotting
plot.m <- select(plot.df, Scaffold, Scaffold_Length, agg.Scaffold_Length, LGs, Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm) %>%
  melt(id = c("agg.Scaffold_Length", "Scaffold", "Scaffold_Length", "LGs"))

# plot
b.p <- ggplot(data = plot.m, aes(x = agg.Scaffold_Length, y = scales::rescale(value, to = c(0,1)))) +
  geom_step(aes(color = variable)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Scaffolds > 10kb, Arranged by Bge Genome Coverage") +
  ylab("Normalized Coverage") +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_discrete(name = "Genome",
                      labels = c("Biomphalaria glabrata", "Bge cell line - Yoshino")) +
  facet_wrap(~LGs, scales = "free_x", nrow = 3)
b.p

ggsave("coverage_normalize.png", b.p, width = 15, height = 12)
ggsave("coverage_normalize.pdf", b.p, width = 15, height = 12)

################################################################################

#### probing the data ####

# histograms of scaffold coverages, grouped by linkage group
hist.df <- select(plot.df, Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm, LGs)
hist.m <- melt(hist.df, id ="LGs")
hist.p <- ggplot(data = hist.m, aes(x = log2(value), y = as.factor(LGs), fill = variable)) +
  geom_joy(alpha = 0.3, scale = 0.9) +
  scale_x_continuous(limits = c(-20, -10))
hist.p

ggsave("total_coverage_distributions.png", all.p, width = 12, height = 12)

# statistics
test.df <- select(plot.df, Scaffold, LGs, Snail_Scaffold_Coverage_Norm, Bge_Scaffold_Coverage_Norm) %>%
  gather(key = variable, value = value, -LGs, -Scaffold) %>%
  group_by(LGs, variable) %>%
  summarise(value = list(value)) %>%
  spread(variable, value) %>%
  group_by(LGs) %>%
  mutate(p_value = ks.test(unlist(Snail_Scaffold_Coverage_Norm), unlist(Bge_Scaffold_Coverage_Norm))$p.value,
         t_value = ks.test(unlist(Snail_Scaffold_Coverage_Norm), unlist(Bge_Scaffold_Coverage_Norm))$statistic)


hist.p <- hist.p + 
  annotate("text", size = 3, x = -10.5, y = 1.5:18.5, label = paste("p-value =", plyr::round_any(test.df$p_value, 0.001, ceiling))) +
  xlab("log2(Normalized Coverage)") +
  ylab("Density per Linkage Group") + 
  scale_fill_discrete(name = "Genome",
                      labels = c("Biomphalaria glabrata", "Bge cell line - Yoshino"))

ggsave("LG_coverage_distributions.png", hist.p, width = 12, height = 12)
ggsave("LG_coverage_distributions.pdf", hist.p, width = 12, height = 12)


all.m <- mutate(hist.m, N = rep("Group"))
all.p <- ggplot(data = all.m, aes(x = log2(value), fill = variable)) +
  geom_density(alpha = 0.3) +
  scale_x_continuous(limits = c(-20, -10)) +
  xlab("log2(Normalized Coverage)") +
  ylab("Density per Genome") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_fill_discrete(name = "Genome",
                      labels = c("Biomphalaria glabrata", "Bge cell line - Yoshino"))
all.p

ggsave("total_coverage_distributions.pdf", all.p, width = 6, height = 6)
ggsave("total_coverage_distributions.png", all.p, width = 6, height = 6)

t.df <- select(plot.df, Scaffold, Scaffold_Length, LGs:Bge_Scaffold_Coverage_Norm) %>%
  mutate(Difference = Bge_Scaffold_Coverage_Norm - Snail_Scaffold_Coverage_Norm) %>%
  arrange(desc(Scaffold_Length)) %>%
  mutate(N = seq.int(nrow(.))) %>%
  group_by(LGs)

quantiles <- t.df %>%
  group_by(LGs) %>% 
  summarise(lowerq = quantile(Difference)[2], upperq = quantile(Difference)[4]) %>%
  mutate(iqr = upperq - lowerq, mild.threshold.upper = (iqr * 1.5) + upperq, mild.threshold.lower = lowerq - (iqr * 1.5))

outliers <- t.df %>%
  select(-N) %>%
  group_by(LGs) %>%
  filter(Difference > ((quantile(Difference)[4] - quantile(Difference)[2]) * 1.5 + quantile(Difference)[4]) | Difference < (quantile(Difference)[2] - (quantile(Difference)[4] - quantile(Difference)[2]) * 1.5)) %>%
  mutate(Outlier = "Outlier")

t.df <- left_join(t.df, outliers)

t.p <- ggplot(data = t.df, aes(x = N, y = log(Difference))) +
  geom_point(size=0.1, alpha=0.2) +
  geom_boxplot(aes(group = LGs), outlier.color = "red") +
  theme_bw() + 
  xlab("Contigs > 10kb, Arranged by Length") +
  ylab("log(Bge coverage - Snail coverage)") +
  theme(axis.text.x=element_blank(),
             axis.ticks.x=element_blank()) +
  facet_wrap(~LGs, scales = "free_x", nrow = 3) +
  geom_text(data = subset(t.df, Outlier == "Outlier"), aes(label = Scaffold), na.rm = TRUE)
t.p



ggsave("diff.pdf", t.p, width = 12, height = 12)








