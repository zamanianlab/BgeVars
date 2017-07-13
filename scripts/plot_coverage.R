library(tidyverse)
library(data.table)


# read in bedgraph files
snail.file = "~/data/Bgla/BAM/snail.bga.bedgraph"
bge.file = "~/data/Bge/BGE.bga.bedgraph"

snail.data <- fread(snail.file)
colnames(snail.data) <- c("Contig", "Start", "Stop", "Coverage")
head(snail.data)
bge.data <- fread(bge.file)
colnames(bge.data) <- c("Contig", "Start", "Stop", "Coverage")
head(bge.data)

# summarize by taking a weighted average of coverage across the contig and getting the contig length
snail.summ.data <- snail.data %>%
  mutate(L_len = Start-Stop+1) %>%
  mutate(C_weight = Coverage*L_len) %>%
  group_by(Contig) %>%
  summarise(C_ave = sum(C_weight)/sum(L_len), C_length = max(Stop)) %>%
  arrange(desc(C_length))

bge.summ.data <- bge.data %>%
  mutate(L_len = Start-Stop+1) %>%
  mutate(C_weight = Coverage*L_len) %>%
  group_by(Contig) %>%
  summarise(C_ave = sum(C_weight)/sum(L_len), C_length = max(Stop)) %>%
  arrange(desc(C_length))

write.table(bge.summ.data, "~/data/Bge/bge.summ.data.csv", sep = "\t")
write.table(snail.summ.data, "~/data/Bgla/BAM/snail.summ.data.csv", sep = "\t")

################################################################################

# plotting

setwd("~/Box Sync/ZamanianLab/LabMembers/Nic/temp/")
snail.summ.data <- read.table("~/Box Sync/ZamanianLab/LabMembers/Nic/temp/snail.summ.data.csv", sep = "\t")
bge.summ.data <- read.table("~/Box Sync/ZamanianLab/LabMembers/Nic/temp/bge.summ.data.csv", sep = "\t")

# arrange by descending coverage
snail.df <- snail.summ.data #%>%
  #arrange(desc(C_ave))
bge.df <- bge.summ.data %>%
  arrange(desc(C_ave))

# rename columns to prepare for merging
colnames(snail.df)[2] <- "S.C_ave"
colnames(bge.df)[2] <- "B.C_ave"
#merge data, reorder columns, filter based on contig size, and normalize coverage measurements
merge.df <- left_join(bge.df, snail.df)
merge.df <- merge.df[c("Contig", "C_length", "S.C_ave", "B.C_ave")] %>%
  filter(C_length > 10000) %>%
  mutate(S.C_ave_norm = S.C_ave * sum(B.C_ave) / sum(S.C_ave)) %>%
  mutate(B.C_ave_norm = B.C_ave * sum(S.C_ave) / sum(B.C_ave))

# add index row  
merge.df <-  mutate(merge.df, N = seq.int(nrow(merge.df)))
# melt for plotting
m.m <- select(merge.df, N, Contig, C_length, S.C_ave_norm, B.C_ave_norm) %>%
  melt(id = c("N", "Contig", "C_length"))



b.p <- ggplot(data = m.m, aes(x = N, y = log10(value+1))) +
  geom_point(aes(color = variable), size=0.1, alpha=0.2) +
  theme_bw() +
  xlab("Contigs > 10kb, Arranged by Bge Genome Coverage") +
  ylab("Coverage") +
  scale_x_continuous(expand = c(0,0)) +
  ggtitle("Coverage Across Reference Contigs")
b.p


ggsave("Bge_snail_coverage.pdf", b.p, width = 12, height = 12)

################################################################################

# probing the data

t.df <- select(merge.df, -N) %>%
  mutate(Diff = B.C_ave_norm - S.C_ave_norm) %>%
  mutate(Group = rep("group", length(Contig))) %>%
  arrange(desc(C_length)) %>%
  mutate(N = seq.int(nrow(t.df)))


t.p <- ggplot(data = t.df, aes(x = N, y = Diff)) +
  geom_point(size=0.1, alpha=0.2) +
  geom_boxplot(aes(group = "Group"), outlier.color = "red") +
  theme_bw() + 
  xlab("Contigs > 10kb, Arranged by Length") +
  ylab("Bge coverage - Snail coverage")
t.p










