library(tidyverse)
library(data.table)

setwd("~/Box Sync/GHdata/Local_BgeVars/")

# read in bedgraph files
snail.file = "snail.bga.bedgraph"
bge.file = "CA301ANXX.bga.bedgraph"
snail.data <- fread(snail.file)
colnames(snail.data) <- c("Contig", "Start", "Stop", "Coverage")
bge.data <- fread(bge.file)
colnames(bge.data) <- c("Contig", "Start", "Stop", "Coverage")


# summarize by taking a weighted average of coverage across the contig and getting the contig length
snail.summ.data <- snail.data %>%
  mutate(C_weight = Coverage*(Stop-Start+1)) %>%
  group_by(Contig) %>%
  summarise(C_ave = sum(C_weight)/sum(Stop-Start+1), C_length = max(Stop)) %>%
  arrange(desc(C_length))

bge.summ.data <- bge.data %>%
  mutate(C_weight = Coverage*(Stop-Start+1)) %>%
  group_by(Contig) %>%
  summarise(C_ave = sum(C_weight)/sum(Stop-Start+1), C_length = max(Stop)) %>%
  arrange(desc(C_length))

save(snail.summ.data, file = "snail.summ.data.rda")
save(bge.summ.data, file = "bge.summ.data.rda")

################################################################################

load("snail.summ.data.rda")
load("bge.summ.data.rda")

snail.df <- snail.summ.data
bge.df <- bge.summ.data

# rename columns to prepare for merging
colnames(snail.df)[2] <- "S.C_ave"
colnames(bge.df)[2] <- "B.C_ave"
# merge data into a single data frame
merge.df <- left_join(bge.df, snail.df)
# add linkage group (LG) information as a new column, prepare by renaming Scaffolds to match with scaffold_mapping.txt
name <- strsplit(as.character(merge.df$Contig), split = "_")
name <- matrix(unlist(name), ncol = 3, byrow = TRUE)
name <- data.frame(name)
merge.df <- mutate(merge.df, new_Contig = as.factor(name$X3))
LG_mapping <- read.table("scaffold_mapping.txt", sep = " ", header = TRUE)
colnames(LG_mapping)[1] <- "new_Contig"
all.df <- left_join(merge.df, LG_mapping)
all.df[is.na(all.df)] <- "Other"
# Add new group, "Ambiguous", as those contigs that are in multiple LGs
all.df$LGs <- as.character(all.df$LGs)
all.df$LGs <- gsub(".*,.*", "Ambiguous", all.df$LGs)
all.df[is.na(all.df)] <- "Other"

save(all.df, file = "all.df.rda")

load("all.df.rda")

# filter based on contig size, and normalize coverage measurements
scale.factor <- (sum(all.df$S.C_ave) + sum(all.df$B.C_ave) / 2)
all.df.norm <- select(all.df, Contig, C_length, S.C_ave, B.C_ave, LGs) %>%
  filter(S.C_ave >= 0, B.C_ave >= 0) %>% #why are there NAs
  mutate(S.C_ave_norm = S.C_ave * C_length / sum(S.C_ave*C_length)) %>%
  mutate(B.C_ave_norm = B.C_ave * C_length / sum(B.C_ave*C_length))


b.p <- ggplot(data = all.df.norm) +
  geom_histogram(aes(x=B.C_ave_norm, color = LGs)) +
  geom_histogram(aes(x=S.C_ave_norm, color = LGs)) +
  theme_bw() +
  xlab("Contigs > 10kb, Arranged by Bge Genome Coverage") +
  ylab("Coverage") +
  scale_x_continuous(expand = c(0,0)) +
  facet_wrap( ~ LGs, scales = "free_y")
  ggtitle("Coverage Across Reference Contigs")# +
#facet_wrap(~m.m$LGs, scales = "free_x")
b.p


temp1 <- sum(all.df.norm$S.C_ave * all.df.norm$C_length)
temp2 <- sum(all.df.norm$B.C_ave * all.df.norm$C_length) 

mutate(S.C_ave_norm= S.C_ave_norm * (sum(S.C_ave*C_length) + sum(B.C_ave*C_length))/ 2) %>%
mutate(B.C_ave_norm= B.C_ave_norm * (sum(S.C_ave*C_length) + sum(B.C_ave*C_length))/ 2)

s.test2 <- sum(all.df.norm$S.C_ave_norm)
b.test2 <- sum(all.df.norm$B.C_ave_norm)


%>%
  filter(C_length > 10000) 
  



# get percentage of contigs and length lost when filtering by C_length

# add index row  
all.df <- all.df %>%
  arrange(LGs, desc(B.C_ave)) %>%
  mutate(N = seq.int(nrow(all.df)))
# melt for plotting
m.m <- select(all.df, N, Contig, C_length, LGs, S.C_ave_norm, B.C_ave_norm) %>%
  melt(id = c("N", "Contig", "C_length", "LGs"))



b.p <- ggplot(data = m.m, aes(x = N, y = log10(value+1))) +
  geom_point(aes(color = LGs, shape = variable), size=0.5, alpha=0.2) +
  theme_bw() +
  xlab("Contigs > 10kb, Arranged by Bge Genome Coverage") +
  ylab("Coverage") +
  scale_x_continuous(expand = c(0,0)) +
  ggtitle("Coverage Across Reference Contigs")# +
  #facet_wrap(~m.m$LGs, scales = "free_x")
b.p


ggsave("all_coverage.pdf", b.p, width = 12, height = 12)

################################################################################

# probing the data

t.df <- select(merge.df, -N) %>%
  mutate(Diff = B.C_ave_norm - S.C_ave_norm) %>%
  mutate(Group = rep("group", length(Contig))) %>%
  arrange(desc(C_length))


t.df <-  mutate(t.df, N = seq.int(nrow(t.df)))


t.p <- ggplot(data = t.df, aes(x = N, y = Diff)) +
  geom_point(size=0.1, alpha=0.2) +
  geom_boxplot(aes(group = "Group"), outlier.color = "red") +
  theme_bw() + 
  xlab("Contigs > 10kb, Arranged by Length") +
  ylab("Bge coverage - Snail coverage")
t.p

ggsave("diff.pdf", t.p, width = 12, height = 12)








