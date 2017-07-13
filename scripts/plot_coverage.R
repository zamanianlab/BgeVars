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

# snail.summ.data <- fread("~/data/Bgla/BAM/snail.summ.data.csv") %>%
#   arrange(desc(C_length))
# bge.summ.data <- fread("~/data/Bge/bge.summ.data.csv")
# 
# snail.summ.data <- read.table("~/Box Sync/ZamanianLab/LabMembers/Nic/temp/snail.summ.data.csv", sep = "\t") %>%
#   arrange(desc(C_length))
# bge.summ.data <- read.table("~/Box Sync/ZamanianLab/LabMembers/Nic/temp/bge.summ.data.csv", sep = "\t")
# 
# colnames(snail.summ.data) <- c("Contig", "C_ave", "C_length")
# colnames(bge.summ.data) <- c("Contig", "C_ave", "C_length")
# 
# snail.summ.data$N <- seq.int(nrow(bge.summ.data))
# 
# b.p <- ggplot(data = bge.summ.data, aes(x = N, y = C_ave)) +
#   stat_bin_hex(bins = 500) +
#   theme_bw() +
#   xlab("Contigs Arranged by Size") +
#   ylab("Coverage") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(limits = c(0, 5000)) +
#   ggtitle("Coverage Across Reference")
# 
# ggsave("mapping_all.pdf", p, width = 12, height = 12)



