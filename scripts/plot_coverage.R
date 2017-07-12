library(tidyverse)
library(data.table)

inputFile = "~/data/Bge/BGE.bga.bedgraph"

all.data <- fread(inputFile)

# filtered.data <- all.data %>%
#   mutate(C_weight = V4*(V3-V2+1))
  
# summ.data <- filtered.data %>%  
#   group_by(V1) %>%
#   #mutate(C_length = max(V3)) %>%
#   summarise(C_ave = mean(C_weight), C_length = max(V3)) %>%
#   arrange(desc(C_ave))
  
summ.data <- all.data %>%
  mutate(L_len = V3-V2+1) %>%
  mutate(C_weight = V4*L_len) %>%
  group_by(V1) %>%
  summarise(C_ave = sum(C_weight)/sum(L_len), C_length = max(V3)) %>%
  arrange(desc(C_length))

write.table(summ.data, "~/data/Bge/summ_data.csv", sep = "\t")

summ.data <- fread("~/Box Sync/GHdata/Local_BgeVars/summ_data.csv") %>%
  arrange(desc(C_length))

summ.data$N <- seq.int(nrow(summ.data))

p <- ggplot(data = summ.data, aes(x = N, y = C_ave)) +
  stat_bin_hex(bins = 500) + 
  theme_bw() +
  xlab("Contigs Arranged by Size") +
  ylab("Coverage") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5000)) +
  ggtitle("Coverage Across Reference")

ggsave("mapping_all.pdf", p, width = 12, height = 12)


