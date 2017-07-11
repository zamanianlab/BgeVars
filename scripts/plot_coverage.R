library(tidyverse)

inputFile = "~/data/Bge/BGE.bga.bedgraph"

all.data <- read.csv(inputFile, header = FALSE, sep = "\t")

filtered.data <- all.data %>%
  mutate(C_weight = V4*(V3-V2+1))
  
summ.data <- filtered.data %>%  
  group_by(V1) %>%
  #mutate(C_length = max(V3)) %>%
  summarise(C_ave = mean(C_weight), C_length = max(V3)) %>%
  arrange(desc(C_ave))
  







