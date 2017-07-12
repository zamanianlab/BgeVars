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
  arrange(desc(C_ave))


write.table(summ.data, "~/data/Bge/summ_data.csv", sep = "\t")






