library(tidyverse)

data = read.csv('/Users/cmdb/qbb2024-answers/week10/output.txt')

ggplot(data, aes(x = factor(Gene), y = nascentRNA)) +
  geom_violin() +
  labs(title = "Nascent mRNA expression for knockouts",
       x = "Gene knockout",
       y = "Nascent mRNA levels") +
  scale_color_manual(values = c("red", "blue", "green", "yellow"))
