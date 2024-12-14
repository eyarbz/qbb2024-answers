library(tidyverse)

data = read.csv('/Users/cmdb/qbb2024-answers/week10/output.txt')

ggplot(data, aes(x = factor(Gene), y = nascentRNA)) +
  geom_violin(fill = "darkgreen", color = "black") +
  labs(title = "Nascent mRNA expression for knockouts",
       x = "Gene knockout",
       y = "Nascent mRNA levels")

ggplot(data, aes(x = factor(Gene), y = PCNA)) +
  geom_violin(fill = "darkblue", color = "black") +
  labs(title = "PCNA expression for knockouts",
       x = "Gene knockout",
       y = "PCNA levels")

ggplot(data, aes(x = factor(Gene), y = ratio)) +
  geom_violin(fill = "darkred", color = "black") +
  labs(title = "PCNA expression for knockouts",
       x = "Gene knockout",
       y = "log2(RNA/PCNA)")

