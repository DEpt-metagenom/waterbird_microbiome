library(ggplot2)
library(viridis)
library(cowplot)
getwd()
df<-read.csv("merged_recount_pozK.tsv", sep="\t")

p_species<-ggplot(df, aes(fill=species, y=abundance, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_hue(l=30, c=200) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_species
ggsave("species_prop_recount_pozK.png", width=20, units="in", limitsize=F)

df$error<-(0.12-df$abundance)-(2*(0.12-df$abundance))

p_err<-ggplot(df, aes(x = sample, y = error)) +
  geom_boxplot(aes()) +
  geom_jitter(aes(col = species, stroke=2), size = 3, shape=21) +
  labs(x = "Sample", y = "Difference (12% - abundance)", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
p_err
ggsave("species_prop_err_pozK.png")
