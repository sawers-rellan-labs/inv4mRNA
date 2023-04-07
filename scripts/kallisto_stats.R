library(dplyr)
library(ggplot2)

kalqc <- read.table("~/Desktop/kallisto_qc.tab", sep = "\t", quote ="", header = FALSE, skip =1)
colnames(kalqc) <- c("sample","pseudoaligned_pct","pseudoaligned", "processed")

quartz()
kalqc %>%
  tidyr::pivot_longer(pseudoaligned_pct:processed, names_to = "var", values_to = "value" ) %>%
  ggplot2::ggplot(aes(x=value))+
  ggtitle("kallisto pseudoalignment") +
  ylab("# tissue libraries") +
  xlab("read stat") +
  geom_histogram()+
  facet_wrap( ~ var, scales = "free")+
  ggpubr::theme_classic2()

kalqc %>% filter(pseudoaligned_pct < 60)


  