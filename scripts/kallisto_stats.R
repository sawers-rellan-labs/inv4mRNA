library(dplyr)
library(ggplot2)

rna_ctn <- read.table("~/Desktop/tuberna_ctn.tab", sep = "\t", quote ="", header = TRUE)
sample_info <- read.csv("~/Desktop/PENN-PHO22_TUBES.csv", quote ="", header =TRUE)
plots <- read.csv("~/Desktop/PENN_PHO22_PLOTS.csv")

metadata <- rna_ctn %>%
  dplyr::inner_join(sample_info, by= "tube") %>%
  dplyr::inner_join(plots, by= c(row ="PHO22"))

metadata$leaf_tissue <- sub(".*_L","", metadata$side_tag, perl =TRUE) 

metadata$leaf_tissue <- sub("L","", metadata$leaf_tissue, perl =TRUE) 
metadata$leaf_tissue <- sub("R","", metadata$leaf_tissue, perl =TRUE) 
metadata$leaf_tissue <- as.integer(metadata$leaf_tissue)
metadata$leaf_tissue


metadata$side_tag

metadata$TIME<- sub("1:","13:",metadata$TIME)
# search gsl sample sheet

write.csv(metadata, file="inv4mRNAseq_metadata.csv")

write.csv()
kalqc <- read.table("~/Desktop/kallisto_qc.tab", sep = "\t", quote ="", header = FALSE, skip =1)
colnames(kalqc) <- c("gsl_sample","pseudoaligned_pct","pseudoaligned", "processed")
kalqc$tube <- substr(kalqc$gsl_sample,1,3)

kalqc

colnames(sample_info)


#sum over genes
# make gene _sample table

library(dplyr)
library(tidyr)
samples <- dir("./quant_out")

all_exp <-lapply(samples, function(x){
  print(x)
  sample_exp <- read.table(file.path("quant_out", x,"abundance.tsv"), sep = "\t", quote ="", header = TRUE)
  sample_exp$est_counts <-as.integer(sample_exp$est_counts)
  sample_exp$gene = factor(sub("_T.+$","",sample_exp$target_id, perl =TRUE))

  out <- sample_exp %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(counts = sum(est_counts))
  out$sample <- x
  out
}) %>% dplyr::bind_rows()



gene_sample_exp <- all_exp %>% pivot_wider( names_from = "sample", values_from ="counts" )


write.csv(gene_sample_exp,"gene_sample_exp.csv")











kalqc %>% filter(pseudoaligned_pct < 60) 
  
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
