library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

rna_ctn <- read.table("~/Desktop/tuberna_ctn.tab", sep = "\t", quote ="", header = TRUE)
sample_info <- read.csv("~/Desktop/PENN-PHO22-TUBES.csv", quote ="", header =TRUE)
plots <- read.csv("~/Desktop/PENN_PHO22_PLOTS.csv")

colnames(metadata)

metadata <- rna_ctn %>%
  dplyr::inner_join(sample_info, by= c(tube= "top_tag")) %>%
  dplyr::inner_join(plots, by= c(row ="PHO22")) 

Treatment

# correct R33 and # R34 they have the leaves swapped
metadata$leaf_tissue[metadata$tube == 'R33'] = 2
metadata$leaf_tissue[metadata$tube == 'R34'] = 1

metadata <- within(metadata,{
TIME = sub("1:","13:",metadata$TIME)
TIME = hm(TIME)
decimal_time <- hour(TIME)  + minute(TIME)/60 + second(TIME) / 3600
}
)

metadata


write.csv(
  metadata %>%
    dplyr::select(
      tube:NCSU_RNA_plant,
      Treatment, genotype, leaf_tissue, side_tag,
      TIME, decimal_time, everything()
      ), 
  file="../data/inv4mRNAseq_metadata.csv", row.names = FALSE)

kalqc <- read.table("~/Desktop/kallisto_qc.tab", sep = "\t", quote ="", header = FALSE, skip =1)
colnames(kalqc) <- c("gsl_sample","pseudoaligned_pct","pseudoaligned", "processed")
kalqc$tube <- substr(kalqc$gsl_sample,1,3)

kalqc

#here you can decide to merge samples!!!

colnames(sample_info)

# check leaf number
ggplot(sampleInfo, aes(y=leaf_number, x =genotype)) + 
  ggbeeswarm::geom_quasirandom() + facet_wrap(~Treatment)
anova(lm(data=sampleInfo,leaf_number ~ Treatment*genotype))


#sum over genes
# make gene _sample table

# if you are thinikg in adding readcounts from poor qc libraries together
# sampleInfo[sampleInfo$genotype=="CTRL" & sampleInfo$Treatment=="Low_P" & sampleInfo$leaf_tissue ==1, c("tube","row","Rep")]

# this is for runing in the server

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
