# DEG characterization.
library(dplyr)
source("~/Projects/NCSU/06_GEA/GEA/sorghum/fastman/fastman.R")

myGFF <- "/System/Volumes/Data/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gff3"
newGFF <- rtracklayer::import(myGFF)

genes = as.data.frame(newGFF)  %>%
        dplyr::filter(type=="gene") 
genes$start
colnames(genes)



#### P effect ####

# DEG for GO terms
trt_deg <- names(get_significant_results(results$`2`$m.Vem_ed,thresh = 1e-6))
length(trt_deg)
cat(trt_deg)
de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)
write.table(trt_deg, file= "P_trt_DEG.txt",
            row.names = FALSE, col.names = FALSE,
            quote= FALSE)
lsfr <- as.data.frame(get_lfsr(results$`2`$m.Vem_ed))
p <- pmin( lsfr$V1,lsfr$V2,lsfr$V3,lsfr$V4)
sum(p<0.05)
length(trt_deg)

str(p)
to_plot <- data.frame( gene_id = names(p), P = p) %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)



quartz(height = 4, width = 12)
fastman (to_plot,
         cex.axis = 1.3, cex.lab=1.3, main = "DEG Phosphorus Treatment",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))

## inv4m effect ####
trt_deg <- names(get_significant_results(results$`3`$m.Vem_ed))
de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)

#:) this look good
p <- get_lfsr(results$`3`$m.Vem_ed)[,3]
str(p)
to_plot <- data.frame( gene_id = names(p), P = p) %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)

quartz(height = 4, width = 12)

fastman (to_plot,
         cex.axis = 1.3, cex.lab=1.3, main = "DEG INV4m effect",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))


# genotype x treatment effect

trt_deg <- names(get_significant_results(results$`4`$m.Vem_ed))
de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)

#:) this look good
to_plot <- as.data.frame(get_lfsr(results$`4`$m.Vem_ed)) %>%
  mutate(P=pmin(V1,V2,V3,V4)) %>%
  tibble::rownames_to_column("gene_id") %>%
  inner_join(genes %>%
               dplyr::select(gene_id, CHR = "seqnames",BP = "start")
  ) %>%
  dplyr::select(gene_id, CHR, BP, P)
hist(to_plot$P)
quartz(height = 4, width = 12)
fastman (to_plot,
         cex.axis = 1.3, cex.lab=1.3, main = "DEG INV4m x Phosphorus interaction effect",
         suggestiveline = -log10(0.05),
         genomewideline = -log10(0.05))




