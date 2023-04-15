# DEG characterization.
library(dplyr)
myGFF <- "/System/Volumes/Data/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gff3"
newGFF <- rtracklayer::import(myGFF)

genes = as.data.frame(newGFF)  %>%
        dplyr::filter(type=="gene") 

genes
# treatment effect
trt_deg <- names(get_significant_results(m.Vem_c_ed))
  
de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)
hist(as.integer(de_gr$seqnames), breaks =  0:10)
hist(de_gr[de_gr$seqnames==4,"start"]/1e6)

#:) this look good

# inv4m effect



# genottype x treatment effect




