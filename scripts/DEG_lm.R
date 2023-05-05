library(edgeR)
library(data.table)
library(ggplot2)
library(limma)
library(ashr)
library(mashr)

counts = fread('../data/inv4mRNAseq_gene_sample_exp.csv',data.table=F)
sampleInfo = fread('../data/inv4mRNAseq_metadata.csv',data.table=F)
# counts = fread('inv4mRNAseq_gene_sample_exp.csv',data.table=F)


genes = data.frame(gene = counts[,2])
counts = as.matrix(counts[,-c(1:2)])
rownames(counts) = genes$gene

sampleNames = colnames(counts)
sampleNames = sapply(sampleNames,function(x) strsplit(x,'_')[[1]][1])
sampleNames %in% sampleInfo$tube
sampleInfo = sampleInfo[match(sampleNames,sampleInfo$tube),]

y = DGEList(counts = counts,samples = sampleInfo)
y$group = interaction(y$samples$Treatment,y$samples$genotype)
keep = filterByExpr(y,group = y$group)

y$samples$lowCount = y$samples$lib.size < 2e7

y_filtered = y[keep,]

mds = plotMDS(y_filtered,pch=21,label = y$samples$side_tag,bg = (as.factor(y$samples$lowCount)==T)+1)
quartz()
plot(mds$x,mds$y,pch=21,bg = (as.factor(y_filtered$samples$lowCount)==T)+1, 
     main = "Expression MDS ")

y_filtered_bySample = y_filtered[,!y_filtered$samples$lowCount]

y_filtered_bySample$samples
table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$leaf_tissue)
table(y_filtered_bySample$samples$genotype,
      y_filtered_bySample$samples$leaf_tissue)

table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$genotype,
      y_filtered_bySample$samples$leaf_tissue)




d = y_filtered_bySample$samples
d$x = mds2$x
d$y = mds2$y

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Treatment))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = decimal_time ))



quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = COLLECTOR))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(
  aes(color = as.factor(row),
      shape = Treatment)) 

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = genotype))

quartz()
d$Treatment <- factor(d$Treatment,levels = c("Low_P","High_P"))

ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = factor(leaf_tissue),shape = Treatment),size=3) 

plot(mds2$x,mds2$y,pch=21,bg = factor(y_filtered_bySample$samples$Treatment),col=0)
plot(mds2$x,mds2$y,col = factor(y_filtered_bySample$samples$genotype))


d$Treatment <- factor(d$Treatment,levels = c("High_P","Low_P"))
design = model.matrix(~leaf_tissue + genotype + Treatment:genotype,d)
y_filtered_bySample = calcNormFactors(y_filtered_bySample)
voomR = voom(y_filtered_bySample,design=design,plot=T)

fit = lmFit(voomR)
ebfit = eBayes(fit)

head(ebfit$coefficients)


design_interaction = model.matrix(~leaf_tissue*Treatment*genotype,d)
fit = lmFit(voomR,design_interaction)
ebfit = eBayes(fit)
head(ebfit$coefficients)
topTable(ebfit,coef = 4)
topTable(ebfit,coef = 5) 
topTable(ebfit,coef = 7) 
topTable(ebfit,coef = 8) 
cat(row.names(topTable(ebfit,coef = 2,number = 5000)))

#gene='Zm00001eb053360'
gene='Zm00001eb003820' # PILNCR1-miR399
gene='Zm00001eb191650' # PHOS2
d$y = voomR$E[gene,]
d$counts = y_filtered_bySample$counts[gene,]
# ggplot(d,aes(x=Treatment,y=y)) + geom_point(aes(color = genotype,group = interaction(Treatment,genotype)),position = position_jitterdodge()) + facet_wrap(~leaf_tissue)

quartz()
ggplot(d,aes(x=Treatment,y=y)) +
  ggtitle(gene) +
  geom_boxplot(aes(color = genotype,group = interaction(Treatment,genotype))) + 
  facet_wrap(~leaf_tissue) 

# ggplot(d,aes(x=Treatment,y=counts)) + geom_boxplot(aes(color = genotype,group = interaction(Treatment,genotype))) + facet_wrap(~leaf_tissue)
ggplot(d,aes(x=genotype,y=y))  +
  ggtitle(gene) +
  geom_boxplot(aes(color = Treatment,group = interaction(Treatment,genotype))) + 
  facet_wrap(~leaf_tissue)




########################################################
# by leaf analysis #####################################
########################################################

# make design in each = ~Treatment*genotype ############
# check effects for interaction





# Akaike information content calculation
aic <- function(x){
  # x is a mashr model object
  # I am assuming x$fitted_g$pi, the number of misture components,
  # is the number of parameters 
  2*length(x$fitted_g$pi) - 2*max(get_loglik(x))
}


## ------------------------------------------------------------------
# Step 4: Run mash for different models
# using just the simple canonical covariances as in the initial introductory vignette.
# correcting for correlations between samples inside a "condition"
# 



tb <- data.frame()
results <- list()

for(coef in 2:4){

  
  effects = SEs = matrix(NA,nrow = nrow(y_filtered_bySample),ncol = 4)
  rownames(effects) = rownames(SEs) = rownames(y_filtered_bySample)
  
  for(x in 1:4) {
    y_filtered_by_leaf = y_filtered_bySample[,y_filtered_bySample$samples$leaf_tissue==x]
    y_filtered_by_leaf$group = interaction(y_filtered_by_leaf$samples$genotype,y_filtered_by_leaf$samples$Treatment)
    # filter_genes = filterByExpr(y_filtered_by_leaf,group = y_filtered_by_leaf$group)
    # y_filtered_by_leaf = y_filtered_by_leaf[filter_genes,]
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~Treatment*genotype,d)
    y_filtered_by_leaf = calcNormFactors(y_filtered_by_leaf)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction,plot=T)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    tt <- topTable(ebfit,coef = coef,sort.by = 'none',n=Inf)
    effects[,x] = tt$logFC
    SEs[,x] = tt$logFC/tt$t
  }
  

  
  # Step 2: Obtain initial data-driven covariance matrices
  
  data = mash_set_data(as.matrix(effects), as.matrix(SEs))
  U.c = cov_canonical(data)
  
  # Step 2: Obtain initial data-driven covariance matrices
  
  # select strong signals
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1,0.05)
  if(length(strong)<20) {
    strong = order(apply(m.1by1$result$lfsr,1,min))[1:min(20,nrow(effects))]
  }
  
  # Perform PCA on data and return list of candidate covariance matrices
  U.pca = cov_pca(data,npc=ncol(effects),subset=strong)
  # npc:	the number of PCs to use, should be less than or equal to n_conditions(data)
  # subset: indices of the subset of data to use (set to NULL for all data)
  # print(names(U.pca))
  
  ## ----------------------------------------------------------------
  # Step 3: prepare canonical/data-driven covariance matrices
  # Perform "extreme deconvolution" (Bovy et al) on a subset of the data
  U.ed = cov_ed(data, U.pca, subset=strong)
  # subset: a subset of data to be used when ED is run (set to NULL for all the data)
  # The function cov_ed is used to apply the ED algorithm from a specified initialization
  # (here U.pca) and to a specified subset of signals.
  
  
  
  # Why U.c and U.ed?  They are supposed to be different ways of estimating the conditions covariance matrix
  V.em_c_ed = mash_estimate_corr_em(data, c(U.c,U.ed), details = TRUE)
  m.Vem_c_ed = V.em_c_ed$mash.model
  m.Vem_c_ed$result$NAs = is.na(effects)
  m.Vem_c_ed$V = V.em_c_ed$V
  
  print(get_loglik(m.Vem_c_ed),digits=10) 
  length(m.Vem_c_ed$fitted_g$pi)
  aic(m.Vem_c_ed)
  
  results[[as.character(coef)]][["m.Vem_c_ed"]] <- m.Vem_c_ed
  
  
  ###
  V.em_ed = mash_estimate_corr_em(data, U.ed, details = TRUE)
  m.Vem_ed = V.em_ed$mash.model
  m.Vem_ed$result$NAs = is.na(effects)
  
  quartz()
  mash_plot_meta(m.Vem_ed,2)
  
  m.Vem_ed$V = V.em_ed$V
  
  print(get_loglik(m.Vem_ed),digits=10) 
  length(m.Vem_ed$fitted_g$pi)
  aic(m.Vem_ed)
  
  results[[as.character(coef)]][["m.Vem_ed"]] <- m.Vem_ed
  
  ###
  V.em_c = mash_estimate_corr_em(data, U.c, details = TRUE)
  m.Vem_c = V.em_c$mash.model
  m.Vem_c$result$NAs = is.na(effects)
  m.Vem_c$V = V.em_c$V
  
  print(get_loglik(m.Vem_c),digits=10) 
  length(m.Vem_c$fitted_g$pi)
  aic(m.Vem_c)
  
  results[[as.character(coef)]][["m.Vem_c"]] <- m.Vem_c
  ###
  m.c = mash(data,U.c)
  # Fitting model with 217 mixture components
  # -57582
  
  print(get_loglik(m.c),digits=10)
  length(m.c$fitted_g$pi)
  aic(m.c)
  
  
  ###
  m.c_ed  = mash(data, c(U.c,U.ed))
  
  # Fitting model with 337 mixture components.
  
  print(get_loglik(m.c_ed),digits=10)
  length(m.c_ed$fitted_g$pi)
  aic(m.c_ed)
  
  results[[as.character(coef)]][["m.c_ed"]] <- m.c_ed
  ####
  # Compare likelihood of the  models
  
  ## Compare model fit
  
  
  model =  c('canonical',
             'canonical + ed',
             'canonical + Vem',
             'ed + Vem', 
             'canonical + ed + Vem')
  
  significant = c(length(get_significant_results(m.c)), 
                  length(get_significant_results(m.c_ed)),
                  length(get_significant_results(m.Vem_c)),
                  length(get_significant_results(m.Vem_ed)),
                  length(get_significant_results(m.Vem_c_ed)))
  
  loglike = c(
    max(get_loglik(m.c)),
    max(get_loglik(m.c_ed)),
    max(get_loglik(m.Vem_c)),
    max(get_loglik(m.Vem_ed)),
    max(get_loglik(m.Vem_c_ed))
  )
  
  k = c(length(m.c$fitted_g$pi), 
        length(m.c_ed$fitted_g$pi), 
        length(m.Vem_c$fitted_g$pi),
        length(m.Vem_ed$fitted_g$pi), 
        length(m.Vem_c_ed$fitted_g$pi))
  
  AIC = c(aic(m.c), aic(m.c_ed), aic(m.Vem_c), aic(m.Vem_ed), aic(m.Vem_c_ed))
  
  tb <- rbind ( tb,
    data.frame(
      coef = coef,
      model = model,
      significant = significant,
      loglike = loglike,
      k = k,
      AIC = AIC)
    )
}
coef_names <- colnames(head(ebfit$coefficients))
tb$coef_name <- coef_names[tb$coef]
library(dplyr)
tb %>%
  dplyr::select(coef, coef_name, everything())

tb %>%
  dplyr::select(coef, coef_name, everything()) %>%
  dplyr::group_by(coef) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1)



length(get_significant_results(results$`4`$m.Vem_ed))
results$`4`$m.Vem_ed
get_lfsr(results$`4`$m.Vem_ed)


deg_idx <- get_significant_results(results$`4`$m.Vem_ed)
deg <- names(deg_idx)
deg_1 <- deg
length(deg)
paste(deg, collapse =  ",")
deg_lfsr <- get_lfsr(results$`4`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,4]),]

head(get_pm(results$`4`$m.Vem_ed)[deg_idx,])

barplot(get_estimated_pi(results$`4`$m.Vem_ed),las = 2)

print(get_pairwise_sharing(results$`4`$m.Vem_ed)) 



deg_idx <- get_significant_results(results$`2`$m.Vem_ed)
deg <- names(deg_idx)
dd <- deg[deg %in% deg_1]

paste(dd, collapse =  ",")
deg_lfsr <- get_lfsr(results$`2`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]


barplot(get_estimated_pi(results$`3`$m.Vem_ed),las = 2)
print(get_pairwise_sharing(results$`2`$m.Vem_ed))



