library(edgeR)
library(data.table)
library(ggplot2)
library(limma)

counts = fread('~/Desktop/inv4mRNAseq_gene_sample_exp.csv',data.table=F)
sampleInfo = fread('~/Desktop/inv4mRNAseq_metadata.csv',data.table=F)

ggplot(sampleInfo, aes(y=leaf_number, x =genotype)) + 
  ggbeeswarm::geom_quasirandom() + facet_wrap(~Treatment)
anova(lm(data=sampleInfo,leaf_number ~ Treatment*genotype))



genes = data.frame(gene = counts[,2])
counts = as.matrix(counts[,-c(1:2)])
rownames(counts) = genes$gene

sampleNames = colnames(counts)
sampleNames = sapply(sampleNames,function(x) strsplit(x,'_')[[1]][1])
sampleNames %in% sampleInfo$tube
sampleInfo = sampleInfo[match(sampleNames,sampleInfo$tube),]

sampleInfo$leaf_tissue[sampleInfo$tube == 'R33'] = 2
sampleInfo$leaf_tissue[sampleInfo$tube == 'R34'] = 1

y = DGEList(counts = counts,samples = sampleInfo)
y$group = interaction(y$samples$Treatment,y$samples$genotype)
keep = filterByExpr(y,group = y$group)

y$samples$lowCount = y$samples$lib.size < 2e7

y_filtered = y[keep,]
mds = plotMDS(y_filtered,pch=21,label = y$samples$side_tag,bg = (as.factor(y$samples$lowCount)==T)+1)
plot(mds$x,mds$y,pch=21,bg = (as.factor(y_filtered$samples$lowCount)==T)+1)

y_filtered_bySample = y_filtered[,!y_filtered$samples$lowCount]

y_filtered_bySample$samples
mds2 = plotMDS(y_filtered_bySample, label = y$samples[!y$samples$lowCount,"side_tag"])

d = y_filtered_bySample$samples
d$x = mds2$x
d$y = mds2$y
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Treatment))
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = as.numeric(sub(':','.',d$TIME,fixed=T))))
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = COLLECTOR))
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = factor(leaf_tissue),shape = Treatment),size=3) #+ geom_text(aes(label = tube))

plot(mds2$x,mds2$y,pch=21,bg = factor(y_filtered_bySample$samples$Treatment),col=0)
plot(mds2$x,mds2$y,col = factor(y_filtered_bySample$samples$genotype))


design = model.matrix(~factor(leaf_tissue) + genotype + Treatment:genotype,d)
voomR = voom(y_filtered_bySample,design=design)

fit = lmFit(voomR)
ebfit = eBayes(fit)

head(ebfit$coefficients)

topTable(ebfit,coef = 6)
topTable(ebfit,coef = 7)

design_interaction = model.matrix(~factor(leaf_tissue) + Treatment*genotype,d)
fit = lmFit(voomR,design_interaction)
ebfit = eBayes(fit)

head(ebfit$coefficients)
topTable(ebfit,coef = 6)

d$y = voomR$E['Zm00001eb233320',]
# ggplot(d,aes(x=Treatment,y=y)) + geom_point(aes(color = genotype,group = interaction(Treatment,genotype)),position = position_jitterdodge()) + facet_wrap(~leaf_tissue)
ggplot(d,aes(x=Treatment,y=2^y)) + geom_boxplot(aes(color = genotype,group = interaction(Treatment,genotype))) + facet_wrap(~leaf_tissue)
ggplot(d,aes(x=genotype,y=2^y)) + geom_boxplot(aes(color = Treatment,group = interaction(Treatment,genotype))) + facet_wrap(~leaf_tissue)

# by leaf analysis #####################################
# 1 split y_filtered_bSample by leaf_tissue
y_filtered_bySample$samples$group <- y_filtered_bySample$samples$leaf_tissue

by_leaf <- splitIntoGroups(y_filtered_bySample)

class(by_leaf$`1`)

# make design in each = ~genotype*

get_top_per_leaf <- function(coef=4){
  # voom -> lmFit -> eBayes -> topTable
  top_genes <- lapply(1:4, function(x){
     y_filtered_by_leaf = y_filtered[,y_filtered$samples$leaf_tissue==x]
     d = y_filtered_by_leaf$samples
     design_interaction = model.matrix(~Treatment*genotype,d)
     voomR = voom(y_filtered_by_leaf ,design=design_interaction)
     fit = lmFit(voomR)
     ebfit = eBayes(fit)
     print(colnames(ebfit$coefficients))
     tt <- topTable(ebfit,coef = coef)
     tt$leaf_tissue = x
     tt
  })
  print(coef_names)
  top_genes 
}

get_top_per_leaf(coef=4)


y_filtered$samples$leaf_bin <- ifelse(y_filtered$samples$leaf_tissue >2,2,1)


get_top_per_leaf <- function(coef=4){
  # voom -> lmFit -> eBayes -> topTable
  top_genes <- lapply(1:2, function(x){
    y_filtered_by_leaf = y_filtered[,y_filtered$samples$leaf_bin==x]
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~Treatment*genotype,d)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    print(colnames(ebfit$coefficients))
    tt <- topTable(ebfit,coef = coef)
    tt$leaf_bin = x
    tt
  })
  print(coef_names)
  top_genes 
}

get_top_per_leaf(coef=4)
