### Example of where to do lfcshrink
# FS12b.glom is a phyloseq object

FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'F')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')

resultsNames(FS12.de)

tmpres <- results(FS12.de, name = 'shed_low_vs_high', cooksCutoff = FALSE)
tmpres <- lfcShrink(FS12.de, res=tmpres, coef = 'shed_low_vs_high', type = 'apeglm')

tmpres[tmpres$padj < 0.1,]



#### example of how to plot NMDS centroids ####

library(funfuns)

library(tidyverse)

OTU <- read_tsv('data/stability.outdoubletons.abund.opti_mcc.0.03.subsample.shared')
META <- read_csv('data/FS9_metadata.csv')
OTU <- OTU %>% column_to_rownames(var = 'Group')
OTU <- OTU[,-c(1:2)]
rownames(OTU) %in% META$Sample
META <- META[META$Sample %in% rownames(OTU),]
OTU <- OTU[rownames(OTU) %in% META$Sample,]
META <- META[match(rownames(OTU), META$Sample),]
rownames(OTU) == META$Sample
META <- META %>% column_to_rownames(var = 'Sample')
META$Set <- paste(META$Treatment, META$Day, sep='_')
rownames(OTU) == rownames(META)

library(funfuns)


metaNMDS <- NMDS_ellipse(metadata = META, OTU_table = OTU, grouping_set = 'Set')


meta_with_NMDS <- metaNMDS[[1]]

# technically this will plot each centroid too many times, 
# one point on top of the other, for each observation in the 'Set'

meta_with_NMDS %>%
  ggplot(aes(x=MDS1, y=MDS2)) + geom_point(color='grey') + 
  geom_point(aes(x=centroidX, y=centroidY, color=Set))


CENTROIDS <- meta_with_NMDS %>%
  select(Set, Treatment, Day, centroidX, centroidY) %>% 
  unique()

#this is probably a better way:
meta_with_NMDS %>%
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(color='grey') + 
  geom_point(data=CENTROIDS, aes(x=centroidX, y=centroidY, color=Treatment, shape=Day), size=4)

