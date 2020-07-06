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
