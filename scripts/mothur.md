# **Mothur**

## **Notes**
Used Silva release version 132 taxonomy and fasta files.
* silva.nr_v132.tax
* silva.nr_v132.pcr.align (renamed as silva.v4.fasta)

## Commands used for mothur (version 1.40.1)
```
make.file(inputdir=/project/fastq/, type=fastq, prefix=stability)

make.contigs(file=stability.files, processors=40)

summary.seqs(fasta=stability.trim.contigs.fasta)

screen.seqs(fasta=stability.trim.contigs.fasta, maxambig=0, maxlength=275, maxhomop=6, group=stability.contigs.groups)

summary.seqs(fasta=current)

unique.seqs(fasta=stability.trim.contigs.good.fasta)

count.seqs(name=current, group=stability.contigs.good.groups)

summary.seqs(count=current, fasta=current, processors=39)

pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F)

rename.file(input=silva.nr_v132.pcr.align, new=silva.v4.fasta)

summary.seqs(fasta=silva.v4.fasta)

align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta, flip=T)

summary.seqs(fasta=current, count=stability.trim.contigs.good.count_table)

screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=6)

summary.seqs(fasta=current, count=current)

filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

pre.cluster(fasta=current, count=current, diffs=2)

chimera.vsearch(fasta=current, count=current, dereplicate=t)

remove.seqs(fasta=current, accnos=current)

summary.seqs(fasta=current, count=current)

classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=silva.v4.fasta, taxonomy=silva.nr_v132.tax, cutoff=80)

remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

summary.tax(taxonomy=current, count=current)

get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=Mock5)

seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, reference=mock_community_16S3.fna, aligned=F)

#Rename these files to use for split.abund command (remove sequences with 2 or fewer occurrences)
#rename.file(input=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, new=stability.outdoubletons.fasta)
#rename.file(input=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, new=stability.outdoubletons.count_table)

#Remove Mock group
remove.groups(count=stability.outdoubletons.count_table, fasta=stability.outdoubletons.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, groups=Mock5)

#Split out doubletons from fasta and count table data (they'll be put in the "rare" files)
split.abund(fasta=stability.outdoubletons.fasta, count=stability.outdoubletons.count_table, cutoff=2, accnos=true)

dist.seqs(fasta=stability.outdoubletons.abund.fasta, cutoff=0.03)

cluster(column=current, count=stability.outdoubletons.abund.count_table)

make.shared(list=current, count=stability.outdoubletons.abund.count_table, label=0.03)

classify.otu(list=stability.outdoubletons.abund.opti_mcc.list, count=stability.outdoubletons.abund.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label=0.03)

count.groups(shared=stability.outdoubletons.abund.opti_mcc.shared)

sub.sample(shared=stability.outdoubletons.abund.opti_mcc.shared, size=1300)
```
