#prepare file with variants for filtering with 4 columns, one row per variant with interval to query (we did variant +/- 200bps), one file per proband, named freebayes_input
#Chromosome, start, end, proband ID

#prepare trio sample dictionary with 3 columns, one row per person, named <proband ID>.bed
#sample cram location, sample ID, proband ID

#launch freebayes for each sample which must first pipe cram file using variant regions for proband to make a new vcf for each sample
#samtools view sample.cram proband.bed -T referenceGenome.fasta -b | freebayes -stdin -f referenceGenome.fasta --vcf sampleID.fb.vdf

#repeatmask annotate freebayes VCFs

#filter by comparing with GATK output, 7 columns with one line per variant, make one file per proband
#proband ID, chromosome, position, ref, alt, tier (we did variant quality), unique variant ID

#run comparison R code on trio VCFs and gatk variant file
#srun trio_vcf_filter.R proband.vcf mom.vcf dad.vcf gatk.variant.file

