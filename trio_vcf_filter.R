#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(GenomicRanges)
library(gtools)
options(scipen=999) #default is 0

args <- commandArgs(trailingOnly = TRUE)
#1 is proband, 2 is mom, 3rd argument is dad, where the arguments are the VCF files to merge
#arg 4 is the full path to the gatk file.

proband_file <- args[1]
proband <- read.delim(proband_file, stringsAsFactors=F, header=F, comment.char="#", colClasses=c("V1"="character"))
cmd <- paste("tabix -H", proband_file, "| tail -1")
proband_columns <- system(cmd, intern=T)
proband_columns <- unlist(str_split(proband_columns, "\\t"))
cmd <- paste("tabix -H", proband_file)
proband_header <- system(cmd, intern=T)
proband_ID <- proband_columns[length(proband_columns)]
proband_columns[1] <- "X.CHROM"
proband_columns[length(proband_columns)] <- "proband"
names(proband) <- proband_columns

mom_file <- args[2]
mom <- read.delim(mom_file, stringsAsFactors=F, header=F, comment.char="#", colClasses=c("V1"="character"))
cmd <- paste("tabix -H", mom_file, "| tail -1")
mom_columns <- system(cmd, intern=T)
mom_columns <- unlist(str_split(mom_columns, "\\t"))
mom_ID <- mom_columns[length(mom_columns)] 
mom_columns[1] <- "X.CHROM"
mom_columns[length(mom_columns)] <- "mom"
names(mom) <- mom_columns

dad_file <- args[3]
dad <- read.delim(dad_file, stringsAsFactors=F, header=F, comment.char="#", colClasses=c("V1"="character"))
cmd <- paste("tabix -H", dad_file, "| tail -1")
dad_columns <- system(cmd, intern=T)
dad_columns <- unlist(str_split(dad_columns, "\\t"))
dad_ID <- dad_columns[length(dad_columns)] 
dad_columns[1] <- "X.CHROM"
dad_columns[length(dad_columns)] <- "dad"
names(dad) <- dad_columns

#read in GATK table file
gatk_file <- args[4]
gatk <- read.delim(gatk_file, stringsAsFactors=F, colClasses=c("X.CHROM"="character"))

#update the header for the final VCF output
proband_header[length(proband_header)] <- paste(proband_header[length(proband_header)], "\t", mom_ID, "\t", dad_ID, sep="")

joined_vcf <- left_join(proband, mom,
                        by=c("X.CHROM", "POS", "REF", "ALT"))
joined_vcf <- left_join(joined_vcf, dad,
                        by=c("X.CHROM", "POS", "REF", "ALT"))
vcf <- joined_vcf %>% select(X.CHROM, POS, ID.x, REF, ALT, QUAL.x, FILTER.x, INFO.x, FORMAT.x, proband, mom, dad)
names(vcf) <- c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "proband", "mom", "dad")
proband_gt <- str_split(vcf$proband, ":") %>% sapply(., "[[", 1)
proband_gt <- str_split(proband_gt, "[/|]")
proband_gt <- lapply(proband_gt, function(x) {x[x=="."] <- "0"; return(x)})
proband_gt <- lapply(proband_gt, as.integer)
proband_gt <- sapply(proband_gt, sum)
vcf$proband2 <- proband_gt
vcf$mom[is.na(vcf$mom)] <- "."
vcf$dad[is.na(vcf$dad)] <- "."

vcf_denovo <- vcf %>% filter(proband2 > 0  & mom == "." & dad == ".")

#######
#start filtering VCF file
#filter on simple repeats and low complexity
message("Filtering.")
simple_repeat <- grepl("Simple_repeat", vcf_denovo$INFO)
low_complexity <- grepl("Low_complexity", vcf_denovo$INFO)
remove <- simple_repeat | low_complexity

vcf_simple_repeat <- vcf_denovo[simple_repeat, ]
vcf_low_complexity <- vcf_denovo[low_complexity, ]

vcf_denovo_filtered <- vcf_denovo[!remove, ]

#futher filter on QUAL < 1
vcf_low_qual <- vcf_denovo_filtered %>% filter(QUAL < 1)
vcf_denovo_filtered <- vcf_denovo_filtered %>% filter(QUAL >= 1)

#prepare for assessing overlapping intervals
AnnotateBEGIN_END <- function(POSITION, REF, ALT) {
  nchar_REF <- nchar(REF)
  nchar_ALT <- nchar(ALT)
  if (nchar_REF == 1 & nchar_ALT == 1) {
    POS1 <- POSITION
    POS2 <- POSITION
  } else if (nchar_REF > nchar_ALT) {
    POS1 <- POSITION
    POS2 <- POSITION + nchar_REF -1
  } else {
    POS1 <- POSITION
    POS2 <- POSITION
  }
  return(data.frame(BEGIN=POS1, END=POS2))
}
AnnotateVarType <- function(REF, ALT) {
  nchar_REF <- nchar(REF)
  nchar_ALT <- nchar(ALT)
  if (nchar_REF == 1 & nchar_ALT == 1) {
    return(data.frame(VarType = "SNP", stringsAsFactors=F))
  } else {
    return(data.frame(VarType = "INDEL", stringsAsFactors=F))
  }
}

positions <- mom %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
mom <- cbind(mom, positions)
vartype <- mom %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
mom <- cbind(mom, vartype)
mom$ID2 <- 1:nrow(mom)

positions <- dad %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
dad <- cbind(dad, positions)
vartype <- dad %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
dad <- cbind(dad, vartype)
dad$ID2 <- 1:nrow(dad)

positions <- vcf_denovo_filtered %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
if (invalid(positions)) {positions <- data.frame(BEGIN=integer(), END=integer())}
vcf_denovo_filtered_intervals <- cbind(vcf_denovo_filtered, positions)

vartype <- vcf_denovo_filtered %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
if (invalid(vartype)) {vartype <- data.frame(VarType=character(), stringsAsFactors=F)}
vcf_denovo_filtered_intervals <- cbind(vcf_denovo_filtered_intervals, vartype)

mom_ranges <- with(mom,
                   GRanges(X.CHROM,
                           IRanges(BEGIN, END),
                           ID2=ID2))
mom_indel_ranges <- with(mom %>% filter(VarType == "INDEL"),
                         GRanges(X.CHROM,
                                 IRanges(BEGIN, END),
                                 ID2=ID2))

dad_ranges <- with(dad,
                   GRanges(X.CHROM,
                           IRanges(BEGIN, END),
                           ID2=ID2))
dad_indel_ranges <- with(dad %>% filter(VarType == "INDEL"),
                         GRanges(X.CHROM,
                                 IRanges(BEGIN, END),
                                 ID2=ID2))

#for proband indels (NOT SNPS), if there is an overlapping parental indel 1 bases either side, then consider it inherited.
RemoveOverlappingIndels <- function(chr, BEGIN, END, vartype="x") {
  if (vartype == "SNP") {return(data.frame(ret=FALSE))}
  BEGIN <- as.integer(BEGIN)
  END <- as.integer(END)
  POS1 <- BEGIN - 2
  POS2 <- END + 2
  test_ranges <- GRanges(chr, IRanges(POS1, POS2))
  
  #perform intersection
  .intersection.mom <- suppressWarnings(intersect(test_ranges, mom_indel_ranges))
  .intersection.dad <- suppressWarnings(intersect(test_ranges, dad_indel_ranges))
  if (length(.intersection.dad) + length(.intersection.mom) > 0) {
    return(data.frame(ret=TRUE))
  } else {
    return(data.frame(ret=FALSE))
  }
} 
remove_overlapping_indels <- vcf_denovo_filtered_intervals %>% rowwise %>% do(RemoveOverlappingIndels(.$X.CHROM, .$BEGIN, .$END, .$VarType)) %>% unlist
if (invalid(remove_overlapping_indels)) {remove_overlapping_indels <- logical()}
vcf_overlapping_indel <- vcf_denovo_filtered_intervals[remove_overlapping_indels, ]

#remove nearby poor quality indel in parents (within 5 bases)
RemoveNearbyIndel <- function(chr, BEGIN, END) {
  BEGIN <- as.integer(BEGIN)
  END <- as.integer(END)
  POS1 <- BEGIN - 5
  POS2 <- END + 5
  test_ranges <- GRanges(chr, IRanges(POS1, POS2))
  
  #intersect with mom:
  .intersection <- suppressWarnings(intersect(test_ranges, mom_indel_ranges))
  if (length(.intersection) == 0) {
    mom_vars <- 0
  } else {
    .overlaps <- findOverlaps(test_ranges, mom_indel_ranges) #query is test_ranges, subject is mom_ranges
    mom_IDs_overlap <- unlist(mom_indel_ranges[subjectHits(.overlaps)]$ID2)
    mom_vars <- mom %>% filter(ID2 %in% mom_IDs_overlap) %>% filter(QUAL < 0.1) %>% nrow
  }
  #intersect with dad:
  .intersection <- suppressWarnings(intersect(test_ranges, dad_indel_ranges))
  if (length(.intersection) == 0) {
    dad_vars <- 0
  } else {
    .overlaps <- findOverlaps(test_ranges, dad_indel_ranges) #query is test_ranges, subject is dad_ranges
    dad_IDs_overlap <- unlist(dad_indel_ranges[subjectHits(.overlaps)]$ID2)
    dad_vars <- dad %>% filter(ID2 %in% dad_IDs_overlap) %>% filter(QUAL < 0.1) %>% nrow
  }
  total_bad_vars <- mom_vars + dad_vars
  if (total_bad_vars > 0) {
    return(data.frame(ret=TRUE))
  } else {
    return(data.frame(ret=FALSE))
  }
}
remove_nearby_indel <- vcf_denovo_filtered_intervals %>% rowwise %>% do(RemoveNearbyIndel(.$X.CHROM, .$BEGIN, .$END)) %>% unlist
if (invalid(remove_nearby_indel)) {remove_nearby_indel <- logical()}
vcf_nearby_indel <- vcf_denovo_filtered_intervals[remove_nearby_indel, ]

#remove any obvious homopolymer runs in REF or ALT of 5 or more
ref_homopolymer <- grepl("([A]){5}", vcf_denovo_filtered_intervals$REF, fixed=FALSE, perl=FALSE) |
  grepl("([C]){5}", vcf_denovo_filtered_intervals$REF, fixed=FALSE, perl=FALSE) | 
  grepl("([T]){5}", vcf_denovo_filtered_intervals$REF, fixed=FALSE, perl=FALSE) |
  grepl("([G]){5}", vcf_denovo_filtered_intervals$REF, fixed=FALSE, perl=FALSE)

alt_homopolymer <- grepl("([A]){5}", vcf_denovo_filtered_intervals$ALT, fixed=FALSE, perl=FALSE) |
  grepl("([C]){5}", vcf_denovo_filtered_intervals$ALT, fixed=FALSE, perl=FALSE) | 
  grepl("([T]){5}", vcf_denovo_filtered_intervals$ALT, fixed=FALSE, perl=FALSE) |
  grepl("([G]){5}", vcf_denovo_filtered_intervals$ALT, fixed=FALSE, perl=FALSE)

homopolymer <- ref_homopolymer | alt_homopolymer
vcf_homopolymer <- vcf_denovo_filtered_intervals[homopolymer, ]

#remove any dinucleotide or trinucleotide repeats of 3 or more in REF or ALT (dtnuc_repeat)
ref_dtnuc_repeat <- grepl("((AA){3}|(AC){3}|(AT){3}|(AG){3}|(CA){3}|(CC){3}|(CT){3}|(CG){3}|(TA){3}|(TC){3}|(TT){3}|(TG){3}|(GA){3}|(GC){3}|(GT){3}|(GG){3})", vcf_denovo_filtered_intervals$REF) | 
  grepl("(AAA){3}|(AAC){3}|(AAT){3}|(AAG){3}|(ACA){3}|(ACC){3}|(ACT){3}|(ACG){3}|(ATA){3}|(ATC){3}|(ATT){3}|(ATG){3}|(AGA){3}|(AGC){3}|(AGT){3}|(AGG){3}|(CAA){3}|(CAC){3}|(CAT){3}|(CAG){3}|(CCA){3}|(CCC){3}|(CCT){3}|(CCG){3}|(CTA){3}|(CTC){3}|(CTT){3}|(CTG){3}|(CGA){3}|(CGC){3}|(CGT){3}|(CGG){3}|(TAA){3}|(TAC){3}|(TAT){3}|(TAG){3}|(TCA){3}|(TCC){3}|(TCT){3}|(TCG){3}|(TTA){3}|(TTC){3}|(TTT){3}|(TTG){3}|(TGA){3}|(TGC){3}|(TGT){3}|(TGG){3}|(GAA){3}|(GAC){3}|(GAT){3}|(GAG){3}|(GCA){3}|(GCC){3}|(GCT){3}|(GCG){3}|(GTA){3}|(GTC){3}|(GTT){3}|(GTG){3}|(GGA){3}|(GGC){3}|(GGT){3}|(GGG){3}", vcf_denovo_filtered_intervals$REF)

alt_dtnuc_repeat <- grepl("((AA){3}|(AC){3}|(AT){3}|(AG){3}|(CA){3}|(CC){3}|(CT){3}|(CG){3}|(TA){3}|(TC){3}|(TT){3}|(TG){3}|(GA){3}|(GC){3}|(GT){3}|(GG){3})", vcf_denovo_filtered_intervals$ALT) | 
  grepl("(AAA){3}|(AAC){3}|(AAT){3}|(AAG){3}|(ACA){3}|(ACC){3}|(ACT){3}|(ACG){3}|(ATA){3}|(ATC){3}|(ATT){3}|(ATG){3}|(AGA){3}|(AGC){3}|(AGT){3}|(AGG){3}|(CAA){3}|(CAC){3}|(CAT){3}|(CAG){3}|(CCA){3}|(CCC){3}|(CCT){3}|(CCG){3}|(CTA){3}|(CTC){3}|(CTT){3}|(CTG){3}|(CGA){3}|(CGC){3}|(CGT){3}|(CGG){3}|(TAA){3}|(TAC){3}|(TAT){3}|(TAG){3}|(TCA){3}|(TCC){3}|(TCT){3}|(TCG){3}|(TTA){3}|(TTC){3}|(TTT){3}|(TTG){3}|(TGA){3}|(TGC){3}|(TGT){3}|(TGG){3}|(GAA){3}|(GAC){3}|(GAT){3}|(GAG){3}|(GCA){3}|(GCC){3}|(GCT){3}|(GCG){3}|(GTA){3}|(GTC){3}|(GTT){3}|(GTG){3}|(GGA){3}|(GGC){3}|(GGT){3}|(GGG){3}", vcf_denovo_filtered_intervals$ALT)

dtnuc_repeat <- ref_dtnuc_repeat | alt_dtnuc_repeat
vcf_dtnuc_repeat <- vcf_denovo_filtered_intervals[dtnuc_repeat, ]

#Now remove the additional filters:
remove_filters <- remove_overlapping_indels | remove_nearby_indel | homopolymer | dtnuc_repeat
vcf_denovo_filtered <- vcf_denovo_filtered[!remove_filters, ]

#keep track of those removed
columns_to_keep <- c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "proband", "mom", "dad")
vcf_simple_repeat <- vcf_simple_repeat[, columns_to_keep]
vcf_low_complexity <- vcf_low_complexity[, columns_to_keep]
vcf_low_qual <- vcf_low_qual[, columns_to_keep]
vcf_overlapping_indel <- vcf_overlapping_indel[, columns_to_keep]
vcf_nearby_indel <- vcf_nearby_indel[, columns_to_keep]
vcf_homopolymer <- vcf_homopolymer[, columns_to_keep]
vcf_dtnuc_repeat <- vcf_dtnuc_repeat[, columns_to_keep]

vcf_simple_repeat$FreeBayesRemoveReason <- rep("Simple_Repeat", nrow(vcf_simple_repeat))
vcf_low_complexity$FreeBayesRemoveReason <- rep("Low_Complexity", nrow(vcf_low_complexity))
vcf_low_qual$FreeBayesRemoveReason <- rep("Low_QUAL", nrow(vcf_low_qual))
vcf_overlapping_indel$FreeBayesRemoveReason <- rep("Overlapping_INDEL", nrow(vcf_overlapping_indel)) 
vcf_nearby_indel$FreeBayesRemoveReason <- rep("Nearby_Indel", nrow(vcf_nearby_indel))
vcf_homopolymer$FreeBayesRemoveReason <- rep("Homopolymer_Indel", nrow(vcf_homopolymer))
vcf_dtnuc_repeat$FreeBayesRemoveReason <- rep("dtNUC_REPEAT", nrow(vcf_dtnuc_repeat))

#prepare files for output
vcf_denovo <- vcf_denovo %>% select(-proband2)
vcf_denovo_filtered <- vcf_denovo_filtered %>% select(-proband2)
vars_removed <- rbind(vcf_simple_repeat, vcf_low_complexity,
                      vcf_low_qual, vcf_overlapping_indel,
                      vcf_nearby_indel,
                      vcf_homopolymer, vcf_dtnuc_repeat)
vars_removed <- vars_removed %>% group_by(X.CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, proband, mom, dad) %>% summarise(FreeBayesRemoveReason=paste(FreeBayesRemoveReason, collapse=";")) %>% as.data.frame

#compare with the gatk file, compare SNPs and INDELs separately
message("Comparing calls to GATK file.")
positions <- vcf_denovo_filtered %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
if (invalid(positions)) {positions <- data.frame(BEGIN=integer(), END=integer())}
vcf_denovo_filtered_ann <- cbind(vcf_denovo_filtered, positions)
vartype <- vcf_denovo_filtered_ann %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
if (invalid(vartype)) {vartype <- data.frame(VarType=character(), stringsAsFactors=F)}
vcf_denovo_filtered_ann <- cbind(vcf_denovo_filtered_ann, vartype)
#vcf_denovo_filtered_ann$ID2 <- 1:nrow(vcf_denovo_filtered_ann)
vcf_denovo_filtered_ann_SNPS <- vcf_denovo_filtered_ann %>% filter(VarType == "SNP")
vcf_denovo_filtered_ann_indel <- vcf_denovo_filtered_ann %>% filter(VarType == "INDEL")

positions <- vars_removed %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
if (invalid(positions)) {positions <- data.frame(BEGIN=integer(), END=integer())}
vars_removed_ann <- cbind(vars_removed, positions)
vartype <- vars_removed_ann %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
if (invalid(vartype)) {vartype <- data.frame(VarType=character(), stringsAsFactors=F)}
vars_removed_ann <- cbind(vars_removed_ann, vartype)
#vars_removed_ann$ID2 <- 1:nrow(vars_removed_ann)
vars_removed_ann_SNPS <- vars_removed_ann %>% filter(VarType == "SNP")
vars_removed_ann_indel <- vars_removed_ann %>% filter(VarType == "INDEL")

positions <- gatk %>% rowwise %>% do(AnnotateBEGIN_END(.$POS, .$REF, .$ALT))
if (invalid(positions)) {positions <- data.frame(BEGIN=integer(), END=integer())}
gatk_ann <- cbind(gatk, positions)
vartype <- gatk_ann %>% rowwise %>% do(AnnotateVarType(.$REF, .$ALT))
if (invalid(vartype)) {vartype <- data.frame(VarType=character(), stringsAsFactors=F)}
gatk_ann <- cbind(gatk_ann, vartype)
if (!invalid(gatk_ann)) {
  gatk_ann$ID2 <- 1:nrow(gatk_ann)
} else {
  gatk_ann$ID2 <- integer()
}

gatk_ann_SNPS <- gatk_ann %>% filter(VarType == "SNP")
gatk_ann_indel <- gatk_ann %>% filter(VarType == "INDEL")

################first do the SNPS
vcf_denovo_filtered_ann_SNPS$FreeBayesRemoveReason <- rep(NA, nrow(vcf_denovo_filtered_ann_SNPS))
free_bayes_SNPs <- rbind(vcf_denovo_filtered_ann_SNPS, vars_removed_ann_SNPS)
if (!invalid(free_bayes_SNPs)) {
  free_bayes_SNPs$ID2 <- 1:nrow(free_bayes_SNPs)
} else {
  free_bayes_SNPs$ID2 <- integer()
}
SNP_compare <- full_join(gatk_ann_SNPS,
                         free_bayes_SNPs,
                         by=c("X.CHROM", "POS", "REF", "ALT"))
#anything that is NA for VarType.y is a remove=NotCalled
SNP_compare$FreeBayesRemoveReason[is.na(SNP_compare$VarType.y)] <- rep("NotCalled", length(SNP_compare$FreeBayesRemoveReason[is.na(SNP_compare$VarType.y)]))

SNP_compare$FreeBayesCompare <- NA
#anything that is NA for uniqueVarID but NA for FreeBayesRemoveReason is an ADD
SNP_compare$FreeBayesCompare[is.na(SNP_compare$uniqueVarID) & is.na(SNP_compare$FreeBayesRemoveReason)] <- rep("ADD", length(SNP_compare$FreeBayesCompare[is.na(SNP_compare$uniqueVarID) & is.na(SNP_compare$FreeBayesRemoveReason)]))
#anything NA for uniqueVarID but not NA for FreeBayesRemove is something only found by FreeBayes but then removed
SNP_compare$FreeBayesCompare[is.na(SNP_compare$uniqueVarID) & !is.na(SNP_compare$FreeBayesRemoveReason)] <- rep("FreeBayesUniqueRemoved", length(SNP_compare$FreeBayesCompare[is.na(SNP_compare$uniqueVarID) & !is.na(SNP_compare$FreeBayesRemoveReason)]))
#anything not NA for uniqueVarID and not NA for FreeBayesRemoveReason was removed by FreeBayes
SNP_compare$FreeBayesCompare[!is.na(SNP_compare$uniqueVarID) & !is.na(SNP_compare$FreeBayesRemoveReason)] <- rep("Removed", length(SNP_compare$FreeBayesCompare[!is.na(SNP_compare$uniqueVarID) & !is.na(SNP_compare$FreeBayesRemoveReason)]))
#any remaining NA for FreeBayesCompare was kept the same
SNP_compare$FreeBayesCompare[is.na(SNP_compare$FreeBayesCompare)] <- rep("Keep", length(SNP_compare$FreeBayesCompare[is.na(SNP_compare$FreeBayesCompare)]))

################now do the Indels
#We will need this consolidating function later, to consolidate rows that consisted of one gatk variant and one freebayes variant, with slightly different pos/ref/alt, but overlapping indels. I consider those the same, so consolidate the two rows down to one.
ConsolidateRowsNA <- function(dfr) {
  #consolidates down to one row
  dfr_2 <- dfr[1, ]
  dfr_2[] <- lapply(dfr, 
                function(x) {
                  x <- unique(x)
                  if (all(is.na(x))) {return(x)}
                  x <- x[!is.na(x)]
                  if (length(x) != 1) {
                    stop("Problem with supplied dataframe. Must have only 1 unique, non-NA value.")
                  } else {
                    return(x)
                  }
                })
  return(dfr_2)
}

vcf_denovo_filtered_ann_indel$FreeBayesRemoveReason <- rep(NA, nrow(vcf_denovo_filtered_ann_indel))
free_bayes_indels <- rbind(vcf_denovo_filtered_ann_indel, vars_removed_ann_indel)
if (!invalid(free_bayes_indels)) {
  free_bayes_indels$ID2 <- 1:nrow(free_bayes_indels)
} else {
  free_bayes_indels$ID2 <- integer()
}

INDEL_compare <- full_join(gatk_ann_indel,
                           free_bayes_indels,
                           by=c("X.CHROM", "POS", "REF", "ALT")) %>% arrange(X.CHROM, POS)
if (nrow(INDEL_compare)>0) {
  INDEL_compare$FreeBayesCompare <- NA
  #any present X.FAM and present FreeBayesRemoveReason is a reject
  INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$FreeBayesRemoveReason)] <- rep("Removed", length(INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$FreeBayesRemoveReason)]))
  #any present X.FAM and present VarType.y but NA FreeBayesRemoveReason is an accept
  INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$VarType.y) & is.na(INDEL_compare$FreeBayesRemoveReason)] <- rep("Keep", length(INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$VarType.y) & is.na(INDEL_compare$FreeBayesRemoveReason)]))
} else {
  INDEL_compare <- INDEL_compare
}

#match the remaining indels by finding intersections.
#remove the matched ones first
INDEL_compare_matched <- INDEL_compare %>% filter(!is.na(FreeBayesCompare))
INDEL_compare_unmatched <- INDEL_compare %>% filter(is.na(FreeBayesCompare))

INDEL_compare_unmatched_gatk <- INDEL_compare_unmatched %>% filter(!is.na(X.FAM))
INDEL_compare_unmatched_fb <- INDEL_compare_unmatched %>% filter(is.na(X.FAM))

gatk_indel_ranges <- with(INDEL_compare_unmatched_gatk,
                          GRanges(X.CHROM,
                                  IRanges(BEGIN.x, END.x),
                                  ID2=ID2.x))
fb_indel_ranges <- with(INDEL_compare_unmatched_fb,
                        GRanges(X.CHROM,
                                IRanges(BEGIN.y - 2, END.y + 2),
                                ID2=ID2.y))

intersection_indel <- suppressWarnings(intersect(gatk_indel_ranges, fb_indel_ranges))

if (length(intersection_indel) == 0) {
  #No matching overlapping indels. The INDEL_compare object is complete.
  #any present X.FAM, absent ID2.y is a NotCalled
  INDEL_compare$FreeBayesRemoveReason[!is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$ID2.y)] <- rep("NotCalled", length(INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$ID2.y)]))
  INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$ID2.y)] <- rep("Removed", length(INDEL_compare$FreeBayesCompare[!is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$ID2.y)]))
  
  #any NA X.FAM and present FreeBayesRemove Reason is something only found by FreeBayes but then removed
  INDEL_compare$FreeBayesCompare[is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$FreeBayesRemoveReason)] <- rep("FreeBayesUniqueRemoved", length(INDEL_compare$FreeBayesCompare[is.na(INDEL_compare$X.FAM) & !is.na(INDEL_compare$FreeBayesRemoveReason)]))
  
  #any NA X.FAM and NA FreeBayesRemove Reason is something added by FreeBayes
  INDEL_compare$FreeBayesCompare[is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$FreeBayesRemoveReason)] <- rep("ADD", length(INDEL_compare$FreeBayesCompare[is.na(INDEL_compare$X.FAM) & is.na(INDEL_compare$FreeBayesRemoveReason)]))
} else {
  #Overlapping indels present. Take the overlaps as the same. Allow only 1 overlap for now.
  non_overlap <- NULL
  overlapping_rows <- NULL
  for (i in 1:nrow(INDEL_compare_unmatched_gatk)) {
    one_row <- INDEL_compare_unmatched_gatk[i, , drop=FALSE]
    gatk_indel_one <- with(one_row,
                           GRanges(X.CHROM,
                                   IRanges(BEGIN.x, END.x),
                                   ID2=ID2.x))
    intersection_one <- suppressWarnings(intersect(gatk_indel_one, fb_indel_ranges))
    if (length(intersection_one) == 0) {
      #this is not a row that has an overlapping indel
      non_overlap <- rbind(non_overlap, one_row)
    } else {
      #an overlap is present. take the first one.
      overlaps <- findOverlaps(gatk_indel_one, fb_indel_ranges)
      fb_ID2_overlap <- unlist(fb_indel_ranges[subjectHits(overlaps)]$ID2)[1]
      #consolidate the fb row with the gatk row
      fb_row <- INDEL_compare_unmatched_fb %>% filter(ID2.y %in% fb_ID2_overlap)
      fb_pos_ref_alt <- fb_row %>% select(POS, REF, ALT)
      names(fb_pos_ref_alt) <- paste(names(fb_pos_ref_alt), "fb", sep=".")
      one_row_pos_ref_alt <- one_row %>% select(POS, REF, ALT)
      
      fb_row <- fb_row %>% select(-POS, -REF, -ALT)
      one_row <- one_row %>% select(-POS, -REF, -ALT)
      consolidated_row <- ConsolidateRowsNA(rbind(one_row, fb_row))
      consolidated_row <- cbind(consolidated_row, one_row_pos_ref_alt, fb_pos_ref_alt)
      overlapping_rows <- rbind(overlapping_rows, consolidated_row)
      #Now remove the overlapping row from the set of fb possbilities
      INDEL_compare_unmatched_fb <- INDEL_compare_unmatched_fb %>% filter(!ID2.y %in% fb_ID2_overlap)
      fb_indel_ranges <- with(INDEL_compare_unmatched_fb,
                              GRanges(X.CHROM,
                                      IRanges(BEGIN.y - 2, END.y + 2),
                                      ID2=ID2.y))
    }
  }
  #now rbind together the non_overlap, the overlap, and the remaining freebayes rows
  INDEL_compare_unmatched_new <- plyr::rbind.fill(non_overlap, overlapping_rows, INDEL_compare_unmatched_fb)
  #Add the annotations
  #any present X.FAM, present POS.fb, and NA FreeBayesRemoveReason is a keep for FreeBayesCompare
  INDEL_compare_unmatched_new$FreeBayesCompare[!is.na(INDEL_compare_unmatched_new$X.FAM) & !is.na(INDEL_compare_unmatched_new$POS.fb) & is.na(INDEL_compare_unmatched_new$FreeBayesRemoveReason)] <- "Keep"
  #any present X.FAM, present POS.gb, and present FreeBayesRemoveReason a Remove for FreeBayesCompare
  INDEL_compare_unmatched_new$FreeBayesCompare[!is.na(INDEL_compare_unmatched_new$X.FAM) & !is.na(INDEL_compare_unmatched_new$POS.fb) & !is.na(INDEL_compare_unmatched_new$FreeBayesRemoveReason)] <- "Removed"
  #any present X.FAM, absent ID2.y is a NotCalled
  INDEL_compare_unmatched_new$FreeBayesRemoveReason[!is.na(INDEL_compare_unmatched_new$X.FAM) & is.na(INDEL_compare_unmatched_new$ID2.y)] <- "NotCalled"
  INDEL_compare_unmatched_new$FreeBayesCompare[!is.na(INDEL_compare_unmatched_new$X.FAM) & is.na(INDEL_compare_unmatched_new$ID2.y)] <- "Removed"
  #any NA X.FAM and present FreeBayesRemove Reason is something only found by FreeBayes but then removed
  INDEL_compare_unmatched_new$FreeBayesCompare[is.na(INDEL_compare_unmatched_new$X.FAM) & !is.na(INDEL_compare_unmatched_new$FreeBayesRemoveReason)] <- "FreeBayesUniqueRemoved"
  
  #any NA X.FAM and NA FreeBayesRemove Reason is something added by FreeBayes
  INDEL_compare_unmatched_new$FreeBayesCompare[is.na(INDEL_compare_unmatched_new$X.FAM) & is.na(INDEL_compare_unmatched_new$FreeBayesRemoveReason)] <- "ADD"
  
  #the new INDEL_compare is an rbind of INDEL_compare_matched, INDEL_compare_unmatched_new
  INDEL_compare <- plyr::rbind.fill(INDEL_compare_matched, INDEL_compare_unmatched_new)
}

#combine the INDEL and SNP compare files:
denovo_gatk_compare <- plyr::rbind.fill(SNP_compare, INDEL_compare) %>% arrange(X.CHROM, POS)
#deal with X.FAM and the column names
denovo_gatk_compare$X.FAM[is.na(denovo_gatk_compare$X.FAM)] <- proband_ID
denovo_gatk_compare <- denovo_gatk_compare %>% select(-starts_with("BEGIN"))
denovo_gatk_compare <- denovo_gatk_compare %>% select(-starts_with("END"))
denovo_gatk_compare <- denovo_gatk_compare %>% select(-starts_with("ID2"))

uniquefyVarType <- function(x, y) {
  #there is always at least one non-NA.
  z <- c(x, y)
  z <- z[!is.na(z)]
  return(z[1])
}
denovo_gatk_compare$VarType <- mapply(uniquefyVarType, denovo_gatk_compare$VarType.x, denovo_gatk_compare$VarType.y)
denovo_gatk_compare <- denovo_gatk_compare %>% select(-VarType.y)
denovo_gatk_compare <- denovo_gatk_compare %>% select(-VarType.x)
names(denovo_gatk_compare) <- gsub("(^.*)(\\.x$)", "\\1\\.gatk", names(denovo_gatk_compare))
names(denovo_gatk_compare) <- gsub("(^.*)(\\.y$)", "\\1\\.fb", names(denovo_gatk_compare))
#reorder the column names
x <- c("FreeBayesCompare", "FreeBayesRemoveReason")
y <- names(denovo_gatk_compare)[!names(denovo_gatk_compare) %in% x]
denovo_gatk_compare <- denovo_gatk_compare[, c(x, y)]

#Write the outputs.
message("Outputting data.")
#set the output filenames
output_denovo_vcf <- paste(proband_ID, ".fb.RM.trio.vcf", sep="")
output_denovo_filtered_vcf <- paste(proband_ID, ".fb.RM.trio.f.vcf", sep="")
output_vars_removed <- paste(proband_ID, ".vars_removed.txt", sep="")
output_gatk_compare <- paste(proband_ID, ".gatk.fb.compare.txt", sep="")

#output the headers
sink(file=output_denovo_vcf)
cat(paste(proband_header, collapse="\n"))
cat("\n")
sink()

sink(file=output_denovo_filtered_vcf)
cat(paste(proband_header, collapse="\n"))
cat("\n")
sink()

#output the data
write.table(vcf_denovo, file=output_denovo_vcf, 
            col.names=F, row.names=F, quote=F, sep="\t", append=T)

write.table(vcf_denovo_filtered, file=output_denovo_filtered_vcf, 
            col.names=F, row.names=F, quote=F, sep="\t", append=T)

write.table(vars_removed, file=output_vars_removed, 
            col.names=T, row.names=F, quote=F, sep="\t")

write.table(denovo_gatk_compare, file=output_gatk_compare, 
            col.names=T, row.names=F, quote=F, sep="\t")

#bgzip and tabix index the VCF files
cmd <- paste("bgzip", output_denovo_vcf)
system(cmd)
cmd <- paste("bgzip", output_denovo_filtered_vcf)
system(cmd)
cmd <- paste("tabix -p vcf", paste(output_denovo_vcf, ".gz", sep=""))
system(cmd)
cmd <- paste("tabix -p vcf", paste(output_denovo_filtered_vcf, ".gz", sep=""))
system(cmd)

#output simple stats
#count SNPs and Indels
Count_Vartype <- function(dfr, vartype="SNP") {
  ref <- nchar(dfr$REF)
  alt <- nchar(dfr$ALT)
  if (vartype != "SNP") {
    return(sum(ref > 1 | alt > 1))
  } else {
    return(sum(ref == 1 & alt ==1))
  }
}
denovo_gatk_compare_SNP <- denovo_gatk_compare %>% filter(VarType == "SNP")
denovo_gatk_compare_indel <- denovo_gatk_compare %>% filter(VarType == "INDEL")
stats_table <- data.frame(proband=proband_ID,
                          
                          before_filter=nrow(vcf_denovo),
                          before_SNP=Count_Vartype(vcf_denovo, vartype="SNP"),
                          before_indel=Count_Vartype(vcf_denovo, vartype="indel"),
                          
                          after_filter=nrow(vcf_denovo_filtered),
                          after_SNP=Count_Vartype(vcf_denovo_filtered, vartype="SNP"),
                          after_indel=Count_Vartype(vcf_denovo_filtered, vartype="indel"),
                          
                          gatk_total=nrow(gatk),
                          gatk_SNP=Count_Vartype(gatk, vartype="SNP"),
                          gatk_indel=Count_Vartype(gatk, vartype="indel"),
                          
                          total_kept=as.numeric(table(denovo_gatk_compare$FreeBayesCompare)["Keep"]),
                          total_removed=as.numeric(table(denovo_gatk_compare$FreeBayesCompare)["Removed"]),
                          total_ADD=as.numeric(table(denovo_gatk_compare$FreeBayesCompare)["ADD"]),
                          
                          total_kept_SNP=as.numeric(table(denovo_gatk_compare_SNP$FreeBayesCompare)["Keep"]),
                          total_removed_SNP=as.numeric(table(denovo_gatk_compare_SNP$FreeBayesCompare)["Removed"]),
                          total_ADD_SNP=as.numeric(table(denovo_gatk_compare_SNP$FreeBayesCompare)["ADD"]),
                          
                          total_kept_indel=as.numeric(table(denovo_gatk_compare_indel$FreeBayesCompare)["Keep"]),
                          total_removed_indel=as.numeric(table(denovo_gatk_compare_indel$FreeBayesCompare)["Removed"]),
                          total_ADD_indel=as.numeric(table(denovo_gatk_compare_indel$FreeBayesCompare)["ADD"]),
                          
                          stringsAsFactors=F)

stats_table[] <- lapply(stats_table, function(x) {x[is.na(x)] <- 0; return(x)})

output_stats <- paste(proband_ID, ".stats.txt", sep="")
write.table(stats_table, file=output_stats, 
            col.names=T, row.names=F, quote=F, sep="\t")







