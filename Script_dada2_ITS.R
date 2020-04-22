library(dada2)
library(ShortRead)
library(Biostrings)

path <- "~/Documentos/Bioinfo_analise/dada2_r/trimmed"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#identificar primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#filtrar bases ambiguas (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#contar numero de vezes que os primers aparecem 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#remover os primers

cutadapt <- "/home/ju/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


#DADA2####
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Quality
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[3:4])


#Filter####

#deignar nomes

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)
out

out2_5 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 5), 
                        truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  


errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)


#Dereplicate identical reads ----
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference ----
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
print(dadaFs)

#inspecionar 
dadaFs[[1]]
head(getSequences(dadaFs[[1]]))

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#view mergers
seqtab <- makeSequenceTable(mergers)
## The sequences being tabled vary in length.
dim(seqtab)


# chimeras 
merger1 <- mergers
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)


#Construct Sequence Table ----
seqtab_filter_quimeras <- makeSequenceTable(merger1.nochim)
dim(seqtab_filter_quimeras)

#seqtab <- makeSequenceTable(mergers)
#dim(seqtab)

table(nchar(getSequences(seqtab_filter_quimeras)))

hist(nchar(getSequences(seqtab_filter_quimeras)), main="Distribution of sequence lengths")
#inspecionar numero de leituras

getN <- function(x) sum(getUniques(x))
track <- cbind(out2_5, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                          getN), rowSums(seqtab_filter_quimeras))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)
track
write.xlsx(track, "./track_2_5.xlsx")


#Atribuir taxonomia ####
unite.ref <- "./sh_general_release_dynamic_04.02.2020.fasta"
taxa <- assignTaxonomy(seqtab_filter_quimeras, unite.ref, multithread = TRUE, tryRC = TRUE) 

write.xlsx(taxa,"./taxa.xlsx")

#inspecionando 
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.xlsx(taxa.print,"./taxa.print.xlsx")

### Phyloseq ----

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library("readr") 
library("tidyverse")
library(microbiome); packageVersion("microbiome") 
library("readxl") 
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}
#importar dados
samples.out <- rownames(seqtab_filter_quimeras)
metadata <- read_tsv("metadado.tsv")
taxa <- taxa.print

seqtab.nochim <- seqtab_filter_quimeras


#extraindo dados do dada2 -----
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa.print
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


#(Re) -Lendo nossos dados

write.table(asv_tab, "ASVs_counts.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax, "ASVs_taxonomy.tsv",
            sep="\t", quote=F, col.names=NA)


count_tab <- read.table("ASVs_counts2.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")[ ]

tax_tab <- as.matrix(read.table("ASVs_taxonomy2.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("metadado_3.tsv", header=T,
                              sep="\t")


