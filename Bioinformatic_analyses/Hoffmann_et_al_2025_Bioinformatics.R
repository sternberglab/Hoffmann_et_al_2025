#!/usr/bin/env Rscript
#
# ============================================================================
# Automated gRNA Target Search for dCas12f Proteins
# ============================================================================
#
# Description:
#   This script performs automated computational identification and analysis 
#   of guide RNA (gRNA) sequences and their putative genomic targets for 
#   domesticated Cas12f (dCas12f) proteins. The analysis pipeline includes:
#   
#   1. Downloading and processing genomes encoding dCas12f homologs
#   2. Identifying gRNA scaffold sequences using covariance models (CM)
#   3. BLASTing guide sequences against their host genomes to find targets
#   4. Filtering putative targets to intergenic regions
#   5. Annotating target loci using Bakta and PFAM HMM searches
#   6. Identifying ncRNA elements associated with target loci
#   7. Categorizing targets by predicted function
#   8. TAM (Target Adjacent Motif) determination
#
# Dependencies:
#   R packages: tidyverse, ggplot2, Biostrings, biomartr, DECIPHER, rentrez,
#               rBLAST, GenomicRanges, ggsci, gggenes, ggbreak, RColorBrewer
#   
#   External tools: Infernal (cmsearch, cmbuild, cmcalibrate), BLAST+, 
#                   MMseqs2, MAFFT (linsi), LocARNA, HMMER3, Bakta, weblogo
#
# Usage:
#   1. Set the working directory to your project folder
#   2. Ensure all input files are in place (see INPUT FILES section below)
#   3. Run sections sequentially - each major section clears the environment
#   4. External command-line tools (commented out) should be run separately
#
# Author: [Author name]
# Date: [Date]
# Associated publication: [Citation]
#
# ============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(biomartr)
library(DECIPHER)
library(rentrez)
library(rBLAST)
library(GenomicRanges)
library(ggsci)
library(gggenes)
library(ggbreak)
library(RColorBrewer)

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR ENVIRONMENT
# ============================================================================

# Set your working directory
# setwd("/path/to/your/working/directory")

# Path to Pfam database (required for PFAM HMM searches)
PFAM_DB_PATH <- "/path/to/Pfam/Pfam-A.hmm"

# Path to Pfam clans file (download from https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/)
PFAM_CLANS_PATH <- "/path/to/Pfam/Pfam-A.clans.tsv"

# NCBI Entrez API key (optional but recommended for faster downloads)
# Get your API key at: https://www.ncbi.nlm.nih.gov/account/settings/
# set_entrez_key("YOUR_API_KEY_HERE")

###############################################################
# Automated Target Search for dCas12f proteins
###############################################################

rm(list = ls())

################################################################
# Build CM from scaffold sequences identified in RIP-seq data
################################################################

# Manually extracted scaffold sequences from RIP-seq reads for 8 loci
# were aligned and used to build a covariance model

# Command-line steps (run outside R):
# -----------------------------------------------------------------
# Align scaffold sequences:
#   mlocarna -tgtdir gRNA.RIP.locarna -keep-sequence-order \
#     -consensus-structure alifold -stockholm \
#     -iterate -v -threads=2 Scaffolds_from_RIP.fna
#
# Build covariance model:
#   cmbuild -F gRNA.RIP.cm gRNA.RIP.locarna/results/result.stk
#
# Calibrate CM:
#   cmcalibrate gRNA.RIP.cm
# -----------------------------------------------------------------


################################################################
# Download dCas12f-encoding genomes with complete assemblies
################################################################

# Import dCas12f homolog information
# INPUT: FASTA file of unique dCas12f protein sequences associated with RpoE
prot <- readAAStringSet("input_data/Unique_dCas12f_assoc_w_RpoE.faa")

# INPUT: TSV file with homolog metadata
df <- read_tsv("input_data/final_homologs.v2.tsv")

# INPUT: GenBank bacterial genome metadata
# This file can be downloaded from NCBI or generated from assembly_summary files
gb <- read_tsv("input_data/Bacterial_genome_metadata.genbank.tsv")

# Check how many assemblies are in the GenBank data
table(df$assembly %in% gb$assembly_accession)
gb <- gb %>% filter(assembly_accession %in% df$assembly)

# Check assembly quality levels
table(paste0(gb$assembly_level, " | ", gb$genome_rep))

# Filter to complete/chromosome-level assemblies only
# 106 are sequenced at "Full" and assembled at either "Chromosome" (n = 5) 
# or "Complete Genome" (n = 101) levels
gb <- gb %>% filter(assembly_level %in% c("Chromosome", "Complete Genome"), 
                    genome_rep == "Full")
df <- df %>% filter(assembly %in% gb$assembly_accession)

# Create output directory for genomes
dir.create("filtered_dCas12f_genomes/", showWarnings = FALSE)

# Get unique genomic accession numbers
gen.acc <- unique(df$gen.acc)

# Initialize empty FASTA file for genome sequences
gen <- DNAStringSet()
writeXStringSet(gen, "filtered_dCas12f_genomes/dCas12f.genomic.fna")

# Download genomes using Entrez in batches
batch.size <- 106
for(i in seq(1, length(gen.acc), by = batch.size)){
  
  # Set max index for this batch
  q <- i + (batch.size - 1)
  if(q > length(gen.acc)) { q <- length(gen.acc)}
  
  # Post accession numbers to NCBI history server
  search.token <- entrez_post(db = "nucleotide", id = gen.acc[i:q])
  Sys.sleep(1)
  
  if(!search.token$QueryKey == "Empty ID list; Nothing to store") {  
    
    # Fetch FASTA sequences
    v <- entrez_fetch(db = "nucleotide", rettype = "fasta", web_history = search.token)
    Sys.sleep(1)
    
    # Write to temporary file
    write(v, "temp.txt")
    
    # Import as DNAStringSet
    z <- readDNAStringSet("temp.txt")
    
    # Append to growing FASTA file and clean up
    writeXStringSet(z, "filtered_dCas12f_genomes/dCas12f.genomic.fna", append = TRUE)
    rm(z)
    system("rm temp.txt")
    
  }
  
  # Print progress
  message(paste0(round((q/length(gen.acc))*100, 2), "%..."))
  rm(search.token)
  
}

# Verify downloaded genomes
gen <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.genomic.fna")
names(gen) <- word(names(gen))
table(gen.acc %in% names(gen)) # Should show 106 = TRUE

# Export metadata files
write_tsv(df, "dCas12f.info.v1.tsv")
write_tsv(gb, "dCas12f.gen_assembly.info.tsv")


################################################################
# Search for gRNA sequences using covariance model
################################################################

rm(list = ls())

# Command-line step (run outside R):
# -----------------------------------------------------------------
# Search genomes with CM:
#   cmsearch --tblout gRNA.RIP.dCas12f.genomic.cm_out.tbl \
#     -A gRNA.RIP.dCas12f.genomic.cm_out.sto \
#     gRNA.RIP.cm filtered_dCas12f_genomes/dCas12f.genomic.fna \
#     > gRNA.RIP.dCas12f.genomic.cm_out.txt
# -----------------------------------------------------------------

# Import CMsearch results
colz <- c("hit", "rem1", "cm", "rm2", "type", "cm.start", "cm.end", 
          "hit.start", "hit.end", "strand", "trunc", "pass", "gc", 
          "bias", "score", "evalue", "inc", "description")
cm <- read_table("gRNA.RIP.dCas12f.genomic.cm_out.tbl", col_names = colz, comment = "#")

# Count hits that passed inclusion threshold per genome
cm.count <- cm %>% filter(inc == "!") %>% count(hit)
table(cm.count$n)
# Distribution of gRNA counts per genome:
#  1  2   3  4  5  6  7  8  18 
# 72  16  6  3  1  1  4  1   1 

barplot(table(cm.count$n), 
        xlab = "Number of gRNAs per genome", 
        ylab = "Number of genomes")


################################################################
# Extract guide sequences from CM hits
################################################################

# Import genomes
gen <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.genomic.fna")
names(gen) <- word(names(gen))

# Initialize columns for sequence extraction
cm$scaffold.nt <- NA
cm$guide.nt <- NA
cm$guide.start <- NA
cm$guide.end <- NA

# Guide length parameter (nucleotides downstream of scaffold)
guide.len <- 20

# Extract scaffold and guide sequences for each CM hit
for(i in 1:nrow(cm)){
  
  if(cm$strand[i] == "+") {
    
    # Extract sequences on plus strand
    cm$scaffold.nt[i] <- unname(as.character(narrow(
      start = cm$hit.start[i],
      end = cm$hit.end[i],
      gen[names(gen) == cm$hit[i]])))
    
    cm$guide.nt[i] <- unname(as.character(narrow(
      start = cm$hit.end[i] + 1,
      end = cm$hit.end[i] + guide.len,
      gen[names(gen) == cm$hit[i]])))
    
    cm$guide.start[i] <- cm$hit.end[i] + 1
    cm$guide.end[i] <- cm$hit.end[i] + guide.len
    
  } else {
    
    # Extract and reverse complement sequences on minus strand
    cm$scaffold.nt[i] <- unname(as.character(reverseComplement(narrow(
      start = cm$hit.end[i],
      end = cm$hit.start[i],
      gen[names(gen) == cm$hit[i]]))))
    
    cm$guide.nt[i] <- unname(as.character(reverseComplement(narrow(
      start = cm$hit.end[i] - guide.len,
      end = cm$hit.end[i] - 1,
      gen[names(gen) == cm$hit[i]]))))
    
    cm$guide.start[i] <- cm$hit.end[i] - guide.len
    cm$guide.end[i] <- cm$hit.end[i] - 1
    
  }
  
}

# Export unique guides for weblogo visualization
length(unique(cm$guide.nt)) # 176 unique sequences

x <- cm %>% filter(inc == "!") %>% select(guide.nt) %>% distinct()
guide.uni <- DNAStringSet(x$guide.nt)
names(guide.uni) <- paste0("A", 1:nrow(x))
writeXStringSet(guide.uni, "Unique_20nt_guides_for_weblogo.fna")

# Export results
write_tsv(cm, "dCas12f_putative_gRNA.info.v1.tsv")
write_tsv(cm.count, "dCas12f_putative_gRNA.counts.v1.tsv")


################################################################
# BLAST guides against their host genomes to find targets
################################################################

dir.create("Target_BLAST", showWarnings = FALSE)

rm(list = ls())

# Import gRNA info
df <- read_tsv("dCas12f_putative_gRNA.info.v1.tsv")

# Import genomes
gen <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.genomic.fna")
names(gen) <- word(names(gen))

# Add local IDs to each gRNA
# "A" prefix = passed inclusion threshold, "B" = did not pass
x <- df %>% filter(inc == "!")
y <- df %>% filter(!inc == "!")
df <- rbind(
  x %>% mutate(local.id = paste0("A", 1:nrow(x))), 
  y %>% mutate(local.id = paste0("B", 1:nrow(y)))
)
rm(x, y)

# BLAST parameters
gen.id <- names(gen)
guide.len <- 14  # Use first 14 nt of guide for BLAST search

# Iterate through each genome
for(i in 1:length(gen.id)){
  
  # Extract genome sequence
  this.gen <- gen[names(gen) == gen.id[i]]
  
  # Export genome as temporary BLAST database
  writeXStringSet(this.gen, "Target_BLAST/temp.db.fna")
  
  # Create FASTA of guide sequences from this genome
  x <- df %>% filter(hit == gen.id[i])
  this.guide <- narrow(DNAStringSet(x$guide.nt), start = 1, width = guide.len)
  names(this.guide) <- paste0(x$hit, "__", x$local.id)
  writeXStringSet(this.guide, "Target_BLAST/temp.query.fna")
  
  # Build BLAST database
  system("makeblastdb -in Target_BLAST/temp.db.fna -dbtype nucl", 
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Run BLASTn search with relaxed parameters for short sequences
  bl <- blast(db = "./Target_BLAST/temp.db.fna")
  cl <- predict(bl, this.guide, 
                custom_format = "qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch",
                BLAST_args = "-evalue 100 -word_size 7 -reward 1 -penalty -3 -gapopen 5 -gapextend 2")
  
  # Export results
  write_tsv(cl, paste0("Target_BLAST/", gen.id[i], ".blast.txt"))
  
  # Clean up temporary files
  system("rm Target_BLAST/temp*")
  message(paste0(i, " of ", length(gen.id)))
}

# Export updated gRNA info
write_tsv(df, "dCas12f_putative_gRNA.info.v2.tsv")


################################################################
# Extract information about putative targets
################################################################

rm(list = ls())

# Import gRNA info
df <- read_tsv("dCas12f_putative_gRNA.info.v2.tsv")

# Import genomes
gen <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.genomic.fna")
names(gen) <- word(names(gen))

# Import and combine all BLAST results
blast <- tibble()
filelist <- list.files("Target_BLAST")
for(i in 1:length(filelist)) { 
  blast <- rbind(blast, read_tsv(paste0("Target_BLAST/", filelist[i]))) 
}

# Initialize columns for target analysis
blast$guide <- NA
blast$target <- NA
blast$duplex <- NA
blast$duplex.mm <- NA
blast$TAM <- NA
blast$seed.match <- NA

df$blast.id <- paste0(df$hit, "__", df$local.id)

# Process each BLAST hit to extract target and TAM sequences
for(i in 1:nrow(blast)){
  
  # Get guide sequence
  blast$guide[i] <- str_sub(df$guide.nt[df$blast.id == blast$qseqid[i]], 
                            start = 1, end = blast$qlen[i])
  
  # Extract target sequence based on strand
  x <- gen[names(gen) == blast$sseqid[i]]
  
  if(blast$sstrand[i] == "plus") { 
    blast$target[i] <- unname(as.character(narrow(
      x, 
      start = blast$sstart[i] - (blast$qstart[i] - 1),
      end = blast$send[i] + (width(blast$guide[i]) - blast$qend[i])))) 
  }
  
  if(blast$sstrand[i] == "minus") { 
    blast$target[i] <- unname(as.character(reverseComplement(narrow(
      x, 
      start = blast$send[i] - (width(blast$guide[i]) - blast$qend[i]),
      end = blast$sstart[i] + (blast$qstart[i] - 1))))) 
  }
  
  # Determine guide-target duplex (match/mismatch pattern)
  this.guide <- str_split(blast$guide[i], "")[[1]]
  this.target <- str_split(blast$target[i], "")[[1]]
  x <- this.guide == this.target
  x[x == FALSE] <- "X"  # Mismatch
  x[x == TRUE] <- "|"   # Match
  blast$duplex[i] <- paste0(x, collapse = "")
  
  # Count mismatches
  blast$duplex.mm[i] <- sum(!this.guide == this.target)
  
  # Extract TAM (Target Adjacent Motif) - 5 nt upstream of target
  x <- gen[names(gen) == blast$sseqid[i]]
  
  if(blast$sstrand[i] == "plus") {
    this.start <- blast$sstart[i] - (blast$qstart[i] - 1) - 5
    if(this.start < 1) { this.start <- 1 }
    blast$TAM[i] <- unname(as.character(narrow(
      x, 
      start = this.start,
      end = blast$sstart[i] - (blast$qstart[i] - 1) - 1))) 
  }
  
  if(blast$sstrand[i] == "minus") { 
    this.end <- blast$sstart[i] + (blast$qstart[i] - 1) + 5
    if(this.end > width(x)) { this.end <- width(x) }
    blast$TAM[i] <- unname(as.character(reverseComplement(narrow(
      x, 
      start = blast$sstart[i] + (blast$qstart[i] - 1) + 1,
      end = this.end)))) 
  }
  
  # Check if seed region (first 6 nt) matches perfectly
  if(str_sub(blast$guide[i], 1, 6) == str_sub(blast$target[i], 1, 6)) { 
    blast$seed.match[i] <- TRUE 
  } else { 
    blast$seed.match[i] <- FALSE 
  }
  
}

# Summary of seed matches
table(blast$seed.match)
# FALSE  TRUE 
# 33682  9031

# Export results
write_tsv(df, "dCas12f_putative_gRNA.info.v3.tsv")
write_tsv(blast, "dCas12f.potential_targets.v1.tsv")


################################################################
# Filter targets to remove self-hits (guide targeting itself)
################################################################

rm(list = ls())

# Import gRNA and potential target info
df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v1.tsv")

# Keep only hits with perfect seed match
hits <- hits.all %>% filter(seed.match == TRUE) # 42,713 -> 9,031 hits

# Convert strand notation
hits$sstrand[hits$sstrand == "plus"] <- "+"
hits$sstrand[hits$sstrand == "minus"] <- "-"

# Identify hits that are actually the guide sequence itself
rm.these <- c()

for(i in 1:nrow(hits)){
  
  # Check if hit coordinates match guide coordinates
  if(hits$sstrand[i] == "+" & 
     nrow(df %>% filter(blast.id == hits$qseqid[i] & guide.start == hits$sstart[i])) == 1) { 
    rm.these <- c(rm.these, i) 
  }
  
  if(hits$sstrand[i] == "-" & 
     nrow(df %>% filter(blast.id == hits$qseqid[i] & guide.end == hits$sstart[i])) == 1) { 
    rm.these <- c(rm.these, i) 
  }
  
}  

# Verify we identified all guide self-hits
length(unique(hits$qseqid)) == length(rm.these) # Should be TRUE

# Remove self-hits
hits <- hits[-rm.these,] # 9,031 -> 8,818 hits

# Export filtered targets
write_tsv(hits, "dCas12f.potential_targets.v2.tsv")


################################################################
# Restrict putative targets to intergenic regions
################################################################

rm(list = ls())

# Import gRNA and potential target info
df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v2.tsv")

# Note: Genome annotation with Bakta should be performed externally
# Command-line steps (run outside R):
# -----------------------------------------------------------------
# Activate Bakta environment:
#   conda activate bakta
#   export BAKTA_DB=/path/to/bakta_db/db
#
# Annotate each genome (use accompanying R script for batch processing):
#   ./filtered_dCas12f_genomes/Annotate_dCas12f_genomes.R
# -----------------------------------------------------------------

# Calculate genomic coordinates for target regions
hits.all <- hits.all %>% mutate(guide.start = NA, guide.end = NA)

# Plus strand hits
hits.all$guide.start[hits.all$sstrand == "+"] <- hits.all$sstart[hits.all$sstrand == "+"]
hits.all$guide.end[hits.all$sstrand == "+"] <- hits.all$guide.start[hits.all$sstrand == "+"] + 13

# Minus strand hits
# Note: For minus strand, guide.start refers to the 3' end genomic position
# and guide.end refers to the 5' end genomic position
hits.all$guide.start[hits.all$sstrand == "-"] <- hits.all$sstart[hits.all$sstrand == "-"] - 13
hits.all$guide.end[hits.all$sstrand == "-"] <- hits.all$sstart[hits.all$sstrand == "-"]

# Filter targets to intergenic regions using GFF annotations
gen.acc <- unique(df$hit)
hits <- tibble()

for(i in 1:length(gen.acc)){
  
  # Extract hits for this genome
  these.hits <- hits.all %>% filter(sseqid == gen.acc[i])
  
  # Import GFF file from Bakta annotation
  annot <- suppressWarnings(read_gff(
    paste0("filtered_dCas12f_genomes/dCas12f_bakta/", gen.acc[i], "/temp.gen.gff3")))
  
  # Keep only CDS annotations for overlap checking
  annot <- annot %>% filter(type == "CDS")
  
  # Create IRanges objects for efficient overlap detection
  x.annot <- IRanges(start = annot$start, end = annot$end)
  x.hits <- IRanges(start = these.hits$guide.start, end = these.hits$guide.end)
  
  # Find overlaps between targets and CDS
  x.over <- findOverlaps(x.hits, x.annot)
  x.over <- unique(queryHits(x.over))
  x.non.over <- setdiff(seq_along(x.hits), x.over)
  
  # Keep only non-overlapping (intergenic) hits
  hits <- rbind(hits, these.hits[x.non.over,])
  
}

# Assign unique IDs to putative targets
hits <- hits %>% mutate(put_target.id = paste0("tar", 1:nrow(hits)))

# Export filtered targets (1,370 intergenic hits)
write_tsv(hits, "dCas12f.potential_targets.v3.tsv")


################################################################
# Extract putative target loci sequences with flanking regions
################################################################

rm(list = ls())

# Import gRNA and potential target info
df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v3.tsv")

# Import genomes
gen <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.genomic.fna")

# Define flank sizes for locus extraction
upstream.flank.size <- 500    # bp upstream of target
downstream.flank.size <- 10000 # bp downstream of target

# Initialize storage
loci <- DNAStringSet()
hits <- hits %>% transform(
  gen.length = NA,
  locus.start = NA,
  locus.end = NA,
  locus.length = NA,
  target.locus.start = NA,
  target.locus.end = NA
)

# Extract locus sequences for each putative target
for(i in 1:nrow(hits)){
  
  this.hit <- hits[i,]
  this.gen <- gen[word(names(gen)) == this.hit$sseqid]
  
  # Calculate locus boundaries based on strand
  if(this.hit$sstrand == "+") {
    locus.start <- this.hit$guide.start - upstream.flank.size
    locus.end <- this.hit$guide.end + downstream.flank.size 
  }
  
  if(this.hit$sstrand == "-") {
    locus.start <- this.hit$guide.start - downstream.flank.size
    locus.end <- this.hit$guide.end + upstream.flank.size 
  }
  
  # Adjust boundaries if they extend beyond contig ends
  shift.start <- 0
  shift.end <- 0
  
  if(locus.start < 1) { 
    old.start <- locus.start
    locus.start <- 1
    shift.start <- locus.start - old.start
    message(paste0("Locus ", i, " adjusted: start value < 1 [", this.hit$sstrand, " strand]"))
  }
  
  if(locus.end > width(this.gen)) { 
    old.end <- locus.end
    locus.end <- width(this.gen) 
    shift.end <- old.end - locus.end
    message(paste0("Locus ", i, " adjusted: end > contig length [", this.hit$sstrand, " strand]"))
  }
  
  # Extract locus sequence
  this.locus <- narrow(this.gen, start = locus.start, end = locus.end)
  
  # Reverse complement minus strand loci for consistent orientation
  if(this.hit$sstrand == "-") { 
    this.locus <- reverseComplement(this.locus) 
  }
  
  # Calculate relative target coordinates within locus
  if(this.hit$sstrand == "+") {
    rel.start <- (upstream.flank.size - shift.start) + 1
    rel.end <- (upstream.flank.size - shift.start) + 14
  }
  
  if(this.hit$sstrand == "-") { 
    rel.start <- (upstream.flank.size - shift.end) + 1
    rel.end <- (upstream.flank.size - shift.end) + 14
  }
  
  # Store locus
  names(this.locus) <- this.hit$put_target.id
  loci <- c(loci, this.locus)
  
  # Record metadata
  hits$gen.length[i] <- width(this.gen)
  hits$locus.start[i] <- locus.start
  hits$locus.end[i] <- locus.end
  hits$locus.length[i] <- width(this.locus)
  hits$target.locus.start[i] <- rel.start
  hits$target.locus.end[i] <- rel.end
  
}

# Export loci sequences and updated target info
write_tsv(hits, "dCas12f.potential_targets.v4.tsv")
writeXStringSet(loci, paste0("filtered_dCas12f_genomes/dCas12f.potential_targets.", 
                             upstream.flank.size, "_", downstream.flank.size, ".loci.fna"))


################################################################
# Cluster loci sequences to identify redundant targets
################################################################

rm(list = ls())

# Import data
df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v4.tsv")

# Command-line step (run outside R):
# -----------------------------------------------------------------
# Cluster loci at 50% nucleotide identity:
#   mmseqs easy-linclust dCas12f.potential_targets.500_10000.loci.fna \
#     dCas12f.potential_targets.500_10000.loci.mm \
#     tmp1 --min-seq-id 0.5 -c 0.5
# -----------------------------------------------------------------

# Import clustering results
cls <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.mm_cluster.tsv",
                col_names = c("rep.id", "put_target.id"))

# Assign cluster IDs
x <- cls %>% select(rep.id) %>% distinct() 
x <- x %>% mutate(cluster.id = paste0("clus", 1:nrow(x)))
cls <- cls %>% left_join(x)

# Count cluster sizes
cls.count <- cls %>% count(cluster.id)

# Visualize cluster size distribution
ggplot() + 
  geom_bar(data = cls.count, aes(x = n)) + 
  theme_bw() + 
  xlab("Cluster members") + 
  ylab("Number of clusters")

# Add cluster info to hits dataframe
hits <- hits %>% left_join(cls %>% select(put_target.id, cluster.id))


################################################################
# Extract and annotate ORFs from target loci
################################################################

rm(list = ls())

# Import data
df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v4.tsv")
loci <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.fna")

# Function to flip GFF coordinates for reverse complement sequences
flip_gff <- function(gff_df, genome_length) {
  
  if (!all(c("start", "end", "strand") %in% colnames(gff_df))) {
    stop("GFF data must contain 'start', 'end', and 'strand' columns.")
  }
  
  # Transform coordinates
  new_start <- genome_length - gff_df$end + 1
  new_end <- genome_length - gff_df$start + 1
  
  gff_df <- gff_df %>%
    dplyr::mutate(start = pmin(new_start, new_end), 
                  end = pmax(new_start, new_end))
  
  # Flip strand
  gff_df$strand <- ifelse(gff_df$strand == "+", "-", 
                          ifelse(gff_df$strand == "-", "+", gff_df$strand))
  
  return(gff_df)
}

# Initialize storage
loci.annot <- tibble()
orfs.in.loci <- AAStringSet()

# Process each target locus
for(i in 1:nrow(hits)){
  
  this.hit <- hits[i,]
  this.locus <- loci[names(loci) == this.hit$put_target.id]
  
  # Import Bakta annotation files
  annot <- suppressWarnings(read_gff(
    paste0("filtered_dCas12f_genomes/dCas12f_bakta/", this.hit$sseqid, "/temp.gen.gff3")))
  annot.cds <- readAAStringSet(
    paste0("filtered_dCas12f_genomes/dCas12f_bakta/", this.hit$sseqid, "/temp.gen.ffn"))
  annot.prot <- readAAStringSet(
    paste0("filtered_dCas12f_genomes/dCas12f_bakta/", this.hit$sseqid, "/temp.gen.faa"))
  
  # Filter annotations
  annot <- annot %>% 
    filter(!type == "region" & !is.na(type)) %>% 
    filter(!is.na(start)) %>% 
    filter(!is.na(end))
  annot$strand[annot$type == "CRISPR"] <- "+"
  annot$strand[!annot$strand %in% c("+", "-")] <- "+"
  
  # Add local IDs
  annot$annot.id <- paste0("annot", 1:nrow(annot))
  
  # Find annotations overlapping the locus
  x.annot <- IRanges(start = annot$start[!is.na(annot$start) & !is.na(annot$end)],
                     end = annot$end[!is.na(annot$start) & !is.na(annot$end)])
  names(x.annot) <- annot$annot.id
  
  x.locus <- IRanges(start = this.hit$locus.start, end = this.hit$locus.end)
  
  x.over <- findOverlaps(x.annot, x.locus)
  x.over <- unique(queryHits(x.over))
  overlapping <- x.annot[x.over]
  
  annot <- annot %>% filter(annot.id %in% names(overlapping))
  
  # Adjust coordinates based on strand
  if(this.hit$sstrand == "+") {
    annot$start <- annot$start - this.hit$locus.start
    annot$end <- annot$end - this.hit$locus.start
  }
  
  if(this.hit$sstrand == "-") { 
    # Add temporary locus row for coordinate transformation
    annot <- rbind(annot, tibble(
      seqid = "locus", source = NA, type = NA,
      start = this.hit$locus.start, end = this.hit$locus.end,
      score = NA, strand = "+", phase = NA,
      attribute = NA, annot.id = NA))
    
    annot <- flip_gff(annot, this.hit$gen.length)
    
    annot.adj <- annot$start[annot$seqid == "locus"]
    annot$start <- annot$start - annot.adj + 1
    annot$end <- annot$end - annot.adj + 1
    
    annot <- annot %>% filter(!seqid == "locus")
  }
  
  # Handle annotations at locus boundaries
  annot$attribute[annot$start < 1] <- paste0(annot$attribute[annot$start < 1], ";partial=10")
  annot$attribute[annot$start > 1] <- paste0(annot$attribute[annot$start > 1], ";partial=00")
  annot$start[annot$start < 1] <- 1
  
  annot$attribute[annot$end > width(this.locus)] <- paste0(
    annot$attribute[annot$end > width(this.locus)], ";partial=01")
  annot$end[annot$end > width(this.locus)] <- width(this.locus)
  
  # Add target annotation
  annot <- rbind(annot, tibble(
    seqid = "Target",
    source = "computational",
    type = "Target",
    start = this.hit$target.locus.start,
    end = this.hit$target.locus.end,
    score = ".",
    strand = "+",
    phase = ".",
    attribute = paste0("ID=", this.hit$put_target.id, ";Name=Putative target"),
    annot.id = this.hit$put_target.id))
  
  annot <- annot %>% arrange(start)
  annot$seqid <- this.hit$put_target.id
  
  # Extract Bakta gene IDs
  annot$bakta.id <- gsub("ID\\=", "", word(annot$attribute, sep = "\\;"))
  
  # Get sequences for annotated features
  annot.cds <- annot.cds[word(names(annot.cds)) %in% annot$bakta.id]
  annot.prot <- annot.prot[word(names(annot.prot)) %in% annot$bakta.id]
  
  x.cds <- tibble(
    bakta.id = word(names(annot.cds)),
    bakta.label = word(names(annot.cds), start = 2, end = -1),
    orf.nt = as.character(annot.cds))
  
  x.prot <- tibble(
    bakta.id = word(names(annot.prot)),
    orf.aa = as.character(annot.prot))
  
  x.cds <- x.cds %>% left_join(x.prot, by = "bakta.id")
  annot <- annot %>% left_join(x.cds, by = "bakta.id")
  
  if(length(x.prot) > 0) {
    x <- AAStringSet(x.prot$orf.aa)
    names(x) <- x.prot$bakta.id
    orfs.in.loci <- c(orfs.in.loci, x)
  }
  
  loci.annot <- rbind(loci.annot, annot)
  
  rm(list = setdiff(ls(), c("df", "hits", "loci", "flip_gff", "loci.annot", "orfs.in.loci", "i")))
  
  message(paste0(i, " of ", nrow(hits), "..."))
  
}

# Remove redundant ORF sequences
orf.df <- tibble(id = names(orfs.in.loci), seq = unname(as.character(orfs.in.loci))) %>% distinct()
orfs.in.loci <- AAStringSet(orf.df$seq)
names(orfs.in.loci) <- orf.df$id

# Export results
writeXStringSet(orfs.in.loci, "filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.ORFs.faa")
write_tsv(loci.annot, "filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v1.tsv")


################################################################
# Cluster translated ORFs from target loci
################################################################

rm(list = ls())

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v4.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v1.tsv")

# Command-line step (run outside R):
# -----------------------------------------------------------------
# Cluster translated ORFs at 40% amino acid identity:
#   mmseqs linclust dCas12f.potential_targets.500_10000.loci.ORFs.faa \
#     dCas12f.potential_targets.500_10000.loci.ORFs.mm \
#     tmp1 --min-seq-id 0.4 -c 0.5
# -----------------------------------------------------------------

# Import clustering results
cls <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.ORFs.mm_cluster.tsv",
                col_names = c("rep.id", "bakta.id"))

# Assign cluster IDs and count sizes
x <- cls %>% select(rep.id) %>% distinct() 
x <- x %>% mutate(cluster.id = paste0("clus", 1:nrow(x)))
cls <- cls %>% left_join(x)

cls.count <- cls %>% count(cluster.id)
cls <- cls %>% left_join(cls.count %>% select(cluster.id, cluster.size = n))

# Visualize cluster distribution
ggplot() + 
  geom_bar(data = cls.count, aes(x = n)) + 
  theme_bw() + 
  xlab("Cluster members") + 
  ylab("Number of clusters")

# Add cluster info to annotations
annot <- annot %>% left_join(cls %>% select(bakta.id, cluster.id, cluster.size))


################################################################
# Annotate loci with dCas12f HMM
################################################################

rm(list = ls())

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v4.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v1.tsv")

# Command-line steps (run outside R):
# -----------------------------------------------------------------
# Build HMM from dCas12f sequences:
#   hmmbuild dCas12f_sequences_in_HAT.unique.hmm \
#     dCas12f_sequences_in_HAT.unique.afa
#
# Search loci ORFs with dCas12f HMM:
#   hmmsearch --tblout dCas12f.potential_targets.500_10000.loci.ORFs.dCas12f_sequences_in_HAT.tbl \
#     --domtblout dCas12f.potential_targets.500_10000.loci.ORFs.dCas12f_sequences_in_HAT.domtbl \
#     dCas12f_sequences_in_HAT.unique.hmm \
#     filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.ORFs.faa \
#     > dCas12f.potential_targets.500_10000.loci.ORFs.dCas12f_sequences_in_HAT.txt
# -----------------------------------------------------------------

# Function to parse HMMER domain table output
parse_hmmscan_domtbl <- function(file_path, collapse_hits = FALSE) {
  
  raw_lines <- readLines(file_path)
  data_lines <- raw_lines[!grepl('^#', raw_lines)]
  
  parsed_lines <- lapply(data_lines, function(line) {
    fields <- strsplit(gsub(" +", " ", trimws(line)), " ")[[1]]
    if (length(fields) < 22) {
      fields <- c(fields, rep(NA, 22 - length(fields)))
    }
    description <- if (length(fields) > 22) paste(fields[23:length(fields)], collapse = " ") else ""
    c(fields[1:22], description)
  })
  
  domtbl_df <- do.call(rbind, parsed_lines)
  domtbl_df <- as.data.frame(domtbl_df, stringsAsFactors = FALSE)
  
  colnames(domtbl_df) <- c(
    "target_name", "accession", "tlen", "query_name", "query_accession", "qlen", 
    "E_value", "score", "bias", "dom_num", "of", "c_Evalue", "i_Evalue", 
    "dom_score", "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", 
    "env_from", "env_to", "acc", "description"
  )
  
  domtbl_df$env_from <- as.numeric(gsub(" .*$", "", domtbl_df$env_from))
  domtbl_df$env_to <- as.numeric(gsub("^.* ", "", domtbl_df$env_to))
  
  numeric_columns <- c("tlen", "qlen", "E_value", "score", "bias", "dom_num", "of", 
                       "c_Evalue", "i_Evalue", "dom_score", "dom_bias", "hmm_from", 
                       "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc")
  
  domtbl_df[numeric_columns] <- lapply(domtbl_df[numeric_columns], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  
  if (collapse_hits) {
    domtbl_df <- collapse_hmm_hits(domtbl_df)
  }
  
  return(domtbl_df)
}

# Function to collapse overlapping HMM hits
collapse_hmm_hits <- function(df) {
  
  required_columns <- c('query_name', 'target_name', 'i_Evalue', 'env_from', 'env_to')
  if (!all(required_columns %in% colnames(df))) {
    stop("Missing one or more required columns in the input dataframe")
  }
  
  collapse_for_query <- function(query_df) {
    
    query_df <- query_df %>% arrange(env_from, env_to)
    
    collapsed_hits <- list()
    current_region <- NULL
    
    for (i in seq_len(nrow(query_df))) {
      row <- query_df[i, ]
      
      if (is.null(current_region)) {
        current_region <- list(
          env_start = row$env_from,
          env_end = row$env_to,
          best_hit = row
        )
      } else {
        if (row$env_from <= current_region$env_end) { 
          current_region$env_end <- max(current_region$env_end, row$env_to) 
          
          if (row$i_Evalue < current_region$best_hit$i_Evalue) {
            current_region$best_hit <- row
          }
        } else {
          collapsed_hits <- append(collapsed_hits, list(current_region$best_hit))
          current_region <- list(
            env_start = row$env_from,
            env_end = row$env_to,
            best_hit = row
          )
        }
      }
    }
    
    if (!is.null(current_region)) {
      collapsed_hits <- append(collapsed_hits, list(current_region$best_hit))
    }
    
    return(do.call(rbind, collapsed_hits))
  }
  
  collapsed_df <- df %>% 
    group_by(query_name) %>% 
    group_split() %>% 
    lapply(collapse_for_query) %>% 
    bind_rows()
  
  return(collapsed_df)
}

# Import HMMsearch results
hmm <- parse_hmmscan_domtbl(
  "filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.ORFs.dCas12f_sequences_in_HAT.domtbl")

# Label dCas12f-encoding ORFs
annot$bakta.label[annot$bakta.id %in% hmm$target_name] <- "dCas12f [custom HMM]"

# Export updated annotations
write_tsv(annot, "filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v2.tsv")


################################################################
# Search for ncRNA elements in target loci
################################################################

dir.create("ncRNA_search", showWarnings = FALSE)

rm(list = ls())

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits <- read_tsv("dCas12f.potential_targets.v4.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v2.tsv")

# Command-line step (run outside R):
# -----------------------------------------------------------------
# Search target loci with nhmmer using reference ncRNA sequence:
#   nhmmer --dfamtblout ncRNA_search/Ata_HTH_RpoE_intergenic.dCas12f.potential_targets.500_10000.dfamtbl \
#     --tblout ncRNA_search/Ata_HTH_RpoE_intergenic.dCas12f.potential_targets.500_10000.tbl \
#     ncRNA_search/Ata_HTH_RpoE_intergenic.fna \
#     filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.fna \
#     > ncRNA_search/Ata_HTH_RpoE_intergenic.dCas12f.potential_targets.500_10000.txt
# -----------------------------------------------------------------

# Import nhmmer results
colz <- c("put_target.id", "rem1", "query", "rem2", "hmm.from", "hmm.to", 
          "ali.from", "ali.to", "env.from", "env.to", "sq.len", "strand", 
          "evalue", "score", "bias", "description")
hmm <- read_table("ncRNA_search/Ata_HTH_RpoE_intergenic.dCas12f.potential_targets.500_10000.tbl",
                  col_names = colz, comment = "#")

# Import loci for sequence extraction
loci <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.fna")

# Extract ncRNA sequences from hits
hmm$ncrna.nt <- NA

for(i in 1:nrow(hmm)) {
  this.locus <- loci[names(loci) == hmm$put_target.id[i]]
  hmm$ncrna.nt[i] <- unname(as.character(narrow(
    this.locus, 
    start = hmm$ali.from[i], 
    end = hmm$ali.to[i])))
}

# Export unique ncRNA sequences
ncrna.nt <- DNAStringSet(hmm$ncrna.nt)
names(ncrna.nt) <- hmm$put_target.id
writeXStringSet(ncrna.nt,
                "ncRNA_search/Ata_HTH_RpoE_intergenic.dCas12f.potential_targets.500_10000.unique.fna")

# Iterative ncRNA refinement continues in subsequent sections...
# (See full workflow in original manuscript methods)


################################################################
# PFAM annotation of target loci ORFs
################################################################

rm(list = ls())

# Command-line step (run outside R):
# -----------------------------------------------------------------
# Scan ORFs with PFAM:
#   hmmscan --tblout dCas12f.potential_targets.500_10000.loci.ORFs.PFAM.cutGA.tbl \
#     --domtblout dCas12f.potential_targets.500_10000.loci.ORFs.PFAM.cutGA.domtbl \
#     --cut_ga /path/to/Pfam/Pfam-A.hmm \
#     dCas12f.potential_targets.500_10000.loci.ORFs.faa \
#     > dCas12f.potential_targets.500_10000.loci.ORFs.PFAM.cutGA.txt
# -----------------------------------------------------------------

# Re-define parsing functions (see previous section)
# [parse_hmmscan_domtbl and collapse_hmm_hits functions defined above]

# Import PFAM results
# hmm <- parse_hmmscan_domtbl(
#   "filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.ORFs.PFAM.cutGA.domtbl",
#   collapse_hits = TRUE)

# Import PFAM Clan information
# Note: Download from https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
colz <- c("hmm", "clan", "clan.name", "hmm.name", "pfam.description")
# clans <- read_tsv(PFAM_CLANS_PATH, col_names = colz)

# The rest of the annotation pipeline continues with categorization of
# target loci based on PFAM domain content...


################################################################
# Visualization function for target loci
################################################################

# Color palette for gene categories
color_pal <- c(
  "HTH" = "#FFFFB3",
  "dCas12f" = "#FB8072",
  "Target" = "#FB8072",
  "ncRNA_RIP" = "#BC80BD",
  "gRNA" = "#FFED6F",
  "SusC" = "#BEBADA",
  "SusD" = "#80B1D3",
  "Peptidase" = "#FDB462",
  "Thioredoxin" = "#FCCDE5",
  "tRNA" = "#B3DE69",
  "rRNA" = "#8DD3C7",
  "ncRNA" = "#CCEBC5",
  "CDS" = "#D9D9D9",
  "regulatory_region" = "#D9D9D9",
  "gap" = "#D9D9D9",
  "CRISPR" = "#D9D9D9",
  "RpoE" = "#B6D4EF",
  "crispr-repeat" = "#D9D9D9",
  "crispr-spacer" = "#D9D9D9"
)

# Function to generate locus maps
graph.loci <- function(annotations_dataframe, output_file_prefix, annot.out = TRUE, notes.out = TRUE) {
  
  pdf(paste0(output_file_prefix, ".pdf"), onefile = TRUE, width = 13.333, height = 7.5)
  
  x1 <- annotations_dataframe
  x1$seqid <- factor(x1$seqid, levels = unique(x1$seqid))
  
  x.1 <- unique(x1$seqid)
  
  if(length(x.1) > 10) {
    
    start.ind <- seq(1, length(x.1), by = 10)
    end.ind <- c(start.ind[2:length(start.ind)] - 1, length(x.1))
    
    for(i in 1:length(start.ind)){
      
      message(paste0(start.ind[i]:end.ind[i], collapse = " "))
      
      g <- ggplot(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]] & 
                                         !type %in% c("Target", "Guide", "Scaffold")),
                  aes(xmin = start, xmax = end, y = seqid, forward = strand.true,
                      fill = gene.cat, label = gene.label)) +
        
        geom_feature_label(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]]) %>% 
                             filter(type == "Target"),
                           aes(x = (start + end) / 2, y = seqid, label = gene.label, forward = NA)) +
        
        geom_feature_label(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]]) %>% 
                             filter(type == "Guide"),
                           aes(x = (start + end) / 2, y = seqid, label = gene.label, forward = NA)) +
        
        geom_segment(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]]) %>% 
                       filter(type == "Target"),
                     aes(x = start, y = seqid, xend = end, yend = seqid), 
                     color = "orange", linewidth = 8) +
        
        geom_segment(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]]) %>% 
                       filter(type == "Guide"),
                     aes(x = start, y = seqid, xend = end, yend = seqid), 
                     color = "purple", linewidth = 8, alpha = 0.8) +
        
        geom_segment(data = x1 %>% filter(seqid %in% x.1[start.ind[i]:end.ind[i]]) %>% 
                       filter(type == "Scaffold"),
                     aes(x = start, y = seqid, xend = end, yend = seqid), 
                     color = "tan", linewidth = 8, alpha = 0.8) +
        
        geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
        geom_gene_label(align = "left") +
        facet_wrap(~ seqid, scales = "free", ncol = 1) +
        scale_fill_manual(values = color_pal) +
        theme_genes() +
        xlab("Position (bp)") +
        theme(axis.title.y = element_blank())
      
      print(g)
    }
    
  } else {
    # Single page for fewer loci
    g <- ggplot(data = x1 %>% filter(!type %in% c("Target", "Guide", "Scaffold")),
                aes(xmin = start, xmax = end, y = seqid, forward = strand.true,
                    fill = gene.cat, label = gene.label)) +
      
      geom_feature_label(data = x1 %>% filter(type == "Target"),
                         aes(x = (start + end) / 2, y = seqid, label = gene.label, forward = NA)) +
      
      geom_segment(data = x1 %>% filter(type == "Target"),
                   aes(x = start, y = seqid, xend = end, yend = seqid), 
                   color = "orange", linewidth = 8) +
      
      geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +
      facet_wrap(~ seqid, scales = "free", ncol = 1) +
      scale_fill_manual(values = color_pal) +
      theme_genes() +
      xlab("Position (bp)") +
      theme(axis.title.y = element_blank())
    
    print(g)
  }
  
  dev.off()
  
  # Export annotation tables if requested
  if(annot.out) {
    write_tsv(annotations_dataframe, paste0(output_file_prefix, ".annot.tsv"))
  }
  
  if(notes.out) {
    notes.df <- annotations_dataframe %>% 
      select(seqid) %>% 
      distinct() %>% 
      mutate(notes = "", keep = "")
    write_tsv(notes.df, paste0(output_file_prefix, ".notes.tsv"))
  }
  
}


################################################################
# Categorize targets by function
################################################################

rm(list = setdiff(ls(), c("graph.loci", "color_pal")))

# Four main functional categories identified for RpoE-associated dCas12f targets:
# Category 1: Regulation of transmembrane/lipoproteins (transporter systems)
# Category 2: Regulation of TCS components (transcriptional regulators)
# Category 3: Auto-regulation
# Category 4: Trans-regulation (multiple dCas12fs regulate each other)

# Import manual scoring notes (generated from visual inspection of locus maps)
# INPUT: Manual scoring file with columns: seqid, ata-like, notes
these <- read_tsv("Loci_maps/All_loci_w_ncRNA.cluster_reps_only.notes.man_score.txt")

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v5.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v7.graphing.tsv")
cls <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.mm95_cluster.tsv",
                col_names = c("rep.id", "put_target.id"))

these$seqid <- gsub("\\*", "", these$seqid)

# Function to count total targets including cluster members
total.these <- function(representative_targets){
  
  x.hits <- hits.all %>% filter(put_target.id %in% representative_targets)
  
  # Add cluster members
  x.hits <- rbind(x.hits, hits.all %>% 
                    filter(put_target.id %in% cls$put_target.id[cls$rep.id %in% representative_targets] &
                             !put_target.id %in% x.hits$put_target.id))
  
  count.targets <- length(unique(x.hits$put_target.id))
  counts.grna <- length(unique(x.hits$qseqid))
  counts.genome <- length(unique(x.hits$sseqid))
  
  message(paste0(count.targets, " targets in ", counts.grna, " gRNAs from ", counts.genome, " genomes"))
}

# Count targets per category
message("Category 1 (Ata-like transporters):")
total.these(these$seqid[these$`ata-like` == "x"])

message("Category 2 (TCS regulation):")
total.these(these$seqid[grep("TCS", these$notes)])

message("Category 3 (Auto-regulation):")
total.these(these$seqid[grep("auto-reg", these$notes)])

message("Category 4 (Trans-regulation):")
total.these(these$seqid[grep("trans-reg", these$notes)])

# Generate locus maps for each category
dir.create("Loci_maps/categorized_v1/", showWarnings = FALSE)

graph.loci(annotations_dataframe = annot %>% filter(seqid %in% these$seqid[these$`ata-like` == "x"]),
           output_file_prefix = "Loci_maps/categorized_v1/Cat1_Ata_like")

graph.loci(annotations_dataframe = annot %>% filter(seqid %in% these$seqid[grep("TCS", these$notes)]),
           output_file_prefix = "Loci_maps/categorized_v1/Cat2_TCS_reg")

graph.loci(annotations_dataframe = annot %>% filter(seqid %in% these$seqid[grep("auto-reg", these$notes)]),
           output_file_prefix = "Loci_maps/categorized_v1/Cat3_auto_reg")

graph.loci(annotations_dataframe = annot %>% filter(seqid %in% these$seqid[grep("trans-reg", these$notes)]),
           output_file_prefix = "Loci_maps/categorized_v1/Cat4_trans_reg")


################################################################
# Extract ncRNAs with flanking regions
################################################################

rm(list = ls())

dir.create("ncRNA_search/target_adjacent_ncRNA", showWarnings = FALSE)

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v5.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v7.graphing.tsv")
cls <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.mm95_cluster.tsv",
                col_names = c("rep.id", "put_target.id"))
loci <- readDNAStringSet("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.fna")

# Fix any coordinate issues
x <- annot$end[annot$start > annot$end]
annot$end[annot$start > annot$end] <- annot$start[annot$start > annot$end]
annot$start[annot$start > annot$end] <- x

# Extract ncRNA sequences
x <- annot %>% filter(ncrna == TRUE)
other.annot <- tibble()
x$ncRNA.flanks <- NA

for(i in 1:nrow(x)){
  
  this.locus <- loci[names(loci) == x$seqid[i]]
  
  # Define flanking region (250 bp each side)
  flank.start <- x$start[i] - 250
  flank.end <- x$end[i] + 250
  
  if(flank.start < 1) { flank.start <- 1 }
  if(flank.end > width(this.locus)) { flank.end <- width(this.locus) }
  
  # Get overlapping annotations
  these.annot <- annot %>% filter(seqid == x$seqid[i])
  
  x.annot <- IRanges(start = these.annot$start, end = these.annot$end)
  x.hits <- IRanges(start = flank.start, end = flank.end)
  
  x.over <- findOverlaps(x.annot, x.hits)
  x.over <- unique(queryHits(x.over))
  these.annot <- these.annot[x.over,]
  
  # Update coordinates
  these.annot$start <- these.annot$start - flank.start + 1
  these.annot$end <- these.annot$end - flank.start + 1
  these.annot$start[these.annot$start < 1] <- 1
  these.annot$end[these.annot$end > width(this.locus)] <- width(this.locus)
  
  other.annot <- rbind(other.annot, these.annot)
  
  # Extract sequence
  x$ncRNA.flanks[i] <- as.character(narrow(this.locus, start = flank.start, end = flank.end))
  
  # Reverse complement if on minus strand
  if(x$strand[i] == "-") { 
    x$ncRNA.flanks[i] <- as.character(reverseComplement(DNAString(x$ncRNA.flanks[i]))) 
  }
}

# Export ncRNA sequences
y <- DNAStringSet(x$ncRNA.flanks)
names(y) <- x$seqid
writeXStringSet(y, "ncRNA_search/target_adjacent_ncRNA/target_adjacent_ncRNA.fna")
writeXStringSet(unique(y), "ncRNA_search/target_adjacent_ncRNA/target_adjacent_ncRNA.unique.fna")

write_tsv(x, "ncRNA_search/target_adjacent_ncRNA/target_adjacent_ncRNA.tsv")
write_tsv(other.annot[,1:9], "ncRNA_search/target_adjacent_ncRNA/target_adjacent_ncRNA.annot.gff", col_names = FALSE)


################################################################
# Match gRNAs to accompanying dCas12f proteins
################################################################

rm(list = ls())

df <- read_tsv("dCas12f_putative_gRNA.info.v3.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v5.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v7.graphing.tsv")

# Add gRNA and genome info to annotations
annot <- annot %>% left_join(
  hits.all %>% select(seqid = put_target.id, gRNA.id = qseqid, gen.acc = sseqid) %>% distinct(), 
  by = "seqid")

# Import dCas12f protein metadata
# INPUT: TSV with dCas12f homolog information including coordinates and sequences
meta <- read_tsv("input_data/final_homologs.v2.tsv") %>% 
  filter(gen.acc %in% df$hit) %>% 
  distinct()

# Assign unique IDs to dCas12f proteins
meta$dcas12f_id <- paste0("dCas12f_", 1:nrow(meta))

# Match each gRNA to its nearest dCas12f
df$assoc.dcas12 <- NA
df$assoc.dcas12.dist <- NA

for(i in 1:nrow(df)){
  
  this.df <- df[i,]
  this.meta <- meta %>% filter(gen.acc == this.df$hit)
  
  if(nrow(this.meta) == 0) { next }
  if(nrow(this.meta) > 1) {
    # Find closest dCas12f by coordinate distance
    this.meta$gRNA.dist <- abs((this.meta$start + this.meta$end)/2 - 
                                 (this.df$hit.start + this.df$guide.end)/2)
    this.meta <- this.meta %>% filter(gRNA.dist == min(gRNA.dist))
  }
  
  df$assoc.dcas12[i] <- this.meta$dcas12f_id
  
  # Calculate distance between dCas12f end and gRNA start
  if(this.df$strand == "+" & this.df$hit.start > this.meta$end) {
    df$assoc.dcas12.dist[i] <- this.df$hit.start - this.meta$end - 1 
  }
  if(this.df$strand == "-" & this.meta$start > this.df$hit.start) { 
    df$assoc.dcas12.dist[i] <- this.meta$start - this.df$hit.start - 1 
  }
  if(this.df$strand == "+" & this.df$hit.start < this.meta$end) {
    df$assoc.dcas12.dist[i] <- this.df$hit.end - this.meta$start - 1 
  }
  if(this.df$strand == "-" & this.meta$start < this.df$hit.start) { 
    df$assoc.dcas12.dist[i] <- this.meta$end - this.df$hit.end - 1 
  }
}

# Filter to gRNAs within 5 kb of a dCas12f
table(abs(df$assoc.dcas12.dist) < 5000) # 200 = TRUE, 14 = FALSE
df <- df %>% filter(abs(df$assoc.dcas12.dist) < 5000)

# Import manual quality scores
# INPUT: Manual scoring file from locus map inspection
these3 <- read_tsv("Loci_maps/categorized_v2/Cat1_Ata_like.notes.man_score.txt") %>% 
  arrange(desc(score)) %>% 
  filter(score > 3)

# Filter annotations to high-confidence targets
x <- unique(annot$seqid[annot$ncrna == TRUE])
annot <- annot %>% filter(seqid %in% c(these3$seqid, x))

# Final target set
hits <- hits.all %>% 
  filter(qseqid %in% df$blast.id) %>% 
  filter(put_target.id %in% annot$seqid)

length(unique(df$assoc.dcas12[df$blast.id %in% hits$qseqid])) # Unique dCas12f proteins

# Export results
dir.create("TAM_determination", showWarnings = FALSE)

hits <- hits %>% left_join(
  df %>% select(qseqid = blast.id, assoc.dcas12, assoc.dcas12.dist), 
  by = "qseqid")

write_tsv(hits, "TAM_determination/High_confidence_targets.tsv")
write_tsv(meta, "TAM_determination/Accompanying_dCas12f_info.tsv")
write_tsv(df, "TAM_determination/gRNAs_accompanied_by_dCas12f.tsv")


################################################################
# Cross-reference dCas12f proteins with phylogenetic tree
################################################################

rm(list = ls())

hits <- read_tsv("TAM_determination/High_confidence_targets.tsv")
prot.df <- read_tsv("TAM_determination/Accompanying_dCas12f_info.tsv") %>% 
  filter(dcas12f_id %in% hits$assoc.dcas12)

# INPUT: Information on dCas12f proteins in phylogenetic subtree
tree <- read_tsv("input_data/Cas12f_Cas12_from_HAT_subtree2.info.v3.tsv")

# Create FASTA of dCas12f proteins with high-confidence targets
prot <- AAStringSet(prot.df$seq)
names(prot) <- prot.df$dcas12f_id

# BLAST dCas12f proteins against reference tree sequences
db.path <- "input_data/Cas12f_Cas12_from_HAT_subtree2.all.faa"

# Build BLAST database
system(paste0("makeblastdb -in ", db.path, " -dbtype prot"), 
       ignore.stdout = TRUE, ignore.stderr = TRUE)

bl <- blast(db = db.path, type = "blastp")
cl <- predict(bl, prot, 
              custom_format = "qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch")

# Keep best hit per query
cl <- cl %>% 
  group_by(qseqid) %>% 
  filter(evalue == min(evalue)) %>%  
  filter(pident == max(pident)) %>% 
  filter(length == max(length)) %>% 
  ungroup()

matches <- cl %>% 
  group_by(qseqid) %>% 
  summarize(match = dplyr::first(sseqid), pident = dplyr::first(pident))

# Count targets per dCas12f
x <- hits %>% count(assoc.dcas12)
matches <- matches %>% left_join(x %>% select(count = n, qseqid = assoc.dcas12))

# Assign colors for ITOL visualization
matches$color <- "#ED1C24"
matches$color[matches$count == 2] <- "#F9ED32"
matches$color[matches$count > 2] <- "#00A651"

write_tsv(matches, "TAM_determination/Distribution_of_dCas12f_w_high_targets.ITOL.tsv")


################################################################
# Generate TAM weblogos for each dCas12f protein
################################################################

tam.path <- "TAM_determination/TAM_logos/" 
dir.create(tam.path, showWarnings = FALSE)

# Sort by target count
matches <- matches %>% arrange(count)

# Generate weblogos (requires weblogo in PATH)
for(i in 1:nrow(matches)){
  
  these.hits <- hits %>% filter(assoc.dcas12 == matches$qseqid[i])
  
  x <- DNAStringSet(these.hits$TAM)
  names(x) <- 1:nrow(these.hits)
  writeXStringSet(x, "temp.fna")
  
  system(paste0("weblogo -F eps --errorbars no -c classic --logo-font Arial -A dna -f temp.fna -o ", 
                tam.path, matches$qseqid[i], "__", matches$match[i], "__", matches$count[i], "_targets.eps"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  system("rm temp.fna")
}


################################################################
# Analyze distance from targets to downstream start codons
################################################################

rm(list = ls())

hits <- read_tsv("TAM_determination/High_confidence_targets.tsv")
hits.all <- read_tsv("dCas12f.potential_targets.v5.tsv")
annot <- read_tsv("filtered_dCas12f_genomes/dCas12f.potential_targets.500_10000.loci.annot.v7.graphing.tsv")

hits$dist.to.ds.cds <- NA

for(i in 1:nrow(hits)){
  
  this.annot <- annot %>% filter(seqid == hits$put_target.id[i])
  
  # Find target annotation
  this.tar <- this.annot %>% filter(seqid == annot.id)
  if(nrow(this.tar) == 0) { this.tar <- this.annot %>% filter(type == "Target") }
  
  # Find downstream CDS
  this.annot <- this.annot %>% filter(type == "CDS", start > this.tar$end)
  
  # Calculate distance
  hits$dist.to.ds.cds[i] <- this.annot$start[1] - this.tar$start
}

# Visualize distance distribution
ggplot() +
  geom_histogram(data = hits, aes(x = dist.to.ds.cds), bins = 30) +
  theme_classic() +
  xlab("Distance from target to downstream CDS start (bp)") +
  ylab("Count")

mean(hits$dist.to.ds.cds, na.rm = TRUE)

# Export data
write_tsv(hits, "TAM_determination/High_confidence_targets_with_distances.tsv")


################################################################
# Analyze guide-target mismatch positions
################################################################

hits <- read_tsv("TAM_determination/High_confidence_targets.tsv")

# Count mismatches at each position
df <- tibble(pos = 1:14, mm = 0)

for(i in 1:nrow(df)){
  df$mm[i] <- sum(str_sub(hits$duplex, start = i, end = i) == "X")
}

df$mm.per <- round(df$mm / nrow(hits) * 100, 2)

# Visualize mismatch distribution
ggplot(df) + 
  geom_tile(aes(x = pos, y = 0, fill = mm.per), color = "black") +
  scale_fill_gradient2(low = "white", high = "red", limits = c(0, 100), name = "% mismatched") +
  coord_fixed(xlim = c(0.5, 14.5)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  scale_x_continuous(breaks = 1:14) +
  xlab("gRNA-target duplex position") +
  geom_text(aes(x = pos, y = 0, label = round(mm.per, 1)), size = 3)


################################################################
# Session information
################################################################

# Print session info for reproducibility
# sessionInfo()


# ============================================================================
# END OF SCRIPT
# ============================================================================
