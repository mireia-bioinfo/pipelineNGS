# pipelineNGS
`pipelineNGS` is a package for processing epigenomic high-throughput data, specifically histone mark ChIP-seq and ATAC-seq.

## Getting started

### Programs you need to have installed and included in your path:

As this package is a wrapper for some command line tools, you need to have this programs in your `$PATH`. If they are not in your `$PATH`, you can also provide the path to the binary files using the appropriate arguments. 

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Samtools](http://www.htslib.org/)
- [MACS2](https://github.com/taoliu/MACS)

### Files you need to have in your local machine:

Additionally, you will need to download reference files to perform the different steps in the pre-processing pipeline:

- Reference genome indexed with Bowtie2. See more information [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome). You can download the fasta files for the different genome builds from the [UCSC download site](https://hgdownload.soe.ucsc.edu/downloads.html#human).
- Chromosome sizes (`gen_sizes`). You can download this file also from the UCSC: [hg38](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) or [hg19](https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes)
- ENCODE Blacklist (`blacklist`). They can be downloaded from [here](https://sites.google.com/site/anshulkundaje/projects/blacklists).

## Installing pipelineNGS in your local machine

Open your R session, install the `devtools` package if it is not already in your machine and type the following:

```
# Install pipelineNGS package
devtools::install_github("mireia-bioinfo/pipelineNGS")

# Load pipelineNGS package
library(pipelineNGS)
```

## Quick start

### Running pipelineNGS with H3K27ac ChIP-seq data

#### Single end

Create a file called `process_samples.R` containing the following code (change `your_fastq_folder` and the suffix of your files to resemble your data):

```
library(pipelineNGS)

## List of fastq files to analyze
fastq_files <- list.files("your_fastq_folder/", pattern="fastq.gz", full.names=TRUE)

## Name of the output files
output_names <- gsub(".fastq.gz", "", basename(fastq_files))

## Proces fastq files
process_epigenome(fastq_files = fastq_files,
                  run_fastqc = FALSE,
                  seq_type = "CHIP",
                  type = "SE",
                  out_name = output_names,
                  index = "/biodata/indices/species/Hsapiens/hg38",
                  blacklist = "/imppc/labs/lplab/share/reference/blacklists/hg38-blacklist.v2.bed",
                  cores=10)
                  
## Get and save stats
stats <- getStats(raw_bam=list.files("BAM/", pattern=".raw.bam$", full.names=TRUE),
                  path_logs = "Logs/",
                  path_peaks = "Peaks",
                  peak_type = "broadPeak") 
save(stats, file="Logs/stats_samples.rda")
```

From the terminal run the following commands: 

```
$ Rscript process_samples.R
```

#### Paired end

```
library(pipelineNGS)

## List of fastq files to analyze
fastq_files <- list.files("your_fastq_folder", pattern="fastq.gz", full.names=TRUE)

## Group files from the same sample
names <- gsub("_read[[:digit:]].fastq.gz", "", basename(fastq_files)) # extract basename
fastq_files <- split(fastq_files, names) # create list

## Name of the output files
output_names <- unique(names)

## Proces fastq files
process_epigenome(fastq_files = fastq_files,
                  run_fastqc = FALSE,
                  seq_type = "CHIP",
                  type = "PE",
                  out_name = output_names,
                  index = "/biodata/indices/species/Hsapiens/hg38",
                  cores=10)
                  
## Get and save stats
stats <- getStats(raw_bam=list.files("BAM/", pattern=".raw.bam$", full.names=TRUE),
                  path_logs = "Logs/",
                  path_peaks = "Peaks",
                  peak_type = "broadPeak") 
save(stats, file="Logs/stats_samples.rda")
```

From the terminal run the following commands: 

```
$ Rscript process_samples.R
```

### Running pipelineNGS with ATAC-seq data

Create a file called `process_samples.R` containing the following code (change `your_fastq_folder` and the suffix of your files to resemble your data):

#### Single end

```
library(pipelineNGS)

## List of fastq files to analyze
fastq_files <- list.files("your_fastq_folder/", pattern="fastq.gz", full.names=TRUE)

## Name of the output files
output_names <- gsub(".fastq.gz", "", basename(fastq_files))

## Proces fastq files
process_epigenome(fastq_files = fastq_files,
                  run_fastqc = FALSE,
                  seq_type = "ATAC",
                  type = "SE",
                  out_name = output_names,
                  index = "/biodata/indices/species/Hsapiens/hg38",
                  blacklist = "/imppc/labs/lplab/share/reference/blacklists/hg38-blacklist.v2.bed",
                  gen_sizes = "/imppc/labs/lplab/share/reference/chrom_sizes/hg38.len",
                  cores=10)

## Get and save stats
stats <- getStats(raw_bam=list.files("BAM/", pattern=".raw.bam$", full.names=TRUE),
                  path_logs = "Logs/",
                  path_peaks = "Peaks",
                  peak_type = "narrowPeak") 
save(stats, file="Logs/stats_samples.rda")
```

From the terminal run the following commands: 

```
$ Rscript process_samples.R
```

#### Paired end


Create a file called `process_samples.R` containing the following code (change `your_fastq_folder` and the suffix of your files to resemble your data):

#### Single end

```
library(pipelineNGS)

## List of fastq files to analyze
fastq_files <- list.files("your_fastq_folder", pattern="fastq.gz", full.names=TRUE)

## Group files from the same sample
names <- gsub("_read[[:digit:]].fastq.gz", "", basename(fastq_files)) # extract basename
fastq_files <- split(fastq_files, names) # create list

## Name of the output files
output_names <- unique(names)

## Proces fastq files
process_epigenome(fastq_files = fastq_files,
                  run_fastqc = FALSE,
                  seq_type = "ATAC",
                  type = "PE",
                  out_name = output_names,
                  index = "/biodata/indices/species/Hsapiens/hg38",
                  blacklist = "/imppc/labs/lplab/share/reference/blacklists/hg38-blacklist.v2.bed",
                  gen_sizes = "/imppc/labs/lplab/share/reference/chrom_sizes/hg38.len",
                  cores=10)

## Get and save stats
stats <- getStats(raw_bam=list.files("BAM/", pattern=".raw.bam$", full.names=TRUE),
                  path_logs = "Logs/",
                  path_peaks = "Peaks",
                  peak_type = "narrowPeak") 
save(stats, file="Logs/stats_samples.rda")
```

From the terminal run the following commands: 

```
$ Rscript process_samples.R
```
