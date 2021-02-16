#' Automated processing of CUT&TAG samples
#'
#' @param fastq_files List containing the different FastQ file pairs (one list element for each sample).
#' @param out_name Sample name to use for each of the samples.
#' @param type Character indicating if reads are paired end (PE) or single end (SE).
#' @param run_fastqc Logical indicating whether to run FastQC or not.
#' @param path_fastqc Path for the FastQC results.
#' @param path_bam Path for the aligment and filtering results.
#' @param path_peaks Path for the peak files.
#' @param path_logs Path for the logs.
#' @param index Path to the index needed by Bowtie2.
#' @param remove List of chromosomes to remove from final BAM file.
#' @param blacklist List of blaclisted regions to remove from BAM file.
#' @param chr_sizes Chromosome sizes to use for bedgraph conversion.
#' @param cores Number of threads to use for the analysis.
#' @param bedtools_bamtobed Path or alias of the bedtools bamtobed utility.
#' @param bedtools_genomecov Path or alias of the bedtools genomecov utility.
#' @param seacr Path or alias of the SEACR utilty.
#' @export
processCUTTAG <- function(fastq_files,
                          out_name,
                          type = "PE",
                          run_fastqc=TRUE,
                          # Directories
                          path_fastqc="FastQC/",
                          path_bam="BAM/",
                          path_peaks="Peaks/",
                          path_logs="Logs/",
                          # 1) Alignment params
                          index="/biodata/indices/species/Hsapiens/ucsc.hg19",
                          # 2) FiltOut params
                          remove=c("chrM", "chrUn", "_random", "_hap", "_gl", "EBV"),
                          blacklist="~/data/consensusBlacklist.bed",
                          # 3) Peak calling
                          chr_sizes = "",
                          # 4) Other params
                          cores=8,
                          # Program paths
                          bedtools_bamtobed = "bedtools bamtobed",
                          bedtools_genomecov = "bedtools genomecov",
                          seacr = "SEACR_1.3.sh") {

  if (type=="PE" & !is(fastq_files, "list")) fastq_files <- list(fastq_files)

  ## 0) Create directory tree -----------------------------
  dir.create(path_fastqc, F)
  dir.create(path_bam, F)
  dir.create(path_peaks, F)
  dir.create(path_logs, F)

  ## 1) Quality control with FastQC -----------------------
  if (run_fastqc) {
    fastqc(files=unlist(fastq_files),
           out_dir=path_fastqc)
  }

  ## 2) Alignment with Bowtie2 ----------------------------
  files <- mapply(alignmentBowtie2,
                  file=fastq_files,
                  out_name=out_name,
                  MoreArgs = list(type=type,
                                  index=index,
                                  out_dir=path_bam,
                                  path_logs=path_logs,
                                  cores=cores,
                                  "--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700")
  )

  ## 3) Post-processing with Samtools ----------------------
  files <- lapply(files,
                  filtOutBAM,
                  path_logs=path_logs,
                  type=type,
                  remove=remove,
                  blacklist=blacklist,
                  cores=cores)

  ## 4) Peak calling with SEACR ----------------------------
  files <- lapply(files,
                  peakCallingSEACR,
                  path_peaks = path_peaks,
                  chr_sizes = chr_sizes,
                  cores = cores,
                  bedtools_bamtobed = bedtools_bamtobed,
                  bedtools_genomecov = bedtools_genomecov,
                  seacr = seacr)
}


#' CUT&TAG peak calling with SEACR
#'
#' @inheritParams processCUTTAG
peakCallingSEACR <- function(bam_file,
                             path_peaks,
                             chr_sizes,
                             cores = 5,
                             bedtools_bamtobed = "bedtools bamtobed",
                             bedtools_genomecov = "bedtools genomecov",
                             seacr = "SEACR_1.3.sh") {
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Starting Peak Calling with SEACR"))
  name <- getNameFromPath(bam_file, suffix=".bam")
  bam_dir <- dirname(bam_file)

  ## 1) Convert bam files to bedgraph
  cmd <- paste("samtools sort -n", bam_file, "-@", cores-1, "-m 2G", "-o -",
               "|", bedtools_bamtobed, "-bedpe -i stdin",
               "| awk '$1==$4 && $6-$2 < 1000 {print $0}'",
               "| cut -f 1,2,6",
               "| sort -dsk1,1 -k2n,2 -k3nr,3 >", file.path(bam_dir, paste0(name, ".bed")))
  system(cmd)
  cmd <- paste(bedtools_genomecov, "-bg -i", file.path(bam_dir, paste0(name, ".bed")),
               "-g", chr_sizes,
               ">", file.path(bam_dir, paste0(name, ".bedgraph")))
  system(cmd)

  ## Remove temporary files
  file.remove(file.path(bam_dir, paste0(name, ".bed")))

  ## 2) Call peaks with SEACR
  cmd <- paste(seacr,
               file.path(bam_dir, paste0(name, ".bedgraph")),
               "0.01 norm stringent",
               file.path(path_peaks, paste0(name, ".peaks")))
  message(paste("\t", cmd))
  system(cmd)
}
