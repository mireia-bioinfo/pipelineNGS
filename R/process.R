#' Automated processing of ChIP-seq and ATAC-seq samples
#'
#' This function performs all necessary steps in the ChIP-seq processing pipeline.
#' @param fastq_files Character string (single-end) or character vector of
#' length 2 (paired-end) with the file names of the samples to be analysed.
#' @param out_name Character vector, with the same length as \code{fastq_files},
#'  indicating the output filenames.
#' @param seq_type Experiment type, either "ATAC" (default) or "CHIP".
#' @param type Sequence type, one of "SE" (single end) or "PE" (paired end).
#' @param suffix_fastq Character indcating the suffix for the fastq files, for
#' example ".fq.gz" or ".fastq.gz".
#' @param path_fastqc Character indicating the output directory for the FastQC reports.
#' @param path_bam Character indicating the output directory for the bam files.
#' @param path_peaks Character indicating the output directory for the peak files.
#' @param path_logs Character indicating the output directory for the logs.
#' @param run_fastqc Logical indcating whether to run (TRUE) or not (FALSE) FastQC. Default: TRUE.
#' @param index Character indicating the location and basename for the
#' [Bowtie2 index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).
#' @param extra_bowtie2 Character containing additional arguments to be passed to
#' bowtie2 alignment call.
#' @param remove Character vector with chr that will be filtered out. Any chromosome
#' name containing matches for these characters will be removed.
#' @param blacklist Character indicating the file containing blacklist regions in
#' bed format. Any reads overlapping these regions will be discarded.
#' @param type_peak Character indicating the type of peak to be called with
#' MACS2, either "narrow" or "broad".
#' @param shift Logical indicating whether the reads should be shifted -100bp
#' and extended to 200bp (TRUE) or not (FALSE, default).
#' @param cores Number of threads to use for the analysis.
#' @param gen_sizes Character string indicating the path where the file with
#' chromosome name and sizes can be found. This argument is necessary only when `type="SE"`.
#' @param chunk Size of the chunk to load into memory for ATAC-seq read offset.
#' This argument is necessary only when `type="SE"`.
#' @export
#' @return Creates the folders `path_fastqc`, `path_bam`, `path_peaks`, `path_logs`,
#' by default in your working directory, containing the output files from de different
#' analyses.
#' @details This function ocesses ATAC-seq or ChIP-seq from FastQ files using
#' the following pipeline:
#' 1. Quality Control (FastQC).
#' 2. Alignment to reference genome (Bowtie2).
#' 3. Post-processing (Samtools), including removing duplicates, blacklisted
#' regions and non-reference chromosomes.
#' 4. (only for ATAC-seq) Offset correction (Samtools).
#' 5. Peak calling (MACS2).
#'
#' This function can process paired and single end FastQ files:
#' -  Single end files. The argument `fastq_files` should be a character vector with the name
#' of each file.
#' - Paired end files. The argument `fastq_files` should be a list, where each element is a
#' vector of size 1, where the first one is the R1 and the second one is the R2.
#'
#' @examples
#' \dontrun{
#' process_epigenome(fastq_files=c("path/to/file.fastq.gz", "path/to/file2.fastq.gz"),
#'                   seq_type="ATAC",
#'                   out_name=c("sample1", "sample2"),
#'                   type="SE",
#'                   cores=8)
#' }
#'
process_epigenome <- function(fastq_files,
                              out_name=NULL,
                              seq_type=c("ATAC", "CHIP", "CT"),
                              type="SE", # or PE
                              suffix_fastq=NULL,
                              cores=8,
                              # Directories
                              path_fastqc="FastQC/",
                              path_bam="BAM/",
                              path_peaks="Peaks/",
                              path_logs="Logs/",
                              # 0) FastQC
                              run_fastqc=TRUE,
                              # 1) Alignment params
                              index="/vault/refs/indexes/hg38",
                              extra_bowtie2="",
                              # 2) FiltOut params
                              remove=c("chrM", "chrUn", "_random", "_hap", "_gl", "EBVls"),
                              blacklist="/vault/refs/hg38-blacklist.v2.bed",
                              # 4) Peak calling params
                              type_peak=c("narrow", "broad"),
                              shift=c(TRUE, FALSE),
                              # 5) ATAC-seq offset
                              chunk=1e7,
                              gen_sizes="/vault/refs/hg38.chromSizes.txt") {

  if(type=="PE" & !is(fastq_files, "list")) stop("Character vectors not allowed as fastq_files for PE sequencing. Your input should be a a list of character vectors.")
  if(type=="PE" & unique(sapply(fastq_files, length))!=2) stop("Your list contains more than two files per element. In PE mode you should have a list, where each element contains R1 and R2.")

  # if (type=="PE" & !is(fastq_files, "list")) fastq_files <- list(fastq_files)

  ## 0) Create directory tree -----------------------------
  dir.create(path_fastqc, F)
  dir.create(path_bam, F)
  dir.create(path_peaks, F)
  dir.create(path_logs, F)

  ## 0) Update necessary parameters
  param_options <- paramdata[paramdata$Parameter==seq_type,]
  if (extra_bowtie2!="") param_options$Extra_bowtie2 <- extra_bowtie2
  if(length(type_peak)==1) param_options$peak_type <- type_peak
  if(length(shift)==1) param_options$shift <- shift

  ## 1) Quality control with fastqc -----------------------
  if (run_fastqc) fastqc(files = unlist(fastq_files),
                         path_fastqc = path_fastqc,
                         cores = cores)

  ## 2) Alignment with Bowtie2 ----------------------------
    files <- mapply(alignmentBowtie2,
                    file=fastq_files,
                    out_name=out_name,
                    MoreArgs = list(suffix_fastq=suffix_fastq,
                                    type=type,
                                    index=index,
                                    path_bam=path_bam,
                                    path_logs=path_logs,
                                    cores=cores,
                                    extra_bowtie2=param_options$Extra_bowtie2)
                    )
    files <- unlist(files)


  ## 3) Post-processing with Samtools ----------------------
  files <- lapply(files,
                  filtOutBAM,
                  path_logs=path_logs,
                  type = type,
                  remove=remove,
                  blacklist=blacklist,
                  cores=cores)

  ## 4) Offset correction for ATAC-seq ---------------------
  if (param_options$Offset_correction & type=="SE") {
    files <- lapply(files,
                    offsetATACSE,
                    chunk=chunk)
  } else if (param_options$Offset_correction & type=="PE") {
    files <- lapply(files,
                    offsetATAC,
                    gen_sizes=gen_sizes,
                    cores=cores)
  }

  ## 5) Peak calling with MACS2 ----------------------------
  files <- lapply(files,
                  callPeak,
                  path_peaks=path_peaks,
                  path_logs=path_logs,
                  type_peak=param_options$peak_type,
                  shift=param_options$shift,
                  type=type)
}
