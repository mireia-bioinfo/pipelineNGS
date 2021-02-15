#' Automated processing of ChIP-seq and ATAC-seq samples
#'
#' This function performs all necessary steps in the ChIP-seq processing pipeline.
#' @param fastq_files Path and filenames for the fastq files to be analyzed.
#' @param seq_type Type of sequencing, either ATAC (default) or CHIP.
#' @param path_fastqc Path for the FastQC output.
#' @param path_bam Path for the output BAM files.
#' @param chunk Size of the chunk to lead into memory for ATAC-seq read offset.
#' @export
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
                              seq_type="ATAC", # or "CHIP"
                              run_fastqc=TRUE,
                              # Directories
                              path_fastqc="FastQC/",
                              path_bam="BAM/",
                              path_peaks="Peaks/",
                              path_logs="Logs/",
                              # 1) Alignment params
                              out_name=NULL,
                              suffix_fastq=NULL,
                              type="SE",
                              index="/biodata/indices/species/Hsapiens/ucsc.hg19",
                              # 2) FiltOut params
                              remove=c("chrM", "chrUn", "_random", "_hap", "_gl", "EBVls"),
                              blacklist="~/data/consensusBlacklist.bed",
                              # 4) Peak calling params
                              peak_type=NULL,
                              shift=NULL,
                              cores=8,
                              chunk=1e7) {
  ## 0) Create directory tree -----------------------------
  dir.create(path_fastqc, F)
  dir.create(path_bam, F)
  dir.create(path_peaks, F)
  dir.create(path_logs, F)

  ## 1) Quality control with fastqc -----------------------
  if (run_fastqc) {
    fastqc(files=fastq_files,
           out_dir=path_fastqc)
  }

  ## 2) Alignment with Bowtie2 ----------------------------
    files <- mapply(alignmentBowtie2,
                    file=fastq_files,
                    out_name=out_name,
                    MoreArgs = list(suffix_fastq=suffix_fastq,
                                    type=type,
                                    index=index,
                                    out_dir=path_bam,
                                    path_logs=path_logs,
                                    cores=cores)
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
  if (seq_type=="ATAC") {
    files <- lapply(files,
                    offsetATACSE,
                    chunk=chunk)
  }

  ## 5) Peak calling with MACS2 ----------------------------

  if(seq_type=="ATAC" & is.null(peak_type)) { peak_type <- "narrow" }
  if(seq_type=="CHIP" & is.null(peak_type)) { peak_type <- "broad" }
  if(seq_type=="ATAC" & is.null(shift)) { shift <- TRUE }
  if(seq_type=="CHIP" & is.null(shift)) { shift <- FALSE }

  files <- lapply(files,
                  callPeak,
                  out_dir=path_peaks,
                  path_logs=path_logs,
                  type=peak_type,
                  shift=shift)
}
