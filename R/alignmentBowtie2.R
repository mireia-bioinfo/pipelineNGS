#' Alignment with Bowtie 2
#'
#' Aligns FastQ files to a reference genome using Bowtie2, converts the output to BAM, sorts it
#' according to genomic position and indexes the final BAM file.
#' @param file Character string (single-end) or character vector (paired-end) with the path and
#' filename of the sample to be analyzed.
#' @param out_name Character string indicating the name for the output file. If not provided, will use
#' suffix_fastq to remove the suffix and use the same name as fastq file.
#' @param suffix_fastq Character vector indicating the suffix to remove from fastq file in order to
#' generate output sam file. Only necessary if out_name is not provided. If this is a paired end sample
#' it should also include the name of read 1 (example: "_read1.fastq.gz").
#' @param type Character string with the type of approach used for the sequencing, either "SE" (default) or "PE".
#' @param out_dir Character string indicating the path where the final bam files should be stored.
#' @param cores Integer indicating the number of cores to use for the analysis.
#' @param index Character string with the path and name for the Bowtie2 index as specified in Bowtie2 manual.
#' @param path_logs Character string with the path for the logs to be stored.
#' @param run Logical indicating whether to run the alignment (for testing purposes). Default: TRUE
#' @param ... Further arguments to be passed to Bowtie2 call.
#' @return Writes "out_name.raw.bam" in out_dir. Also generates a log for the alignment (sample.alignment.log) in your path_logs.
#' @export
#' @examples
#' \dontrun{
#' alignmentBowtie2(file=paste0(path, "fastq/NL1_h3k27ac_sample.fastq.gz"),
#'                  suffix_fastq=".fastq.gz",
#'                  type="SE",
#'                  index=" /biodata/indices/species/Hsapiens/ucsc.hg19",
#'                  out_dir=paste0(path, "bam/"),
#'                  path_logs=paste0(path, "logs/"),
#'                  cores=6,
#'                  run=FALSE)
#' }
alignmentBowtie2 <- function(file,
                             out_name=NULL,
                             suffix_fastq=NULL,
                             type="SE",
                             out_dir="bam/",
                             cores=6,
                             index="/biodata/indices/species/Hsapiens/ucsc.hg19",
                             path_logs="logs/",
                             run=TRUE,
                             ...) {
  ## Create folders out_dir and path_logs
  .al <- out_dir; dir.create(.al, F) # Create BAM directory
  .log <- path_logs; dir.create(.log, F) # Create directory for logs

  if (!file.exists(file)) stop(paste("Input file", file, "does not exist"))

  ## Get names for aligned files
  if (is.null(out_name)) {
    name <- getNameFromPath(file[1], suffix=suffix_fastq)
  } else {
    name <- out_name
  }

  ## Check single-end or paired-end input
  if (type=="SE") { ## Align single end file
    if (length(file)>1) stop("Incorrect number of files, should be 1 for single-end alignment.")
    al_fastq <- paste("-U", file)

  } else if (type=="PE") { ## Align paired end file
    if (length(file)>2) stop("Incorrect number of files, should be 2 for paired-end alignment.")
    al_fastq <- paste("-1", file[1],
                      "-2", file[2], "-X2000")
  }

  ## Start alignment
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Alignment with Bowtie 2: ", name))
  cmd <- paste("bowtie2 -t",
               "-x", index,
               al_fastq, # Different parameter depending on SE or PE
               "-p", cores,
               ...,
               "2>", paste0(.log, name, ".alignment.log"),
               "| samtools view - -b -@", cores-1,
               "| samtools sort - -@", cores-1, "-m 2G",
               "-o", paste0(.al, name, ".raw.bam"),
               "; samtools index", paste0(.al, name, ".raw.bam"), "-@", cores-1
               )
  message(paste("\t", cmd))
  if(run) system(cmd)

  return(invisible(paste0(.al, name, ".raw.bam"))) # Returns path and filename
}
