#' Alignment with Bowtie 2
#'
#' Aligns FastQ files to a reference genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
#' converts the output to BAM, sorts it according to genomic position and
#' indexes the final BAM file.
#' @param file Character string (single-end) or character vector of
#' length 2 (paired-end) with the file names of the samples to be analysed.
#' @param run Logical indicating whether to run the alignment (for testing purposes). Default: TRUE
#' @inheritParams process_epigenome
#' @return Writes `out_name.raw.bam` in `path_bam`. Also generates a log for the alignment (`sample.alignment.log`) in your `path_logs`.
#' @export
#' @examples
#' \dontrun{
#' alignmentBowtie2(file=paste0(path, "fastq/NL1_h3k27ac_sample.fastq.gz"),
#'                  suffix_fastq=".fastq.gz",
#'                  type="SE",
#'                  index="/vault/refs/indexes/hg38",
#'                  path_bam=paste0(path, "bam/"),
#'                  path_logs=paste0(path, "logs/"),
#'                  cores=6,
#'                  run=FALSE)
#' }
alignmentBowtie2 <- function(file,
                             out_name=NULL,
                             type="SE",
                             suffix_fastq=NULL,
                             cores=6,
                             index="/vault/refs/indexes/hg38",
                             path_bam="bam/",
                             path_logs="logs/",
                             extra_bowtie2="",
                             run=TRUE) {
  ## Create folders path_bam and path_logs
  .al <- path_bam; dir.create(.al, F) # Create BAM directory
  .log <- path_logs; dir.create(.log, F) # Create directory for logs

  if (!all(file.exists(unlist(file)))) stop(paste("Input file", file, "does not exist"))

  ## Get names for aligned files
  if (is.null(out_name)) {
    name <- getNameFromPath(file[1], suffix=suffix_fastq)
  } else {
    name <- out_name
  }

  ## Check single-end or paired-end input
  if (type=="SE") { ## Align single end file
    al_fastq <- paste("-U", paste0(file, collapse=" "))
  } else if (type=="PE") { ## Align paired end file
    if (length(file)>2) stop("Incorrect number of files, should be 2 for paired-end alignment.")
    al_fastq <- paste("-1", paste0(unlist(file[1]), collapse=" "),
                      "-2", paste0(unlist(file[2]), collapse=" "))
  }

  ## Start alignment
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Alignment with Bowtie 2: ", name))
  cmd <- paste("bowtie2 -t --local",
               "-x", index,
               al_fastq, # Different parameter depending on SE or PE
               "-p", cores,
               extra_bowtie2,
               "2>", paste0(.log, name, ".alignment.log"),
               "| samtools view - -b -@", cores-1,
               "| samtools sort - -@", cores-1, "-m 2G", "-o", paste0(.al, name, ".raw.bam"),
               "; samtools index", paste0(.al, name, ".raw.bam"), "-@", cores-1,
               "; samtools idxstats", paste0(.al, name, ".raw.bam"), ">", file.path(.log, paste0(name, ".raw.idxstats.log"))
               )
  message(paste("\t", cmd))
  if(run) system(cmd)

  return(invisible(paste0(.al, name, ".raw.bam"))) # Returns path and filename
}
