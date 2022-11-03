#' Generate new fastqc files with the adapters trimmed using trim_galore method
#'
#' trim_galore() is a function that will automatically
#' @param experiment_type Character string detailing the type of experiment (i.e. if the sequencing is paired or single end)
#' @param adapter Character string with the type of adapter used when sequencing
#' @param read1 Character string with the full path of the first read of the sample
#' @param read2 Character string with the full path of the second read of the sample
#' @param output_directory Character string with the path of the output directory
#' @param return the FastQC report once the adapters are trimmed as well as a trimming report and the reads trimmed
#' @examples
#' \dontrun{
#' trim_galore(experiment_type ="paired", adapter ="illumina", read1 = r1.fastq.gz, read2 = r2.fastq.gz, output_directory = "~/username/folder"))
#' }

trim_galore <- function(experiment_type, adapter,
                        read1, read2, output_directory) {
  
  cp <- paste(paste0("trim_galore --",
                     experiment_type, " --",
                     adapter, " ",
                     paste0(read1), " ", paste0(read2),
                     " --fastqc", "
                     -o ", output_directory))
  system(cp)
  
}