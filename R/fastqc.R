#' Quality Control with FASTQC
#'
#' Runs FASTQC for quality control analysis of your FastQ reads.
#' @param files String or character vector with the path and filename of the fastq samples.
#' @param out_dir Character string with the name of the folder where you want to save FastQC results.
#' @param cores Integer indicating the number of cores to use for the analysis.
#' @return It creates a folder "out_dir" that contains the html and other files that are generated as a result of running FASTQC over your samples.
#' @export
#' @examples
#' \dontrun{
#' fastqc(samples="path/to/fastq/file.fastq.gz", out_dir="fastqc/", cores=6)
#' }
fastqc <- function(files,
                   out_dir="fastqc/",
                   cores=6) {
  message(">> Quality Control with FASTQC")
  dir.create(out_dir, F)

  for (i in files) {
    message(paste("----", i))
    cmd <- paste("fastqc",
                 i,
                 "-o", out_dir,
                 "-t", cores)
    system(cmd)
  }
}
