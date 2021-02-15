#' Filter out reads from BAM file
#'
#' Function that filters out the following reads: 1) aligned no non-cannonical chromosomes, 2) not aligned, 3) aligned to ENCODE black-listed regions and 3) duplicated reads.
#' @param file Filenme and path for the BAM file to be filter out.
#' @param path_logs Path for the output logs.
#' @param remove List of regular expressions to filter out from chromosome names.
#' @param blacklist Path for ENCODE blacklist in bed format.
#' @param cores Number of cores to use for the analysis.
#' @return File without the ".raw" preffix containing the final reads, along with its index (.bai).
#' Also returns invisibly the name of the output file.
#' @export
#' @examples
#' \dontrun{
#' file <- filtOutBAM(file=file,
#'                    path_logs=path_logs,
#'                    cores=8)
#' }
filtOutBAM <- function(file,
                       path_logs,
                       remove=c("chrM", "chrUn", "_random", "_hap", "_gl"),
                       blacklist="~/data/consensusBlacklist.bed",
                       cores=6) {

  ## Select aligned, non-blacklisted and cannonical chr reads.
  message(paste0("[", format(Sys.time(), "%X"), "] ",
                 ">> Removing unaligned and blacklisted reads: ", file))
  cmd <- paste("samtools idxstats", file,
               "| cut -f 1 |", paste("grep -v", remove, collapse=" | "),
               "| xargs samtools view -b", file, "-@", cores-1, "-F 4", # select only mapped reads
               "| samtools view - -b -L", blacklist, "-U", gsub(".raw", ".tmp", file, fixed=TRUE),
               "-o trash.bam"
               )
  system(cmd)

  ## Remove duplicates
  message(paste0("[", format(Sys.time(), "%X"), "] ",
                 ">> Removing duplicates: ", file))
  cmd <- paste("samtools markdup",
               gsub(".raw", ".tmp", file, fixed=TRUE), # input
               gsub(".raw", "", file, fixed=TRUE), # output
               "-r -s 2>",
               paste0(path_logs, getNameFromPath(file, suffix=".raw.bam"), ".rmdup.log"),
               "; samtools index", gsub(".raw", "", file, fixed=TRUE), "-@", cores-1,
               "; samtools idxstats", gsub(".raw", "", file, fixed=TRUE), ">", file.path(path_logs, paste0(getNameFromPath(file, suffix=".raw.bam"), ".idxstats.log"))
               )
  system(cmd)

  ## Remove temporary files
  file.remove("trash.bam",
              gsub(".raw", ".tmp", file, fixed=TRUE))

  return(invisible(gsub(".raw", "", file, fixed=TRUE)))
}

