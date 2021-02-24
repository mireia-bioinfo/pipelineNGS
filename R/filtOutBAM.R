#' Filter out reads from BAM file
#'
#' Function that filters out the following reads: 1) aligned no non-cannonical chromosomes, 2) not aligned, 3) aligned to ENCODE black-listed regions and 3) duplicated reads.
#' @param file Filenme and path for the BAM file to be filter out.
#' @param path_logs Path for the output logs.
#' @param type Either "SE" (single end) or "PE" (paired end).
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
                       type = "SE",
                       remove=c("chrM", "chrUn", "_random", "_hap", "_gl", "EBV"),
                       blacklist="~/data/consensusBlacklist.bed",
                       cores=6) {

  ## Remove blacklisted regions if necessary
  if (file.exists(blacklist)) {
    rm_blacklist <- paste("| samtools view - -b -L", blacklist, "-U", gsub(".raw", ".tmp", file, fixed=TRUE), "-o trash.bam")
  } else {
    rm_blacklist <- paste("-o", gsub(".raw", ".tmp", file, fixed=TRUE))
  }

  ## Select aligned, non-blacklisted and cannonical chr reads.
  message(paste0("[", format(Sys.time(), "%X"), "] ",
                 ">> Removing unaligned and blacklisted reads: ", file))
  cmd <- paste("samtools idxstats", file,
               "| cut -f 1 |", paste("grep -v", remove, collapse=" | "),
               "| xargs samtools view -b", file, "-@", cores-1, "-F 4", # select only mapped reads
               rm_blacklist
               )
  message(paste("\t", cmd))
  system(cmd)

  ## Fixmate if type is PE
  if(type=="PE") {
    message(paste0("[", format(Sys.time(), "%X"), "] ",
                   ">> Fixing mate: ", file))
    cmd_fixmate <- paste("samtools sort -n", gsub(".raw", ".tmp", file, fixed=TRUE),
                         "-@", cores-1, "-m 2G -o - | samtools fixmate -m - -",
                         "| samtools sort - -@", cores-1, "-m 2G -o", gsub(".raw", ".fixmate", file, fixed=TRUE))
    message(paste("\t", cmd_fixmate))
    system(cmd_fixmate)
    input_rmdup <- gsub(".raw", ".fixmate", file, fixed=TRUE)
  } else {
    input_rmdup <- gsub(".raw", ".tmp", file, fixed=TRUE)
  }

  ## Remove duplicates
  message(paste0("[", format(Sys.time(), "%X"), "] ",
                 ">> Removing duplicates: ", file))
  cmd <- paste("samtools markdup",
               input_rmdup, # input
               gsub(".raw", "", file, fixed=TRUE), # output
               "-r -s 2>",
               paste0(path_logs, getNameFromPath(file, suffix=".raw.bam"), ".rmdup.log"),
               "; samtools index", gsub(".raw", "", file, fixed=TRUE), "-@", cores-1,
               "; samtools idxstats", gsub(".raw", "", file, fixed=TRUE), ">", file.path(path_logs, paste0(getNameFromPath(file, suffix=".raw.bam"), ".idxstats.log"))
               )
  message(paste("\t", cmd))
  system(cmd)

  ## Remove temporary files
  rm_files <- c("trash.bam",
                gsub(".raw", ".fixmate", file, fixed=TRUE),
                gsub(".raw", ".tmp", file, fixed=TRUE))
  file.remove(rm_files[file.exists(rm_files)])

  return(invisible(gsub(".raw", "", file, fixed=TRUE)))
}

