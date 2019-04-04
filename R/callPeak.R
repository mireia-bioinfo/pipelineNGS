#' Peak calling with MACS2
#'
#' Calls the peaks present in your sample bam files and returns bed files with the location of
#' the enriched regions.
#' @param file Character string with the filename for the BAM file.
#' @param out_dir Character string with the path for the peak calling result to be saved.
#' @param path_logs Character string with the path for the logs to be stored.
#' @param type Character string indicating if the resulting peaks should be "narrow" (default) or "broad".
#' @param shift Logical indicating wether the reads should be shifted -100bp and extended to 200bp or not (default).
#' @return Bed files in out_dir with the location of the enriched regions.
#' @export
#' @examples
#' \dontrun{
#' callPeak()
#' }
callPeak <- function(file,
                     out_dir,
                     path_logs,
                     type="narrow",
                     shift=FALSE) {

  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Starting Peak Calling with MACS2"))

  dir.create(out_dir, F)

  ## Create parameters
  if (type=="narrow") type_spec <- "-q 0.05 --nomodel"
  else if (type=="broad") type_spec <- "--broad --broad-cutoff 0.1 --nomodel"
  if (shift==TRUE) type_spec <- paste(type_spec, "--shift -100 --extsize 200")

  ## Temporary folder
  tmp <- paste0(out_dir, "tmp/"); dir.create(tmp, F)

  ## Get name from bam path
  name <- getNameFromPath(file, suffix=".bam")
  name <- gsub(".offset", "", name) # In case is ATAC-seq data

  ## Run peak calling
  message(paste0("[", format(Sys.time(), "%X"), "] ", "---- ", file))
  cp <- paste("macs2 callpeak",
              "-f BAM -t", file,
              "-g hs --outdir", out_dir, "-n", name,
              "--tempdir", tmp,
              type_spec,
              "2>", paste0(path_logs, name, ".macs2.log")
              )
  print(cp)
  system(cp)

  # system(paste0("rm -r", tmp))
}
