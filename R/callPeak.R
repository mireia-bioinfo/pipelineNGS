#' Peak calling with MACS2
#'
#' Calls the peaks present in your sample bam files and returns bed files with the location of
#' the enriched regions.
#' @param file Character string with the filename for the BAM file.
#' @inheritParams process_epigenome
#' @return Bed files in `path_peaks` with the location of the enriched regions.
#' @export
#' @examples
#' \dontrun{
#' callPeak()
#' }
callPeak <- function(file,
                     path_peaks,
                     path_logs,
                     type_peak="narrow",
                     shift=FALSE,
                     type = "SE" # or PE
                     ) {

  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Starting Peak Calling with MACS2"))

  dir.create(path_peaks, F)

  ## Create parameters
  if (type_peak=="narrow") type_spec <- "-q 0.05 --nomodel"
  else if (type_peak=="broad") type_spec <- "--broad --broad-cutoff 0.1 --nomodel"

  if (shift==TRUE) type_spec <- paste(type_spec, "--shift -100 --extsize 200")

  if (type == "SE") format_spec <- "BAM"
  else if (type == "PE") format_spec <- "BAMPE"

  ## Temporary folder
  tmp <- paste0(path_peaks, "tmp/"); dir.create(tmp, F)

  ## Get name from bam path
  name <- getNameFromPath(file, suffix=".bam")
  name <- gsub(".offset", "", name) # In case is ATAC-seq data

  ## Run peak calling
  message(paste0("[", format(Sys.time(), "%X"), "] ", "---- ", file))
  cp <- paste("macs2 callpeak",
              "-f", format_spec,
              "-t", file,
              "-g hs --outdir", path_peaks, "-n", name,
              "--tempdir", tmp,
              type_spec,
              "2>", paste0(path_logs, name, ".macs2.log")
              )
  print(cp)
  system(cp)

  system(paste0("rm -r", tmp))
}
