#' Get name from path
#'
#' Obtains a file name from a path.
#' @param path Character string or vector with the path.
#' @param prefix Character string with the prefix you want to remove from the filename. For example: "ATAC_"
#' @param suffix Characer string with the suffix you want to remove from the filename. For example: ".bam"
#' @return A character string or vector with the filename.
#' @export
#' @examples
#' \dontrun{
#' getNameFromPath("test/data/cov_file.txt", suffix=".txt", preffix="cov_")
#' }
getNameFromPath <- function(path, prefix=NULL, suffix=NULL) {
  names <- basename(path)
  if (!is.null(suffix)) names <- gsub(suffix, "", names, fixed=TRUE)
  if (!is.null(prefix)) names <- gsub(prefix, "", names, fixed=TRUE)

  return(names)
}
