#' Get name from path
#'
#' Obtains a file name from a path.
#' @param path Character string or vector with the path.
#' @param suffix Characer string with the suffix you want to remove from the filename. For example: ".bam"
#' @param prefix Character string with the prefix you want to remove from the filename. For example: "ATAC_"
#' @return A character string or vector with the filename.
#' @export
#' @examples
#' \dontrun{
#' getNameFromPath("test/data/cov_file.txt", suffix=".txt", preffix="cov_")
#' }
getNameFromPath <- function(path, suffix=NULL, prefix=NULL) {
  names <- sapply(path, .getNameFromPathSingle, suffix=suffix, prefix=prefix)
  names(names) <- NULL
  return(names)
}

.getNameFromPathSingle <- function(path,
                                   suffix=NULL,
                                   prefix=NULL) {
  split <- strsplit(path, "/")[[1]]
  name <- split[length(split)]

  if (!is.null(suffix)) name <- gsub(suffix, "", name, fixed=TRUE)
  if (!is.null(prefix)) name <- gsub(prefix, "", name, fixed=TRUE)

  return(name)
}
