#' Merge BAM files
#'
#' Merges several BAM files in one.
#' @param files Name for the bam files you want to merge.
#' @param path_bam Character string indicating the path where the input bam files are.
#' @param path_out Character string indicating the path where the final merged bam files should be stored.
#' @param out_name Character string with the name for the merged BAM file.
#' @param cores Integer specifying the number of cores to use.
#' @return Bed files in path_pc with the location of the enriched regions.
#' @export
#' @examples
#' \dontrun{callPeak()}

mergeBAMs <- function(files, out_name, path_bam, path_out, cores=6) {
  merge <- paste("samtools merge -f -@", cores-1,
                 paste0(path_out, out_name, ".bam"),
                 paste0(path_bam, files, collapse=" "))
  system(merge)
}
