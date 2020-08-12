#' Conversion from BAM to visualization formats
#'
#' Converts post-processed BAM files to visualization formats such as BedGraph (compressed) and BigWig.
#' @param samples Character string or character vector with the name of the samples (whitout any specific file format suffix).
#' @param summaryStats Character string with the location of the summary statistics .rda file.
#' @param scales Number or vector with the scales for the samples. Default is using those in summaryStats.
#' @param path_bam Character string indicating the path where the final bam files should be stored.
#' @param path_vis Character string indicating the path where you want to save visualization files (it will be created if it does not exist).
#' @param cores Number of cores to use.
#' @param perc_mem Percentage of memory to use.
#' @param gen_sizes Character string indicating the path where the file with chromosome name and sizes can be found. Default value: "~/data/hg19.len"
#' @param path_bedGraphToBigWig Character string indicating the path where the program converter from bedGraph to BigWig can be found. Default value: "~/tools/bedGraphToBigWig".
#' @param path_bedtools_genomecov Character string indicating the path where the bedtools genomecov binary can be found.
#' @param unc Logical indicating whether to compress or not the resulting bigWig file.
#' @return Two visualization files: bedgraph (.bdg) and bigwig (.bw) in your path_vis.
#' @export
#' @examples
#' \dontrun{
#' bamToVis("thyroid_sample", summaryStats="testData/logs/summary.rda", path_bam="testData/bam/",
#'          path_vis="testData/vis/")
#' }
bamToVis <- function(samples, summaryStats="", path_bam, path_vis, scales="",
                     gen_sizes="~/refs/hg19.chromSizes.txt",
                     cores = 8,
                     perc_mem = "50%",
                     path_bedGraphToBigWig="bedGraphToBigWig",
                     path_betools_genomecov="bedtools genomecov",
                     unc = FALSE) {
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Converting to BedGraph and Bigwig"))
  dir.create(path_vis, F)
  if (summaryStats!="") load(summaryStats)

  c = 1
  for (i in samples) {
    message(paste0("[", format(Sys.time(), "%X"), "] ", "---- ", i))

    if (scales=="") {
      scale <- sum.df[grep(i, sum.df$samples, fixed=TRUE),16]
    } else {
      scale <- scales[c] }

    c = c + 1
    message(paste0("[", format(Sys.time(), "%X"), "] ", "------ Converting to BedGraph"))
    bdg <- paste(path_betools_genomecov, "-bg -split",
                 "-scale", scale,
                 "-ibam", file.path(path_bam, paste0(i, ".bam")),
                 "| LC_COLLATE=C sort -k 1,1 -k2,2n --parallel", cores, "-S", perc_mem,
                 ">", file.path(path_vis, paste0(i, ".bedgraph"))
                 )
    system(bdg)
    message(paste0("[", format(Sys.time(), "%X"), "] ", "------ Converting to BigWig"))
    bigWig <- paste(path_bedGraphToBigWig,
                    file.path(path_vis, paste0(i, ".bedgraph")),
                    gen_sizes,
                    file.path(path_vis, paste0(i, ".bw"))
                    )
    if (unc) paste(bigWig, "-unc")

    system(bigWig)

    message(paste0("[", format(Sys.time(), "%X"), "] ", "------ Compressing BedGraph"))
        comp.bdg <- paste("bgzip -f",
                      file.path(path_vis, paste0(i, ".bedgraph")))
    system(comp.bdg)
  }
}
