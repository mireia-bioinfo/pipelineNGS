#' Conversion from BAM to visualization formats
#'
#' Converts post-processed BAM files to visualization formats such as BedGraph (compressed) and BigWig.
#' @param bamfile Path and name of the bam file to convert.
#' @param scaling_factor Factor to use for scaling.
#' @param path_bw Character string indicating the path where you want to save visualization files (it will be created if it does not exist).
#' @param chr_sizes Character string indicating the path where the file with chromosome name and sizes can be found.
#' @param bedtools_genomecov Character string indicating the path where the program converter from bedGraph to BigWig can be found. Default assumes it is in your PATH.
#' @param bdg_to_bw Character string indicating the path where the bedtools genomecov binary can be found. Default assumes it is in your PATH.
#' @param uncompressed_bw Logical indicating whether to compress or not the resulting bigWig file.
#' @param cores Number of cores to use.
#' @param perc_mem Percentage of memory to use.
#' @return Two visualization files: bedgraph (.bdg) and bigwig (.bw) in your path_bw.
#' @export
bamToBigWig <- function(bamfile,
                        scaling_factor=1,
                        path_bw,
                        chr_sizes,
                        bedtools_genomecov = "bedtools genomecov",
                        bdg_to_bw = "bedGraphToBigWig",
                        uncompressed_bw = FALSE,
                        cores = 5,
                        perc_mem = "50%") {
  ## Create output directory
  dir.create(path_bw, F)

  ## Obtain filename
  filename <- gsub(".bam", "", basename(bamfile))

  ## Make bedgraph with bedtools
  bdg <- paste(bedtools_genomecov, "-bg -split",
               "-scale", scaling_factor,
               "-ibam", bamfile,
               "| LC_COLLATE=C sort -k 1,1 -k2,2n --parallel", cores, "-S", perc_mem,
               ">", file.path(path_bw, paste0(filename, ".bedgraph"))
  )
  system(bdg)

  ## Convert to bigWig with UCSC tools
  message(paste0("[", format(Sys.time(), "%X"), "] ", "------ Converting to BigWig"))
  bigWig <- paste(bdg_to_bw,
                  file.path(path_bw, paste0(filename, ".bedgraph")),
                  chr_sizes,
                  file.path(path_bw, paste0(filename, ".bw"))
  )
  if (uncompressed_bw) paste(bigWig, "-unc")

  system(bigWig)

  ## Compress bedgraph
  message(paste0("[", format(Sys.time(), "%X"), "] ", "------ Compressing BedGraph"))
  comp.bdg <- paste("bgzip -f",
                    file.path(path_bw, paste0(filename, ".bedgraph")))
  system(comp.bdg)
}
