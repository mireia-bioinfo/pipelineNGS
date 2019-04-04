#' Annotate regions to nearest gene
#'
#' Annotates genomic regions to nearest feature in a GRanges object.
#' @param regions GRanges including the query regions
#' @param genes GRanges including data for genes to be annotated to. Should contain strand info.
#' @param promLimit Limit for annotating a region as a promoter. Default: 2kb.
#' @return Returns an annotated data.frame including info of the regions and genes, along with distance to the gene and annotation.
#' @export
annotateToNearest <- function(regions,
                              genes,
                              promLimit=2000) {
  genes.p <- promoters(genes, upstream=1, downstream=0, use.names=TRUE)

  hits <- distanceToNearest(regions, genes.p, ignore.strand=TRUE)

  reg.anno <- cbind(data.frame(regions)[queryHits(hits), -c(4:5)],
                    data.frame(genes)[subjectHits(hits), -c(4:5)],
                    mcols(hits)$distance)
  colnames(reg.anno)[ncol(reg.anno)] <- "distanceToTSS"
  reg.anno$annotation <- "Distal RE"
  reg.anno$annotation[abs(reg.anno$distanceToTSS)<promLimit] <- "Promoter"
  return(reg.anno)
}
