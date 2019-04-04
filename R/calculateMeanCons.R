#' Calculates mean conservation score in a set of regions
#'
#' Using a set of regions as input, extends them \code{scope} from the center (upstream
#' and downstream) and bins them in \code{bins} bp. Then obtains the conservation score from
#' a phastCons bigWig file, and calculates the mean per bin in all regions.
#'
#' @param regions Set of regions as a \code{GRanges} object to use for the analysis. The first mcols column should either
#' contain a unique identifier for the region or be empty.
#' @param scope Number of bp to extend the regions from the peak center, both upstream
#' and downstream (Length of region = scope*2). Default=500bp.
#' @param bin Number of bp for binning the regions. Default=50bp.
#' @param phastConsBW Path and name of the bigWig file containing the phastCons conservation
#' scores (or any other bigWig for calculating mean score).
#' @import GenomicRanges
#' @return A data.frame object containing the bin position (from the peak center) and
#' the mean conservation score.
#' @export

calculateMeanCons <- function(regions,
                              scope=500,
                              bin=50,
                              phastConsBW="~/data/phastCons_46_placentalMammals/phastCons46way.placental.bw") {
  ## Generate a unique identifier for each query region (if not present)
  if(ncol(mcols(regions))==0 |
     length(unique(mcols(regions)[,1]))!=length(regions)) {
    regions$GeneID <- paste0("region_", 1:length(regions))
  } else {
    colnames(mcols(regions))[1] <- "GeneID"
  }

  ## Extend regions to scope*2
  regions.ext <- resize(regions, width=scope*2, fix="center")
  regions.ext$center <- start(ranges(regions.ext)) + width(regions)/2

  ## Bin regions
  regions.bin <- tile(regions.ext, width=bin)
  n <- unique(sapply(regions.bin, length))
  regions.bin <- unlist(regions.bin)
  regions.bin$GeneID <- rep(regions.ext$GeneID, each=n) ## add identifier
  regions.bin$center <- rep(regions.ext$center, each=n) ## add center
  regions.bin$pos <- start(regions.bin) - regions.bin$center ## add pos relative to center
  # unique(regions.bin$pos)

  ## Import conservation scores
  cons <- rtracklayer::import(phastConsBW,
                 which=regions.ext)

  ## Overlap data from regions and scores
  hits <- findOverlaps(regions.bin, cons)
  hits.df <- data.frame(hits)

  ## Calculate mean of all scores overlapping query regions
  cons_list <- split(score(cons[subjectHits(hits)]), queryHits(hits))
  cons_mean <- sapply(cons_list, mean, na.rm=T)

  ## Add mean for every region and position
  regions.cons <- regions.bin[unique(queryHits(hits)),]
  regions.cons$meanCons <- cons_mean

  ## Calculate mean for each position
  sel.cons.df <- data.frame(regions.cons)[,c(8,9)]

  mean_pos <- sapply(unique(sel.cons.df$pos),
                     function(x) mean(sel.cons.df$meanCons[sel.cons.df$pos==x]))

  ## Create output data.frame
  df <- data.frame("position"=unique(sel.cons.df$pos),
                   "meanCons"=mean_pos)
  return(df)
}

