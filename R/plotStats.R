#' Plot alignment statistics
#'
#' @param stats Table obtained from getStats.
#' @param threshold Threshold for drawing the line (currently 80 percent). Indicates that all samples should have > 80 percent uniquely mapped reads.
#' @return Ggplot with the alignment statistics per sample.
#' @import ggplot2
#' @export
plotAlignment <- function(stats, threshold=80) {
  ## Make data.frame for plotting
  df_al <- stats[,c("aligned_reads", "multi_reads")]
  df_al$unaligned <- stats$total_reads - stats$aligned_reads
  df_al$unique <- stats$aligned_reads - stats$multi_reads
  df_al <- df_al[,-grep("aligned_reads", colnames(df_al))]
  df_al <- df_al / stats$total_reads * 100
  rownames(df_al) <- stats$sampleID
  df_al <- reshape2::melt(as.matrix(df_al))
  colnames(df_al) <- c("sampleID", "Type", "Percentage")
  df_al$Type <- factor(df_al$Type,
                       levels=rev(c("unique", "multi_reads", "unaligned")),
                       labels=rev(c("Unique", "Multi-mapping", "Not aligned")))

  ## Plot alignment results
  ggplot(df_al, aes(sampleID, Percentage)) +
    geom_bar(aes(fill=Type), stat="identity", color="black") +
    geom_hline(yintercept=threshold, lty=2, color="dark red") +
    coord_flip() +
    theme(legend.position = "top")
}
