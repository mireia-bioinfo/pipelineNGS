#' Assigning counts to transcripts with HTseqCount
#' 
#' Using an input BAM file, it assigns the counts to each transcript annotated in the transcriptome of interest.
#' @param samples Single element or character vector with the name of the samples (whitout any specific file format suffix).
#' @param path_bam Character string indicating the path where the final bam files should be stored.
#' @param transcript Transcriptome of reference to use for assigning the counts.
#' @return A file ".htseq.stats" with the general statistics and a file ".htseq.counts" with the number of counts assigned to each transcript.
#' @export
#' @examples 
#' \dontrun{
#' HTseqCount(samples=c("sample1", "sample2"), path_bam="~/test/samples/bam/",
#'                 transcript="~/data/gencode.v18.annotation.gtf")
#' }


HTseqCount <- function(samples, bam_path, transcript) {
  message(">> Assigning counts to known transcripts (HTseq-count)")
  for (i in samples) {
    message(paste("----", i))
    cmd <- paste("htseq-count",
                 "-s no",
                 "-r pos",
                 "-f bam", paste0(.al, i, ".bam"),
                 transcript,
                 ">", paste0(.al, i, "-htseq.txt")
    )
    system(cmd)
    
    cmd <- paste("tail -n 5", paste0(.al, i, "-htseq.txt"), ">", paste0(.al, i, ".htseq.stats"))
    system(cmd)
    cmd <- paste("head -n -5", paste0(.al, i, "-htseq.txt"), ">", paste0(.al, i, ".htseq.counts"))
    system(cmd)
  }
}
