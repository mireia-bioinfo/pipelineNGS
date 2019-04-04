#' Alignment with TopHat
#' 
#' Aligns FastQ files to a reference genome using TopHat (RNA-seq data).
#' @param samples Single element or character vector with the name of the samples (whitout any specific file format suffix).
#' @param type Character string with the type of approach used for the sequencing, either "SE" (default) or "PE".
#' @param path_fastq Character string indicating the path which contains the fastq files.
#' @param suffix_fastq Character string with the suffix of the input fastq files (ex. "fq.gz", "fastq.gz").
#' @param path_bam Character string indicating the path where the final bam files should be stored.
#' @param cores Integer indicating the number of cores to use for the analysis.
#' @param trans_index Indexed gene annotation (in our case gencode v18) by TopHat.
#' @param index Path with the path and name for the Bowtie2 as specified in Bowtie2 manual.
#' @return Returns a BAM file with the algimnent.
#' @export
#' @examples
#' \dontrun{
#' alignmentTophat(samples=c("sample1", "sample2"), type="SE", index="/biodata/index/ucsc.hg19",
#'                 path_fastq="~/test/samples/raw/", suffix_fastq="fq.gz", path_bam="~/test/samples/bam/",
#'                 cores=6, trans_index="~/data/gencode.v18.annotation.tophat2index")
#'}

alignmentTophat <- function(samples, type="SE", path_fastq, suffix_fastq, path_bam, cores,
                            trans_index, index) {
  .al <- path_bam; dir.create(.al, F)
  .tmp <- paste0(path_bam, "/tmp/"); dir.create(.tmp, F)
  
  message(">> Alignment with TopHat")
  if (type == "SE") {
    for (i in samples) {
      message(paste("----", i))
      cmd <- paste("tophat2",
                   "-r 100",
                   "-o", paste0(.al, i),
                   "-p", cores,
                   "--no-coverage-search",
                   paste0("--transcriptome-index=", trans_index),
                   index,
                   paste0(path_fastq, i, ".", suffix_fastq),
                   "2>", paste0(.log, i, ".alignment.log")
                   )
      message(paste("\t", cmd))
      system(cmd)
      
      # Rename files and move them to root folder
      cmd <- paste("cp", paste0(.al, i, "/accepted_hits.bam"), paste0(.al, i, ".bam"))
      system(cmd)
    }
  } else if (type=="PE") {
    for (i in samples) {
      message(paste("----f", i))
      cmd <- paste("tophat2",
                   "-r 100",
                   "-o", paste0(.al, i),
                   "--no-coverage-search",
                   "-p", cores,
                   paste0("--transcriptome-index=", params$parameters$trans.index),
                   params$parameters$index,
                   paste0(path_fastq, i, "_1.", suffix_fastq),
                   paste0(path_fastq, i, "_2.", suffix_fastq),
                   "2>", paste0(.log, i, ".alignment.log")
                   )
      message(paste("\t", cmd))
      system(cmd)
      
      # Rename files and move them to root folder
      cmd <- paste("cp", paste0(.al, i, "/accepted_hits.bam"), paste0(.al, i, ".bam"))
      system(cmd)
    }
  } else {
    stop("Invalidad type of experiment - should be either SE or PE")
  }
}