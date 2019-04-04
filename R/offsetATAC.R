#' Offset correction for ATAC-seq data
#'
#' It performs the offset correction needed for ATAC-seq data: +4bp in forward strand; -5bp in reverse strand.
#' @param file Character string with the filename and path for the BAM file.
#' @param cores Number of cores to use for the processing.
#' @param gen_sizes Character string indicating the path where the file with chromosome name and sizes can be found. Default value: "~/data/hg19.len"
#' @return Returns the final BAM file, sorted and indexed, in path_bam.
#' @export
#' @examples
#' \dontrun{
#' offsetATAC(file)
#' }
offsetATAC <- function(file,
                       cores=6,
                       gen_sizes="~/data/hg19.len") {
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Offset correction Tn5: ", file))

  ## Offset correction
  cmd <- paste("bedtools bamtobed -i", file,
               "| awk 'BEGIN {OFS = \"\t\"} ; {if ($6 ==\"+\") print $1, $2 + 4, $3 +4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}'",
               "| awk -v OFS='\t' '$2<0 {$2=0} 1'",
               "| bedtools bedtobam -g", gen_sizes,
               "| samtools sort -", "-o", gsub(".bam", ".tmp.bam", file), "-m 1G -@", cores-1
               )
  system(cmd)

  ## Rehead BAM file
  cmd <- paste("samtools view -H", file,
               "| samtools reheader -", gsub(".bam", ".tmp.bam", file),
               ">", gsub(".bam", ".offset.bam", file),
               "; samtools index", gsub(".bam", ".offset.bam", file), "-@", cores-1)
  system(cmd)
  file.remove(gsub(".bam", ".tmp.bam", file))

  return(invisible(gsub(".bam", ".offset.bam", file)))
}
