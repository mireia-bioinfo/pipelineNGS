#' Offset correction for ATAC-seq data
#'
#' It performs the offset correction needed for ATAC-seq data: +4bp in forward strand; -5bp in reverse strand.
#' @param file Character string with the filename and path for the BAM file.
#' @inheritParams process_epigenome
#' @return Returns the final BAM file, sorted and indexed, in path_bam.
#' @export
#' @examples
#' \dontrun{
#' offsetATAC(file)
#' }
offsetATAC <- function(file,
                       cores=6,
                       gen_sizes="/vault/refs/hg38.chromSizes.txt") {
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

#' Offset correction for SE ATAC-seq data
#'
#' It performs offset corrrection on single end ATAC-seq files.
#' @param file Path for the BAM file.
#' @param positive Positive offset: 4.
#' @param negative Negative offset: 5.
#' @inheritParams process_epigenome
#' @importClassesFrom Rsamtools BamFile
#' @importMethodsFrom S4Vectors metadata
#' @importFrom rtracklayer export
#' @export
offsetATACSE <- function(file,
                         chunk=1e6,
                         positive = 4L,
                         negative = 5L) {
  message(paste0("[", format(Sys.time(), "%X"), "] ", ">> Offset correction Tn5: ", file))

  outbam <- gsub(".bam", ".offset.bam", file)

  ## Obtain possible tags to read
  possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2",
                                  "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                  "TC", "UQ"),
                      "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                    "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                    "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                    "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                    "U2"))
  bamTop100 <- Rsamtools::scanBam(Rsamtools::BamFile(file, yieldSize = 100),
                                  param = Rsamtools::ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]

  ## Load BAM file
  gal <- ATACseqQC::readBamFile(file, tag=tags, asMates=FALSE, bigFile=TRUE)
  meta <- metadata(gal)

  ## Shift alignments using chunk sizes
  index <- ifelse(length(meta$index)>0, meta$index, meta$file)
  bamfile <- Rsamtools::BamFile(meta$file, index=index,
                                yieldSize=chunk, asMates = meta$asMates)
  outfile <- NULL
  open(bamfile)
  on.exit(close(bamfile))

  while (length(chunk0 <- GenomicAlignments::readGAlignments(bamfile, param=meta$param))) {
    gal1 <- ATACseqQC::shiftGAlignments(chunk0,
                                        positive = positive,
                                        negative = negative)
    outfile <- c(tempfile(fileext = ".bam"), outfile)
    ATACseqQC:::exportBamFile(gal1, outfile[1])
    rm(gal1)
  }
  close(bamfile)
  on.exit()

  ## Create merged files from outfiles
  if(length(outfile)>1){
    mergedfile <- Rsamtools::mergeBam(outfile,
                                      destination=tempfile(fileext = ".bam"),
                                      indexDestination=TRUE,
                                      header=meta$file)
    unlink(outfile)
    unlink(paste0(outfile, ".bai"))
  }else{
    if(length(outfile)==1){
      mergedfile <- outfile
    }else{
      stop("Can not get any proper mapped reads from  your inputs.")
    }
  }

  ## Rename file to match outbam name
  if(!missing(outbam)){
    file.copy(from=mergedfile, to=outbam)
    file.copy(from=paste0(mergedfile, ".bai"),
                to=paste0(outbam, ".bai"))
    file.remove(mergedfile)
    file.remove(paste0(mergedfile, ".bai"))
  }else{
    gal1 <- GenomicAlignments::readGAlignments(mergedfile, param = meta$param)
    mcols(gal1)$MD <- NULL
    names(gal1) <- mcols(gal1)$qname
    gal1 <- gal1[order(names(gal1))]
    unlink(mergedfile)
    unlink(paste0(mergedfile, ".bai"))
  }
  return(outbam)
}
