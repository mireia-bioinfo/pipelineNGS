#############################################################
## Create the required phastcons46way file for pipelineNGS ##
#############################################################

out_dir <- "~/data/phastCons_46_placentalMammals"
wigToBigWig <- "~/bin/wigToBigWig"
chrom_info <- "~/refs/hg19_chrom_info.txt"

## List necessary files -------------------------------------------------------

url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/"
filenames = RCurl::getURL(url, verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
files <- XML::getHTMLLinks(filenames)

chr <- paste0("chr", c(1:22, "X", "Y"), ".")
files <- files[sapply(chr, grep, files, fixed=TRUE)]
files <- paste0(url, files)

## Download files -------------------------------------------------------------

dir.create(out_dir, recursive=TRUE)

for (i in files) {
  download.file(i, destfile=file.path(out_dir,
                                      basename(i)))
}

## Convert to bigWig ----------------------------------------------------------

wigs <- file.path(out_dir, basename(files))

cmds <- paste(wigToBigWig, wigs, chrom_info, gsub(".wigFix.gz", ".bw", wigs))

for (c in cmds) {
  message("Running", c)
  system2(c)
}

## Merge bigWig files ---------------------------------------------------------

cmd <- "bigWigMerge "
