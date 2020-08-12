#############################################################
## Create the required phastcons46way file for pipelineNGS ##
#############################################################

out_dir <- "~/data/phastCons_46_placentalMammals"
wigToBigWig <- "wigToBigWig"
chrom_info <- "~/refs/hg19.chromSizes.txt"

## List necessary files -------------------------------------------------------

url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/"
filenames = RCurl::getURL(url, verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
files <- XML::getHTMLLinks(filenames)

chr <- paste0("chr", c(1:22, "X", "Y"), ".")
files <- files[sapply(chr, grep, files, fixed=TRUE)]
files <- paste0(url, files)

## Download files -------------------------------------------------------------

dir.create(out_dir, recursive=TRUE)

# for (i in files) {
#   download.file(i, destfile=file.path(out_dir,
#                                       basename(i)))
# }

## Convert to bigWig ----------------------------------------------------------

wigs <- file.path(out_dir, basename(files))
wigs <- wigs[order(wigs)]
bw <- gsub(".wigFix.gz", ".bw", wigs)

cmds <- paste(wigToBigWig, wigs, chrom_info, bw)

for (c in cmds) {
#  message("Running >>>", c)
#  system(c)
}

## Merge bigWigs ---------------------
cmd <- paste("bigWigMerge", paste0(bw, collapse=" "), file.path(out_dir, "placental_mammals.bedGraph"))

print(cmd)
system(cmd)

## Convert to bigWig ----------------
# Sort
cmd <- paste("sort -k1,1 -k2,2n --buffer-size=24G", file.path(out_dir, "placental_mammals.bedGraph"),
	     ">", file.path(out_dir, "placental_mammals.srtd.bedGraph"))
system(cmd)

# Convert
cmd <- paste("bedGraphToBigWig", file.path(out_dir, "placental_mammals.srtd.bedGraph"), chrom_info, 
	     file.path(out_dir, "placental_mammals.bw"))
system(cmd)
#cmd <- paste("zcat", paste0(wigs, collapse=" "), "|", wigToBigWig, "-clip stdin", 
#	     chrom_info, file.path(out_dir, "placental_mammals.bw"))
#print(cmd)
#system(cmd)
