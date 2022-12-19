# pipelineNGS 0.1.1

-   Fixed error `In file.rename(from = paste0(mergedfile, ".bai"), to = paste0(outbam,  :  cannot rename file '/tmp/RtmpkINXQ0/file2548651671fb07.bam.bai' to 'BAM/EC04.offset.bam.bai', reason 'Invalid cross-device link'`. Now `offsetATACSE()` uses `file.copy()` instead of `file.rename()`.

# pipelineNGS 0.1.0

-   Fixed PE stats to total reads, not pairs.

# pipelineNGS 0.1.0

-   Added a `NEWS.md` file to track changes to the package.
-   Added `pkgdown` site
-   Current experiments: ATAC, CHIP, CT (CUT&TAG) and CR (CUT&RUN).
-   Default alignment is `--local` to avoid trimming reads.
-   Function `process_CUTTAG()` is deprecated. Now all analysis can be performed with calls to `process_epigenome()`.
