# pipelineNGS 0.1.0

-   Added a `NEWS.md` file to track changes to the package.
-   Added `pkgdown` site
-   Current experiments: ATAC, CHIP, CT (CUT&TAG) and CR (CUT&RUN).
-   Default alignment is `--local` to avoid trimming reads.
-   Function `process_CUTTAG()` is deprecated. Now all analysis can be performed with calls to `process_epigenome()`.
