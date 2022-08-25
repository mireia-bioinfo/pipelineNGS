#' Data parameters for different experiments
#'
#' A dataset containing the different parameters used for the pipelines
#' analyzing different types of experiments.
#'
#' @format A data frame with 4 rows and 6 columns
#' \describe{
#'   \item{Parameter}{Name of the parameter used as input (`seq_type`) to the function [process_epigenome()].}
#'   \item{Experiment}{Complete name of the experiment}
#'   \item{Extra_bowtie2}{Extra parameters to be passed to the bowtie2 aligment.}
#'   \item{Offset_correction}{Whether to perform (`TRUE`) or not (`FALSE`) transposase insertion offset correction.}
#'   \item{peak_type}{Whether to call `narrow` or `broad` peaks in MACS2.}
#'   \item{shift}{Whether to shift the coordinates of the peaks in MACS2.}
#' }
"paramdata"
