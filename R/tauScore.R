#' Calculate tissue-specificity score
#' 
#' Given a list of counts for a gene/open chromatin region in different tissues, it 
#' returns a t score that goes from 0 (housekeeping) to 1 (tissue-specific).
#' @param counts Integer vector of counts of a specific gene/region in different tissue samples.
#' @return Tau score for tissue specificity of the specified gene/region in the specified tissues.
#' @export
#' @examples 
#' 

tauScore <- function(counts) {
  if(any(is.na(counts))) stop('NA\'s need to be 0.')
  if(any(counts<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-counts/max(counts))/(length(counts)-1)
} 