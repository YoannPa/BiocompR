
#' Filters irrelevant warnings matching a regular expression.
#'
#' @param pattern A \code{character} to be processed as a regular expression.
#' @return A \code{character} list of warning messages without the warnings
#'         matching the pattern.
#' @author Yoann Pageaud.
#' @export

warn.handle <- function(pattern, ...){
  withCallingHandlers(..., warning = function(warning){
    condition <- conditionMessage(warning)
    if(grepl(pattern, condition)){ invokeRestart("muffleWarning") }
  })
}
