
##IMPORTS

##FUNCTIONS

# reduce_IDs ###################################################################

#' Simplify samples IDs.
#' 
#' @param samples.names  A \code{character} vector of sample IDs.
#' @return a \code{character} vector of simplified sample IDs.
#' @author Yoann Pageaud.

reduce_IDs<-function(samples.names){
  #hsc --> h & mut--> m & mu-chip-wt --> wt
  smpl.names<-strsplit(unlist(strsplit(unlist(lapply(strsplit(sample.names,
                                                              split = "sc"),
                                                     paste0,collapse="")),
                                       split = "u-chip-mut")),
                       split = "mu-chip-")
  #mpp --> m & wt_1 --> w
  smpl.names<-strsplit(unlist(lapply(strsplit(unlist(lapply(smpl.names,paste0,
                                                            collapse = "")),
                                              split = "pp"),paste0,
                                     collapse = "")),"t_1")
  smpl.names<-toupper(unlist(lapply(strsplit(unlist(smpl.names),split = "_"),
                                    paste0,collapse="")))
  
  return(smpl.names)
}