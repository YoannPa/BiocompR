
#' Keeps, removes or imputes missing values in a matrix or a data.frame based on
#' sample groups.
#'
#' @param data   A \code{matrix} or \code{data.frame} with column names.
#' @param method A \code{character} which defines the method to use for managing
#'               NAs:
#'               \itemize{
#'                \item{'keep' does nothing.}
#'                \item{'remove' removes all rows containing at least 1 NA.}
#'                \item{'impute' removes rows containing only NAs and impute
#'                      missing values in remaining rows. The imputation uses
#'                      the groups defined in the parameter 'groups': e.g. in
#'                      row 1 if the value of sample_1 is NA, the missing value
#'                      will be imputed using the median of values from other
#'                      samples of the same group than sample_1.}
#'               }
#' @param groups A vector of length ncol(data) specifying the groups to which
#'               samples belong. Order of groups in the vector has to match the
#'               order of column names in data to properly associate samples and
#'               groups.
#' @value A \code{matrix}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references
#' @keywords internal

manage.na<-function(data, method = "remove", groups){
  if(method != "keep"){ #If NA should not be kept
    if(any(is.na(data))) { #If any NA
      if(method == "remove"){ #Remove all rows containing any NA
        data<-data[complete.cases(data),]

      } else if(method == "impute"){ #Impute NAs with the median value by group
        #Remove rows containing only NAs
        data<- data[!apply(X = is.na(data),MARGIN = 1,FUN = all), ,drop = FALSE]

        #Get groups of samples from sample conditions
        grp_tbl<-data.frame(samples = colnames(data), groups = groups)
        sample_grps<-unique(groups)

        #Get median by cell type/line
        apply(X = head(DMRs), MARGIN = 1, FUN = function(row){
          if(anyNA(data[row,])){
            row_vec<-sapply(sample_grps, function(grp){
              #List sample names matching the group
              samples.in.grp<-grp_tbl[groups == grp,]$samples
              #List Methylation values of the group on the row
              meth_vals<-data[row,samples.in.grp]
              #If more than 1 sample in the group get median
              if (length(samples.in.grp) > 1) { median(meth_vals,na.rm = T)}
              else { meth_vals } #Else get value of sample
            })
            #If no NA in group medians keep row for manipulation; Else remove it
            if (anyNA(row_vec) == F) {
              #If no NA in full row do nothing
              if(anyNA(data[row,])){
                #Replace NA value by median of its group if some
                data[row,]<-sapply(X = names(data[row,]), FUN = function(smpl){
                  #Get methylation value
                  valmeth<-data[row, smpl]
                  if (is.na(valmeth)){
                    row_vec[grp_tbl[grp_tbl$samples == smpl,]$groups]
                  } else { valmeth } #Else keep value as is
                })
              }
            } else { stop("Some group medians have for value NA.")
              #TODO: Handle this case if the issue arises.
            }
          }
        })
      }
    }
  }
  return(data)
}
