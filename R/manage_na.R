
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
#' @param ncores An \code{integer} to specify the number of cores/threads to be
#'               used to parallel-process the matrix (Default: ncores = 1).
#' @return A \code{matrix}.
#' @author Yoann Pageaud.
#' @examples
#' # Creating a dataframe of 4 columns with 50 values in each:
#' df <- data.frame("col1" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col2" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col3" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col4" = rnorm(n = 50, mean = 25, sd = 4))
#' # col1 and col4 are missing data
#' df$col1[1:10] <- NA
#' df$col4[40:50] <- NA
#' # Create 2 groups: grp1 & grp2.
#' # col1 and col2 are in grp1, col3 and col4 are in grp2.
#' grps <- c("grp1","grp1","grp2","grp2")
#' # Remove all rows containing at least 1 missing value
#' manage.na(data = test, method = "remove", groups = grps)
#' # Impute by group the missing data (and remove rows missing all values)
#' manage.na(data = test, method = "impute", groups = grps)
#' @export
#' @references
#' @keywords internal

manage.na <- function(data, method = "remove", groups, ncores = 1){
  if(is.data.frame(data)){ data <- as.matrix(data) }
  if(method != "keep"){ #If NA should not be kept
    if(any(is.na(data))) { #If any NA
      if(method == "remove"){ #Remove all rows containing any NA
        data <- data[complete.cases(data),]
      } else if(method == "impute"){ #Impute NAs with the median value by group
        #Parallel-remove rows containing only NAs
        data <- data[
          !unlist(parallel::mclapply(
            X = seq(nrow(data)), mc.cores = ncores, FUN = function(r){
                             all(is.na(data[r, ])) })), , drop = FALSE]
        #Get groups of samples from sample conditions
        grp_tbl <- data.frame(samples = colnames(data), groups = groups)
        sample_grps <- unique(groups)
        #Get median by groups
        #TODO: Make parallel apply when it will be possible.
        invisible(lapply(X = seq(nrow(data)), FUN = function(r){
          invisible(lapply(X = sample_grps, FUN = function(grp){
            #List sample names matching the group
            samples.in.grp <- grp_tbl[groups == grp, ]$samples
            #List Methylation values of the group on the row
            if(all(is.na(data[r, samples.in.grp]))){
              #TODO: Handle this case if the issue arises.
              stop("The group as NA for all its values.")
            } else if(anyNA(data[r, samples.in.grp])){
              data[r, samples.in.grp][is.na(data[r, samples.in.grp])] <<-
                median(x = data[r, samples.in.grp], na.rm = TRUE)
            }
          }))
        }))
      }
    }
  }
  return(data)
}
