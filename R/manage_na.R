
#' Imputes missing at random (MAR) data in a matrix row.
#'
#' @param data         A \code{matrix}.
#' @param miss.mat     A logical \code{matrix} where cells are TRUE when there
#'                     is a numeric value in data, or FALSE when there are NAs
#'                     in data. miss.mat is only used when vinterpolate is TRUE
#'                     (Default: miss.mat = NULL).
#' @param grp_tbl      A \code{dataframe} with 2 columns:
#'                     \itemize{
#'                      \item{samples - which contains all column names from
#'                      data.}
#'                      \item{groups - which contains the groups to which each
#'                      column belong to.}
#'                     }
#' @param nr           An \code{integer} specifying the total number of rows in
#'                     the matrix data.
#' @param r            An \code{integer} specifying the position of the row of
#'                     interest for imputation in the matrix data.
#' @param grp.fun      A \code{character} specifying the function to be used on
#'                     values from a same group
#'                     (Default: impute.fun = "median").
#' @param impute.fun   A \code{character} specifying the function to apply to
#'                     all groups values obtained using 'grp.fun' for imputation
#'                     (Default: impute.fun = "median").
#' @param vinterpolate A \code{logical} specifying whether vertical
#'                     interpolation should be applied on missing values, based
#'                     on previous and next available values in a column
#'                     (vinterpolate = TRUE) or not
#'                     (Default: vinterpolate = FALSE). It is combined with the
#'                     function use for imputation by group. vinterpolate is not
#'                     use if method is not set to 'impute'.
#' @return A \code{numeric} vector containing the imputed values of the row for
#'         each groups of columns.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

row.impute.na <- function(
  data, miss.mat = NULL, grp_tbl, nr = NULL, r, grp.fun = "median",
  impute.fun = "median", vinterpolate = FALSE){
  sub.data <- data[r, ]
  # Apply group function to each group to get 1 value per group
  grp_vals <- unlist(lapply(X = unique(grp_tbl$groups), FUN = function(grp){
    samples.in.grp <- grp_tbl[grp_tbl$groups == grp, ]$samples
    grp.data <- sub.data[samples.in.grp]
    if(all(is.na(grp.data))){
      message(paste(
        "grp.fun failed on row", r,"in group", grp, "because of missing data."))
      grp_val <- NA
      names(grp_val) <- grp
      grp_val
    } else {
      grp_val <- eval(parse(text = paste0(
        grp.fun, "(x = grp.data, na.rm = TRUE)")))
      names(grp_val) <- grp
      grp_val
    }
  }))
  # Impute missing data
  if(vinterpolate){
    impute.vals <- unlist(lapply(
      X = names(sub.data[is.na(sub.data)]), FUN = function(i){
        grp.val <- eval(parse(text = paste0(
          impute.fun, "(x = grp_vals, na.rm = TRUE)")))
        lead.val <- rev(data[which(miss.mat[seq(r - 1), i]), i])[1]
        if(r < nr){
          lag.val <- data[r + which(miss.mat[(r + 1):nr, i]), i][1]
        } else { lag.val <- NA }
        if(is.na(lead.val)){
          grp_tbl[grp_tbl$samples == i, ]$groups
          impute.val <- eval(parse(text = paste0(
            impute.fun, "(c(grp.val, lag.val), na.rm = TRUE)")))
        } else if(is.na(lag.val)){
          impute.val <- eval(parse(text = paste0(
            impute.fun, "(c(grp.val, lead.val), na.rm = TRUE)")))
        } else {
          impute.val <- eval(parse(text = paste0(
            impute.fun, "(c(grp.val, lead.val, lag.val), na.rm = TRUE)")))
        }
        names(impute.val) <- i
        impute.val
      }))
  } else {
    impute.vals <- rep(eval(parse(text = paste0(
      impute.fun, "(x = grp_vals, na.rm = TRUE)"))), length(
        data[r, ][is.na(data[r, ])]))
    names(impute.vals) <- names(data[r, ][is.na(data[r, ])])
  }
  return(impute.vals)
}


#' Keeps, removes or imputes missing at random (MAR) values in a matrix or a
#' data.frame accounting for group size.
#'
#' @param data         A \code{matrix} or \code{data.frame} with column names.
#' @param method       A \code{character} which defines the method to use for
#'                     managing NAs:
#'                     \itemize{
#'                      \item{'keep' does nothing.}
#'                      \item{'remove' removes all rows containing at least 1
#'                      NA.}
#'                      \item{'impute' removes rows containing only NAs and
#'                      imputes missing at random values (MARs) in remaining
#'                      rows. The imputation uses the groups defined in the
#'                      parameter 'groups': e.g. in row 1 if the value of
#'                      sample_1 is NA, the missing value will be imputed
#'                      calculating the median of all sample groups' medians,
#'                      including the sample group in which sample_1 is.}
#'                     }
#' @param groups       A vector of length ncol(data) specifying the groups to
#'                     which samples belong. Order of groups in the vector has
#'                     to match the order of column names in data to properly
#'                     associate samples and groups.
#' @param grp.fun      A \code{character} specifying the function to be used on
#'                     values from a same group
#'                     (Default: impute.fun = "median").
#' @param impute.fun   A \code{character} specifying the function to apply to
#'                     all groups values obtained using 'grp.fun' for imputation
#'                     (Default: impute.fun = "median").
#' @param vinterpolate A \code{logical} specifying whether vertical
#'                     interpolation should be applied on missing values, based
#'                     on previous and next available values in a column
#'                     (vinterpolate = TRUE) or not
#'                     (Default: vinterpolate = FALSE). It is combined with the
#'                     function use for imputation by group. vinterpolate is not
#'                     use if method is not set to 'impute'.
#' @param row.na.cut   An \code{integer} specifying the cut-off for imputation
#'                     on rows missing values. if NULL, all rows with missing
#'                     values will be processed. if row.na.cut equals an
#'                     integer, only rows containing a number of NAs smaller or
#'                     equal to the cut-off will be processed
#'                     (Default: row.na.cut = NULL).
#' @param ncores       An \code{integer} to specify the number of cores/threads
#'                     to be used to parallel-process the matrix
#'                     (Default: ncores = 1).
#' @return A \code{matrix}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Creating a data.frame of 4 columns with 50 values in each:
#' df <- data.frame("col1" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col2" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col3" = rnorm(n = 50, mean = 25, sd = 4),
#'                  "col4" = rnorm(n = 50, mean = 25, sd = 4))
#' # col1 and col4 contain missing data
#' df$col1[1:10] <- NA
#' df$col4[40:50] <- NA
#' df
#' # Create 2 groups: grp1 & grp2.
#' # col1 and col2 are in grp1, col3 and col4 are in grp2.
#' grps <- c("grp1","grp1","grp2","grp2")
#' # Remove all rows containing at least 1 missing value
#' manage.na(data = df, method = "remove", groups = grps)
#' # Impute by group the missing data (and remove rows missing all values)
#' manage.na(data = df, method = "impute", groups = grps)
#' @references \href{https://stefvanbuuren.name/fimd/sec-MCAR.html}{Concepts of MCAR, MAR and MNAR - Stef van Buuren}

manage.na <- function(
  data, method = "remove", groups, grp.fun = "median", impute.fun = "median",
  vinterpolate = FALSE, row.na.cut = NULL, ncores = 1){
  if(is.data.frame(data)){ data <- as.matrix(data) }
  if(method != "keep"){ #If NA should not be kept
    if(any(is.na(data))) { #If any NA
      if(method == "remove"){ #Remove all rows containing any NA
        data_rname <- rownames(data)[
          apply(X = data, MARGIN = 1, FUN = anyNA) == FALSE]
        data <- data[stats::complete.cases(data), ]
        if(is.vector(data)){
          data <- matrix(
            data, nrow = 1, dimnames = list(data_rname, names(data)))
        }
      } else if(method == "impute"){ #Impute NAs with the median value by group
        #Parallel-remove rows containing only NAs
        data <- data[!unlist(parallel::mclapply(
          X = seq(nrow(data)), mc.cores = ncores, FUN = function(r){
            all(is.na(data[r, ])) })), , drop = FALSE]
        #Create logical matrix of missing value for vinterpolate
        if(vinterpolate){
          miss.mat <- !is.na(x = data)
          nr <- nrow(data)
        } else {
          miss.mat <- NULL
          nr <- NULL
        }
        # Precompute rowSums(is.na(data))
        pos.cov <- as.integer(rowSums(is.na(data)))
        if(!is.null(row.na.cut)){
          rows.to.impute <- which(pos.cov <= row.na.cut & pos.cov > 0)
        } else { rows.to.impute <- which(pos.cov > 0) }
        #Get groups of samples from sample conditions
        grp_tbl <- data.frame(samples = colnames(data), groups = groups)
        #TODO: Make parallel apply when it will be possible.
        invisible(lapply(X = rows.to.impute, FUN = function(r){
          grp.imp <- row.impute.na(
            data = data, miss.mat = miss.mat, grp_tbl = grp_tbl, nr = nr, r = r,
            grp.fun = grp.fun, impute.fun = impute.fun,
            vinterpolate = vinterpolate)
          base::`<<-` (data[r, names(grp.imp)], grp.imp)
          cat(paste0(round((r/max(rows.to.impute))*100, 3), "%\n"))
        }))
        if(vinterpolate){ rm(miss.mat) }
      }
    }
  }
  return(data)
}
