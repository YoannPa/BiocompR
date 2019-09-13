
##IMPORTS
Imports = c("data.table")
lapply(Imports, library, character.only = T)

##FUNCTIONS

# get.ks.stat ##################################################################

#' @description Function description.
#' 
#' @param df.ks.tests A \code{data.frame} containing raw results from multiple
#'                    Kolmogorov-Smirnov tests.
#' @param statistic   A \code{character} specifying the type of statistic to
#'                    retrieve from the test
#'                    (Supported: statistic = c('n','stat','p')).
#' @value a \code{matrix} containing the KS test's statistic values of interest.
#' @author Yoann Pageaud.

get.ks.stat<-function(table_combinations, df.ks.tests, statistic){
  ks.res<-df.ks.tests[[statistic]]
  tbl.ks<-cbind(table_combinations,ks.res)
  #Cast a molten tables into a matrix
  mat.ks<-dcast(data=tbl.ks, formula=Var1 ~ Var2, value.var="ks.res")
  rownames(mat.ks)<-mat.ks$Var1
  mat<-as.matrix(mat.ks[,-1])
  return(mat)
}

# pairwise.ks ##################################################################

#' @description Compute a Kolmogrov-Smirnov test between all columns of a
#'              data.frame.
#' 
#' @param data      A \code{data.frame} of numerical values.
#' @param statistic A \code{character} specifying the type of statistic to
#'                  retrieve from the test
#'                  (Supported: statistic = c('n','stat','p')).
#' @param ncores    An \code{integer} to specify the number of cores/threads to
#'                  be used to parallel-run tests.
#' @value A \code{list} containing a matrix of the statistic values, a
#'        dataframe of the pairwise KS raw test results, and a table of pairwise
#'        combinations.
#' @author Yoann Pageaud.

pairwise.ks<-function(data, statistic, ncores){
  table_combinations<-expand.grid(colnames(data),colnames(data))
  List_ks.tests<-mclapply(
    seq(nrow(table_combinations)),mc.cores = ncores,function(i){
      #Compute KS test
      ks.res<-ks.test(x = data[,table_combinations[i,1]],
                      y = data[,table_combinations[i,2]])
      #Create table with all statistics of the KS test
      data.frame('n' = nrow(data[complete.cases(
        data[,c(table_combinations[i,1], table_combinations[i,2])]),]),
        'stat' = ks.res$statistic, 'p' = ks.res$p.value, row.names = NULL)
    })
  #Create Table
  df.ks.tests<-do.call("rbind",List_ks.tests)
  #Get matrix of the stat of interest
  mat<-get.ks.stat(table_combinations = table_combinations,
                   df.ks.tests = df.ks.tests, statistic = statistic)
  return(list("res.statistic" = mat, "res.test" = df.ks.tests,
              "table_combinations" = table_combinations))
}
