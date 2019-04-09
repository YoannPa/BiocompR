
##IMPORTS
Imports = c('data.table','ggplot2','Hmisc','fastcluster','parallelDist',
            'fastcluster','RColorBrewer')
lapply(Imports, library, character.only = T)

##FUNCTION

# ggclust.heatmap ##############################################################

#' Create a clustered heatmap using GGplot.
#' 
#' @param data             A \code{matrix} or a \code{data.frame} with samples
#'                         by columns and regions by rows. Row names have to be
#'                         of the format "chr:start-end" and matching at least
#'                         some regions clustered.  
#' @param regions.clusters A \code{data.frame} of BED format with the cluster
#'                         IDs (chr, start, end, clusterID).
#' @param order.type       A \code{character} specifying how regions should be
#'                         ordered within a cluster. 2 possible methods:
#'                         SD - ordered by increasing order of standard
#'                         deviation
#'                         Hclust - ordered following a hierarchical clustering
#'                         method.
#' @param distance.method  A \code{character} specifying the distance
#'                         calculation method if regions should be ordered
#'                         following a hierarchical clustering method.
#'                         Distance methods for continuous input variables:
#'                         bhjattacharyya, bray, canberra, chord, divergence,
#'                         dtw, euclidean, fJaccard, geodesic, hellinger,
#'                         kullback, mahalanobis, manhattan, maximum, minkowski,
#'                         podani, soergel, wave, whittaker.
#'                         Distance methods for binary input variables:
#'                         binary, braun-blanquet, cosine, dice, fager, faith,
#'                         hamman, hamming, kulczynski1, kulczynski2, michael,
#'                         mountford, mozley, ochiai, phi, russel,
#'                         simple matching,  simpson, stiles, tanimoto, yule,
#'                         yule2.
#' @param select.clust     A \code{character} vector specifying the clusters to
#'                         keep for the plot.
#' @param select.sample    A \code{character} vector specifying the samples to
#'                         keep for the plot (must match the columns names of
#'                         the matrix). 
#' @param groups           A \code{character} vector specifying the name of the
#'                         groups for samples. The length of the vector must
#'                         match the number of columns in the matrix.
#' @param groups.order     A \code{character} vector specifying the names of the
#'                         groups (must be unique) in the wanted order.
#' @param sampling         A \code{logical} to specify if a random sampling
#'                         should be applied on each cluster.
#' @param sampling.size    An \code{integer} specifying the number of regions to
#'                         randomly select from each cluster if sampling = TRUE.
#'                         If one of the clusters contains less regions than the
#'                         amount specified, all regions of the cluster will be
#'                         selected.
#' @param na.rm            A \code{logical} to specify if NAs should be removed
#'                         or not from the data. If order.type = "Hclust", NAs
#'                         are removed anyway.
#' @param title            A \code{character} to be used as the plot title.
#' @param palette          A \code{character} vector to be used as a color
#'                         palette to plot the data.
#' @param na.col           A \code{character} matching a color to used for NAs
#'                         on the plot.
#' @param margins          A \code{double} vector of length 4 to give the
#'                         thickness of each margin of the plot.
#' @return A \code{gg} plot.
#' @author Yoann Pageaud.

ggclust.heatmap<-function(data, regions.clusters, order.type="Hclust",
                          distance.method, select.clust = NULL,
                          select.sample = NULL,
                          groups=NULL,groups.order=NULL,
                          sampling = F,sampling.size = 500, na.rm = F, title="",
                          palette = c("steelblue", "gray95", "darkorange"),
                          na.col = "black", margins = c(0.5,0,0,0)){
  #Concatenate regions chr + start + end
  regions.clusters$V2<-apply(regions.clusters[,c(2,3)],1,paste0,collapse = "-")
  regions.clusters$V1<-as.integer(gsub("[^0-9.]",replacement ="",
                                       regions.clusters$V1))
  regions.clusters<-regions.clusters[order(regions.clusters$V1,
                                           regions.clusters$V3),]
  regions.clusters$V1<-as.character(regions.clusters$V1)
  regions.clusters$V2<-apply(regions.clusters[,c(1,2)],1,paste0,collapse = ":")
  regions.clusters<-regions.clusters[,c(2,4)]
  
  #Only Consider the regions that have been clustered.
  data<-as.data.frame(data[rownames(data) %in% regions.clusters$V2,])
  data<-setDT(data, keep.rownames = TRUE)
  data[, Clusters := regions.clusters$V4]
  
  #Remove all NAs
  if(na.rm){
    data<-na.omit(data)
  } else {
    if(order.type == "Hclust"){ #Hclust need NAs to be removed anyway
      na.rm <- T
      data<-na.omit(data)
    }
  }
  
  #Random Sampling by Clusters if sampling TRUE
  if(sampling){
    data<-data[,.SD[sample(.N, min(sampling.size,.N))],by = Clusters]
  }
  #Keep Clusters of interest
  if(is.null(select.clust) == F) {
    data<-data[Clusters %in% select.clust]  
  }
  #Keep Samples of interest
  if(is.null(select.sample) == F){
    cols<-c("rn","Clusters",select.sample)
    data<-data[,..cols]  
  }
  
  #Order Regions in clusters
  if(order.type == "SD"){ #Order regions by cluster and then by SD
    data[, sd := apply(data[, 3:(ncol(data))], 1, sd, na.rm = T)]
    data<-data[order(as.numeric(Clusters), -sd)]
    data<-data[,sd:=NULL]
    row_count<-nrow(data)
    distance.method<-"none"
    melt_data<-melt.data.table(data,id.vars = c("rn","Clusters"))
    
  } else {#Order regions by cluster and then by hierarchical cluster
    row_count<-nrow(data)
    data_list<-split(data,f = data$Clusters)
    #Remove Clusters and and sd columns from all dataframes in list
    data_list_names<-lapply(data_list, function(x) x[["rn"]])
    # data_list<-lapply(data_list, function(x) x[,c("rn","Clusters"):=NULL])
    lapply(data_list, function(x) x[,c("rn","Clusters"):=NULL])
    data_list<-lapply(data_list,as.matrix)
    row_dist<-lapply(data_list,parDist,method = distance.method,threads = 7)
    row_hclust<-lapply(row_dist, hclust)
    Dend_list<-lapply(row_hclust, as.dendrogram)
    row.order<-lapply(Dend_list, order.dendrogram)
    namesdata<-names(data_list)
    data_list<-lapply(1:length(data_list), function(i){
      xi = data_list[[i]]
      xi = xi[row.order[[i]],, drop = FALSE]
      xi
    })
    names(data_list)<-namesdata
    data_list<-lapply(data_list,as.data.frame)
    data_list<-lapply(data_list,setDT, keep.rownames = TRUE)
    data_list<-suppressWarnings(lapply(data_list,melt.data.table))
    data_list<-Map(f = cbind, data_list,Clusters = as.integer(names(data_list)))
    melt_data<-do.call(rbind,data_list)
  }
  #Add groups
  if(is.null(groups) == F){
    split_melt<-split(melt_data,f = melt_data$variable)
    split_melt<-Map(cbind,split_melt,Groups = groups)
    melt_data<-do.call(rbind,split_melt)
    melt_data$Groups<-factor(melt_data$Groups,levels = groups.order)  
  }
  
  #Plot Heatmap Cluster
  if(na.rm){
    Clustheatmap<- ggplot(data = melt_data, aes(x = variable, y = rn,
                                                fill = value)) +
      scale_fill_gradientn(colors = palette)
  } else {
    Clustheatmap<- ggplot(data = melt_data, aes(x = variable, y = rn,
                                                fill = value,colour="")) +
      scale_fill_gradientn(colors = palette,na.value = na.col) +
      scale_colour_manual(values=NA) +
      guides(colour=guide_legend("N.A.", override.aes=list(fill=na.col)))
  }
  Clustheatmap<-Clustheatmap +
    geom_tile() +
    theme(axis.title = element_text(size = 13), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11, angle = 90, hjust = 1,
                                     vjust = 0.5, face = "bold"),
          plot.margin=unit(margins,"cm"),
          legend.text=element_text(size= 11),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 12,face = "bold"),
          strip.text.y = element_text(angle = 180),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste(title,"\n",row_count,"regions ordered by",order.type,
                       "using",capitalize(distance.method),"distance."),
         x = "Samples", y = "Clusters")
  
  if(is.null(groups) == F){
    Clustheatmap <- Clustheatmap +
      facet_grid(Clusters ~ Groups, scales = "free", switch="y")
  } else {
    Clustheatmap <- Clustheatmap +
      facet_grid(Clusters ~ ., scales = "free", switch="y")
  }
  Clustheatmap
}
