##IMPORTS
Imports = c("ggplot2","plyr","IRanges")
lapply(Imports, library, character.only = T)


##DATA 
pal_sunset<-c("red","gold","blue4")
pal_westworld<-c("sienna1","lightgoldenrod","skyblue3")
pal_startrek<-c("red","goldenrod1","dodgerblue")
pal_margesimpson<-c("lightgreen","gold","dodgerblue2")
pal_eighties<-c("darkviolet","deeppink","blue4")
pal_insta<-c("deeppink","red","goldenrod1")

##FUNCTIONS

# sunset #######################################################################

#' @description draw a sunset plot showing the completeness of a dataset.
#' 
#' @param pos.cov          An \code{integer} vector containing the sum of
#'                         non-0/non-NA values by rows of the dataset.
#' @param title            A \code{character} to specify the title of your plot.
#' @param col.pal          A \code{character} vector of length 3 matching R
#'                         colors to used as a palette for the plot.
#'                         Alternatively col.pal can use the palettes provides
#'                         within the script of the function:
#'                         pal_sunset (default), pal_westworld, pal_startrek,
#'                         pal_margesimpson, pal_eighties, pal_insta. 
#' @param horizontal       A \code{logical} to specify whether the plot should
#'                         be drawn vertically (default) or horizontally.
#' @param keep_2nd_ticks   A \code{logical} to specify whether the ticks of the
#'                         secondary axis should be displayed or not (default).
#' @param n.grad           An \code{integer} specifying the number of graduation
#'                         to use on the Y-Axis (Default: n.grad = 15).
#' @param display.cutoff   A \code{double} specifying the the minimum height of
#'                         a bin to display the label of the number of samples
#'                         on it.
#' @param display.num.smpl A \code{double} specifying the minimum height of a
#'                         bin to display the label of the number of samples on
#'                         its side.
#' @param display.sep      A \code{double} specifying the linewidth of the
#'                         separations between the bins.
#' @param lgd.pos          A \code{character} specifying the position of the
#'                         legend (for more information about the possible
#'                         values, check the theme() function from the ggplot2
#'                         package for the parameter "legend.position").
#'  
#' @value a \code{gg} plot of the sunset.
#' @author Yoann Pageaud.

sunset<-function(pos.cov, title, col.pal = pal_sunset, horizontal = F,
                 keep_2nd_ticks = F, n.grad = 15, display.cutoff = 0.03,
                 display.num.smpl = 0.01, display.sep = 0.005,
                 lgd.pos = "bottom",reverse=F){
  
  N <- max(pos.cov)
  x <- tabulate(pos.cov + 1L, N + 1L)
  smpl_string<-paste0(c("NONE",seq(N))," (",round((x/sum(x))*100),"%",")")
  Sample.Amt <-c(0,seq(N))
  
  if(reverse){
    x<-rev(x)
    smpl_string<-rev(smpl_string)
    Sample.Amt<-rev(Sample.Amt)
  }
  
  #Create dataframes
  dframe<-data.frame(Sample.Amt = Sample.Amt,
                     CpG.Covered = x,
                     label.pos = cumsum(x)-0.5*x,
                     percent = smpl_string,
                     white_lines = cumsum(x),
                     diff_bins = x/sum(x),
                     cumulated = cumsum(x),
                     stringsAsFactors = F)
  
  df_right_y<-data.frame(CpG.Covered = x,
                         cumulated_smpl = cumsum(x)-0.5*x,
                         right_Y = smpl_string,
                         diff_bins = x/sum(x),
                         stringsAsFactors = F)
  
  dframe[dframe$CpG.Covered < display.cutoff * max(dframe$CpG.Covered),
         ]$percent<-" "
  dframe[dframe$CpG.Covered/sum(x) < display.sep,]$white_lines<-NA
  df_right_y[df_right_y$CpG.Covered > display.cutoff * max(dframe$CpG.Covered),
             ]$right_Y<-" "
  df_right_y[df_right_y$diff_bins < display.num.smpl,]$right_Y<-" "
  
  #Get positions of strings that will not be displayed
  pos_str<-match(smpl_string[-(sort(match(c(df_right_y$right_Y,dframe$percent),
                                          smpl_string)))],smpl_string)
  #Get groups of following samples
  smpl_intervals<-IRanges()
  for(i in pos_str){
    if(i+1 - i == 1){
      if(!is.na(match(i+1,pos_str))){
        smpl_intervals<-c(smpl_intervals,IRanges(start = i, end = i+1))
      }
    }
  }
  smpl_intervals<-reduce(smpl_intervals)
  #Get groups average label.pos cumulative diff_bins
  if(length(smpl_intervals) != 0){
    for (i in seq(length(smpl_intervals))) {
      df_grp<-dframe[c(start(smpl_intervals[i]):end(smpl_intervals[i])),
                     c(1:3,6,7)]
      cpg_cov_grp<-sum(df_grp$CpG.Covered)
      per_grp<-sum(df_grp$diff_bins)
      str_grp<-paste0(df_grp[as.character(start(smpl_intervals[i])),1]," to ",
                      df_grp[as.character(end(smpl_intervals[i])),1]," (",
                      round(per_grp*100),"%)")
      if(rownames(head(df_grp,n=1L))=="1"){
        pos_lab_grp<-(0 + tail(df_grp$cumulated,n=1L))/2
      } else {
        pos_lab_grp<-(dframe[as.integer(rownames(head(df_grp,n=1L)))-1,
                             ]$cumulated + tail(df_grp$cumulated,n=1L))/2  
      }
      if(per_grp > display.num.smpl){
        df_right_y<-rbind(df_right_y,data.frame("CpG.Covered" = cpg_cov_grp,
                                                "cumulated_smpl" = pos_lab_grp,
                                                "right_Y" = str_grp,
                                                "diff_bins" = per_grp))
      }
    }  
  }
  
  if(keep_2nd_ticks == F) {
    df_right_y<-df_right_y[df_right_y$right_Y != " ",]  
  }
  
  #Vector for white lines
  white_lines<-dframe[!is.na(dframe$white_lines),]$white_lines
  
  
  #"Sunset" Plot of the Amount of CpGs Covered by Number of Samples
  Sunset<-ggplot(data=dframe,aes(x = 0,y=CpG.Covered, fill = Sample.Amt)) +
    theme_gray() +
    theme(legend.title.align=0.5,
          legend.text=element_text(size=12)) +
    geom_bar(stat = "identity") +
    geom_text(aes(y=label.pos, label=percent), vjust=0.25, hjust=0.5,
              colour = "white", size = 5) +
    scale_fill_gradient2(low=col.pal[1],mid=col.pal[2],high=col.pal[3],
                         midpoint = round(N/2)) +
    scale_y_continuous(expand = c(0,0),
                       breaks =
                         seq(0, sum(dframe$CpG.Covered),
                             by = round_any(sum(dframe$CpG.Covered)/n.grad,
                                            10^(nchar(as.character(
                                              as.integer(round(
                                                sum(dframe$CpG.Covered)/n.grad,
                                                1))))-1))),
                       sec.axis = sec_axis(trans = ~.,
                                           breaks = df_right_y$cumulated_smpl,
                                           labels = df_right_y$right_Y)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_hline(yintercept=white_lines,color="white")
  
  if(horizontal == T) {
    Sunset <- Sunset +
      ggtitle(title) +
      theme(plot.title = element_text(title,
                                      size = 15,hjust = 0.5),
            axis.title = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 12,face = "bold"),
            legend.title = element_text(size = 14,vjust = 0.8),
            legend.position=lgd.pos,
            plot.margin = margin(0,0,0,20)) +
      labs(fill = "Number of Samples") +
      coord_flip()
  } else {
    Sunset <- Sunset +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12,face = "bold"),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 14),
            legend.position=lgd.pos) +
      ylab(title) +
      labs(fill = "Number\nof\nSamples")
  }
  Sunset
}