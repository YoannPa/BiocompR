##CpG Coverage Summary

# setwd("/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Hematopoiesis_DNMT1_Hypomorph")
setwd("/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Prostate Cancer Projects/")

## L I B R A R I E S ###########################################################
library(ggplot2)
library(plyr)
library(IRanges)
## G L O B A L S ###############################################################

## File storing the summary objects
# FILE.SUMMARY <- paste0("data/","CpG_Table_Merged_DNMT1.RDS")
FILE.SUMMARY <- paste0("CpG_Sites_Table/","PC_CpG_Methylation_Table.RDS")

## Load the statistics table
tbl <- readRDS(FILE.SUMMARY)
tbl$pos.cov <- as.integer(tbl$pos.cov)
tbl$pos.cov.ten <- as.integer(tbl$pos.cov.ten)
# av.cov: average coverage
# pos.cov: number of samples with positive coverage
# av.pos.cov: average coverage amongth samples with pos. coverage
# min.cov: minimum coverage
# fith.perc.cov: fifth percentile coverage
# fiftyth.perc.cov: median coverage
# ninetyfith.perc.cov: 95th percentile coverage
# pos.cov.ten: number of samples with coverage >= 10
# av.meth: average methylation (among samples with positive coverage)
# sd.meth: SD of methylation (among samples with positive coverage)
rm(FILE.SUMMARY)

#Choose horizontal format (TRUE or FALSE) ?
horizontal <-F
#Keep empty Y axis ticks (TRUE or FALSE) ?
keep_y_ticks <-F
#Set an integer for the number of graduations on Y axis
n.graduation<-15
#Set a cutoff for the display of the percentage
display.cutoff<-0.03
#Set a cutoff for displaying the number of samples on the right Y-axis
display_num_smpl<- 0.01
#Set a cutoff for displaying horizontal white separation
display_sep<-0.005
#Set legend position ("top","left","bottom","right")
lgd_pos<-"bottom"
  
#Palettes available
pal_sunset<-c("red","gold","blue4")
pal_westworld<-c("sienna1","lightgoldenrod","skyblue3")
pal_startrek<-c("red","goldenrod1","dodgerblue")
pal_margesimpson<-c("lightgreen","gold","dodgerblue2")
pal_eighties<-c("darkviolet","deeppink","blue4")
pal_instagram<-c("deeppink","red","goldenrod1")

used_pal<-pal_sunset


## Summary of which CpGs are covered
N <- max(tbl$pos.cov)
x <- tabulate(tbl$pos.cov + 1L, N + 1L)

#Strings to print on the plot
smpl_string<-paste0(c("NONE",seq(N-1),"ALL")," (",round((x/sum(x))*100),"%",")")


#Create dataframes
dframe<-data.frame(Sample.Amt = c(0,seq(N)),
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

dframe[dframe$CpG.Covered < display.cutoff * max(dframe$CpG.Covered),]$percent<-" "
dframe[dframe$CpG.Covered/sum(x) < display_sep,]$white_lines<-NA
df_right_y[df_right_y$CpG.Covered > display.cutoff * max(dframe$CpG.Covered),]$right_Y<-" "
df_right_y[df_right_y$diff_bins < display_num_smpl,]$right_Y<-" "

#Get positions of strings that will not be displayed
pos_str<-match(smpl_string[-(sort(match(c(df_right_y$right_Y,dframe$percent),smpl_string)))],smpl_string)
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
for (i in seq(length(smpl_intervals))) {
  df_grp<-dframe[c(start(smpl_intervals[i]):end(smpl_intervals[i])),c(1:3,6,7)]
  cpg_cov_grp<-sum(df_grp$CpG.Covered)
  per_grp<-sum(df_grp$diff_bins)
  str_grp<-paste0(df_grp[as.character(start(smpl_intervals[i])),1]," to ",
                  df_grp[as.character(end(smpl_intervals[i])),1]," (",
                  round(per_grp*100),"%)")
  pos_lab_grp<-(dframe[as.integer(rownames(head(df_grp,n=1L)))-1,7] +
                  tail(df_grp$cumulated,n=1L))/2
    if(per_grp > display_num_smpl){
      df_right_y<-rbind(df_right_y,data.frame("CpG.Covered" = cpg_cov_grp,
                                              "cumulated_smpl" = pos_lab_grp,
                                              "right_Y" = str_grp,
                                              "diff_bins" = per_grp))
      }
}


if(keep_y_ticks == F) {
  df_right_y<-df_right_y[df_right_y$right_Y != " ",]  
}

#Vector for white lines
white_lines<-dframe[!is.na(dframe$white_lines),]$white_lines


#"Sunset" Plot of the Amount of CpGs Covered by Number of Samples
Sunset<-ggplot(data=dframe, aes(x = 0,y=CpG.Covered, fill = Sample.Amt)) +
  theme_gray() +
  theme(legend.title.align=0.5,
        legend.text=element_text(size=12)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y=label.pos, label=percent), vjust=0.25, hjust=0.5, colour = "white", size = 5) +
  scale_fill_gradient2(low=used_pal[1],mid=used_pal[2],high=used_pal[3],midpoint = round(N/2)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, sum(dframe$CpG.Covered),
                                  by = round_any(sum(dframe$CpG.Covered)/n.graduation,
                                                 10^(nchar(as.character(as.integer(round(sum(dframe$CpG.Covered)/n.graduation,
                                                                                         1))))-1))),
                     sec.axis = sec_axis(trans = ~.,breaks = df_right_y$cumulated_smpl,labels = df_right_y$right_Y)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept=white_lines,color="white")

if(horizontal == T) {
  Sunset <- Sunset +
    ggtitle("Amount of CpGs Covered in Samples") +
    theme(plot.title = element_text("Amount of CpGs Covered in Samples",
                                    size = 15,hjust = 0.5),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12,face = "bold"),
          legend.title = element_text(size = 14,vjust = 0.8),
          legend.position=lgd_pos,
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
          legend.position=lgd_pos) +
    ylab("Amount of CpGs Covered in Samples") +
    labs(fill = "Number\nof\nSamples")
}

Sunset
