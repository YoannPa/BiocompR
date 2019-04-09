#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
################################_Plot_Heatmaps_#################################
Version = '1.0'
Date = '2018-07-26'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 3.5.2 (2018-12-20)',
'RStudio version 1.1.453 – © 2009-2018','circlize','reshape2','parallelDist',
'fastcluster','ggdendro','ggplot2','grid','gridExtra','Hmisc','gtable','ggsci')
Description = 'Plot Methylation Heatmaps with various options:
Missing values imputation; Several distance hierarchical clustering methods;
Row ordering; Ranking or Overviewing DMRs/DMLs; NA handling;
Plotting dendrograms.'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##Imports
Imports = c("circlize","reshape2","parallelDist","fastcluster","ggdendro",
            "ggplot2","grid","gridExtra","Hmisc","gtable","ggsci","RColorBrewer")
lapply(Imports, library, character.only = T)

##Parameters
#Load Matrix (ALL)
DMRs<-readRDS('/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Hematopoiesis_DNMT1_Hypomorph/data/DML_Filt_Merged_HSPC_DMRs_NoSNPs_DMRCov30_Dens3_DMRpercent02_P-val005.RDS')

#Specify which samples to display
# samples.selected<-seq(ncol(DMRs)) #Load Full Matrix
# samples.selected<-c(1,2,4,5,7,8,10,11,13,15) #Load Matrix (Mut & Normal)
# samples.selected<-c(2,3,5,6,8,9,11,12,14,16) #Load Matrix (Mut & WT)
samples.selected<-c(1,3,4,6,7,9,10,12:16) #Load Matrix (Normal & WT)

#Choose how NAs should be handled if there are some ("keep","impute","remove")
#"keep": keep all NAs;
#"impute": assign missing value the median of its group. if no value remove row;
#"remove": remove all rows containing any NA;
handle_na<-"remove"

#Select a distance measure method for columns
dist_methods.cols<-c("euclidean","maximum","manhattan",
                "canberra","binary","minkowski")
selected_dist.cols<-dist_methods.cols[3]

#Select a distance measure method for columns
dist_methods.rows<-c("euclidean","maximum","manhattan",
                     "canberra","binary","minkowski")
selected_dist.rows<-dist_methods.rows[3]

#Select ordering method for columns of the heatmap
col_order <- c("distance method","cell type","genotype")
selected_col_order <- col_order[1]

#Select ordering method for rows of the heatmap
row_order <- c("variability across samples","distance method","default")
selected_row_order <- row_order[2]

#Select heatmap type ("Rank","Overview")
htmp_type<-"Overview"

#Set number of DMRs to display
dmr_rank<-50

#Should the dendrogram of columns be displayed (TRUE, FALSE)
display_dend_col<-T

#Should the dendrogram of rows be displayed (TRUE, FALSE)
display_dend_row<-T

#Genotype Legend Order
# legend.genotype.order<-c(2,3,1)
legend.genotype.order<-c(1,2,3)

#Select the ggsci palette of your choice with
# palette.genotypes<-pal_npg("nrc",alpha = 0.9)(10)[c(1,2,3)] # Nature Reviews Cancer
palette.genotypes<-pal_npg("nrc",alpha = 0.9)(10)[c(3,1,2)] # Nature Reviews Cancer

#palette.genotypes<-pal_lancet("lanonc",alpha = 0.7)(9) # Lancet Oncology
#palette.genotypes<-pal_jco("default",alpha = 0.9)(10) # Journal of Clinical Oncology

#Select Palette for values
# palette.heatmap<-c("steelblue", "gray95", "darkorange") #Default palette
palette.heatmap<-rev(brewer.pal(n = 11,name = "RdBu"))

##_Extract Cell types and Genotypes_############################################

#Get cell types/lines
samples_split<-sapply(colnames(DMRs), strsplit,split="_")
Cell_types<-unlist(lapply(samples_split,head,n=1))

#Get genotypes
Genotypes<-lapply(samples_split, tail, n=2L)
Genotypes<-unlist(lapply(Genotypes, paste0, collapse=" "))
Genotypes[!grepl("chip", Genotypes, ignore.case=FALSE)] <- "C57BL_6J"
Genotypes[grepl("wt", Genotypes, ignore.case=FALSE)] <- "129P2_OlaHsd_WT"
Genotypes[grepl("mut", Genotypes, ignore.case=FALSE)] <- "129P2_OlaHsd_Mut"

Geno_short<-unlist(lapply(sapply(Genotypes, strsplit, "_"),tail,1))
Geno_short[grepl("6J", Geno_short, ignore.case=FALSE)] <- "Normal"
#Set samples conditions
Sample_cond<-paste(Cell_types,Genotypes,sep = "_")

#Create Sample Table
sample_tbl<-data.frame("Genotypes" = Genotypes, "Short.genotypes" = Geno_short,
					   "Cell.types" = Cell_types,
					   "Sample.Conditions" = Sample_cond,
					   "Order" = seq(ncol(DMRs)))

##_Create genotype-color table_#################################################

#Initialize table
color_table<-data.frame("Genotypes" = unique(sample_tbl$Short.genotypes)[legend.genotype.order],
                        "Colors" = palette.genotypes,stringsAsFactors = F)

#subset samples from table
subset.smpl<-colnames(DMRs)[samples.selected]

#subset genotypes from table
subset.genotypes<-unique(sample_tbl[subset.smpl,"Short.genotypes"])

#Update color table
color_table<-color_table[color_table$Genotypes %in% subset.genotypes,]

#Subset sample table by the samples selected
sample_tbl<-sample_tbl[samples.selected,]

#Subset Matrix with selected samples
DMRs<-DMRs[,samples.selected]

##_Handle NAs_##################################################################
if(handle_na != "keep"){ #If NA should not be kept
	
	if(any(is.na(DMRs))) { #If any NA
		if(handle_na == "remove"){ #Remove all rows containing any NA
			DMRs<-DMRs[complete.cases(DMRs),]
			
		} else { #Impute NAs with the median value by group 
			#Remove rows containing only NAs
			DMRs<- DMRs[!apply(is.na(DMRs), 1, all), ,
										drop = FALSE]
			#Get conditions
			conditions <- sample_tbl$Sample.Conditions
			
			#Get groups of samples from sample conditions
			sample_grps<-unique(conditions)
			
			#Get median by cell type/line
			list_val_rows<-list()
			for (row in seq(nrow(DMRs))) {
				
				row_vec<-c()
				for (grp in sample_grps){
					
					#List sample names matching the group
					samples.in.grp<-rownames(sample_tbl[conditions == grp,])
					
					#List Methylation values of the group on the row
					meth_vals<-DMRs[row,samples.in.grp]
					
					# If more than 1 sample in the group get median
					if (length(samples.in.grp) > 1) {
						row_vec<-c(row_vec,median(meth_vals,na.rm = T))
					} else { #Else get value of sample
						row_vec<-c(row_vec,meth_vals)
					}
				}
				df_row_med<-data.frame(sample_grps,row_vec)
				#If no NA in group medians keep row for manipulation; Else remove it
				if (anyNA(row_vec) == F) {
					final_row <- c()
					#If no NA in full row do nothing
					if(anyNA(DMRs[row,]) == F) {
						# print(DMRs[row,])
						final_row <- c(DMRs[row,])
						#Add to list of valid rows
						list_val_rows[[row]]<-matrix(final_row, nrow = 1,
													 dimnames = list(rownames(DMRs)[row],
													 				colnames(DMRs)))
					} else { #Else replace NA value by median of its group if some
						
						for (smpl in names(DMRs[row,])) {
							
							#Get methylation value
							valmeth<-DMRs[row,
												  match(smpl,
												  	  rownames(sample_tbl))]
							
							#If value is a NA, replace by median of its group
							if (is.na(valmeth) == T) {
								final_row <- c(final_row,
											   df_row_med[df_row_med$sample_grps ==
											   		   	sample_tbl[rownames(sample_tbl) ==
											   		   			   	smpl,
											   		   			   ]$Sample.Conditions,
											   		   ]$row_vec)
							} else { #Else keep value as is
								final_row <- c(final_row, valmeth)
							}
						}
						#Add to list of valid rows
						list_val_rows[[row]]<-matrix(final_row, nrow = 1,
													 dimnames = list(rownames(DMRs)[row],
													 				colnames(DMRs)))
					}
				}
			}
			#Rbind all rows in the list into a matrix
			DMRs<-do.call("rbind",list_val_rows)
			
		}
	}  
}

##_Create Dendrogram_###########################################################

#Remove NAs if some for dendrogram matrix
dend_mat<- DMRs[complete.cases(DMRs),]
#Create Columns Dendrogram
ddgr <- as.dendrogram(hclust(parDist(t(dend_mat),method = selected_dist.cols,
                                     threads = 7)))

#Create Rows Dendrogram
row_dist<-parDist(dend_mat,method = selected_dist.rows,threads = 7)
row_hclust<-hclust(row_dist)
rm(row_dist)
rowclust<-as.dendrogram(row_hclust)

ddgr_dat<-dendro_data(ddgr) #Dendrogram data

#Function for creating dendograms
ggdend <- function(df,orientation) {
	ddplot<- ggplot() +
		geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
		theme_minimal() +
		theme(axis.title = element_blank(),
			  axis.text = element_blank(),
			  axis.ticks = element_blank(),
			  panel.grid = element_blank()) +
		expand_limits(x = c(0.4,max(df$x)+0.6), y = 0) +
		scale_x_continuous(expand = c(0, 0))
	
	if(orientation == "h"){
		ddplot <- ddplot +
			theme(plot.margin=unit(c(0.5,0,0,0),"cm")) +
			scale_y_continuous(expand = c(0, 0))
		
	} else {
		ddplot <- ddplot +
			theme(plot.margin=unit(c(0.02,0,0.02,0),"cm")) +
			scale_y_reverse(expand = c(0,0)) +
			coord_flip()
	}
	ddplot
}

if (selected_col_order == "distance method"){
	#Get dendrogram segments for columns
	if(display_dend_col == T){
		ddgr_seg_col <- ggdend(ddgr_dat$segments, orientation = "h") 
	}
}

#Get dendrogram segments for rows
if(display_dend_row == T){
	if(any(is.na(DMRs)) == F) {
		ddgr_seg_row <- ggdend(dendro_data(rowclust)$segments, orientation = "v") 
	} else {
		stop("Cannot create dendrogram of a matrix with rows containing NAs! Set display_dend_row to FALSE.")
	}
}

##_Create Heatmap with Ggplot2_#################################################

#Order of samples for the heatmap
if (selected_col_order == "distance method"){ #Order samples by distance method
	column.order<-order.dendrogram(ddgr)
	
} else if(selected_col_order == "cell type"){ #Order samples by cell types
	column.order<- sample_tbl[order(sample_tbl$Cell.types,
									sample_tbl$Short.genotypes),"Order"]
	
} else { #Order samples by genotypes
	column.order<- sample_tbl[order(sample_tbl$Genotypes),"Order"]
}

if(selected_row_order == "distance method") {
  row.order<-order.dendrogram(rowclust) #Get groups of close rows
}

if (selected_row_order == "variability across samples") {
  dframe <- DMRs[order(apply(DMRs, 1, sd, na.rm = T),
                               decreasing = TRUE),
                         column.order, drop = FALSE]
  if(htmp_type == "Rank") { #Keep Top N DMRs with N given by "dmr_rank"
    melted_mat <- melt(head(dframe,n=dmr_rank)) #Melt Matrix into a dataframe
  } else{ #Overview Heatmap Layout; Keep all DMRs
    melted_mat <- melt(dframe) #Melt Matrix into a dataframe
  }
  melted_mat$Var1<-factor(melted_mat$Var1,
                          levels(melted_mat$Var1)[length(levels(melted_mat$Var1)):1])
  colnames(melted_mat)[3]<-"Methylation"

} else if (selected_row_order == "distance method") {
  if (any(is.na(DMRs)) == F) {
    dframe<-DMRs[row.order,column.order, drop = FALSE]
    melted_mat <- melt(dframe) #Melt Matrix into a dataframe
    colnames(melted_mat)[3]<-"Methylation"
  } else {
    stop("Cannot apply distance method on rows containing NAs!")
  }
  
} else {
  dframe<-DMRs[,column.order, drop = FALSE]
  melted_mat <- melt(dframe) #Melt Matrix into a dataframe
  colnames(melted_mat)[3]<-"Methylation"
}

#Check if any(is.na()) is True or False and adapt ggplot accordingly
if (any(is.na(DMRs)) == F) { #No NA
  #Plot Heatmap
  htmp <- ggplot(data = melted_mat, aes(x = Var2, y = Var1, fill = Methylation)) +
    scale_fill_gradientn(colours = palette.heatmap)
  
} else { #NA present
  htmp <- ggplot(data = melted_mat,
                 aes(x = Var2, y = Var1, fill = Methylation, colour="")) +
    scale_fill_gradientn(colours = palette.heatmap, na.value = "black") +
    scale_colour_manual(values=NA) +
    guides(colour=guide_legend("N.A.", override.aes=list(colour="black")))
}

if (htmp_type == "Rank") {
  htmp <- htmp +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm"))
  
  Y_labels<-ggplot(data = melted_mat, aes(x = Var2, y = Var1, fill = Methylation)) +
    scale_x_discrete(expand = c(0,0),limits=c(0,0)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          plot.margin=unit(c(0,-15,-0.23,0.2),"cm"),
          panel.grid = element_blank(),
          panel.border = element_blank()) +
    ylab("DMR Locations")
  
  X_labels<-ggplot(data = melted_mat, aes(x = Var2, y = Var1, fill = Methylation)) +
    scale_y_discrete(expand = c(0,0),limits=c(0,0)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 11, angle = 90, hjust = 1,
                                     vjust = 0.5, face = "bold"),
          plot.margin=unit(c(0,0,0.2,-0.07),"cm"),
          panel.grid = element_blank(),
          panel.border = element_blank()) +
    xlab("Samples")
  
  
} else{ #For Distance Method or None
  if(display_dend_row == F){
    htmp <- htmp +
      theme(axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 11, angle = 90, hjust = 1,
                                       vjust = 0.5, face = "bold"),
            plot.margin=unit(c(0,0,0.2,0),"cm"))  
  } else {
    htmp <- htmp +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm"))
    
    X_labels<-ggplot(data = melted_mat, aes(x = Var2, y = Var1, fill = Methylation)) +
      scale_y_discrete(expand = c(0,0),limits=c(0,0)) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 11, angle = 90, hjust = 1,
                                       vjust = 0.5, face = "bold"),
            plot.margin=unit(c(0,0,0.2,-0.07),"cm"),
            panel.grid = element_blank(),
            panel.border = element_blank()) +
      xlab("Samples")
  }
}

htmp <- htmp +
  geom_tile() +
  theme(legend.text=element_text(size= 11),
        legend.title = element_text(size = 12)) +
  xlab("Samples")

##_Create Color Side Bar_#######################################################
Samples<-colnames(dframe)
#Create Genotypes column
Genotypes<-sample_tbl$Short.genotypes[column.order]
Genotypes<-factor(Genotypes,levels = color_table$Genotypes)
Genotype_table<-data.frame("Samples" = factor(Samples,levels = Samples),
                           "Genotypes" = Genotypes,row.names = NULL)

#Plot Color Sidebar
col_sidebar<-ggplot(Genotype_table, aes(Samples,"Genotypes")) +
  theme_minimal() +
  geom_tile(aes(fill = Genotypes)) +
  ylab(label = "Genotypes") +
  theme(legend.position = "right",
        legend.text=element_text(size= 11),
        legend.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm")) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values=color_table$Colors)

##_Extract Legend_##############################################################
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

htmp_legend<-g_legend(htmp)
sidebar_legend<-g_legend(col_sidebar)

##_Combine Dendrogram with Color Sidebar and Heatmap_###########################
#TODO check for margins in grid.arrange
if (htmp_type == "Rank") {
  grid.arrange(top = textGrob(paste0("Top",dmr_rank," DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs"),
                              gp = gpar(fontsize = 15, font = 1)),
               textGrob(""),
               textGrob(paste0("(", capitalize(selected_dist.cols), " distance; ",
                               "Rows ordered by ",selected_row_order,")"),
                        gp = gpar(fontsize=12, fontface=3L), hjust = 0.5),
               textGrob(""),
               textGrob(""),
               ddgr_seg_col,
               textGrob(""),
               textGrob(""),
               col_sidebar + theme(legend.position="none"),
               textGrob(""),
               Y_labels,
               htmp + theme(legend.position="none"),
               arrangeGrob(textGrob(""),
                           sidebar_legend, htmp_legend,
                           textGrob(""),ncol = 1),
               textGrob(""),
               X_labels,
               textGrob(""),
               ncol = 3, widths = c(2.5,8,2), heights = c(1,5,1,30,8))
  
} else {
  if(display_dend_col == T & display_dend_row == F){ #Display Columns Dendrogram only
    grid.arrange(top = textGrob("DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                                gp = gpar(fontsize = 15, font = 1)),
                 textGrob(paste0("(", capitalize(selected_dist.cols), " distance; ",
                                 "Rows ordered by ",selected_row_order,"; ",
                                 nrow(DMRs)," DMRs)"),
                          gp = gpar(fontsize=12, fontface=3L), hjust = 0.43),
                 textGrob(""),
                 arrangeGrob(ddgr_seg_col,
                             col_sidebar + theme(legend.position="none"),
                             htmp  + theme(legend.position="none"),
                             ncol = 1, heights = c(4,1,32)),
                 arrangeGrob(textGrob(""),
                             sidebar_legend, htmp_legend,
                             textGrob(""), textGrob(""), ncol = 1),
                 ncol = 2, widths = c(15,2), heights = c(1,30))
    
  } else if(display_dend_col == T & display_dend_row == T){ #Display Both Dendrogram
    grid.arrange(top = textGrob("DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                                gp = gpar(fontsize = 15, font = 1)),
                 textGrob(""),
                 textGrob(paste0("(", capitalize(selected_dist.cols), " distance; ",
                                 "Rows ordered by ",selected_row_order,"; ",
                                 nrow(DMRs)," DMRs)"),
                          gp = gpar(fontsize=12, fontface=3L), hjust = 0.43),
                 textGrob(""),
                 textGrob(""),
                 ddgr_seg_col,
                 textGrob(""),
                 textGrob(""),
                 col_sidebar + theme(legend.position="none"),
                 textGrob(""),
                 ddgr_seg_row,
                 htmp  + theme(legend.position="none"),
                 arrangeGrob(textGrob(""),
                             sidebar_legend, htmp_legend,
                             textGrob(""), textGrob(""), ncol = 1),
                 textGrob(""),
                 X_labels,
                 textGrob(""),
                 ncol = 3, widths = c(1,6,1), heights = c(1,5,1,30,8))
    
  } else if(display_dend_col == F & display_dend_row == T){ #Display Row Dendrogram only
    grid.arrange(top = textGrob("DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                                gp = gpar(fontsize = 15, font = 1)),
                 textGrob(""),
                 textGrob(paste0("(", capitalize(selected_dist.cols), " distance; ",
                                 "Rows ordered by ",selected_row_order,"; ",
                                 nrow(DMRs)," DMRs)"),
                          gp = gpar(fontsize=12, fontface=3L), hjust = 0.43),
                 textGrob(""),
                 textGrob(""),
                 col_sidebar + theme(legend.position="none"),
                 textGrob(""),
                 ddgr_seg_row,
                 htmp  + theme(legend.position="none"),
                 arrangeGrob(textGrob(""),
                             sidebar_legend, htmp_legend,
                             textGrob(""), textGrob(""), ncol = 1),
                 textGrob(""),
                 X_labels,
                 textGrob(""),
                 ncol = 3, widths = c(1,6,1), heights = c(1,1,30,8))
    
  } else { #Hide Dendrograms
    grid.arrange(top = textGrob("DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                                gp = gpar(fontsize = 15, font = 1)),
                 textGrob(paste0("(", capitalize(selected_dist.cols), " distance; ",
                                 "Rows ordered by ",selected_row_order,"; ",
                                 nrow(DMRs)," DMRs)"),
                          gp = gpar(fontsize=12, fontface=3L), hjust = 0.43),
                 textGrob(""),
                 arrangeGrob(col_sidebar + theme(legend.position="none"),
                             htmp  + theme(legend.position="none"),
                             ncol = 1, heights = c(1,36)),
                 arrangeGrob(textGrob(""),
                             sidebar_legend, htmp_legend,
                             textGrob(""), textGrob(""), ncol = 1),
                 ncol = 2, widths = c(15,2), heights = c(1,30))
  }
}
