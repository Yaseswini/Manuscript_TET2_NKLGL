## +++++++++++++++++++++++++++++++++ ## 
## 	 Notes : hg19 reference			 ## 
# main question : 
# consequences of TET2 mutations
# TET2 was done on PBMC dan 
# lowVAF is 18% 
# proposed mechanisms : altered methylation can activate STAT3 signaling
# 
# TET2 mutant
# 1444       9/26/2012          purified NK
# 1511      10/29/2014        purified CD94+
# 1820       7/16/2013          purified CD94+


# TET2 wild-type
# 1856      8/5/2014             PBMC (78% leukemic burden)
# 1866       9/30/2014           purified CD94+
# 1791      4/16/2013           PBMC (93% leukemic burden)
#
## +++++++++++++++++++++++++++++++++ ## 
.libPaths( c( .libPaths() , "/home/yn9w/R/x86_64-pc-linux-gnu-library/3.6" ) )
options( stringsAsFactors = F )
options("scipen"=100, "digits"=4)
set.seed(47) 

library(tidyverse)
library(readr)
library(data.table)
library(ggfortify)
library(pheatmap)
library(methylKit)
library(ggrepel)

projectDIR = "/Users/yaseswini/Documents/GitHub/Manuscript_TET2_NKLGL"

source( file.path( projectDIR , "scripts/TET2_NKLGL_HelperFunctions.R" ) )

gffRDA = file.path( projectDIR , "data" , "gencode_v19_2020-08-12.rds" )
gffData = readRDS( gffRDA )

rrbsData = <rrbs_data_dir>

rrbsFiles = list.files( rrbsData , full.names = T , pattern = "*.bed" )
rrbsMeta = read_tsv( file.path( rrbsData , "metadata.txt" ) ) %>% 
				mutate( SampleID = gsub( "^X", "" , Sample_ID ) , 
							BedFilePath = file.path( rrbsData , gsub(".gz","",Sample_Filename) ) , 
							CpGFilePath = file.path( rrbsData , "FilesForMethylKit" , gsub(".bed.gz",".txt",Sample_Filename) ) ) %>% data.frame

errbsYN_analysisDIR = "<output_directory>"


## ========  Reformatting files for differential methylation & unsupervised analysis ====== ## 
rrbsMeta = read_tsv( file.path( rrbsData , "metadata.txt" ) ) %>% 
				mutate( SampleID = gsub( "^X", "" , Sample_ID ) , 
						BedFilePath = file.path( rrbsData , gsub(".gz","",Sample_Filename) ) ) %>% data.frame
fileList_new = list()
for( i in 1:nrow(rrbsMeta) )
{	
	dat = read_tsv( rrbsMeta[i,"BedFilePath"] , col_names = c("chr","start","end","freqC","numC","numT") )
	# required column : chrBase , chr , base , strand , coverage , freqC , freqT 
	dat1 = dat %>% 
				mutate( chrBase = paste0(chr,".",start) ,
						strand = "+" , 
						base = start , 
						coverage = numC + numT ,
						freqT = (numT/(numC+numT))*100 ) %>% 
				dplyr::select( chrBase , chr , base , strand , coverage , freqC, freqT )

	write.table( dat1 , file = file.path( rrbsData , "FilesForMethylKit" , paste0( rrbsMeta[i,"SampleID"] , ".txt" ) ) , sep = "\t" , row.names = F , col.names = T , quote = F )
	fileList_new[[ rrbsMeta[i,"SampleID"] ]] = file.path( rrbsData , "FilesForMethylKit" , paste0( rrbsMeta[i,"SampleID"] , ".txt" ) )

	## filter for coverage : 
	dat_filtered = dat1 %>% filter( coverage >= 10 )
	write.table( dat_filtered , file = file.path( rrbsData , "FilesForMethylKit" , paste0( rrbsMeta[i,"SampleID"] , ".mincov10.txt" ) ) , sep = "\t" , row.names = F , col.names = T , quote = F )
}


## ==== Generate a table of methylation values @ 0x coverage and @ 10x coverage 
rrbsFiles = list.files( file.path( <rrbsData> , "FilesForMethylKit" ) , full.names = T , pattern = "*.txt" )
names(rrbsFiles) = basename(rrbsFiles)
rrbsFiles_0x = rrbsFiles[ ! grepl("mincov10",names(rrbsFiles)) ]
rrbsFiles_10x = rrbsFiles[ grepl("mincov10",names(rrbsFiles)) ]
names(rrbsFiles_10x) = gsub(".mincov10","",names(rrbsFiles_10x))
rrbsList_0x = list() ; rrbsList_10x = list()

for( fn in rrbsFiles_0x )
{
	dt = fread( fn )
	dt = dt[, c("chrBase","freqC"), with=FALSE]
	setnames( dt , "freqC", gsub(".txt","",basename(fn)) )
	setkey( dt , chrBase )
	rrbsList_0x[[ basename(fn) ]] = dt
}
methTab_0x = Reduce( function(...) merge(...,all=T) , rrbsList_0x )
methTab_0x = methTab_0x %>% as_tibble() %>% column_to_rownames( var = "chrBase" )
methTab_0x_noNAs = na.omit( methTab_0x  )

for( fn in rrbsFiles_10x )
{
	dt = fread( fn )
	dt = dt[, c("chrBase","freqC"), with=FALSE]
	setnames( dt , "freqC", gsub(".txt","",basename(fn)) )
	setkey( dt , chrBase )
	rrbsList_10x[[ basename(fn) ]] = dt
}
methTab_10x = Reduce( function(...) merge(...,all=T) , rrbsList_10x )
methTab_10x = methTab_10x %>% as_tibble() %>% column_to_rownames( var = "chrBase" )
methTab_10x_noNAs = na.omit( methTab_10x  )
colnames(methTab_10x) = gsub(".mincov10","",colnames(methTab_10x)) 
colnames(methTab_10x_noNAs) = gsub(".mincov10","",colnames(methTab_10x_noNAs)) 

write.table( methTab_0x , file = file.path( errbsYN_analysisDIR , "methTab_0x.txt" ) , sep = "\t" , quote = F , row.names = T , col.names = NA )
write.table( methTab_0x_noNAs , file = file.path( errbsYN_analysisDIR , "methTab_0x_naOmit.txt" ) , sep = "\t" , quote = F , row.names = T , col.names = NA )
write.table( methTab_10x , file = file.path( errbsYN_analysisDIR , "methTab_10x.txt" ) , sep = "\t" , quote = F , row.names = T , col.names = NA )
write.table( methTab_10x_noNAs , file = file.path( errbsYN_analysisDIR , "methTab_10x_naOmit.txt" ) , sep = "\t" , quote = F , row.names = T , col.names = NA )

## ==== Reading the above files === ##
methTab_0x = read_tsv( file.path( errbsYN_analysisDIR , "methTab_0x.txt" ) ) %>% column_to_rownames( var = "X1" )
methTab_10x = read_tsv( file.path( errbsYN_analysisDIR , "methTab_10x.txt" ) ) %>% column_to_rownames( var = "X1" )
methTab_0x_noNAs = read_tsv( file.path( errbsYN_analysisDIR , "methTab_0x_naOmit.txt" ) ) %>% column_to_rownames( var = "X1" )
methTab_10x_noNAs = read_tsv( file.path( errbsYN_analysisDIR , "methTab_10x_naOmit.txt" ) ) %>% column_to_rownames( var = "X1" )

### == Annotating cpg sites === ##
allCpGs = methTab_0x %>% rownames_to_column( var = "chrBase" ) %>% mutate( chrBase1 = chrBase ) %>% separate( chrBase1 , c("chr","start") ,sep = "[.]")	
allCpGs$end = allCpGs$start
allCpGs.AnnotationPromoter = annotate_sites( loci = allCpGs[ , c("chr","start","end") ] , 
									  		 referenceAnnot = gffData ,
                				  	  		 preference = c("Promoter") ,
                				 	 		 overlap.type = "within" )
allCpGs.AnnotationPromoter = merge_annotationRes( allCpGs.AnnotationPromoter ) %>% dplyr::rename( "chrBase_3coord" = "chrBase" ) %>% mutate( chrBase = gsub( "^(chr.*)[.].*", "\\1" , chrBase_3coord ) ) 
allCpGs.AnnotationPromoterTab = allCpGs %>% left_join( allCpGs.AnnotationPromoter ) %>% filter( Region == "Promoter" )


### ==== Identifying highvariance genes ==== ## 
quantileThresh = 0.90
cpg_iqr_0x = apply( methTab_0x_noNAs , 1 , IQR )
iqr_0x_quantiles = quantile( cpg_iqr_0x , quantileThresh )
highvarBasedonIQR_cpgs_0x = cpg_iqr_0x[ cpg_iqr_0x > iqr_0x_quantiles ]
highvarBasedonIQR.cpgs_0xTab = methTab_0x_noNAs[ names(highvarBasedonIQR_cpgs_0x) , ]

cpgs_10x = rownames(methTab_10x_noNAs)
cpgs_10x_iqr = cpg_iqr_0x[ names(cpg_iqr_0x) %in% cpgs_10x ]
iqr_10x_quantiles = quantile( cpgs_10x_iqr , quantileThresh )
highvarBasedonIQR_cpgs_10x = cpgs_10x_iqr[ cpgs_10x_iqr > iqr_10x_quantiles ]
highvarBasedonIQR_cpgs_10xTab = methTab_10x_noNAs[ names(highvarBasedonIQR_cpgs_10x) , ]

# cpg_iqr_10x = apply( methTab_10x_noNAs , 1 , IQR )
# iqr_10x_quantiles = quantile( cpg_iqr_10x , quantileThresh )
# highvarBasedonIQR_cpgs_10x = cpg_iqr_10x[ cpg_iqr_10x > iqr_10x_quantiles ]
# highvarBasedonIQR_cpgs_10xTab = methTab_10x_noNAs[ names(highvarBasedonIQR_cpgs_10x) , ]

write.table( highvarBasedonIQR.cpgs_0xTab , file = file.path( errbsYN_analysisDIR , "methTab_0x_highvariance_90IQR.txt" ) , row.names = T , col.names = NA , quote = F , sep = "\t" )
write.table( highvarBasedonIQR_cpgs_10xTab , file = file.path( errbsYN_analysisDIR , "methTab_10x_highvariance_90IQR.txt" ) , row.names = T , col.names = NA , quote = F , sep = "\t" ) 

## =================================================================== ## 
## 			Map CpGs in the promoter regions of genes of interest	   ## 
## =================================================================== ## 
genes = c( "TET2",paste0("PIAS",c(1:4)),
		    paste0("SOCS",c(1:7)), 
		    paste0("SHP",c(1:2)) , 
		    paste0("PTPN",LETTERS), 
		    paste0("PTPR",LETTERS), 
		    paste("CISH") , 
		    paste0("PTPN",c(1:30)),"PTPN20CP" )

dir.create( file.path( errbsYN_analysisDIR , "MethylationPlotter" ) )
dir.create( file.path( errbsYN_analysisDIR , "Heatmaps" ) )

##  PTPRN and PTPRD ## 
for( g1 in genes )
{
	g1_CpGs_promoter = allCpGs.AnnotationPromoter %>% dplyr::filter( gene_name == g1 ) %>% .$chrBase %>% unique

	g1_promoterCpGs_0x = methTab_0x[ g1_CpGs_promoter , ]
	g1_promoterCpGs_0x_noNAs = methTab_0x_noNAs[ g1_CpGs_promoter , ]

	g1_promoterCpGs_10x = methTab_10x[ which( rownames(methTab_10x) %in% g1_CpGs_promoter ) , ]
	g1_promoterCpGs_10x_noNAs = methTab_10x_noNAs[ which( rownames(methTab_10x_noNAs) %in% g1_CpGs_promoter ) , ]

	## File for methylation plotter 
	datForMethylationPlotter_10x = t( g1_promoterCpGs_10x ) %>% data.frame %>% rownames_to_column( var = "SampleID" ) %>% left_join( rrbsMeta %>% dplyr::select( SampleID , TET2Status1 ) ) %>% arrange( TET2Status1 )
	datForMethylationPlotter_0x = t( g1_promoterCpGs_0x ) %>% data.frame %>% rownames_to_column( var = "SampleID" ) %>% left_join( rrbsMeta %>% dplyr::select( SampleID , TET2Status1 ) )  %>% arrange( TET2Status1 )

	write.table( datForMethylationPlotter_10x , file = file.path( errbsYN_analysisDIR , "MethylationPlotter" , paste0( g1 , "_allCpGs_at_10xCoverage.txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T )
	write.table( datForMethylationPlotter_0x , file = file.path( errbsYN_analysisDIR , "MethylationPlotter" , paste0( g1 , "_allCpGs_at_0xCoverage.txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T  )

	## Heatmaps for both 10x and 0x 
	pdf(  file.path( errbsYN_analysisDIR , "Heatmaps" , paste0( g1 , "_allCpGs_at_0x_and_10x.pdf" ) ) )
	heatmapAnnotation.df = rrbsMeta %>% 
							dplyr::select( SampleID , TET2Status , TET2Status1 ) %>% 
							mutate( TET2Status1 = factor(TET2Status1,levels=c("Normal","WT","Mut")) ) %>%
							arrange( TET2Status1 ) %>%
							column_to_rownames( var = "SampleID" ) %>% 
							data.frame

	heatmapAnnotation.colors = list( TET2Status1 = c("Mut"="#440154FF","WT"="#21908CFF","Normal"="#5DC963FF") , 
								 	 TET2Status = c("Q1523*" = "#F7F7F7" , "X1394_splice" = "#B2ABD2" , "R1452*" = "#5E3C99" , "Normal" = "#5DC963FF" , "WT" = "#21908CFF") )

	heatmapAnnotation.obj = HeatmapAnnotation( df = heatmapAnnotation.df  , 
								 		   col = heatmapAnnotation.colors , 
								  		   gap = unit(3, "points") , 
								 		   border = TRUE , 
								  		   show_legend = TRUE , 
								  		   simple_anno_size_adjust = TRUE , 
								  		   simple_anno_size = unit(0.3, "cm") )

	if( nrow(g1_promoterCpGs_0x) > 0 )
	{
		g1_promoterCpGs_0x_htData = g1_promoterCpGs_0x[ , rownames(heatmapAnnotation.df) ]
		heatmap.breaks = seq(0,100,by=5)
		heatmap.cols = colorRampPalette( c("blue" ,"white","orange") )( length(heatmap.breaks)-1 )
		#heatmap.colFun = colorRamp2( heatmap.breaks , rnaseq.cols )

		ht = Heatmap( as.matrix( g1_promoterCpGs_0x_htData ) , 
			  col = heatmap.cols , 
			  name = "% CpG methylation", 
			  cluster_rows = FALSE , cluster_columns = FALSE , 
			  clustering_distance_rows = "pearson" , clustering_distance_columns = "pearson" , 
			  clustering_method_rows = "ward" , clustering_method_columns = "ward" , 
			  show_row_names  = FALSE , use_raster = FALSE , top_annotation = heatmapAnnotation.obj , 
			  column_names_gp = gpar(fontsize = 3) , 
			  column_title = paste0( "All CpGs ( 0x ) in " , g1 , " promoters" , "\n" , "Number of CpGs : " , nrow(g1_promoterCpGs_0x) ) )
		draw( ht )
		cat( g1 , "\n" )
	}
	if( nrow(g1_promoterCpGs_10x) > 0 )
	{
		g1_promoterCpGs_10x_htData = g1_promoterCpGs_10x[ , rownames(heatmapAnnotation.df) ]
		heatmap.breaks = seq(0,100,by=5)
		heatmap.cols = colorRampPalette( c("blue" ,"white","orange") )( length(heatmap.breaks)-1 )
		#heatmap.colFun = colorRamp2( heatmap.breaks , rnaseq.cols )

		ht = Heatmap( as.matrix(g1_promoterCpGs_10x_htData) , 
			  col = heatmap.cols , 
			  name = "% CpG methylation", 
			  cluster_rows = FALSE , cluster_columns = FALSE , 
			  clustering_distance_rows = "pearson" , clustering_distance_columns = "pearson" , 
			  clustering_method_rows = "ward" , clustering_method_columns = "ward" , 
			  show_row_names  = FALSE , use_raster = FALSE , top_annotation = heatmapAnnotation.obj , 
			  column_names_gp = gpar(fontsize = 3) , 
			  column_title = paste0( "All CpGs ( 10x ) in " , g1 , " promoters\nNumber of CpGs : " , nrow(g1_promoterCpGs_10x_htData) ) )
		draw( ht )
		cat( g1 , "\n" )
	}
	dev.off()
}


## ========================================================== ## 
## 												PCA 							  								##
## ========================================================== ## 
# rrbsMeta = read_tsv( file.path( rrbsData , "metadata.txt" ) ) %>% 
# 				mutate( SampleID = gsub( "^X", "" , Sample_ID ) , BedFilePath = file.path( rrbsData , Sample_Filename ) , 
# 					    TET2Status = if_else( is.na(TET2Status) , "Normal" , TET2Status ) , 
# 					    TET2Status1 = if_else( TET2Status != "Normal" & TET2Status != "WT" , "Mut" , TET2Status ), 
# 					    pcaLabel = paste0( SampleID , "( " , TET2Status , " )" ) ) %>% 
# 				mutate( TET2Status1 = factor( TET2Status1 , levels = c("Normal","WT","Mut")) )

rrbsMeta = read_tsv( file.path( rrbsData , "metadata.txt" ) ) %>% 
 				mutate( SampleID = gsub( "^X", "" , Sample_ID ) , BedFilePath = file.path( rrbsData , Sample_Filename ) , 
 					    TET2Status = if_else( is.na(TET2Status) , "Normal" , TET2Status ) , 
 					    TET2Status1 = if_else( TET2Status != "Normal" & TET2Status != "WT" , "Mut" , TET2Status ), 
 					    pcaGroup = if_else( TET2Status != "Normal" & TET2Status != "WT" , "TET2 mutant" , TET2Status ) ) %>% 
 				mutate( TET2Status1 = factor( TET2Status1 , levels = c("Normal","WT","Mut")) )


sampleCol = c("TET2 mutant"="#440154FF","WT"="#21908CFF","Normal"="#5DC963FF")

# highvarBasedonIQR.cpgs_0xTab = read_tsv( file.path( errbsYN_analysisDIR , "methTab_0x_highvariance_90IQR.txt" ) ) %>% column_to_rownames( var = "X1" )
# highvarBasedonIQR_cpgs_10xTab = read_tsv( file.path( errbsYN_analysisDIR , "methTab_10x_highvariance_90IQR.txt" ) ) %>% column_to_rownames( var = "X1" )

### === 0x coverage 
pcaBasedonIQR = prcomp( t(  highvarBasedonIQR.cpgs_0xTab  ) , center = T )
pcaBasedonIQR_scores = pcaBasedonIQR$x
pcaBasedonIQR_variance = unclass( summary(pcaBasedonIQR) )$importance["Proportion of Variance",]
datPCA = data.frame( SampleID = rownames(pcaBasedonIQR_scores), PC1=pcaBasedonIQR_scores[,1] , PC2=pcaBasedonIQR_scores[,2] )
rownames(datPCA) = NULL 
xmax = ymax = ceiling( max( datPCA[ , c("PC1","PC2") ] ) )
xmin = ymin = floor( min( datPCA[ , c("PC1","PC2") ] ) )
if( abs(xmin) > xmax ) { xmax = -xmin } else if( abs(xmin) < xmax ) { xmin = -xmax }  # xmin is always negative so -xmin will make it positive
datPCA = datPCA %>% left_join( rrbsMeta )

plot.title = paste( paste0("Highvariance CpGs @ " , quantileThresh , " percentile of IQR") , "\nNumber of data points = " , ncol( t( highvarBasedonIQR.cpgs_0xTab ) )  , sep = "" )
xlab = paste("PC1" , "( Variance = " , (as.numeric(pcaBasedonIQR_variance["PC1"]) * 100 ) , "% ) " , sep = "" )
ylab = paste("PC2" , " ( Variance = " , (as.numeric(pcaBasedonIQR_variance["PC2"]) * 100 ) , "% ) "  , sep = "" )
  
p = ggplot( datPCA , aes(x=PC1,y=PC2,color=pcaGroup,label=SampleID) ) 
p = p + geom_point( stat='identity',size = 3 ) 
p = p + scale_color_manual(values = sampleCol) + theme_bw() + theme( legend.position="bottom" , legend.text = element_text(size = 7),legend.title=element_blank() , plot.title = element_text( size = 9.5 , face = "bold" ) , axis.text = element_text( size = 7 ) , axis.title = element_text( size = 9 ) ) + scale_x_continuous( limits = c(xmin , xmax) ) + scale_y_continuous( limits = c(xmin,xmax) ) 
p = p + labs( x = xlab , y = ylab , title=plot.title ) + geom_text( size = 2 , hjust = 1 , vjust = 1 )
ggsave( p , file = file.path( errbsYN_analysisDIR , "TET2_NKLGL_HighvarianceCpGs_at0xcoverage_PCA_15Feb2021.pdf" ) )


### === 10x coverage 
pcaBasedonIQR = prcomp( t(  highvarBasedonIQR_cpgs_10xTab  ) , center = T )
pcaBasedonIQR_scores = pcaBasedonIQR$x
pcaBasedonIQR_variance = unclass( summary(pcaBasedonIQR) )$importance["Proportion of Variance",]
datPCA = data.frame( SampleID = rownames(pcaBasedonIQR_scores), PC1=pcaBasedonIQR_scores[,1] , PC2=pcaBasedonIQR_scores[,2] )
rownames(datPCA) = NULL 
xmax = ymax = ceiling( max( datPCA[ , c("PC1","PC2") ] ) )
xmin = ymin = floor( min( datPCA[ , c("PC1","PC2") ] ) )
if( abs(xmin) > xmax ) { xmax = -xmin } else if( abs(xmin) < xmax ) { xmin = -xmax }  # xmin is always negative so -xmin will make it positive
datPCA = datPCA %>% left_join( rrbsMeta )

plot.title = paste( paste0("Highvariance CpGs @ " , quantileThresh , " percentile of IQR") , "\nNumber of data points = " , ncol( t( highvarBasedonIQR_cpgs_10xTab ) )  , sep = "" )
xlab = paste("PC1" , "( Variance = " , (as.numeric(pcaBasedonIQR_variance["PC1"]) * 100 ) , "% ) " , sep = "" )
ylab = paste("PC2" , " ( Variance = " , (as.numeric(pcaBasedonIQR_variance["PC2"]) * 100 ) , "% ) "  , sep = "" )
  
p = ggplot( datPCA , aes(x=PC1,y=PC2,color=pcaGroup,label=SampleID) ) 
p = p + geom_point( stat='identity',size = 3 ) 
p = p + scale_color_manual(values = sampleCol) + theme_bw() + theme( legend.position="bottom" , legend.text = element_text(size = 7),legend.title=element_blank() , plot.title = element_text( size = 9.5 , face = "bold" ) , axis.text = element_text( size = 7 ) , axis.title = element_text( size = 9 ) ) + scale_x_continuous( limits = c(xmin , xmax) ) + scale_y_continuous( limits = c(xmin,xmax) ) 
p = p + labs( x = xlab , y = ylab , title=plot.title ) + geom_text_repel( )
ggsave( p , file = file.path( errbsYN_analysisDIR , "TET2_NKLGL_HighvarianceCpGs_at10xcoverage_PCA_16Feb2021.pdf" ) )



 