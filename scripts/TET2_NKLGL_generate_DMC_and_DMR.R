## TET2_NKLGL Differentially Methylated Sites & Regions ## 

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


### === Calculating differentially methylated sites for tet2mut vs tet2wt === ### 
rrbsMeta = rrbsMeta %>% mutate( CpGFilePath = file.path( rrbsData , "FilesForMethylKit" , gsub(".bed.gz",".txt",Sample_Filename) ) ) %>% data.frame
methDiff_thresh = 20
qvalue_thresh = 0.01

### === tet2mut vs normals 
tet2Mut_vs_normals = as.list( c( rrbsMeta %>% filter( TET2Status1 == "Mut" ) %>% .$CpGFilePath , 
					   			 rrbsMeta %>% filter( TET2Status1 == "Normal" ) %>% .$CpGFilePath ) )
names(tet2Mut_vs_normals) = c( rrbsMeta %>% filter( TET2Status1 == "Mut" ) %>% .$SampleID , 
							   rrbsMeta %>% filter( TET2Status1 == "Normal" ) %>% .$SampleID )

outLabel = "TET2Mut_vs_Normals"
fileList = as.list(tet2Mut_vs_normals)
trt = c( 1,1,1,0,0,0,0,0 )

tet2Mut_vs_normals_methylkitRes  = run_methylkit( file.list = tet2Mut_vs_normals , 
			   		   		  					  sample.id = as.list( names(tet2Mut_vs_normals) ) , 
			   		  		  					  assembly = "hg19" , 
			  		    	  					  treatment = trt , 
			   		   		  					 outLabel = file.path( errbsYN_analysisDIR , "TET2Mut_vs_normals"  ) )

tet2mut_vs_normals_DMCs = tet2Mut_vs_normals_methylkitRes %>% 
										data.frame %>% 
										filter( qvalue < qvalue_thresh & abs(meth.diff) > methDiff_thresh ) %>% 
										mutate(  chrBase = paste( chr , start , end , sep = "." ) ) %>% 
										dplyr::select( chrBase , meth.diff ) %>% 
										dplyr::rename( "TET2mut_vs_Normals" = "meth.diff" )

### ==== Calculating differentially methylated sites for Tet2 wt vs normals ==== ### 
tet2WT_vs_normals = as.list( c( rrbsMeta %>% filter( TET2Status1 == "WT" ) %>% .$CpGFilePath , 
					   			rrbsMeta %>% filter( TET2Status1 == "Normal" ) %>% .$CpGFilePath ) )
names(tet2WT_vs_normals) = c( rrbsMeta %>% filter( TET2Status1 == "WT" ) %>% .$SampleID , 
							  rrbsMeta %>% filter( TET2Status1 == "Normal" ) %>% .$SampleID )

outLabel = "TET2WT_vs_Normals"
fileList = as.list(tet2WT_vs_normals)
trt = c( 1,1,1,0,0,0,0,0 )

tet2WT_vs_normals_methylkitRes = run_methylkit( file.list = tet2WT_vs_normals , 
			   		   		  sample.id = as.list( names(tet2WT_vs_normals) ) , 
			   		  		  assembly = "hg19" , 
			  		    	  treatment = trt , 
			   		   		  outLabel = file.path( errbsYN_analysisDIR  , "TET2WT_vs_Normals"  ) )

tet2WT_vs_normals_DMCs = tet2WT_vs_normals_methylkitRes %>% 
										data.frame %>% 
										filter( qvalue < qvalue_thresh & abs(meth.diff) > methDiff_thresh ) %>% 
										mutate(  chrBase = paste( chr , start , end , sep = "." ) ) %>% 
										dplyr::select( chrBase , meth.diff ) %>% 
										dplyr::rename( "TET2WT_vs_Normals" = "meth.diff" )


dmcTab = tet2mut_vs_normals_DMCs %>% left_join( tet2WT_vs_normals_DMCs )
dmc.barplot = dmc_barplot( dmcTab  , colOrder = c("TET2WT_vs_Normals","TET2mut_vs_Normals") )
dmc.barplot = dmc.barplot + coord_fixed( 30000 )

ggsave( dmc.barplot , file = file.path( errbsYN_analysisDIR , paste0( "TET2_NKLGL_DifferentiallyMethylatedSites_" , Sys.Date() , ".png" ) ) , width = 5 , height = 5 )

### ==== Annotating dmcs to genomic sites ==== ## 
dmcTab = dmcTab %>% mutate( chrBase1 = chrBase ) %>% separate( chrBase1 , c("chr","start","end") , sep = "[.]" )
dmcTab.annotation = annotate_sites( loci = dmcTab[ , c("chr","start","end") ] ,
	                				  referenceAnnot = gffData ,
	                				  preference = c("Promoter","3UTR","5UTR","Exon","Intron","Intergenic") ,
	                				  overlap.type = "within" )
dmcTab.annotation = merge_annotationRes( dmcTab.annotation )
dmcTab.annotation = dmcTab %>% left_join( dmcTab.annotation )

dmcTab.annotation_toWrite = dmcTab.annotation %>% 
						mutate( status = case_when( is.na(TET2WT_vs_Normals) & ! is.na(TET2mut_vs_Normals)  ~ "Unique TET2 Mut" ) ) %>% 
						dplyr::select( -chr , -start , -end , -tx_id , -tx_name , -transcript_type , -exon_id , -exon_id1 , -exon_rank) %>% 
						unique
dmcTab.annotationPromoter = dmcTab.annotation_toWrite %>% filter( Region == "Promoter" )

write.table( dmcTab.annotation , file = file.path( errbsYN_analysisDIR , paste0( "TET2_NKLGL_DifferentiallyMethylatedSites_" , Sys.Date() , ".txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T )
write.table( dmcTab.annotationPromoter , file = file.path( errbsYN_analysisDIR , paste0( "TET2_NKLGL_DifferentiallyMethylatedSites_In_Promoter_" , Sys.Date() , ".txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T )

## === Extracting unique CpGs in TET2mut vs NBM
TET2mut_vs_Normals_uniqueHyper = dmcTab.annotation %>% filter( is.na(TET2WT_vs_Normals) & ! is.na(TET2mut_vs_Normals) ) %>% filter( TET2mut_vs_Normals > 0 )
TET2mut_vs_Normals_uniqueHypo = dmcTab.annotation %>% filter( is.na(TET2WT_vs_Normals) & ! is.na(TET2mut_vs_Normals) ) %>% filter( TET2mut_vs_Normals < 0 )

TET2mut_vs_Normals_uniqueHyperBed = TET2mut_vs_Normals_uniqueHyper %>% 
										dplyr::select( chr , start , end ) %>% 
										unique %>% filter( ! grepl("^chrM",chr) )

TET2mut_vs_Normals_uniqueHypoBed = TET2mut_vs_Normals_uniqueHypo %>% 
										dplyr::select( chr , start , end ) %>% 
										unique %>% filter( ! grepl("^chrM",chr) )


write.table( TET2mut_vs_Normals_uniqueHyper , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHyper.txt" ) , sep = "\t" , quote = F , row.names = F , col.names = T )
write.table( TET2mut_vs_Normals_uniqueHyperBed , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHyper.bed" ) , sep = "\t" , quote = F , row.names = F , col.names = F )
write.table( TET2mut_vs_Normals_uniqueHypo , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHypo.txt" ) , sep = "\t" , quote = F , row.names = F , col.names = T )
write.table( TET2mut_vs_Normals_uniqueHypoBed , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHypo.bed" ) , sep = "\t" , quote = F , row.names = F , col.names = F )

### ==== background cpgs for great tool
tet2WT_vs_normals_methylkitRes = get(load(file.path( errbsYN_analysisDIR , "TET2WT_vs_Normals_myDiff_2021-02-04.rda" ))) %>% dplyr::select( chr , start , end )
tet2Mut_vs_normals_methylkitRes = get(load(file.path( errbsYN_analysisDIR , "TET2Mut_vs_normals_myDiff_2021-02-04.rda" )))  %>% dplyr::select( chr , start , end )

background_cpgs = rbind( tet2WT_vs_normals_methylkitRes , tet2Mut_vs_normals_methylkitRes ) %>% unique
write.table( background_cpgs , file = file.path( errbsYN_analysisDIR , "backgroundCpGs_for_GREAT.bed" ) , quote = F , sep = "\t" , row.names = F , col.names = F )

## ============== Edmr ======= ## 
run_edmr( tet2WT_vs_normals_methylkitRes , outfile.label = file.path( errbsYN_analysisDIR , "TET2WT_vs_Normals" ) )
run_edmr( tet2Mut_vs_normals_methylkitRes , outfile.label = file.path( errbsYN_analysisDIR , "TET2Mut_vs_Normals" ) )

tet2WT_vs_normals_dmr = get( load( file.path( errbsYN_analysisDIR , "TET2WT_vs_Normals_myeDMR_2021-02-05.rda" ) ) )  %>% data.frame %>% 
								filter( DMR.qvalue < 0.01 & abs(mean.meth.diff) > 25 ) %>% 
								mutate( "TET2WT_vs_Normals" = mean.meth.diff , "chrBase" = paste( seqnames , start , end , sep = ".") ) %>% 
								dplyr::select( chrBase , TET2WT_vs_Normals )

tet2Mut_vs_normals_dmr = get( load( file.path( errbsYN_analysisDIR , "TET2Mut_vs_Normals_myeDMR_2021-02-05.rda" ) ) ) %>% data.frame %>% 
								filter( DMR.qvalue < 0.01 & abs(mean.meth.diff) > 25 )  %>% 
								mutate( "TET2mut_vs_Normals" = mean.meth.diff , "chrBase" = paste( seqnames , start , end , sep = ".") ) %>% 
								dplyr::select( chrBase , TET2mut_vs_Normals )

dmrTab = tet2Mut_vs_normals_dmr %>% left_join( tet2WT_vs_normals_dmr )
dmr.barplot = dmc_barplot( dmrTab  , colOrder = c("TET2WT_vs_Normals","TET2mut_vs_Normals") )
dmr.barplot = dmr.barplot + coord_fixed( 2000 )
ggsave( dmr.barplot , file = file.path( errbsYN_analysisDIR , paste0( "TET2_NKLGL_DifferentiallyMethylatedRegions_" , Sys.Date() , ".pdf" ) ) , width = 5 , height = 5 )

dmrTab = dmrTab %>% mutate( chrBase1 = chrBase ) %>% separate( chrBase1 , c("chr","start","end") , sep = "[.]" )
dmrTab.annotation = annotate_sites( loci = dmrTab[ , c("chr","start","end") ] ,
                				  referenceAnnot = gffData ,
                				  preference = c("Promoter","3UTR","5UTR","Exon","Intron","Intergenic") ,
                				  overlap.type = "within" )
dmrTab.annotation = merge_annotationRes( dmrTab.annotation )
dmrTab.annotation = dmrTab %>% left_join( dmrTab.annotation )

dmrTab.annotation_toWrite = dmrTab.annotation %>% 
						mutate( status = case_when( is.na(TET2WT_vs_Normals) & ! is.na(TET2mut_vs_Normals)  ~ "Unique TET2 Mut" ) ) %>% 
						dplyr::select( -chr , -start , -end , -tx_id , -tx_name , -transcript_type , -exon_id , -exon_id1 , -exon_rank) %>% 
						unique
dmrTab.annotationPromoter = dmrTab.annotation_toWrite %>% filter( Region == "Promoter" )
write.table( dmrTab.annotation_toWrite , file = file.path( errbsYN_analysisDIR , paste0("TET2_NKLGL_DifferentiallyMethylatedRegions_" , Sys.Date() , ".txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T )

hyperDMR_uniqueTET2mut = dmrTab.annotation %>% filter( ! grepl("^chrM",chr) ) %>% filter( is.na(TET2WT_vs_Normals) & TET2mut_vs_Normals > 0 ) %>% dplyr::select( chr , start , end ) %>% unique
hypoDMR_uniqueTET2mut = dmrTab.annotation %>% filter( ! grepl("^chrM",chr) ) %>% filter( is.na(TET2WT_vs_Normals) & TET2mut_vs_Normals < 0 )  %>% dplyr::select( chr , start , end ) %>% unique

backgroundDMRs = rbind( get( load( file.path( errbsYN_analysisDIR , "TET2WT_vs_Normals_myeDMR_2021-02-05.rda" ) ) ) %>% data.frame %>% dplyr::select( seqnames , start , end ) , 
						get( load( file.path( errbsYN_analysisDIR , "TET2Mut_vs_Normals_myeDMR_2021-02-05.rda" ) ) ) %>% data.frame %>% dplyr::select( seqnames , start , end ) )
backgroundDMRs = backgroundDMRs %>% filter( ! grepl("^chrM",seqnames) )


write.table( hyperDMR_uniqueTET2mut , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHyperDMRs.bed" ) , sep = "\t" , row.names = F , col.names = F , quote = F ) 
write.table( hypoDMR_uniqueTET2mut , file = file.path( errbsYN_analysisDIR , "TET2mut_vs_Normals_uniqueHypoDMRs.bed" ) , sep = "\t" , row.names = F , col.names = F , quote = F )
write.table( backgroundDMRs , file = file.path( errbsYN_analysisDIR , "backgroundDMRs_GREATtool.bed" ) , quote = F , sep = "\t" , row.names = F ,col.names = F )


