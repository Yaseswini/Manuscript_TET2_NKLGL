## Helper functions ## 


source( "/sfs/qumulo/qproject/gblab/users/yn9w/Scripts/errbs/errbs_analysis.R" )
source( "/sfs/qumulo/qproject/gblab/users/yn9w/Scripts/annotate_loci_old.R" )


## Run methylkit ## 
## ------------------------------------------------------------ ## 
##                               methylkit                                                      ## 
## ------------------------------------------------------------ ## 
## Description : run methylkit                  
## Input : file.list -> list of methylcall files
##                 sample.id -> list of sample names 
##                 assembly -> assembly ( hg19 etc )
##                 treatment -> vector of 1s and 0s indicating treatment / control
##                 context -> CpG 
##                 outLabel -> label of the output file
## Produces 3 output files : one with output , args and sessionInfo
run_methylkit <- function( file.list , sample.id , assembly , treatment , context = "CpG" , min.per.group = 2L , num.cores = 10 , outLabel )
{
    if( !is.list(file.list) & missing(file.list) & !is.list(sample.id) & missing(sample.id) ) { stop( "Specify files and samplenames " ) }
    if( missing(treatment) ) { stop( "Specify treatment" ) }

    library(methylKit)
    #cat( "*** Args *** \n" , file.list , "\n----\n" , sample.id ,"\n-----\n" , assembly , "\n---\n", treatment,"\n----\n" ,context,"\n---\n",num.cores,"\n")
    myobj = methRead(file.list , sample.id = sample.id , assembly=assembly , treatment=treatment , context=context)
    print( "myobj created" )

    filtered.myobj = filterByCoverage(myobj, lo.count=10,lo.perc=NULL, hi.count=500,hi.perc=NULL)
    print( "filtered myobj" )

    filtered.myobj = normalizeCoverage(filtered.myobj)
    print( "filtered myobj - normalized for Coverage" )

    meth=unite(filtered.myobj,destrand=FALSE,min.per.group=min.per.group)
    print( "meth obj created" )

    myDiff=calculateDiffMeth(meth,mc.cores=num.cores)
    print( "myDiff calculated" )

    outfile = paste( outLabel , "_myDiff_" , Sys.Date() , ".rda" , sep = "" )
    save( myDiff , file = outfile )

    # print a log file with commands 
    # modified on 10th march , 2017
    args.outfile = paste( file_path_sans_ext(outfile) , "_methylkitArgs_" , Sys.Date() , ".txt" , sep = "" )

    dat = data.frame( "file.list" = unlist( file.list ) , "sample.id" = unlist( sample.id ) , "treatment" = treatment )
    write.table( dat , file = args.outfile , sep = "\t" , row.names = F , col.names = T , quote = F )

    additional.args = paste( "assembly = " , assembly , "\n" ,
                             "context = " , context , "\n" ,
                             "min.per.group = " , min.per.group , "\n" ,
                             "num.cores = " , num.cores , "\n" ,
                             "outfile = " , outfile , "\n"  , sep = "" )

    write( additional.args , file = args.outfile , sep = "\t" , append = T )

    # print sessionInfo 
    writeLines( capture.output( sessionInfo() ) , paste( outfile , "_sessionInfo" , "_" , Sys.Date() , sep = "" ) )

    return(myDiff)
}


## 
## ------------------------------------------------------------ ## 
##                               edmr                                                           ## 
## ------------------------------------------------------------ ## 
## Description : run edmr               
## Input : myDiff -> output of methykit 
##                 outfile.label -> label of the output file
run_edmr <- function( myDiff , DMC.qvalue = 0.01 , DMC.methdiff = 25 , outfile.label )
{
   # packages 
   library(edmr)
     # Function from Sheng's github :- edmr.R
    # edmr=function(myDiff, step=100, dist="none", DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, plot=FALSE, main="", mode=1, ACF=TRUE, fuzzypval=1){
    #myDiff=as.data.frame(myDiff)
    #if(mode==1){
    #  myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, 
    #                 DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="both", ACF=ACF)
    #} else if (mode==2){
    #  hyper.myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs,
    #                       DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="hyper", ACF=ACF)
    #  hypo.myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs,
    #                      DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="hypo", ACF=ACF)
    #  myDMR=c(hyper.myDMR, hypo.myDMR)
    #}
    #else {
    #  stop ("mode = 1 or 2")
    #}
    #
  # ..................................

  myMixmdl = myDiff.to.mixmdl(myDiff)   # Fits a mixture model 
  # Background DMRs 
   bkgd.DMC.qvalue = 1;
   bkgd.DMC.methdiff = 0;
   bkgd.num.DMCs = 0;
   bkgd.num.CpGs = 3;
   bkgd.DMR.methDiff = 0;
   bkgd.plot = FALSE;
   bkgd.mode = 2;
   bkgd.ACF = T;
   bkgd.fuzzypval = 1

   myDMR.background = edmr( myDiff = as.data.frame(myDiff) ,
                            step = 100 ,
                            dist = "none" ,
                            DMC.qvalue = bkgd.DMC.qvalue ,
                            DMC.methdiff = bkgd.DMC.methdiff ,
                            num.DMCs = bkgd.num.DMCs ,
                            num.CpGs = bkgd.num.CpGs ,
                            DMR.methdiff = bkgd.DMR.methDiff ,
                            mode = 2 ,
                            ACF = T ,
                            fuzzypval = bkgd.fuzzypval )

   save( myDMR.background , file = paste( outfile.label ,"_background_", Sys.Date(), ".rda" , sep="" ) )
   # edmr
   num.CpGs = 3 ; DMR.methDiff = 1 ;
   myDMR = edmr( myDiff = myDiff , step = 100 , dist = "none" ,
                 DMC.qvalue = DMC.qvalue ,
                 DMC.methdiff = DMC.methdiff ,
                 num.DMCs = 3 ,
                 num.CpGs = num.CpGs ,
                 DMR.methdiff = DMR.methDiff ,
                 plot = FALSE ,
                 mode = 2 ,
                 ACF = TRUE ,
                 fuzzypval = 0.1 )

    save( myDMR , file = paste( outfile.label, "_myeDMR_" , Sys.Date() , ".rda",sep="") )

      # arguments 
     backgroundArgs = paste( "# Background DMR args : " , sep = "\n" )
     backgroundArgs = paste( backgroundArgs , "\n" , "DMC.qvalue = " , bkgd.DMC.qvalue , sep = "" )
     backgroundArgs = paste( backgroundArgs , "\n" , "DMC.methdiff = " , bkgd.DMC.methdiff , sep = "" )
     backgroundArgs = paste( backgroundArgs , "\n" , "num.DMCs = " , bkgd.num.DMCs , sep  = "" )
     backgroundArgs = paste( backgroundArgs , "\n" , "num.CpGs = " , bkgd.num.CpGs , sep = "" )
     backgroundArgs = paste( backgroundArgs , "\n" , "DMR.methdiff = " , bkgd.DMR.methDiff , sep = "" )
     backgroundArgs = paste( backgroundArgs , "\n" , "mode = " , bkgd.mode , "\n" , "ACF = T" , "\n" , "fuzzypval = " , bkgd.fuzzypval , sep = "" )

     edmrArgs = paste( "# edmr args : " , sep = "\n" )
     edmrArgs = paste( edmrArgs , "\n" , "DMC.qvalue = " , DMC.qvalue , sep = "" )
     edmrArgs = paste( edmrArgs , "\n" , "DMC.methdiff = " , DMC.methdiff , sep = "" )
     edmrArgs = paste( edmrArgs , "\n" , "num.DMCs = 3" , "\n" , "plot = FALSE", "\n" , "mode = 2" , "\n" , "ACF = TRUE" , "\n" , "fuzzypval = 0.1" )

     args = c( backgroundArgs , edmrArgs )

     cat( args , file = paste( outfile.label , "_edmrArgs.txt" , sep = "" ) , append = F )

      argsDf = data.frame( sampleID = myDiff@sample.ids , treatment = myDiff@treatment , assembly = rep( myDiff@assembly , length(myDiff@sample.ids) ) )
      write.table( argsDf , file = paste( outfile.label , "_edmrArgs.txt" , sep = "" ) , append = T , sep = "\t" , quote = F , row.names = F , col.names = T )


writeLines( capture.output( sessionInfo() ) , paste( outfile.label , "_sessionInfo_" , Sys.Date() , sep = "" ) )

}



## Prepare data for methylation plot 
prepare_data_dmc_barplot <- function( methDat , colOrder )
{  
  tab = methDat
  if( is.element("chrBase",colnames(methDat)) ) {
    rownames(tab) = tab$chrBase ;
    tab$chrBase = NULL
  }
  temp1 = tab[ , colOrder ]
  temp1$NACount = apply( tab , 1 , function(x) { return( length( x[is.na(x)] ) ) } )
  tab = temp1[ which( temp1$NACount < length(colOrder) ) , ]
  #
  temp = setNames( data.frame( rownames(tab) , tab[ , colOrder[[1]] ] ) , c("chrBase" , colOrder[[1]]) )
  temp$effect = "1Hypo"
  temp[ which(temp[[ colOrder[[1]] ]] > 0 ) , "effect" ] = "2Hyper"
  temp[ which( is.na( temp[[ colOrder[[1]] ]] ) ) , "effect" ] = "3NA"
  temp[ , colOrder[[1]] ] = NULL
  #altunaDat = replaceColName("effect",colOrder[[1]],temp)
  plotDat = temp
  colnames( plotDat )[ colnames( plotDat ) == "effect" ] = colOrder[[1]]
  
  for( i in 2:length(colOrder) )
  {   
    temp = setNames( data.frame( rownames(tab) , tab[ , colOrder[[i]] ] ) , c("chrBase" , colOrder[[i]]) )
    temp$effect = "1Hypo"
    temp[ which(temp[[ colOrder[[i]] ]] > 0 ) , "effect" ] = "2Hyper"
    temp[ which( is.na( temp[[ colOrder[[i]] ]] ) ) , "effect" ] = "3NA"
    temp[ , colOrder[[i]] ] = NULL
    #temp = replaceColName("effect",colOrder[[i]],temp)
    colnames( temp )[ colnames( temp ) == "effect" ] = colOrder[[i]]
    
    plotDat = merge( plotDat , temp , by = "chrBase" , all = T )
  }
  
  #altunaDat$chrBase = NULL 
  return(plotDat)
}

dmc_barplot <- function( tab , colOrder , xlabName = "" , plotTitle = "" )
{
  plotDat = prepare_data_dmc_barplot( methDat = tab , colOrder )
  rownames(plotDat) = plotDat$chrBase
  plotDat$chrBase = NULL

  cat( nrow(tab) , "\t" , nrow(plotDat) , "\n" )

  plotDat = plotDat[ do.call( order , as.list( plotDat ) ) , ]
  plotDat[ plotDat == "1Hypo" ] = "Hypo"
  plotDat[ plotDat == "2Hyper" ] = "Hyper"
  plotDat[ plotDat == "3NA" ] = "NA"
  plotDat[ , "id" ] = c( 1:nrow( plotDat ) )
  plotDat = melt( plotDat , "id" )
  plotDat$variable = factor( plotDat$variable , levels = colOrder )

  p = ggplot(plotDat,aes(x=id,y=variable,fill=value))+geom_tile()+
    scale_fill_manual(values=c("Hyper"="gold2","Hypo"="blue3","NA"="gray86")) +
    theme_bw() +
    labs(y="",fill="") +
    theme( plot.title = element_text(size=10) ,
           panel.grid.major = element_blank() ,
           panel.grid.minor = element_blank() ,
           axis.text = element_text(colour = "black",size = 6) ,
           axis.title = element_text( size = 6.5 ) ,
           legend.text=element_text(size = 6.5,colour="black"),legend.position="bottom")

  p = p + labs( x = xlabName , title = plotTitle )

  return(p)
}