#########################################################
## Bayes support functions
## Andy Pohl - 30058005
## University of Calgary
## March 2017
#########################################################

## A collection of support functions for bayesian analysis
## Some gathered from various authors as indiciated in the comments

##

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  # Source: Krushke - https://github.com/boboppie/kruschke-doing_bayesian_data_analysis
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

summarizePost = function( paramSampleVec , 
                          compVal=NULL , ROPE=NULL , credMass=0.95 ) {
  # Summarises the posterior of a parameter from a provided MCMC sample vector
  # Arguments:
  #   paramSampleVec
  #     is a vector of representative values from a probability distribution for a parameter.
  #   compVal
  #      is a scalar of a comparative Value to examine
  #   ROPE
  #      is a vector of length 2 containing the lower and upper bounds for the ROPE
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   Returns a named vector with the posterior summary
  # Source: Krushke - https://github.com/boboppie/kruschke-doing_bayesian_data_analysis
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round( effectiveSize( paramSampleVec ) , 1 )
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                    / length( paramSampleVec ) )
  } else {
    compVal=NA
    pcgtCompVal=NA
  }
  if ( !is.null(ROPE) ) {
    pcltRope = ( 100 * sum( paramSampleVec < ROPE[1] ) 
                 / length( paramSampleVec ) )
    pcgtRope = ( 100 * sum( paramSampleVec > ROPE[2] ) 
                 / length( paramSampleVec ) )
    pcinRope = 100-(pcltRope+pcgtRope)
  } else { 
    ROPE = c(NA,NA)
    pcltRope=NA 
    pcgtRope=NA 
    pcinRope=NA 
  }  
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam , 
             ESS=mcmcEffSz ,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] , 
             CompVal=compVal , PcntGtCompVal=pcgtCompVal , 
             ROPElow=ROPE[1] , ROPEhigh=ROPE[2] ,
             PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) )
}


openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
                        ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file=paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
}


plotPost = function( paramSampleVec , cenTend=c("mode","median","mean")[1],
                     dispCenTend = FALSE, dispHDI = TRUE, disbCRL = FALSE,
                     compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=FALSE , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if ( class(paramSampleVec) == "mcmc.list" ) {
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks=c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border="skyblue"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
    densCurve = density( paramSampleVec , adjust=2 )
    lines(densCurve$x, densCurve$y, type = "l", lwd = 5, col = 'royalblue', bty = 'n')
    lines( c(0,0) , c(1.1*max(histinfo$density),min(histinfo$density)-0.04) , 
           lty="solid" , lwd=1 , col="black" )
    #plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
     #     xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
      #    main=main , cex=cex , cex.lab=cex.lab , ... )
    #polygon(densCurve$x, densCurve$y, col = 'lightgrey', border = col)
  }
  
  # Dispaly values
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display central tendency:
  mn = mean(paramSampleVec)
  med = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  if(dispCenTend){
    if ( cenTend=="mode" ){ 
      #text( mo , cenTendHt ,
      #      bquote(mode==.(signif(mo,3))) , adj=c(.5,0) , cex=cex )
      lines( c(mo,mo) , c(1.1*max(histinfo$density),0) , 
             lty="dotted" , lwd=2 , col="red" )
    }
    if ( cenTend=="median" ){ 
      text( med , cenTendHt ,
            bquote(median==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
    }
    if ( cenTend=="mean" ){ 
      #text( mn , cenTendHt ,
       #     bquote(mean==.(signif(mn,3))) , adj=c(.5,0) , cex=cex )
      lines( c(mo,mo) , c(1.1*max(histinfo$density),0) , 
             lty="dotted" , lwd=2 , col="red" )
    }
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal = sum( paramSampleVec > compVal ) / length( paramSampleVec ) 
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) , 
           lty="dotted" , lwd=2 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(round(100*pLtCompVal,1)) * "% < " *
                    .(signif(compVal,3)) * " < " * 
                    .(round(100*pGtCompVal,1)) * "%" ) ,
          adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"compVal"] = compVal
    postSummary[,"pGtCompVal"] = pGtCompVal
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE = ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE = ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " * 
                    .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " * 
                    .(round(100*pGtROPE,1)) * "%" ) ,
          adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pLtROPE"]=pLtROPE
    postSummary[,"pInROPE"]=pInROPE
    postSummary[,"pGtROPE"]=pGtROPE
  }
  # Display the HDI.
  if(dispHDI){
    lines( c() , c(0,0) , lwd=4 , lend=1 )
    text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
          adj=c(.5,-1.7) , cex=cex )
    text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
          adj=c(HDItextPlace,-0.5) , cex=cex )
    text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
          adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  }
  
  if(disbCRL){
    CRL = quantile(paramSampleVec, prob = c(0.025, 0.975))
    lines( CRL , c(0,0) , lwd=4 , lend=1 )
    text( CRL[1] , 0 , bquote(.(signif(CRL[1],3))) ,
          adj=c(HDItextPlace,-0.5) , cex=cex )
    text( CRL[2] , 0 , bquote(.(signif(CRL[2],3))) ,
          adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  }
  par(xpd=F)
  #
  return( postSummary )
}

# Function for plotting a single horizontal violin within a plot:
plotViolin = function( mcmcSamp , plotHt=0 , dataVal=NULL , dataN=500 ,  
                       cenTend=c("mode","median", "mean")[1] ,
                       cexMin=0.5 , cexMax=1.0 , threshN=700 , slopeN=6/1000 ,
                       topOnly=FALSE , violinWd=c(0.90,1.50)[topOnly+1], 
                       origInt = NULL, axisLim =NULL) {
  # Plot density curve:
  densCurve = density( mcmcSamp )
  polygon( densCurve$x , plotHt+densCurve$y/max(densCurve$y)*(violinWd/2) , 
           col="skyblue" , border="skyblue" )
  if ( !topOnly ) {
    polygon( densCurve$x , plotHt-densCurve$y/max(densCurve$y)*(violinWd/2) , 
             col="skyblue" , border="skyblue" )
  }

  # Plot central tendency of study sIdx:
  mcmcDensity = density(mcmcSamp)
  # 
   if ( cenTend=="mode" ) { 
     mcmcMode = mcmcDensity$x[which.max(mcmcDensity$y)]
     mcmcCenTend = mcmcMode 
     # Plot HDI:
     hdiLim = HDIofMCMC( mcmcSamp )
     lines( hdiLim , rep(plotHt,2) , lend=1 , lwd=4, col = 'black' )
     
  # # Display numerical details of mode and HDI:
  text( par('usr')[2] , plotHt ,
        labels=paste0( format(round(mcmcCenTend,2),nsmall=2)
                       ," (", format(round((hdiLim[1]),2),nsmall=2)
                       ,",", format(round((hdiLim[2]),2),nsmall=2) ,")" ) ,
        adj=c(1.1,-0.5) , cex=0.75 )
  }
  #
  if ( cenTend=="mean" ) { 
    mcmcCenTend = mean(mcmcSamp) 
  #Plot CRI
  CRL = quantile(mcmcSamp, probs = c(0.025, 0.975))
  lines( CRL , rep(plotHt,2) , lend=1 , lwd=4 )
  text( par('usr')[2] , plotHt , 
        labels=paste0( format(round(mcmcCenTend,2),nsmall=2)
                       ," (", format(round((CRL[1]),2),nsmall=2)
                       ,",", format(round((CRL[2]),2),nsmall=2) ,")" ) , 
        adj=c(1.1,-0.5) , cex=0.75 )}
  lines( rep(mcmcCenTend,2) , c(plotHt,plotHt+(violinWd/2)) , 
         lwd=2 , lend=1 )

  
  
  if ( !topOnly ) {
    lines( rep(mcmcCenTend,2) , c(plotHt-(violinWd/2),plotHt) , 
           lwd=2 , lend=1 )
  }

  # Plot data point of study sIdx:
  if (!is.null(origInt)) {

    if(round(origInt[1],5)< axisLim[1] & round(origInt[2],5)> axisLim[2]){
      points(round(axisLim[1],5), plotHt -.1, 
             pch=-9668 , col="black" , lwd=1)
      points(round(axisLim[2],5), plotHt -.1, 
             pch=-9658 , col="black" , lwd=1)
      lines(round(axisLim,5)-.1 , rep(plotHt-.1,2) , lend=1 , lwd=2, lty = 2 )
    }else if(round(origInt[1],5)< axisLim[1]){
      points(round(axisLim[1],5), plotHt -.1, 
             pch=-9668 , col="black" , lwd=1)
      points(origInt[2], plotHt -.1, 
             pch="|" , col="black" , lwd=1)
      lines(c(round(axisLim[1],5)-.1, origInt[2]-.1) , rep(plotHt-.1,2) , lend=1 , lwd=2, lty = 2 )
    } else if(round(origInt[2],5)> axisLim[2]){
      points(origInt[1], plotHt -.1, 
             pch="|" , col="black" , lwd=1)
      points(round(axisLim[2],5), plotHt -.1, 
             pch=-9658 , col="black" , lwd=1)
      lines(c(origInt[1], round(axisLim[2],5)-.1) , rep(plotHt-.1,2) , lend=1 , lwd=2, lty = 2 )
    } else {
      points(origInt[1], plotHt -.1, 
             pch="|" , col="black" , lwd=1)
      points(origInt[2], plotHt -.1, 
             pch="|" , col="black" , lwd=1)
      lines(origInt , rep(plotHt-.1,2) , lend=1 , lwd=2, lty = 2 )
    }
    
  }
  if(!is.null(dataVal)){
  pChar=15
  points(dataVal, plotHt -.1, 
         pch=pChar , col="grey30" , lwd=1 , bg="grey" ,
         cex = cexMin+(cexMax-cexMin)/(1+exp(-slopeN*(dataN-threshN))) )
  }
}
