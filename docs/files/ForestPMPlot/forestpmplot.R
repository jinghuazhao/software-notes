
# This R code is for ploting Forest plots and PM-Plots with
# genome-wide association meta analysis result from MetaSoft.

# Copyright (C) 2014  Eun Yong Kang (eunyong.kang@gmail.com)

#The program is free for academic use. Please contact Eun Yong Kang
#<eunyong.kang@gmail.com> if you are interested in using the software for
#commercial purposes.

#The software must not be modified and distributed without prior
#permission of the author.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



################################################################################
#
#  This R code requires 7 inputs
#
#  1. filename containing selected rows of metasoft output file
#  2. filename containing selected rows of metasoft input file
#  3. filename containing rsids and candidate Gene name for selected SNPs
#  4. filename continaing studynames.
#  5. filename continaing display index of the studynames.
#  6. maximum -log10(p-value) for PM-plot
#  7. height of the output file in inches

library(plotrix)
## library(maptools)

meta.colors<-function(all.elements,box="black",lines="gray",summary="black",zero="lightgray",
                      mirror="lightblue",text="black", axes="black",background=NA){

    if (missing(all.elements)){
        return(list(box=box, lines=lines, summary=summary,
             zero=zero, mirror=mirror, text=text,
             axes=axes, background=background))
    }

    if (is.null(all.elements))
        all.elements<-par("fg")

    return(list(box=all.elements, lines=all.elements,
             summary=all.elements, zero=all.elements,
             mirror=all.elements, text=all.elements,
             axes=all.elements, background=NA))

}



metaplot <- function( mn, se, nn=NULL, labels=NULL, conf.level = .95,
		      xlab = "Odds ratio", ylab = "Study Reference",
		       xlim = NULL, summn = NULL,
		      sumse = NULL, sumnn = NULL,
		      summlabel = "Summary", logeffect = FALSE,
		      lwd = 2, boxsize = 1,
		      zero = as.numeric(logeffect),
                      colors=meta.colors(), xaxt="s", logticks=TRUE,
		      ... ) {
    nth<-function(x,i){
        x[ (i-1) %% length(x) +1]
    }
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    ok <- is.finite( mn + se )
    if ( is.null( xlim ) )
        xlim <- c( min( mn[ok] - ci.value * se[ok], na.rm = TRUE )-0.00,
      		     max( mn[ok] + ci.value * se[ok], na.rm = TRUE ) )
    ##par( pty="s" )
    n <- length( mn )
    if ( logeffect ) {
        xlog <- "x"
        nxlim <- exp( xlim )
    }
    else {
        xlog <- ""
        nxlim <- xlim
    }

    leftedge<-nxlim[1]
    lalim = nxlim

    if ( !is.null( labels ) ) {
        if ( logeffect )
            nxlim[1] <- nxlim[1] / sqrt( nxlim[2] / nxlim[1] )
        else
          #nxlim[1] <- nxlim[1] - 0.5 * ( nxlim[2] - nxlim[1] )
          #lalim[1] <- nxlim[1] - 0.6 * ( nxlim[2] - nxlim[1] )
          lalim[1] <- nxlim[1] - 0.33 * ( nxlim[2] - nxlim[1] )
          nxlim[1] <- nxlim[1] - 0.7 * ( nxlim[2] - nxlim[1] )

        labels<-as.character(labels)

    }
    par( xaxt = "n",yaxt = "n", bg=colors$background )
    plot( nxlim,c( 1,-n-2-3 * !is.null( summn ) ),
          type = "n", bty = "n", xaxt = "n", yaxt = "n",
          log = xlog, xlab=xlab, ylab=ylab,..., col.lab=colors$axes )

    par( xaxt = "s" )
    if (xaxt=="s"){
        if (logeffect) {
            if (logticks){
                ats<-round( 10 ^ pretty( log( exp( xlim ),10), 8,min.n=6  ), 2 )
                ats<-ats[ats> exp(xlim[1]) & ats< 10^(par("usr")[2])]
                axis( 1, at = ats, col= colors$axes, col.axis= colors$axes)
            } else {
                ats<-pretty(exp(xlim),8, min.n=6)
                ats<-ats[ats> exp(xlim[1]) & ats <10^(par("usr")[2])]
                axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
            }
        }  else {
            ats<-pretty(xlim, 6)
            ##ats<-ats[ats> xlim[1] & ats <xlim[2]]
            axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
        }
    }

    if ( !is.null( zero )&& zero>leftedge )
        abline( v = zero, lty = 2, lwd = 2 ,col=colors$zero)

    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    lower <- mn - ci.value * se
    upper <- mn + ci.value * se
    if ( logeffect ){
        lower <- exp( lower )
        upper <- exp( upper )
    }
    for ( i in 1:n ){
        if ( is.na( lower[i]+upper[i] ) )
            next
        lines( c( lower[i], upper[i] ), c( -i, -i ), lwd = lwd, col=nth(colors$lines,i),... )
    }

    if ( !is.null( labels ) )
        #text( rep( nxlim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0 )
        #text( rep( lalim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=-0.4 )
        #text( rep( lalim[1]+0.4, n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0 )
        ## added by yurang.park   	
        if(max(nchar(studies))>5){
            mxstrlenlabel <- max(strwidth(labels))-0.4    		    
        }else{
        		mxstrlenlabel <- 0         
        }        
        text( rep( lalim[1]-mxstrlenlabel, n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0)            
        

    if ( is.null( nn ) )
        nn <- se ^ -2
    yscale <- 0.3 * boxsize / max( sqrt( nn ), na.rm = TRUE )

    if ( logeffect ) {
        scale <- ( nxlim[2] / nxlim[1] ) ^ ( yscale / ( 4 + n ) )
        xl <- exp( mn ) * ( scale ^ -sqrt( nn ) )
        xr <- exp( mn ) * ( scale ^ sqrt( nn ) )
    }
    else {
        scale <- yscale * ( nxlim[2] - nxlim[1] ) / ( 4 + n )
        xl <- mn - scale * sqrt( nn )
        xr <- mn + scale * sqrt( nn )
    }
    yb <- ( 1:n ) - yscale * sqrt( nn )
    yt <- ( 1:n ) + yscale * sqrt( nn )
    for ( i in 1:n ) {
        if ( !is.finite( mn[i] ) )
            next
        rect( xl[i], -yb[i], xr[i], -yt[i], col = nth(colors$box,i),border=nth(colors$box,i))
    }
    if ( !is.null( summn ) ) {
        if ( logeffect ) {
            x0 <- exp( summn )
            xl <- exp( summn - ci.value * sumse )
            xr <- exp( summn + ci.value * sumse )
        }
        else{
            x0 <- summn
            xl <- summn - ci.value * sumse
            xr <- summn + ci.value * sumse
        }
        y0 <- n + 3
        yb <- n + 3 - sqrt( sumnn ) * yscale
        yt <- n + 3 + sqrt( sumnn ) * yscale
        polygon( c( xl, x0, xr, x0 ), -c( y0, yt, y0, yb ),
    	         col = colors$summary, border = colors$summary )
        #text( nxlim[1], -y0, labels = summlabel, adj = 0,col=colors$text )
        text( lalim[1], -y0, labels = summlabel, adj = 0,col=colors$text )
    }
}



funnelplot.default<-function(x,se,size=1/se,summ=NULL,xlab="Effect",ylab="Size",
		colors=meta.colors(),
		conf.level=0.95,plot.conf=FALSE,zero=NULL,mirror=FALSE,...)
{
   finite<-function(x) x[is.finite(x)]

   if (mirror && is.null(summ))
	stop("Can't do a mirror plot without a summary value")

   if (plot.conf){
	ci<--qnorm((1-conf.level)/2)
	xlim<-range(finite(c(zero,x-ci*se,x+ci*se)))
	if (mirror)
	  xlim<-range(finite(c(xlim,2*summ-x-ci*se,2*summ-x+ci*se)))
   }else{
      xlim<-range(finite(c(zero,x)))
      if (mirror)
	xlim<-range(finite(c(xlim,2*summ-x,2*summ+x)))
   }
   plot(x,size,ylim=c(0,max(size)*1.1),xlim=xlim,xlab=xlab,
	ylab=ylab,col=if(is.null(colors$points)) par("fg") else colors$points)

   if (plot.conf)
       segments(x-ci*se,size,x+ci*se,size,col=if(is.null(colors$conf)) par("fg") else colors$conf,lwd=2)
   if (!is.null(summ))
	   abline(v=summ,col=if(is.null(colors$summary)) par("fg") else colors$summary,lty=2,lwd=2)

   if(!is.null(zero))
	abline(v=zero,col=if(is.null(colors$zero)) par("fg") else colors$zero,lwd=2)

   if(mirror){
	points(2*summ-x,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror)
       if (plot.conf)
	segments(2*summ-x-ci*se,size,2*summ-x+ci*se,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror,lwd=2)
  }
}



extract_output <- function(outtab) {
    dim(outtab) <- c(1, length(outtab))
    ncol = dim(outtab)[2]
    nsnp = dim(outtab)[1]
    mat <- as.matrix(outtab[,17:ncol]) # p-m data
    if (nsnp == 1) {
        dim(mat) <- c(1, length(mat))
    }
    ncol <- ncol - 16
    p <- as.numeric(mat[,1:(ncol/2)])
    h <- as.numeric(mat[,(ncol/2+1):ncol])
    ptab = as.numeric(outtab[, 9])
    summary <- as.numeric(outtab[,7]) # random effects model summary
    summaryse <- as.numeric(outtab[,8])
    return(list(p=p, h=h, ptab=ptab, summary=summary, summaryse=summaryse))
}



extract_input <- function(intab) {
    dim(intab) <- c(1, length(intab))
    # read data for forest plot
    #data <- intab[sel, ]
    ncol = dim(intab)[2]
    nsnp = dim(intab)[1]
    mat <- as.matrix(intab[,2:ncol])
    if (nsnp) {
        dim(mat) <- c(1, length(mat))
        effects <- mat[,c(TRUE,FALSE)]
        stderrs <- mat[,c(FALSE,TRUE)]
        dim(effects) <- c(1, length(effects))
        dim(stderrs) <- c(1, length(stderrs))
    }
    else {
        effects <- mat[,c(TRUE,FALSE)]
        stderrs <- mat[,c(FALSE,TRUE)]
    }
    return(list(effects=effects, stderrs=stderrs))
}





pval <- function(p) {
    if (is.na(p)) {
        return("")
    }
    if (p < 0.0001) {
        pp = as.numeric(p)
        ps = as.character(pp)
        pa = paste(substr(ps, 1, 4), "x")
        po = substr(ps, 10, 12)
        po = as.character(as.numeric(po) *(-1))
        #rsid <- bquote(bold(.(pa)) ~ bold("10")^bold(.(po)) )
        rsid <- bquote(.(pa) ~ "10"^.(po) )
        return(rsid)
    }
    else {
        pa = as.character(p)
        pa = paste("  ", substr(pa, 1, 6), sep="")
        #rsid <- bquote(bold(.(pa)))
        rsid <- bquote(.(pa))
        return(rsid)
    }
}

## New function to spread text to not overlap ## Buhm 2/5/14
spread <- function(pos, mingap) {
  ord <- order(pos)
  pos.spread <- pos
  previous <- -99999999
  for (i in 1:length(pos)) {
    pos.spread[i] <- max(previous + mingap, pos[ord][i])
    previous <- pos.spread[i]
  }
  return(pos.spread[order(ord)])
}

pmplot <- function(maxylim=10, p, h, ptab, summary, summaryse, effects, stderrs, newrsid, studies, change, height, pthreshold, pmplot_file) {

    i = 1
    nsnps = length(ptab)
    flip <- rep(1, nsnps)
    studies.short <- 1:length(studies)
    cstudies = studies
    #change = 1:length(studies)
    #change <- c(1, 10, 2, 11, 3, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16, 9, 17)
    for (k in 1:length(cstudies)) {
        cstudies[k] = paste(which(change == k), ".", cstudies[k], sep="")
    }

    id <- 1:length(studies)
    for (k in 1:length(cstudies)) {
        id[k] = which(change == k)
    }
    print(pmplot_file)
    
    # Forest plot
    ## added by yurang.park
    maxsl<-max(nchar(studies))*0.2
    pdf(pmplot_file, width=11.4+maxsl, height=height)
    
    pp = min(as.numeric(ptab[i]))
    ps = as.character(pp)
    pa = paste(substr(ps, 1, 4), "x")
    po = substr(ps, (nchar(ps)-1), nchar(ps))
    po = as.character(as.numeric(po) *(-1))
    if (is.na(newrsid[i, 2])) {
        rsid <- bquote(bold(.(newrsid[i,1])) ~ bold(" (Meta ") ~ bolditalic(P) ~ bold(.(" = ")) ~ bold(.(pa)) ~ bold("10")^bold(.(po)) ~ bold(")") )
    }
    else {
        genename = newrsid[i,2]
        rsid <- bquote(bold(.(newrsid[i,1])) ~ bold(" (Meta ") ~ bolditalic(P) ~ bold(.(" = ")) ~ bold(.(pa)) ~ bold("10")^bold(.(po)) ~ bold(")") )
    }
    #layout(matrix(c(1,2,3,4),1,2,byrow=TRUE),widths=c(6.3,5.5))
    layout(matrix(c(1,2,3,4),1,2,byrow=TRUE),widths=c(7.8+maxsl,5.9))
    par(mar=c(5,0,4,1)+0.1)
    metaplot(as.numeric(effects[i,change])*flip[i], as.numeric(stderrs[i,change]), 
             labels=cstudies[change], 
             xlab="Log odds ratio", ylab="",
             summn=as.numeric(summary[i])*flip[i],sumse=as.numeric(summaryse[i]),
             sumnn=1/(as.numeric(summaryse[i])^2),
             summlabel="RE Summary",
             main=rsid)
    off = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0)
    if (!is.na(newrsid[i, 2])) {
        genename = newrsid[i,2]
        mtext(bquote(bold(.("Gene : ")) ~ bolditalic(.(genename))),line=0, cex=1.20, adj=(0.118+off[i])-0.0685+0.192)
    }
    # color dots
    # compute x-location
    ## added by yurang.park
    if(max(nchar(studies))>5){
    	mxstrlenlabel <- max(strwidth(studies))-0.35  
    }else{
    	mxstrlenlabel <- 0
    } 
    mn = as.numeric(effects[i,change])*flip[i]
    se = as.numeric(stderrs[i,change])
    conf.level = 0.95
    ci.value <- -qnorm((1 - conf.level)/2)
    ok <- is.finite(mn + se)
    xlim <- c(min(mn[ok] - ci.value * se[ok], na.rm = TRUE), max(mn[ok] + ci.value * se[ok], na.rm = TRUE))
    nxlim <- xlim
    
    x_loc = (nxlim[1] - 0.3962 * (nxlim[2] - nxlim[1]))- mxstrlenlabel
    inds <- which(h[change] < 0.9 & h[change] > 0.1)
    inds <- inds*-1 + 0.1 - 0.08
    points(rep(x_loc, length(inds)),inds,
         col="limegreen", cex=1.5, pch=19)
    inds <- which(h[change] >= 0.9)
    inds <- inds*-1 + 0.1 - 0.08
    points(rep(x_loc, length(inds)),inds,
         col="#f31515", cex=1.5, pch=19)
    inds <- which(h[change] <= 0.1)
    inds <- inds*-1 + 0.1 - 0.08
    points(rep(x_loc, length(inds)),inds,
         col="blue", cex=1.5, pch=19)
    #s_x_loc = nxlim[1] - 0.2942 * (nxlim[2] - nxlim[1])
    #s_x_loc = nxlim[1] - 0.2942 * (nxlim[2] - nxlim[1])
    s_x_loc = (nxlim[1] - 0.1242 * (nxlim[2] - nxlim[1])) - (mxstrlenlabel+0.1)
    text(s_x_loc, 0, expression(bold("Study Name")))
    
    pv=bquote(bolditalic("P") ~ bold("-value"))
    p_x_loc = (nxlim[1] - 0.5360 * (nxlim[2] - nxlim[1])) - (mxstrlenlabel-0.05)
    text(p_x_loc, 0.06, pv, cex=1.1)

    p_x_loc = (nxlim[1] - 0.5642 * (nxlim[2] - nxlim[1]))- (mxstrlenlabel-0.05)
    for (k in 1:length(change)) {
        y_loc = k*-1.00 + 0.1 - 0.08
        text(p_x_loc, y_loc, pval(p[change[k]]), cex=.8)
    }


    # P-M plot
    par(mar=c(5,3,4,2.5)+0.1)
 
    plot(0, 0,
         type="n",xlim=c(0,1), ylim=c(0,maxylim),yaxt='s',
         #xlab="m-value",ylab="-log(p-value)", bty="n",
         xlab="m-value",ylab="", bty="n",
         main = "PM-Plot")
    ##modified by yurang.park
    ##COL.darkgray = "#E0E0E0"
    COL.lightgray = "#F3F3F3"
    rect(0, 0, 0.1, maxylim, col=COL.lightgray, border=NA)
    rect(0.1, 0, 0.9, maxylim, col=COL.lightgray, border=NA)
    rect(0.9, 0, 1, maxylim, col=COL.lightgray, border=NA)
    abline(v=axTicks(1), col="white")
    abline(h=axTicks(2), col="white")
    abline(v=c(0.1, 0.9), col="white")

    ## variable dot size
    max.CEX = 3
    min.CEX = 1
    approxN = as.numeric(stderrs[i,])^-2
    scale.CEX = approxN / max(approxN[!is.na(approxN)]) * max.CEX
    scale.CEX[scale.CEX < min.CEX] = min.CEX

    ## DOTS!!
    x.TEXTGAP = 1/10
    y.TEXTGAP = maxylim/15
    line.color = "#808080"
    ## green dots
    inds = 0
    tinds <- which(h[i,] < 0.9 & h[i,] > 0.1)
    counts = 1
    for (w in 1:length(tinds)) {
        if (length(which(tinds[w] == change)) != 0) {
            inds[counts] = tinds[w]
            counts = counts + 1
        }
    }
    if (length(inds) > 1 || (inds != 0 && length(inds) > 0)) {
      inds <- inds[order(scale.CEX[inds], decreasing=T)]
      points(h[i,inds],-log10(p[i,inds]),
             bg="limegreen", cex=scale.CEX[inds], pch=21, col=line.color)
      thigmophobe.labels(x=h[i,inds],y=-log10(p[i,inds]),labels=id[inds], cex = 0.75)
    }
    ## red dots
    inds = 0
    tinds <- which(h[i,] >= 0.9)
    counts = 1
    for (w in 1:length(tinds)) {
        if (length(which(tinds[w] == change)) != 0) {
            inds[counts] = tinds[w]
            counts = counts + 1
        }
    }
    if (length(inds) > 1 || (inds != 0 && length(inds) > 0)) {
      inds <- inds[order(scale.CEX[inds], decreasing=T)]
      xpos <- h[i, inds]
      ypos <- -log10(p[i,inds])
      ypos.spread <- spread(ypos, maxylim/50)
      for (j in 1:length(inds)) {
        lines(c(xpos[j], xpos[j]-x.TEXTGAP), c(ypos[j], ypos.spread[j]), col=line.color)
        points(xpos[j], ypos[j],
               bg="#f31515", cex=scale.CEX[inds[j]], pch=21, col=line.color)
        text(x=xpos[j]-x.TEXTGAP, y=ypos.spread[j], labels=id[inds[j]], cex=0.75, pos=2)
      }
    }
    ## blue dots
    inds = 0
    tinds <- which(h[i,] <= 0.1)
    counts = 1
    for (w in 1:length(tinds)) {
        if (length(which(tinds[w] == change)) != 0) {
            inds[counts] = tinds[w]
            counts = counts + 1
        }
    }
    if (length(inds) > 1 || (inds != 0 && length(inds) > 0)) {
      inds <- inds[order(scale.CEX[inds], decreasing=T)]
      xpos <- h[i, inds]
      ypos <- -log10(p[i,inds])
      xpos.spread <- spread(xpos, 1/20)
      for (j in 1:length(inds)) {
        lines(c(xpos[j], xpos.spread[j]), c(ypos[j], ypos[j]+y.TEXTGAP), col=line.color)
        points(xpos[j], ypos[j],
               bg="blue", cex=scale.CEX[inds[j]], pch=21, col=line.color)
        text(x=xpos.spread[j], y=ypos[j]+y.TEXTGAP, labels=id[inds[j]], cex=0.75, pos=3)
      }
    }

    abline(h=-log10(pthreshold),lty="dashed")

    title(ylab=expression(-log[10](p)),mgp=c(2,1,0), xpd=F)
    #title(mgp=c(3,1,0), xpd=F)

    ## pthreshold = 1E-6


    # key explation modified by yurang.park
    #st = 8
    #ste = st -0.03
    #points(0.08, st, col="#f31515", cex=1.5, pch=19)
    #text(0.360, ste, "Study has an effect (m > .9)", cex = 0.8)
    #points(0.08, st-0.5, col="blue", cex=1.5, pch=19)
    #text(0.445, ste-0.5, "Study does not have an effect (m < .1)", cex=0.8)
    #points(0.08, st-1, col="limegreen", cex=1.5, pch=19)
    #text(0.443, ste-1, "Study's effect is uncertain (.1< m < .9)", cex=0.8)

    dev.off()
}



forest_pm_plot <- function(outtab, intab, rsids, genename, studies, change, maxylim=10, height, pthreshold=1E-6, pmplot_file) {

    outs = extract_output(outtab)
    p = outs$p
    h = outs$h
    dim(h) <- c(1, length(h))
    dim(p) <- c(1, length(p))
    ptab = outs$ptab
    summary = outs$summary
    summaryse = outs$summaryse


    ins = extract_input(intab)
    effects = ins$effects
    stderrs= ins$stderrs

    newrsid <- c(rsids, genename)
    #print(newrsid)
    dim(newrsid) <- c(1, 2)


    pmplot(maxylim, p, h, ptab, summary, summaryse, effects, stderrs, newrsid, studies, change, height, pthreshold, pmplot_file)

}



args=(commandArgs(TRUE))
output_file = args[1]
input_file = args[2]
gene_file = args[3]
studynames_file = args[4]
studyorder_file = args[5]
pmplot_file = args[6]
height= as.numeric(args[7])



outtab = as.matrix(read.table(output_file))
intab = as.matrix(read.table(input_file))
rsids = outtab[1]
genes = as.matrix(read.table(gene_file))
studies <- as.matrix(read.table(studynames_file))
change <- scan(studyorder_file)

forest_pm_plot(outtab, intab, rsids, genes, studies, change, maxylim=10, height, pthreshold=1E-6, pmplot_file)










