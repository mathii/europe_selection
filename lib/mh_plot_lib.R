## Very sorry about these hardcoded paths.
## if(Sys.info()[4]=="dhcp179.well.ox.ac.uk"){
##   source("~/mathii/gwas/whr-strat-t2d/meta_hits.R")
## }else{
##   source("~/gwas/whr-strat-t2d/meta_hits.R")
## }


MH.plot <- function(data, flip=FALSE, add=FALSE, chr.labels=TRUE, shift = 0, cap = 100, color.loci=data.frame(),range.shift=2, ...){

  ## Deal with X chromosome
  if(typeof(data$CHR)=="character"){
    cat("Converting CHR field to numeric\n")
    data[data$CHR=="X","CHR"] <- "23"
    data[data$CHR=="Y","CHR"] <- "24"
    data$CHR <- as.numeric(data$CHR)
  }
             
  data<- data[data$PVAL>-1,]
  
  obspval <- (data$PVAL)
  chr <- (data$CHR)
  pos <- (data$POS)
  obsmax <- min(cap, trunc(max(-log10(obspval))))+range.shift
  
  sort.ind <- order(chr, pos)  # sorting index
  chr <- chr[sort.ind]
  pos <- pos[sort.ind]
  obspval <- obspval[sort.ind]
  
  x <- get.chr.starts(chr,pos)

  locX = trunc(pos/100) + x[chr]
  locY = -log10(obspval)
  plot.sym <- ifelse( locY>cap, 17, 20 )

  ## Set colours - alternating per chromosome and red for off the scale
  col1=rgb(0,0,108,maxColorValue=255)
  col2=rgb(100,149,237,maxColorValue=255)
  col3=rgb(205,50,50,maxColorValue=255)
  col4 <- ifelse (chr%%2==0, col1, col2)
  curcol <- ifelse (locY>cap, col3, col4) 
  pcl <- rep(FALSE, length(curcol))

  if(NROW(color.loci)>0){
    for(i in 1:NROW(color.loci)){
      pos.from<-color.loci[i,"pos.from"]
      pos.to <- color.loci[i,"pos.to"]
      this.chr <- color.loci[i,"chr"]
      curcol[chr==this.chr & pos>=pos.from & pos<=pos.to] <- color.loci[i,"col"]
      pcl[chr==this.chr & pos>=pos.from & pos<=pos.to] <- TRUE
    }
  }
  
  locY = pmin(locY,cap)
  if(flip){locY <- -locY}
  locY= locY+shift

  if( add ){
    points(locX,locY,pch=plot.sym,col=curcol,cex=0.8, ...)
  }
  else{
    plot(locX,locY,pch=plot.sym,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(-obsmax/20,obsmax),cex=0.8, ...)
    axis(2,las=1)
    mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
  }

  points(locX[pcl],locY[pcl],pch=plot.sym[pcl],col=curcol[pcl],cex=0.8, ...)
  
  if(NROW(color.loci)>0){
    
  }

  put.text.at <- -obsmax/20
  if( chr.labels ){
    for (i in 1:21)
      {
        labpos = (x[i+1] + x[i]) / 2
        text(labpos,put.text.at,i,cex=0.8)
      }
    labpos = (x[22]+((x[22]-x[21])/2))
    text(labpos,put.text.at,"22",cex=0.8)
    labpos = (x[23]+750000)            #Midpoint of the X chromosome
    text(labpos,put.text.at,"X",cex=0.8)
  }
}

MH.over.under <- function( data.up, data.down, caption.up, caption.down, cap = 100, top.hits = 0, up.hits = 0, down.hits=0, close = FALSE, highlight = c(), color.loci=data.frame(), range.shift=2, ... ){
  range = min( ceiling( -log10( min( c( data.up$PVAL, data.down$PVAL ) ) ) ), cap )+range.shift 
  ylimit = c(-range-1,range)
  xlimit = c(min(c(data.up$POS, data.down$POS)), 32000000 )
  laba = c(abs((-range):0),0:range)
  labpos=(-range-1):range
  plot(100, 100, col = "white", xlim=xlimit, ylim=ylimit, yaxt = "n", xaxt = "n", bty = "n", xlab = "", ylab = "-log10(P-value)", ... )
  axis(2, at = labpos, labels = laba)
  mtext(caption.down, 1 )
  mtext(caption.up, 3 )
  message( "Plotting upper dataset" )
  MH.plot(data.up, flip=FALSE, add=TRUE, chr.labels=TRUE, cap=cap, color.loci=color.loci, ... )
  message( "Plotting lower dataset" )
  MH.plot(data.down, flip=TRUE, add = TRUE, chr.labels=FALSE, shift = -1, cap=cap, color.loci=color.loci,... )

##   if( top.hits!=0 & up.hits==0 & down.hits==0 ){
##     up.hits=top.hits
##     down.hits=top.hits
##   }

  if((typeof(up.hits)!="list")){
    up.hits <- top.n.hits(data.up, n=top.hits, close=close)
  }

  if( dim(up.hits)[1] > 0 ){
    print( up.hits )
    add.top.hits(data.up, up.hits, cap=cap, highlight=highlight)
  }

  if(typeof(down.hits)!="list"){       #Not data frame of hits
    down.hits <- top.n.hits(data.down, n=top.hits, close=close)
  }

  if( dim(down.hits)[1] > 0 ){
    print( down.hits )
    add.top.hits(data.down, down.hits, flip=TRUE, shift=-1, cap=cap, highlight=highlight)
  }


}

## Add a list of top hits to a manhatten plot
add.top.hits <- function(data, top.hits, flip=FALSE, shift=0, cap = 100, highlight = c()){
  x <- get.chr.starts(data$CHR, data$POS)

  text.text <- c()
  text.height <- c()
  text.xpos <- c()
  text.textcex <- c()
  text.col <- c()
  
  counter = 0
  for( i in 1:(dim(top.hits)[1])){
    xpos = x[top.hits[i,"CHR"]]+top.hits[i,"POS"]/100
    height = min(-log10(top.hits[i,"PVAL"]), cap )
    if( top.hits[i,"SNP"] %in% highlight ){
      textcol="darkgreen"
      textcex=0.75
    }
    else{
      textcol="black"
      textcex=0.75
    }
##     lines(c(xpos,xpos),c(0+shift, height * (1-2*flip)+shift), col="red", lty=3)

    text.text <- c(text.text, top.hits[i,"SNP"])
    ## Hacks
    text.height <- c(text.height, height* (1-2*flip)+shift)

    text.xpos <- c( text.xpos, xpos)
    text.textcex <- c(text.textcex, textcex)
    text.col <- c(text.col, textcol)
  }

  text.df <- data.frame(text=as.character(text.text), height=text.height, xpos=text.xpos, cex=text.textcex, col=as.character(text.col), stringsAsFactors=FALSE )
  text.df <- text.df[order(text.df$xpos),]
  text.df$xpos.new <- spread(text.df$xpos, 5e5 )

  for( i in 1:nrow(text.df) ){
    if(text.df$text[i]=="rs5015480"){text.df$xpos.new[i]=text.df$xpos.new[i]+1000000}
    if(text.df$text[i]=="rs11196175"){text.df$xpos.new[i]=text.df$xpos.new[i]+1000000}
  }  
  
  ##   text( text.df$xpos.new, 11*(1-2*flip)+shift, text.df$text, cex=text.df$cex, col=text.df$col, srt=90 )
  text( text.df$xpos.new, pmin(pmax(text.df$height + 3*(1-2*flip),-12.5),11) , text.df$text, cex=text.df$cex, col=text.df$col, srt=90 )

  print(text.df)
  
  for( i in 1:nrow(text.df) ){
    lines(c(text.df[i,"xpos"],text.df[i,"xpos.new"]), c(text.df[i,"height"],min(max(text.df[i,"height"]+2.3*(1-1.15*flip)+shift,-10.8),10.3) + shift), col="red", lty=3)
  }   
}

QQ.plot <- function(data, add=FALSE, pt.col="black", ...){
  data<- data[data$PVAL>-1,]
  obspval <- sort(data$PVAL)
  logobspval <- -(log10(obspval))
  exppval <- c(1:length(obspval))
  logexppval <- -(log10( (exppval-0.5)/length(exppval)))
  obsmax <- trunc(max(logobspval))+1
  expmax <- trunc(max(logexppval))+1
  if(!add){
  plot(c(0,expmax), c(0,expmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l", ...)
}
  points(logexppval, logobspval, pch=23, cex=.4, bg=pt.col)

}

genomic.inflation <- function(pvals){
  return( median(qchisq(pvals, 1, lower.tail=FALSE)) / qchisq(0.5,1) )
}

draw.qq.ci <- function(N, fill.col="lightgrey", border.col="grey"){
    lo <- -log10(sapply(1:N,function(x) qbeta(.025,x,N-x+1)))
    up <- -log10(sapply(1:N,function(x) qbeta(.975,x,N-x+1)))
    e = -log10(((1:N)-0.5)/N)
    polygon(c(e,rev(e)),c(lo,rev(up)),col=fill.col,border=border.col)


}

qqPlotOfPValues <- function( pvals, gi = NULL, add=FALSE, BW=FALSE, jitter=FALSE, col="blue", linecol="red", ci=TRUE,  ... ){
##   p <- sort(unique(pvals))
  p <- sort(pvals)
  e = -log10(((1:length(p))-0.5)/length( p ))
  o = -log10(p)

  if(jitter){o <-  jitter(o)}

  if(add){
    points( e, o, pch= 20, bty = "n", cex=0.3, col=col, ... )
  }
  else{

    plot( e, o, xlab = "Expected -log10(p-values)", ylab = "Observed -log10(p-values)", pch= 20, bty = "n", cex=0.3, type="n", ... )
    if(!BW & ci){
        draw.qq.ci(length(p))
    }
    abline(0,1,col=linecol)
    points( e, o, pch= 20, cex=0.3, col=col, ... )
  }

  if(BW){
    lines(e,lo,lty=2)
    lines(e,up,lty=2)
    abline(0,1,col="black")
  }

  if( !is.null(gi)){
    medx = median(e)
    medy = median(o)
    offset = 0.05 * max(e)
    lines(c(medx,medx-0.5*offset),c(medy,medy+offset),col="red",lty=3)
    text(medx-0.5*offset,medy+offset+0.1,substitute(lambda == gi, list(gi=round(gi,2)) ) )
  }
}

## A function by Greg Snow from http://www.mail-archive.com/r-help@r-project.org/msg15983.html
## Spreads out the labels so they are reasonably spread out. 
spread <- function(x, mindiff) {
  df <- x[-1] - x[-length(x)]
  i <- 1
  while (any(df < mindiff)) {
    x[c(df < mindiff, FALSE)] <- x[c(df < mindiff, FALSE)] - mindiff/10
    x[c(FALSE, df < mindiff)] <- x[c(FALSE, df < mindiff)] + mindiff/10
    df <- x[-1] - x[-length(x)]
    i <- i + 1
    if (i > 100) {
      break
    }
  }
  x
}

## Get the start positions of each chromosome, for plotting purposes
get.chr.starts <- function(chr,pos){
  x <- 1:24
  for (i in 1:24)
    {
      if (!any(chr==i)){
        x[i] <- 0
      }
      else{
        curchr=which(chr==i)
        x[i] <- trunc((max(pos[curchr]))/100)
      }
    }
  for (i in 2:24)
    {
      x[i] <- x[i] + x[i-1] + 100
    }
  x <- c(0,x)
  return( x )
}
