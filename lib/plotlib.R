## Functions for plotting alelle counts etc.
library(plotrix)

################################################################################
## Add allele counts to a plot
## Plot as individual points, or a pie chart. 
################################################################################

plot.ac <- function(a, b, cnt, pie=(a+b>20)){
    xpos=list("Germany"=10, "England"=0, "Greece"=20, "Spain"=-5, "Hungary"=20, "Luxembourg"=5, "Finland"=25, "Italy"=12, "Sweeden"=15, "Samara"=50, "Karelia"=27)
    ypos=list("Germany"=50, "England"=53, "Greece"=40, "Spain"=40, "Hungary"=45, "Luxembourg"=50, "Finland"=65, "Italy"=42, "Sweeden"=58, "Samara"=53, "Karelia"=62)
    x=xpos[[cnt]]
    y=ypos[[cnt]]
    if(pie){
        pv <- c(a,b)
        if(a==0){pv <- b}
        if(b==0){pv <- a}
        if(any(pv>0)){floating.pie(x,y, pv, col=c("blue", "red"), radius=4 ) }
    } else {
        if(a>0){
            points(x+(1:a),rep(y,a), col="blue", pch=16)}
        if(b>0){
            points(x+(1:b),rep(y-1,b), col="red", pch=16)}
    }
}
