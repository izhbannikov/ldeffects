# Plotting funtions #

getPlotLD <- function(m, t) {

    pld <- ggplot(data=data.frame(m), aes(t)) + 
           geom_line(aes(y = ld,color='D'),cex=2) + 
           geom_line(aes(y = r2,colour = 'r2'),cex=2) + 
           theme_bw() +
           xlab("t") + ylab("D(t), r2(t)") + 
           theme(axis.text=element_text(size=16), 
                 axis.title=element_text(size=22,face="bold")) +
           theme(legend.justification=c(1,0), 
                 legend.position=c(1, .5), 
                 legend.title=element_blank(), 
                 legend.text = element_text(size = 16))

    return(pld)
}


getPlotPS <- function(m, t) {
  
    ps <- ggplot(data=data.frame(m), aes(t)) + theme_bw() +
          geom_line(aes(y = S1carr,color='S1carr'),cex=2) +
          geom_line(aes(y = S1non, color="S1non"), linetype="dashed",cex=2) +
          xlab("t") + ylab("S(t)") + 
          theme(axis.text=element_text(size=16), 
                axis.title=element_text(size=22,face="bold")) +
          theme(legend.justification=c(1,0), 
                legend.position=c(1,0.5), 
                legend.title=element_blank(), 
                legend.text = element_text(size = 16)) +
          scale_colour_manual(values=c("red","blue4"))
  
    return(ps)
}

getPlotPMIJ <- function(m, t) {
    
    pmij <- ggplot(data=data.frame(m), aes(t)) + theme_bw() +
            geom_line(aes(y = m00,color='m00', size=Mij),cex=2) +
            geom_line(aes(y = m01,color='m01', size=Mij),cex=2) +
            geom_line(aes(y = m10,color='m10', size=Mij),cex=2) +
            geom_line(aes(y = m11,color='m11', size=Mij),cex=2) +
            xlab("t") + ylab("mij(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
            theme(legend.position=c(0.9, .5), legend.title=element_blank(), legend.text = element_text(size = 16))
  
    return(pmij)
}