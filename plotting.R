# Plotting funtions #

###### For mPlot #####

getPlotPM1 <- function(m, t) {
    
    pm1 <- ggplot(data=data.frame(m), aes(t)) + 
           geom_line(aes(y = m1t,color='m1'),cex=2) + 
           theme_bw() +
           xlab("t") + ylab("m1,m2(t)") + 
           theme(axis.text=element_text(size=16), 
                 axis.title=element_text(size=22,face="bold")) +
           theme(legend.justification=c(1,0), 
                 legend.position=c(1,0.5), 
                 legend.title=element_blank(), 
                 legend.text = element_text(size = 16)) + 
           geom_line(aes(y = m2t,color='m2', color='m2'),cex=2) + 
           scale_colour_manual(values=c("red","blue4"))
    
    return(pm1)
}


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

getTitleMPlot <- function(pvv1, pvv2, m, input, dd) {
  
    title <- paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
               "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
               ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
               "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
               ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
               ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
               ifelse(input$gomp_mu00==FALSE, 
                      paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), 
               ifelse(input$epistasis==T, paste("\nc = ",round(input$c,3), "; R1 = ", round(input$R1,3), "; R2 = ", round(input$R2,3), sep=""), ""),
               sep="")

    return(title)  
}

####### muPlot ###########
getPlotMu <- function(mu,t) {
    pmu <- ggplot(data=data.frame(mu), aes(t)) + geom_line(aes(y = mu1,color="mu1"),cex=2) + theme_bw() +
           xlab("t") + ylab("mu1(t),mu0(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
           theme(legend.justification=c(1,0), legend.position=c(1,0.5), legend.title=element_blank(), legend.text = element_text(size = 16)) + 
           scale_colour_manual(values=c("blue4","red")) +
           geom_line(aes(y = mu0, color="mu0"), linetype="dashed", cex=2)
  return(pmu)
}

getTitleMuPlot <- function(pvv1, pvv2, m, input, dd) {
  
  title <- paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                 "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                 ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                 "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                 ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                 ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                 ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""),
                 ifelse(input$epistasis==T, paste("\nc = ",round(input$c,3), "; R1 = ", round(input$R1,3), "; R2 = ", round(input$R2,3), sep=""), ""),
                 sep="")
           
  return(title)  
}

######## muPlotHap ########
getPlotMuIJ <- function(mu, t) {
    pmu.ij=ggplot(data=data.frame(mu), aes(t)) + theme_bw() +
           geom_line(aes(y = mu00,color='mu00', size=Mij),cex=2) +
           geom_line(aes(y = mu01,color='mu01', size=Mij),cex=2) +
           geom_line(aes(y = mu10,color='mu10', size=Mij),cex=2) +
           geom_line(aes(y = mu11,color='mu11', size=Mij),cex=2) +
           xlab("t") + ylab("muij(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
           theme(legend.position=c(0.9, .5), legend.title=element_blank(), legend.text = element_text(size = 16))
    return(pmu.ij)
}

####### mafPlot ########
getPlotPMIJ2 <- function(m, t) {
    pmij <- ggplot(data=data.frame(m), aes(t)) + theme_bw() +
            geom_line(aes(y = m00,color='m00', size=Mij),cex=2) +
            geom_line(aes(y = m01,color='m01', size=Mij),cex=2) +
            geom_line(aes(y = m10,color='m10', size=Mij),cex=2) +
            geom_line(aes(y = m11,color='m11', size=Mij),cex=2) +
            xlab("t") + ylab("mij(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
            theme(legend.position=c(0.9, .5), legend.title=element_blank(), legend.text = element_text(size = 16))
    
    return(pmij)
}

getMafPlotTitle <- function(pvv1, pvv2, m, input, dd) {
    title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
               "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
               ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
               "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
               ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
               ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
               ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), 
               ifelse(input$epistasis==T, paste("\nc = ",round(input$c,3), "; R1 = ", round(input$R1,3), "; R2 = ", round(input$R2,3), sep=""), ""),
               sep="")
    return(title)
}

###### ldPlot ######
getLDPlot <- function(m,t) {
  pld <- ggplot(data=data.frame(m), aes(t)) + 
    geom_line(aes(y = ld,color='D'),cex=2) + 
    geom_line(aes(y = r2,colour = 'r2'),cex=2) + 
    theme_bw() +
    xlab("t") + ylab("D(t), r2(t)") + 
    theme(axis.text=element_text(size=16), 
          axis.title=element_text(size=22,face="bold")) +
    theme(legend.justification=c(1,0), 
          legend.position=c(1, .5), 
          legend.title=element_blank(), legend.text = element_text(size = 16))
  return(pld)
}