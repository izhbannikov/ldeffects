library(shiny)
library(ggplot2)
library(grid)
library(extrafont)
source("multiplot.R")
source("functions.R")
source("utils.R")
source("plotting.R")

oldState<-NULL
newState<-NULL
t0 <- NULL
pvv1 <- NULL
pvv2 <- NULL
m00 <- NULL
m10 <- NULL
m01 <- NULL

# Define server logic #
shinyServer(function(input, output, session) {
  
  observe({
    newState<<-getState(input)
    #cat(which(oldState-newState != 0),"\n")
    if(length(which(oldState-newState != 0)) == 1) {
      i<-which(oldState-newState != 0)[1]
      
      if(!is.na(i)){
        if(input$constrained == F) {
          rem <- 1-newState[i]
          a <- sum(newState[-i])
          
          if(a==0) newState[-i] <<- rem/length(newState[-i])
          else newState[-i] <<- rem*(newState[-i]/a)
        
          for(j in 1:length(newState))
            if(j!=i)
              updateSliderInput(session, paste0("slider", j), value=newState[j])
        }
    }
    
    #---------------#
  }
  oldState<<-newState
  
    
  })
  
  dataset <- reactive({
    saved <- loadData()
    len.sav <- dim(saved)[1]
    if(len.sav) {
      updateSliderInput(session, "time", value=c(saved[len.sav,1],saved[len.sav,2]))
      updateSliderInput(session, "slider1", value=saved[len.sav,3])
      updateSliderInput(session, "slider2", value=saved[len.sav,4])
      updateSliderInput(session, "slider3", value=saved[len.sav,5])
      updateSliderInput(session, "slider4", value=saved[len.sav,6])
      updateSliderInput(session, "mu00", value=saved[len.sav,7])
      updateSliderInput(session, "H1", value=saved[len.sav,8])
      updateSliderInput(session, "H2", value=saved[len.sav,9])
      updateSliderInput(session, "D1", value=saved[len.sav,10])
      updateSliderInput(session, "D2", value=saved[len.sav,11])
      updateSliderInput(session, "epistasis", value=saved[len.sav,12])
      updateSliderInput(session, "c", value=saved[len.sav,13])
      updateSliderInput(session, "R1", value=saved[len.sav,14])
      updateSliderInput(session, "R2", value=saved[len.sav,15])
    }
  })
  
  data <- reactive({
    
    vals<-getState(input)
    m00 <<- vals[1]
    m01 <<- vals[2]
    m10 <<- vals[3]
    m11 <<- vals[4]
    
    if(input$gomp_mu00 == TRUE) {
      
      amu00 <- input$a_mu00
      bmu00 <- input$b_mu00
    
      t1 <- input$time[1]
      t2 <- input$time[2]
      if(input$constrained==F) {
        pars <- list(m00=m00, m01=m01, m10=m10, m11=m11, amu00=amu00, bmu00=bmu00, t1=t1, t2=t2, dcase=input$dcase, D1=input$D1, D2=input$D2, H1=input$H1, H2=input$H2, epistasis=input$epistasis, c=input$c, R1=input$R1, R2=input$R2)
        dd <- calc_gompertz(pars)
        pvv1 <<- m10 + m11
        pvv2 <<- m01 + m11
      } else {
        m00 <<- 1 + m11 - pvv1 - pvv2
        m10 <<- pvv1 - m11
        m01 <<- pvv2 - m11
        pars <- list(m00=m00, m01=m01, m10=m10, m11=m11, amu00=amu00, bmu00=bmu00, t1=t1, t2=t2, dcase=input$dcase, D1=input$D1, D2=input$D2, H1=input$H1, H2=input$H2, epistasis=input$epistasis, c=input$c, R1=input$R1, R2=input$R2)
        dd <- calc_gompertz(pars)
      }
      
    } else {
      mu00 <- input$mu00
      
      if(input$dcase == T) {
        mu10 <- mu00*(1+input$D1)
        mu01 <- mu00*(1+input$D2)
        mu11 <- mu00*(1+input$D1 + input$D2)
      } else {
        mu10 <- mu00*input$H1
        mu01 <- mu00*input$H2
        mu11 <- mu00*input$H1*input$H2
      }
      
      if(input$epistasis) {
        mu10 <- mu00 + input$R2
        mu01 <- mu00 + input$R1
        mu11 <- mu00 + input$R1 + input$R2 + input$c*input$R1*input$R2
      }
    
      t1 <- input$time[1]
      t2 <- input$time[2]
      #-----------------#
    
      if(input$constrained==F) {
        pars <- list(m00=m00, m01=m01, m10=m10, m11=m11, mu10=mu10, mu00=mu00, mu11=mu11, mu01=mu01, t1=t1, t2=t2)
        dd <- calc(pars)
        pvv1 <<- m10 + m11
        pvv2 <<- m01 + m11
      } else {
        m00 <<- 1 + m11 - pvv1 - pvv2
        m10 <<- pvv1 - m11
        m01 <<- pvv2 - m11
        pars <- list(m00=m00, m01=m01, m10=m10, m11=m11, mu10=mu10, mu00=mu00, mu11=mu11, mu01=mu01, t1=t1, t2=t2)
        dd <- calc(pars)
      }
    
      dd$mu00 <- mu00
      dd$mu10 <- mu10
      dd$mu01 <- mu01
      dd$mu11 <- mu11
    }
    dd[["pars"]] <- pars
    dd
  })
  
  
  mPlot <- function(cols=1, save=F){
    
    dd <- data()
    m <- dd$m
    
    pm1 <- ggplot(data=data.frame(m), aes(t)) + geom_line(aes(y = m1t,color='m1'),cex=2) + theme_bw() +
      xlab("t") + ylab("m1,m2(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
      theme(legend.justification=c(1,0), legend.position=c(1,0.5), legend.title=element_blank(), legend.text = element_text(size = 16)) + 
      geom_line(aes(y = m2t,color='m2', color='m2'),cex=2) + 
      scale_colour_manual(values=c("red","blue4"))
    
    # LD and r2:
    pld <- getPlotLD(m, t)
    ps <- getPlotPS(m, t)
    pmij <- getPlotPMIJ(m, t)
    
      
    if(save==F) {
      multiplot(pm1, pld, ps, pmij, cols=cols, 
              title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                          "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                          ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                          "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                          ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                          ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                          ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
              titlesize=12,titlefont="Courier", titleface=2)
    } else {
      if(input$notitle_main==TRUE) { 
             multiplot(pm1, pld, ps, pmij, cols=cols)
      } 
      else {
             multiplot(pm1, pld, ps, pmij, cols=cols, 
                       title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                                    "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                                    ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                                    "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                                    ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                                    ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                                    ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                       titlesize=12,titlefont="Courier", titleface=2)
      }
      
      
    }
  }
  
  muPlot <- function(save=F, cols=1) {
    dd <- data()
    t1 <- input$time[1]
    t2 <- input$time[2]
    m <- dd$m
    
    m0t <- 1 - m[,"m1t"]
    mu1 <- dd$mu10*m[,"m10"]/m[,"m1t"] + dd$mu11*m[,"m11"]/m[,"m1t"]
    mu0 <- dd$mu00*m[,"m00"]/m0t + dd$mu01*m[,"m01"]/m0t
    if(input$mu.log) {
      mu1 <- log(mu1)
      mu0 <- log(mu0)
    } 
    mu <- cbind(t=t1:t2, mu1=mu1, mu0=mu0)
    
    pmu <- ggplot(data=data.frame(mu), aes(t)) + geom_line(aes(y = mu1,color="mu1"),cex=2) + theme_bw() +
      xlab("t") + ylab("mu1(t),mu0(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
      theme(legend.justification=c(1,0), legend.position=c(1,0.5), legend.title=element_blank(), legend.text = element_text(size = 16)) + 
      scale_colour_manual(values=c("blue4","red")) +
      geom_line(aes(y = mu0, color="mu0"), linetype="dashed", cex=2)
      
    
    if(save==F) {
      pmu <- pmu + ggtitle(paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                     "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                     ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                     "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                     ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                     ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                     ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep="")) +
                  theme(plot.title = element_text(face="bold", family="Courier", size = 12))
      pmu
    } else {
      
      if(input$notitle_mortality==TRUE) { 
        pmu
      } 
      else {
        pmu <- pmu + ggtitle(paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                                   "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                                   ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                                   "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                                   ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                                   ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                                   ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep="")) +
          theme(plot.title = element_text(face="bold", family="Courier", size = 12))
      }
    }
  }
  
  muPlotHap <- function(save=F, cols=1) {
    dd <- data()
    t1 <- input$time[1]
    t2 <- input$time[2]
    m <- dd$m
    a <- input$a_mu00
    b <- input$b_mu00
    D1 <- input$D1
    D2 <- input$D2
    H1 <- input$H1
    H2 <- input$H2
    dcase <- input$dcase
    epistasis <- input$epistasis
    c <- input$c
    R1 <- input$R1
    R2 <- input$R2
    pars <- dd$pars
    pars$t <- t1:t2
    
    if(input$gomp_mu00) {
      mu00 <- mu00t(pars)
      mu10 <- mu10t(pars)
      mu01 <- mu01t(pars)
      mu11 <- mu11t(pars)
    } else {
      mu00 <- dd$mu00
      mu10 <- dd$mu10
      mu01 <- dd$mu01
      mu11 <- dd$mu11
    }
    
    if(input$mu.log) {
      mu00 <- log(mu00)
      mu10 <- log(mu10)
      mu01 <- log(mu01)
      mu11 <- log(mu11)
    } 
    
    mu <- cbind(t=t1:t2, mu00=mu00, mu10=mu10, mu01=mu01, mu11=mu11)
    
    pmu.ij=ggplot(data=data.frame(mu), aes(t)) + theme_bw() +
      geom_line(aes(y = mu00,color='mu00', size=Mij),cex=2) +
      geom_line(aes(y = mu01,color='mu01', size=Mij),cex=2) +
      geom_line(aes(y = mu10,color='mu10', size=Mij),cex=2) +
      geom_line(aes(y = mu11,color='mu11', size=Mij),cex=2) +
      xlab("t") + ylab("muij(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
      theme(legend.position=c(0.9, .5), legend.title=element_blank(), legend.text = element_text(size = 16))
    
    if(save==F) {
      pmu.ij <- pmu.ij + ggtitle(paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                                 "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                                 ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                                 "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                                 ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                                 ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                                 ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep="")) +
        theme(plot.title = element_text(face="bold", family="Courier", size = 12))
      pmu.ij
    } else {
      
      if(input$notitle_mortality_hap==TRUE) { 
        pmu.ij
      } 
      else {
        pmu.ij <- pmu.ij + ggtitle(paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                                   "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                                   ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                                   "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                                   ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                                   ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                                   ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep="")) +
          theme(plot.title = element_text(face="bold", family="Courier", size = 12))
      }
    }
  }
  
  mafPlot <- function(cols=1, save=F){
    
    dd <- data()
    m <- dd$m
    
    pmij=ggplot(data=data.frame(m), aes(t)) + theme_bw() +
      geom_line(aes(y = m00,color='m00', size=Mij),cex=2) +
      geom_line(aes(y = m01,color='m01', size=Mij),cex=2) +
      geom_line(aes(y = m10,color='m10', size=Mij),cex=2) +
      geom_line(aes(y = m11,color='m11', size=Mij),cex=2) +
      xlab("t") + ylab("mij(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
      theme(legend.position=c(0.9, .5), legend.title=element_blank(), legend.text = element_text(size = 16))
    
    if(save==F) {
      multiplot(pmij, cols=cols, 
                title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                             "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                             ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                             "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                             ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                             ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                             ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                titlesize=12,titlefont="Courier", titleface=2)
    } else {
      if(input$notitle_maf==TRUE) { 
        multiplot(pmij, cols=cols)
      } 
      else {
        multiplot(pmij, cols=cols, 
                  title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                               "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                               ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                               "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                               ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                               ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                               ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                  titlesize=12,titlefont="Courier", titleface=2)
      }
      
      
    }
  }
  
  ldPlot <- function(cols=1, save=F){
    
    dd <- data()
    m <- dd$m
    
    # LD and r2:
    
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
    
    if(save==F) {
      multiplot(pld, cols=cols, 
                title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                             "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                             ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                             "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                             ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                             ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                             ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                titlesize=12,titlefont="Courier", titleface=2)
    } else {
      if(input$notitle_ld==TRUE) { 
        multiplot(pld, cols=cols)
      } 
      else {
        multiplot(pld, cols=cols, 
                  title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                               "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                               ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                               "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                               ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                               ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                               ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                  titlesize=12,titlefont="Courier", titleface=2)
      }
    }
  }
  
  survPlot <- function(cols=1, save=F){
    
    dd <- data()
    m <- dd$m
    
    # Survival:
    ps <- ggplot(data=data.frame(m), aes(t)) + theme_bw() +
      geom_line(aes(y = S1carr,color='S1carr'),cex=2) +
      geom_line(aes(y = S1non, color="S1non"), linetype="dashed",cex=2) +
      xlab("t") + ylab("S(t)") + theme(axis.text=element_text(size=16), axis.title=element_text(size=22,face="bold")) +
      theme(legend.justification=c(1,0), legend.position=c(1,0.5), legend.title=element_blank(), legend.text = element_text(size = 16)) +
      scale_colour_manual(values=c("red","blue4"))
    
    if(save==F) {
      multiplot(ps, cols=cols, 
                title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                             "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                             ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                             "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                             ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                             ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                             ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                titlesize=12,titlefont="Courier", titleface=2)
    } else {
      if(input$notitle_surv==TRUE) { 
        multiplot(ps, cols=cols)
      } 
      else {
        multiplot(ps, cols=cols, 
                  title= paste("P(V1=1)=",round(pvv1,3), "; P(V2=1)=",round(pvv2,3),
                               "; m1(t0) = ", round(m[,"m1t"][1],3), "; m2(t0) = ", round(m[,"m2t"][1],3), "; LD(t0) = ", round(m[,"ld"][1],3),
                               ";\nm00(t0) = ", round(m[,"m00"][1],3), "; m01(t0) = ", round(m[,"m01"][1],3),
                               "; m10(t0) = ", round(m[,"m10"][1],3), "; m11(t0) = ", round(m[,"m11"][1],3),
                               ifelse(input$dcase==F, paste(";\nH1 =", round(input$H1,3)), paste(";\nD1 =", round(input$D1,3))), 
                               ifelse(input$dcase==F, paste("; H2 =", round(input$H2,3)), paste("; D2 =",round(input$D2,3))),
                               ifelse(input$gomp_mu00==FALSE, paste(";\nmu00 = ", round(dd$mu00,3), "; mu10 = ", round(dd$mu10,3), "; mu01 = ", round(dd$mu01,3), "; mu11 = ", round(dd$mu11,3), sep=""), ""), sep=""),
                  titlesize=12,titlefont="Courier", titleface=2)
      }
    }
  }
  
  getPVV <- function() {
    dd <- data()
    dd$pvv
  }
  
  
  output$distPlot <- renderPlot({
    dataset()
    print(mPlot())
  })
  
  output$mortalityPlot <- renderPlot({
    print(muPlot())
  })
  
  output$mortalityPlotHap <- renderPlot({
    print(muPlotHap())
  })
  
  output$mafPlot <- renderPlot({
    print(mafPlot())
  })
  
  output$ldPlot <- renderPlot({
    print(ldPlot())
  })
  
  output$survPlot <- renderPlot({
    print(survPlot())
  })
  
  output$downloadPlot <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      png(file, width = 1920, height = 1024)
      print(mPlot(2,T))
      dev.off()
    }
  )    
  
  output$downloadPlotMu <- downloadHandler(
    filename = "plot_mu.png",
    content = function(file) {
      png(file, width = 640, height = 480)
      print(muPlot(T,1))
      dev.off()
    }
  ) 
  
  output$downloadPlotMuHap <- downloadHandler(
    filename = "plot_mu_hap.png",
    content = function(file) {
      png(file, width = 640, height = 480)
      print(muPlotHap(T,1))
      dev.off()
    }
  ) 
  
  output$downloadPlotMAF <- downloadHandler(
    filename = "plot_maf.png",
    content = function(file) {
      png(file, width = 640, height = 480)
      print(mafPlot(T,1))
      dev.off()
    }
  )
  
  output$downloadPlotLD <- downloadHandler(
    filename = "plot_ld.png",
    content = function(file) {
      png(file, width = 640, height = 480)
      print(ldPlot(T,1))
      dev.off()
    }
  )
  
  output$downloadPlotSurvival <- downloadHandler(
    filename = "plot_surv.png",
    content = function(file) {
      png(file, width = 640, height = 480)
      print(survPlot(T,1))
      dev.off()
    }
  )
  
  session$onSessionEnded(function() {
    isolate({
      # This will get executed when a session exits
      saveData(c(input$time[1], input$time[2], 
                 input$slider1, input$slider2, input$slider3, input$slider4,
                 input$mu00, input$H1, input$H2, input$D1, input$D2))
    })
  })
  
  output$downloadDocs <- downloadHandler(
       filename = function() {
         fname <- "Documentation.pdf"
         fname
       },
       content = function(file) {
         file.copy('Documentation.pdf', file)
       }
   )
  
})