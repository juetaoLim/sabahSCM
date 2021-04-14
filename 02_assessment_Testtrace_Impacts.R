rm(list=ls())
library(readxl)
library(reshape2)
df <- read_excel("~/nBox/NatExpSabah/data/testing/testing.xlsx")

#DiD Counterfactual#
#DiD Counterfactual#
#DiD Counterfactual#
#DiD Counterfactual#

#create time trend for data
tt <- seq(1,nrow(df))
df1 <- cbind(df,tt)
df1$`Number of Tests (PCR)` <- NULL
df1$`Number of Tests (AG)` <- NULL
df1$Total <- NULL
#remove regions adjacent, high flight to sabah-population ratio
df1$Johor <- NULL
df1$KualaLumpur <- NULL
df1$Putrajaya <- NULL
df1$Selangor <- NULL
df1$Sarawak <- NULL
df1$Labuan <- NULL
df1$Kelantan <- NULL
#collapse data
df1 <- melt(df1,id.vars=c("tt","Date","Treatment"))

#create treatment only for sabah region
df1$Treatment[which(df1$variable!="Sabah")] <- 0 
df1 <- data.frame(df1,did=df1$Treatment)
#recode variables, treatment is treatment period for Sabah, DiD is treatment x T after election
ind <- which(df1$variable=="Sabah")
df1$did[ind] <- df1$did[ind] * df1$tt[ind]
glm1 <- glm(value~as.factor(variable) + tt + Treatment+did,data=df1)
summary(glm1)

#Plot Test Positive Rates, Nationally and for each region#
#Plot Test Positive Rates, Nationally and for each region#
#Plot Test Positive Rates, Nationally and for each region#
#Plot Test Positive Rates, Nationally and for each region#

pcrRatio <- df$Total/df$`Number of Tests (PCR)`
antigenRatio <- df$Total/df$`Number of Tests (AG)`
totalRatio <- df$Total/(df$`Number of Tests (AG)`+df$`Number of Tests (PCR)`)
dates <- df$Date
treatment <- df$Treatment
plotter1<- function(p=pcrRatio, #pcrratio
         a=antigenRatio, #antigen ratio
         pa=totalRatio, #total ratio
         d=dates, #dates
         treat=treatment){ #treatment
  
  
  plot(x=c(0,length(p)),
       y=c(0,max(p,a)),
       col="white",
       xlab="",
       ylab="",
       xaxt='n')
  #get very smooth lines for antigen 
  tt <- seq(1,length(a))
  ntt <- data.frame(tt=seq(1,length(a),length.out=10000))
  gamA <- gam(a~s(tt,10))
  gamA <- predict(gamA,newdata=ntt)
  lines(y=gamA,x=ntt$tt,col="darkgreen")  
  points(a,col="seagreen3",pch=16)
  
  #get very smooth lines for pcr 
  tt <- seq(1,length(p))
  ntt <- data.frame(tt=seq(1,length(p),length.out=10000))
  gamA <- gam(p~s(tt,10))
  gamA <- predict(gamA,newdata=ntt)
  lines(y=gamA,x=ntt$tt,col="blue")  
  points(p,col="navy",pch=16)
  
  #get very smooth lines for total 
  tt <- seq(1,length(pa))
  ntt <- data.frame(tt=seq(1,length(pa),length.out=10000))
  gamA <- gam(pa~s(tt,10))
  gamA <- predict(gamA,newdata=ntt)
  lines(y=gamA,x=ntt$tt,col="red4")  
  points(pa,col="pink",pch=16)
  
  #plot treatment
  abline(v=min(which(treat==1)),lty=2)
  
  #legend 
  legend(x="topleft",legend=c("Antigen",
                              "PCR",
                              "Total Test",
                              "Election"), cex=0.9,
         
         pch=c(16,16,16,NA),lty=c(NA,NA,NA,2),bty='n',col=c("darkgreen",
                                                            "navy",
                                                            "red4",
                                                            "black"))
  d <- as.Date(d)
  axis(side=1,at=seq(1,length(p),by=7),labels=d[seq(1,length(p),by=7)])
  mtext("C: Positive Percentage",adj=0,padj=-0.1, cex=0.9)
}

plotter2<- function(a=df$Total, #total case
                    b=df$Sabah, #sabah case
                    d=dates, #dates
                    treat=treatment){ #treatment
  
  
  plot(x=c(0,length(a)),
       y=c(0,max(a,b)),
       col="white",
       ylab="",
       xlab="",
       xaxt='n')
  #get very smooth lines for TOTAL 
  tt <- seq(1,length(a))
  ntt <- data.frame(tt=seq(1,length(a),length.out=10000))
  gamA <- gam(a~s(tt,20))
  gamA <- predict(gamA,newdata=ntt)
  
  xInt <- c(0,ntt$tt,max(ntt$tt),0)
  yInt <- c(0,gamA,0,0)
  polygon(x=xInt,y=yInt,col="pink",border="pink")
  
  #get very smooth lines for sabah 
  tt <- seq(1,length(b))
  ntt <- data.frame(tt=seq(1,length(b),length.out=10000))
  gamA <- gam(b~s(tt,20))
  gamA <- predict(gamA,newdata=ntt)
  xInt <- c(0,ntt$tt,max(ntt$tt),0)
  yInt <- c(0,gamA,0,0)
  polygon(x=xInt,y=yInt,col="steelblue1",border="steelblue1")
  
  
  #plot treatment
  abline(v=min(which(treat==1)),lty=2)
  
  #legend 
  legend(x="topleft",legend=c("Total",
                              "Sabah","Election"),
         cex=0.9,
         pch=c(16,16,NA),lty=c(NA,NA,2),bty='n',col=c("pink",
                                                      "steelblue1",
                                                      "black"))
  d <- as.Date(d)
  axis(side=1,at=seq(1,length(d),by=7),labels=d[seq(1,length(d),by=7)])
  mtext("A: Case Counts",adj=0,padj=-0.1, cex=0.9)
}


plotter3<- function(a=df$`Number of Tests (PCR)`, #total case
                    b=df$`Number of Tests (AG)`, #sabah case
                    d=dates, #dates
                    treat=treatment){ #treatment
  
  
  plot(x=c(0,length(a)),
       y=c(0,max(a,b)),
       col="white",
       ylab="",
       xlab="",
       xaxt='n')
  #get very smooth lines for TOTAL 
  tt <- seq(1,length(a))
  ntt <- data.frame(tt=seq(1,length(a),length.out=10000))
  gamA <- gam(a~s(tt,20))
  gamA <- predict(gamA,newdata=ntt)
  
  xInt <- c(0,ntt$tt,max(ntt$tt),0)
  yInt <- c(0,gamA,0,0)
  polygon(x=xInt,y=yInt,col="khaki",border="khaki")
  
  #get very smooth lines for sabah 
  tt <- seq(1,length(b))
  ntt <- data.frame(tt=seq(1,length(b),length.out=10000))
  gamA <- gam(b~s(tt,20))
  gamA <- predict(gamA,newdata=ntt)
  xInt <- c(0,ntt$tt,max(ntt$tt),0)
  yInt <- c(0,gamA,0,0)
  polygon(x=xInt,y=yInt,col="indianred1",border="indianred1")
  
  
  #plot treatment
  abline(v=min(which(treat==1)),lty=2)
  
  #legend 
  legend(x="topleft",legend=c("Antigen",
                              "PCR","Election"),
         cex=0.9,
         pch=c(16,16,NA),lty=c(NA,NA,2),bty='n',col=c("khaki",
                                                      "indianred1",
                                                      "black"))
  d <- as.Date(d)
  axis(side=1,at=seq(1,length(d),by=7),labels=d[seq(1,length(d),by=7)])
  mtext("B: Number of Tests",adj=0,padj=-0.1, cex=0.9)
}

dev.off()
pdf("C:\\Users\\user\\Documents\\nBox\\NatExpSabah\\plots\\testing.pdf",width=8,height=3)
par(mfrow=c(1,3),las=1,pty="s")
plotter2()
plotter3()
plotter1()
dev.off()

