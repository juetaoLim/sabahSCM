rm(list=ls())
library(readxl)
cumulative <- read_excel("~/nBox/NatExpSabah/data/counts_city/raw/cumulative.xlsx")
roll <- read_excel("~/nBox/NatExpSabah/data/counts_city/raw/14dayRoll.xlsx")

name1 <- cumulative$`District/City`
name2 <- cumulative$State

cumulative$`District/City` <- NULL
cumulative$State <- NULL

counts <- apply(cumulative,MARGIN=1,diff)
counts[which(counts<0)] <- 0
counts <- t(counts)
#recover raw case count from 14 day rolling sum
roll$`District/City` <- NULL
roll$State          <- NULL
store <- list()
k<-1
for (i in 1:nrow(roll)){
  
temp1 <- round(as.numeric(counts[i,]))  #indexes district
temp2 <- round(as.numeric(roll[i,]))    #indexes district
  
out <- c(temp1,temp2)
#adjust 14 day rolling mean portion 
indAdj_start <- length(temp1) + 1
indAdj_end   <- length(out)

for (j in indAdj_start:indAdj_end){
  #take last 14 days of raw counts for first set
  ind <- c((j-13):(j-1))    #past 13 days of counts
  value <- out[j]           #current 14 day rolling
  value <- value - sum(out[ind])
  if (value<0) value <-  0  #correct data errors
  out[j] <- value
 
}
store[[k]] <- out
k <- k+1
  
}

store <- Reduce(rbind,store)
store <- data.frame(name1,name2,store)
save(store,file="~/nBox/NatExpSabah/data/counts_city/out/out.RData")#original file
