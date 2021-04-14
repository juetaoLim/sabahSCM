rm(list=ls())
load("~/nBox/HY_flightdata/flightforalexlab.RData")
df <- subset(flightforalexlab,flightforalexlab$DEPCOUNTRY=="Malaysia")
df <- subset(df,df$ARRCOUNTRY=="Malaysia")
#2013 Oct Numbers
df <- subset(df,df$Month==10)
df <- subset(df,df$Year==2013)
df <- subset(df,df$ARRNAME=="Sabah")