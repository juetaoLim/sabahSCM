# Data accessed 11/11/2020
# Code written 11/11/2020

# This code estimates treatment effects with spillover effects (Cao and Dowd, 2019)
rm(list=ls())
##### CHANGE HERE ##########################################################
set.path <- "~\\nBox\\NatExpSabah\\"
############################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(LowRankQP, readxl)

# load functions from Abadie & L'Hour (2019)
setwd(paste0(set.path,"\\R functions"))
source("regsynth.r")
source("wsoll1.r")
source("TZero.r")

sensMatch <- function(CUM){
# prepare data
setwd(paste0(set.path,"\\data\\counts_city\\out\\"))
store <- read_excel("out.xlsx")

start.date <- "2020-03-22"
end.date <- "2020-11-8"
analysis.start <- "2020-06-10" # start of RMCO
analysis.end <- "2020-10-12"   #travel restrictions from Sabah start
treat.start <- "2020-09-26"

dates <- seq(as.Date(start.date), as.Date(end.date), "days") # check if data is updated
names(store) <- c("City", "State", as.character(dates))
select <- (which(dates == analysis.start)+2):(which(dates == analysis.end)+2)
#select <- (which(dates == analysis.start)+2):ncol(store)
subset.store <- cbind(store[,1:2], store[,select])
if (CUM){
  cum_sum <- t(apply(store[,select],MARGIN=1,function(x) cumsum(x)))
  subset.store <- cbind(store[,1:2], cum_sum)
}
index <- subset.store$State == "Sabah"
sabah.store <- subset.store[index,]
control.store <- subset.store[!index,]

#add covariates to control dataframe
#add covariates to control dataframe
#add covariates to control dataframe
dfcov <- read_excel(paste0(set.path,"\\data\\covariates\\citylevel4.xlsx"))
dfcov <- data.frame(City=dfcov$name1,
                    Pop=dfcov$`population count`,
                    AccessM=dfcov$accessibility_mean,
                    AccessR=dfcov$accessibility_range,
                    VehicleStateA = dfcov$`State Number Vehicles on Road 2015 (Active)`,
                    VehicleStateNA = dfcov$`State Number Vehicles on Road 2015 (NonActive)`,
                    StateHosp=dfcov$`State MOH Hospital 2016`,                                    
                    StateMed=dfcov$`State MOH Specialised medical institute 2016`,
                    StateHealth1=dfcov$`State NUMBER OF ASSISTANT MEDICAL OFFICERS IN MOH 2018`,
                    StateHealth2=dfcov$`State NUMBER OF MEDICAL REHABILITATION OFFICERS IN MOH 2018`,
                    StateHealth3=dfcov$`State NUMBER OF MEDICAL OFFICERS AND SPECIALIST IN MOH 2018`,
                    Density=dfcov$`State Pop. density`,                                         
                    UrbanPop=dfcov$`State Urban pop(%)`,                                        
                    Ethnicity1=dfcov$`State Bumiputra (%)`,                                        
                    Ethnicity2=dfcov$`State Chinese (%)`)
dfcov[,-1] <- data.frame(apply(dfcov[,-1], 2, function(x) as.numeric(as.character(x))))
#randomly add in covariates
rand <- round(runif(1,min=1,max=ncol(dfcov)))
rand <- round(runif(rand,min=3,max=ncol(dfcov)))
dfcov <- dfcov[,c(1,2,rand)]
#order covariates in same way as case counts
#aggregate store and controls
sabah.store <- merge(y=sabah.store, x=dfcov, by.x="City", by.y="City", sort=F)
control.store <- merge(y=control.store, x=dfcov, by.x="City", by.y="City", sort=F)

covariate.index <- which(names(sabah.store) == analysis.start) - 1

sabah.store[,(covariate.index+1):ncol(sabah.store)] <- sabah.store[,(covariate.index+1):ncol(sabah.store)]/(sabah.store$Pop/100)
control.store[,(covariate.index+1):ncol(control.store)] <- control.store[,(covariate.index+1):ncol(control.store)]/(control.store$Pop/100)

# prepare preintervention data
select.train <- which(colnames(sabah.store) == treat.start)
X0 <- t(as.matrix(control.store[,-c(1,covariate.index, select.train:ncol(sabah.store))])) #untreated 
colnames(X0) <- control.store$City
X1 <- t(as.matrix(sabah.store[,-c(1,covariate.index, select.train:ncol(sabah.store))])) #treated 
colnames(X1) <- sabah.store$City
X <- cbind(X1,X0)

# post-intervention outcomes
Y0T <- t(as.matrix(control.store[,-c(1:(select.train-1))])) #untreated
colnames(Y0T) <- control.store$City
Y1T <- t(as.matrix(sabah.store[,-c(1:(select.train-1))]))
colnames(Y1T) <- sabah.store$City
YT <- cbind(Y1T,Y0T)

T0 <- nrow(t(as.matrix(control.store[,-c(1:covariate.index, select.train:ncol(control.store))])))
T1 <- nrow(Y0T)
C0 <- ncol(X0) # number of untreated units
C1 <- ncol(X1) # number of treated units
Tot <- ncol(X)


B <- matrix(NA, Tot, Tot)
for (i in 1:Tot) {
  temp.X0 <- X[,-i]
  pseudo.X1 <- as.matrix(X[,i])
  W <- regsynth(temp.X0, pseudo.X1, pen=0, tol=1e-6, parallel=FALSE)
  if (i == 1) {
    B[,i] <- c(0,W) 
  }  else {
    B[,i] <- append(W, 0, after = W[i-1])    
  } 
}

# Assume treatment spillovers based on flight numbers, previous election season
cov2 <- read_excel(paste0(set.path,"\\data\\covariates\\state.xlsx"))
# Second stage to estimate spillover effects
treat.spill <- diag(C1) 
#resolve flight and state information with city level
dist <- cov2$CityToSabahFlightDistanceKM
norm.dist <- (dist - min(dist))/ (max(dist) - min(dist))
geom.dist <- exp(-norm.dist)
dist.merge <- merge(x=data.frame(State=cov2$State, geom.dist), y=subset.store[,c("City","State")], by.x="State", by.y="State",sort=F)
#get appropriate flight distance to sabah from city/state
#get ordered naming for pre-treatment, non=intervention, spillover group
ind1 <- colnames(X0)
#flight data with names
ind2 <- dist.merge$City
#match ordered naming to flight data for Avec spillover matrix
order <- match(ind1,ind2)
Avec <- dist.merge$geom.dist[order]

#create spillover matrix
A2 <- cbind(matrix(0, C0, C1), Avec)
A1 <- cbind(treat.spill, rep(0, C1))
A <- rbind(A1, A2)
# create auxilary matrix M
I <- diag(Tot)
M <- t(I - B) %*% (I - B)
AMA_inv <- solve(t(A)%*%M%*%A)
gamma <- matrix(NA, dim(A)[2], T1)
for (t in 1:T1) gamma[,t] <- AMA_inv %*% t(A) %*% M %*% YT[t,]
alpha <- A %*% gamma # full effect vector on all districts

# Test for no-treatment effect, this is a joint test
# Null hypothesis: spillover for Manjung = 0 & ... & spillover for Hulu Selangor = 0
# Alt hypothesis: not null
# collect outcomes
testY0 <- t(as.matrix(control.store[,-c(1:covariate.index, ncol(control.store))])) #untreated
colnames(testY0) <- control.store$City
testY1 <- t(as.matrix(sabah.store[,-c(1:covariate.index, ncol(control.store))]))
colnames(testY1) <- sabah.store$City
testY <- cbind(testY1,testY0)

C <- cbind(matrix(0, C0, C1), diag(C0))
d <- 0
G_hat <-A %*% (AMA_inv %*% t(A)) %*% t(I - B)

p_val <- numeric(T1)

for ( t in 1:T1) {
  alpha_t <- alpha[,t] 
  P <- t(C %*% alpha_t - d) %*% (C %*% alpha_t - d)
  P_t <- numeric(T0+t-1)
  
  for (j in 1:(T0+t-1)) { 
    error <- testY[j,]-(B %*% testY[j,])
    P_t[j] <- t(C %*% G_hat %*% error) %*% C %*% G_hat %*% error
  }
  
  p_val[t] <- mean(P[1,1] < P_t)
}

# Test for individual treatment effect/spillover effect (element-wise testing of alpha), this is a single test
# Null hypothesis: effect for \alpha_i = 0 
# Alt hypothesis: \alpha_i \neq 0
# Hypothesis
path <- list()
G_hat <-A %*% (AMA_inv %*% t(A)) %*% t(I - B)
alpha_sig <- 0.05
for (i in 1:Tot) {
  
  C <- rep(0,Tot)
  C[i] <- 1
  d <- 0
  
  p_val <- numeric(T1)
  lb <- numeric(T1)
  ub <- numeric(T1)
  
  for ( t in 1:T1) {
    alpha_t <- alpha[,t] 
    P <- t(C %*% alpha_t - d) %*% (C %*% alpha_t - d)
    P_t <- numeric(T0+t-1)
    u_hat <- numeric(T0+t-1)
    
    for (j in 1:(T0+t-1)) { 
      error <- testY[j,]-(B %*% testY[j,])
      u_hat[j] <- C %*% G_hat %*% error
      P_t[j] <- t(u_hat) %*% u_hat
    }
    p_val[t] <- mean(P[1,1] < P_t)
    se <- sd(u_hat)
    lb[t] <- C %*% alpha_t + qnorm(alpha_sig/2,0,se)
    ub[t] <- C %*% alpha_t + qnorm(1-alpha_sig/2,0,se)
  }
  
  out <- rbind(ub, alpha[i,], lb, p_val)
  path[[i]] <- out  
  
}


# Test for average treatment effect for Sabah, this is a joint test
# Null hypothesis: effect for 1/24 sum_i=1^24 \alpha_i = 0 
# Alt hypothesis: not null
# Hypothesis
C <- rep(0,Tot)
C[1:C1] <- 1/C1
d <- 0
alpha_sig <- 0.05
G_hat <-A %*% (AMA_inv %*% t(A)) %*% t(I - B)

p_val <- numeric(T1)
lb <- numeric(T1)
ub <- numeric(T1)

for ( t in 1:T1) {
  alpha_t <- alpha[,t] 
  P <- t(C %*% alpha_t - d) %*% (C %*% alpha_t - d)
  P_t <- numeric(T0+t-1)
  u_hat <- numeric(T0+t-1)
  
  for (j in 1:(T0+t-1)) { 
    error <- testY[j,]-(B %*% testY[j,])
    u_hat[j] <- C %*% G_hat %*% error
    P_t[j] <- t(C %*% G_hat %*% error) %*% C %*% G_hat %*% error
  }
  
  p_val[t] <- mean(P[1,1] < P_t)
  se <- sd(u_hat)
  lb[t] <- C %*% alpha_t + qnorm(alpha_sig/2,0,se)
  ub[t] <- C %*% alpha_t + qnorm(1-alpha_sig/2,0,se)
}

out <- rbind(ub, colMeans(alpha[1:C1,]), lb, p_val)


#########get treatment effects for each region for map plotting###############
#########get treatment effects for each region for map plotting###############
#########get treatment effects for each region for map plotting###############
names <- c(sabah.store$City,control.store$City)
#cumulative treatment over treatment period for sabah per pop/100

if (!CUM){
  cumulativeTreat <- rowSums(alpha) * c(sabah.store$Pop,control.store$Pop)/100
  
  sab <- cumulativeTreat[1:length(sabah.store$Pop)]
  spillover <- cumulativeTreat[-c(1:length(sabah.store$Pop))]
}

if (CUM){
  
  cumulativeTreat <- alpha[,ncol(alpha)] * c(sabah.store$Pop,control.store$Pop)/100
  # #CIs
  sab <- cumulativeTreat[1:length(sabah.store$Pop)]
  spillover <- cumulativeTreat[-c(1:length(sabah.store$Pop))]
  
}


out1 <- list(treat=sum(sab),spill=sum(spillover),pather=path,outer=out,names=names)
return(out1)
}

storer <- list()
for (tt in 1:10000){
  try(storer[[tt]] <- sensMatch(CUM=T))
  try(print(c(storer[[tt]]$treat,storer[[tt]]$spill)))
}
save(storer,file="~\\nBox\\NatExpSabah\\submission\\placeboEvents\\out2.Rdata")
