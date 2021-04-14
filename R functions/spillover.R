#---------------------------------------------
#     Synthetic Controls with Spillovers
#---------------------------------------------
# Based on Cao & Dowd (2019)
# www.codowd.com for links to paper.
# 
# Code ported from Matlab by Connor Dowd, Original 
# code by Jianfei Cao which can be found at:
# https://voices.uchicago.edu/jianfeicao/research/
# 

#pretreatment period only!
scm_weights = function(Y,treat.col=1) {
  # Takes a Matrix of outcomes, Y,
  # with columns representing units and rows representing times
  # Where (by default) the first column is the treated unit,
  # and (by necessity), none of the rows represents treatment yet.
  # Returns weights and intercept.

  #Checking input
  if (!is.matrix(Y)) Y = as.matrix(Y)
  #Grabbing dimensions
  t = nrow(Y)
  n = ncol(Y)
  
  #Splitting into treatment, control, as well as demeaning
  Y_treated = Y[,treat.col]
  Y_untreated = Y[,-treat.col]
  
  means = colMeans(Y)
  demeaned = t(apply(Y,1,function(x) x-means))
  Y_demeaned = demeaned[,treat.col]
  X_demeaned = demeaned[,-treat.col]
  
  #Optimization
  #
  #Initial guess -- diff-in-diff weights (ish)
  b_init = rep(1,n-1)/n 
  #Criterion function
  Q = function(b) sum((Y_demeaned - (X_demeaned%*%b)[,1])^2)
  #Constraints
  A_eq = rbind(matrix(-1,ncol=n-1),diag(n-1)) #ui
  B_eq = c(-1,rep(0,n-1)) #ci
  
  opt = constrOptim(b_init,Q,grad=NULL,ui=A_eq,ci=B_eq)
  if (opt$convergence != 0) warning("convergence may be in doubt.")
  b_hat = rep(0,n)
  b_hat[-treat.col] = opt$par
  a_hat = means[treat.col]-means%*%b_hat
  return(list(a_hat=a_hat[1,1],b_hat=b_hat))
}


scm_batch_weights = function(Y) {
  #Takes a matrix Y
  #Where columns are units
  #and rows are time periods
  #Runs SCM, taking each column as the column to be predicted in turn.
  #Returns vector of intercepts and matrix of weights (0 on diag).
  n = ncol(Y)
  #Apply scm to each column.
  listed = lapply(1:n,
            function(i) scm_weights(Y,treat.col=i))
  #Unlist outputs. Make each row the weights for one unit.
  B_hat = t(sapply(listed,function(x) x$b_hat))
  A_hat = sapply(listed,function(x) x$a_hat)
  return(list(A_hat=A_hat,B_hat=B_hat))
}

sp_andrews = function(Y0,Y1,A,loo=F,C=diag(nrow(A))[,1],d=0,alpha_sig=0.05) {
  #Y0 is pretreatment outcomes for all units
  #Y1 is post-treatment outcomes for all units
  #C and d define the null hypothesis: C*alpha = d
  #A is the matrix of spillovers/treatments for gamma
  #alpha_sig is the significance level.
  #loo is a boolean controlling use of the slower leave-one-out estimation, which has better finite sample properties.
  #
  n = ncol(Y0)
  t = nrow(Y0)
  
  #Estimate weights, intercepts
  batch = scm_batch_weights(Y0)
  A_hat = batch$A_hat
  B_hat = batch$B_hat
  
  #Invertibility condition:
  i_less_B = diag(n)-B_hat
  M_hat = t(i_less_B)%*%(i_less_B)
  AMA_inv = solve(t(A)%*%M_hat%*%A)
  
  #Estimate TEs
  gamma_hat = AMA_inv%*%(t(A)%*%t(i_less_B)%*%(i_less_B%*%Y1-A_hat))
  alpha_hat = A%*%gamma_hat
  
  #Test Statistic
  P = t(C%*%alpha_hat - d)%*%(C%*%alpha_hat-d)
  G_hat = A%*%(AMA_inv%*%t(A))%*%t(i_less_B)
  
  P_t = sapply(1:t,function(j){
    if (loo) {
      newbatch = scm_batch_weights(Y0[-j,])
      a_t = newbatch$A_hat
      b_t = newbatch$B_hat
    } else {
      a_t = A_hat
      b_t = B_hat
    }
    error = Y0[j,]-(a_t+b_t%*%Y0[j,])
    u_hat = C%*%G_hat%*%error
    out = t(u_hat)%*%u_hat
    out[1,1]
  })
  pval = mean(P[1,1]<P_t)
  test = pval <= alpha_sig
  list(Estimates=gamma_hat,pval_joint = pval,test_reject = test)
}

sp_andrews_te = function(Y0,Y1,A,loo=F,Normal=F,alpha_sig=0.05) {
  #Performs estimation and inference for null: alpha_1 = 0.
  # i.e. fixing C, d, such that inverting 
  # the test is straightforward.
  #Functionally much the same as sp_andrews
  #
  #Normal is a boolean controlling use of the normal distribution for inferential CI / pval. 
  #
  n = ncol(Y0)
  t = nrow(Y0)
  C = matrix(rep(0,n),nrow=1)
  C[,which(!!A[,1])] = 1
  d = 0
  
  #Estimate weights, intercepts
  batch = scm_batch_weights(Y0)
  A_hat = batch$A_hat
  B_hat = batch$B_hat
  
  #Invertibility condition:
  i_less_B = diag(n)-B_hat
  M_hat = t(i_less_B)%*%(i_less_B)
  AMA_inv = solve(t(A)%*%M_hat%*%A)
  
  #Estimate TEs
  gamma_hat = AMA_inv%*%(t(A)%*%t(i_less_B)%*%(i_less_B%*%Y1-A_hat))
  alpha_hat = A%*%gamma_hat
  
  #Test Statistic
  P = t(C%*%alpha_hat - d)%*%(C%*%alpha_hat-d)
  G_hat = A%*%(AMA_inv%*%t(A))%*%t(i_less_B)
  
  errors_across_t = sapply(1:t,function(j){
    if (loo) {
      newbatch = scm_batch_weights(Y0[-j,])
      a_t = newbatch$A_hat
      b_t = newbatch$B_hat
    } else {
      a_t = A_hat
      b_t = B_hat
    }
    error = Y0[j,]-(a_t+b_t%*%Y0[j,])
    u_hat = C%*%G_hat%*%error
    out = t(u_hat)%*%u_hat
    c(u_hat[1,1],out[1,1])
  })
  P_t = errors_across_t[2,]
  u_hat = errors_across_t[1,]
  
  if (Normal) {
    pval = pnorm(P[1,1],0,sd(P_t))
    pval = min(pval,1-pval)*2
    
    se = sd(u_hat)
    lb = alpha_hat[1]+qnorm(alpha_sig/2,0,se)
    ub = alpha_hat[1]+qnorm(1-alpha_sig/2,0,se)
  } else {
    pval = mean(P[1,1]<P_t)
    
    lb = alpha_hat[1]+quantile(u_hat,alpha_sig/2)
    ub = alpha_hat[1]+quantile(u_hat,1-alpha_sig/2)
  }
  test = pval <= alpha_sig
  
  return(list(Estimates=gamma_hat,pval = pval, 
              test_reject=test,CI = c(lb,ub),
              level=alpha_sig))
}

