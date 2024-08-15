####Compute RMSPE and its standard deviation

library(MASS)
library(tictoc)
library(BDgraph)


#Parameters
n<-100;p<-150;s<-10
snr<-10
a<-5
RMSP<-0
std<-0

for (i in 1:10){

threshold_Gamma<-0.5

set.seed(28173)


#Sigma
Prime <- list(
  c(1:10),c(11:15))


Sigma <- diag(p)
for (i in Prime) {
  Sigma[i, i] <- 0.1
}


for (i in 1:15){
  Sigma[i,i]=1
}



###Construction of G (dependencies about response variables)
G <- matrix(0, s, s)

Prime2 <- list(
  c(1:2),
  c(3:6),
  c(7:10)
)
Res <- Prime2

for (i in Prime2) {
  G[i, i] <- 1
}



#Prise en compte des correlations des x dans gamma

for (i in 1:length(Prime2)) {
  Sigma[1:3 +(i+3)*10,1:3 +(i+3)*10] <- 0.1
  Sigma[1:4 +(i+1)*30,1:4 +(i+1)*30] <- 0.1
}


for (i in 1:p){
  Sigma[i,i]=1
}


####Ajouter correlations entre genes qui fonctionnent ensemble mais pas pour la maladie dans Sigma
#Sigma0=Sigma
#Sigma0[146:150,146:150]<-0.1
#diag(Sigma0)<-1



####On simule alors X avec Sigma 0 mais le reste reste avce Sigma (mrfG par ex). Prior knowledge based on Sigma/Sigma0
###Ensuite, on peut supprimer des infos de Sigma (prior knowledge we don't know yet)




##############Simulation of X
mu<-c(rep(0,p))

#X
#X<-mvrnorm(n, mu=mu,Sigma=Sigma0)
set.seed(1234+a)
X.test<-mvrnorm(n, mu=mu,Sigma=Sigma)


###Construction of B
set.seed(4321)
B.test<-matrix(rnorm(p*s,mean = 0,sd = 1),nrow = p)





###Construction of Gamma
Gamma <- matrix(0, p, s)


#Correlation present in Sigma
Gamma[1:10,c(1:2,8:10)] <- 1
Gamma[11:15,6:10] <- 1

#Response variable that works together (G)
for (i in 1:length(Prime2)){
  Gamma[1:3 +(i+3)*10,Prime2[[i]]]<-1
  Gamma[1:4 +(i+1)*30,Prime2[[i]]]<-1
  
}

#Additional variable to be selected (random here), wide need to be 1 (1 x selected only)
Gamma[70,1:6]<-1
Gamma[85,1:7]<-1
Gamma[97,8:10]<-1
Gamma[133,2:4]<-1

sparsity_Gamma<- sum(Gamma)/(s*p)
#0.1093333



#We select the B coefficients (one selected) with Gamma
for (i in 1:p){
  for (j in 1:s){
    B.test[i,j]<-B.test[i,j]*(Gamma[i,j]==1)
    
  }
}


###XB
XB.test<-X.test%*%B.test

###Finally, Y

###Construction of E
set.seed(1357+a)
E_tilde.test<-matrix(rnorm(n*s,mean=0,sd=1), nrow=n, ncol=s)

#M
M <- matrix(0.9, s, s)
diag(M)<-rep(1,s)

#P
P <- BDgraph::rgwish(n = 1, adj = G, b = 3, D=M)

var<-solve(P)



###Maintenant, reglage du snr (signal-noise ratio)


factor <- 2
factor_min <- 0.001
factor_max <- 1000
count <- 0
maxit <- 10000

factor_prev <- 1


repeat{
  
  ### Sample the error 
  var <- var / factor * factor_prev
  cVar <- chol(as.matrix(var))
  E <- E_tilde.test %*% cVar
  S<-diag(diag(t(cVar)))
  sigma<-S*S
  
  ### Sample Y
  Y.test <- XB.test + E
  
  S<-diag(diag(cVar))
  sigma<-S*S
  ### S/N Ratio
  emp_snr <- mean(diag(var(XB.test) %*% solve(sigma)))
  ##############
  
  if (abs(emp_snr - snr) < (snr / 2) | count > maxit) {
    break
  } else {
    if (emp_snr < snr) { # increase factor
      factor_min <- factor
    } else { # decrease factor
      factor_max <- factor
    }
    factor_prev <- factor
    factor <- (factor_min + factor_max) / 2
  }
  count <- count + 1
}
emp_snr
count



#Erreur de test 
RMSE<-sqrt(mean((Y.test-X.test%*%B.hat)^2))
RMSP=RMSP+RMSE
a=a+5

}

RMSP<-RMSP/10
RMSP








#Standard deviation

for (i in 1:10){
  
  threshold_Gamma<-0.5
  
  set.seed(28173)
  
  
  #Sigma
  Prime <- list(
    c(1:10),c(11:15))
  
  
  Sigma <- diag(p)
  for (i in Prime) {
    Sigma[i, i] <- 0.1
  }
  
  
  for (i in 1:15){
    Sigma[i,i]=1
  }
  
  
  
  ###Construction of G (dependencies about response variables)
  G <- matrix(0, s, s)
  
  Prime2 <- list(
    c(1:2),
    c(3:6),
    c(7:10)
  )
  Res <- Prime2
  
  for (i in Prime2) {
    G[i, i] <- 1
  }
  
  
  
  #Prise en compte des correlations des x dans gamma
  
  for (i in 1:length(Prime2)) {
    Sigma[1:3 +(i+3)*10,1:3 +(i+3)*10] <- 0.1
    Sigma[1:4 +(i+1)*30,1:4 +(i+1)*30] <- 0.1
  }
  
  
  for (i in 1:p){
    Sigma[i,i]=1
  }
  
  
  ####Ajouter correlations entre genes qui fonctionnent ensemble mais pas pour la maladie dans Sigma
  #Sigma0=Sigma
  #Sigma0[146:150,146:150]<-0.1
  #diag(Sigma0)<-1
  
  
  
  ####On simule alors X avec Sigma 0 mais le reste reste avce Sigma (mrfG par ex). Prior knowledge based on Sigma/Sigma0
  ###Ensuite, on peut supprimer des infos de Sigma (prior knowledge we don't know yet)
  
  
  
  
  ##############Simulation of X
  mu<-c(rep(0,p))
  
  #X
  #X<-mvrnorm(n, mu=mu,Sigma=Sigma0)
  set.seed(1234+a)
  X.test<-mvrnorm(n, mu=mu,Sigma=Sigma)
  
  
  ###Construction of B
  set.seed(4321)
  B.test<-matrix(rnorm(p*s,mean = 0,sd = 1),nrow = p)
  
  
  
  
  
  ###Construction of Gamma
  Gamma <- matrix(0, p, s)
  
  
  #Correlation present in Sigma
  Gamma[1:10,c(1:2,8:10)] <- 1
  Gamma[11:15,6:10] <- 1
  
  #Response variable that works together (G)
  for (i in 1:length(Prime2)){
    Gamma[1:3 +(i+3)*10,Prime2[[i]]]<-1
    Gamma[1:4 +(i+1)*30,Prime2[[i]]]<-1
    
  }
  
  #Additional variable to be selected (random here), wide need to be 1 (1 x selected only)
  Gamma[70,1:6]<-1
  Gamma[85,1:7]<-1
  Gamma[97,8:10]<-1
  Gamma[133,2:4]<-1
  
  sparsity_Gamma<- sum(Gamma)/(s*p)
  #0.1093333
  
  
  
  #We select the B coefficients (one selected) with Gamma
  for (i in 1:p){
    for (j in 1:s){
      B.test[i,j]<-B.test[i,j]*(Gamma[i,j]==1)
      
    }
  }
  
  
  ###XB
  XB.test<-X.test%*%B.test
  
  ###Finally, Y
  
  ###Construction of E
  set.seed(1357+a)
  E_tilde.test<-matrix(rnorm(n*s,mean=0,sd=1), nrow=n, ncol=s)
  
  #M
  M <- matrix(0.9, s, s)
  diag(M)<-rep(1,s)
  
  #P
  P <- BDgraph::rgwish(n = 1, adj = G, b = 3, D=M)
  
  var<-solve(P)
  
  
  
  ###signal-noise ratio
  
  
  factor <- 2
  factor_min <- 0.001
  factor_max <- 1000
  count <- 0
  maxit <- 10000
  
  factor_prev <- 1
  
  
  repeat{
    
    ### Sample the error 
    var <- var / factor * factor_prev
    cVar <- chol(as.matrix(var))
    E <- E_tilde.test %*% cVar
    S<-diag(diag(t(cVar)))
    sigma<-S*S
    
    ### Sample Y
    Y.test <- XB.test + E
    
    S<-diag(diag(cVar))
    sigma<-S*S
    ### S/N Ratio
    emp_snr <- mean(diag(var(XB.test) %*% solve(sigma)))
    ##############
    
    if (abs(emp_snr - snr) < (snr / 2) | count > maxit) {
      break
    } else {
      if (emp_snr < snr) { # increase factor
        factor_min <- factor
      } else { # decrease factor
        factor_max <- factor
      }
      factor_prev <- factor
      factor <- (factor_min + factor_max) / 2
    }
    count <- count + 1
  }
  emp_snr
  count
  
  
  
  #Test error
  RMSE<-sqrt(mean((Y.test-X.test%*%B.hat)^2))
  std=std+(RMSE-RMSP)**2
  
  a=a+5
  
}

std<-sqrt(std/10)
std 


