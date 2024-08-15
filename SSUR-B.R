###################################
##Simulation of the variables
rm(list=objects())
library(MASS)
#library(BayesSUR,lib.loc="/data/grdegoix/Rlib")
library(tictoc)
library(BDgraph)
library(fields)
#library(fields,lib.loc="/data/grdegoix/Rlib")

#Parameters
n<-100;p<-150;s<-10
snr<-10

#BayesSUR parameters
covariancePrior<-"IW"
gammaPrior<-"hierarchical"

threshold_Gamma<-0.5

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


##############Simulation of X
mu<-c(rep(0,p))

#X
set.seed(1234)
#X<-mvrnorm(n, mu=mu,Sigma=Sigma0)
X<-mvrnorm(n, mu=mu,Sigma=Sigma)


###Construction of B
set.seed(4321)
B<-matrix(rnorm(p*s,mean = 0,sd = 1),nrow = p)


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




#We create the final B matrix xith the latent indicator matrix Gamma
for (i in 1:p){
  for (j in 1:s){
    B[i,j]<-B[i,j]*(Gamma[i,j]==1)
    
  }
}


###XB
XB<-X%*%B

###Finally, Y

###Construction of E
set.seed(1357)
E_tilde<-matrix(rnorm(n*s,mean=0,sd=1), nrow=n, ncol=s)

#M
M <- matrix(0.9, s, s)
diag(M)<-rep(1,s)

#P
P <- BDgraph::rgwish(n = 1, adj = G, b = 3, D=M) #b=3 because our G graph is decomposible with 3 components (the one on the diagonal)

var<-solve(P)



###signal-noise ratio setting to 1


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
  E <- E_tilde %*% cVar
  S<-diag(diag(t(cVar)))
  sigma<-S*S
  
  ### Sample Y
  Y <- XB + E
  
  S<-diag(diag(cVar))
  sigma<-S*S
  ### S/N Ratio
  emp_snr <- mean(diag(var(XB) %*% solve(sigma)))
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



#Get the simulation
pdf("Simulation.pdf",width=12,height=4)
bw_palette <- colorRampPalette(c("white", "black"))(256)

neg_values <- seq(min(B), 0, length.out = 128)
pos_values <- seq(0.00001, max(B), length.out = 128)
color_breaks <- c(neg_values, pos_values)
beta_palette <- colorRampPalette(c("blue", "white", "red"))(length(color_breaks) - 1)
layout(matrix(1:3, ncol = 3))

image.plot(1:p, 1:s,B,main="B",col=beta_palette,legend.width=1.5,breaks=color_breaks) 
image.plot(1:p, 1:s,Gamma,main="Gamma",col=bw_palette,legend.width=1.5) 
image.plot(1:s, 1:s,G,main="G",col=bw_palette,legend.width=1.5) 

dev.off()



##############################################################################################################
#####################################Fit a BayesSUR model#####################################################
##############################################################################################################

hyperpar <- list(a_omega=1.0,b_omega=1.0,a_eta=4,b_eta=6,a_tau=10,b_tau=10)
#Run the model
set.seed(28173)
tic("Timeofmodelfitting")
fit<-BayesSUR(Y=Y,X=X,outFilePath="SSURB",nIter=200000,
              nChains=3, burnin=100000,
              covariancePrior=covariancePrior,
              gammaPrior=gammaPrior, mrfG=mrfG,hyperpar=hyperpar,
              output_CPO = TRUE)
toc()


##############################################################################################################
#####################################Get the results#####################################################
##############################################################################################################

#MCMC diagnostics
pdf("SSURB.MCMC_diagnostic.pdf")
plot(fit, estimator = "logP", type = "diagnostics")
dev.off()


#Estimators
pdf("SSURB.pdf",width=12,height=4)
plotEstimator(fit, estimator = c("beta", "gamma","Gy"),Pmax = 0.5, beta.type = "conditional")
dev.off()


#RMSE
B.hat<-getEstimator(fit, estimator = "beta", Pmax = 0.5, beta.type = "conditional") 
RMSE<-sqrt(mean((Y-X%*%B.hat)^2))




#Confusion matrix


Gamma_hat<-getEstimator(fit, estimator = "gamma")

accuracy<-sum((Gamma_hat>=threshold_Gamma & Gamma)|(Gamma<threshold_Gamma & !Gamma ))/(p*s)
sensitivity<-sum(Gamma_hat>=threshold_Gamma & Gamma)/(sum(Gamma_hat>=threshold_Gamma & Gamma)+sum(Gamma_hat<threshold_Gamma & Gamma))
specificity<-sum(Gamma_hat<threshold_Gamma & !Gamma)/(sum(Gamma_hat<threshold_Gamma & !Gamma)+sum(Gamma_hat>=threshold_Gamma & !Gamma))


#ROC curve

library(pROC)
pdf("SSURB.ROC_curve.pdf")
roc_score = roc(as.vector(Gamma),as.vector(Gamma_hat))  # AUC score
plot.roc(roc_score, main="ROC curve ", legacy.axes=TRUE)
dev.off()



###################################################################
##############   AFFICHAGE DES RESULTATS   ###############################
######################################################################


# Charger les packages nécessaires
library("readxl")
library("writexl",lib.loc="/data/grdegoix/Rlib")
library(purrr)

# Lire le fichier Excel existant
existing_data <- read_excel("Results_SSURB.xlsx", sheet = "Sheet1")

# Afficher les données existantes
print(existing_data)

# Créer des nouvelles données à ajouter
hyperpar=flatten_dbl(hyperpar)
new_data <- data.frame(
  a_tau=hyperpar[5],
  b_tau=hyperpar[6],
  a_eta=hyperpar[3],
  b_eta=hyperpar[4],
  a_omega = hyperpar[1],
  b_omega = hyperpar[2],
  roc_score$auc,
  WAIC = elpd(fit, method="WAIC"),
  LOO = elpd(fit, method="LOO"),
  RMSE = RMSE,
  accuracy = accuracy,
  sensitivity = sensitivity,
  specificity = specificity
)

# Combiner les données existantes et nouvelles
combined_data <- rbind(existing_data, new_data)

# Sauvegarder les données combinées dans un fichier Excel
write_xlsx(combined_data, path = "Results_SSURB.xlsx")

# Message de confirmation
cat("Les données ont été ajoutées au fichier Results_SSURB.xlsx")

