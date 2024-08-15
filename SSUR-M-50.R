###################################
##Simulation of the variables
rm(list=objects())
library(MASS)
library(BayesSUR,lib.loc="/data/grdegoix/Rlib")
library(tictoc)
library(BDgraph)
#library(fields)
library(fields,lib.loc="/data/grdegoix/Rlib")

#Parameters
n<-100;p<-150;s<-10
snr<-10
#Tester ensuite pour snr =5,3 10 pour les diff methodes

#BayesSUR parameters
covariancePrior<-"HIW"
gammaPrior<-"MRF"

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

################################## Sigma50 #####################################

Prime <- list(
  c(6:10),c(13:15))


Sigma50 <- diag(p)
for (i in Prime) {
  Sigma50[i, i] <- 0.1
}


for (i in 1:15){
  Sigma50[i,i]=1
}




#Prise en compte des correlations des x dans gamma

for (i in 1:length(Prime2)) {
  Sigma50[1:3 +(i+3)*10,1:3 +(i+3)*10] <- 0.1
  Sigma50[1:4 +(i+1)*30,1:4 +(i+1)*30] <- 0.1
}


for (i in 1:p){
  Sigma50[i,i]=1
}
############################################################################











####Ajouter correlations entre genes qui fonctionnent ensemble mais pas pour la maladie dans Sigma
#Sigma0=Sigma
#Sigma0[146:150,146:150]<-0.1
#diag(Sigma0)<-1



####On simule alors X avec Sigma 0 mais le reste reste avce Sigma (mrfG par ex). Prior knowledge based on Sigma/Sigma0
###Ensuite, on peut supprimer des infos de Sigma (prior knowledge we don't know yet)




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




######Construction de mrfG
adjacency_matrix <- function(matrix) {
  adj_matrix <- ifelse(matrix != 0, 1, 0)
  return(adj_matrix)
}
#Omega<-adjacency_matrix(Sigma)
Omega<-adjacency_matrix(solve(Sigma50))
#mrfG<-kronecker(G-diag(s),Omega) #On enlève la diagonal de s car sinon pb dans correlation sur mrfG

###Change in the way to see the problem. No more kronecker product

#########Correlations from X1 à X10
#Omega1_10
Omega1_10<-Omega
Omega1_10[-(6:10), ] <- 0  
Omega1_10[, -(6:10)] <- 0 
#A.1_10
A.1_10<-diag(s)
A.1_10[3:7,3:7]<-0
#C.1_10
comb1_2<-combn(1:2,2)
comb8_10<-combn(8:10,2)
C.1_10<-matrix(rep(0,s*s),nrow=s,ncol=s)

for (i in 1:length(comb1_2[1,])){
  C.1_10[comb1_2[1,i],comb1_2[2,i]]<-1
  C.1_10[comb1_2[2,i],comb1_2[1,i]]<-1
  
}

for (i in 1:length(comb8_10[1,])){
  C.1_10[comb8_10[1,i],comb8_10[2,i]]<-1
  C.1_10[comb8_10[2,i],comb8_10[1,i]]<-1
  
}

mrfG<-kronecker(C.1_10,Omega1_10)+
  kronecker(A.1_10,Omega1_10-diag(c(rep(0,5),rep(1,5),rep(0,p-10)))) #In the second part, we remove the diagonal because it is include b choosing d



#On rajoute les corrélations de X11 à X15
Omega11_15<-Omega
Omega11_15[-(13:15), ] <- 0  
Omega11_15[, -(13:15)] <- 0 
#A.1_10
A.11_15<-diag(s)
A.11_15[1:5,1:5]<-0
#C.11_15
comb6_10<-combn(6:10,2)
C.11_15<-matrix(rep(0,s*s),nrow=s,ncol=s)

for (i in 1:length(comb6_10[1,])){
  C.11_15[comb6_10[1,i],comb6_10[2,i]]<-1
  C.11_15[comb6_10[2,i],comb6_10[1,i]]<-1
  
}

mrfG<-mrfG+kronecker(C.11_15,Omega11_15)+
  kronecker(A.11_15,Omega11_15-diag(c(rep(0,12),rep(1,3),rep(0,p-15))))



#On rajoute les corrélations de X41 à X43 avec Y1 et Y2
Omega41_43<-Omega
Omega41_43[-43, ] <- 0  
Omega41_43[, -43] <- 0 
#A.41_43
A.41_43=diag(10)
A.41_43[3:10,3:10]<-0
#C.41_43
C.41_43<-matrix(rep(0,s*s),nrow=s,ncol=s)

for (i in 1:length(comb1_2[1,])){
  C.41_43[comb1_2[1,i],comb1_2[2,i]]<-1
  C.41_43[comb1_2[2,i],comb1_2[1,i]]<-1
  
}


mrfG<-mrfG+kronecker(C.41_43,Omega41_43)+
  kronecker(A.41_43,Omega41_43-diag(c(rep(0,42),rep(1,1),rep(0,p-43))))


#On rajoute les corrélations de X61 à X64 avec Y1 et Y2
Omega61_64<-Omega
Omega61_64[-(63:64), ] <- 0  
Omega61_64[, -(63:64)] <- 0 
#A.61_64
A.61_64=diag(10)
A.61_64[3:10,3:10]<-0
#C.61_64
C.61_64<-matrix(rep(0,s*s),nrow=s,ncol=s)

for (i in 1:length(comb1_2[1,])){
  C.61_64[comb1_2[1,i],comb1_2[2,i]]<-1
  C.61_64[comb1_2[2,i],comb1_2[1,i]]<-1
  
}
mrfG<-mrfG+kronecker(C.61_64,Omega61_64)+
  kronecker(A.61_64,Omega61_64-diag(c(rep(0,62),rep(1,2),rep(0,p-64))))


#Corrélations de X51 à X53 avec Y3,4,5,6
Omega51_53<-Omega
Omega51_53[-(52:53), ] <- 0  
Omega51_53[, -(52:53)] <- 0 
#A.51_53
A.51_53=diag(10)
A.51_53[c(1:2,7:10),c(1:2,7:10)]<-0
#C.51_53
C.51_53<-matrix(rep(0,s*s),nrow=s,ncol=s)
comb3_6<-combn(3:6,2)
for (i in 1:length(comb3_6[1,])){
  C.51_53[comb3_6[1,i],comb3_6[2,i]]<-1
  C.51_53[comb3_6[2,i],comb3_6[1,i]]<-1
  
}

mrfG<-mrfG+kronecker(C.51_53,Omega51_53)+
  kronecker(A.51_53,Omega51_53-diag(c(rep(0,51),rep(1,2),rep(0,p-53))))

#Corrélations de X91 à X94 avec Y3,4,5,6
Omega91_94<-Omega
Omega91_94[-(93:94), ] <- 0  
Omega91_94[, -(93:94)] <- 0 
#A.91_94
A.91_94=diag(10)
A.91_94[c(1:2,7:10),c(1:2,7:10)]<-0
#C.91_94
C.91_94<-matrix(rep(0,s*s),nrow=s,ncol=s)
for (i in 1:length(comb3_6[1,])){
  C.91_94[comb3_6[1,i],comb3_6[2,i]]<-1
  C.91_94[comb3_6[2,i],comb3_6[1,i]]<-1
  
}


mrfG<-mrfG+kronecker(C.91_94,Omega91_94)+
  kronecker(A.91_94,Omega91_94-diag(c(rep(0,92),rep(1,2),rep(0,p-94))))


#Corrélations de X61 à X63 avec Y7,8,9,10
Omega61_63<-Omega
Omega61_63[-63, ] <- 0  
Omega61_63[, -63] <- 0 
#A.61_63
A.61_63=diag(10)
A.61_63[1:6,1:6]<-0
#C.61_63
C.61_63<-matrix(rep(0,s*s),nrow=s,ncol=s)
comb7_10<-combn(7:10,2)
for (i in 1:length(comb7_10[1,])){
  C.61_63[comb7_10[1,i],comb7_10[2,i]]<-1
  C.61_63[comb7_10[2,i],comb7_10[1,i]]<-1
  
}


mrfG<-mrfG+kronecker(C.61_63,Omega61_63)+
  kronecker(A.61_63,Omega61_63-diag(c(rep(0,62),rep(1,1),rep(0,p-63))))


#Corrélations de X121 à X124 avec Y7,8,9,10
Omega121_124<-Omega
Omega121_124[-(123:124), ] <- 0  
Omega121_124[, -(123:124)] <- 0 
#A.61_63
A.121_124=diag(10)
A.121_124[1:6,1:6]<-0
#C.121_124
C.121_124<-matrix(rep(0,s*s),nrow=s,ncol=s)
for (i in 1:length(comb7_10[1,])){
  C.121_124[comb7_10[1,i],comb7_10[2,i]]<-1
  C.121_124[comb7_10[2,i],comb7_10[1,i]]<-1
  
}


mrfG<-mrfG+kronecker(C.121_124,Omega121_124)+
  kronecker(A.121_124,Omega121_124-diag(c(rep(0,122),rep(1,2),rep(0,p-124))))

#Variation of the matrix (else we have the true one, this is non representative of the reality)
#mrfG[,c(80:300)]<-0
#mrfG[c(1200:1400),]<-0



##################################################################################
mrfG<-which(mrfG==1, arr.ind=TRUE)



#We select the B coefficients (one selected) with Gamma
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



#Get the real ones (for comparison)
pdf("dSURM.Real.pdf",width=12,height=4)
bw_palette <- colorRampPalette(c("white", "black"))(256)
beta_palette <- colorRampPalette(c("blue", "white", "red"))(256) #needed because of the negative values
layout(matrix(1:3, ncol = 3))

image.plot(1:p, 1:s,B,main="B",col=beta_palette,legend.width=1.5) 
image.plot(1:p, 1:s,Gamma,main="Gamma",col=bw_palette,legend.width=1.5) 
image.plot(1:s, 1:s,G,main="G",col=bw_palette,legend.width=1.5) 

dev.off()



##############################################################################################################
#####################################Fit a BayesSUR model#####################################################
##############################################################################################################

hyperpar <- list(a_eta=3.5,b_eta=6.5,mrf_d=-1.8,mrf_e=10)  

#Run the model
set.seed(28173)
tic("Timeofmodelfitting")
fit<-BayesSUR(Y=Y,X=X,outFilePath="SSURM54",nIter=200000,
              nChains=3, burnin=100000,
              covariancePrior=covariancePrior,
              gammaPrior=gammaPrior, mrfG=mrfG,hyperpar=hyperpar,
              output_CPO = TRUE)
toc()


##############################################################################################################
#####################################Get the results#####################################################
##############################################################################################################

#MCMC diagnostics
pdf("SSURM50.MCMC_diagnostic.pdf")
plot(fit, estimator = "logP", type = "diagnostics")
dev.off()


#Get the Estimators
pdf("SSURM50.Estimators.pdf",width=12,height=4)
plotEstimator(fit, estimator = c("beta", "gamma","Gy"),Pmax = 0.5, beta.type = "conditional")
dev.off()


#MSE
B.hat<-getEstimator(fit, estimator = "beta", Pmax = 0.5, beta.type = "conditional") 
RMSE<-sqrt(mean((Y-X%*%B.hat)^2))




#Confusion matrix


Gamma_hat<-getEstimator(fit, estimator = "gamma")

accuracy<-sum((Gamma_hat>=threshold_Gamma & Gamma)|(Gamma<threshold_Gamma & !Gamma ))/(p*s)
sensitivity<-sum(Gamma_hat>=threshold_Gamma & Gamma)/(sum(Gamma_hat>=threshold_Gamma & Gamma)+sum(Gamma_hat<threshold_Gamma & Gamma))
specificity<-sum(Gamma_hat<threshold_Gamma & !Gamma)/(sum(Gamma_hat<threshold_Gamma & !Gamma)+sum(Gamma_hat>=threshold_Gamma & !Gamma))

#G metrics
G_hat<-getEstimator(fit, estimator = "Gy")

G.accuracy<-sum((G_hat>=threshold_Gamma & G)|(G<threshold_Gamma & !G ))/prod(dim(G))
G.sensitivity<-sum(G_hat>=threshold_Gamma & G)/(sum(G_hat>=threshold_Gamma & G)+sum(G_hat<threshold_Gamma & G))
G.specificity<-sum(G_hat<threshold_Gamma & !G)/(sum(G_hat<threshold_Gamma & !G)+sum(G_hat>=threshold_Gamma & !G))


G.sensitivity<-sum(G_hat>=threshold_Gamma & G)/sum(G)


#ROC curve

library(pROC)
pdf("SSURM50.ROC_curve.pdf")
roc_score = roc(as.vector(Gamma),as.vector(Gamma_hat))  # AUC score
plot.roc(roc_score, main="ROC curve ", legacy.axes=TRUE)
dev.off()

roc_score


###################################################################
##############   AFFICHAGE DES RESULTATS   ###############################
######################################################################


# Charger les packages nécessaires
library("readxl")
library("writexl",lib.loc="/data/grdegoix/Rlib")
library(purrr)

# Lire le fichier Excel existant
existing_data <- read_excel("Results_SSURM50.xlsx", sheet = "Sheet1")

# Afficher les données existantes
print(existing_data)

# Créer des nouvelles données à ajouter
hyperpar=flatten_dbl(hyperpar)
new_data <- data.frame(
  a_eta = hyperpar[1],
  b_eta = hyperpar[2],
  mrf_d = hyperpar[3],
  mrf_e = hyperpar[4],
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
write_xlsx(combined_data, path = "Results_SSURM50.xlsx")

# Message de confirmation
cat("Les données ont été ajoutées au fichier Results_SSURM50.xlsx")




