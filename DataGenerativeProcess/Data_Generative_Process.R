##############################################################################################
# https://github.com/salleuska/IRT_raters
# ------------------------------ DATA GENERATING PROCESS -------------------------------------
rm(list=ls())
library(MASS)

N               = 100                                                           # Subjects sample size
I               = 5                                                             # Number of items
R               = 5                                                             # Number of raters
                                                                        
K               = 4                                                             # Number of categories
C               = 1                                                             # Number of students' clusters
                                                          
sp              = 3                                                             # number of subjects per rater
tot             = sp * N * I

# ------------------------------------------------------------------------------


generateData <- function( mu, sigma, c_p_s_s, N, R, tot, PPi, II, RRi, ARi, lambda, beta, delta, R_tau, R_alpha) {
  
  y              = matrix(0, nrow=tot, ncol=1)
  c_s            = matrix(0, nrow=R, ncol=1)
  eta            = matrix(0, nrow=R, ncol=1)
  pic            = matrix(0, nrow=tot, ncol=K)
  pi             = matrix(0, nrow=tot, ncol=K)
  ppi            = matrix(0, nrow=tot, ncol=K)
  
  # Subjects random allocation
  for(j in 1:N){
    c_s[j]      = sample(length(c_p_s_s), 1, prob=c_p_s_s)
    eta[j]      = rnorm(1, mu[c_s[j]], sqrt(sigma[c_s[j]]))
  }
  
  # Observations
  
  for ( n in 1:tot ) {
   
     if(ARi[n]==0){
      
      for(k in 1:K){
        pic[n,k]          <-  exp( R_alpha[RRi[n]] * lambda[II[n]] * (eta[PPi[n]] + beta[II[n]] + delta[k] + R_tau[RRi[n],1]  ) )   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
        pi[n,k]           <-  sum(pic[n,1:k])     
      }
      
    }else{
      
      for(k in 1:K){
        pic[n,k]          <-  exp( R_alpha[RRi[n]] * lambda[II[n]] * (eta[PPi[n]] + beta[II[n]] + delta[k] + R_tau[RRi[n],1] + R_tau[RRi[n],2]*eta[ARi[n]] ) )   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
        pi[n,k]           <-  sum(pic[n,1:k])     
      }
      
    }
    
    ppi[n,1:K]      <- pi[n,1:K]/sum(pi[n,1:K])
    
    y[n]   = sample(1:K,1,prob=ppi[n,1:K])
    
    }
  
  
  list(y=y, eta=eta, c_s=c_s, ppi=ppi)
}

################################################################################
################################################################################

# ------------ Subjects' Latent Ability Distribution ---------------------------
# Let's assign a distribution to the subjects' latent ability ??? 

Grid_s              = seq(-5, 5, len = 200)
mu                  = c(-1,1)
sigma               = c(0.4,0.3)
c_p_s_s             = c(0.5,0.5)
par(mfrow=c(1, 1))

# plot
curve(c_p_s_s[1]*dnorm(x,mu[1],sqrt(sigma[1])) 
      + c_p_s_s[2]*dnorm(x,mu[2],sqrt(sigma[2])),from=min(Grid_s),to=max(Grid_s),col="red",lwd=2,
      xlab=quote(eta),ylab="")

# ------------ Raters' Features Distribution ---------------------------
# Rater's features: systematic bias, autoregressive coefficient, consistency, respectively. 

Omega = rbind(                                                                   # Correlation matrix
              c(1, 0.1, 0.3),
              c(0.1, 1, 0.1),
              c(0.3, 0.1, 1))

sigma = c(1,0.1,0.1)                                                            # standard deviations vector

Sigma = diag(sigma) %*% Omega %*% diag(sigma)

mu = c(0, 1, 0)                                                                   # MVN mean 

R_features   = mvrnorm( R-1, mu, Sigma) 


R_tau        = rbind(R_features[,1:2],c(-sum(R_features[,1]), 
                                        -sum(R_features[,2]) ))                 #Constraints on one Rater's features
R_alpha      = c(exp(R_features[,3]),exp(-sum(R_features[,3])))

# R_alpha      = rep(1,I) # ???

par(mfrow=c(1, 3))
plot(R_features[,1],xlab="raters",ylab="tau",cex=2,pch=17,col="red",
     main="Raters' systematic biases",xaxt = "n")                               # raters' systematic biases
axis(1,at=1:R)
plot(R_features[,2],xlab="raters",ylab="rho",cex=2,pch=18,col="red",
     main="Raters' autoregressive parameters",xaxt = "n")                       # raters' systematic biases
axis(1,at=1:R)

plot(R_alpha,xlab="raters",ylab="alpha",cex=2,pch=8,col="red",
     main="Raters' consistency",xaxt = "n")                                     # raters' systematic biases
axis(1,at=1:R)

# plot(tau[,1],tau[,2])

# ------------ Items' Parameters Distribution ---------------------------
l_lambda    = rep(NA,I)
lambda      = rep(NA,I)
beta        = rep(NA,I)
delta       = rep(NA,K)

l_lambda    = rnorm(I-1,0,1)
beta        = rnorm(I-1,0,1)

l_lambda[I] = -sum(l_lambda[1:(I-1)])

lambda      = exp(l_lambda)                                                     # constraints on items parameters
beta[I]     = -sum(beta[1:(I-1)])


delta[2:K]  = rgamma(length(2:K),1,3)                                   
delta[1]    = 0

par(mfrow=c(1, 3))
plot(lambda,pch=18,cex=2,xlab="Items",main="Items' discimination",xaxt = "n")
axis(1,at=1:I)
plot(beta,pch=20,cex=2,xlab="Items", main="Items' easiness",xaxt = "n")
axis(1,at=1:I)
plot(delta[2:K],pch=20,cex=2,xlab="Categories",main="Categories' steps",xaxt = "n")
axis(1,at=1:(K-1),labels=2:K)
################################################################################
# Let's assign subcjets to raters 

Q_mat = matrix(0,N,R)                                                            # maps the subjetcs into the raters
gg = rep(1/R,R)
for(j in 1:N){
  if(sp == R){
    ind          = 1:R 
  }else{
    ind          = sample(1:R,sp,replace = FALSE,prob = exp(gg)) 
  }
  
  Q_mat[j,ind] = 1
  
}
# Let's Check it out!

apply(Q_mat,2,sum)
apply(Q_mat,1,sum)==sp

# Let's create the vectors mapping the observations into: the students (PP), the items (II) and the raters (RR)
# We also map each student to the previous rated by the same rater, for the AutoRegressive path (AR)
PP  = which(Q_mat[,1]==1)
RR  = rep(1,length(PP))
AR  = c(0,which(Q_mat[,1]==1)[-length(which(Q_mat[,1]==1))])

for(i in 2:I){
  PP = c(PP,which(Q_mat[,i]==1)) 
  RR = c(RR,rep(i,length(which(Q_mat[,i]==1))))
  AR = c(AR,c(0,which(Q_mat[,i]==1)[-length(which(Q_mat[,i]==1))]))
}


# Let's Check them out!
length(which(PP==1))==sp
length(which(RR==1))==apply(Q_mat,2,sum)[1]
length(AR) == length(PP)

# Let's link them to the items IDs
PPi = rep(PP[1],I)
RRi = rep(RR[1],I)
ARi = rep(AR[1],I)

for(n in 2:length(PP)){
  PPi = c(PPi, rep(PP[n],I))
  RRi = c(RRi, rep(RR[n],I))
  ARi = c(ARi, rep(AR[n],I))
}
II = rep(1:I,length(PP))

length(II)==length(PPi)
length(II)==length(ARi)
# Let's take a look
(cbind(II,PPi,RRi,ARi))[1:20,]

################################################################################



Data = generateData(mu,sigma ,c_p_s_s, N, 
                    R, length(PPi), PPi, II, RRi, ARi, 
                    lambda, beta, delta, R_tau,R_alpha) 

# Check on the generated data:
# Generated Etas
par(mfrow=c(2, 1))
curve(c_p_s_s[1]*dnorm(x,mu[1],sqrt(sigma[1])) 
      + c_p_s_s[2]*dnorm(x,mu[2],sqrt(sigma[2])),from=min(Grid_s),to=max(Grid_s),col="red",lwd=2,
      xlab=quote(eta),ylab="",main="Eta distribution")
plot(density(Data$eta),main="Empirical Eta distribution")

# Just to give an example of observed scores: 
plot(table(Data$y[which(PPi==19)]),ylab="frequency",main="Subject 19 - observed scores")
plot(table(Data$y[which(RRi==1)]),ylab="frequency",main="Rater 1 - observed scores")





