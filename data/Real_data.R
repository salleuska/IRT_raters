rm(list=ls())
library(nimble)
setwd("C:/Users/39388/Dropbox/Il mio PC (LAPTOP-NO4UO9GH)/Desktop/Bocconi/Sally")

Data  <- read.csv("real_data.csv")

P         = 30                                                                  # number of subjects

working   = which(Data[,1+2]!="NA")   
scores    = Data[working,1+2]       
II        = Data[working,1]
RRi       = Data[working,2]
PPi       = rep(1,length(working))

for(i in 2:P){
  working = which(Data[,i+2]!="NA")   
  
  scores  = c(scores, Data[working,i+2])       
  II      = c(II, Data[working,1])
  RRi     = c(RRi, Data[working,2])
  PPi     = c(PPi, rep(i,length(working)))
  
}

K   = length(table(scores))
N   = length(scores)
R   = length(table(RRi))
I   = length(table(II))

vec = list(scores,II,RRi,PPi,K,N,R,I)

write.csv(vec, file="OCSE_Long.csv")


ARi = rep(NA,N)
  
  for(i in 1:R){
    working                                   = which(RRi==i)
    names_p                                   = c("0", names(table(PPi[working])))
    
    for(j in 2:length(names_p)){
      ARi[which(PPi==names_p[j] & RRi == i)]    = names_p[j-1]
    }
  
    }

prova=cbind( PPi[which(RRi==1)], RRi[which(RRi==1)], ARi[which(RRi==1)])
prova[20:600,]

vec = list(scores,II,RRi,PPi,ARi,K,N,R,I)


write.csv(vec, file="OCSE_Long.csv")







