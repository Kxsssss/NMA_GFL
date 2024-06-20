rm(list=ls())
#install.packages("R2jags")
# load packages
library(R2jags)

set.seed(42)

N.sim=1000

### helper function
# Given a group list G, return the same list with group names.
rename_G = function(G){
  if(1 %in% sapply(G, min)){
    group1_id <- which(sapply(G, min) == 1 )
    names(G)[group1_id] <- "group1"
    if(length(G) > 1) names(G)[-group1_id] <- paste0("group", 2:length(G))
  }else{
    if(length(G) > 0) names(G) <- paste0("group", 2:(length(G)+1))
  }
  singles <- setdiff(1:N.treat, unlist(G))
  if(length(singles) > 0) G <- c(G, singles = list(singles))
  G
}

#####################
model3="
model{
 
  for(i in 1:NS){
    D1[i]<-d[t2[i]]-d[t1[i]]
    prec[i]<-1/(SE[i]*SE[i])
    y[i]~dnorm(D1[i],prec[i])}
    
  for (i in 1:NT){
    for (j in i:NT){
      D[j,i]<-d[j]-d[i]}}
  for (i in 2:NT){
    dref[i]<-d[i]-d[1]}
  #TreatmeNT hierarchy
  order[1:NT]<- NT+1- rank(d[1:NT])
  for(k in 1:NT) {
    # this is when the outcome is positive - omit  'NT+1-' when the outcome is negative
    most.effective[k]<-equals(order[k],1)
    for(j in 1:NT) {
      effectiveness[k,j]<- equals(order[k],j)}}
  for(k in 1:NT) {
    for(j in 1:NT) {
      cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])}}
  #SUCRAS#
  for(k in 1:NT) {
    SUCRA[k]<- sum(cumeffectiveness[k,1:(NT-1)]) /(NT-1)}
"

group_1_bugs <- "
### group 1
for(i in 1:length(group1)){
 d[group1[i]] <- dtemp[i]
}
dtemp[1:length(group1)] ~ dmnorm(md1[1:length(group1)],prec1[1:length(group1),1:length(group1)])
for (i in 1:length(group1)) {md1[i]<-md0 }
for(k in 1:length(group1)){taud1[k,k]<-tau.sqd}
for (i in 1:(length(group1) - 1)){
  for (k in (i+1):length(group1)){
    taud1[i,k]<-0.5*tau.sqd
    taud1[k,i]<-taud1[i,k]}}
prec1[1:length(group1),1:length(group1)]<-inverse(taud1[1:length(group1),1:length(group1)])
md0~dnorm(0,0.1)
taud ~ dunif(0,2)                               
tau.sqd<- pow(taud,2)
"

other_groups_bugs <- function(id){
  gsub("%s", id, "
  ### group %s
  for(i in group%s){
    d[i]~dnorm(md%s,sd%s)} 
  md%s~dnorm(0,0.1)
  sd%s<-1/(td%s*td%s)
  td%s~dunif(0,2)
  ")
}

single_bugs <- "
  ### Singletons
  for(i in singles){
    d[i]~dnorm(0,0.0001)}
  "
ref_single_bugs <- "
  ### Singletons with ref
  d[1]<-0
  for(i in 2:length(singles)){
    d[singles[i]]~dnorm(0,0.0001)}
  "

setwd("~/Desktop/codes for examples/ex7 Parkinson's diff")
load("H2.RData")
setwd("~/Desktop/the_dark_side_of_the_force-master/simulation study 2/ExchangedGFL/H2")

N.treat = 10

for (sim.id in 1:N.sim) {
  G = group[[sim.id]]
  G = rename_G(G)
  code_bugs <- model3
  # If a group includes the reference
  if ("group1" %in% names(G)){
    code_bugs <- paste0(c(code_bugs,
                          group_1_bugs))
    # if there are singles
    if("singles" %in% names(G)){
      # If there are other groups
      if(length(G) > 2){ 
        num_other_group = 2:(length(G)-1)
        code_bugs <- paste0(c(code_bugs,
                              sapply(num_other_group, other_groups_bugs)))}
      code_bugs <- paste0(c(code_bugs, single_bugs))
    }else{# no singles
      # If there are other groups
      if(length(G) > 1){ 
        num_other_group = 2:length(G)
        code_bugs <- paste0(c(code_bugs,
                              sapply(num_other_group, other_groups_bugs)))}
    }
    
  }else{ # if the reference is a singleton
    
    # no other groups
    code_bugs <- paste0(c(code_bugs, ref_single_bugs))

    # If there are other groups as well
    if(length(G) > 1){ code_bugs <- paste0(c(model3,
                                             sapply(2:length(G), 
                                                    other_groups_bugs), 
                                                     ref_single_bugs))}
  }
  code_bugs <- paste0(c(code_bugs, "}"))  
  sink(paste("tempH2", sim.id, ".txt"))
  cat(code_bugs)
  sink()
}



################################
######## some definitions
################################
MAE=c()
bias=c()
SSresults=c()
difference.best.worst=c()
sd.best.worst=c()
D=list()
jags.m=list()
difference.best.reference=c()
sd.best.ref=c()
coverage=c()

######################

params=c() 
for (i in 1:(N.treat-1)){
  for (j in (i+1):N.treat){
    params=c(params, paste("D[",j,",",i,"]",sep=""))
  }}
for (i in 2:(N.treat)){
  params=c(params, paste("dref[",i,"]",sep=""))
}
for (i in 1:(N.treat)){
  params=c(params, paste("SUCRA[",i,"]",sep=""))
}

#number of D parameters
no.D=N.treat*(N.treat-1)/2

TE2=c(TE[9],TE[1:8])


setwd("~/Desktop/codes for examples/ex7 Parkinson's diff")
load("H2.RData")
setwd("~/Desktop/the_dark_side_of_the_force-master/simulation study 2/ExchangedGFL/H2")


for (i in 1:N.sim){
  #model1.spec<-textConnection(model1.string) 
  initialval = NULL
  model_file = paste("tempH2", i, ".txt")
  G = rename_G(group[[i]])
  
  data2 <- c(list(y = data1[[i]]$TE,SE=data1[[i]]$seTE, 
                  NS=length(data1[[i]]$studlab), t1=data1[[i]]$t1,t2=data1[[i]]$t2, 
                  NT=N.treat), G)
  
  jags.m[[i]] <- jags(data=data2,initialval,parameters.to.save = params, 
                               n.chains = 2, n.iter = 15000, n.thin=1, n.burnin = 5000, 
                               DIC=F, model.file = model_file)
  print(i)
  
  ## bias and MAE of basic parameters  
  bias[i]=(mean(jags.m[[i]]$BUGSoutput$summary[(no.D+N.treat+1):(no.D+2*N.treat-1),1]-TE))
  MAE[i]=mean(abs(jags.m[[i]]$BUGSoutput$summary[(no.D+N.treat+1):(no.D+2*N.treat-1),1]-TE))
  ### 95% CrI does not include 0
  SSresults[i]=sum(jags.m[[i]]$BUGSoutput$summary[1:no.D,3]>0|jags.m[[i]]$BUGSoutput$summary[1:no.D,7]<0) 
 
  ## best and worst treatment
  best.treat=which.max(jags.m[[i]]$BUGSoutput$summary[(no.D+1):(no.D+N.treat),1])
  best.treat=substr(names(best.treat),7,nchar(names(best.treat))-1)
  worst.treat=which.min(jags.m[[i]]$BUGSoutput$summary[(no.D+1):(no.D+N.treat),1])
  worst.treat=substr(names(worst.treat),7,nchar(names(worst.treat))-1)
  index1.difference.best.worst= (as.numeric(best.treat)>
                                   as.numeric(worst.treat))*which(rownames(jags.m[[i]]$BUGSoutput$summary)==
                                                                    paste("D[",best.treat,",",worst.treat,"]",sep=""))
  index2.difference.best.worst= (as.numeric(best.treat)<
                                   as.numeric(worst.treat))*which(rownames(jags.m[[i]]$BUGSoutput$summary)==
                                                                    paste("D[",worst.treat,",",best.treat,"]",sep=""))
  index.difference.best.worst=max(index1.difference.best.worst,index2.difference.best.worst)
  difference.best.worst[i]=abs(jags.m[[i]]$BUGSoutput$summary[index.difference.best.worst,1])
  sd.best.worst[i]=abs(jags.m[[i]]$BUGSoutput$summary[index.difference.best.worst,2])
  
  ## false positives/negatives
  D[[i]]=jags.m[[i]]$BUGSoutput$summary[1:(N.treat*(N.treat-1)/2),c(3,7)]
  
  
  index.difference.best.ref=which(rownames(jags.m[[i]]$BUGSoutput$summary)==paste("D[",best.treat,",",1,"]",sep=""))
  
  if(length(index.difference.best.ref)!=0){
    difference.best.reference[i]=abs(jags.m[[i]]$BUGSoutput$summary[index.difference.best.ref,1])
    sd.best.ref[i]=(jags.m[[i]]$BUGSoutput$summary[index.difference.best.ref,2])} else {
      difference.best.reference[i]=0
      sd.best.ref[i]=0
    }
  
  coverage[i]=(mean(jags.m[[i]]$BUGSoutput$summary[(no.D+N.treat+1):(no.D+2*N.treat-1),3]<
                      TE2&jags.m[[i]]$BUGSoutput$summary[(no.D+N.treat+1):(no.D+2*N.treat-1),7]>TE2))
  
  ## delete jags from memory
  jags.m[[i]]=NULL 
} 

save(MAE, 
     bias, 
     SSresults, 
     difference.best.worst, 
     sd.best.worst, 
     D,
     jags.m,
     difference.best.reference, 
     sd.best.ref, 
     coverage,file = "E-GFL_H2.RData")





