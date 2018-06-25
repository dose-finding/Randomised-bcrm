library("bcrm")
library("mvtnorm")

truep<-c(0.08,0.10,0.12,0.15,0.25,0.40,0.45,0.47) # Scenario 1
p.tox0<-c(0.08,0.25,0.35,0.45,0.55,0.65,0.70,0.75) # Skeleton 
target.tox<-0.25 #Target Toxicity
nsims<-10000     #Number of Simulations
Ntotal<-48       #Sample size
threshold<-0.05  #The different in the ET to be tested

ff<-"logit2" #The model to be used
f<-model<-function(dose,alpha){1/(1+exp(-log(alpha[,1])-alpha[,2]*dose))} # The model

### Specifying the prior ###
s1<-1.5
s2<-0.75
mu<-c(log(0.08/0.92)-s1/2,1-s2/2)
Sigma<-rbind(c(s1,0),c(0,s2))
prior.alpha<-list(4,mu,Sigma)



### Run L2R###
set.seed(101)
bcrmodel<-bcrm(stop=list(nmax=Ntotal),p.tox0=p.tox0
               ,ff="logit2",prior.alpha=list(4,mu,Sigma),cohort=3,target.tox=target.tox,cohort.control=1,
               ,constrain=FALSE,pointest="mean",method="rjags",start=1,simulate=TRUE,nsims=nsims,truep=truep)

### Saving the results (intermediate step) ###
result.tox<-mat.or.vec(nsims,length(truep))
result.notox<-mat.or.vec(nsims,length(truep))
for (i in 1:nsims){
  result.tox[i,]<-bcrmodel[[i]]$tox
  result.notox[i,]<-bcrmodel[[i]]$notox
}

#### Computing the probability of the ET being greater than alpha=0.90 ####
all.estimates<-mat.or.vec(nsims,2)
burnin.itr<-2000
production.itr<-10000
alpha.prior.plug<-exp(prior.alpha[[2]]+diag(prior.alpha[[3]])/2)
sdose<-find.x(ff,p.tox0,alpha=alpha.prior.plug)
result.probability<-mat.or.vec(nsims,1)
all.prob.estimates<-mat.or.vec(nsims,length(truep))
for (t in 1:nsims){
  post<-Posterior.rjags(result.tox[t,],result.notox[t,],sdose,ff,prior.alpha,burnin.itr,production.itr)
  samples.sdose<-sapply(sdose,function(x){model(x,post)})
  all.prob.estimates[t,]<-probab.means<-apply(samples.sdose,2,mean)
  rec<-which.min(probab.means-target.tox)^2
  mtd.distribution<-model(rec,post)
  control.distribution<-post[,1]+threshold
  final.prob<-mean(mtd.distribution>control.distribution)
  result.probability[t]<-sum(final.prob)
  cat(t,"out of",nsims,"\n")
}
###########

print.bcrm.sim(bcrmodel) # printing the results

mean(result.probability>0.9) # the probability to detect the ET

set.seed(101)
p.opt<-benchmark(P=truep,N=Ntotal,target=target.tox,sim=nsims) # computing the estimates for the benchmark

error<-mat.or.vec(nsims,1)
for (i in 1:nsims){
  error[i]<-sqrt(sum((all.prob.estimates[i,]-truep)^2)/sum((p.opt[i,]-truep)^2))
}
NMSE<-mean(error)
NMSE  # printing NMSE

