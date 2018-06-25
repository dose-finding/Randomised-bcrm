
truep<-c(0.08,0.10,0.12,0.15,0.25,0.40,0.45,0.47) # Scenario 1
p.tox0<-c(0.08,0.25,0.35,0.45,0.55,0.65,0.70,0.75) # Skeleton 
target.tox<-0.25 #Target Toxicity
nsims<-10000     #Number of Simulations
Ntotal<-48       #Sample size
threshold<-0.05  #The different in the ET to be tested
ff<-"emax"       #The model to be used
f<-emax.model<-function(dose,alpha) {alpha[,1]+ (dose^alpha[,2]*alpha[,3])/(dose^alpha[,2]+alpha[,4])}   # The model

set.seed(101)
bcrmodel<-bcrm(stop=list(nmax=Ntotal),p.tox0=p.tox0
               ,ff="emax",prior.alpha=list(4,0,0),cohort=3,target.tox=target.tox, cohort.control=1,
               ,constrain=FALSE,pointest="mean",method="rjags",start=2,simulate=TRUE,nsims=nsims,truep=truep) # running the CRM with Emax model and randomization

### Saving the results (intermediate step) ###
result.tox<-mat.or.vec(nsims,length(truep))
result.notox<-mat.or.vec(nsims,length(truep))
for (i in 1:nsims){
  result.tox[i,]<-bcrmodel[[i]]$tox
  result.notox[i,]<-bcrmodel[[i]]$notox
}

#### Computing the probability of the ET being greater than alpha=0.90 ####
all.estimates<-mat.or.vec(nsims,4)
burnin.itr<-2000
production.itr<-10000
alpha.prior.plug<-c(0.08,1,0.92,1)
sdose<-find.x(ff,p.tox0,alpha=alpha.prior.plug)
prior.alpha<-list(4,0,0)
result.probability<-mat.or.vec(nsims,1)
all.prob.estimates<-mat.or.vec(nsims,length(truep))
for (t in 1:nsims){
  post<-Posterior.rjags(result.tox[t,],result.notox[t,],sdose,ff,prior.alpha,burnin.itr,production.itr)
  samples.sdose<-sapply(sdose,function(x){emax.model(x,post)})
  all.prob.estimates[t,]<-probab.means<-apply(samples.sdose,2,mean)
  rec<-which.min(probab.means-target.tox)^2
  mtd.distribution<-emax.model(rec,post)
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
