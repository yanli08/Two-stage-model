# library(rjags)
library(R2jags)

###################### input data
# replace 0 with NA to avoid JAGS error-node inconsis ent with paraents
#######################
sim.yr=1995:2016 #modeled year time serires
n.yr=length(sim.yr) #number of modeled years

# catch data
dd.catch=read.table('input_catch_data.txt',header=T,sep='\t',row.names=NULL) #106 fish
dd.catch=dd.catch[which(dd.catch$Year %in% sim.yr),] #pick up time series of catch data that will be used in the model
dd.catch=dd.catch[,-1] #remove 1st col year

# abundance index data
dd.cpue=read.table('input_cpue.txt',header=T,sep='\t',row.names=NULL) #standardized, fish per tow
dd.cpue=dd.cpue[which(dd.cpue$Year %in% sim.yr),] #pick up time series of index data that will be used in the model
dd.cpue=dd.cpue[,-1] #remove 1st col year
colnum.sep=c(3,6,8,10,12,14) #col number in the index data that is surveyed in september, except spawner cpue because sep spCPUE_y calibrate sp_y
dd.cpue[2:n.yr,colnum.sep]=dd.cpue[1:(n.yr-1),colnum.sep] #sep cpue_y calibrate N_y+1 assuming most harvest occur May-Oct;spwanning occur september; computationally easier to modify here than the jags model file; when plotting, need to match back to original cpue data

# leave.out=12 #col number of indices to leave out, SEAMAP spawnerout
# dd.cpue=dd.cpue[,-(leave.out)]

##################### input parameters
n.catch=ncol(dd.catch) #number of catch data/coln in catch data
n.cpue=ncol(dd.cpue) #number of index/coln in cpue data
colnum.RM=1:3 #col number in the index data that is for recurit male 
colnum.RF=4:6 #col number in the index data that is for recurit female
colnum.NM=7:10 #col number in the index data that is for fully recruited male
colnum.NF=11:14 #col number in the index data that is for fully recruited female
colnum.SP=15:16 #col number in the index data that is for spawner

sex.ratioF=0.5 #proportion of female in recruits sex raio
seleN=1 #selectivirty for fully recruited male and female
matR=0.044 #proportion of mature female among female recruits
matN=0.9 #proportion of mature female among female fully recriuted
spM=0 #fraction for spawner mortality until spawning;1 indicate spawn at the end of time step

#################### input prior for ninital popn size,catchability and proportion of mature female in female at each stage
priorN=read.table('input_priorN.txt',header=T,sep='\t',row.names=NULL)
priorQ=read.table('input_priorQ_const.txt',header=T,sep='\t',row.names=NULL)
# priorQ=priorQ[,-(leave.out)]

#################### parameter name to track
Ry.name=rep(NA,n.yr) #total recruits, 106
for(y in 1:n.yr){
	Ry.name[y]=paste('Ry[',y,']',sep='')
}

N.name=rep(NA,n.catch*n.yr) #popn size by stage, 106
for(r in 1:n.catch){
	for(y in 1:n.yr){
		N.name[(r-1)*n.yr+y]=paste('N[',y,',',r,']',sep='')
	}
}

NySp.name=rep(NA,n.yr) #spawner
for(y in 1:n.yr){
	NySp.name[y]=paste('NySp[',y,']',sep='')
}

sele.name=rep(NA,n.catch) #selectivity, match dd.catch col
for(r in 1:n.catch){
	sele.name[r]=paste('sele[',r,']',sep='')
}

muCatch.name=rep(NA,n.catch*n.yr) #logCatch, match dd.catch col
for(r in 1:n.catch){
	for(y in 1:n.yr){
		muCatch.name[(r-1)*n.yr+y]=paste('muCatch[',y,',',r,']',sep='')
	}
}

muCPUE.name=rep(NA,n.cpue*n.yr) #logCPUE, match dd.cpue col
for(r in 1:n.cpue){
	for(y in 1:n.yr){
		muCPUE.name[(r-1)*n.yr+y]=paste('muCPUE[',y,',',r,']',sep='')
	}
}

qq.name=rep(NA,n.cpue) #catchability for index, match dd.cpue col
for(r in 1:n.cpue){
	qq.name[r]=paste('qq[',r,']',sep='')
}

Fy.name=rep(NA,n.yr) #fishing mortality F
for(y in 1:n.yr){
	Fy.name[y]=paste('Fy[',y,']',sep='')
}

sigmaCatch.name=rep(NA,n.catch) #observation error for catch data, match dd.catch col
for(r in 1:n.catch){
	sigmaCatch.name[r]=paste('sigmaCatch[',r,']',sep='')
}

sigmaCPUE.name=rep(NA,n.cpue) #observation error for index data, match dd.cpue col
for(r in 1:n.cpue){
	sigmaCPUE.name[r]=paste('sigmaCPUE[',r,']',sep='')
}

##################################### model
model=function(){
	#natrual mortality
	for(i in 1:4){
		m[i]~dlnorm(log(mbar),tau.m) %_% T(0.4,3)
	}		
	tau.m=pow(sigmaM,-2)
	
	#N matrix for popn size by stage:1st col-male recruits;2nd col-female recruits;3rd col-male fully recruit;4th col-female fully recruit	
	#recruit 106, age-1 recuitment, so one-year delay, ssn[y] gives recruit[y+1]
	for(y in 1:n.yr){
		Ry[y]=Rybar*exp(errorR[y]) #total recruits
		errorR[y]~dnorm(0,tau.R)
	}
	tau.R=pow(sigmaR,-2)
	
	for(y in 1:n.yr){
		N[y,1]=Ry[y]*(1-sex.ratioF) #male recruits
		N[y,2]=Ry[y]*sex.ratioF #female rectuis
	}
	
	#fully recruited 106	
	for(i in 3:4){
		for(y in 1:(n.yr-1)){
			N[y+1,i]=(N[y,i]*exp(-m[i]-FySele[y,i])+N[y,i-2]*exp(-m[i-2]-FySele[y,i-2]))*exp(errorN[y+1,i-2])
		}
	}
	
	for(i in 1:2){
		for(y in 1:n.yr){
			errorN[y,i]~dnorm(0,tau.N)
		}
	}
		
	tau.N=pow(sigmaN,-2)
	
	for(i in 1:n.catch){
		for(y in 1:n.yr){
			FySele[y,i]=Fy[y]*sele[i]
		}
	}
	
	for(y in 1:n.yr){
		NySp[y]=N[y,2]*matR*exp(spM*(-m[2]-FySele[y,2]))+N[y,4]*matN*exp(spM*(-m[4]-FySele[y,4])) #spawner
	}
	
	#catch data 106
	for(i in 1:n.catch){
		for(y in 1:n.yr){
			dd.catch[y,i]~dlnorm(log(muCatch[y,i]),tau.catch[i])
			muCatch[y,i]=(FySele[y,i]/(FySele[y,i]+m[i]))*(1-exp(-m[i]-FySele[y,i]))*N[y,i]+0.00001 #add a small constant if needed to avoid error during intial:inconsisten node with parents
		}
		tau.catch[i]=pow(sigmaCatch[i],-2)
	}
			
	#cpue, fish per tow
	#make a new N matrix for calculating cpue,match dd.cpue col
	for(i in colnum.RM){
		N.new[1:n.yr,i]=N[1:n.yr,1] #male recruits
	}
	for(i in colnum.RF){
		N.new[1:n.yr,i]=N[1:n.yr,2] #female recruits
	}
	for(i in colnum.NM){
		N.new[1:n.yr,i]=N[1:n.yr,3] #male fully recruited
	}
	for(i in colnum.NF){
		N.new[1:n.yr,i]=N[1:n.yr,4] #female fully recruited
	}
	for(i in colnum.SP){
		N.new[1:n.yr,i]=NySp[1:n.yr] #spawner
	}
	
	for(i in 1:n.cpue){
		for(y in 1:n.yr){
			dd.cpue[y,i]~dlnorm(log(muCPUE[y,i]),tau.cpue[i])
			muCPUE[y,i]=N.new[y,i]*qq[i]+0.00001 #add a small constant if needed to avoid error during intial:inconsisten node with parents
		}
		tau.cpue[i]=pow(sigmaCPUE[i],-2)
	}	
		
	#priors
	mbar~dunif(0.5,2)
	sigmaM~dunif(0.001,1)
	
	Rybar~dunif(priorN[1,1],priorN[2,1]) #106 fish
	sigmaR~dunif(0.001,10)
	
	N[1,3]~dunif(priorN[1,2],priorN[2,2]) #106 fish
	N[1,4]~dunif(priorN[1,3],priorN[2,3]) #106 fish
	sigmaN~dunif(0.001,10)
	
	for(y in 1:n.yr){
		Fy[y]~dunif(0.001,3)
	}
	
	for(i in 1:n.catch){
		sigmaCatch[i]~dunif(0.001,10)
	}
	
	sele[1]~dunif(0,0.6)
	sele[2]~dunif(0,0.6)
	sele[3]=seleN
	sele[4]=seleN
	
	for(i in 1:n.cpue){
		qq[i]~dunif(priorQ[1,i],priorQ[2,i])
		sigmaCPUE[i]~dunif(0.2,10)
	}	
} #end model

datafit=list('dd.catch'=dd.catch,'dd.cpue'=dd.cpue,'n.yr'=n.yr,'n.catch'=n.catch,'n.cpue'=n.cpue,'colnum.RM'=colnum.RM,'colnum.RF'=colnum.RF,'colnum.NM'=colnum.NM,'colnum.NF'=colnum.NF,'colnum.SP'=colnum.SP,'sex.ratioF'=sex.ratioF,'priorN'=priorN,'priorQ'=priorQ,'matR'=matR,'matN'=matN,'seleN'=seleN,'spM'=spM)
para=c('Rybar','sigmaR','sigmaN','mbar','sigmaM','m[1]','m[2]','m[3]','m[4]','sele[1]','sele[2]',Ry.name,N.name,NySp.name,muCatch.name,muCPUE.name,qq.name,Fy.name,sigmaCatch.name,sigmaCPUE.name)

# initial=function(){
# list('Linf'=runif(1,538,800),'k'=runif(1,0,1),'t0'=runif(1,-2,0),'sigma.age'=runif(1,0,10))
# }

# initial=list(list('Linf'=600,'k'=0.1,'t0'=-1.4,'sigma.age'=38),list('Linf'=630,'k'=0.2,'t0'=-1.5,'sigma.age'=40),list('Linf'=700,'k'=0.15,'t0'=-1.6,'sigma.age'=35))

fit=jags(data=datafit,inits=NULL,parameters.to.save=para,model.file=model,n.chains=3,n.iter=500000,n.burnin=470000,n.thin=10,DIC=T)

# #display output
# plot(fit)
# print(fit,intervals=c(0.025,0.5,0.975))
# traceplot(fit)

# ########### gelman-rubin (1992) diagnostic, potential scale reduction factor, target distribution is normal, upper CI =1 indicate converge
# fit.mcmc=as.mcmc(fit)
# gr=gelman.diag(fit.mcmc)
# write.table(gr$psrf,'gelman_psrf.txt',row.names=T,sep='\t',quote=F)
# write.table(gr$mpsrf,'gelman_mpsrf.txt',row.names=T,sep='\t',quote=F)

# xyplot(fit.mcmc)
# densityplot(fit.mcmc)

###########extract para values after burn-in and thin
output.array=fit$BUGSoutput$sims.array #[,3,448]
output=rbind(output.array[,1,],output.array[,2,],output.array[,3,]) #convert to matrix; the matrix directly from fit output (fit$BUGSoutput$sims.matrix) does not sort by chain

output=output[,c(para,'deviance')] #re-organize the parameter order in the matrix

write.table(output,'output_par.txt',row.names=F,sep='\t',quote=F)

##################trace plot and density plot
n.keep=nrow(output)/3 #num of posterior samples kept at each chain after burn-in and thin, 3 chains
# n.keep=fit$BUGSoutput$n.keep #num of posterior samples kept at each chain after burn-in and thin, 3 chains

n.col=5 #num of coln per plot
n.row=10 #num of row per plot
n.pannel=n.col*n.row #num of pannels per plot
n.plot=ceiling(ncol(output)/n.pannel)

# #trace plot
# for(hh in 1:n.plot){
	# pdf(paste('trace_',hh,'.pdf',sep=''),width=8,height=12)
	
	# if(hh<n.plot){
		# n.fig=n.pannel #number of pannels per plot
		# layout(matrix(1:n.fig,nrow=n.row,ncol=n.col,byrow=T))
		# par(mar=c(2,4,1,1),oma=c(2,2,1,1))	
	# }
	
	# if(hh==n.plot){
		# n.fig=ncol(output)-(n.plot-1)*n.pannel
		# layout(matrix(c(1:n.fig,rep(0,(n.pannel-n.fig))),nrow=n.row,ncol=n.col,byrow=T))
		# par(mar=c(2,4,1,1),oma=c(2,2,1,1))
	# }
	
	# for(pp in 1:n.fig){
		# y.min=min(output[,(hh-1)*n.pannel+pp],na.rm=T)
		# y.max=max(output[,(hh-1)*n.pannel+pp],na.rm=T)
		# plot(x=1:n.keep,y=output[1:n.keep,(hh-1)*n.pannel+pp],type='l',lwd=1,col='red',ylim=c(y.min,y.max),xlab='',ylab=colnames(output)[(hh-1)*n.pannel+pp],main=NULL)
		# lines(x=1:n.keep,y=output[(n.keep+1):(n.keep*2),(hh-1)*n.pannel+pp],type='l',lwd=1,col='green')
		# lines(x=1:n.keep,y=output[(2*n.keep+1):(n.keep*3),(hh-1)*n.pannel+pp],type='l',lwd=1,col='blue') #3 chains in total
	# }
	# mtext('Iteration',1,outer=T)

# dev.off()
# }

#density plot
for(hh in 1:n.plot){
	pdf(paste('density_',hh,'.pdf',sep=''),width=8,height=12)
	
	if(hh<n.plot){
		n.fig=n.pannel #number of pannels per plot
		layout(matrix(1:n.fig,nrow=n.row,ncol=n.col,byrow=T))
		par(mar=c(4,2,1,1),oma=c(2,2,1,1))	
	}
	
	if(hh==n.plot){
		n.fig=ncol(output)-(n.plot-1)*n.pannel
		layout(matrix(c(1:n.fig,rep(0,(n.pannel-n.fig))),nrow=n.row,ncol=n.col,byrow=T))
		par(mar=c(4,2,1,1),oma=c(2,2,1,1))
	}
	
	for(pp in 1:n.fig){
		den.1=density(output[1:n.keep,(hh-1)*n.pannel+pp])
		den.2=density(output[(n.keep+1):(n.keep*2),(hh-1)*n.pannel+pp])
		den.3=density(output[(2*n.keep+1):(n.keep*3),(hh-1)*n.pannel+pp])
		
		x.min=min(output[,(hh-1)*n.pannel+pp],na.rm=T)
		x.max=max(output[,(hh-1)*n.pannel+pp],na.rm=T)
		y.max=max(max(den.1$y),max(den.2$y),max(den.3$y))
		
		plot(den.1,type='l',lwd=1,col='red',xlim=c(x.min,x.max),ylim=c(0,y.max),ylab='',xlab=colnames(output)[(hh-1)*n.pannel+pp],main='')
		lines(den.2,type='l',lwd=1,col='green')
		lines(den.3,type='l',lwd=1,col='blue') #3 chains in total
	}
	mtext('Density',2,outer=T)

dev.off()
}

###############summary values for para
# output=as.matrix(read.table('output_par.txt',header=T,row.names=NULL,sep='\t',check.names=F))
####mean
para.est=matrix(NA,nrow=ncol(output)+2,ncol=4)
rownames(para.est)=c(colnames(output),'pD','DIC')
colnames(para.est)=c('mean','lower','upper','cv')

for(i in 1:ncol(output)){
para.est[i,'mean']=round(mean(output[,i],na.rm=T),6)
para.est[i,'upper']=round(quantile(output[,i],probs=0.975,na.rm=T),6)
para.est[i,'lower']=round(quantile(output[,i],probs=0.025,na.rm=T),6)
para.est[i,'cv']=round(sd(output[,i],na.rm=T)/abs(para.est[i,'mean']),2)
}

para.est['pD',1]=round(fit$BUGSoutput$pD,3)
para.est['DIC',1]=round(fit$BUGSoutput$DIC,3)

write.table(para.est,'summary_par_mean.txt',row.names=T,sep='\t',quote=F) #NOTE:the catch and cpue value in summary file is not log-scale anymore
