setwd("C:/Users/Graham/Dropbox/Manuscripts/VI Boas/MANUSCRIPT/New MS/Simulations/MetaSim")


#logistic growth
LG=function(N,r,K,T){
	plot(0,0,type="n",xlim=c(0,T),ylim=c(0,2*K),xlab="generation",ylab="N")
	for(i in 1:T){
		N=N+r*N*(K-N)
		points(i,N)
		print(N)
}}

SMM=function(n,N,theta,T,P){
	mu=theta/(4*N)
	for(i in 1:T){
		mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
		P=sample(P,2*N,TRUE)
		P=P+mutations
		}
	cbind(sample(P,n),sample(P,n))
}

# make genepop file
pop=function(L=20,n=20,N=1000,theta=4,T=10000,P=rep(50,2000)){
	X=SMM(n,N,theta,T,P)
	for(i in 2:L){
		X=cbind(X,SMM(n,N,theta,T,P))}
	X
	}

bottle=function(n=20,N=100,theta=4,T=1000,P=rep(50,200),time=rexp(1000,1/100),b=8,r=.002,K=100){
	mu=theta/(4*N)
	for(i in 1:T){
		mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
		P=sample(P,2*N,TRUE)
		P=P+mutations
		}
	N=b
	P=sample(P,2*N)
	for(i in 1:time){#if time=0, you get two generations
		N=N+r*N*(K-N)
		N=replace(N,N<1,1)
		mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
		P=sample(P,2*N,TRUE)
		P=P+mutations
		}
	if(N>n){cbind(sample(P,n),sample(P,n))}
	else{matrix(P,ncol=2)}
	}
bottle.L=function(L=20,n=20,N=100,theta=4,T=1000,P=rep(50,200),time=rexp(1000,1/100),b=8,r=.002,K=100){
	X=bottle(n,N,theta,T,P,time,b,r,K)
	for(i in 2:L){
		X=cbind(X,bottle(n,N,theta,T,P,time,b,r,K))
		}
	X
	}

# simulate a stepping stone metapopulation
smd=function(N,mu,P){#one generation of stepwise mutation (no drift)
	mutations=sample(c(-1,1),size=2*N,replace=TRUE)*sample(c(0,1),size=2*N,replace=TRUE,prob=c(1-mu,mu))
	P+mutations
	}
migrate=function(pops,N,m){#one generation of migration among pops according to migration matrix m
	#pops is a matrix of k populations and L allele frequencies
	#N causes each population's contribution of migrants to be proportional to its size
	#npops is then the expected number of gametes with each allele in each population
	k=nrow(pops)
	npops=pops
	for(i in 1:k){
		npops[i,]=t(m[i,]*N)%*%pops
		}
	npops
	}
Fst=function(x){
	# assume pops is a matrix with population membership in the first column and two columns of alleles
	# calculates pairwise Fst with no bias correction
	tab=table(pop=rep(x[,1],2),allele=c(x[,2],x[,3]))
	k=nrow(tab)
	n=rowSums(tab)
	tab=tab/n
	f=matrix(NA,k,k)
	for(i in 1:k){
		for(j in 1:k){
		ptab=tab[c(i,j),]
		Ht=1-sum(colMeans(ptab)^2)
		Hs=1-sum(colMeans(ptab^2))
		f[i,j]=(Ht-Hs)/Ht
		}
		}
	f
	}
	
MMD=function(pops,N,m,mu,T,r,K,n){
	#simulate T generations of mutation, migration, and drift in a set of populations
	#this version works for exactly 12 populations
	#pops is a table of allele frequencies
	#N is a vector of initial population sizes
	#r is the intrinsic growth rate
	#K is a vector of carrying capacities
	k=nrow(pops)
	alleles=colnames(pops)
	#plot(rep(0,k),N,xlim=c(0,T),ylim=c(0,2*max(K)),type="n",xlab="generation",ylab="N")
	for(i in 1:T){
		#migration
		popsm=migrate(pops,N,m)
		N=rowSums(popsm)
		#drift+mutation
		N=N+r*N*(K-N)
		N=round(N)
		N=replace(N,N<1,1)
		#points(rep(i,k),N)
		P=list(rep(NA,2*N[1]),rep(NA,2*N[2]),rep(NA,2*N[3]),rep(NA,2*N[4]),rep(NA,2*N[5]),rep(NA,2*N[6]),rep(NA,2*N[7]),rep(NA,2*N[8]),rep(NA,2*N[9]),rep(NA,2*N[10]),rep(NA,2*N[11]),rep(NA,2*N[12]))
		pnames=list(rep(1,2*N[1]),rep(2,2*N[2]),rep(3,2*N[3]),rep(4,2*N[4]),rep(5,2*N[5]),rep(6,2*N[6]),rep(7,2*N[7]),rep(8,2*N[8]),rep(9,2*N[9]),rep(10,2*N[10]),rep(11,2*N[11]),rep(12,2*N[12]))
		for(j in 1:k){
			P[[j]]=smd(N[j],mu,as.numeric(sample(alleles,size=2*N[j],replace=TRUE,prob=popsm[j,])))
			}
		alleles=names(table(unlist(P)))
		pops=table(unlist(pnames),unlist(P))/(2*N)	
		}
	
	p=list(rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n),rep(NA,2*n))
	pnames=list(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n),rep(11,n),rep(12,n))
	for(j in 1:k){
			p[[j]]=sample(alleles,size=2*n,replace=T,prob=pops[j,])
			}
	sam=matrix(as.numeric(unlist(p)),ncol=2,byrow=T)
	f=Fst(cbind(unlist(pnames),sam))
	list(pops=pops, sam=sam, Fst=f, N=N, pnames=unlist(pnames))
		}



# function for focal population model with extinction
FPM=function(r,K,L,m,mu,time,n){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		eqpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 goes extinct
		bpop=eqpop
		bpop[1,]=0
		N2=c(1,rep(K,11))
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r,K=K,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

# function for focal population model with bottleneck
FPMb=function(r,K,L,m,mu,time,n,b){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		bpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 drops to b individuals
		N2=c(b,rep(K,11))
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r,K=K,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

# function for 12 simultaneous bottlenecks and recovery to arbitrary size
Bottle12=function(r,K,L,m,mu,time,n,b,r2,K2){
	pops=table(1:12,rep(50,12))
	N1=rep(K,12)
	Exdat=as.vector(matrix(rep(c("A","B","C","D","E","F","G","H","I","J","K","L"),n),ncol=12,byrow=TRUE))
	pnames=Exdat
	f=matrix(0,ncol=12,nrow=12)
	for(i in 1:L){
		# step 1: simulate for 10N generations to achieve equilibrium
		bpop=MMD(pops=pops,N=N1,m=m,mu=mu,T=10*K,r=r,K=K,n=n)$pops
		# step 2: pop1 drops to b individuals
		N2=rep(b,12)
		#step 3: simulate more generations
		npop=MMD(pops=bpop,N=N2,m=m,mu=mu,T=time,r=r2,K=K2,n=n)$sam
		Exdat=cbind(Exdat,npop)
		f=f+Fst(cbind(pnames,npop))
		}
	list(data=Exdat,Fst=f/L)
	}

# simulate isolated populations	
m=0
x=c(
1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,
.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,
.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,0,
0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,0,
0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,0,
0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,0,
0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,0,
0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,0,
0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,0,
0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,.1*m/2,
.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m,.9*m/2,
.9*m/2,.1*m/2,0,0,0,0,0,0,0,.1*m/2,.9*m/2,1-m)
m=matrix(x,nrow=12,ncol=12)

# set pre and post bottleneck sizes
Ne <- 1000
Nb <- 1000
alpha <- Ne/Nb


#times <- rep(tau*Nb,8) # not clear for Cornuet and Luikart whether this should be Ne or Nb
times <- c(rep(2,20),rep(8,20),rep(16,20),rep(32,20),rep(64,20))
#write.table(data.frame(alpha,times),file="ENM1000.txt")
#write.table(data.frame(alpha,times),file="ENM1000.Fst.txt")

for(i in 1:length(times)){
	B1000=Bottle12(r=.001,K=100,L=9,m=m,mu=1/100,time=times[i],n=16,b=Nb,r2=.001,K2=Nb)
	write.table(B1000$data,file="ENM1000.txt",append=TRUE)
	write.table(B1000$Fst,file="ENM1000.Fst.txt",append=TRUE)
	}
	
