#Generate a million agents
n.Agents=1e6
#Length distn of agents
Agents=round(c(rlnorm(0.6*n.Agents,log(25),0.1),rlnorm(0.4*n.Agents,log(40),0.15))) 
Agents=round(c(rlnorm(0.2*n.Agents,log(20),0.05),rlnorm(0.2*n.Agents,log(30),0.05),rlnorm(0.2*n.Agents,log(40),0.08),rlnorm(0.4*n.Agents,log(50),0.15))) 

Freq=tabulate(Agents)/n.Agents

size.max=max(Agents)
lengths = 0:(size.max+1)-0.5
hist(Agents,breaks=lengths)


#M-H using N(y[i-1],sigma) proposal. Trial runs:
par(mfrow=c(3,3))
n=1000 #sample size
y=rep(NA,n)
sigma=2 #The higher sigma, less the autocorrelation in samples

for(sims in 1:9) {
y[1]=sample(Agents,1)
for(i in 2:n) {
  y.proposed=round(rnorm(1,y[i-1],sigma))
  if(y.proposed<1 | y.proposed>size.max) 
    y.proposed=y[i-1]
  ratio=Freq[y.proposed]/Freq[y[i-1]]
  y[i]=
    ifelse(runif(1)<ratio,y.proposed, y[i-1])
}

hist(y)
}


#M-H using N(y[k],sigma) proposal. Check marginal over n.sims runs.
n=200 #fish per sample (A haul)
sigma=2; Sig=paste0("sigma=",sigma)
y=rep(NA,n)
n.sims=50 # number of samples Number of hauls
n.tot=n.sims*n
y.all=as.numeric()

for(sims in 1:n.sims) {
  if(sims%%1000==0) 
    cat(sims,"...")
  y[1]=sample(Agents,1)
  for(i in 2:n) {
    y.proposed=round(rnorm(1,y[i-1],sigma))
    if(y.proposed<1 | y.proposed>size.max) 
      y.proposed=y[i-1]
    ratio=Freq[y.proposed]/Freq[y[i-1]]
    y[i]=ifelse(runif(1)<ratio,y.proposed,y[i-1])
  }
  y.all=c(y.all,y)
}

par(mfrow=c(2,1))
hist(Agents,breaks=0:(size.max+1)-0.5,main=paste("Agents,",Sig))
hist(y.all,breaks=0:(size.max+1)-0.5,main=paste("Marginal distn of samples,",Sig))

Freq.y=tabulate(y.all,nbins=length(Freq))/n.tot
barplot(rbind(Freq,Freq.y),beside=T,main=Sig)


## redo the above 1000 times
n_rep = 1000
par(mfrow = c(1,1))
plot(0.5:(size.max)-0.5, Freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.1))
for(j in 1:n_rep) {
  y.all=as.numeric()
  
  for(sims in 1:n.sims) {
    if(sims%%1000==0) 
      cat(sims,"...")
    y[1]=sample(Agents,1)
    for(i in 2:n) {
      y.proposed=round(rnorm(1,y[i-1],sigma))
      if(y.proposed<1 | y.proposed>size.max) 
        y.proposed=y[i-1]
      ratio=Freq[y.proposed]/Freq[y[i-1]]
      y[i]=ifelse(runif(1)<ratio,y.proposed,y[i-1])
    }
    y.all=c(y.all,y)
  }
  Freq.y=tabulate(y.all,nbins=length(Freq))/n.tot
  lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
}
lines(0.5:(size.max)-0.5, Freq, type = "l",lwd = 3, col = "black")


n=200 #sample size
sigma=2; Sig=paste0("sigma=",sigma)
y=rep(NA,n)
n.sims=5000
n.tot=n.sims*n
y.all=as.numeric()

for(sims in 1:n.sims) {
  if(sims%%1000==0) 
    cat(sims,"...")
  y[1]=sample(Agents,1)
  for(i in 2:n) {
    y.proposed=round(rnorm(1,y[i-1],sigma))
    if(y.proposed<1 | y.proposed>size.max) 
      y.proposed=y[i-1]
    ratio=Freq[y.proposed]/Freq[y[i-1]]
    y[i]=ifelse(runif(1)<ratio,y.proposed,y[i-1])
  }
  y.all=c(y.all,y)
}

par(mfrow=c(2,1))
hist(Agents,breaks=0:(size.max+1)-0.5,main=paste("Agents,",Sig))
hist(y.all,breaks=0:(size.max+1)-0.5,main=paste("Marginal distn of samples,",Sig))

