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
  y[i]= ifelse(runif(1) < ratio, y.proposed, y[i-1])
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
Freq_est = matrix(NA, nrow = n_rep, ncol = length(0.5:(size.max)-0.5))
par(mfrow = c(1,1))
plot(0.5:(size.max)-0.5, Freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.1))
for(j in 1:n_rep) {
  if(j%%100==0) 
    cat(j,"...")
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
  #lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
  Freq_est[j, ] = Freq.y;
}
lines(0.5:(size.max)-0.5, Freq, type = "l",lwd = 3, col = "black")

## calculate MSE
diff_ = sweep(Freq_est, MARGIN = 2, STATS = Freq, FUN = "-")
dim(diff_)
round(cbind(diff_[1,],Freq_est[1,] - Freq),3)

## plot bias
boxplot(diff_, xlab = "length bin", ylab = expression(hat(theta) - theta), ylim = c(-0.1,0.1), main = paste0("sigma = ", sigma))

## smaller sigma =  1
sigma = 1
n_rep = 1000
Freq_est_sig_1 = matrix(NA, nrow = n_rep, ncol = length(0.5:(size.max)-0.5))
par(mfrow = c(1,1))
plot(0.5:(size.max)-0.5, Freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.1))
for(j in 1:n_rep) {
  if(j%%100==0) 
    cat(j,"...")
  y.all=as.numeric()
  for(sims in 1:n.sims) {
    if(sims%%1000==0) 
      cat(sims,"...")
    y[1]=sample(Agents,1)
    for(i in 2:n) {
      y.proposed=round(rnorm(1,y[i-1],sigma))
      if(y.proposed < 1 | y.proposed > size.max) 
        y.proposed=y[i-1]
      ratio=Freq[y.proposed]/Freq[y[i-1]]
      y[i]=ifelse(runif(1)<ratio,y.proposed,y[i-1])
    }
    y.all=c(y.all,y)
  }
  Freq.y=tabulate(y.all,nbins=length(Freq))/n.tot
  #lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
  Freq_est_sig_1[j, ] = Freq.y;
}
lines(0.5:(size.max)-0.5, Freq, type = "l",lwd = 3, col = "black")

## calculate MSE
diff_ = sweep(Freq_est_sig_1, MARGIN = 2, STATS = Freq, FUN = "-")
dim(diff_)
round(cbind(diff_[1,],Freq_est[1,] - Freq),3)

boxplot(diff_, xlab = "length bin", ylab = expression(hat(theta) - theta), ylim = c(-0.1,0.1), main = paste0("sigma = ", sigma))

######################
## Do the same for AF
######################
logis_ogive = function (X, a50, a95) {
  1/(1 + 19^((a50 - X)/a95))
}

## Params
ages = 1:20
R0 = 1000
M = 0.2
a50 = 3.4
ato95 = 2.4
S_age = logis_ogive(ages, a50, ato95)

A = length(ages)
N_age = vector()
N_age[1] = R0
set.seed(21)
for(age_ndx in 2:A)
  N_age[age_ndx] = N_age[age_ndx - 1] * exp(-M) * exp(rnorm(1,0,0.5))


target_freq = N_age * S_age / sum(N_age * S_age)
barplot(height = target_freq, names = ages, main = "Target Distribution", xlab = "Ages", ylab = "Frequency")

size.max = A
## smaller sigma =  1
sigma = 0.4
n_rep = 1000
Freq_est_sig_1 = matrix(NA, nrow = n_rep, ncol = A)
par(mfrow = c(1,1))
plot(1:A, target_freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.4))
out_bounds = vector()
hist_ratio = as.numeric()
for(j in 1:n_rep) {

  y.all=as.numeric()
  out_bounds = as.numeric()
  hist_ratio = as.numeric()
  for(sims in 1:n.sims) {

    y[1]=sample(1:A, 1)
    for(i in 2:n) {
      y.proposed = round(rnorm(1, y[i-1], sigma))
      if(y.proposed < 1 | y.proposed > size.max) {
        y.proposed=y[i-1]
        out_bounds = c(out_bounds, 1)
      } else {
        out_bounds = c(out_bounds, 0)
      }
      ratio = target_freq[y.proposed] / target_freq[y[i-1]]
      hist_ratio = c(hist_ratio, ratio)
      y[i] = ifelse(runif(1) < ratio, y.proposed, y[i-1])
    }
    y.all=c(y.all,y)

  }
  if(j %% 100==0) {
    cat("\n",j,"...")
    cat(table(out_bounds))
    hist(hist_ratio)
  }
  Freq.y=tabulate(y.all,nbins=A) / (n.sims * n)
  #lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
  Freq_est_sig_1[j, ] = Freq.y;
}
table(out_bounds)
## calculate MSE
diff_ = sweep(Freq_est_sig_1, MARGIN = 2, STATS = target_freq, FUN = "-")
dim(diff_)
round(cbind(diff_[1,],Freq_est_sig_1[1,] - target_freq),3)
boxplot(diff_, xlab = "Age", ylab = expression(hat(theta) - theta), ylim = c(-0.1,0.1), main = paste0("sigma = ", sigma))

## find the breakdown in M-H algorithm what happens with a "smoother" AF controlled by deviations added to exponential decline 
set.seed(21)
for(age_ndx in 2:A)
  N_age[age_ndx] = N_age[age_ndx - 1] * exp(-M) * exp(rnorm(1,0,0.2))
target_freq = N_age * S_age / sum(N_age * S_age)
barplot(height = target_freq, names = ages, main = "Target Distribution", xlab = "Ages", ylab = "Frequency")
## smaller sigma =  1
sigma = 0.4
n_rep = 1000
Freq_est_sig_1 = matrix(NA, nrow = n_rep, ncol = A)
par(mfrow = c(1,1))
plot(1:A, target_freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.4))
out_bounds = vector()
for(j in 1:n_rep) {

  y.all=as.numeric()
  out_bounds = as.numeric()
  for(sims in 1:n.sims) {
    
    y[1]=sample(1:A, 1)
    for(i in 2:n) {
      y.proposed = round(rnorm(1, y[i-1], sigma))
      if(y.proposed < 1 | y.proposed > size.max) {
        y.proposed=y[i-1]
        out_bounds = c(out_bounds, 1)
      } else {
        out_bounds = c(out_bounds, 0)
      }
      ratio = target_freq[y.proposed] / target_freq[y[i-1]]
      y[i] = ifelse(runif(1) < ratio, y.proposed, y[i-1])
    }
    y.all=c(y.all,y)
  }
  if(j %% 100==0) {
    cat("\n",j,"...")
    cat(table(out_bounds))
    hist(hist_ratio)
  }
  Freq.y=tabulate(y.all,nbins=A) / (n.sims * n)
  #lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
  Freq_est_sig_1[j, ] = Freq.y;
}
table(out_bounds)
## calculate MSE
diff_ = sweep(Freq_est_sig_1, MARGIN = 2, STATS = target_freq, FUN = "-")
dim(diff_)
round(cbind(diff_[1,],Freq_est_sig_1[1,] - target_freq),3)
par(mfrow = c(1,1))
boxplot(diff_, xlab = "Age", ylab = expression(hat(theta) - theta), ylim = c(-0.1,0.1), main = paste0("sigma = ", sigma, " n = ", n," n.sim = ", n.sims))

## larger n = 1000
n = 1000
size.max = A
## smaller sigma =  1
n_rep = 1000
Freq_est_sig_1 = matrix(NA, nrow = n_rep, ncol = A)
par(mfrow = c(1,1))
plot(1:A, target_freq, type = "l",lwd = 2, col = "blue", cex = 1.2, xlab = "Lengths", ylab = "Proportion", ylim = c(0,0.4))
out_bounds = vector()
for(j in 1:n_rep) {
  
  y.all=as.numeric()
  out_bounds = as.numeric()
  for(sims in 1:n.sims) {
    
    y[1]=sample(1:A, 1)
    for(i in 2:n) {
      y.proposed = round(rnorm(1, y[i-1], sigma))
      ##print(paste0("y_i + 1 = ", y.proposed, " y_i = ", y[i-1], " sigma = ", sigma))
      if(y.proposed < 1 | y.proposed > size.max) {
        y.proposed=y[i-1]
        out_bounds = c(out_bounds, 1)
      } else {
        out_bounds = c(out_bounds, 0)
      }
      ratio = target_freq[y.proposed] / target_freq[y[i-1]]
      y[i] = ifelse(runif(1) < ratio, y.proposed, y[i-1])
    }
    y.all=c(y.all,y)
  }
  if(j %% 100==0) {
    cat("\n",j,"...")
    cat(table(out_bounds))
    hist(hist_ratio)
  }
  Freq.y=tabulate(y.all,nbins=A) / (n.sims * n)
  #lines(0.5:(size.max)-0.5, Freq.y, col = adjustcolor(col = "red", alpha.f = 0.3), lwd =3)
  Freq_est_sig_1[j, ] = Freq.y;
}
table(out_bounds)
## calculate MSE
diff_ = sweep(Freq_est_sig_1, MARGIN = 2, STATS = target_freq, FUN = "-")
dim(diff_)
round(cbind(diff_[1,],Freq_est_sig_1[1,] - target_freq),3)
par(mfrow = c(1,1))
boxplot(diff_, xlab = "Age", ylab = expression(hat(theta) - theta), ylim = c(-0.1,0.1), main = paste0("sigma = ", sigma, " n = ", n," n.sim = ", n.sims))


## is the small support??
## lengh bin 100, ages = 20 lets redo with 100 ages
## Params
ages = 1:100
R0 = 1000
M = 0.075
a50 = 3.4
ato95 = 2.4
S_age = logis_ogive(ages, a50, ato95)

A = length(ages)
N_age = vector()
N_age[1] = R0
set.seed(21)
for(age_ndx in 2:A)
  N_age[age_ndx] = N_age[age_ndx - 1] * exp(-M) * exp(rnorm(1,0,0.5))
















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

