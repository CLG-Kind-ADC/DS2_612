# hw 3
# question 2
p2l = NULL
for (i in 1:10000){
n = 3
# rank x = 1
r = 1
x = c(-1, 0, 1)
errors = rnorm(n, mean = 0, sd = 1)
y = x + errors
regr = lm(y~x-1)
betahat=regr$coefficients
eta = 0.1
z_eta=qnorm(eta)
# make new point x0, y0
x0 = 1
error0 = rnorm(1,mean=0,sd=1)
y0=error0+x0

c1 = t(x0) %*% solve(t(x) %*% x) %*% x0
C = sqrt(c1)

ncp = -z_eta/C
w = C*qt(0.95, df = n-r,ncp)

s = sqrt( sum((y-regr$fitted.values)^2) / (n-r))

L=x0 * betahat - w * s
p2l[i]=L
}
hist(p2l)
truel = x0 +qnorm(0.1)
truel
mean(p2l)

# try using a normal distribution
p2l.2 = NULL
for (i in 1:10000){
  n = 3
  # rank x = 1
  r = 1
  x = c(-1, 0, 1)
  errors = rnorm(n, mean = 0, sd = 1)
  y = x + errors
  regr = lm(y~x-1)
  betahat=regr$coefficients
  eta = 0.1
  z_eta=qnorm(eta)
  # make new point x0, y0
  x0 = 1
  error0 = rnorm(1,mean=0,sd=1)
  y0=error0+x0
  
  s = sqrt( sum((y-regr$fitted.values)^2) / (n-r))
  # normal quantile for w this time
  w = qnorm(0.95, mean=x0*betahat, sd=s)
  
  L=x0 * betahat - w * s
  p2l.2[i]=L
}
hist(p2l.2)
truel.2 = x0 +qnorm(0.1)
truel.2
mean(p2l.2)
####################################
####################################
####################################
#problem 3

cc=qchisq(0.95,df=1)
power3=1-pchisq(cc,df=1,ncp = 8)
power3

####################################
####################################
####################################
# question 4
snp=c(rep(-1,times=60),rep(0,times=80),rep(1,times=60))
var(snp)
# 0.6030151

# we have trait = snp + error
# assume independence of snp andd trait error (?)
# in the case of indep, var(trait) = var(snp) + var(error)
# and var(snp) = 0.04*var(trait)
# var(trait) = var(snp)/0.04 = 15.07538
# var(error) = 0.96*(var(snp)/0.04)
# = 14.47236
# let's just call the errors normally distributed

n = 200
r = 1
r0= 0

ncp=sum((snp)^2) / 15.07538
CC=qf(0.95,df1=r-r0, df2 = n-r)
pow0_ = 1-pf(CC,r-r0, n-r,ncp=ncp)

pvals_ = NULL
for (i in 1:1000){
  # Power: Probability(reject the null)
  trait = snp + rnorm(length(snp),mean=0,sd = sqrt(0.96*(var(snp)/0.04)))
  
  linmodq4=summary(lm(trait ~ snp))
  fs=linmodq4$fstat
  pvals_[i] <- 1 - pf(fs[1], fs[2], fs[3])
}
#pvals_

pow1_<- mean(pvals_ < .05)  # for later comparison
# "exact" power
pow0_
# simulated power
pow1_
