# Power example
# (Cuthbert Daniel -> Scheffe -> Stapleton)
#  "Suppose 4 different kinds of alloy steel are prepared by varying the
#  method of manufacture. It is expected that the tensile strength
#  will be of order 150,000 psi and the SD of duplicate specimens from
#  the same batch will be about 3000 psi.
#
#  Suppose n=10 specimens of each kind are tested and that
#  mu1=150000, mu2=149000, mu3=148000, mu4=153000
#  What is the power of the resulting alpha=.05 F test?"

#+++++++++++++++++++++++++++++++++++++++++
# Do power by simulation -----------------
#+++++++++++++++++++++++++++++++++++++++++
library(MASS)
# make some data
n <- 10
mu <- rep(c(150,149,148,153),each=n) # units are 1000 psi
sig <- 3
y <- rnorm(4*n, mean=mu, sd=sig)

batch <- rep(letters[1:4], each=n)
#plot(y ~ batch)  # error

batch <- as.factor(rep(letters[1:4], each=n))
plot(y ~ batch)
plot(batch,y)

lm1 <- lm(y ~ batch)
# Can see X matrix with the command:  
model.matrix(lm1)
# (We'll talk about this more shortly.)

summary(lm1)
# I want to get the P value of the F test.  Where is it?
s <- summary(lm1)
names(s)
s$fstat
sf <- s$fstat
args(pf)
1 - pf(sf[1], sf[2], sf[3])
pval <- 1 - pf(sf[1], sf[2], sf[3])
pval
s

# Do a simulation to approximate the power for this alternative,
# say for a test of "size" .05.
#
# Start by putting the above in a loop. 
# Repeat nrep times and record all the P values.
nrep <- 1000
pvals <- numeric(nrep)
batch <- as.factor(rep(letters[1:4], each=n))
for(it in 1:nrep){
  mu <- rep(c(150,149,148,153),each=n) # units are 1000 psi
  y <- rnorm(4*n, mean=mu, sd=sig)
  lm1 <- lm(y~batch)
  s <- summary(lm1)
  sf <- s$fstat
  pvals[it] <- 1 - pf(sf[1], sf[2], sf[3])
}

pvals
truehist(pvals)
mean(pvals < .05)  # <-- power

pow1 <- mean(pvals < .05)  # for later comparison

# A very good sanity check to do: make sure P values are uniform under null:
nullpvals <- numeric(nrep)
for(it in 1:nrep){
  mu <- rep(c(150,150,150,150),each=n) # units are 1000 psi
  y <- rnorm(4*n, mean=mu, sd=sig)
  batch <- as.factor(rep(letters[1:4], rep(n,4)))
  lm1 <- lm(y~batch)
  s <- summary(lm1)
  sf <- s$fstat
  nullpvals[it] <- 1 - pf(sf[1], sf[2], sf[3])
}

#truepvals
truehist(nullpvals)
mean(nullpvals < .05)

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# power via noncentral F distribution formulas -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++

n2 <- function(x)sum(x^2) # function to 
# --- a ---
# including personal recaps
# units are 1000 psi

n <- 10 # number of each 
mu <- rep(c(150,149,148,153),each=n) # observed means thru testing
sig <- 3 # SD is 3000

mu0 <- mean(mu)*rep(1,4*n) # null? all mu's 150 (here, the "mean of the means")
ncp <- n2(mu-mu0)/sig^2 # squared difference between [alt?] and [null?], div by SD
# SD assumed same for both?

crit <- qf(.95, 4-1, 4*n-4) # crit[ical value]
# value that, for our dsn with 3 numdf and 36 dendf, is more extreme than 95% of the data
# (95th percentile)

pow <- 1 - pf(crit, 4-1, 4*n-4, ncp) #~10.5% chance of being 95th percentile or greater.
# 1 - 10.5% means ~89.5% of being lower than 95th percentile.
# this is our power
pow
#^^ compare with pow1 above.
pow1


# E.g. suppose you are planning an experiment and want to look at how power increases with sample size
mypow <- function(n){ #(varying on n: number of obs to be taken FOR EACH manufacturer)
  mu <- rep(c(150,149,148,153),each=n) # units are 1000 psi
  mu0 <- mean(mu)*rep(1,4*n)
  sig <- 3
  ncp <- n2(mu-mu0)/sig^2
  crit <- qf(.95, 4-1, 4*n-4)
  pow <- 1 - pf(crit, 4-1, 4*n-4, ncp)
  return(pow)
}

mypow(10)
ns <- 5:15
pows <- ns
for(i in 1:length(ns))pows[i] <- mypow(ns[i])
plot(ns, pows, type='b')
abline(h=.8, lty=2) #shows the point at which power exceeds 0.8 (80%)

#++++++++++++++++++++++++++++++++++++++++++
# Stapleton's second question -------------
#++++++++++++++++++++++++++++++++++++++++++

# Suppose we want to design a study that will have power at least 0.90 for alpha=0.05
# in the case that any two of the four means differ by 8000 or more.
# How many observations should be taken on each alloy?

# Trial and error; start with 10 observations per alloy.
# evidently, we were able to go down to 6 and stay above 95% power.

# --- b ---
n <- 6
crit <- qf(p = .95, df1 = 4-1, df2 = 4*n-4)
ncp <- n2(c(rep(4,n), rep(-4,n), rep(0,2*n))) / sig^2
pow <- 1 - pf(crit, 4-1, 4*n-4, ncp)
pow




#+++++++++++++++++++++++++++++++++++++++++
# Added variable plots -------------------
#+++++++++++++++++++++++++++++++++++++++++
url <- "http://www.stat.yale.edu/~jtc5/312_612/data/bodymeasurements.txt"
bodyfat <- read.table(url, header=T, sep="\t")[,-1]
bodyfat0 <- bodyfat
set.seed(12345)
n <- 30
bodyfat <- bodyfat0[sample(1:nrow(bodyfat0),n),] #takes random 30 rows from the total possible rows (250)
# (as random as it can be, with that seed set)
attach(bodyfat)
bf <- Pct.BF
ht <- Height
ab <- Abdomen
ab <- ab/2.54 #why divide by 2.54?

lmha <- lm(bf ~ ht + ab)
summary(lmha)

# We're interested in visualizing the effect of ht **in the presence of ab**, 
# that is, does adding ht enhance the explanatory ability of ab?
# We could look at it this way: after we do our best with ab, take
# residuals to get the part of bf that is not explained by ab,
# and then see how well ht can explain those residuals.

    # Aside on treatment of missing values
    lm(bf~ab)$resid
    resid(lm(bf~ab))

    ab.tmp <- ab
    ab.tmp[2] <- NA
    lm1 <- lm(bf~ab.tmp)
    lm1$resid
    resid(lm1)
    lm2 <- lm(bf~ab.tmp, na.action=na.exclude)
    lm2$resid
    resid(lm2) 
    #^^ this way keeps n the same and puts NA's in appropriately in residuals

# --- c ---
bfr <- lm(bf~ab)$resid
plot(bfr ~ ht)

# compare the regression coefficients:
lm(bfr ~ ht)
lm(bf~ht+ab)

# compare the residuals
plot(lm(bf~ht+ab)$resid, lm(bfr ~ ht)$resid)

# Another idea:
htr <- lm(ht ~ ab)$resid
lm(bf ~ htr)
lmha
#^^ same coef of ht!
# How about the residuals for the two regressions?
plot(lmha$residuals, lm(bf ~ ht)$resid) # bad : residuals tend upward
plot(lmha$residuals, lm(bf ~ htr)$resid) # still tending somewhat upward...

# Try residualizing both bf and ht: added variable plot
bfr <- lm(bf~ab)$resid
htr <- lm(ht~ab)$resid
plot(htr,bfr) # little pattern to the residuals (residualized bf, ht on ab)

# compare coefficients:
lm(bfr~htr)
lm(bf~ht+ab)

# compare residuals
plot(lm(bf~ht+ab)$resid, lm(bfr~htr)$resid)
abline(0,1)

summary(lm(bfr~htr))

lm(bfr~htr)
lm(bf~htr)
plot(bfr~htr)
plot(bf~htr)
plot(lm(bfr~htr)$resid, lm(bf~htr)$resid)

