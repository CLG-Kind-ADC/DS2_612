library(MASS)

# From last time about body fat regression ------------
# the old body fat data

url = "http://www.stat.yale.edu/~jtc5/312_612/data/bodymeasurements.txt"

bodyfat0 = read.table(url, header=T, sep="\t")[,-1]
set.seed(1)
n = 30
bodyfat = bodyfat0[sample(1:nrow(bodyfat0),n),]
dim(bodyfat)
names(bodyfat)
attach(bodyfat)
bf = Pct.BF
ht = Height
ab = Abdomen
ab = ab/2.54  

plot(bf~ht)
cor(bf, ht)

lmh = lm(bf~ht)
summary(lmh)
# ^^ not surprising; it's not true that taller people have systematically 
# less or more percent bodyfat...

plot(bf ~ ab)
lma = lm(bf ~ ab)
summary(lma)

lmha = lm(bf ~ ht + ab)
summary(lmha)

# The two variables ab and ht are useful *together* in predicting bf
# even though ht wasn't useful by itself.  
# If we fix a value for ab, then ht becomes useful in predicting bf.

# Do the tests from scratch...
#lmha t tests:
out <- summary(lmha)

bf
ht
ab
X <- cbind(1, ht,ab)
X
y <- bf
bhat <- solve(t(X) %*% X) %*% t(X) %*% y
bhat
out

    # equiv:
    solve(t(X)%*%X, t(X)%*%y)

# covariance matrix of bhat is sig2 (X'X)^(-1)
# estimated cov matrix of bhat is s2 (X'X)^(-1)
#     s2 is sum of (y-yhat)^2 / (n-r)
yhat <- X %*% bhat
yhat
fitted(lmha)

# "Residual standard error" in out
s2 <- sum((y-yhat)^2)/(30-3)
s2
sqrt(s2)
out

# "Std. Error" column in out
s2 * solve(t(X) %*% X) # estimated cov matrix of bhat
diag(s2 * solve(t(X) %*% X))
sqrt(diag(s2 * solve(t(X) %*% X)))
out 
ses <- sqrt(diag(s2 * solve(t(X) %*% X)))
ses

# "t value" column in out:
tstats <- bhat / ses
tstats
out

# "Pr(>|t|)" column in out
2 * (1 - pt(abs(tstats),df=30-3))
out

2 * (1 - pt(abs(tstats),df=28)) 
# Yes, not the same as in out. Can tell degrees of freedom used is indeed 30-3.

# New stuff about body fat ------------------------

# F tests:
ybar <- mean(y)
Fstat <- sum((yhat - ybar)^2)/(3-1) / (sum((y-yhat)^2)/(30-3))
Fstat
s
1 - pf(q=Fstat,df1=3-1,df2=30-3)
1 - pf(q=Fstat,df1=3-1,df2=28)


# multiple R^2
sum((yhat-ybar)^2)/sum((y-ybar)^2)
s
cor(y, yhat)^2
plot(yhat,y)

1 - sum((y-yhat)^2)/sum((y-ybar)^2)

# Adjusted R^2
1 - (sum((y-yhat)^2)/(30-3)/var(y))

# Can be negative:
set.seed(1)
x <- 1:10
yy <- rnorm(10)
y <- lm(yy~x)$res
plot(x,y)
cor(x,y)
summary(lm(y~x))
# Note adjusted R-squared is -0.125

