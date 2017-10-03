gram_schmidt_lindep = function(matrixq1){
  x=sqrt(matrixq1[,1][1]^2+matrixq1[,1][2]^2+matrixq1[,1][3]^2+matrixq1[,1][4]^2)
  u0=matrixq1[,1]/x
  #print(u0)
  y=matrixq1[,2]-sum((matrixq1[,2])*u0)*u0
  u1=y/sqrt(y[1]^2+y[2]^2+y[3]^2+y[4]^2)
  #print(u1)
  return(matrix(c(u0,u1),byrow=F,nrow=4,ncol=2))
}
q1matrix=matrix(c(1,1,1,1,1,0,1,0,1,3,1,3),byrow=F,nrow=4)
orthogbasis = gram_schmidt_lindep(q1matrix)

hatmat=orthogbasis %*% t(orthogbasis)
all.equal(hatmat,t(hatmat))
all.equal(hatmat,hatmat %*% hatmat)

rho = matrix(c(1,-2,1,1),byrow=F,nrow=4)
a=hatmat %*% rho

rho2 = matrix(c(1,1,1,-2),byrow=F,nrow=4)
atilde=hatmat %*% rho2

##################
##################

y=matrix(c(82,79,74,83,80,81,84,81),byrow=F,nrow=8)
x=matrix(c(10,9,9,11,11,10,10,12,15,14,13,15,14,14,16,13),byrow=F,nrow=8,ncol=2)

betahat=solve(t(x) %*% x) %*% t(x) %*% y
# values of betahat_naught and betahat_1

x %*% betahat

yhat = x %*% betahat
s2 = sum((y - yhat)^2)/6 # df 6: # paramaters - rank of the X matrix
# s2 is our estimate of sigmasquared

s2 * solve(t(x) %*% x) # estimated cov matrix of bhat
diag(s2 * solve(t(x) %*% x))
sqrt(diag(s2 * solve(t(x) %*% x)))

ses <- sqrt(diag(s2 * solve(t(x) %*% x)))

lambda1 = matrix(c(1,0),byrow=F, nrow= 2, ncol = 1)
lambda2 = matrix(c(1,1),byrow=F, nrow=2, ncol=1)

# alpha-level two way confidence level of lambda prime beta

confinfbeta1= c(t(lambda1) %*% betahat - qt(1-0.05/2, df = 6) * sqrt(
  s2 * t(lambda1) %*% solve(t(x) %*% x) %*% lambda1),
  
  t(lambda1) %*% betahat + qt(1-0.05/2, df = 6) * sqrt(
    s2 * t(lambda1) %*% solve(t(x) %*% x) %*% lambda1)

)  # 6 = n-6 = 8-2

confinfbeta1plusbeta2 = c(t(lambda2) %*% betahat - qt(1-0.05/2, df = 6) * sqrt(
  s2 * t(lambda2) %*% solve(t(x) %*% x) %*% lambda2),
  
  t(lambda2) %*% betahat + qt(1-0.05/2, df = 6) * sqrt(
    s2 * t(lambda2) %*% solve(t(x) %*% x) %*% lambda2))

# c) Perform an alpha = 0.01 test for H0: B2 = 3
partclambda = matrix(c(0,1),nrow = 2,ncol = 1)
tstatc = (t(partclambda) %*% betahat - 3) / (s2 * t(partclambda) 
                                          %*% solve(t(x) %*% x) %*% partclambda)
# look up significance to the 0.01 level in my t-chart for n-r = 8-2 = 6 df

# d)
# Test 
partdlambda = matrix(c(1,-1), nrow = 2)
tstatd = (t(partdlambda) %*% betahat) / (s2 * t(partdlambda) 
                                            %*% solve(t(x) %*% x) %*% partdlambda)
# probability of being more "extreme" than tstatd