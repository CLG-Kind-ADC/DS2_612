#' HW 1 STAT 612

# problem 1
## regress pressure on boiling point
library(alr4)
plot(x=Forbes$bp,y=Forbes$pres)
forbes=lm(Forbes$pres~Forbes$bp)
abline(a=coef(forbes[1],b=coef(forbes[2])),col="blue")

## regress log of pressure on boiling point
plot(x=Forbes$bp,y=Forbes$lpres)
forbes2=lm(Forbes$lpres ~ Forbes$bp)
forbes2
abline(a=coef(forbes2[1]),b=coef(forbes2[2]),col="red")

## plot both
plot(x=Forbes$bp,y=Forbes$pres,ylim=c(0,160),xlim=c(195,210))
points(x=Forbes$bp,y=Forbes$lpres)
abline(a=coef(forbes[1],b=coef(forbes[2])),col="blue")
abline(a=coef(forbes2[1]),b=coef(forbes2[2]),col="red")

plot(fitted(forbes), residuals(forbes))
plot(fitted(forbes2), residuals(forbes2))
#' Notice: besides the one outlier of which we were warned, the
#' residuals vs fitted of the log of pressure
#' shows relatively little pattern: a good thing.
#' 
#' However, the graph of pressure's residuals shows some pattern.
#' This is a little concerning.

# #5
#' Let o1, ... , or be an orthonormal basis for C(X)
#' O=[o1,...,or]
#' Then OO' is the perpendicular projection operator onto C(X)
#' 

A = matrix(c(2,0,4,1,5,7,1,-5,-3),nrow=3,ncol=3,byrow=T)
B = matrix(c(1,0,0,0,0,1,0,1,0),nrow=3,ncol=3,byrow=T)
C= matrix(c(1,4,1,2,5,1,-3,0,1),nrow=3,byrow=T)

# gram schmidt
#a function?

gram_schmidt_lindep = function(matrix3x3){
  x=sqrt(matrix3x3[,1][1]^2+matrix3x3[,1][2]^2+matrix3x3[,1][3]^2)
  u0=matrix3x3[,1]/x
  print(u0)
  y=matrix3x3[,2]-sum((matrix3x3[,2])*u0)*u0
  u1=y/sqrt(y[1]^2+y[2]^2+y[3]^2)
  print(u1)
  z=matrix3x3[,3]-sum(matrix3x3[,3]*u0)*u0-sum(matrix3x3[,3]*u1)*u1
  u2=z/sqrt(z[1]^2+z[2]^2+z[3]^2)
  print(u2)
  return(matrix(c(u0,u1),byrow=F,nrow=3,ncol=2))
}

gram_schmidt = function(matrix3x3){
  x=sqrt(matrix3x3[,1][1]^2+matrix3x3[,1][2]^2+matrix3x3[,1][3]^2)
  u0=matrix3x3[,1]/x
  print(u0)
  y=matrix3x3[,2]-sum((matrix3x3[,2])*u0)*u0
  u1=y/sqrt(y[1]^2+y[2]^2+y[3]^2)
  print(u1)
  z=matrix3x3[,3]-sum(matrix3x3[,3]*u0)*u0-sum(matrix3x3[,3]*u1)*u1
  u2=z/sqrt(z[1]^2+z[2]^2+z[3]^2)
  print(u2)
  return(matrix(c(u0,u1,u2),byrow=F,nrow=3,ncol=3))
}

# these are orthonormal bases
Aorthonormal = gram_schmidt_lindep(A)
Borthonormal = gram_schmidt(B)
Corthonormal = gram_schmidt_lindep(C)

# projection matrices
projA = Aorthonormal%*%t(Aorthonormal)
projB = Borthonormal%*%t(Borthonormal)
projC = Corthonormal%*%t(Corthonormal)

projA
projB
projC