
#Tipus d'interes
ti1= 0.0144
ti2= 0.0145
ti3= 0.0149
ti4= 0.0151
ti5= 0.0155 
ti6= 0.0159 
ti7= 0.0165
ti8= 0.0169
ti9= 0.0178
ti10= 0.0181

r=c(ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,ti9,ti10) 

#venciments
T=c(1/12, 1/6, 1/4, 0.5, 1,2,3,5,7,10)

#Venciments a interpolar
T1= 1.5/12
T2= 10/12
T3= 6.5
T4= 9

#calcul del polinomi interpolador (Vandermonde)

polyinterp=function(x,y)
{
  if(length(x)!=length(y))
    stop ("La longitud de los vectores x e y debe ser la misma" )
  n=length(x)-1
  vandermonde=rep(1,length(x))
  for(i in 1:n){
    
    xi=x^i
    vandermonde=cbind(vandermonde,xi)
  }
  beta=solve(vandermonde,y)
  names(beta)=NULL
  return(beta)
}

#regla de Horner per avaluar polinomis
horner=function(x,coefs)
{
  y=rep(0,length(x))
  for(i in length(coefs):1)
    y=coefs[i]+x*y
  return(y)
}

p=polyinterp(T,r)
p

#resultats interpolacio amb el polinomi interpolador
r1p=horner(T1,p)
r2p=horner(T2,p)
r3p=horner(T3,p)
r4p=horner(T4,p)

#grafic del polinomi interpolador

library(ggplot2)

#grafic del polinomi interpolador
x=seq(1/12,10,0.5)
y=horner(x,p)

dat <- data.frame(cbind(x, y))

ggplot(dat, aes(x=x, y=y)) + 
  geom_point(size=3, col='blue')+
  stat_function(fun =p,  colour = 'lightblue')

#lagrange

lagrange= function(x,y,a){
  n =length(x)
  if(a < min(x) || max(x)<a)stop ("no s´està interpolant")
  X = matrix(rep(x, times=n), n, n, byrow = T)
  mN= a -X; diag(mN)= 1
  mD= X- t(X); diag(mD)= 1
  Lnk=apply(mN, 1, prod)/apply(mD, 2, prod)
  sum(y*Lnk)
}

T=c(1/12, 1/6, 1/4, 0.5, 1,2,3,5,7,10)
r= c(ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,ti9,ti10)

lagrange(T,r, T1)
lagrange(T,r, T2)
lagrange(T,r, T3)
lagrange(T,r, T4)


library(polynom)
library(PolynomF)

polyLagr=poly.calc(T,r)

polyLagr

x=seq(1/12,10,0.5)
y=horner(x,polyLagr)

dat <- data.frame(cbind(x, y))

ggplot(dat, aes(x=x, y=y)) + 
  geom_point(size=3, col='blue')

library(pracma)

xs = c(T1,T2,T3, T4)

(newtonInterp(T, r, xs))

(lagrangeInterp(T, r, xs))

#interpolacio per trams lineal
#r=m*T+b
linterp=function(x1,y1,x2,y2)
{
  m=(y2-y1)/(x2-x1)
  b=y2-(m*x2)
  return(c(b,m))
}
#resultats interpolacio per trams lineal
p1=linterp(1/12,ti1,1/6,ti2)
r1=p1[[2]]*T1+p1[[1]]

p2=linterp(0.5,ti4,1,ti5)
r2=p2[[2]]*T2+p2[[1]]

p3=linterp(5,ti8,7,ti9)
r3=p3[[2]]*T3+p3[[1]]

p4=linterp(7,ti9,10,ti10)
r4=p4[[2]]*T4+p4[[1]]

