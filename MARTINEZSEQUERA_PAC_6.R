
#parametres opcio
S=105
K=110
r=0.005
sigma=0.4
T=0.5

#alternatives per calcular Phi
#formula de Hastings
te=function(x)
{
  alpha= 0.2316419
  return(1/(1+alpha*x))
}

phi=function(x)
{
  return((1/sqrt(2*pi))*exp(-x^2/2))
}

PhiHpos=function(x)
{
  a1= 0.319381530
  a2= -0.356563782
  a3= 1.781477937
  a4= -1.821255978
  a5= 1.330274429
  t=te(x)
  Phitilde= 1-phi(x)*(a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5)
  return(Phitilde)
}

PhiH=function(x)
{
  if(x>=0) return(PhiHpos(x))
  else return(1-PhiHpos(-x))
}


blackScholes <- function(S, K, r, T, sigma) {  
  d1 <- (log(S/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  call = S*PhiH(d1) - K*exp(-r*T)*PhiH(d2)
  return(call)
}

blackScholes(S, K, r, T, sigma)

#regla dels trapezis
trap=function(f,a,b, m)
{
  x=seq(a,b,length.out=m+1)
  y=f(x)
  p.area=sum((y[2:(m+1)]+y[1:m]))
  p.area=p.area*abs(b-a)/(2*m)
  return(p.area)
}

PhiTpos=function(x)
{
  return(0.5+trap(phi,0,x,m))
}

PhiT=function(x)
{
  if(x>=0) return(PhiTpos(x))
  else return(1-PhiTpos(-x))
}

#regla de Simpson
simp=function(f,a,b,m){
  
  x.ends=seq(a,b,length.out=m+1)
  y.ends=f(x.ends)
  x.mids=(x.ends[2:(m+1)]-x.ends[1:m])/2+x.ends[1:m]
  y.mids=f(x.mids)
  p.area=sum(y.ends[2:(m+1)]+4*y.mids[1:m]+ y.ends[1:m])
  p.area=p.area*abs(b-a)/(6*m)
  return(p.area)
}

PhiSpos=function(x)
{
  return(0.5+simp(phi,0,x,m))
}

PhiS=function(x)
{
  if(x>=0) return(PhiSpos(x))
  else return(1-PhiSpos(-x))
}


#Romberg (si tab=TRUE torna tota la matriu)
romberg=function(f,a,b,m,tab=FALSE){
  R=matrix(NA,nrow=m,ncol=m)
  R[1,1]=trap(f,a,b,m=1)
  for(j in 2:m){
    R[j,1]=trap(f,a,b,m=2^(j-1))
    for (k in 2:j){
      k4=4^(k-1)
      R[j,k]=k4*R[j,k-1]-R[j-1,k-1]
      R[j,k]=R[j,k]/(k4-1)
    }
  }
  if (tab==TRUE) 
    return (R)
  return (R[m,m])
}

PhiRompos=function(x)
{
  return(0.5+romberg(phi,0,x,m,FALSE))
}

PhiRom=function(x)
{
  if(x>=0) return(PhiRompos(x))
  else return(1-PhiRompos(-x))
}

#Monte Carlo
mcint=function(f,a,b,m){
  set.seed(1255)
  x=runif(m,min=a,max=b)
  y.hat=f(x)
  area= (b-a)*sum(y.hat)/m
  return(area)
}

PhiMCpos=function(x)
{
  return(0.5+mcint(phi,0,x,m))
}

PhiMC=function(x)
{
  if(x>=0) return(PhiMCpos(x))
  else return(1-PhiMCpos(-x))
}

m=100
#preu de la call
BScall=function(S, K, r, T, sigma) {
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*PhiMC(d1) - K*exp(-r*T)*PhiMC(d2)
  return(call)
}

MC100 <-BScall(S, K, r, T, sigma)

m=1000
#preu de la call
BScall=function(S, K, r, T, sigma) {
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*PhiMC(d1) - K*exp(-r*T)*PhiMC(d2)
  return(call)
}

MC1000 <- BScall(S, K, r, T, sigma)


m=10000
#preu de la call
BScall=function(S, K, r, T, sigma) {
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*PhiMC(d1) - K*exp(-r*T)*PhiMC(d2)
  return(call)
}

MC <- BScall(S, K, r, T, sigma)

#errors
m=5
#preu de la call
BScall=function(x){
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*x(d1) - K*exp(-r*T)*x(d2)
  return(call)
}

vPhiH=BScall(PhiH) #valor amb la Phi de Hastings
vPhiT=BScall(PhiT) #valor amb la Phi de trapezis
vPhiS=BScall(PhiS) #valor amb la Phi de Simpson
vPhiRom=BScall(PhiRom) #valor amb la Phi de Romberg
vPhiMC=BScall(PhiMC) #valor amb la Phi de Monte Carlo

#errors
abserrT=abs(vPhiT-vPhiH)
relerrT=abserrT/vPhiH

abserrS=abs(vPhiS-vPhiH)
relerrS=abserrS/vPhiH

abserrRom=abs(vPhiRom-vPhiH)
relerrRom=abserrRom/vPhiH

abserrMC=abs(vPhiMC-vPhiH)
relerrMC=abserrMC/vPhiH

Mètode<-c("Hastings", "Trapezis", "Simpson", "Roomberg", "Monte Carlo" )
valor<- c(vPhiH, vPhiT, vPhiS, vPhiRom, vPhiMC)
e.rel<- round(c(0,relerrT, relerrS, relerrRom, relerrMC),4)
e.ab<- round(c(0,abserrT, abserrS, abserrRom, abserrMC), 3)

data.frame(Mètode,valor, e.rel, e.ab)


m=10
#preu de la call
BScall=function(x){
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*x(d1) - K*exp(-r*T)*x(d2)
  return(call)
}

vPhiH=BScall(PhiH) #valor amb la Phi de Hastings
vPhiT=BScall(PhiT) #valor amb la Phi de trapezis
vPhiS=BScall(PhiS) #valor amb la Phi de Simpson
vPhiRom=BScall(PhiRom) #valor amb la Phi de Romberg
vPhiMC=BScall(PhiMC) #valor amb la Phi de Monte Carlo

#errors
abserrT=abs(vPhiT-vPhiH)
relerrT=abserrT/vPhiH

abserrS=abs(vPhiS-vPhiH)
relerrS=abserrS/vPhiH

abserrRom=abs(vPhiRom-vPhiH)
relerrRom=abserrRom/vPhiH

abserrMC=abs(vPhiMC-vPhiH)
relerrMC=abserrMC/vPhiH

Mètode<-c("Hastings", "Trapezis", "Simpson", "Roomberg", "Monte Carlo" )
valor<- c(vPhiH, vPhiT, vPhiS, vPhiRom, vPhiMC)
e.rel<- round(c(0,relerrT, relerrS, relerrRom, relerrMC),4)
e.ab<- round(c(0,abserrT, abserrS, abserrRom, abserrMC), 3)

data.frame(Mètode,valor, e.rel, e.ab)


m=20
#preu de la call
BScall=function(x){
  d1=(log(S/K)+(r+0.5*sigma^2)*(T))/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T)
  call=  S*x(d1) - K*exp(-r*T)*x(d2)
  return(call)
}

vPhiH=BScall(PhiH) #valor amb la Phi de Hastings
vPhiT=BScall(PhiT) #valor amb la Phi de trapezis
vPhiS=BScall(PhiS) #valor amb la Phi de Simpson
vPhiRom=BScall(PhiRom) #valor amb la Phi de Romberg
vPhiMC=BScall(PhiMC) #valor amb la Phi de Monte Carlo

#errors
abserrT=abs(vPhiT-vPhiH)
relerrT=abserrT/vPhiH

abserrS=abs(vPhiS-vPhiH)
relerrS=abserrS/vPhiH

abserrRom=abs(vPhiRom-vPhiH)
relerrRom=abserrRom/vPhiH

abserrMC=abs(vPhiMC-vPhiH)
relerrMC=abserrMC/vPhiH

abserrM =abs(MC-vPhiH)
relerrM =abserrM/vPhiH

Mètode<-c("Hastings", "Trapezis", "Simpson", "Roomberg", "Monte Carlo" , "MC 10000")
valor<- c(vPhiH, vPhiT, vPhiS, vPhiRom, vPhiMC, MC)
e.rel<- round(c(0,relerrT, relerrS, relerrRom, relerrMC, relerrM),4)
e.ab<- round(c(0,abserrT, abserrS, abserrRom, abserrMC, abserrM), 3)

data.frame(Mètode,valor, e.rel, e.ab)

