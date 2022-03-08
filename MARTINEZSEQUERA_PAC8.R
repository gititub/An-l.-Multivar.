myReadData_byDate = function(file, date_ini, date_fin, column){
  
  df = read.csv(file, sep = ',', row.names = 1)
  idx_Dates = as.character( seq(as.Date(date_ini, format = '%d/%m/%Y'), as.Date(date_fin, format = '%d/%m/%Y'), 'days') )
  return( na.omit(df[idx_Dates, column]) )
}

#EXERCICI 1
Y = myReadData_byDate("C:/Users/hp/Downloads/WHO-COVID-19-global-data-SPAIN.csv", '15/07/2020', '15/11/2020', 'New_deaths')

m = length(Y)
X = 1:m # definir la variable que representa als dies
data= data.frame(X,Y)
# Representacio de les dades
library(ggplot2)
p<- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("New deaths") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(0, 120, by =10))

#Model exponencial

lny= log(Y)
exp_model<- nls(lny~ log(a)+ b*X)

sum= summary(exp_model)

cor(lny, predict(exp_model))

a = coefficients(exp_model)[1]
b = coefficients(exp_model)[2]

#Recuperem la fòrmula original:
y1= a*exp(b*X)

p + geom_smooth(aes(x=X, y=y1))

# Equació potència

logy = log10(Y)
logx = log10(X)

pot_model<- nls(logy ~ log10(a2) + b2*logx)

sum2= summary(pot_model)

#Estimem la bondat d´ajust

cor(logy, predict(pot_model))

#Recuperem la fòrmula original
a2 = coefficients(pot_model)[1]
b2 = coefficients(pot_model)[2]

y2= a2*(X^b2)

p + geom_smooth(aes(x=X, y= y2))

# Equació de la taxa de creixement de saturació

#1/y = a0 + a1*1/X
#a0 = 1/a3
#a1 = b3/a3

invy = 1/Y
invx = 1/X

sat_model<- nls(invy ~ 1/a3 + (b3/a3)*invx)

sum3= summary(sat_model)

cor(invy, predict(sat_model))

a3 = coefficients(sat_model)[1]
b3 = coefficients(sat_model)[2]

y3= (a3*X)/(b3+X)

p + geom_smooth(aes(x=X, y=y3))

# 2. Realitzar el mateix ajust mitjançant regressió lineal (n = 1).

# Construccio de la matriu d'equacions normals

myPhi = function(x, n){
  
  Phi = matrix(1, length(x), n+1)
  for (i in 1:n) {
    Phi[, i+1] = x^i
  }
  return(Phi)
  
}

# Funcio que avalua un ajust polinomic definit pel coeficients "c" en els punts "x" (que pot ser un vector de punts)
myeval = function(x, c){
  
  f = 0
  for (i in 1:length(c)) {
    f = f + c[i]*x^(i-1)
  }
  return(f)
}


#grau 1
n=1
lout = lsfit(myPhi(X, n), Y)
ckP = lout[[1]] #coeficients

fP = myeval(X, ckP)

#Alternativa:
model<- lm(Y~X)
coef<- coef(model)
f<- myeval(X, coef)

coefs<- summary(model)$coefficients

summary(model)

p + geom_abline(intercept = coefs[1], 
                slope = coefs[2],
                colour="green")

#Estimem la bondat d´ajust
cor(Y, predict(model))

p + geom_smooth(aes(x=X, y= y2), colour="red")+
  geom_smooth(aes(x=X, y= y1))+
  geom_abline(intercept = coefs[1], 
              slope = coefs[2],
              colour="green")

St = sum( (Y - mean(Y))^2 ) # discrepancies de les dades respecte a la mitjana
Sr_E = sum( (Y - y1)^2 ) # exponencial
Sr_P = sum( (Y - y2)^2 ) # potència
Sr_L = sum( (Y - f)^2 )# sumatori dels residus per a l'ajust lineal

r_E = sqrt((St - Sr_E)/St)# coeficient de correlacio- exponencial
r_P = sqrt((St - Sr_P)/St)# coeficient de correlacio - potència
r_L = sqrt((St - Sr_L)/St)# coeficient de correlacio - lineal



x10 = m+10 # abscissa de la prediccio

#Predicció amb funció potència
y10_P= a2*(x10^b2)

#Predicció amb exponencial
y10_E= a*exp(b*x10)

#Predicció amb regressió lineal
y10_L = coefs[1] + coefs[2]*x10

#EXERCICI 2

Y = myReadData_byDate("C:/Users/hp/Downloads/WHO-COVID-19-global-data-SPAIN.csv", '03/01/2020', '03/09/2021', 'Cumulative_cases')

m = length(Y)
X = 1:m# definir la variable que representa als dies

# Representacio de les dades
data= data.frame(X,Y)

library(ggplot2)
p<- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("cumulative cases") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(10, 610, by =40)) 

seg = 7
color = c('cyan', 'yellow', 'magenta', 'blue', 'red', 'orange', 'purple')
for (i in 1:seg){ 
  ini_segmentos = c(1,75,200,330,370,400,550)
  fin_segmentos = c(75,200,330,370,400,550,610)
  s = ini_segmentos[i]:fin_segmentos[i]
  Xs = X[s]
  Ys = Y[Xs] 
  lm = lm(Ys ~ Xs, data=data)
  ck = coef(lm) 
  f = myeval(Xs, ck)
  p <- p + geom_abline(intercept=ck[1],slope =
                         ck[2],col=color[i])
}


#n= 1
lm.n1<- lm(Y ~ X)
my_coef1<-coef(lm.n1)
# Representacio de les dades
data= data.frame(X,Y)

p1 <- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("cumulative cases n=1") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(10, 610, by =40))


p1<- p1 + geom_abline(intercept = my_coef1[1], 
                      slope = my_coef1[2],
                      colour="green")

#n=2
lm.n2<- lm(Y ~ poly(X, degree = 2, raw = TRUE))
my_coef2<-coef(lm.n2)
f2<- myeval(X, my_coef2)


p2 <- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("cumulative cases n=2") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(10, 610, by =40))


p2<- p2 + stat_smooth(method = "lm",
                      formula = Y ~ poly(X, 2, raw = TRUE),
                      size = 1)
summary(lm.n1)
summary(lm.n2)

#plot

prediccions <- predict(lm.n2, newdata = data, se.fit = TRUE,level = 0.95)

# càlcul del interval de confiança superior i inferior al 95%
interval_conf <-data.frame(inferior= prediccions$fit-                                     1.96*prediccions$se.fit,
                           superior = prediccions$fit+
                             1.96*prediccions$se.fit)

plot(x = X, y = Y, pch = 20, col = "darkgrey")
title("Polinomi de grau 2")
lines(x = X, prediccions$fit, col = "red", lwd = 1)
lines(x = X, interval_conf$inferior, col = "blue", lwd = 1)
lines(x = X, interval_conf$superior, col = "blue", lwd = 1)

#anova
anova(lm.n1,lm.n2)

#punts de tall
x1= 1:74
x2=75:200
x3=200:330
x4=330:370
x5=370:400
x6=400:550
x7=551:610

lm1<- lm(Y[x1]~x1)
lm2<- lm(Y[x2]~x2)
lm3<- lm(Y[x3]~x3)
lm4<- lm(Y[x4]~x4)
lm5<- lm(Y[x5]~x5)
lm6<- lm(Y[x6]~x6)
lm7<- lm(Y[x7]~x7)

a1=coef(lm1)[1]
b1=coef(lm1)[2]
a2=coef(lm2)[1]
b2=coef(lm2)[2]

x=(a1-a2)/(b2-b1)

a1=coef(lm2)[1]
b1=coef(lm2)[2]
a2=coef(lm3)[1]
b2=coef(lm3)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm3)[1]
b1=coef(lm3)[2]
a2=coef(lm4)[1]
b2=coef(lm4)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm4)[1]
b1=coef(lm4)[2]
a2=coef(lm5)[1]
b2=coef(lm5)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm5)[1]
b1=coef(lm5)[2]
a2=coef(lm6)[1]
b2=coef(lm6)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm6)[1]
b1=coef(lm6)[2]
a2=coef(lm7)[1]
b2=coef(lm7)[2]
x=(a1-a2)/(b2-b1)

#paquet segmented de R:

library(segmented)

my.seg <- segmented(lm.n1, 
                    seg.Z = ~ X,
                    psi = NA)
summary(my.seg)

#Punts de tall:
my.seg$psi

#pendents:
slope(my.seg)

#Representació gràfica dels segments:

my.fitted <- fitted(my.seg)
my.model <- data.frame(x=X, y = my.fitted)

ggplot(my.model, aes(x,y)) + geom_line()


#REPETIR EXERCICI CUMULATIVE_DEATHS
Y = myReadData_byDate("C:/Users/hp/Downloads/WHO-COVID-19-global-data-SPAIN.csv", '03/01/2020', '03/09/2021', 'Cumulative_deaths')

m = length(Y)
X = 1:m# definir la variable que representa als dies

# Representacio de les dades
data= data.frame(X,Y)

library(ggplot2)
P<- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("cumulative deaths") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(10, 610, by =40))+
  theme(plot.title = element_text(hjust = 0.5)) 

# 2B. Regressió segmentada de "Cumulative_deaths":

seg = 7
color = c('cyan', 'yellow', 'magenta', 'blue', 'red', 'orange', 'purple')
for (i in 1:seg){ 
  ini_segmentos = c(1,75,220,310,380,430,560)
  fin_segmentos = c(75,220,310,380,430,560,610)
  s = ini_segmentos[i]:fin_segmentos[i]
  Xs = X[s]
  Ys = Y[Xs] 
  lm = lm(Ys ~ Xs, data=data)
  ck = coef(lm) 
  f = myeval(Xs, ck)
  P <- P + geom_abline(intercept=ck[1],slope = ck[2],col=color[i])
  
}

#n= 1

lm.n1<- lm(Y ~ X)
my_coef1<-coef(lm.n1)

# Representacio de les dades
data= data.frame(X,Y)

p1 <- ggplot(data, aes(x = X, y = Y)) +
  geom_point(col = "darkgrey") +
  ggtitle("cumulative cases n=1") +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_line()) +
  scale_x_continuous(breaks = seq(10, 610, by =40))


p1<- p1 + geom_abline(intercept = my_coef1[1], 
                      slope = my_coef1[2],
                      colour="green")

#n=2
lm.n2<- lm(Y ~ poly(X, degree = 2, raw = TRUE))
my_coef2<-coef(lm.n2)
f2<- myeval(X, my_coef2)


prediccions <- predict(lm.n2, newdata = data, se.fit = TRUE,level = 0.95)

# càlcul del interval de confiança superior i inferior al 95%
interval_conf <-data.frame(inferior= prediccions$fit-                                     1.96*prediccions$se.fit,
                           superior = prediccions$fit+
                             1.96*prediccions$se.fit)

plot(x = X, y = Y, pch = 20, col = "darkgrey")
title("Polinomi de grau 2")
lines(x = X, prediccions$fit, col = "red", lwd = 1)
lines(x = X, interval_conf$inferior, col = "blue", lwd = 1)
lines(x = X, interval_conf$superior, col = "blue", lwd = 1)

#punts de tall

ini_segmentos = c(1,75,220,310,380,430,560)
fin_segmentos = c(75,220,310,380,430,560,610)

x1=ini_segmentos[1]:fin_segmentos[1]
x2=ini_segmentos[2]:fin_segmentos[2]
x3=ini_segmentos[3]:fin_segmentos[3]
x4=ini_segmentos[4]:fin_segmentos[4]
x5=ini_segmentos[5]:fin_segmentos[5]
x6=ini_segmentos[6]:fin_segmentos[6]
x7=ini_segmentos[7]:fin_segmentos[7]

lm1<- lm(Y[x1]~x1)
lm2<- lm(Y[x2]~x2)
lm3<- lm(Y[x3]~x3)
lm4<- lm(Y[x4]~x4)
lm5<- lm(Y[x5]~x5)
lm6<- lm(Y[x6]~x6)
lm7<- lm(Y[x7]~x7)

a1=coef(lm1)[1]
b1=coef(lm1)[2]
a2=coef(lm2)[1]
b2=coef(lm2)[2]

x=(a1-a2)/(b2-b1)

a1=coef(lm2)[1]
b1=coef(lm2)[2]
a2=coef(lm3)[1]
b2=coef(lm3)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm3)[1]
b1=coef(lm3)[2]
a2=coef(lm4)[1]
b2=coef(lm4)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm4)[1]
b1=coef(lm4)[2]
a2=coef(lm5)[1]
b2=coef(lm5)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm5)[1]
b1=coef(lm5)[2]
a2=coef(lm6)[1]
b2=coef(lm6)[2]
x=(a1-a2)/(b2-b1)

a1=coef(lm6)[1]
b1=coef(lm6)[2]
a2=coef(lm7)[1]
b2=coef(lm7)[2]
x=(a1-a2)/(b2-b1)

#Step function
modelo.step = lm(Y ~ cut(X, 7))
summary(modelo.step)
pred <- predict(modelo.step, data)
plot(Y ~ X, col="darkgrey")
lines(X, pred, col="red", lwd=2)
title("Step function, cuts = 7")

# REGRES LOCAL
modelo.local <- loess(Y ~ X, span = 0.7)

ggplot(data = data, aes(x = X, y =Y)) +
  geom_point(col = "darkgrey") +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.2, 
              color = "orange", se = F) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.7, 
              color = "brown", se = F) +
  labs(title = "Regressió local") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = "span = 0.2", x = 400, y = 40), size = 5, 
            colour = "orange") +
  geom_text(aes(label = "span = 0.7", x = 400, y = 10000), size = 5, 
            colour = "brown")
