

###########AJUSTE DE DISTRIBUCIONES PARA LA VARIABLE INGRESO, PARETO, GAMA, GAMA GENERAL, BETA GENERAL. 
###########CON RESTRICCION CUENTAS NACIONALES , ES DECIR, RESTICION A QUE LA ESPERANZA DE CADA DISTRIBUCION 
########### SEA IGUAL A EL PROMEDIO REPORTADO POR CUENTAS NACIONALES. 

#### librerias
library(alabama)
library(gamlss)
library(survey)
library(laeken)

#### leer datos

prop<-read.csv("D:\\survey\\proprorciones_pib.csv") # proporciones del pib
totales<-read.csv("D:\\survey\\totales.csv") #totales cuentas nacionales por año 
val<-read.csv("D:\\survey\\val_iniciales_rescate.csv")
factores<-read.csv("D:\\survey\\factores.csv")  # datos de la muestra total de hogares por entidad y por año desde el 2008
años<-c(2004,2006,2008,2010,2012)



### prueba con j=5 osea el año 2012
for (j in 1:length(años))  ## conteo por año 

### prueba con j=5 osea el año 2012
for (j in 1:length(años))  ## conteo por año 

{
  info<-paste ("ingresos",años[j],sep="")

data<-read.csv(paste("D:\\survey\\",info,".csv",sep=""))  # datos de la muestra
#data<-read.csv(paste("D:\\survey\\mcs12",".csv",sep="")) # datos del mocoso

jj<-which(totales[,1]==años[j])
tot<- totales[jj,2]*1000000 # total ingreso por año j cuentas nacionales 

######################## para obtener el total del ingreso

  ent<-data
  muestra<-ent$ing_cor
  fac<-ent$factor_hog
  total_ingreso <-sum(muestra*fac) # total ingreso muestra 
  
  tot_hog<-sum(fac) # total de hogares de la poblacion 
##  QUITAREMOS TODOS LOS VALORES MENORES A 2000
w<-c()
w<-which(muestra<2000)
if(length(w)==0)
{muestra<-muestra} else 
{muestra<-muestra[-w]
 fac<-fac[-w]
}

fquitados<-tot_hog-sum(fac)  ## re-austar los valores de expansión 
coci<-fquitados/length(fac)
nfac<-(fac+coci)
fac<-nfac  ## nuevos factores ajustados 

########## SI LOS DATOS SON SAT 


data<-read.csv("D:\\survey\\resumen.csv",sep="\t")
muestra<-data$v4
w<-which(is.na(muestra))  ## quitando NA
muestra<-muestra[-w]
### NOTA: QUITAR LOS PESOS A LAS FUNCIONES LOGVERO PARA AGILIZAR EL PROCESO.

###########tot_hog<-31559379 3.8 personas por hogar aprox 4
fac<-rep(1,length(muestra))  ## vector de unos 
ci <- 92733.62   ###

######################## para obtener el total del ingreso por entidad
t<-c()
pro<-c()
data<-read.csv(paste("D:\\survey\\mcs12",".csv",sep="")) # datos del mocoso
for (ei in 3:5) ####ciclo por entidad
{
  ## conteo de entidad i prueba con aguascalientes 
  suma_total<-sum(data$ing_cor*data$factor_hog)
  ent<-subset(data, data$entidad==ei) # archivo de la entidad 
  #ent<-data
  muestra<-ent$ing_cor
  
  fac<-ent$factor_hog
  t[ei]<-sum(muestra*fac)
  pro[ei]<-t[ei]/suma_total
}
tot_hog<-factores$tot_hog2012[ei]



for (ei in 20:32) ####ciclo por entidad
{
 ## conteo de entidad i prueba con aguascalientes 
ent<-subset(data, data$entidad==ei) # archivo de la entidad 
muestra<-ent$ing_cor
fac<-ent$factor
tot_hog<-factores$tot_hog2012[ei]
#tot_hog<-sum(fac)
##  QUITAREMOS TODOS LOS VALORES MENORES A 2000
w<-c()
w<-which(muestra<2000)
if(length(w)==0)
 {muestra<-muestra} else 
 {muestra<-muestra[-w]
 fac<-fac[-w]
 }

 fquitados<-tot_hog-sum(fac)  ## re-austar los valores de expansión 
 coci<-fquitados/length(fac)
 nfac<-(fac+coci)
 fac<-nfac


#pro<-prop$X20012p[ei] # valor de proporción del ingreso entidad i
proc<-pro[ei] 
tote<-tot*proc/4  # total por entidad i trimestral CN 
ci <- c(tote/tot_hog)  ##promedio estatal de ingreso segun cuentas nacionaes
ci<-round(ci,2)
#ci


probs = c(0.1, 0.5, 1, 2, 5, 10,20,30,40, 50,60,70,80,90,95,98,99,99.5,99.9)/100  # vector probabilidades para quantiles
s<-seq(min(muestra),max(muestra),length.out = 10000)  ## para calcular valores de distri solo para lognormal y beta II



########################AJUSTES##############################

### Distribución GB 4 parámetros 

hin<-function(x)  # desigualdad todos los parámetros positivos 
{
  h <- rep(NA, 3)
  h[1]<-x[1]
  h[2]<-x[3]
  h[3]<-x[4]
  h
}



heq<- function(x,y,z,w)  #esperanza igual a constante ci, cuentas nacionales
{ h <- rep(NA, 1)
  h[1]<-ci-(x[1]*beta(x[3]+(1/x[2]),x[4]-(1/x[2]))/beta(x[3],x[4]))  
  h
}

# estimadores de la muestra caso: sin restricción 
modbe2<-gamlss(muestra~1, family=GB2,weights=fac,control=gamlss.control(c.crit = 0.001, n.cyc = 50))

modbe2<-refit(modbe2)   ## si que no converge volver a ajustar
#fit <- vglm(muestra ~ 1, genbetaII, trace = TRUE)
## valores iniciales de los estiomadores
theta<-as.vector(c(fitted(modbe2,c("mu"))[1],fitted(modbe2,c("sigma"))[1],fitted(modbe2,c("nu"))[1],fitted(modbe2,c("tau"))[1]))  
betamues<- constrOptim.nl(par=theta, fn=fbeta2, hin=hin)
xGBsin<-betamues$par # estimadores sin restricción
x<-xGBsin # solo para calculos
ExGBsin<-x[1]*beta(x[3]+(1/x[2]),x[4]-(1/x[2]))/beta(x[3],x[4])
vGBsin<-betamues$value  # valor verosimilitud

set.seed(10)
ma<- rGB2(10000, x[1], x[2],x[3],x[4])
giniGBsin<-gini(ma) # gini
#meanGBsin<-mean(ma)  # promedio 
par<-x
# quantiles 
qGBsin<-qGB2(probs, mu=par[1],sigma=par[2],nu=par[3],tau=par[4], lower.tail = TRUE, 
             log.p = FALSE) 


intGBsin<-c()
for ( i in 1: length(s))
{intGBsin[i]<-integrate(function(x) dGB2(x=x, mu=par[1], sigma=par[2],nu=par[3],tau=par[4]), 
                        lower=0, upper=s[i] )$v}



### caso restringido
ansbeta<-c()
ansbeta <- constrOptim.nl(par=theta, fn=fbeta2, heq=heq,hin=hin)
xGBcon<-ansbeta$par #estimadores de parámetros de distribición
x<-ansbeta$par
ExGBcon<-x[1]*beta(x[3]+(1/x[2]),x[4]-(1/x[2]))/beta(x[3],x[4])
ExGBcon<-round(ExGBcon,2)
if ((ci-ExGBcon)<1) 
ExGBcon<-ExGBcon else { asig<-paste("GB",años[j],sep='')
 theta<-as.matrix(na.omit(val[which(names(val)==asig)]))
 ansbeta<- constrOptim.nl(par=theta, fn=fbeta2, heq=heq,hin=hin)
 xGBcon<-ansbeta$par #estimadores de parámetros de distribición
 x<-ansbeta$par
 ExGBcon<-x[1]*beta(x[3]+(1/x[2]),x[4]-(1/x[2]))/beta(x[3],x[4])
      }

vGBcon<-ansbeta$value

set.seed(10)
ma<- rGB2(10000, x[1], x[2],x[3],x[4])
meanGBcon<-mean(ma)
giniGBcon<-gini(ma)
par<-x
qGBcon<-qGB2(probs, mu=par[1],sigma=par[2],nu=par[3],tau=par[4], lower.tail = TRUE, 
             log.p = FALSE)


# integrales para ver desigualdad de 10 y 90 porciento
q10<-qGB2(p=c(.10),mu=par[1],sigma=par[2],nu=par[3],tau=par[4])
q90<-qGB2(p=c(.90),mu=par[1],sigma=par[2],nu=par[3],tau=par[4])

den<-integrate(function(x) x*dGB2(x=x, mu=par[1], sigma=par[2],nu=par[3],tau=par[4]), 
               lower= 0, upper=q10 )
numa<-integrate(function(x) x*dGB2(x=x, mu=par[1], sigma=par[2],nu=par[3],tau=par[4]), 
                lower= 0, upper=q90 )
num<-ci-numa$v
resGB<-num/den$v
numGB<-num
denGB<-den$v


mGB<-max(muestra)
i<-integrate(function(x) x*dGB2(x=x, mu=par[1], sigma=par[2],nu=par[3],tau=par[4]), 
             lower=0, upper=mGB)
cocGB<-ci-i$v


intGBcon<-c()
for ( i in 1: length(s))
{intGBcon[i]<-integrate(function(x) dGB2(x=x, mu=par[1], sigma=par[2],nu=par[3],tau=par[4]), 
                        lower=0, upper=s[i] )$v}

par<-xGBcon
pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
qq<-c()
qq[1]<-0
for ( i in 2:12)
  
  qq[i]<-qGB2(p=pes[i],mu=par[1],sigma=par[2],nu=par[3],tau=par[4])
for ( i in 2:12)
  ii[i]<-integrate(function(x) x*dGB2(x=x, mu=par[1],sigma=par[2],nu=par[3],tau=par[4]), 
                   lower= qq[i-1], upper=qq[i])$v


ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dGB2(x=x, mu=par[1],sigma=par[2],nu=par[3],tau=par[4]), 
                              lower= qq[12], upper=100000000000)$v
ii<-ii[-1]
ii[13]<-sum(ii)
if (ii[13]==ExGBcon) 
  ii[13]<-ii[13] else 
  {ii[12]<-ExGBcon-sum(ii[-c(12,13)])
   ii[13]<-ExGBcon}
iiGB<-ii
perGB<-qq[12]  ## percentil .9999999 
t(t(ii))


## Distribución GAMA dos parámetros 


heq<- function(x)  #esperanza igual a constante 
{ h <- rep(NA, 1)
  h[1]<-ci-x[1]
  h
}

hin<-function(x)  ## desigualdades 
{
  h <- rep(NA, 2)
  h[1]<-x[1]
  h[2]<-x[2]
}


# caso no restringido 
modgaa<-gamlss(muestra~1, family=GA,weights=fac)
mu<-fitted(modgaa,"mu")[1]
sigma<-fitted(modgaa,"sigma")[1]
theta<-c(mu,sigma)
gamamues <- constrOptim.nl(par=theta, fn=fgama2, hin=hin)
xGsin<-gamamues$par
x<-c()
x<-gamamues$par

vGsin<-gamamues$value
ExGsin<-x[1]
set.seed(10)
ma<- rGA(10000, x[1], x[2])
meanGsin<-mean(ma)
giniGsin<-gini(ma)
par<-x
qGsin<-qGA(probs, mu=par[1],sigma=par[2], lower.tail = TRUE, 
           log.p = FALSE)


## caso restringido 
ansgama<-c()
ansgama <- constrOptim.nl(par=xGsin, fn=fgama2, heq=heq,hin=hin)
xGcon<-ansgama$par
x<-ansgama$par

ExGcon<-round(x[1],2)
if ((ci-ExGcon)<1) 
  ExGcon<-ExGcon else {asig<-paste("G",años[j],sep='')
 theta<-as.matrix(na.omit(val[which(names(val)==asig)]))
 ansgama<- constrOptim.nl(par=theta, fn=fgama2, heq=heq,hin=hin)
 xGcon<-ansgama$par
 x<-ansgama$par
 ExGcon<-x[1]}

vGcon<-ansgama$value

set.seed(10)
ma<- rGA(10000, x[1], x[2])
#meanGcon<-mean(ma)
giniGcon<-gini(ma)
par<-x
qGcon<-qGA(probs, mu=par[1],sigma=par[2], lower.tail = TRUE, 
           log.p = FALSE)

q10<-qGA(p=c(.10),mu=par[1],sigma=par[2])
q90<-qGA(p=c(.90),mu=par[1],sigma=par[2])

den<-integrate(function(x) x*dGA(x=x, mu=par[1], sigma=par[2]), 
               lower= 0, upper=q10 )
numa<-integrate(function(x) x*dGA(x=x, mu=par[1], sigma=par[2]), 
                lower= 0, upper=q90 )
num<-ci-numa$v
resG<-num/den$v
numG<-num
denG<-den$v

mG<-max(muestra)
i<-integrate(function(x) x*dGA(x=x, mu=par[1], sigma=par[2]), 
             lower=0, upper=mG)
cocG<-ci-i$v


par<-xGBcon
pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
qq<-c()
qq[1]<-0
for ( i in 2:12)
  
  qq[i]<-qGB2(p=pes[i],mu=par[1],sigma=par[2],nu=par[3],tau=par[4])
for ( i in 2:12)
  ii[i]<-integrate(function(x) x*dGB2(x=x, mu=par[1],sigma=par[2],nu=par[3],tau=par[4]), 
                   lower= qq[i-1], upper=qq[i])$v


ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dGB2(x=x, mu=par[1],sigma=par[2],nu=par[3],tau=par[4]), 
                              lower= qq[12], upper=100000000000)$v
ii<-ii[-1]
ii[13]<-sum(ii)
if (ii[13]==ExGBcon) 
  ii[13]<-ii[13] else 
  {ii[12]<-ExGBcon-sum(ii[-c(12,13)])
   ii[13]<-ExGBcon}
iiGB<-ii
perGB<-qq[12]  ## percentil .9999999 
t(t(ii))






############ distribución log-normal 


heq<- function(x)  #esperanza igual a constante 
{ h <- rep(NA, 1)
  h[1]<-ci-((exp(x[2]^2))^(0.5)*exp(x[1]))
  h
}

hin<-function(x)
{
  h <- rep(NA, 2)
  h[1]<-x[1]
  h[2]<-x[2]
  h
}

#caso no restringido 

log<-gamlss(muestra~1, family=LOGNO,weights=fac)
mu<-fitted(log,"mu")[1]
sigma<-fitted(log,"sigma")[1]
theta<-c(mu,sigma)

logmues <- constrOptim.nl(par=theta, fn=flog2, hin=hin)            
xLNsin<-logmues$par
x<-logmues$par
vLNsin<-logmues$value
ExLNsin<-(exp(x[2]^2))^(0.5)*exp(x[1])
set.seed(10)
ma<- rLOGNO(10000, x[1], x[2])
#meanLNsin<-mean(ma)
giniLNsin<-gini(ma)
par<-x
qLNsin<-qLOGNO(probs, mu=par[1],sigma=par[2], lower.tail = TRUE, 
               log.p = FALSE)

intLNsin<-c()
for ( i in 1: length(s))
  {intLNsin[i]<-integrate(function(x) dLOGNO(x=x, mu=par[1], sigma=par[2]), 
                        lower=0, upper=s[i] )$v}


## caso restringido
#anslog<-c()
anslog <- constrOptim.nl(par=xLNsin, fn=flog2, heq=heq,hin=hin)
xLNcon<-anslog$par
x<-anslog$par
ExLNcon<-(exp(x[2]^2))^(0.5)*exp(x[1])  ## esperanza
ExLNcon<-round(ExLNcon,2)
if ((ci-ExLNcon)<1) 
  ExLNcon<-ExLNcon else 
{asig<-paste("LN",años[j],sep='')
 theta<-as.matrix(na.omit(val[which(names(val)==asig)]))
 anslog<- constrOptim.nl(par=theta, fn=flog2, heq=heq,hin=hin)
 xLNcon<-anslog$par
 x<-anslog$par
 ExLNcon<-(exp(x[2]^2))^(0.5)*exp(x[1])
}

vLNcon<-anslog$value               ### verosimilitud

set.seed(10)
ma<- rLOGNO(10000, x[1], x[2])
#meanLNcon<-mean(ma)
giniLNcon<-gini(ma)
par<-x
qLNcon<-qLOGNO(probs, mu=par[1],sigma=par[2], lower.tail = TRUE, 
               log.p = FALSE)

q10<-qLOGNO(p=c(.10),mu=par[1],sigma=par[2])
q90<-qLOGNO(p=c(.90),mu=par[1],sigma=par[2])

den<-integrate(function(x) x*dLOGNO(x=x, mu=par[1], sigma=par[2]), 
               lower= 0, upper=q10 )
numa<-integrate(function(x) x*dLOGNO(x=x, mu=par[1], sigma=par[2]), 
                lower= 0, upper=q90 )
num<-ci-numa$v
resLN<-num/den$v
numLN<-num
denLN<-den$v

mLN<-max(muestra)
i<-integrate(function(x) x*dLOGNO(x=x, mu=par[1], sigma=par[2]), 
             lower= 0, upper=mLN )
cocLN<-ci-i$v


intLNcon<-c()
for ( i in 1: length(s))
{intLNcon[i]<-integrate(function(x) dLOGNO(x=x, mu=par[1], sigma=par[2]), 
                        lower=0, upper=s[i] )$v}

par<-xLNcon

pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
qq<-c()
qq[1]<-0
for ( i in 2:12)
  qq[i]<-qLOGNO(p=pes[i],mu=par[1],sigma=par[2])
ii<-c()
for ( i in 2:12)
  ii[i]<-integrate(function(x) x*dLOGNO(x=x, mu=par[1],sigma=par[2]), 
                   lower= qq[i-1], upper=qq[i])$v


ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dLOGNO(x=x, mu=par[1],sigma=par[2]), 
                              lower= qq[12], upper=1000000000)$v
ii<-ii[-1]
ii[13]<-sum(ii)
if (ii[13]==ExLNcon) 
  ii[13]<-ii[13] else 
  {ii[12]<-ExLNcon-sum(ii[-c(12,13)])
   ii[13]<-ExLNcon}
iiLN<-ii
perLN<-qq[12]  ## percentil .9999999 
t(t(iiLN))



## distribucion gamma generalizada 3 parametros 


heq<- function(x)  #esperanza igual a constante 
{ 
  h <- rep(NA, 1)
  t<-1/(x[2]^2*x[3]^2)
  h[1]<-ci-(x[1]*gamma(t+(1/x[3]))/(t^(1/x[3])*gamma(t)))
  
  h
  
}

hin<-function(x)
{
  h <- rep(NA, 2)
  h[1]<-x[1]
  h[2]<-x[2]
  h
}

#samples<-sample(muestra,.10*length(muestra))
#caso no restringido
gammues<-gamlss(muestra~1, family=GG,weights=fac)
gmag<-gammues
#refit(gmag)
mu<-fitted(gmag,"mu")[1]
sigma<-fitted(gmag,"sigma")[1]
nu<-fitted(gmag,"nu")[1]
theta<-c(mu,sigma,nu)

gammues<-c()
x<-c()
gammues<-constrOptim.nl(par=theta, fn=fgamG, hin=hin)
xGGsin<-gammues$par
x<-gammues$par
t<-1/(x[2]^2*x[3]^2)
ExGGsin<-(x[1]*gamma(t+(1/x[3])))/(t^(1/x[3])*gamma(t))
vGGsin<-gammues$value

set.seed(10)
ma<- rGG(1000, x[1], x[2],x[3])
#meanGGsin<-mean(ma)
giniGGsin<-gini(ma)
par<-x
qGGsin<-qGG(probs, mu=par[1],sigma=par[2], nu=par[3],lower.tail = TRUE, 
            log.p = FALSE)


### caso restringido 
ansgam<-c()
asig<-paste("GG",años[j],sep='')
theta<-as.matrix(na.omit(val[which(names(val)==asig)]))
ansgam<- constrOptim.nl(par=theta, fn=fgamG, heq=heq,hin=hin)
xGGcon<-ansgam$par
x<-ansgam$par
t<-1/(x[2]^2*x[3]^2)
ExGGcon<-x[1]*gamma(t+(1/x[3]))/(t^(1/x[3])*gamma(t))
if ((ci-ExGGcon)<1) 
  ExGGcon<-ExGGcon  else  
    {asig<-paste("GG",años[j],sep='')
   theta<-as.matrix(na.omit(val[which(names(val)==asig)]))
   ansgam<- constrOptim.nl(par=theta, fn=fgamG, heq=heq,hin=hin)
   xGGcon<-ansgam$par
   x<-ansgam$par
   t<-1/(x[2]^2*x[3]^2)
   ExGGcon<-x[1]*gamma(t+(1/x[3]))/(t^(1/x[3])*gamma(t))
  }


vGGcon<-ansgam$value
set.seed(10)
ma<- rGG(10000, x[1], x[2],x[3])
#meanGGcon<-mean(ma)
giniGGcon<-gini(ma)

par<-x
qGGcon<-qGG(probs, mu=par[1],sigma=par[2], nu=par[3],lower.tail = TRUE, 
            log.p = FALSE)

q10<-qGG(p=c(.10),mu=par[1],sigma=par[2],nu=par[3])
q90<-qGG(p=c(.90),mu=par[1],sigma=par[2],nu=par[3])

den<-integrate(function(x) x*dGG(x=x, mu=par[1], sigma=par[2],nu=par[3]), 
               lower= 0, upper=q10)
numa<-integrate(function(x) x*dGG(x=x, mu=par[1], sigma=par[2],nu=par[3]), 
                lower= 0, upper=q90 )
num<-ci-numa$v
resGG<-num/den$v
numGG<-num
denGG<-den$v

mGG<-max(muestra)

i<-integrate(function(x) x*dGG(x=x, mu=par[1], sigma=par[2],nu=par[3]), 
             lower= 0, upper=mGG )
cocGG<-ci-i$v

par<-xGGcon
pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
qq<-c()
qq[1]<-0
ii<-c()
for ( i in 2:length(pes))
  qq[i]<-qGG(p=pes[i],mu=par[1],sigma=par[2],nu=par[3])

for ( i in 2:12)
  ii[i]<-integrate(function(x) x*dGG(x=x, mu=par[1],sigma=par[2],nu=par[3]), 
                   lower= qq[i-1], upper=qq[i])$v


ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dGG(x=x, mu=par[1],sigma=par[2],nu=par[3]), 
                              lower= qq[12], upper=100000000000)$v
ii<-ii[-1]
ii[13]<-sum(ii)
if (ii[13]==ExGGcon) 
  ii[13]<-ii[13] else 
  {ii[12]<-ExGGcon-sum(ii[-c(12,13)])
   ii[13]<-ExGGcon}
iiGG<-ii
perGG<-qq[12]  ## percentil .9999999 

t(t(iiGG))




#################### distribución pareto #################
########################################################
############################################################

library(actuar)


hin<-function(x)
{
  h <- rep(NA, 2)
  h[1]<-x[1]-1
  h[2]<-x[2]
  h
}


heq<- function(x)  #esperanza igual a constante 
{  h <- rep(NA, 1)
   h[1]<-ci-(mpareto(1,x[1],x[2]))
   h
}


fparet2<-function(x,m=muestra,w=fac) ## funcion de verosi
{ 
  fi<-c()
  fi<- -1*w*(dpareto(m, x[1],x[2],log=TRUE))
  sumaf<-sum(fi)
  return(sumaf)
  
}


parmues <- constrOptim.nl(par=theta, fn=fparet2, hin=hin)            
xPsin<-parmues$par
x<-parmues$par
vPsin<-parmues$value
ExPsin<-mpareto(1,x[1],x[2])
set.seed(10)
ma<- rpareto(10000, x[1], x[2])
#meanLNsin<-mean(ma)
giniPsin<-gini(ma)
par<-x
qPsin<-qpareto(probs, shape=par[1],scale=par[2])

intPsin<-c()
for ( i in 1: length(s))
{intPsin[i]<-integrate(function(x) dpareto(x=x, shape=par[1],scale=par[2]), 
                        lower=0, upper=s[i] )$v}


## caso restringido

anspar<- constrOptim.nl(par=xPsin, fn=fparet2, heq=heq,hin=hin)
xPcon<-anspar$par
x<- xPcon
ExPcon <-mpareto(1,x[1],x[2])  ## esperanza
ExPcon<-round(ExPcon,2)

vPcon<-anspar$value               ### verosimilitud

set.seed(10)
ma<- rpareto(10000, x[1], x[2])
#meanLNcon<-mean(ma)
giniPcon<-gini(ma)
par<-x
qPcon<-qpareto(probs, shape=par[1],scale=par[2])

q10<-qpareto(p=c(.10),shape=par[1],scale=par[2])
q90<-qpareto(p=c(.90),shape=par[1],scale=par[2])

den<-integrate(function(x) x*dpareto(x=x, par[1],par[2]), 
               lower= 0, upper=q10 )
numa<-integrate(function(x) x*dpareto(x=x, par[1], par[2]), 
                lower= 0, upper=q90 )
num<- ExPcon-numa$v
resP<-num/den$v
numP<-num
denP<-den$v

mP<-max(muestra)
i<-integrate(function(x) x*dpareto(x=x, par[1], par[2]), 
             lower= 0, upper=mP )
cocP<- ExPcon-i$v


intPcon<-c()
for ( i in 1: length(s))
  
{intPcon[i]<-integrate(function(x) dpareto(x=x,par[1], par[2]), 
                       lower=0, upper=s[i] )$v}



par<-xPcon
pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
qq<-c()
qq[1]<-0
for ( i in 2:12)
  #q[i]<-qGG(p=pes[i],mu=par[1],sigma=par[2],nu=par[3])
  qq[i]<-qpareto(p=pes[i],par[1],par[2])
ii<-c()
for ( i in 2:12)
  ii[i]<-integrate(function(x) x*dpareto(x=x, par[1], par[2]), 
                   lower= qq[i-1], upper=qq[i])$v


ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dpareto(x=x, par[1], par[2]), 
                              lower= qq[12], upper=1000000000)$v
ii<-ii[-1]
ii[13]<-sum(ii)
if (ii[13]==ExPcon) 
  ii[13]<-ii[13] else 
  {ii[12]<-abs(ExPcon-sum(ii[-c(12,13)]))
   ii[13]<-ExPcon}
iiP<-ii
perP<-qq[12]  ## percentil .9999999 
t(t(iiP))

##### fin de ajuste de distribuciones 
#### datos de la muestra 
if (length(w)==0)
{data.desp <- svydesign(id=~upm,weights=~factor_hog,data=ent)} else
{data.desp <- svydesign(id=~upm[-w],weights=~factor_hog[-w],data=ent)
}
media<-svymean(~muestra,data.desp)
quantdis<- svyquantile(~muestra, data.desp, probs,ci=TRUE)


##### salida de información#############


a1<-c(vPcon,vLNcon,vGGcon,vGBcon)
a2<-c(vPsin,vLNsin,vGGsin,vGBsin)
a3<-c(giniPsin,giniLNsin,giniGGsin,giniGBsin)
#a3<-c(giniGsin$value,giniLNsin$value,giniGGsin$value,giniGBsin$value)/100
a4<-c(as.numeric(ExPsin),as.numeric(ExLNsin),as.numeric(ExGGsin),as.numeric(ExGBsin))
av1<-c(0,0 ,0 ,0)
di<-c("-","-","-","-")
d<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),4,4)
d[1:2,1]<-t(t(xPsin))
d[1:2,2]<-t(t(xLNsin))
d[1:3,3]<-t(t(xGGsin))
d[,4]<-t(t(xGBsin))
av2<-c(0,0 ,0 ,0)
a5<-cbind(qPsin,qLNsin,qGGsin,qGBsin)

salida1<-(rbind(a1,a2,a3,a4,di,d,di,a5))
colnames(salida1)<-c("P","L-N","GG","GB2")
rownames(salida1)<-c("Verosi-restringido","Verosi-no restingido","Gini-SR",
                     "Esperanza-SR","parámetros estimados-SR","mu","sigma","nu","tau","Percentiles-SR","0.1"," 0.5", "1", "2", "5", "10","20","30","40", "50","60","70","80","90","95","98","99","99.5","99.9")
#write.csv(salida1,"salida1.csv")
b3<-c(giniPcon,giniLNcon,giniGGcon,giniGBcon)
#b3<-c(giniGcon$value,giniLNcon$value,giniGGcon$value,giniGBcon$value)/100
b4<-c(as.numeric(ExPcon),as.numeric(ExLNcon),as.numeric(ExGGcon),as.numeric(ExGBcon))
bv1<-c(0,0 ,0 ,0)
d2<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),4,4)
d2[1:2,1]<-t(t(xPcon))
d2[1:2,2]<-t(t(xLNcon))
d2[1:3,3]<-t(t(xGGcon))
d2[1:4,4]<-t(t(xGBcon))

b5<-cbind(qPcon,qLNcon,qGGcon,qGBcon)
b6<-c(numP,numLN,numGG,numGB)
b7<-c(denP,denLN,denGG,denGB)
b8<-c(resP,resLN,resGG,resGB)
b9<-c(cocP,cocLN,cocGG,cocGB)
b10<-b9/b4
b11<-cbind(iiP,iiLN,iiGG,iiGB) ## integrales 
b12<-c(perP,perLN,perGG,perGB)
mG<-max(muestra)
vm<-c(mG,"-","-","-") # valor mayor de la muestra
promi<-c(matrix (media)[1],"-","-","-") # promedio enigh
sdi<-c(as.numeric(cbind(c(data.frame(media)[2]))),"-","-","-") # desviacion del promediop
tem1<-c(tote,"-","-","-")   ## total de ingreso CN trimestre
tem2<-c(tot_hog,"-","-","-")  ## total hogares 
totencu<-c(sum(muestra*fac),"-","-","-")  ## total ingreso segun encuesta 

salida2<-rbind(b3,b4,di,d2,di,b5,b12,b6,b7,b8,di,b9,b10,vm,promi,sdi,tem1,totencu,tem2)

colnames(salida2)<-c("P","L-N","GG","GB2")
rownames(salida2)<-c("Gini-CR",
                     "Esperanza-CR","Parámetros estimados-CR","mu","sigma","nu","tau","Percentiles-CR","0.1"," 0.5", "1", "2", "5", "10","20","30","40", "50","60","70","80","90","95","98","99","99.5","99.9",
                     "99.9999","numerador","denominador","cociente","Relación acumulados mayores al máximo observado a ingresos totales","valor","cociente",
                     "valor máximo de la muestra","promedio ENIGH","sd promedio","tot_SCNM","tot_encuesta","Hogares"
)

### probabilidades mj
c1<-cbind(intGBsin,intGBcon,intLNsin,intLNcon)
colnames(c1)<-c("P(x<mj)GB2_sin_restrición","P(x<mj)GB2_con_restrición","P(x<mj)LN_sin_restrición","P(x<mj)LN_con_restrición")
s<-as.character(s)
rownames(c1)<-s

#write.csv(c1,"pxmj.csv")

## percentiles muestra expandida 

ddd<-cbind(t(quantdis$qua),t(matrix(quantdis$CI,nrow=2)))
#ddd<-cbind(t(t(quantdis)),cbind(t(t(quantdis)),t(t(quantdis))))
ddd<-cbind(ddd,rep(0,19))
colnames(ddd)<-c("percentil","límite inferior","límite superior","-")
rownames(ddd)<-c("0.1","0.5", "1", "2", "5", "10","20","30","40", "50","60","70","80","90","95","98","99","99.5","99.9")

salida<-paste("salida",i,sep="")
integrales<-b11
rownames(integrales)<-c(
  "0---10",
  "10---20",
  "20---30",
  "30---40",
  "40---50",
  "50---60",
  "60---70",
  "70---80",
  "80---90",
  "90---99",
  "99---999999",
  "999999--Inf",
  "suma")


salida3<-(rbind(salida1,c("-","-","-","-"),c("-","-","-","-"),c("P","L-N","GG","GB2"),salida2,c("-","-","-","-"),c("-","-","-","-"),c("percentil","límite inferior","límite superior","-"),ddd,c("-","-","-","-"),c("-","-","-","-"),c("P_GB2_SR","P_GB2_CR","PLN_SR","PLN_CR"),c1,c("-","-","-","-"),c("-","-","-","-"),c("iP","iLN","iGG","iGB"),integrales))

carpeta<-paste("D:","//","survey","//",años[j],"//",ei,sep="")
write.csv(salida3,paste(carpeta,"//","salida_proMCS_datosMCS",ei,".","csv",sep=""))

}




## declaracion de la FUNCIóN GINI
# se puede usar tambien la libreria laeken tiene una funcion "gini"###

gini <- function(x, unbiased = TRUE, na.rm = FALSE){
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind <- is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x <- x[!na.ind]
  n <- length(x)
  mu <- mean(x)
  N <- if (unbiased) n * (n - 1) else n * n
  ox <- x[order(x)] ## quitar la informacion irrelevante DROP 
  ## la función crossprod = productos cruzados 
  dsum <- drop(crossprod(2 * 1:n - n - 1,  ox)) 
  dsum / (mu * N)
}


fbeta2<-function(x,m=muestra,w=fac) ## funciÓn de verosimilitud
{ 
  fi<-c()
  fi<- -1*w*dGB2(m, x[1], x[2],x[3],x[4],log=TRUE)
  sumaf<-sum(fi)
  return(sumaf)
}



fgama2<-function(x,m=muestra,w=fac) ## funcion de verosi
{ 
  fi<-c()
  fi<- -1*w*(dGA(m, x[1], x[2],log=TRUE))
  sumaf<-sum(fi)
  return(sumaf)
 }

flog2<-function(x,m=muestra,w=fac) ## funcion de verosi
{ 
  fi<-c()
  fi<- -1*w*(dLOGNO(m, x[1], x[2],log=TRUE))
  sumaf<-sum(fi)
  return(sumaf)
  
}


fgamG<-function(x,m=muestra,w=fac) ## funcion de verosi
{ 
  fi<-c()
  fi<- -1*w*dGG(m, x[1], x[2],x[3],log=TRUE)
  sumaf<-sum(fi)
  return(sumaf)
}




####visualize the densities
income<-muestra
u=seq(0,2e5,length=251)
svyhist(~income,data.desp,breaks=seq(0,max(muestra)+100000,by=2000),col=rgb(0,0,1,.5),border="white",
        xlim=c(0,2e5),probability=TRUE,main="Histogam of income",xlab="income")

v_b <- dGB2(u, xGBsin[1],xGBsin[2],xGBsin[3],xGBsin[4])
v_br <- dGB2(u, xGBcon[1],xGBcon[2],xGBcon[3],xGBcon[4])
v_ga <- dGA(u, xGsin[1],xGsin[2])
v_gar <- dGA(u, xGcon[1],xGcon[2])
v_ln <- dLNO(u, xLNsin[1],xLNsin[2])
v_lnr <- dLNO(u, xLNcon[1],xLNcon[2])
v_gg <- dGG(u, xGGsin[1],xGGsin[2],xGGsin[3])
lines(u,v_b,col="red",lwd=2)
lines(u,v_ga,col="black",lwd=2)             
lines(u,v_ln,col="blue",lwd=1)   
lines(u,v_gg,col="pink",lwd=3)   
lines(u,v_,col="black",lwd=2)

lines(u,v_br,col="blue",lwd=2)
lines(u,v_gar,col="black",lwd=2)             
lines(u,v_lnr,col="blue",lwd=1)  


cdf.mue<-svycdf(~muestra,data.desp)
plot(cdf.mue,col="black",xlim=c(0,400000))
v_b2 <- pGB2(u, xGBsin[1],xGBsin[2],xGBsin[3],xGBsin[4])
v_br2 <- pGB2(u, xGBcon[1],xGBcon[2],xGBcon[3],xGBcon[4])
lines(u,v_b2,col="red",lwd=2)
lines(u,v_br2,col="blue",lwd=2)
 lines(u,v_ln,col=rgb(1,0,0,.4),lwd=2)


## curvas de lorenz
library("ineq")
plot(Lc(ma))
 lines(Lc.lognorm,param=1.5,col="red")
 lines(Lc.lognorm,param=1.2,col="red")
 lines(Lc.lognorm,param=.8,col="red")



##### obtener algunas integrales especiales

lis<-list()
for ( j in 1:32)
{
  info<-paste ("salida_proMCS_datosMCS",j,sep="")
  carpeta<-paste("D:","//","survey","//",2012,"//",j,"//",info,sep="")
  
  data<-read.csv(paste(carpeta,".csv",sep=""))
  data<-as.matrix(data)
  ii<-c()
  q<-c()
  pes<-c( 0, 0.1,.2,.3,.4,.5,.6,.7,.8,.9,.99,.999999)
  par<-data[36:37,3]
  par<-as.numeric(par)
 # par<-c(10.3217
  #       , 1.23881
         
         
  #) 
  q[1]<-0
  for ( i in 2:12)
    #q[i]<-qGG(p=pes[i],mu=par[1],sigma=par[2],nu=par[3])
 q[i]<-qGB2(p=pes[i],mu=par[1],sigma=par[2],nu=par[3],tau=par[4])
  for ( i in 2:12)
    ii[i]<-integrate(function(x) x*dLOGNO(x=x, mu=par[1], sigma=par[2]), 
                     lower= q[i-1], upper=q[i])$v
  
  
  ii[i+1]<-  ii[i+1]<-integrate(function(x) x*dLOGNO(x=x, mu=par[1], sigma=par[2]), 
                                lower= q[12], upper=1000000000)$v
  ii<-ii[-1]
  ii[13]<-sum(ii)
  t(t(ii))
  lis[[j]]<-t(t(ii))
}


lis<-list()
for ( j in 1:32)
{
  info<-paste ("salida_proMCS_datosMCS",j,sep="")
  carpeta<-paste("D:","//","survey","//",2012,"//",j,"//",info,sep="")
  
  data<-read.csv(paste(carpeta,".csv",sep=""))
  data<-as.matrix(data)
  promedio<-as.numeric(data[68,2])
  restricc<-as.numeric(data[34,3])
  maximo<-as.numeric(data[67,2])
  ginis<-as.numeric(data[33,3])
  cociente<-as.numeric(data[63,3])
 tem<- cbind(promedio,restricc,maximo,ginis,cociente)
   tem<-t(t(tem))
  lis[[j]]<-list(tem)

}
