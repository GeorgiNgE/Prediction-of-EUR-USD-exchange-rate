
# Libraries ---------------------------------------------------------------
library(psych)
library(tsoutliers)
library(sandwich)
library(strucchange)
library(tseries)

# Prelucrare date --------------------------------------------------------------------

date<-read.csv("date2.csv",stringsAsFactors = FALSE)
date[,"Date"]<-as.Date(date[,"Date"],"%m/%d/%Y")
nw <- length(date$Price)
difw <- log(date$Price[-1]/date$Price[-nw])

# Descriere date ----------------------------------------------------------
JarqueBera.test(date$Price) #nu este normal distribuita
JarqueBera.test(difw) #nu este normal distribuita

statsw<-describe(date$Price)
stats<-describe(difw)
stats_af<-as.matrix(c(statsw$mean,statsw$sd,statsw$min,statsw$max,statsw$skew,statsw$kurtosis))
stats_af<-cbind(stats_af,c(stats$mean,stats$sd,stats$min,stats$max,stats$skew,stats$kurtosis))
colnames(stats_af)<-c("Statistici descriptive curs de schimb","Statistici descriptive randamente")
rownames(stats_af)<-c("Medie","Abatere standard","Minim","Maxim","Skewness","Kurtosis")
round(stats_af,2)


# Afisare -----------------------------------------------------------------

#plot-ul cu pretul si trasarea mediei
par(mfrow=c(2,2))
plot(date,type="l",col=4)
abline(h=val$mean,col="dark cyan")
boxplot(date$Price,col=4,ylab = "EUR/USD")
hist(date$Price,xlab="EUR/USD",col=4,main="Histograma EUR/USD")
qqnorm(date$Price,cex=.5)
qqline(date$Price,ylab ="EUR/USD",col=4,lwd=3)

dev.new()
par(mfrow=c(1,2))
hist(difw,xlab="Randamente",col=4,main="Histograma EUR/USD")
plot(date[-1,"Date"],difw,type="l",col=4,xlab="Data",ylab = "Randamente")

# Schimbari structurale ---------------------------------------------------
bp_ts_con<-breakpoints(date$Price~1)
plot(bp_ts_con)
summary(bp_ts_con) 
seg1<-date[1:779,]
seg2<-date[780:2852,]
seg3<-date[2853:4523,]

#grafic al datelor cu schimbarile din structura
dev.new()
par(mfrow=c(1,2))
plot(date,type='l',main="Setul de date cu schimbarile structurale")
lines(fitted(bp_ts_con, breaks = 2) ~date$Date , col = 4, lwd = 3)
plot(bp_ts_con,main="Valoarea criteriilor pentru numarul de segmente")


# MF-DFA ------------------------------------------------------------------
n <- length(date$Price)
difw <- log(date$Price[-1]/date$Price[-n])
#seria initiala
scale=10:904 #bine
q<- -5:5
m<-1
bw<-MFDFA(difw, scale, m, q)
mdw<-max(bw$spec$hq)-min(bw$spec$hq)

dev.new()
par(mai=rep(1, 4))
par(mfrow=c(1,2))
##  q-order Hurst exponent
plot(q, bw$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent", ylim=c(min(bw$Hq),max(bw$Hq)))
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)

## Multifractal spectrum
plot(bw$spec$hq, bw$spec$Dq, col=1, axes=F, pch=16, main="Multifractal spectrum",
     ylab=expression('D'[q]),cex.lab=1.8, cex.axis=1.8,
     xlab=expression('h'[q]))
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
x1=bw$spec$hq
y1=bw$spec$Dq
rr<-poly_fit(x1,y1,4)
mm1<-rr$model1
mm<-rr$polyfit
x2<-seq(0,max(x1)+1,0.01)
curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
lines(x2,curv, col="red", lwd=2)

#histograma hurst
dev.new()
par(mai=rep(1, 4))
par(mfrow=c(1,1))
hist(bw$Hq,xlab="Valorile exponentului Hurst generalizat",col=4, main= "Repartitia exponentului Hurst")

#seg1
n <- length(seg1$Price)
dif <- log(seg1$Price[-1]/seg1$Price[-n])
scale=10:18 # a iesit bine asa
q<- -5:5
m<-1
b1<-MFDFA(dif, scale, m, q)
md1<-max(b1$spec$hq)-min(b1$spec$hq)

#seg2
n <- length(seg2$Price)
dif <- log(seg2$Price[-1]/seg2$Price[-n])
scale=10:415 #binisor 
q<- -5:5
m<-1
b2<-MFDFA(dif, scale, m, q)
md2<-max(b2$spec$hq)-min(b2$spec$hq)

#seg3
n <- length(seg3$Price)
dif <- log(seg3$Price[-1]/seg3$Price[-n])
scale=10:27 #binisor 27
q<- -5:5
m<-1
b3<-MFDFA(dif, scale, m, q)
md3<-max(b3$spec$hq)-min(b3$spec$hq)


dev.new()
par(mai=rep(1, 2))
par(mfrow=c(3,2))
##  q-order Hurst exponent
plot(q, b1$Hq, col=1, axes= F, ylab="h(q)", pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent- segment1", ylim=c(min(b1$Hq),max(b1$Hq)))
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)

## Multifractal spectrum
plot(b1$spec$hq, b1$spec$Dq, col=1, axes=F, pch=16, main="Multifractal spectrum- segment1",
     ylab="D(q)",cex.lab=1.8, cex.axis=1.8,
     xlab="h(q)")
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
x1=b1$spec$hq
y1=b1$spec$Dq
rr<-poly_fit(x1,y1,4)
mm1<-rr$model1
mm<-rr$polyfit
x2<-seq(0,max(x1)+1,0.01)
curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
lines(x2,curv, col="red", lwd=2)

##  q-order Hurst exponent
plot(q, b2$Hq, col=1, axes= F, ylab="h(q)", pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent- segment2", ylim=c(min(b2$Hq),max(b2$Hq)))
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)

## Multifractal spectrum
plot(b2$spec$hq, b2$spec$Dq, col=1, axes=F, pch=16, main="Multifractal spectrum- segment2",
     ylab="D(q)",cex.lab=1.8, cex.axis=1.8,
     xlab="h(q)")
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
x1=b2$spec$hq
y1=b2$spec$Dq
rr<-poly_fit(x1,y1,4)
mm1<-rr$model1
mm<-rr$polyfit
x2<-seq(0,max(x1)+1,0.01)
curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
lines(x2,curv, col="red", lwd=2)

##  q-order Hurst exponent
plot(q, b3$Hq, col=1, axes= F, ylab="h(q)", pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent- segment3", ylim=c(min(b3$Hq),max(b3$Hq)))
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)

## Multifractal spectrum
plot(b3$spec$hq, b3$spec$Dq, col=1, axes=F, pch=16, main="Multifractal spectrum- segment3",
     ylab="D(q)",cex.lab=1.8, cex.axis=1.8,
     xlab="h(q)")
grid(col="midnightblue")
axis(1, cex=4)
axis(2, cex=4)
x1=b3$spec$hq
y1=b3$spec$Dq
rr<-poly_fit(x1,y1,4)
mm1<-rr$model1
mm<-rr$polyfit
x2<-seq(0,max(x1)+1,0.01)
curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
lines(x2,curv, col="red", lwd=2)
rm(vm)

vm<-as.matrix(c(mdw,md1,md2,md3))
colnames(vm)<-"Gradul de multifractalitate"
rownames(vm)<-c("Set intreg","Prima perioada","A doua perioada","A treia perioada")
vm


