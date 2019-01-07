## Package yang dibutuhkan dan User-Defined Function
library(ggplot2)
library(forecast)
library(fBasics)
library(tseries)
library(plotly)
library(DT)
library(car)
library(knitr)
library(lattice)

## Data Import
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
dta=read.csv("D:/Training - StatTalk 20180303/data.csv",sep=",")[,-1]

##Change Data Frame to Time-Series Training Data
Y=ts(dta$Y[-c(89,90)],start=c(2009,2),frequency =12)
dY=diff(ts(Y,start=c(2010,2),frequency =12),1)
head(Y)

X1=ts(dta$X1[-c(89,90)],start=c(2009,2),frequency =12)
dX1=diff(ts(X1,start=c(2010,2),frequency =12),1)
head(X1)

X2=ts(dta$X2[-c(89,90)],start=c(2009,2),frequency =12)
dX2=diff(ts(X2,start=c(2010,2),frequency =12),1)

## Prosedur ARIMAX Model
## 1. Eksplorasi Data
## 2. Regression Model
## 3. Checking Regression Error
## 4. ARIMA Model

## 1. Eksplorasi Data
ggplotly(autoplot(Y)+geom_line(col="#3DA087",size=1.4)+
  geom_point(size=2,col="#E48C92")+theme_classic())
  
## Plot ACF dari variabel Y
acf(Y,plot=T)
pacf(Y,plot=T)

## Membuat Variabel dummy
dummy_2013=c(rep(0,48),rep(1,12),rep(0,30))
dta1=dta[,-5]
dta2=ts(cbind(dta1,dummy_2013),start=c(2010,2),frequency =12)
dta1=ts(cbind(dta1,dummy_2013)[-c(89,90),],frequency =12)
dummy_2013=c(rep(0,48),rep(1,12),rep(0,28))

## Plot Variabel Exogenous
ggplotly(autoplot(X1)+geom_line(col="#3DA087",size=1.4)+
  geom_point(size=2,col="#E48C92")+theme_classic())
  
ggplotly(autoplot(X2)+geom_line(col="#3DA087",size=1.4)+
  geom_point(size=2,col="#E48C92")+theme_classic())

## UJI ADF untuk melihat apakah data stasioner atau tidak
apply(dta,2,adf.test)

## UJI ADF setelah data di difference
apply(apply(dta,2,diff,1),2,adf.test)

## 2. Regression Model
vif(lm(Y~X1+X2+dummy_2013))

## Kemudian model regresinya seperti dibawah ini
printmodel = function(fit_arima,digits=4,se = T) {

  if (length(fit_arima$coef)>0) {
    cat("\nCoefficients : \n")
    coef = round(fit_arima$coef, digits=digits)
    if (se && nrow(fit_arima$var.coef)) {
      ses  = rep(0,length(coef))
      ses[fit_arima$mask] = round(sqrt(diag(fit_arima$var.coef)),digits=digits)
      coef = matrix(coef,1,dimnames=list(NULL,names(coef)))
      coef = rbind(coef,s.e.=ses)
      statt = coef[1,]/ses
      pval = 2*pt(abs(statt),df=length(fit_arima$residual)-1,lower.tail=F)
      coef = rbind(coef,t=round(statt,digits=digits),sign.=round(pval,digits=digits))
      row.names(coef) = c("estimate","standard error","t-value","p-value")
      coef = t(coef)
    }
    d=eval(fit_arima$call[[3]])[2]
    if(d>0){
    y=diff(eval(fit_arima$call[[2]]),d)
    }else{
      y=eval(fit_arima$call[[2]])
    }
    
    SSTo = sum((y-mean(y))^2)
    SSE = sum(fit_arima$residuals^2)
    `R-Squared`= paste(round(100*(1-(SSE/SSTo)),2),'%',sep='')
    
    
    if (fit_arima$call[1] == "arimax()"){
      npar=length(fit_arima$coef)
      nstar= length(fit_arima$residuals)-fit_arima$arma[6]-fit_arima$arma[7]*fit_arima$arma[5]
      bic=fit_arima$aic+npar*log((nstar)-2)
      aicc=fit_arima$aic+2*npar*((nstar)/(nstar-npar-1)-1)
      
      goodness = matrix(c(fit_arima$sigma2,fit_arima$loglik,fit_arima$aic,aicc,bic),5,1,dimnames = 
                          list(c("Standard Error of Model","Log-Likelihood","AIC","AICC","BIC")," "))
      
    } else {
    goodness = matrix(c(fit_arima$sigma2,fit_arima$loglik,fit_arima$aic,fit_arima$aicc,fit_arima$bic),5,1,
                      dimnames = list(c("Standard Error of Model","Log-Likelihood","AIC","AICC","BIC")," "))
    }
    print.default(coef,print.gap = 2)
  
    print.default(rbind(round(goodness,3),`R-Squared`),print.gap = 2,quote = F)
    
  }
  }
  
  fit1=Arima(Y,order=c(0,0,0),include.mean = F,include.drift=F,xreg=dta1[,2:4])
printmodel(fit1)

## Membuat residual dari model regresi
resreg1=residuals(fit1,"regression")

ggljbox=function(res,nlag=length(res)/4) {

pval <- numeric(nlag)
for (i in 1L:nlag)  {
  pval[i]=Box.test(res, i, type = "Ljung-Box")$p.value
}
lb<- data.frame(Lag = 1:nlag, `p value` = pval, Bound=0.05)

colnames(lb) <- c("Lag", "p value", "Bound")

ggplot(data = lb,aes_string(x="Lag")) + 
  geom_point(aes_string(y = "`p value`"),colour="#22A7F0",size=3) +
  geom_hline(yintercept = 0.05,colour="red",size=1.5)+
  scale_y_continuous(limits = c(0, 1)) + 
  ggtitle("p values for Ljung-Box statistic")+theme_classic()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=12),plot.title = 
          element_text(hjust=0.5,face="bold",size=20))
}

ggljbox(resreg1)

# 3. ARIMA Model

## Spesifikasi Model ARIMA
acf(resreg1, plot=T)
pacf(resreg1, plot=T)

acf(diff(resreg1,1), plot=T)
pacf(diff(resreg1,1), plot=T)

nsdiffs(Y,12)

auto.arima(Y,d=1,max.p=2,max.q=3,max.P=1,max.Q =1, ic='bic',trace=T,allowdrift = F,allowmean = F,xreg=dta1[,2:4])

## Estimasi Model ARIMA
fit11=Arima(Y,order=c(0,1,0),seasonal=c(0,0,1),
            include.mean = F,include.drift=F,xreg=dta1[,2:4])
printmodel(fit11)

## Evaluasi Model
ggljbox(residuals(fit11,"innovation"))

## Terakhir adalah  akurasi data testing
f1=forecast(fit11,xreg=dta2[c(89,90),2:4],level = c(95))
accuracy(f1,dta2[c(89,90),1])[2,c(2,3,5)]

fit21=Arima(Y,order=c(0,1,0),seasonal=c(0,0,1),
            include.mean = F,include.drift=F)
printmodel(fit21)

f2=forecast(fit21,level = c(95))
accuracy(f2,dta2[c(89,90),1])[2,c(2,3,5)]