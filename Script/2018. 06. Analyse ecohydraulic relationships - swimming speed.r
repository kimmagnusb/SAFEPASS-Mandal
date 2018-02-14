
rm(list=ls())
options(digits=10)

library(nlme)
library(geoR)

### FUNCTIONS

Transform2Normal.2p<-function(Ori){
T.box = boxcoxfit(Ori, lambda2=TRUE)
lambda = T.box$lambda[1]
lambda2 = T.box$lambda[2]
if(lambda==0){T.norm = log(Ori + lambda2)}
if(lambda!=0){T.norm = ((Ori + lambda2) ^ lambda - 1) / lambda}
return(T.norm)
}

### READ AND CLEAN DATA

DIR<-"M:\\My Documents\\NINA projects\\Safepass\\Data2018\\Processed\\"

D <- read.csv(paste(DIR,"NewSmoltData for LME.2018-02-07.csv",sep=""))
nrow(D)
tf<-colnames(D)%in%c("AngDif.C","SwimSpeed.C","u","v","w","ke","Len","FlowSpeed","FF","Tag")
D<-D[,tf]
D<-na.omit(D)
nrow(D)

D$AngDif.C.norm<-Transform2Normal.2p(D$AngDif.C)
D$SwimSpeed.C.norm<-Transform2Normal.2p(D$SwimSpeed.C)

### MODELS

g.AD1<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+Len, random=~1+FlowSpeed | Tag, method="ML", data=D)
g.AD2<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+Len+SwimSpeed.C, random=~1+FlowSpeed | Tag, method="ML", data=D)
g.AD3<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+FF, random=~1+FlowSpeed | Tag, method="ML", data=D)
g.AD4<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+FF+SwimSpeed.C, random=~1+FlowSpeed | Tag, method="ML", data=D)

g.SS1<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+Len, random=~1+FlowSpeed | Tag, method="ML", data=D, control = list(opt = "optim"))
g.SS2<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+Len+AngDif.C, random=~1+FlowSpeed | Tag, method="ML", data=D, control = list(opt = "optim"))
g.SS3<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+FF, random=~1+FlowSpeed | Tag, method="ML", data=D, control = list(opt = "optim"))
g.SS4<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+FF+AngDif.C.norm, random=~1+FlowSpeed | Tag, method="ML", data=D, control = list(opt = "optim"))



AIC(g.AD1,g.AD2,g.AD3,g.AD4)
#g.AD2 is best model

AIC(g.SS1,g.SS2,g.SS3,g.SS4)
#g.SS4 is best model

summary(g.AD2)$tTable
summary(g.SS4)$tTable

### FOLLOWING SHOWS PREDICTIONS FROM MODEL COEFFICIENT
### FOR INTERCEPT AND FLOW SPEED

g.AD2.coef<-coef(g.AD2)
g.SS4.coef<-coef(g.SS4)
FlowSpeed<-c(0,0.72)
L<-nrow(g.AD2.coef)
Pred.AngDif<-Pred.FlowSpeed<-array(dim=c(L,length(FlowSpeed)))
for (i in 1:L)
{
	Pred.AngDif[i,]<-g.AD2.coef[i,1]+g.AD2.coef[i,2]*FlowSpeed
	Pred.FlowSpeed[i,]<-g.SS4.coef[i,1]+g.SS4.coef[i,2]*FlowSpeed
}

XLAB<-expression(paste(   "Flow speed (m s", {}^-1,")"     ) )

par(las=0)
par(mfrow=c(2,2))
barplot(sort(g.AD2.coef[,2]),ylab="Slope",main="Normalized angular difference",xlab="Tag (sorted by slope")
barplot(sort(g.SS4.coef[,2]),ylab="Slope",main="Normalized swimming speed",xlab="Tag (sorted by slope")
matplot(FlowSpeed,t(Pred.AngDif),type="l",xlab=XLAB,ylab="Normalized angular difference")
matplot(FlowSpeed,t(Pred.FlowSpeed),type="l",xlab=XLAB,ylab="Normalized swimming speed")




############

#g.AD1.autcor<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+Len, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D)
#g.AD2.autcor<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+Len+SwimSpeed.C, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D)
#g.AD3.autcor<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+FF, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D)
#g.AD4.autcor<-lme(AngDif.C.norm~FlowSpeed+u+v+ke+FF+SwimSpeed.C, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D)

#g.SS1.autcor<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+Len, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D, control = list(opt = "optim"))
#g.SS2.autcor<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+Len+AngDif.C, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D, control = list(opt = "optim"))
#g.SS3.autcor<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+FF, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D, control = list(opt = "optim"))
#g.SS4.autcor<-lme(SwimSpeed.C.norm~FlowSpeed+u+v+ke+FF+AngDif.C.norm, random=~1+FlowSpeed | Tag, method="ML", correlation = corAR1(), data=D, control = list(opt = "optim"))


