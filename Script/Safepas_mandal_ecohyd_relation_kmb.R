### READ AND CLEAN DATA
library(RCurl)
library(bbmle)
library(glmmTMB)

x <- getURL("https://raw.githubusercontent.com/kimmagnusb/SAFEPASS-Mandal/master/Data/NewSmoltData%20for%20LME.filt.Nozero.2018-03-23.csv")
D <- read.csv(text=x)


D$Tag <- as.factor(D$Tag)

#Make AngDif.B that can have both negative and positive values
D$AngDif.B.alt<-D$FlowDir-D$GroundDir.B 

D$Sum.posErr<-D$PosErrX+D$PosErrY




tf<-colnames(D)%in%c("AngDif.B.alt","AngDif.B","Sum.posErr","GroundDir.B","GroundSpeed.B","FlowDir","SwimSpeed.B","u","v","w","Pos","ke","Len","FlowSpeed","FF","Tag")
D<-D[,tf]
D<-na.omit(D)
# Devide the data into whats going on in 2D-polygon, and whats going on in the caotic 3D-polygon
D2<-D[D$Pos==2,]
D3<-D[D$Pos==3,]

#

# Remove strange individuals (different individuals in different poly)
#12, 24, 38, 40, 48, 51, 73, 77, 80, 89 (suggestion from henrik)
#12, 40,51, 77, 80, 89 (suggestion from ana, for the 2d)


## For 2d-part
D2_cleaned<-D2[!(D2$Tag=="12"| D2$Tag=="40"| D2$Tag=="51"|D2$Tag=="77"| D2$Tag=="80"| D2$Tag=="89"),]

hist(D2_cleaned$AngDif.B.alt)
hist(D2_cleaned$SwimSpeed.B)

# Angel deviation from flow (binom):
#D2_cleaned$ang.dev.flow<-ifelse(D2_cleaned$AngDif.B>10,"Dev","No_dev")

#table(D2_cleaned$ang.dev.flow)## Hmmm, seems like fish always deviates from flow, strange..

# Smoothing flow trajectories ++ based on running mean ()
library(zoo)
library(data.table)

setDT(D2_cleaned)
D2_cleaned[, MA5.FlowSpeed := rollmeanr(FlowSpeed, 5, na.pad = TRUE), by = Tag]
D2_cleaned[, MA5.FlowDir := rollmeanr(FlowDir, 5, na.pad = TRUE), by = Tag]
D2_cleaned[, MA5.GroundDir.B := rollmeanr(GroundDir.B, 5, na.pad = TRUE), by = Tag]
D2_cleaned[, MA5.u := rollmeanr(u, 5, na.pad = TRUE), by = Tag]
D2_cleaned[, MA5.v := rollmeanr(v, 5, na.pad = TRUE), by = Tag]
D2_cleaned[, SUM5.ke := sum(v, 5, na.pad = TRUE), by = Tag]


# Recalculate angle dif (Richard)
x1=cos(D2_cleaned$MA5.FlowDir*pi/180)
y1=sin(D2_cleaned$MA5.FlowDir*pi/180)
x2=cos(D2_cleaned$MA5.GroundDir.B*pi/180)
y2=sin(D2_cleaned$MA5.GroundDir.B*pi/180)
D2_cleaned$MA5.AngDif.B<-acos((x1*x2 + y1*y2))*180/pi

# Recaluculate angle dif alt (simple approach)
D2_cleaned$MA5.AngDif.B.alt<-D2_cleaned$MA5.FlowDir-D2_cleaned$MA5.GroundDir.B 

D2_cleaned<-na.omit(D2_cleaned)

# make bivariate variable
D2_cleaned$AngDif.B.bin<-ifelse(D2_cleaned$MA5.AngDif.B.alt >= -20 & D2_cleaned$MA5.AngDif.B.alt <= 20, 1, 0)


#Add time-variable for auto-cor.
times=NULL
for(i in levels(factor(D2_cleaned$Tag))){
  temp <- 1:nrow(D2_cleaned[D2_cleaned$Tag==i,])
  times <- as.vector(c(times,temp))
}


D2_cleaned$times <- as.factor(times) 

#####################################
###Using glmmTMB (no running mean)
#Testing different covar-structures (see vignette in package) 

#mod.1TMB<-glmmTMB(AngDif.B.alt ~ toep(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.2TMB<-glmmTMB(AngDif.B.alt ~ ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 

#mod.3TMB<-glmmTMB(AngDif.B.alt ~ cs(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.4TMB<-glmmTMB(AngDif.B.alt ~ diag(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects

## Testing fixed effect structure (always including an int with u and ground dir)

#Nb: these are rather slow for some reason 

mod1.noav.TMB<-glmmTMB(AngDif.B.alt ~ u+ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod2.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod3.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed+ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod4.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod5.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod6.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + v + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod7.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + v*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 

## Big issue: Convergence problems with glmmTMB and not with NLME. Not sure if they are set up poperly though


mod.1TMB$sdr$pdHess # do not converge
summary(mod.1TMB) # Returns: Error in `colnames<-`(`*tmp*`, value = "Corr") : 
#length of 'dimnames' [2] not equal to array extent
#In addition: Warning message:
#  In base::rbind(..., deparse.level = deparse.level) :
#  number of columns of result is not a multiple of vector length (arg 1)
D2_cleaned$Tag<-factor(D2_cleaned$Tag)

# averagred models (dismissed for now)

#mod1.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod2.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed*MA5.u, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod3.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed+MA5.u, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod4.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed*MA5.u + ke, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod5.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed*MA5.u + ke + Len, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod6.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed*MA5.u+ MA5.u*ke + Len, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)

#mod7.lme<-lme(MA5.AngDif.B.alt ~ MA5.FlowSpeed*MA5.u+ MA5.u*ke + Len + SwimSpeed.B, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)



#mod7.lme<-lme(AngDif.B.alt ~ FlowSpeed*u+ u*ke + Len, random=~1|Tag, method="ML", correlation = corAR1(), data=D2_cleaned)




#AICtab(mod5.lme,mod6.lme,mod7.lme)
#plot(mod7.lme) #looks sort of ok
#MuMIn::r.squaredGLMM(mod7.lme) # Biological significans...well 
# R2m        R2c 
# 0.06224640 0.06225344 

#####################
#NLME:Non averaged 
mod1.noav.lme<-lme(AngDif.B.alt ~ u, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod4.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod5.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod6.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B + v + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod7.noav.lme<-lme(AngDif.B.alt ~ u*GroundDir.B + v +v:GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)


AICtab(mod1.noav.lme,mod2.noav.lme,mod3.noav.lme,mod4.noav.lme,mod5.noav.lme,mod6.noav.lme,mod7.noav.lme)
plot(mod6.noav.lme) # looks strange, but ok in theory
MuMIn::r.squaredGLMM(mod6.noav.lme)
#R2m       R2c 
#0.9975408 0.9975873 

ang.dif.lme <-lme(AngDif.B.alt ~ u*GroundDir.B + v + ke, random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)





#Predicting
pred_frame <- with(D2_cleaned, expand.grid(u=seq(min(u),max(u), by=0.05),v=seq(min(v),max(v), by=0.05),GroundDir.B= seq(min(GroundDir.B),max(GroundDir.B), by=10), ke= seq(0,0.0035, by=0.0001),Tag=levels(Tag)))

pred<-predict(ang.dif.lme,pred_frame,
              se.fit = F)

pred_frame$pred<-as.data.frame(pred)[,1]
#pred_frame$upr<-as.data.frame(pred)[,3]
# Plot some results
ggplot(pred_frame, aes(x=u, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=ke, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=v, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=GroundDir.B, y=pred)) +stat_smooth()+theme_bw()




####################
# Modeling grounddir
mod1.2.lme<-lme(GroundDir.B ~ FlowDir , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.2.lme<-lme(GroundDir.B ~ FlowDir + FlowSpeed , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.2.lme<-lme(GroundDir.B ~ FlowDir*FlowSpeed , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
AICtab(mod1.2.lme,mod2.2.lme,mod3.2.lme)
mod4.2.lme<-lme(GroundDir.B ~ FlowDir*FlowSpeed +ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
AIC(mod3.2.lme,mod4.2.lme)

plot(mod3.2.lme) # not perfect 
MuMIn::r.squaredGLMM(mod3.2.lme)
#R2m         R2c 
#0.003157171 0.020087431 

##############
# SwimSpeed 

mod1.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed*AngDif.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod4.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod5.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod6.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod7.noav.speed.lme<-lme(SwimSpeed.B ~ u*GroundDir.B + v +v:GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod8.noav.speed.lme<-lme(SwimSpeed.B ~ 1, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod9.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed*AngDif.B , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)


AICtab(mod1.noav.speed.lme,mod2.noav.speed.lme,mod3.noav.speed.lme,mod4.noav.speed.lme,mod5.noav.speed.lme,mod6.noav.speed.lme,mod7.noav.speed.lme,mod8.noav.speed.lme,mod9.noav.speed.lme)
plot(mod3.noav.speed.lme) #looks almost ok..heteroscedasticity. Probably due to a lack of information in the variables. Also, doubtful if everything is linear. Letting it slip for now

MuMIn::r.squaredGLMM(mod3.noav.speed.lme)
#  R2m       R2c 
#0.3327892 0.6328491 


swimspeed.lme<-lme(SwimSpeed.B ~ FlowSpeed*AngDif.B + ke, random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)

summary(swimspeed.lme)

#Predicting
pred_frame <- with(D2_cleaned, expand.grid(FlowSpeed=seq(min(FlowSpeed),max(FlowSpeed), by=0.05),AngDif.B= seq(min(AngDif.B),max(AngDif.B), by=5), ke= seq(0,0.0035, by=0.0001),Tag=levels(Tag)))

pred<-predict(swimspeed.lme,pred_frame,
              se.fit = F)


pred_frame$pred<-as.data.frame(pred)[,1]
#pred_frame$upr<-as.data.frame(pred)[,3]
#pred_frame$lwr<-as.data.frame(pred)[,2]

# Plot some results
ggplot(pred_frame, aes(x=FlowSpeed, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=ke, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=AngDif.B, y=pred)) +stat_smooth()+theme_bw()





###Modeling bivariate ang dif () 

#TMB
mod.1TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.2TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u + ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.3TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.4TMB<-glmmTMB(AngDif.B.bin ~ (MA5.FlowSpeed+MA5.u)*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.5TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed+MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.6TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u+ Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.7TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u+ MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

mod.8TMB<-glmmTMB(AngDif.B.bin ~ MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects


AICtab(mod.1TMB,mod.2TMB,mod.3TMB,mod.4TMB,mod.5TMB,mod.6TMB,mod.7TMB,mod.8TMB)

summary(mod.7TMB)

confint(mod.7TMB)


#LMER
#mod5.glmer.bin<-glmer(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u + ke + Len +(1|Tag),family="binomial", data=D2_cleaned)
# Do not converge

table(D2_cleaned$AngDif.B.bin)


#Predicting
pred_frame <- with(D2_cleaned, expand.grid(MA5.FlowSpeed=seq(0.15,0.55, by=0.01),MA5.u= seq(-0.5,0.4, by=0.1), ke= seq(0,0.0035, by=0.0001),Len=mean(Len),Tag=levels(Tag)))

pred<-predict(mod.7TMB,pred_frame,
              se.fit = F)

predi<-predict(mod.7TMB,D2_cleaned,
               se.fit = F)

pred_frame$pred<-as.data.frame(pred)[,1]
pred_frame$upr<-as.data.frame(pred)[,3]


#for pred.frame
pred_frame$fit<-round(pred[,1],2)
pred_frame$lwr<-round(pred[,2],2)
pred_frame$upr<-round(pred[,3],2)


#Plotting
ggplot(pred_frame, aes(x=MA5.u, y=pred)) +stat_smooth()+theme_bw()
#theme(legend.background = element_rect(fill="white",colour="white"))+
#facet_grid( .~Type_fugl)+
#ylab("Estimates of bycatch per trip (bycatch only)")+
#xlab("Min. fishing depth (m)")+
#theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold")) 


ggplot(pred_frame, aes(x=MA5.FlowSpeed, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=ke, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=SwimSpeed.B, y=pred)) +stat_smooth()+theme_bw()


