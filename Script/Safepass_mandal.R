### READ AND CLEAN DATA
library(RCurl)
library(bbmle)
library(glmmTMB)

x <- getURL("https://raw.githubusercontent.com/kimmagnusb/SAFEPASS-Mandal/master/Data/NewSmoltData%20for%20LME.filt.Nozero.2018-03-23.csv")
D <- read.csv(text=x)
Mandal.3D <- read.csv("/data/home/kim.barum/R_raw_data/NewSmoltData.3D.csv",header=T,sep=",", na.strings = "")

D$Tag <- as.factor(D$Tag)
Mandal.3D$Tag <- as.factor(Mandal.3D$Tag)
#Make AngDif.B that can have both negative and positive values
D$AngDif.B.alt<-D$FlowDir-D$GroundDir.B 

D$Sum.posErr<-D$PosErrX+D$PosErrY




tf<-colnames(D)%in%c("AngDif.B.alt","AngDif.B","Sum.posErr","SSX.B","SSY.B","GroundDir.B","GroundSpeed.B","FlowDir","SwimSpeed.B","u","v","w","Pos","ke","Len","FlowSpeed","FF","Tag")
D<-D[,tf]
D<-na.omit(D)
#3d-data only
tf<-colnames(Mandal.3D)%in%c("AngDif.B.alt","AngDif.B","Sum.posErr","SSX.B","SSY.B","SSZ.B","GroundDir.B","GroundSpeed.B","FlowDir","SwimSpeed.B","u","v","w","Pos","ke","Len","FlowSpeed","FF","Tag")
Mandal.3D<-Mandal.3D[,tf]

Mandal.3D$GroundDir.B<-as.numeric(as.character(Mandal.3D$GroundDir.B))
Mandal.3D$GroundSpeed.B<-as.numeric(as.character(Mandal.3D$GroundSpeed.B))
Mandal.3D$SwimSpeed.B<-as.numeric(as.character(Mandal.3D$SwimSpeed.B))
Mandal.3D$AngDif.B<-as.numeric(as.character(Mandal.3D$AngDif.B))

Mandal.3D<-Mandal.3D[complete.cases(Mandal.3D), ]

# Devide the data into whats going on in 2D-polygon, and whats going on in the caotic 3D-polygon
D2<-D[D$Pos==2,]
D3<-D[D$Pos==3,]

#

# Remove strange individuals (different individuals in different poly)
#12, 24, 38, 40, 48, 51, 73, 77, 80, 89 (suggestion from henrik)
#12, 40,51, 77, 80, 89 (suggestion from ana, for the 2d)


## For 2d-part
D2_cleaned<-D2[!(D2$Tag=="12"| D2$Tag=="40"| D2$Tag=="51"|D2$Tag=="77"| D2$Tag=="80"| D2$Tag=="89"),]

## For 3d-part
Mandal.3D.cleaned <- Mandal.3D[!(Mandal.3D$Tag=="4"| Mandal.3D$Tag=="24"| Mandal.3D$Tag=="48"|Mandal.3D$Tag=="89"),]


hist(D2_cleaned$AngDif.B.alt)
hist(D2_cleaned$SwimSpeed.B)


hist(Mandal.3D.cleaned$GroundDir.B)
hist(Mandal.3D.cleaned$SwimSpeed.B)
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
#       2d-data
###Using glmmTMB (no running mean)
#Testing different covar-structures (see vignette in package) 

#mod.1TMB<-glmmTMB(AngDif.B.alt ~ toep(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects
D2_cleaned$AngDif.B<-D2_cleaned$AngDif.B/180
D2_cleaned$AngDif.B<-D2_cleaned$AngDif.B*180
mod.2TMB<-glmmTMB(AngDif.B ~ ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.3TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + Len  + ke+ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 

mod.4TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 

mod.11TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),data=D2_cleaned) 
mod.12TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u*SwimSpeed.B + v*SwimSpeed.B + u*v,data=D2_cleaned) 
mod.13TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u+v+ke,data=D2_cleaned) 
mod.14TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u*v,data=D2_cleaned) 
mod.15TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u*SwimSpeed.B + v*SwimSpeed.B + u*v+ke,data=D2_cleaned) 
mod.16TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~ke,data=D2_cleaned) 
mod.17TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u,data=D2_cleaned) 
mod.18TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~v,data=D2_cleaned) 
mod.19TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B,data=D2_cleaned) 
mod.20TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B+ke+u+v,data=D2_cleaned) 
mod.21TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke*v+ke*u+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B+ke+u+v,data=D2_cleaned) 


mod.5TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke + Len+ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.6TMB<-glmmTMB(AngDif.B ~ SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.7TMB<-glmmTMB(AngDif.B ~ u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.8TMB<-glmmTMB(AngDif.B ~ u+v+ ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.9TMB<-glmmTMB(AngDif.B ~ SwimSpeed.B+ ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 
mod.10TMB<-glmmTMB(AngDif.B ~ ke+ ar1(times + 0 | Tag)+(1|Tag),family="beta",data=D2_cleaned) 

AICtab(mod.2TMB,mod.3TMB,mod.4TMB,mod.5TMB,mod.6TMB,mod.7TMB,mod.8TMB,mod.9TMB,mod.10TMB,mod.11TMB,mod.12TMB)

AICtab(mod.12TMB,mod.13TMB,mod.14TMB,mod.15TMB,mod.16TMB)
AICtab(mod.16TMB,mod.17TMB,mod.18TMB,mod.19TMB,mod.20TMB,mod.15TMB,mod.13TMB,mod.21TMB)


#Model validation

sim_nbz = simulate(mod.19TMB, nsim = 10)
sim_nbz = do.call(cbind, sim_nbz)
sim_res_nbz = createDHARMa(simulatedResponse = sim_nbz, 
                           observedResponse = D2_cleaned$AngDif.B,
                           fittedPredictedResponse = predict(mod.19TMB),
                           integerResponse = TRUE)
plotSimulatedResiduals(sim_res_nbz)

simo=simulate(mod.20TMB, seed=1)
Simdat=D2_cleaned
Simdat$AngDif.B=simo[[1]]
Simdat=transform(Simdat,  
                 AngDif.B = AngDif.B, 
                 type="simulated")
D2_cleaned$type = "observed"  
Dat=rbind(D2_cleaned, Simdat) 

ggplot(Dat,  aes(AngDif.B, colour=type))+geom_density()

pred4<-predict(mod.20TMB,
               se.fit = F,type="response")

summary(lm(D2_cleaned$AngDif.B~as.data.frame(pred2)[,1]))
plot(as.data.frame(pred3)[,1],D2_cleaned$AngDif.B)


MuMIn::r.squaredGLMM(mod.4TMB)

#mod.3TMB<-glmmTMB(AngDif.B.alt ~ cs(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.4TMB<-glmmTMB(AngDif.B.alt ~ diag(times + 0 | Tag)+(1|Tag),dispformula=~0,data=D2_cleaned) #dispformula=~0 forces variance into the random effects

## Testing fixed effect structure (always including an int with u and ground dir)

#Nb: these are rather slow for some reason 

mod1.noav.TMB<-glmmTMB(AngDif.B.alt ~ u+ar1(times + 0 | Tag)+(1|Tag),,data=D2_cleaned) 
mod2.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod3.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed+ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod4.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod5.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + FlowSpeed + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod6.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + v + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 
mod7.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + v*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 

## Big issue: Convergence problems with glmmTMB and not with NLME. Not sure if they are set up poperly though
mod2.noav.TMB<-glmmTMB(AngDif.B.alt ~ u*GroundDir.B + ar1(times + 0 | Tag)+(1|Tag),data=D2_cleaned) 


mod.12TMB$sdr$pdHess # do not converge
summary(mod.12TMB) # Returns: Error in `colnames<-`(`*tmp*`, value = "Corr") : 
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



mod1.noav.TMB<-glmmTMB(AngDif.B ~ u+ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod2.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod3.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + FlowSpeed+ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod4.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod5.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + FlowSpeed + ke + ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod6.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + v + ke + ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 
mod7.noav.TMB<-glmmTMB(AngDif.B ~ u*GroundDir.B + v*GroundDir.B + ke + ar1(times + 0 | Tag)+(1|Tag),family=list(family="beta",link="logit"),data=D2_cleaned) 



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
pred_frame <- with(D2_cleaned, expand.grid(u=seq(min(u),max(u), by=0.05),v=seq(min(v),max(v), by=0.05),SwimSpeed.B= seq(min(SwimSpeed.B),max(SwimSpeed.B), by=0.05), ke= seq(0,0.0035, by=0.0001),Tag=levels(Tag),times))

pred<-predict(mod.15TMB,
              se.fit = F,type="response")

pred2<-predict(mod.13TMB,
               se.fit = F,type="response")

pred3<-predict(mod.19TMB,
               se.fit = F,type="response")

summary(lm(as.data.frame(pred2)[,1]~D2_cleaned$AngDif.B))
plot(as.data.frame(pred2)[,1]~D2_cleaned$AngDif.B)

summary(lm(as.data.frame(pred3)[,1]~D2_cleaned$AngDif.B))
plot(as.data.frame(pred3)[,1]~D2_cleaned$AngDif.B)


pred_frame$pred<-as.data.frame(pred)[,1]
#pred_frame$upr<-as.data.frame(pred)[,3]
# Plot some results
ggplot(pred_frame, aes(x=u, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=ke, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=v, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=GroundDir.B, y=pred)) +stat_smooth()+theme_bw()
#####

mod1.noav.lme<-lme(AngDif.B ~ u, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.noav.lme<-lme(AngDif.B ~ u*GroundDir.B + FlowSpeed, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.noav.lme<-lme(AngDif.B ~ u*GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod4.noav.lme<-lme(AngDif.B ~ u*GroundDir.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod5.noav.lme<-lme(AngDif.B ~ u*GroundDir.B + FlowSpeed + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod6.noav.lme<-lme(AngDif.B ~ u*GroundDir.B + v + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod7.noav.lme<-lme(AngDif.B ~ u*GroundDir.B + v +v:GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)

AICtab(mod1.noav.lme,mod2.noav.lme,mod3.noav.lme,mod4.noav.lme,mod5.noav.lme,mod6.noav.lme,mod7.noav.lme)
plot(mod7.noav.lme) # looks strange, but ok in theory
MuMIn::r.squaredGLMM(mod7.noav.lme)

####################
# Modeling grounddir
mod1.2.lme<-lme(GroundDir.B ~ FlowDir , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.2.lme<-lme(GroundDir.B ~ FlowDir + FlowSpeed , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.2.lme<-lme(GroundDir.B ~ FlowDir*FlowSpeed , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
AICtab(mod1.2.lme,mod2.2.lme,mod3.2.lme)
mod4.2.lme<-lme(GroundDir.B ~ FlowDir*FlowSpeed +ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod5.2.lme<-lme(GroundDir.B ~ FlowDir*FlowSpeed +u, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)



mod6.2.lme <- lme(GroundDir.B ~ FlowDir*FlowSpeed + ke * Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod7.2.lme <- lme(GroundDir.B ~ FlowDir*SwimSpeed.B + ke * Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)


AICtab(mod3.2.lme,mod4.2.lme,mod5.2.lme,mod6.2.lme,mod7.2.lme)
AICtab(mod3.2.lme,mod6.2.lme,mod7.2.lme)

plot(mod7.2.lme) # not perfect 
MuMIn::r.squaredGLMM(mod7.2.lme) #...not perfect at all
#R2m        R2c 
#0.02931053 0.04979340 

##############
# SwimSpeed 

#mod1.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
#mod2.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
#mod3.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed*AngDif.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod1.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod2.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod3.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod4.noav.speed.lme<-lme(SwimSpeed.B ~ u*GroundDir.B + v*GroundDir.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod5.noav.speed.lme<-lme(SwimSpeed.B ~ v*AngDif.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod6.noav.speed.lme<-lme(SwimSpeed.B ~ 1, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod7.noav.speed.lme<-lme(SwimSpeed.B ~ u*GroundDir.B + v*GroundDir.B + u*v + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod8.noav.speed.lme<-lme(SwimSpeed.B ~ u*GroundDir.B + v*GroundDir.B + u*v + w + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod9.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + Len  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
mod10.noav.speed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + Len  * ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
#mod9.noav.speed.lme<-lme(SwimSpeed.B ~ FlowSpeed*AngDif.B , random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)


AICtab(mod1.noav.speed.lme,mod2.noav.speed.lme,mod3.noav.speed.lme,mod4.noav.speed.lme,mod5.noav.speed.lme,mod6.noav.speed.lme,mod7.noav.speed.lme,mod8.noav.speed.lme,mod9.noav.speed.lme,mod10.noav.speed.lme)
plot(mod9.noav.speed.lme) #looks almost ok..heteroscedasticity. Probably due to a lack of information in the variables. Also, doubtful if everything is linear. Letting it slip for now
plot(mod2.noav.speed.lme)


MuMIn::r.squaredGLMM(mod9.noav.speed.lme)
#  R2m       R2c 
#0.3327892 0.6328491 


swimspeed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + Len  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed2.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + u*v  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed3.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + u*v  + ke + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed4.lme<-lme(SwimSpeed.B ~  AngDif.B + u*v  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed5.lme<-lme(SwimSpeed.B ~ u*v  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed6.lme<-lme(SwimSpeed.B ~ v*AngDif.B + u*v  + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
swimspeed7.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + u*v  + ke * Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)
AICtab(swimspeed.lme,swimspeed2.lme,swimspeed3.lme,swimspeed4.lme,swimspeed5.lme,swimspeed6.lme,swimspeed7.lme)

coef(swimspeed.lme)

summary(swimspeed.lme)$tTable

tidy(swimspeed.lme, effects = "fixed")

#Predicting
pred_frame <- with(D2_cleaned, expand.grid(u=seq(min(u),max(u), by=0.05),v=seq(min(v),max(v), by=0.05),w=seq(min(w),max(w), by=0.05),AngDif.B= seq(min(AngDif.B),max(AngDif.B), by=10), ke= seq(0,0.03, by=0.001),Len = mean(Len),Tag=levels(Tag)))

pred<-predict(swimspeed2.lme,pred_frame,
              se.fit = F)


pred_frame$pred<-as.data.frame(pred)[,1]
#pred_frame$upr<-as.data.frame(pred)[,3]
#pred_frame$lwr<-as.data.frame(pred)[,2]

# Plot some results
ggplot(pred_frame, aes(x=u, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=ke, y=pred)) +stat_smooth()+theme_bw()
ggplot(pred_frame, aes(x=v, y=pred)) +stat_smooth()+theme_bw()

contourplot(pred~u*AngDif.B, data=pred_frame, label.style="align", xlab="U", ylab="AngDif.B",cuts=50, as.table=T )
trellis.focus("panel", 1, 1)
lpoints(D2_cleaned$u,D2_cleaned$AngDif.B)
trellis.unfocus()

contourplot(pred~v*AngDif.B, data=pred_frame, label.style="align", xlab="V", ylab="AngDif.B",cuts=50, as.table=T )
trellis.focus("panel", 1, 1)
lpoints(D2_cleaned$u,D2_cleaned$AngDif.B)
trellis.unfocus()

contourplot(pred~v*u, data=pred_frame, label.style="align", xlab="V", ylab="U",cuts=50, ,as.table=T )
trellis.focus("panel", 1, 1)
lpoints(D2_cleaned$v,D2_cleaned$u)
trellis.unfocus()

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
ggplot(pred_frame, aes(x=Len, y=pred)) +stat_smooth()+theme_bw()

#####################################
#       3d-data
################
#make time variable for cor-function
times=NULL
for(i in levels(factor(Mandal.3D.cleaned$Tag))){
  temp <- 1:nrow(Mandal.3D.cleaned[Mandal.3D.cleaned$Tag==i,])
  times <- as.vector(c(times,temp))
}

Mandal.3D.cleaned$times_f <- as.factor(times)
Mandal.3D.cleaned$times <- times


## Swimmingspeed (horizontal)
# Using dregde

Global.mod.lme <- lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + w + ke*u*v + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

dredged.mod.lme <- dredge(Global.mod.lme)

subset(dredged.mod.lme, cumsum(weight) <= .95 )
# construcing "best" model based on dregde (lazy) --> aprox. 30% of AIC weights
swimspeed.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + w + ke*u + ke*v + u*v, random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

MuMIn::r.squaredGLMM(swimspeed.lme)

plot(swimspeed.lme)


summary(Mandal.3D.cleaned$u)
summary(Mandal.3D.cleaned$AngDif.B)
summary(Mandal.3D.cleaned$v)
summary(Mandal.3D.cleaned$w)
summary(Mandal.3D.cleaned$ke)

pred_frame <- with(Mandal.3D.cleaned, expand.grid(u=seq(min(u),max(u), by=0.05),v=seq(min(v),max(v), by=0.1),w=seq(min(w),max(w), by=0.05),AngDif.B= seq(min(AngDif.B),max(AngDif.B), by=10), ke= seq(0,0.03, by=0.005),Len = mean(Len),Tag=levels(Tag)))

pred<-predict(swimspeed.lme,pred_frame,
              full = T)

pred_frame$pred<-as.data.frame(pred)[,1]

library(lattice)

contourplot(pred~u*AngDif.B, data=pred_frame, label.style="align", xlab="U", ylab="AngDif.B", as.table=T )
trellis.focus("panel", 1, 1)
lpoints(Mandal.3D.cleaned$u,Mandal.3D.cleaned$AngDif.B)
trellis.unfocus()

contourplot(pred~v*AngDif.B, data=pred_frame, label.style="align", xlab="V", ylab="AngDif.B", as.table=T )
trellis.focus("panel", 1, 1)
lpoints(Mandal.3D.cleaned$v,Mandal.3D.cleaned$AngDif.B)
trellis.unfocus()

contourplot(pred~u*v, data=pred_frame, label.style="align", xlab="U", ylab="V", as.table=T )
trellis.focus("panel",1, 1)
lpoints(Mandal.3D.cleaned$u,Mandal.3D.cleaned$v)
trellis.unfocus()


contourplot(pred~ke*v, data=pred_frame, label.style="align", xlab="ke", ylab="V", as.table=T )
trellis.focus("panel",1, 1)
lpoints(Mandal.3D.cleaned$ke,Mandal.3D.cleaned$v)
trellis.unfocus()

contourplot(pred~ke*u, data=pred_frame, label.style="align", xlab="ke", ylab="U", as.table=T )
trellis.focus("panel",1, 1)
lpoints(Mandal.3D.cleaned$ke,Mandal.3D.cleaned$u)
trellis.unfocus()

ggplot(pred_frame, aes(x=w, y=pred)) +stat_smooth()+theme_bw()

# Model with Flowspeed and AngDif.B in stead of components of flow

Global.mod2.lme <- lme(SwimSpeed.B ~ AngDif.B*FlowSpeed + FlowSpeed*ke + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

dredged.mod2.lme <- dredge(Global.mod2.lme)

subset(dredged.mod2.lme, cumsum(weight) <= .95 )
# construcing "best" model based on dregde (lazy)
swimspeed.2.lme<-lme(SwimSpeed.B ~ AngDif.B*FlowSpeed + FlowSpeed*ke, random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

MuMIn::r.squaredGLMM(swimspeed.2.lme)

plot(swimspeed.2.lme)


pred_frame2 <- with(Mandal.3D.cleaned, expand.grid(FlowSpeed=seq(min(FlowSpeed),max(FlowSpeed), by=0.05),AngDif.B= seq(min(AngDif.B),max(AngDif.B), by=10), ke= seq(0,0.03, by=0.005),Tag=levels(Tag)))

pred2<-predict(swimspeed.2.lme,pred_frame2,
               full = T)

pred_frame2$pred<-as.data.frame(pred2)[,1]

contourplot(pred~FlowSpeed*AngDif.B, data=pred_frame2, label.style="align", xlab="FlowSpeed", ylab="AngDif.B", as.table=T )
trellis.focus("panel", 1, 1)
lpoints(Mandal.3D.cleaned$FlowSpeed,Mandal.3D.cleaned$AngDif.B)
trellis.unfocus()

contourplot(pred~FlowSpeed*ke, data=pred_frame2, label.style="align", xlab="FlowSpeed", ylab="ke", as.table=T )
trellis.focus("panel",1, 1)
lpoints(Mandal.3D.cleaned$FlowSpeed,Mandal.3D.cleaned$ke)
trellis.unfocus()


## Swimmingspeed (3D)
# Using dregde
Mandal.3D.cleaned$SSX.B<-as.numeric(levels(Mandal.3D.cleaned$SSX.B)[Mandal.3D.cleaned$SSX.B])
Mandal.3D.cleaned$SSY.B<-as.numeric(levels(Mandal.3D.cleaned$SSY.B)[Mandal.3D.cleaned$SSY.B])
Mandal.3D.cleaned$SSZ.B<-as.numeric(levels(Mandal.3D.cleaned$SSZ.B)[Mandal.3D.cleaned$SSZ.B])


#Mandal.3D.cleaned$SS.3D.A<-sqrt(Mandal.3D.cleaned$SSX.A^2+Mandal.3D.cleaned$SSY.A^2+Mandal.3D.cleaned$SSZ.A^2)
Mandal.3D.cleaned$SS.3D.B<-sqrt(Mandal.3D.cleaned$SSX.B^2+Mandal.3D.cleaned$SSY.B^2+Mandal.3D.cleaned$SSZ.B^2)
#Mandal.3D.cleaned$SS.3D.C<-sqrt(Mandal.3D.cleaned$SSX.C^2+Mandal.3D.cleaned$SSY.C^2+Mandal.3D.cleaned$SSZ.C^2)
Mandal.3D.cleaned<-Mandal.3D.cleaned[complete.cases(Mandal.3D.cleaned), ]


Global.mod.lme <- lme(SS.3D.B ~ u*AngDif.B + v*AngDif.B + w*u +w*ke +w*v + ke*u+ke*v+u*v + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dredged.mod.lme <- dredge(Global.mod.lme)

subset(dredged.mod.lme, cumsum(weight) <= .95 )

swimspeed.3D.lme<-lme(SS.3D.B ~ AngDif.B + ke + u + v + w + AngDif.B:u + AngDif.B:v + ke:v + ke:w + u:v + v:w , random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

MuMIn::r.squaredGLMM(swimspeed.3D.lme)
MuMIn::r.squaredGLMM(swimspeed.lme)


plot(swimspeed.3D.lme)
summary(swimspeed.3D.lme)


pred_frame3d <- with(Mandal.3D.cleaned, expand.grid(AngDif.B=seq(min(AngDif.B),max(AngDif.B), by=20),u= seq(min(u),max(u), by=0.05),v=seq(min(v),max(v), by=0.05),w=seq(min(w),max(w), by=0.05) ,ke= seq(0,0.03, by=0.005),Tag=levels(Tag)))

pred3d<-predict(swimspeed.3D.lme,pred_frame3d,
                full = T)

pred_frame3d$pred3d<-as.data.frame(pred3d)[,1]

pred_frame3d<-pred_frame3d[complete.cases(pred_frame3d), ]

pred_frame3d[pred_frame3d$pred3d<0,]$pred3d <- 0



contourplot(pred3d~u*ke, data=pred_frame3d, label.style="align", xlab="u", ylab="ke",cuts=30 ,as.table=T )

contourplot(pred3d~v*ke, data=pred_frame3d, label.style="align", xlab="v", ylab="ke",cuts=30 ,as.table=T )

contourplot(pred3d~w*ke, data=pred_frame3d, label.style="align", xlab="w", ylab="ke",cuts=30 ,as.table=T )

# Produce tables

stargazer(,
          #note you have to specify type
          type = "html",
          #note that the argument is "out" not "file"
          out="swimspeed.doc")

stargazer(swimspeed.lme, swimspeed.3D.lme,
          type="html",
          title="Regression Results",
          intercept.bottom = F,
          intercept.top = T,
          ci = F, digits=2,
          notes = "This is a caption.",
          model.names = T,
          single.row = T,
          out="swim speed_2.doc")
#covariate.labels =
#  c
#("Constant","Veriscolor",
# "Virginica", "Petal Width",
# "Versicolor x Petal Width",
)#"Virginica x Petal Width"))

stargazer(mod.20TMB,
          type="html",
          title="Regression Results",
          intercept.bottom = F,
          intercept.top = T,
          ci = F, digits=2,
          notes = "This is a caption.",
          model.names = T,
          single.row = T,
          out="angdif.doc")


stargazer(swimspeed.lme, swimspeed.3D.lme,type="html" ,title="Regression Results",
          dep.var.labels=c("Overall Rating","High Rating"),
          order=c("learning", "privileges"),
          keep.stat="n", ci=TRUE, ci.level=0.90, single.row=TRUE,out="swim speed_2.doc")








sjp.lme = function(fit, type = 'fe', printPlot = FALSE,  showIntercept = F) {
  if(type == "fe" & printPlot == FALSE) {
    confintervals = data.frame(nlme::intervals(fit)$fixed)
    confintervals$name = rownames(confintervals)
    mydf = confintervals[confintervals$name != "(Intercept)", ]
    names(mydf) = c("lower.CI", "OR", "upper.CI", "name")
    pvals = data.frame(nlme:::summary.lme(fit)$tTable)
    pvals$name = rownames(pvals)
    mydf = merge(mydf, pvals[, c("name","p.value")], by = 'name')
    mydf$x = mydf$sorting = mydf$fade = mydf$grp = NA
    ov = mydf$OR
    pv = mydf$p.value
    ps <- rep("", length(ov))
    ps <- sprintf("%.*f", 2, ov)
    for (i in 1:length(pv)) {
      ps[i] <- sjmisc::trim(paste(ps[i], sjPlot:::get_p_stars(pv[i])))
    }
    mydf$p = ps
    rownames(mydf) = mydf$name
    mydf = mydf[, c("OR", "lower.CI", "upper.CI", "p", "grp", "sorting", "x", "fade")]
    invisible(structure(class = "sjplme", list(plot = NULL, plot.list = NULL, 
                                               mydf = mydf)))
  } else {
    stop("Only type=fe and printPlot=F possible now.")
  }
}


sjp.lme(swimspeed.lme,file="swimspeed_linear.doc")

plot_model(swimspeed.lme, type="pred")

## Swimmingspeed (3D)
# Using dregde
Mandal.3D.cleaned$SSX.B<-as.numeric(levels(Mandal.3D.cleaned$SSX.B)[Mandal.3D.cleaned$SSX.B])
Mandal.3D.cleaned$SSY.B<-as.numeric(levels(Mandal.3D.cleaned$SSY.B)[Mandal.3D.cleaned$SSY.B])
Mandal.3D.cleaned$SSZ.B<-as.numeric(levels(Mandal.3D.cleaned$SSZ.B)[Mandal.3D.cleaned$SSZ.B])


#Mandal.3D.cleaned$SS.3D.A<-sqrt(Mandal.3D.cleaned$SSX.A^2+Mandal.3D.cleaned$SSY.A^2+Mandal.3D.cleaned$SSZ.A^2)
Mandal.3D.cleaned$SS.3D.B<-sqrt(Mandal.3D.cleaned$SSX.B^2+Mandal.3D.cleaned$SSY.B^2+Mandal.3D.cleaned$SSZ.B^2)
#Mandal.3D.cleaned$SS.3D.C<-sqrt(Mandal.3D.cleaned$SSX.C^2+Mandal.3D.cleaned$SSY.C^2+Mandal.3D.cleaned$SSZ.C^2)
Mandal.3D.cleaned<-Mandal.3D.cleaned[complete.cases(Mandal.3D.cleaned), ]




## Direction
dir.mod.lme <- lme(GroundDir.B ~ FlowDir, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod2.lme <- lme(GroundDir.B ~ FlowDir*ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod3.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod4.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod5.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + ke, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod6.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + w, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod7.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + ke + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )
dir.mod8.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + ke * Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

AICtab(dir.mod.lme,dir.mod2.lme,dir.mod3.lme,dir.mod4.lme,dir.mod5.lme,dir.mod6.lme,dir.mod7.lme,dir.mod8.lme)

dir.mod.3D.lme <- lme(GroundDir.B ~ FlowDir*SS.3D.B + ke * Len, random=~1|Tag, method="REML", correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

#trying gam
dir.mod.3D.gam <- gamm(GroundDir.B ~ FlowDir*SS.3D.B + ke * Len, random=list(Tag=~1), correlation = corAR1(form=~1|Tag), data=Mandal.3D.cleaned )

AICtab(dir.mod.3D.lme,dir.mod.3D.gam)

MuMIn::r.squaredGLMM(dir.mod.3D.lme)

plot(dir.mod.3D.lme)

pred_frame3dg <- with(Mandal.3D.cleaned, expand.grid(FlowDir=seq(min(FlowDir),max(FlowDir), by=20),SS.3D.B= seq(min(SS.3D.B),max(SS.3D.B), by=0.05),Len=seq(min(Len),max(Len), by=10),ke= seq(0,0.03, by=0.005),Tag=levels(Tag)))

pred3dg<-predict(dir.mod.3D.lme,pred_frame3dg,
                 full = T)

pred_frame3dg$pred3dg<-as.data.frame(pred3dg)[,1]

contourplot(pred3dg~FlowDir*SS.3D.B, data=pred_frame3dg, label.style="align", xlab="FlowDir", ylab="SS.3D.B",cuts=30 ,as.table=T )

contourplot(pred3dg~v*ke, data=pred_frame3dg, label.style="align", xlab="v", ylab="ke",cuts=30 ,as.table=T )

contourplot(pred3dg~w*ke, data=pred_frame3dg, label.style="align", xlab="w", ylab="ke",cuts=30 ,as.table=T )

###############
Mandal.3D.cleaned$AngDif.B<-Mandal.3D.cleaned$AngDif.B/180
mod3d.1TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + u*v + w*v +w*u + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned) 
mod3d.2TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + u*v + w*v +w*u + ke+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned) 
mod3d.3TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + u+v + w + ke+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned) 
mod3d.4TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u+v + w + ke+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned) 
mod3d.5TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + u*v + w + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned) 
mod3d.6TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned)
mod3d.7TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned)
mod3d.8TMB<-glmmTMB(AngDif.B ~ SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned)
mod3d.9TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned)
mod3d.10TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned)
mod3d.11TMB<-glmmTMB(AngDif.B ~  v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+u+v+w,data=Mandal.3D.cleaned)
mod3d.12TMB<-glmmTMB(AngDif.B ~  v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+u+ke+v+w,data=Mandal.3D.cleaned)



mod3d.7.1TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B  + w+ u*v + ke,data=Mandal.3D.cleaned)
mod3d.7.2TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B,data=Mandal.3D.cleaned)
mod3d.7.3TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke,data=Mandal.3D.cleaned)
mod3d.7.4TMB<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+w,data=Mandal.3D.cleaned)



AICctab(mod3d.1TMB,mod3d.2TMB,mod3d.3TMB,mod3d.4TMB,mod3d.5TMB,mod3d.6TMB,mod3d.7TMB,mod3d.8TMB,mod3d.9TMB,mod3d.10TMB,mod3d.11TMB,mod3d.12TMB)
AICctab(mod3d.7TMB,mod3d.7.1TMB,mod3d.7.2TMB,mod3d.7.3TMB,mod3d.7.4TMB)


simo=simulate(mod3d.7TMB, seed=1)
Simdat=Mandal.3D.cleaned
Simdat$AngDif.B=simo[[1]]
Simdat=transform(Simdat,  
                 AngDif.B = AngDif.B, 
                 type="simulated")
Mandal.3D.cleaned$type = "observed"  
Dat=rbind(Mandal.3D.cleaned, Simdat) 

ggplot(Dat,  aes(AngDif.B, colour=type))+geom_density()

pred5<-predict(mod3d.7TMB,
               se.fit = F,type="response")

summary(lm(Mandal.3D.cleaned$AngDif.B~as.data.frame(pred5)[,1]))
plot(as.data.frame(pred5)[,1],Mandal.3D.cleaned$AngDif.B)


summary(mod3d.6TMB)





####
# testing best model
###
# slicing data set by randomly choosing 50% of the individuals

(ids <- sample(unique(Mandal.3D.cleaned$Tag), 20))

#
Mandal.3D.cleaned_50.1 <- Mandal.3D.cleaned[Mandal.3D.cleaned$Tag %in% ids, ]
Mandal.3D.cleaned_50.2<-Mandal.3D.cleaned[!Mandal.3D.cleaned$Tag %in% ids, ]

testmod<-glmmTMB(AngDif.B ~ u*SS.3D.B + v*SS.3D.B + w*SS.3D.B+ u*v + ke*v+ke*u+ar1(times_f + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SS.3D.B+ke+u+v+w,data=Mandal.3D.cleaned_50.1)
pred.test<-predict(testmod,newdata = Mandal.3D.cleaned_50.2 ,allow.new.levels =T,
                   se.fit = F,type="response")

summary(lm(Mandal.3D.cleaned_50.2$AngDif.B~as.data.frame(pred.test)[,1]))
plot(as.data.frame(pred.test)[,1],Mandal.3D.cleaned_50.2$AngDif.B)


#########################
#      brt-model        #
#########################

#' Genneral information about BRT-models (Boosted Regression Trees)
#' Brt-models are a type of regression model technique which is very flexible. It thus has some benefits when fitting ecological data, which is usually non-linear and so on.
#' Given care in the model-fittingg, brt can give predictive advantages over methods as e.g. glm or gam.
#' The following script uses the dismo and gmb packages to optimise a brt model for the given species.
#' Analytically, BRT regularization involves jointly optimizing the number of  trees, learning rate and tree complexity.
#' Here optimizing numbers of trees is done through the gbm.step function, and the script includes a function that tries to aid in the optimazation prosess for learning rate and tree complexity (the get.train.diganostic.func). 
#' The next step is to fit the actual model used for predictions (brt_mod)
#'  Its is recomended to read through "A working guide to boosted regression trees" (Elith, et al 2009), before atempting your first go.




###############################
# part 1: load and filter data
###############################
library(stringr)
library(dismo)
library(dplyr)
library(gbm)
library(doParallel)



Mandal.3D.cleaned<-as.data.frame(Mandal.3D.cleaned) # convert to data.frame - needed for gbm.step input


###############################
# part 2: Make the brt model
###############################


#It is encuraged to do this with paralell computing speeds the prosess up to some extent.
#Identify cores on current system
cores<-detectCores(all.tests = FALSE, logical = FALSE)
cores

#Create training function for gbm.step
get.train.diganostic.func=function(tree.com,learn){
  #set seed for reproducibility
  k1<-gbm.step(data=Mandal.3D.cleaned, 
               gbm.x = c("Tag", "u","v","w","ke","FlowDir","FlowSpeed","Len"), #Include variables at will here
               gbm.y = "GroundDir.B",
               family = "gaussian", 
               tree.complexity = tree.com,
               learning.rate = learn,
               bag.fraction =0.7 ,
               prev.stratify=TRUE,
               n.folds=10,
               n.trees=500,
               step.size=100,
               silent=TRUE,
               plot.main = FALSE,
               n.cores=cores)
  
  k.out=list(interaction.depth=k1$interaction.depth,
             shrinkage=k1$shrinkage,
             n.trees=k1$n.trees,
             AUC=k1$self.statistics$discrimination,
             cv.AUC=k1$cv.statistics$discrimination.mean,
             deviance=k1$self.statistics$mean.resid,
             cv.deviance=k1$cv.statistics$deviance.mean)  
  return(k.out)
}

#define complexity and learning rate
tree.complexity<-c(1:9)
learning.rate<-c(0.01,0.025,0.005,0.001)

#setup parallel backend to use n processors
cl<-makeCluster(cores)
registerDoParallel(cl)

#Run the actual function
foreach(i = tree.complexity) %do% {
  foreach(j = learning.rate) %do% {
    nam=paste0("gbm_tc",i,"lr",j)
    assign(nam,get.train.diganostic.func(tree.com=i,learn=j))
    
  }
}

#Stop parallel
stopCluster(cl)
registerDoSEQ()

#Find all item in workspace that contain "gbm_tc"
train.all<-ls(pattern="gbm_tc")

#cbind each list that contains "gbm_tc"
train.results<-list(do.call(cbind,mget(train.all)))

#Place in a data frame
train.results<- do.call(rbind, lapply(train.results, rbind))
train.results <- data.frame(matrix(unlist(train.results),ncol=7 , byrow=T))

#Change column names
colnames(train.results)<-c("TC","LR","n.trees", "AUC", "cv.AUC", "dev", "cv.dev")

#Round 4:7 down to 3 digits
train.results[,4:7]<-round(train.results[,4:7],digits=3)

#Sort by cv.dev, cv.AUC, AUC
train.results<-train.results[order(train.results$cv.dev,-train.results$cv.AUC, -train.results$AUC),]

train.results #Includes a dataframe with ordered (numbered) choice based on AUC cv.dev and cv.AUC, be aware that there are mutiple ways of judging the models...

# Use best parametrization from train.results 

brt_mod<-gbm.fixed(data=analyse.df, gbm.x = c( ), gbm.y = "introduced",family = "bernoulli",tree.complexity = 1, learning.rate = 0.01,bag.fraction = 1,n.trees=1100)
names(brt_mod$gbm.call)[1] <- "dataframe"

predictors<-gbm.simplify(brt_mod,n.folds = 10, n.drops = "auto", alpha = 1, prev.stratify = TRUE, 
                         eval.data = NULL, plot = TRUE)



