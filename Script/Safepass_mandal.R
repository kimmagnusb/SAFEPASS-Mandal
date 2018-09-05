### READ AND CLEAN DATA
library(RCurl)
library(bbmle)
library(glmmTMB)

#2-d data
x <- getURL("https://raw.githubusercontent.com/kimmagnusb/SAFEPASS-Mandal/master/Data/NewSmoltData%20for%20LME.filt.Nozero.2018-03-23.csv")
D <- read.csv(text=x)

D$Tag <- as.factor(D$Tag)
#Make AngDif.B that can have both negative and positive values
D$AngDif.B.alt<-D$FlowDir-D$GroundDir.B 

D$Sum.posErr<-D$PosErrX+D$PosErrY




tf<-colnames(D)%in%c("AngDif.B.alt","AngDif.B","Sum.posErr","SSX.B","SSY.B","GroundDir.B","GroundSpeed.B","FlowDir","SwimSpeed.B","u","v","w","Pos","ke","Len","FlowSpeed","FF","Tag")
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


hist(Mandal.3D.cleaned$GroundDir.B)
hist(Mandal.3D.cleaned$SwimSpeed.B)
# Angel deviation from flow (binom):

## The following section can be omitted
#####################################
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


#################################

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

mod.15TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~u*SwimSpeed.B + v*SwimSpeed.B + u*v+ke,data=D2_cleaned) 

mod.19TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B,data=D2_cleaned) 
mod.20TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B+ke+u+v,data=D2_cleaned) 
mod.21TMB<-glmmTMB(AngDif.B ~ u*SwimSpeed.B + v*SwimSpeed.B + u*v  + ke*v+ke*u+ar1(times + 0 | Tag)+(1|Tag),family=beta_family(),dispformula=~SwimSpeed.B+ke+u+v,data=D2_cleaned) 






#Predicting
#Include orkladata here
pred<-predict(mod.15TMB,newdata= #Orkladata
              se.fit = F,type="response")

#For alternative model setup
pred2<-predict(mod.20TMB,
               se.fit = F,type="response")

pred3<-predict(mod.21TMB,
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

##############
# SwimSpeed 


swimspeed3.lme<-lme(SwimSpeed.B ~ u*AngDif.B + v*AngDif.B + u*v  + ke + Len, random=~1|Tag, method="ML", correlation = corAR1(form=~1|Tag), data=D2_cleaned)

pred<-predict(swimspeed3.lme,#orkladata,
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
#mod.1TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.2TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u + ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.3TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.4TMB<-glmmTMB(AngDif.B.bin ~ (MA5.FlowSpeed+MA5.u)*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.5TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed+MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.6TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u+ Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.7TMB<-glmmTMB(AngDif.B.bin ~ MA5.FlowSpeed*MA5.u+ MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects

#mod.8TMB<-glmmTMB(AngDif.B.bin ~ MA5.u*ke + Len +(1|Tag),family="binomial",data=D2_cleaned) #dispformula=~0 forces variance into the random effects




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



