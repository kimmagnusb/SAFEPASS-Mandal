
############
# CREATED: 2018-02-04
############

rm(list=ls()) # CLEAR ALL DATA

library(geoR)
library(nlme)
library(lattice)
library(colorRamps)

DIR<-"M:\\My Documents\\NINA projects\\Safepass\\Data2018\\Processed\\"

# CORVIF FUNCTION IS FROM ZUUR ET AL.
#Library files for courses provided by: Highland Statistics Ltd.
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
#Copyright Highland Statistics LTD.
#VIF FUNCTION.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  #cat("Correlations of the variables\n\n")
  #tmp_cor <- cor(dataz,use="complete.obs")
  #print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS

############

Transform2Normal.2p<-function(Ori){
T.box = boxcoxfit(Ori, lambda2=TRUE)
lambda = T.box$lambda[1]
lambda2 = T.box$lambda[2]
if(lambda==0){T.norm = log(Ori + lambda2)}
if(lambda!=0){T.norm = ((Ori + lambda2) ^ lambda - 1) / lambda}
return(T.norm)
}

D <- read.csv(paste(DIR,"NewSmoltData for LME.2018-02-04.csv",sep=""))


D$PosErr<-apply(cbind(D$PosErrX,D$PosErrY),1,mean,na.rm=TRUE)
D$Wt<-1/D$PosErr


######## BOXPLOTS

D$FlowSpeed.c<-">0.5"
tf<-D$FlowSpeed<0.5; D$FlowSpeed.c[tf]<-"0.4-0.5"
tf<-D$FlowSpeed<0.4; D$FlowSpeed.c[tf]<-"0.3-0.4"
tf<-D$FlowSpeed<0.3; D$FlowSpeed.c[tf]<-"0.2-0.3"
tf<-D$FlowSpeed<0.2; D$FlowSpeed.c[tf]<-"0.1-0.2"
tf<-D$FlowSpeed<0.1; D$FlowSpeed.c[tf]<-"0.0-0.1"
D$FlowSpeed.c<-factor(D$FlowSpeed.c,levels=c("0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",">0.5"))

par(mfrow=c(1,3))
x<-tapply(D$AngDif,list(D$Tag,D$FlowSpeed.c),mean)
boxplot(x,ylab="Angular difference",xlab="Flow speed",main="All data")
tf<-D$Pos==2
x<-tapply(D$AngDif[tf],list(D$Tag[tf],D$FlowSpeed.c[tf]),mean)
boxplot(x,ylab="Angular difference",xlab="Flow speed",main="2D grid")
tf<-D$Pos==3
x<-tapply(D$AngDif[tf],list(D$Tag[tf],D$FlowSpeed.c[tf]),mean)
boxplot(x,ylab="Angular difference",xlab="Flow speed",main="3D grid")


### LME

tf<-!is.na(D$AngDif); D<-D[tf,]

D$AngDif.norm<-Transform2Normal.2p(D$AngDif)

tf<-colnames(D)%in%c("FlowSpeed","ke","Len","FF")
Pred<-D[,tf]
corvif(Pred)

# FIT ALL PREDICTORS
g1<-lme(AngDif.norm~FlowSpeed+ke+Len+FF, weights=~Wt, random=~1 | Tag, method="ML", data=D)
summary(g1)$tTable

# REMOVE FF BECAUSE IT IS NOT STRONGLY SIGNIFICANT
g2<-lme(AngDif.norm~FlowSpeed+ke+Len, weights=~Wt, random=~1 | Tag, method="ML", data=D)
summary(g2)$tTable

AIC(g1,g2)

# EXAMINE INTERACTION BETWEEN KE AND LEN
g3<-lme(AngDif.norm~FlowSpeed+ke*Len, weights=~Wt, random=~1 | Tag, method="ML", data=D)
summary(g3)$tTable

AIC(g1,g2,g3)


# PREDICT FROM MODEL
g<-lme(AngDif.norm~FlowSpeed+FlowSpeed:ke, weights=~Wt, random=~1 | Tag, method="ML", data=D)
summary(g)

ke<-seq(0,0.02,0.0001)
FlowSpeed<-seq(0,0.72,0.001)
new.dat <- expand.grid(FlowSpeed=FlowSpeed,ke=ke)
gg<-predict(g,newdata=new.dat,level=0)
new.dat<-data.frame(new.dat)
new.dat$gg<-gg

par(mfrow=c(1,1))
plot(new.dat$FlowSpeed,new.dat$ke,pch=15,cex=2,col=new.dat$colvec)
points(D$FlowSpeed,D$ke,pch=1,cex=0.5)

colvec<-blue2green2red(100); colvec=colvec[1:99]
levelplot(gg ~ FlowSpeed * ke, data = new.dat,cuts = 98, region = TRUE,col.regions=colvec)




###################
# FOLLOWING EXAMPLES INCLUDE AUTOCORRELATION
# BUT WILL ONLY WORK IF YOU LOTS OF RAM
###################

gg1<-lme(AngDif.norm~FlowSpeed+ke+Len+FF, weights=~Wt, random=~1 | Tag, correlation = corAR1(), method="ML", data=D)
gg2<-lme(AngDif.norm~FlowSpeed+ke+Len, weights=~Wt, random=~1 | Tag, correlation = corAR1(), method="ML", data=D)
gg3<-lme(AngDif.norm~FlowSpeed+ke*Len, weights=~Wt, random=~1 | Tag, correlation = corAR1(), method="ML", data=D)

AIC(gg1,gg2,gg3)

