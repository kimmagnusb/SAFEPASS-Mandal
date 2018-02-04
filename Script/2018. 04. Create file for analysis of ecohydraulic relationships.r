
############
# CREATED: 2018-02-04
############

rm(list=ls()) # CLEAR ALL DATA

DIR<-"M:\\My Documents\\NINA projects\\Safepass\\Data2018\\Processed\\"



geopolar2=function(x,y)
{
L=length(x)
thetarad=thetadeg=rep(NA,L)
true=y>0 & x>=0 
true[is.na(true)]=FALSE
thetarad[true] = atan(x[true]/y[true])
true=y<0 & x>= 0
true[is.na(true)]=FALSE
thetarad[true] = atan(x[true]/y[true]) + pi
true=x<=0 & y<0 
true[is.na(true)]=FALSE
thetarad[true] = atan(x[true]/y[true]) + pi
true=x<=0 & y>0 
true[is.na(true)]=FALSE
thetarad[true] = atan(x[true]/y[true]) + 2*pi
true=x>0 & y==0
true[is.na(true)]=FALSE
thetarad[true]=0.5*pi
true=x==0 & y<0
true[is.na(true)]=FALSE
thetarad[true]=pi
true=x<0 & y==0
true[is.na(true)]=FALSE
thetarad[true]=1.5*pi
true=x==0 & y>0
true[is.na(true)]=FALSE
thetarad[true]=2*pi
thetadeg=thetarad*180/pi
#thetadeg=round(thetadeg,0)
dist=sqrt(x^2+y^2)
out=cbind(thetarad,thetadeg,dist)
return(out)
}

### READ DATA

DD <- read.csv(paste(DIR,"NewSmoltSSIIM.SBG.Env.2018-02-04.csv",sep=""))


### CREATE NEW ARRAY

D<-data.frame(Tag=DD$tag,DateTime=DD$time,X=DD$fish_x,Y=DD$fish_y,Z=DD$fish_z,DistX=DD$dist_x,DistY=DD$dist_y)
D$DateTime<-strptime(D$DateTime,format="%Y-%m-%d %H:%M:%S")

### CALCULATE GROUND SPEED

TagU<-unique(D$Tag)
TagUL<-length(TagU)

GSX.ls<-GSY.ls<-GSZ.ls<-list()
XD.ls<-YD.ls<-list()
	
for (i in 1:TagUL)
{
	tf<-D$Tag==TagU[i]
	Sub<-D[tf,]
	SubL<-nrow(Sub)

	D1<-Sub$DateTime[1:(SubL-2)]
	D2<-Sub$DateTime[3:(SubL)]
	TimeDif<-as.numeric(difftime(D2,D1,units="secs"))
	XDist<-Sub$X[3:(SubL)]-Sub$X[1:(SubL-2)]
	YDist<-Sub$Y[3:(SubL)]-Sub$Y[1:(SubL-2)]
	ZDist<-Sub$Z[3:(SubL)]-Sub$Z[1:(SubL-2)]
	XVel<-XDist/TimeDif
	YVel<-YDist/TimeDif
	ZVel<-ZDist/TimeDif

	GSX.ls[[i]]<-c(NA,XVel,NA)
	GSY.ls[[i]]<-c(NA,YVel,NA)
	GSZ.ls[[i]]<-c(NA,ZVel,NA)

	XD.ls[[i]]<-c(NA,diff(Sub$X))
	YD.ls[[i]]<-c(NA,diff(Sub$Y))

	}

GSX<-unlist(GSX.ls)
GSY<-unlist(GSY.ls)
GSZ<-unlist(GSZ.ls)
XD<-unlist(XD.ls)
YD<-unlist(YD.ls)


D$PosErrX<-DD$pos_err_x
D$PosErrY<-DD$pos_err_y
D$PosEryZ<-DD$pos_err_z

D$GSX<-GSX
D$GSY<-GSY
D$GSZ<-GSZ
D$XD<-XD
D$YD<-YD
D$u<-DD$u
D$v<-DD$v
D$w<-DD$w
D$ke<-DD$ke
D$Pos<-DD$positioning
D$Len<-DD$len
D$Wei<-DD$wei
D$SI<-DD$SI
D$FF<-DD$FF

temp<-geopolar2(D$DistX,D$DistY)
D$GroundSpeed<-temp[,3]
D$GroundDir<-temp[,2]
D$FlowDir<-geopolar2(D$u,D$v)[,2]

x1=cos(D$FlowDir*pi/180)
y1=sin(D$FlowDir*pi/180)
x2=cos(D$GroundDir*pi/180)
y2=sin(D$GroundDir*pi/180)
D$AngDif<-acos((x1*x2 + y1*y2))*180/pi
D$FlowSpeed<-sqrt(D$u^2+D$v^2)

write.csv(file=paste(DIR,"NewSmoltData for LME.2018-02-04.csv",sep=""),D,row.names=FALSE)
