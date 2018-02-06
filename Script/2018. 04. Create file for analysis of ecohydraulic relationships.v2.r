
############
# CREATED: 2018-02-06
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

D<-data.frame(Tag=DD$tag,DateTime=DD$time,X=DD$fish_x,Y=DD$fish_y,Z=DD$fish_z)
D$DateTime<-strptime(D$DateTime,format="%Y-%m-%d %H:%M:%S")

### CALCULATE GROUND SPEED

TagU<-unique(D$Tag)
TagUL<-length(TagU)

GSX.A.ls<-GSY.A.ls<-GSZ.A.ls<-list()
GSX.B.ls<-GSY.B.ls<-GSZ.B.ls<-list()
GSX.C.ls<-GSY.C.ls<-GSZ.C.ls<-list()

XD.ls<-YD.ls<-list()
	
for (i in 1:TagUL)
{
	tf<-D$Tag==TagU[i]
	Sub<-D[tf,]
	SubL<-nrow(Sub)

# METHOD A	
	D1<-Sub$DateTime[1:(SubL-1)]
	D2<-Sub$DateTime[2:(SubL)]
	TimeDif.A<-as.numeric(difftime(D2,D1,units="secs"))
	XDist.A<-Sub$X[2:(SubL)]-Sub$X[1:(SubL-1)]
	YDist.A<-Sub$Y[2:(SubL)]-Sub$Y[1:(SubL-1)]
	ZDist.A<-Sub$Z[2:(SubL)]-Sub$Z[1:(SubL-1)]
	XVel.A<-XDist.A/TimeDif.A
	YVel.A<-YDist.A/TimeDif.A
	ZVel.A<-ZDist.A/TimeDif.A
	GSX.A.ls[[i]]<-c(XVel.A,NA)
	GSY.A.ls[[i]]<-c(YVel.A,NA)
	GSZ.A.ls[[i]]<-c(ZVel.A,NA)
	
# METHOD B	
	D1<-Sub$DateTime[1:(SubL-2)]
	D2<-Sub$DateTime[3:(SubL)]
	TimeDif.B<-as.numeric(difftime(D2,D1,units="secs"))
	XDist.B<-Sub$X[3:(SubL)]-Sub$X[1:(SubL-2)]
	YDist.B<-Sub$Y[3:(SubL)]-Sub$Y[1:(SubL-2)]
	ZDist.B<-Sub$Z[3:(SubL)]-Sub$Z[1:(SubL-2)]
	XVel.B<-XDist.B/TimeDif.B
	YVel.B<-YDist.B/TimeDif.B
	ZVel.B<-ZDist.B/TimeDif.B
	GSX.B.ls[[i]]<-c(NA,XVel.B,NA)
	GSY.B.ls[[i]]<-c(NA,YVel.B,NA)
	GSZ.B.ls[[i]]<-c(NA,ZVel.B,NA)

# METHOD A	
	D1<-Sub$DateTime[1:(SubL-1)]
	D2<-Sub$DateTime[2:(SubL)]
	TimeDif.C<-as.numeric(difftime(D2,D1,units="secs"))
	XDist.C<-Sub$X[2:(SubL)]-Sub$X[1:(SubL-1)]
	YDist.C<-Sub$Y[2:(SubL)]-Sub$Y[1:(SubL-1)]
	ZDist.C<-Sub$Z[2:(SubL)]-Sub$Z[1:(SubL-1)]
	XVel.C<-XDist.A/TimeDif.C
	YVel.C<-YDist.A/TimeDif.C
	ZVel.C<-ZDist.A/TimeDif.C
	GSX.C.ls[[i]]<-c(NA,XVel.C)
	GSY.C.ls[[i]]<-c(NA,YVel.C)
	GSZ.C.ls[[i]]<-c(NA,ZVel.C)	
	
	XD.ls[[i]]<-c(NA,diff(Sub$X))
	YD.ls[[i]]<-c(NA,diff(Sub$Y))

	}

GSX.A<-unlist(GSX.A.ls)
GSX.B<-unlist(GSX.B.ls)
GSX.C<-unlist(GSX.C.ls)

GSY.A<-unlist(GSY.A.ls)
GSY.B<-unlist(GSY.B.ls)
GSY.C<-unlist(GSY.C.ls)

GSZ.A<-unlist(GSZ.A.ls)
GSZ.B<-unlist(GSZ.B.ls)
GSZ.C<-unlist(GSZ.C.ls)

XD<-unlist(XD.ls)
YD<-unlist(YD.ls)


D$PosErrX<-DD$pos_err_x
D$PosErrY<-DD$pos_err_y
D$PosEryZ<-DD$pos_err_z

D$GSX.A<-GSX.A
D$GSX.B<-GSX.B
D$GSX.C<-GSX.C

D$GSY.A<-GSY.A
D$GSY.B<-GSY.B
D$GSY.C<-GSY.C

D$GSZ.A<-GSZ.A
D$GSZ.B<-GSZ.B
D$GSZ.C<-GSZ.C

D$XD<-XD
D$YD<-YD
D$u<-DD$u
D$v<-DD$v
D$w<-DD$w
D$ke<-DD$ke
D$eps<-DD$eps
D$Pos<-DD$positioning
D$Len<-DD$len
D$Wei<-DD$wei
D$SI<-DD$SI
D$FF<-DD$FF

temp<-geopolar2(D$GSX.A,D$GSY.A)
D$GroundSpeed.A<-temp[,3]
D$GroundDir.A<-temp[,2]

temp<-geopolar2(D$GSX.B,D$GSY.B)
D$GroundSpeed.B<-temp[,3]
D$GroundDir.B<-temp[,2]

temp<-geopolar2(D$GSX.C,D$GSY.C)
D$GroundSpeed.C<-temp[,3]
D$GroundDir.C<-temp[,2]

D$FlowDir<-geopolar2(D$u,D$v)[,2]
D$FlowSpeed<-sqrt(D$u^2+D$v^2)

x1=cos(D$FlowDir*pi/180)
y1=sin(D$FlowDir*pi/180)
x2=cos(D$GroundDir.A*pi/180)
y2=sin(D$GroundDir.A*pi/180)
D$AngDif.A<-acos((x1*x2 + y1*y2))*180/pi
x2=cos(D$GroundDir.B*pi/180)
y2=sin(D$GroundDir.B*pi/180)
D$AngDif.B<-acos((x1*x2 + y1*y2))*180/pi
x2=cos(D$GroundDir.C*pi/180)
y2=sin(D$GroundDir.C*pi/180)
D$AngDif.C<-acos((x1*x2 + y1*y2))*180/pi



write.csv(file=paste(DIR,"NewSmoltData for LME.2018-02-06.csv",sep=""),D,row.names=FALSE)
