library(ggplot2)
library(mgcv)
library(nlme)
library(caTools)
library(plyr)
dat_all <- read.table('../Data/NewSmoltData for LME.2018-02-07.csv', sep=",", header=TRUE)
dat_all$Tag <- as.factor(dat_all$Tag)
dat_all$DateTime <- as.POSIXct(dat_all$DateTime, tz="UTC")
str(dat_all)

dat_all$diff <- dat_all$FlowDir - dat_all$GroundDir.B
idx <- which(dat_all$diff < -180)
dat_all$diff[idx] <- dat_all$diff[idx] + 360


dat <- dat_all[c('Tag', 'DateTime', 'X','Y','Z','PosErrX','PosErrY', 'FlowDir', 'SwimDir.A', 'SwimDir.B', 'SwimDir.C', 'GroundDir.B','AngDif.B','diff','ke','u','v','SwimSpeed.B','FlowSpeed', 'GroundSpeed.B')]
# dat <- dat[dat$Y >= 6462600 & dat$Y <= 6462660, ]
dat <- dat[dat$Y >= 6462575 & dat$Y <= 6462675, ]
dat <- dat[!dat$Tag %in% c(9,15,53,57,59,63,77,80,97,98),]
# dat <- dat[!dat$Tag %in% c(12, 24, 38, 40, 48, 51, 73, 77, 80, 89),]
# dat <- dat[(dat$PosErrX + dat$PosErrY) <= 3, ]
# dat <- na.omit(dat)
dat$speed_diff <- dat$FlowSpeed - dat$SwimSpeed.B

plot(GroundSpeed.B ~ SwimSpeed.B + FlowSpeed, data=dat)

par(mfrow=c(3,3), pch=".")
plot(AngDif.B~ ke, data=dat)
plot(AngDif.B~ log(ke), data=dat)
plot(AngDif.B~ FlowSpeed, data=dat_all)
plot(AngDif.B~ FlowSpeed, data=dat)
plot(diff~ FlowSpeed, data=dat)
plot(AngDif.B~ u, data=dat)
plot(AngDif.B~ v, data=dat)
plot(AngDif.B~ X, data=dat)
plot(GroundDir.B ~ FlowDir, data=dat)

dat$angdif_runmean <- unlist(by(dat, dat$Tag, function(k) runmean(k$AngDif.B, k=1), simplify=TRUE))

tags <- unique(dat$Tag)
for(i in 1:length(tags)){
	tag <- tags[i]
	dat_i <- dat[dat$Tag == tag, ]
	plot(runmean(-1*diff(dat_i$X) + diff(dat_i$Y), k=9))
	
	plot(rollapply(dat_i, 5, function(k) weighted.mean(as.numeric(k[, 'AngDif.B']), as.numeric(k[,'GroundSpeed.B'])), by.column=FALSE))
	
	Xs <- apply(dat_i, 1, function(k) as.numeric(k['X']) + rnorm(25, 0, sd=as.numeric(k['PosErrX'])))
	Ys <- apply(dat_i, 1, function(k) as.numeric(k['Y']) + rnorm(25, 0, sd=as.numeric(k['PosErrY'])))
	
	diffs <- (-1*apply(Xs, 1, diff) + apply(Ys, 1, diff))
	matplot(t(apply(diffs, 1, quantile, probs=c(0.25,0.5,0.75))), type="l")
	
}





pdf('tracks_angDif.pdf')
for(i in 1:length(tags)){
	tag <- tags[i]
	p <- ggplot(data=dat[dat$Tag == tag, ], aes(x=X, y=Y), title=tag) + geom_point(aes(col=(angdif_runmean > 90))) + coord_fixed(xlim=c(413625, 413775), ylim=c(6462575, 6462675)) + scale_colour_discrete() + geom_path() + ggtitle(tag)
	print(p)
}
dev.off()

ggplot(data=dat, aes(x=X, y=Y, col=angdif_runmean )) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10)) + facet_wrap(~Tag)






plot(runmean(dat$AngDif.B, k=5))


m <- gam(GroundDir.B ~ s(FlowDir), data=dat)
plot(predict(m) ~ dat$FlowDir)


ggplot(data=dat, aes(x=X, y=AngDif.B, col=Tag)) + geom_point() + facet_wrap(~Tag)
ggplot(data=dat, aes(x=X, y=AngDif.B, col=Tag)) + geom_point() + facet_wrap(~Tag)
# ggplot(data=dat, aes(x=DateTime, y=diff, col=Tag)) + geom_point() + facet_wrap(~Tag)


par(mfrow=c(3,2), pch=".")
plot(SwimSpeed.B~ ke, data=dat)
plot(SwimSpeed.B~ log(ke), data=dat)
plot(SwimSpeed.B~ FlowSpeed, data=dat)
plot(SwimSpeed.B~ u, data=dat)
plot(SwimSpeed.B~ v, data=dat)



plot(u~ FlowSpeed, data=dat); abline(0,1)
plot(v~ FlowSpeed, data=dat); abline(0,1)
plot(u~ log(ke), data=dat)
plot(v~ log(ke), data=dat)
ggplot(data=dat, aes(x=u, y=v, colour=FlowSpeed)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))

m <- gam(FlowSpeed ~ ti(u,v), data=dat)
vis.gam(m, theta=195)

plot(dat$ke)

dat$ke_std <- (dat$ke - mean(dat$ke)) / (2*sd(dat$ke))
dat$u_std <- (dat$u - mean(dat$u)) / (2*sd(dat$u))
dat$v_std <- (dat$v - mean(dat$v)) / (2*sd(dat$v))

dat_sub <- dat[which(!is.na(dat$SwimSpeed.B) & !is.na(dat$ke_std) & !is.na(dat$u_std) & !is.na(dat$v_std)), ]
m <- lme(SwimSpeed.B~ ke_std + u_std*v_std, data=dat_sub, random=~1|Tag)
plot(predict(m) ~ dat_sub$SwimSpeed.B); abline(0,1)
plot(predict(m) ~ dat_sub$ke_std)
plot(predict(m) ~ dat_sub$u)
plot(predict(m) ~ dat_sub$v)
plot(predict(m) ~ dat_sub$FlowSpeed)


dat_sub <- dat[which(!is.na(dat$SwimSpeed.B) & !is.na(dat$AngDif.B) & !is.na(dat$ke_std) & !is.na(dat$u_std) & !is.na(dat$v_std)), ]
m <- lme(abs(diff)~ ke_std + u_std*v_std, data=dat_sub, random=~1|Tag)
# m <- lme(AngDif.B~  v_std, data=dat_sub, random=~1|Tag)
summary(m)
plot(predict(m) ~ dat_sub$AngDif.B); abline(0,1)
plot(predict(m) ~ dat_sub$ke_std)
plot(predict(m) ~ dat_sub$u)
plot(predict(m) ~ dat_sub$v)
plot(predict(m) ~ dat_sub$FlowSpeed)

plot(resid(m, type="normalized"))
plot(resid(m, type="normalized") ~ dat_sub$ke_std)
plot(resid(m, type="normalized") ~ dat_sub$u_std)
plot(resid(m, type="normalized") ~ dat_sub$v_std)
plot(resid(m, type="normalized") ~ dat_sub$FlowSpeed)
plot(resid(m, type="normalized") ~ dat_sub$Tag)
plot(resid(m, type="normalized") ~ dat_sub$Y)



head(dat)
str(dat)
plot((PosErrY + PosErrX) ~ Y, data=dat)

ggplot(data=dat, aes(x=X, y=Y, colour=log(PosErrX + PosErrY))) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10)) + geom_hline(aes(yintercept=6462600)) + geom_hline(aes(yintercept=6462660))
ggplot(data=dat, aes(x=X, y=Y, colour=abs(diff))) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10)) + geom_hline(aes(yintercept=6462600)) + geom_hline(aes(yintercept=6462660))
ggplot(data=dat, aes(x=X, y=Y, colour=FlowDir)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=SwimDir.B)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=GroundDir.B)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=log(ke))) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=u)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=v)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=FlowSpeed)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=SwimSpeed.B)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))
ggplot(data=dat, aes(x=X, y=Y, colour=speed_diff)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = heat.colors(10))

ggplot(data=dat, aes(x=X, y=diff, col=Tag)) + geom_line() 


plot(SwimSpeed.B ~ FlowSpeed, data=dat)
plot(FlowSpeed ~ u, data=dat)
plot(SwimSpeed.B ~ u, data=dat)
plot(SwimSpeed.B ~ v, data=dat)
plot(SwimSpeed.B ~ log(ke), data=dat)
plot(u ~ log(ke), data=dat)
plot(FlowSpeed ~ log(ke), data=dat)
plot(speed_diff ~ (ke), data=dat)

plot(abs(diff) ~ log(ke), data=dat)



m <- gam(abs(diff) ~ FlowSpeed + te(ke, X) , data=dat)
vis.gam(m, theta=45, view=c("ke","X"))
plot(m)
plot(ke~ X, data=dat)




plot(abs(diff) ~ FlowSpeed, data=dat)
plot(abs(diff) ~ log(ke), data=dat)
m <- lme(abs(diff) ~ FlowSpeed, random=~1|Tag, data=dat, na.action=na.omit)
summary(m)



plot(speed_diff ~ abs(diff), data=dat)



par(mfrow=c(2,1))
plot(SwimDir.B ~ FlowDir, data=dat); abline(a=0,b=1)
plot(GroundDir.B ~ FlowDir, data=dat); abline(a=0,b=1)

plot(SwimSpeed.B ~ FlowSpeed, data=dat)

plot(SwimDir.A ~ SwimDir.B, data=dat[dat$Tag == 10, ]); abline(a=0,b=1)
plot(dat[dat$Tag == 10, ]$SwimDir.B)

plot(GroundDir.B ~ FlowDir, data=dat); abline(0,1)

ggplot(data=dat, aes(x=X, y=ke, col=abs(diff))) + geom_point()+ scale_colour_gradientn(colours = topo.colors(10))
ggplot(data=dat, aes(x=X, col=ke, y=abs(diff))) + geom_point()+ scale_colour_gradientn(colours = topo.colors(10))

plot(abs(diff) ~ ke, data=dat)


m <- gam(abs(diff) ~ te(X, (ke)), data=dat)
vis.gam(m)


plot(dat$diff ~ dat$X, col=dat$Tag)
m <- gam(abs(diff) ~ s(X,Y), data=dat)
vis.gam(m)

ggplot(data=dat, aes(x=X, y=Y, colour=abs(diff))) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10)) + facet_wrap(~Tag)
ggplot(data=dat, aes(x=X, y=Y, colour=AngDif.B)) + geom_point() + coord_fixed() + scale_colour_gradientn(colours = topo.colors(10)) + facet_wrap(~Tag)



tag <- 10
dat_tag <- dat[dat$Tag == tag, ]

plot(SwimDir.B ~ GroundDir.B, data=dat)
plot(dat_tag$SwimDir.B)
plot(dat_tag$GroundDir.B)
plot(Y~X, data=dat_tag, asp=1)
plot(SwimDir.B ~ Y, data=dat_tag)


plot(abs(dat$FlowDir - dat$GroundDir.B) ~ dat$AngDif.B)

plot(AngDif.B ~ X, data=dat)

library(mgcv)
m <- gam(AngDif.B ~ te(X,Y), data=dat)
vis.gam(m)


head(dat)

head(dat_tag)
diff(dat_tag$Y)


gnu <- pi*1.5
sin(gnu) + cos(gnu)


