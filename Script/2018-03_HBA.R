library(ggplot2)
library(mgcv)
dat_all <- read.table('../Data/NewSmoltData for LME.2018-02-07.csv', sep=",", header=TRUE)
dat_all$Tag <- as.factor(dat_all$Tag)
dat_all$DateTime <- as.POSIXct(dat_all$DateTime, tz="UTC")
str(dat_all)

dat_all$diff <- dat_all$FlowDir - dat_all$GroundDir.B
idx <- which(dat_all$diff < -180)
dat_all$diff[idx] <- dat_all$diff[idx] + 360

dat <- dat_all[c('Tag', 'DateTime', 'X','Y','Z','PosErrX','PosErrY', 'FlowDir', 'SwimDir.A', 'SwimDir.B', 'SwimDir.C', 'GroundDir.B','AngDif.B','diff','ke','u','v','SwimSpeed.B','FlowSpeed')]
dat <- dat[dat$Y >= 6462600 & dat$Y <= 6462660, ]
# dat <- dat[dat$Y >= 6462575 & dat$Y <= 6462675, ]
# dat <- dat[!dat$Tag %in% c(12, 24, 38, 40, 48, 51, 73, 77, 80, 89),]
dat <- dat[(dat$PosErrX + dat$PosErrY) <= 3, ]
# dat <- na.omit(dat)
dat$speed_diff <- dat$FlowSpeed - dat$SwimSpeed.B

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
plot(speed_diff ~ (ke), data=dat)

plot(abs(diff) ~ log(ke), data=dat)



m <- gam(abs(diff) ~ FlowSpeed + te(ke, X) , data=dat)
vis.gam(m, theta=45, view=c("ke","X"))
plot(m)
plot(ke~ X, data=dat)






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


