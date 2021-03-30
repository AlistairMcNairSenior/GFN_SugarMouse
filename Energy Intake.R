


rm(list=ls())

# Load the packages
library(plyr)
library(fields)
library(mgcv)
library(sp)


# Load the data
data<-read.csv("FileName.csv")
head(data)
str(data)




data$in.P<-as.factor(data$in.P)




plot(data$in.Fru, data$energy.in, col=data$in.P, ylim=c(10, 70), xlim=c(-0.5, 8), pch=1, cex=0.4, xlab="Fructose (kJ/g)", ylab="Energy Intake (kJ/day)")


legend(6.2, 25, as.character(sort(unique(data$in.P))), col=c(1:3), pch=1, title = "Protein (kJ/g)")


model<-gam(energy.in ~ s(in.Fru, by=in.P, k=4) + in.P, data=data)
summary(model)




x<-fitted(model, type="response")
y<-resid(model, type="pearson")
test<-gam(y ~s(x))
print(summary(test))


new.data<-data.frame(in.P=as.factor(1.43), in.Fru=seq(0, 7, 1))
out<-predict(model, newdata=new.data, se.fit=T)


lines(new.data$in.Fru, out$fit, col=1, lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=1)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=1)



new.data<-data.frame(in.P=as.factor(2.86), in.Fru=seq(0, 6, 1))
out<-predict(model, newdata=new.data, se.fit=T)


lines(new.data$in.Fru, out$fit, col=2, lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=2)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=2)


new.data<-data.frame(in.P=as.factor(4.29), in.Fru=seq(0, 5, 1))
out<-predict(model, newdata=new.data, se.fit=T)


lines(new.data$in.Fru, out$fit, col=3, lwd=2)

lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=3)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=3)



