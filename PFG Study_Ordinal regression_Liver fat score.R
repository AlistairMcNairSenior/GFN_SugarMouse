#######################################################
#### Mouse VGAM - vector generalized additive model####
#######################################################

####----
#setwd
dat=read.csv("/Users/jibran/Desktop/PFG Study_Mastersheet_Histology_Scores.csv")

#data (dealing with missings)
#selecting variables to be used 
usedvars <- c("eaten.P", "eaten.Fru", "eaten.Glu", "Fat.score" )
#Subset used
subset <- dat [usedvars]
#Excluding missing data from subset
dat <- na.omit(subset)

dat
fix(dat)
attach(dat)
####----

#install.packages("VGAM")
require(VGAM)

#macronutrient eaten
vgam2<-vgam(ordered(Fat.score)~s(eaten.P)+s(eaten.Fru)+s(eaten.Glu),data=dat,family=propodds())


opr<-par(mfrow=c(1,3), mar=c(5,5,7,1))
len<-51        # not sure what this number means just yet

#adjust font sizes here
par(cex=0.8)
par(cex.axis=1.8,cex.lab=2.5, cex.main=1, lwd=3)

#Protein
mdf1<-expand.grid(eaten.P=seq(min(dat$eaten.P),max(dat$eaten.P),len=len),
  eaten.Fru=median(dat$eaten.Fru),eaten.Glu=median(dat$eaten.Glu))
pred.P<-predict(vgam2,newdata=mdf1,type="response")

#plot(range(mdf1$eaten.P),c(0,1),type="n",xlab=labs1[1],ylab="Fitted probabilities")
plot(range(mdf1$eaten.P),c(0,1),type="n",xlab="Protein eaten", ylab="Fitted probabilities")
#plot(range(mdf1$eaten.P),c(0,1),xlab="Protein eaten",ylab="Fitted probabilities")

grid(col=8)
for(i in 1:4){ lines(mdf1$eaten.P,pred.P[,i],col=i+1) }

legend(4.5, 1, legend=0:3,lty=1,col=2:5,cex=1.8,bg="white")

#graphics.off()

#Fructose
mdf2<-expand.grid(eaten.Fru=seq(min(dat$eaten.Fru),max(dat$eaten.Fru),len=len),eaten.P=median(dat$eaten.P),eaten.Glu=median(dat$eaten.Glu))
pred.Fru<-predict(vgam2,newdata=mdf2,type="response")
#plot(range(mdf2$eaten.Fru),c(0,1),type="n",xlab=labs2[2],ylab="Fitted probabilities")
plot(range(mdf2$eaten.Fru),c(0,1),type="n",xlab="Fructose eaten", ylab="Fitted probabilities")
grid(col=8)
for(i in 1:4){lines(mdf2$eaten.Fru,pred.Fru[,i],col=i+1) }



#graphics.off()

#Glucose
mdf3<-expand.grid(eaten.Glu=seq(min(dat$eaten.Glu),max(dat$eaten.Glu),len=len),eaten.P=median(dat$eaten.P),eaten.Fru=median(dat$eaten.Fru))
pred.Glu<-predict(vgam2,newdata=mdf3,type="response")
#plot(range(mdf3$eaten.Glu),c(0,1),type="n",xlab=labs2[3],ylab="Fitted probabilities")
plot(range(mdf3$eaten.Glu),c(0,1),type="n",xlab="Glucose eaten", ylab="Fitted probabilities")
grid(col=8)
for(i in 1:4){lines(mdf3$eaten.Glu,pred.Glu[,i],col=i+1)}

summary(mdf1)
summary(vgam2)
#legend("top",legend=0:3,lty=1,col=2:5, title="Liver Fat score",cex=1.0,horiz=TRUE,bg="white")

#graphics.off()
