

# Clear the environment
rm(list=ls())

# Load the packages
library(plyr)
library(fields)
library(mgcv)
library(sp)


# Load the Drosophila data
data<-read.csv("/Users/jibran/Desktop/PFG Study_FoodWaterIntakeByCage_Wk14.csv")
head(data)
str(data)

# I will illustrate using this fly data - at different % P and C - think of as % sucrose and starch
# These data are not actually isocaloric 0 they are at different concentrations (conc), but here we'll treat 'conc' as your % proteins

# We need to make sure conc is set as a factor (i.e. a category)
# JW changed conc to Protein in diet.
data$in.P<-as.factor(data$in.P)

# The column X.P is the percent protein - I will use this as my outcome - you could use % Sucrose or Starch
# JW changed X.P to Fructose in diet.
# The outcome is lifespan, but you can use intake obviously
# JW changed lifespan to energy intake


# Plot lifespan as function of % P in the diet - colored by concentration - in your case it will be %Sucrose
# JW plotted energy intake as a function of fructose in diet.
plot(data$in.Fru, data$energy.in, col=data$in.P, ylim=c(10, 70), xlim=c(-0.5, 8), pch=1, cex=0.4, xlab="Fructose (kJ/g)", ylab="Energy Intake (kJ/day)")

# Legend if you want
legend(6.2, 25, as.character(sort(unique(data$in.P))), col=c(1:3), pch=1, title = "Protein (kJ/g)")

# Lets fit a models
# Use k=4 as it fits the best model for this type of data.
model<-gam(energy.in ~ s(in.Fru, by=in.P, k=4) + in.P, data=data)
summary(model)


# Family: gaussian 
# Link function: identity 

# Formula:
# lifespan ~ s(X.P, by = conc, k = 5) + conc

# Parametric coefficients:


# Approximate significance of smooth terms:


# The parametric coefficients are the differences in lifespan between 45g/l conc and each other group. 
# For JW, parametric coefficients are differences in energy intake at three levels of dietary protein.
# The smooth terms tell us the effect of protein in each concentration
# For JW, smooth terms show the effect of fructose on energy intake at each level of dietary protein content

# We can use the same test we used from the surfaces GAMs

x<-fitted(model, type="response")
y<-resid(model, type="pearson")
test<-gam(y ~s(x))
print(summary(test))

# Looks good!

# If we want to add a curve to the plot we need to use predictions from the model - one prediction for each conc

# Here is 45
# For JW, this is protein content in 10% protein diet.
new.data<-data.frame(in.P=as.factor(1.43), in.Fru=seq(0, 7, 1))
out<-predict(model, newdata=new.data, se.fit=T)

# Add a line
lines(new.data$in.Fru, out$fit, col=1, lwd=2)
# And the SE
lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=1)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=1)


# Lets do 360 - you could do all if you want
# For JW, this is 20% protein diet.
new.data<-data.frame(in.P=as.factor(2.86), in.Fru=seq(0, 6, 1))
out<-predict(model, newdata=new.data, se.fit=T)

# Add a line
lines(new.data$in.Fru, out$fit, col=2, lwd=2)
# And the SE
lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=2)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=2)

# Lets do 360 - you could do all if you want
# For JW, this is 30% protein diet.
new.data<-data.frame(in.P=as.factor(4.29), in.Fru=seq(0, 5, 1))
out<-predict(model, newdata=new.data, se.fit=T)

# Add a line
lines(new.data$in.Fru, out$fit, col=3, lwd=2)
# And the SE
lines(new.data$in.Fru, out$fit + out$se.fit, lty=3, col=3)
lines(new.data$in.Fru, out$fit - out$se.fit, lty=3, col=3)



