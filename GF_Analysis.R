

rm(list=ls())

# Load the packages. 
library(mixexp)
library(plyr)
library(sp)
library(arm)
library(mgcv)
library(geometry)

# Function for generatig colors for surfaces
source("0. Header_Functions.R")


data<-read.csv("FileName.csv")
head(data)



# Set the resolution of the surface
surface.resolution<-501

# How many values to round surface
round.surf<-3

# This specifies the color scheme for surface - it is actually a function that returns a function
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# How many different colours should we use on the plot
no.cols<-256

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# How many levels should there be on the surface 
nlev<-5

		
# The traits we are interested in 

traits<-c("Body.Weight")

titles<-c("Body Weight")

# Order to plot the nutrients
nutrient.order<-c("eaten.Fru", "eaten.Glu", "eaten.P")

# Values to predict over
x.limits<-c(0, 25)
y.limits<-c(0, 25)

# Decide which P values to use
z.vals<-round(quantile(data$eaten.P)[2:4])
fit.resolution=101

# Labels
x.label<-"Fructose eaten"
y.label<-"Glucose eaten"
z.label<-"Protein eaten"


x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=fit.resolution)
y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=fit.resolution)
z.new<-z.vals
Fitted.List<-as.data.frame(expand.grid(x.new, y.new, z.new))
names(Fitted.List)<-nutrient.order
in.poly<-as.numeric(inhull(Fitted.List[,c(1:3)], data[,names(Fitted.List)]) != -1)


pdf("GAMs.pdf", height=5.5, width=length(z.new) * 5)

# Set the layout
par(mfrow=c(1,length(z.new)), mar=c(5,5,5,1))	

# csv file to write results to
csv.file<-"GAM_results.csv"


write.table("", file=csv.file, sep=",", row.names=F, col.names=F)

for(k in 1:length(traits)){

	# Find the right outcome
	data$this.outcome<-data[,traits[k]]
	
	
	GAM<-gam(this.outcome ~ s(eaten.P, eaten.Fru, eaten.Glu, k = 11), data=data, method="REML")
	
	
	x<-fitted(GAM, type="response")
	y<-resid(GAM, type="pearson")
	test<-gam(y ~s(x))
	print(summary(test))
	
	# Save the model's results
	write.table(traits[k], file=csv.file, sep=",", append=T, row.names=F, col.names=F)
	write.table(t(c("", colnames(summary(GAM)$p.table))), file=csv.file, sep=",", append=T, col.names=F, row.names=F)
	write.table(summary(GAM)$p.table, file=csv.file, sep=",", append=T, col.names=F)
	write.table(t(c("", colnames(summary(GAM)$s.table))), file=csv.file, sep=",", append=T, col.names=F, row.names=F)
	write.table(summary(GAM)$s.table, file=csv.file, sep=",", append=T, col.names=F)
	write.table(summary(GAM)$dev.expl * 100, file=csv.file, sep=",", append=T, row.names="Dev. Expl.", col.names=F)
	write.table("", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	
	# Do the predictions
	predictions<-predict(GAM, newdata=Fitted.List, type="response")
	

	predictions[which(in.poly == 0)]<-NA

	# Find the min and max values
	mn<-min(predictions, na.rm=T)
	mx<-max(predictions, na.rm=T)
	
	# Do the 3 quantiles
	for(i in 1:length(z.new)){
					
		# Subset for the ith quantile
		ith_Quantile<-predictions[which(Fitted.List[, nutrient.order[3]] == z.new[i])]
								
		surf<-matrix(ith_Quantile, nrow=fit.resolution)
					
		locs<-round((range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols)
		
		#Axis labels
		image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab=x.label, ylab=y.label, cex.lab=2.0, axes=FALSE)
		
  	#Protein Eaten on GF plots written within "( )"
		mtext(paste0( "(", z.label, " = ", round(z.new[i], 2), ")"), line=-2, cex=1.3)
		
		#Font size of the numbers on x and y axis
		axis(1, cex.axis=1.5)
		axis(2, cex.axis=1.5)
	
		#GF contours
		contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn, mx), nlev), labcex=1.0)
	
		#Plot title
		if(i == 2){mtext(titles[k], line=2, cex=0.75)}
		
	


		

					
	}

}


dev.off()
