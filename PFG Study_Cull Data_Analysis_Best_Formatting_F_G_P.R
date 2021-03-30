
#Some column titles may require re-typing. 
# Clear the environment
rm(list=ls())

# Load the packages. Install package 'geometry' if script is not working.
library(mixexp)
library(plyr)
library(sp)
library(arm)
library(mgcv)
library(geometry)

# Function for generatig colors for surfaces
source("0. Header_Functions.R")

# Read in the means and SDs on each diet
data<-read.csv("/Users/jibran/Desktop/PFG Study_Mastersheet_Cull Data.csv")
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
#I changed it from 3 to 5. Value of 8 looked too busy. This number is based on personal preference.
nlev<-5

		
########### PART 2 INTAKE GAMs ################# Part 1 was for RMTs and deleted from this script to make it simpler.

# The traits we are interested in 
# Note that the script does not work for Plasma Creatinine (gives an error message) probably because all the mice have the value of 18 as their result
traits<-c("Body.Weight.Cull","Body.Length","Blood.Glucose.Cull","Liver.Weight","Heart.Weight","Kidney.Weight","Spleen.Weight","Pancreas.Weight","Gonadal.Weight","Subcut.Weight","BAT.Weight","Gastroc.Weight","Quad.Weight","Colon.Length","Ratio.Visc.Subcut","Liver.Weight.Normalised","Heart.Weight.Normalised","Kidney.Weight.Normalised","Spleen.Weight.Normalised","Pancreas.Weight.Normalised","Gonadal.Weight.Normalised","Subcut.Weight.Normalised","BAT.Weight.Normalised","Gastroc.Weight.Normalised","Quad.Weight.Normalised","Liver.Triglyceride","Plasma.FGF21","Plasma.Sodium","Plasma.Chloride","Plasma.Urea","Plasma.Bilirubin","Plasma.Albumin","Plasma.Protein","Plasma.GGT","Plasma.ALT","Plasma.AST","Plasma.Uric.Acid","Plasma.Cholesterol","Plasma.Triglyceride")

titles<-c("Body Weight_Cull","Body Length", "BGL Cull", "Liver Weight", "Heart Weight", "Kidney Weight","Spleen Weight", "Pancreas Weight", "Gonadal Fat Weight", "Subcut Fat Weight", "BAT Weight","Gastroc Weight", "Quad Weight", "Colon Length","Ratio Visc Subcut","Liver Weight Normalised", "Heart Weight Normalised", "Kidney Weight Normalised", "Spleen Weight Normalised", "Pancreas Weight Normalised", "Gonadal Weight Normalised", "Subcut Weight Normalised", "BAT Weight Normalised", "Gastroc Weight Normalised", "Quad Weight Normalised", "Liver TG", "Plasma FGF21","Plasma Sodium","Plasma Chloride","Plasma Urea","Plasma Bilirubin","Plasma Albumin","Plasma Protein","Plasma GGT","Plasma ALT","Plasma AST","Plasma Uric Acid","Plasma Cholesterol","Plasma TG")

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


# Fitted list to hold some results for later
x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=fit.resolution)
y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=fit.resolution)
z.new<-z.vals
Fitted.List<-as.data.frame(expand.grid(x.new, y.new, z.new))
names(Fitted.List)<-nutrient.order
in.poly<-as.numeric(inhull(Fitted.List[,c(1:3)], data[,names(Fitted.List)]) != -1)

# Open the pdf file for plotting
pdf("GAMs.pdf", height=5.5, width=length(z.new) * 5)

# Set the layout
par(mfrow=c(1,length(z.new)), mar=c(5,5,5,1))	

# csv file to write results to
csv.file<-"GAM_results.csv"

# Open the file for results
write.table("", file=csv.file, sep=",", row.names=F, col.names=F)

for(k in 1:length(traits)){

	# Find the right outcome
	data$this.outcome<-data[,traits[k]]
	
	# Fit the model	
	# Start with the lowest and increase it till you find the k-value that works for most of the results. 
	# Only increase k-value further for those outcomes that give significant p value in the test for k-value
	# The script works for all k-values between 1 and 10 but it gives this message: "There were 12 warnings (use warnings() to see them)". This message disappears at k value of 11. Check with Alistair.
	
	GAM<-gam(this.outcome ~ s(eaten.P, eaten.Fru, eaten.Glu, k = 11), data=data, method="REML")
	
	# Add the test for k-value here, the test results should be non-significant if k value is suitable for analysis. 
	# These values are not the numbers for statistical significance of my results, these numbers do not need to go in the paper.
	# Only use the statistics saved to excel sheet for the paper.
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
	
	# Edit out based on the marker list
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
