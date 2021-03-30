
# A function created to find the outer perimeter over which the surface should be fitted for proportional data
findConvex.prop<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	x.new<-seq(0,1,len=res)
	y.new<-seq(0,1,len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}



inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) { 
	
	# https://tolstoy.newcastle.edu.au/R/e8/help/09/12/8784.html
   calpts <- as.matrix(calpts) 
   testpts <- as.matrix(testpts) 
   p <- dim(calpts)[2] 
	cx <- dim(testpts)[1] # rows in testpts
	nt <- dim(hull)[1] # number of simplexes in hull 
   nrmls <- matrix(NA, nt, p)

   degenflag <- matrix(TRUE, nt, 1) 
   
	for (i in 1:nt) { 
    		
		nullsp<-t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
  
		if (dim(nullsp)[1] == 1){
			nrmls[i,]<-nullsp
       		degenflag[i]<-FALSE
		}
	}
	
	 
   
	if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
	nrmls <- nrmls[!degenflag,] 
   	nt <- dim(nrmls)[1] 

   center = apply(calpts, 2, mean) 
   a <- calpts[hull[!degenflag,1],] 


   	nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)

	dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
	nrmls <- nrmls*matrix(dp, nt, p)

   aN <- diag(a %*% t(nrmls)) 
   val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 

   val[abs(val) < tol] <- 0 
   as.integer(sign(val)) 
}

