require("RcppArmadillo")
require("Rcpp")

stsreg_C<-function(x,y,xout=FALSE,outfun=out,iter=10,sc=pbvar,varfun=pbvar,
corfun=pbcor,plotit=FALSE,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		temp1<-.Call("stsregp1_C", X=x,Y=y)
		coef<-temp1$coef
		res<-temp1$res
	}
	if(ncol(x)>1){
		temp1<-.Call("stsreg_for", X=x,Y=y, IT=as.integer(iter))
		coef<-c(temp1$alpha, temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	e.pow<-varfun(yhat)/varfun(y)
	if(!is.na(e.pow)){
		if(e.pow>=1)
			e.pow<-corfun(yhat,y)$cor^2
		e.pow=as.numeric(e.pow)
		stre=sqrt(e.pow)
	}
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

tstsreg_C<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,plotit=FALSE,...){
	#
	# Compute a modified Theil-Sen regression estimator.
	# Use s-type initial estimate, eliminate points with
	# outlying residuals, then do regular Theil-Sen
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	res=stsreg_C(x,y)$res
	chk<-abs(res-median(res))/mad(res)
	xx<-x[chk<=2,]
	yy<-y[chk<=2]
	temp<-tsreg(xx,yy)
	list(coef=temp$coef,residuals=temp$res)
}

tshdreg_C<- function(x,y,HD=TRUE,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,tol=.0001,...){
	#
	#  Compute Theil-Sen regression estimator
	#
	#  Use back-fitting
	#  when there is more than one predictor
	#  and estimate intercept using Harrel-Davis estimator
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		coef<-.Call("tshd_C", X=x, Y=y, hd=as.integer(HD))
		res<-y-coef[2]*x-coef[1]
	}
	if(ncol(x)>1){
		temp1<-.Call("tshdreg_for", X=x,Y=y, IT=as.integer(iter), TOL=tol)
		coef<-c(temp1$alpha, temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	temp=varfun(y)
	if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
	e.pow=NULL
	if(temp>0){
		e.pow<-varfun(yhat)/varfun(y)
		if(!is.na(e.pow)){
			if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
			e.pow=as.numeric(e.pow)
			stre=sqrt(e.pow)
		}
	}
	res=NULL
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}


tsreg_C<-function(x,y,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,WARN=TRUE,HD=FALSE,...){
	#
	#  Compute Theil-Sen regression estimator
	#
	#  Use Gauss-Seidel algorithm
	#  when there is more than one predictor
	#
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		temp1<-.Call("tsp1reg_C", X=x, Y=y, HD=as.integer(HD))
		coef<-temp1$coef
		res<-temp1$res
	} 
	if(ncol(x)>1){
		temp1<-.Call("tsreg_for", X=x, Y=y, IT=as.integer(iter), HD=as.integer(HD))
		coef<-c(temp1$alpha,temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	temp=varfun(y)
	if(temp==0){
		if(WARN)print("Warning: When computing strength of association, measure of variation=0")
	}
	e.pow=NULL
	if(temp>0){
		e.pow<-varfun(yhat)/varfun(y)
		if(!is.na(e.pow)){
			if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
			e.pow=as.numeric(e.pow)
			stre=sqrt(e.pow)
		}
	}
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

outmgv_C<-function(x,y=NULL,plotit=TRUE,outfun=outbox,se=TRUE,op=1,
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,STAND=FALSE,...){
	#
	# Check for outliers using mgv method
	#
	# NOTE: if columns of the input matrix are reordered, this can
	# have an effect on the results due to rounding error when calling
	# the R function eigen.
	#
	#  (Argument STAND is included simply to avoid programming issues when outmgv is called by other 
	#  functions.)
	#
	require("RcppArmadillo")
	require("Rcpp")

	if(is.null(y[1]))m<-x
	if(!is.null(y[1]))m<-cbind(x,y)
	m=elimna(m)
	m=as.matrix(m)
	nv=nrow(m)
	temp<-mgvar_C(m,se=se,op=op,cov.fun=cov.fun,SEED=SEED)
	#if(fast)temp<-mgvdep.for(m,se=se)$distance
	temp[is.na(temp)]<-0
	if(ncol(m)==1){
		temp2=outpro(m)
		nout=temp2$n.out
		keep=temp2$keep
		temp2=temp2$out.id
	}
	if(ncol(m)>1){
		if(ncol(m)==2)temp2<-outfun(temp,...)$out.id
		if(ncol(m)>2)temp2<-outbox(temp,mbox=TRUE,gval=sqrt(qchisq(.975,ncol(m))))$out.id
		vec<-c(1:nrow(m))
		flag<-rep(T,nrow(m))
		flag[temp2]<-F
		vec<-vec[flag]
		vals<-c(1:nrow(m))
		keep<-vals[flag]
		if(plotit && ncol(m)==2){
			x<-m[,1]
			y<-m[,2]
			plot(x,y,type="n",xlab=xlab,ylab=ylab)
			flag<-rep(T,length(y))
			flag[temp2]<-F
			points(x[flag],y[flag],pch="*")
			points(x[temp2],y[temp2],pch="o")
		} else if(plotit && ncol(m)==3){
			require("scatterplot3d")
			x1<-m[,1]
			x2<-m[,2]
			y<-m[,3]
			flag<-rep(T,length(y))
			flag[temp2]<-F
			scatter3d<-scatterplot3d(x1, x2, y, type="n")
 			scatter3d$points3d(x1[flag], x2[flag], y[flag], pch="*")
 			scatter3d$points3d(x1[temp2], x2[temp2], y[temp2], pch="+", col="red")
		}
		nout=0
		if(!is.na(temp2[1]))nout=length(temp2)
	}
	list(n=nv,n.out=nout,out.id=temp2,keep=keep)
}

mgvar_C<-function(m,se=FALSE,op=0,cov.fun=covmve,SEED=TRUE){
#
# Find the center of a scatterplot, add point that
# increases the generalized variance by smallest amount
# continue for all points
# return the generalized variance
#  values corresponding to each point.
# The central values and point(s) closest to it get NA
#
# op=0 find central points using pairwise differences
# op!=0 find central points using measure of location
# used by cov.fun
#
# choices for cov.fun include
# covmve
# covmcd
# tbs (Rocke's measures of location
# rmba (Olive's median ball algorithm)
#
require("RcppArmadillo")
require("Rcpp")

	if(op==0)temp<-apgdis(m,se=se)$distance
	if(op!=0)temp<-out(m,cov.fun=cov.fun,plotit=FALSE,SEED=SEED)$dis
	flag<-(temp!=min(temp))
	temp2<-temp
	temp2[!flag]<-max(temp)
	flag2<-(temp2!=min(temp2))
	flag[!flag2]<-F
	if(sum(flag)>0)
		varvec<-.Call("mgvar_while", X=as.numeric(flag), M=m)
	else varvec<-NA
}


fdepthv2_C<-function(m,pts=NA,plotit=TRUE){
#
# Determine depth of points in pts relative to
# points in m
#
# Draw a line between each pair of distinct points
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# This function is slower than fdepth and requires
# space for a nc by nc matrix, nc=(n^2-n)/2.
# But it allows
# data to have a singular covariance matrix
# and it provides a more accurate approximation of
# halfspace depth.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  When plotting,
#  center is marked with a cross, +.
#
require("RcppArmadillo")
require("Rcpp")
	m<-elimna(m) # Remove missing values
	if(!is.na(pts[1]))remm<-m
	if(!is.matrix(m))dep<-unidepth(m)
	if(is.matrix(m)){
		nm<-nrow(m)
		nt<-nm
		nm1<-nm+1
	    if(!is.na(pts[1])){
    		if(ncol(m)!=ncol(pts))
    			stop("Number of columns of m is not equal to number of columns for pts")
			nt<-nm+nrow(pts)
			}
		}
	    if(ncol(m)==1)depth<-unidepth(m)
    	if(ncol(m)>1){
			m<-elimna(m) # Remove missing values
			nc<-(nrow(m)^2-nrow(m))/2
		  #  if(is.na(pts[1]))mdep <- matrix(0,nrow=nc,ncol=nrow(m))
		   # if(!is.na(pts[1])){
			#	mdep <- matrix(0,nrow=nc,ncol=nrow(pts))
			#}
		#ic<-0
		if(is.na(pts[1])) pts=matrix(, 2,2)
		mdep<-t(.Call("fdepthv2_for", M=m, PTS=pts))
		dep<-apply(mdep,2,min)
	}
    	if(ncol(m)==2 &&is.na(pts[1])){
			flag<-chull(m)
			dep[flag]<-min(dep)
		}
	    if(ncol(m)==2){
    		if(is.na(pts[1]) && plotit){
				plot(m, pch="+", cex=0.7)
				x<-m
				temp<-dep
				flag<-(temp>=median(temp))
				xx<-x[flag,]
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				lines(xx[temp,], col="red")
				lines(xx[c(temp[1],temp[length(temp)]),], , col="red")
			}
		}
		dep
}


outpro_C<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=FALSE,tr=.2,q=.5,pr=TRUE,...){
#
# Detect outliers using a modification of the
# Stahel-Donoho  projection method.
#
# Determine center of data cloud, for each point,
# connect it with center, project points onto this line
# and use distances between projected points to detect
# outliers. A boxplot method is used on the
# projected distances.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data.
#
# op=T
# means the .5 depth contour is plotted
# based on data with outliers removed.
#
# op=F
# means .5 depth contour is plotted without removing outliers.
#
#  MM=F  Use interquatile range when checking for outliers
#  MM=T  uses MAD.
#
#  If value for center is not specified,
#  there are four options for computing the center of the
#  cloud of points when computing projections:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
#
#  args q and tr having are not used by this function. They are included to deal
#  with situations where smoothers have optional arguments for q and tr
#
#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
	m<-as.matrix(m)
	if(pr){
		if(!STAND){
			#if(ncol(m)>1)cat("STAND=FALSE. If measures are on different scales,", 
			#					"might want to use STAND=TRUE\n")
		}
	}
	library(MASS)
	m=elimna(m)
	m<-as.matrix(m)
	nv=nrow(m)
	if(ncol(m)==1){
		dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
		dis<-sqrt(dis)
		dis[is.na(dis)]=0
		crit<-sqrt(qchisq(.975,1))
		chk<-ifelse(dis>crit,1,0)
		vec<-c(1:nrow(m))
		outid<-vec[chk==1]
		keep<-vec[chk==0]
	}
	if(ncol(m)>1){
		if(STAND)m=standm(m,est=median,scat=mad)
		if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
		if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
		if(cop==1 && is.na(center[1])){
		if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
		if(ncol(m)==2){
			tempd<-NA
			for(i in 1:nrow(m))
			tempd[i]<-depth(m[i,1],m[i,2],m)
			mdep<-max(tempd)
			flag<-(tempd==mdep)
			if(sum(flag)==1)center<-m[flag,]
			if(sum(flag)>1)center<-apply(m[flag,],2,mean)
		}
	}
	if(cop==2 && is.na(center[1])){
		center<-cov.mcd(m)$center
	}
	if(cop==4 && is.na(center[1])){
		center<-cov.mve(m)$center
	}
	if(cop==3 && is.na(center[1])){
		center<-apply(m,2,median)
	}
	if(cop==5 && is.na(center[1])){
		center<-tbs(m)$center
	}
	if(cop==6 && is.na(center[1])){
		center<-rmba(m)$center
	}
	if(cop==7 && is.na(center[1])){
		center<-spat(m)
	}
	outid.flag<-.Call("outpro_for", 
				 M=m, 
				 GVAL=gval, 
				 CENTER=center, 
				 MM=MM
				 )
	idv<-1:nrow(m)
	outid<-idv[outid.flag]
	keep<-idv[!outid.flag]
	if(ncol(m)==2){
		if(plotit){
			plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
			points(m[keep,1],m[keep,2],pch="*")
			if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
				if(op){
					tempd<-NA
					keep<-keep[!is.na(keep)]
					mm<-m[keep,]
					for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
					mdep<-max(tempd)
					flag<-(tempd==mdep)
					if(sum(flag)==1)center<-mm[flag,]
					if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
					m<-mm
				}
				points(center[1],center[2],pch="+")
				x<-m
				temp<-fdepth(m,plotit=FALSE)
				flag<-(temp>=median(temp))
				xx<-x[flag,]
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				lines(xx[temp,])
				lines(xx[c(temp[1],temp[length(temp)]),])
			}
		}
	}
	list(n=nv,n.out=length(outid),out.id=outid,keep=keep)
}



skip_C<-function(m,cop=6,MM=FALSE,op=3,mgv.op=0,outpro.cop=3,...){
#
# m is an n by p matrix
#
# Compute skipped location and covariance matrix
#
# op=1:
# Eliminate outliers using a projection method
# That is, first determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
#
# op=2 use mgv (function outmgv) method to eliminate outliers
#
# Eliminate any outliers and compute means
#  using remaining data.
# mgv.op=0, mgv uses all pairwise distances to determine center of the data
# mgv.op=1 uses MVE
# mgv.op=2 uses MCD
#
temp<-NA
m<-elimna(m)
if(op==1)temp<-outpro(m,plotit=FALSE,MM=MM,cop=outpro.cop,...)$keep
if(op==2)temp<-outmgv(m,plotit=FALSE,op=mgv.op,...)$keep
if(op==3)temp<-outpro_C(m,plotit=FALSE,MM=MM,cop=outpro.cop,...)$keep
if(op==4)temp<-outmgv_C(m,plotit=FALSE,op=mgv.op,...)$keep
val<-var(m[temp,])
loc<-apply(m[temp,],2,mean)
list(center=loc,cov=val)
}

addedOutliers <- function( outliers = NULL, p = NULL, n = 20,  nout = 1, rho = 0.9, cutoff = .01,
          					  niter = 500, seed = TRUE, FUN = outmgv_test, loc_scale = c(0, 2), 
          					  scale = TRUE, center = TRUE, ... ){
     require( parallel ) 
     if( is.null(p) && is.null(outliers) )
         stop( "'outliers' and 'p' cannot both be NULL. Please specify one of them." )
     if(seed)
     	set.seed(1)
     simdata <- split( data.frame( rmul( n = ( n - nout )*niter, p = p, rho = rho, ...)), f = rep(1:niter, each = n - nout))
     if( is.null(outliers) ) {
         if(seed) 
         	set.seed(2)
         outliers <- rmul( nout*niter, p, mean = loc_scale[1], sd = loc_scale[2] )
         outliers <- split( as.data.frame( outliers ) , f = rep(1:niter, each = nout) ) 
         temp  <- mapply(function( x, y ){ 
         					x <- rbind(as.matrix(x), as.matrix(y))
                            FUN( x, plotit = FALSE, cutoff = cutoff, ...)[2:3]}, x = outliers, y = simdata)
         #temp2  <- mapply(function( x ){ FUN( x , plotit=FALSE,cutoff=cutoff, ...)[2:3] }, x = simdata )
     } else {
         p <- length(outliers)
         outliers <- matrix( rep(outliers, nout), nout, byrow = TRUE )
         temp  <- mapply(function( x, y ){ 
         					FUN( rbind(outliers, as.matrix(x)), plotit = FALSE, cutoff = cutoff, ...)[2:3] }, x = simdata)
     }
     temp.totalout <-  unlist( temp[1, ])
     #temp.outp <- sum(unlist( mapply( function( x ){ if( length(x) == 1 && is.na(x)  ) 0 
     #                                                  else if( sum( x <= nout ) == nout ) 1 
     #                                                  else 0 }, x = temp[2,]) ))/niter
     temp.outp2 <- unlist( mapply( function( x ){ if( length(x) == 1 && is.na(x)  ) 0 
                                                       else sum( x <= nout )}, x = temp[2,]) )
     list(`Avg. total number of identified observations:` = mean(temp.totalout),
          `Avg. total number of of non-outliers identified:` = mean(temp.totalout - temp.outp2),
          `Avg. total number of of outliers identified:` = mean(temp.outp2),
          `Avg. rate of added outliers identified:` = mean(temp.outp2)/nout)
 }



