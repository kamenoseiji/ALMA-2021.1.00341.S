#-------- Functions to be used in GBT reduction
#---- baselinfit( dataframe, weight ) : A function to subtact baseline from a given spectrum. The baseline is calculated by 3rd order spline
baselinefit <- function( filename, weight, Skip=3, velRef=0.0){
	#  filename : Input text data, output from GBTIDL
	#  weight : weight to select baseline range
	tempframe <- read.table(filename, skip=Skip, header=F)
	names(tempframe) <- c("velocity", "rawflux")
	baseline <- predict(smooth.spline((1:nrow(tempframe)), tempframe$rawflux, all.knots=TRUE, df=8, w=weight), (1:nrow(tempframe)))$y
	dataframe <- data.frame(velocity = tempframe$velocity - velRef, flux = tempframe$rawflux - baseline, baseline = baseline)
	return(dataframe)
}

#---- fitdata( dataframe, veloc_range, rms_range) : A function to clip given data within a given velocity range, to prepare for spectral fitting. 
fitdata <- function( dataframe, veloc_range = c(950, 2050), rms_range = c(1000, 1200) ){
	#  datafram : Input data frame
	#  veloc_range : Velocity range to clip the input data
	#  rms_range : Velocity range to calculate standard deviation of the spectrum. Line-free channels should be used.
	tmpdata <- subset(dataframe, ((velocity > veloc_range[1]) & (velocity < veloc_range[2])))
	rms <- sd(subset(dataframe, ((velocity > rms_range[1]) & (velocity < rms_range[2])))$flux)
	return( data.frame(tmpdata, err = rep(rms, length(tmpdata$velocity))) )
}
	

#---- compplot( velocity, gausscomp ) : A function to plot Gaussian components
compplot <- function( velocity, gausscomp ){
	#  velocity : a vector of velocities to plot
	#  gausscomp : Gaussian components … num_comp x 3 (mean velocity, width, and peak intensity) x 4 (estimate, sd, t-value, p-value)
	num_comp <- length( gausscomp ) / 12   # 12 parameters per component
	for( comp_index in 1:num_comp ){
		gausspeak <- gausscomp[3* comp_index ]
		gaussmean <- gausscomp[3* comp_index - 2]
		gausswidth <- gausscomp[3*comp_index - 1]
		lines(velocity, gausspeak * exp(-0.5* ((velocity - gaussmean) / gausswidth)^2), col=comp_index, lty='dashed')
		cat(sprintf("Comp%d : VLSR=%6.1f km/s Sigma=%6.2f Peak=%6.2f mJy\n", comp_index, gaussmean, gausswidth, 1e3*gausspeak))
	}
}

compdata <- function(comp_index, mjd, gausscomp1, gausscomp2, gausscomp3, gausscomp4, gausscomp5, gausscomp6){
	velocity <- c(gausscomp1[3* comp_index - 2], gausscomp2[3* comp_index - 2], gausscomp3[3* comp_index - 2], gausscomp4[3* comp_index - 2], gausscomp5[3* comp_index - 2], gausscomp6[3* comp_index - 2])
	width    <- c(gausscomp1[3* comp_index - 1], gausscomp2[3* comp_index - 1], gausscomp3[3* comp_index - 1], gausscomp4[3* comp_index - 1], gausscomp5[3* comp_index - 1], gausscomp6[3* comp_index - 1])
	peak     <-  c(gausscomp1[3* comp_index], gausscomp2[3* comp_index], gausscomp3[3* comp_index], gausscomp4[3* comp_index], gausscomp5[3* comp_index], gausscomp6[3* comp_index])
	velocity_err <- c(gausscomp1[6* comp_index - 2], gausscomp2[6* comp_index - 2], gausscomp3[6* comp_index - 2], gausscomp4[6* comp_index - 2], gausscomp5[6* comp_index - 2], gausscomp6[6* comp_index - 2])
	width_err    <- c(gausscomp1[6* comp_index - 1], gausscomp2[6* comp_index - 1], gausscomp3[6* comp_index - 1], gausscomp4[6* comp_index - 1], gausscomp5[6* comp_index - 1], gausscomp6[6* comp_index - 1])
	peak_err     <-  c(gausscomp1[6* comp_index], gausscomp2[6* comp_index], gausscomp3[6* comp_index], gausscomp4[6* comp_index], gausscomp5[6* comp_index], gausscomp6[6* comp_index])
	return(data.frame(mjd, velocity, width, peak, velocity_err, width_err, peak_err))	
}


#---- compplot( velocity, gausscomp ) : A function to plot Gaussian components
compLorentzPlot <- function( velocity, comp ){
	#  velocity : a vector of velocities to plot
	#  comp : Lorentzian Function components … num_comp x 3 (mean velocity, width, and peak intensity) x 4 (estimate, sd, t-value, p-value)
	num_comp <- length( comp ) / 12   # 12 parameters per component
	for( comp_index in 1:num_comp ){
		peak <- comp[3* comp_index ]
		mean <- comp[3* comp_index - 2]
		width <- comp[3*comp_index - 1]
		lines(velocity, peak / (((velocity - mean)/width)^2 + 1 ), col=comp_index, lty='dashed')
		cat(sprintf("Comp%d : VLSR=%6.1f km/s FWHM=%6.2f Peak=%6.2f mJy\n", comp_index, mean, width, 1e3*peak))
	}
}
