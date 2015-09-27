## This script contains functions for NC Lichen FD project

####################################################################
### FUNCTIONS

# A function that calculate topographic openness (after Yokoyama et al 2002)
calc_openness = function(x, dem, L){
	
	# Make profile lines from focal point in cardinal directions (0,45,90, etc..)
	Ds = seq(0, 315, 45)
	endpts = destPoint(x, Ds, L)
	profiles = sapply(1:nrow(endpts), function(i){
		y = endpts[i,]
		Dname = Ds[i]
		lines = Lines(Line(rbind(x,y)), ID=Dname)
	})
	profiles = SpatialLines(profiles, proj4string=CRS('+proj=longlat'))
	
	# Extract cell numbers and values from dem along lines
	elev_profile = extract(dem, profiles)
	focal_elev = extract(dem, x)
	elev_profile = lapply(elev_profile, function(p) p - focal_elev)
	
	# Calculate distances to focal cell along profiles
	dists = distanceFromPoints(dem, x)
	dist_profile = extract(dists, profiles)

	# Calculate angle of elevation along each profile
	angles = sapply(1:8, function(p) atan(elev_profile[[p]]/dist_profile[[p]])/(2*pi)*360)
	
	# Drop 0s that arise from profiles including the focal cell
	angles = lapply(angles, function(p){
		new_p = p
		zeros = which(p==0)
		if(zeros %in% c(1, length(p))) new_p = p[-zeros]
		new_p
	})	
	
	# Find maximum and minimum angle along each profile
	maxs = lapply(angles, max)
	mins = lapply(angles, min)
	
	# Calculate zeniths and nadirs
	zeniths = 90-unlist(maxs)
	nadirs = 90+unlist(mins)
	
	# Calculate openness
	open_pos = mean(zeniths)
	open_neg = mean(nadirs)

	data.frame(OpenPos=open_pos, OpenNeg=open_neg)
}


## A function that matches a lichen to its morphotype
# x is a vector with names of traits that are included in the morphs table
# morphs is a dataframe of traits corresponding to morphotypes. Rownames should be MorphIDs.

match_morpho = function(x, morphs){
	traitnames = colnames(morphs)
	use_x = data.frame(x[traitnames])

	morphs[is.na(morphs)] = 'NA'
	use_x[is.na(use_x)] = 'NA'
	
	keep_morphs = morphs
	for(i in 1:length(traitnames)){
		keep_morphs = subset(keep_morphs, keep_morphs[,i] == use_x[,i])
	}
	

	#matches = sapply(1:nrow(morphs), function(i){
	#	sum(morphs[i,]==use_x)
	#})
	
	this_morph = ifelse(nrow(keep_morphs)==1, rownames(keep_morphs), NA)

	this_morph
}

# A function that calculates Rao's Q (or mpd)- returns 0 for one taxon
# dmat : a distance matrix between species
# abun : a relative abundance matrix (rows sum to 1)
# dmat and abun must be in the same order
calc_rao = function(abun, dmat){
	N = nrow(dmat)
	UP = upper.tri(matrix(1, nrow=N, ncol=N), diag=F)
	
	sapply(1:nrow(abun), function(i){
		sum(abun[i,]%*%t(abun[i,])*dmat*UP)
	})
}

# A function that calculates the percentile of an observation x in null distribution null
calc_p = function(x, null){
	# Remove any na values
	null = null[!is.na(null)]

	min(sum(x<=null), sum(x>=null))/length(null)
}

# A function that calculates the z-score of an observation x in null distribution null
calc_z = function(x, null){
	# Remove any na values
	null = null[!is.na(null)]

	(x - mean(null))/sqrt(var(null))
}


# A function that calculated the fixed effects variance of a lmm model
calc_fixedvar = function(x){
	var(as.numeric(as.vector(fixef(x)) %*% t(model.matrix(x))))
}

# A function that extracts all variance components from a lmm model
get_allvar = function(x){
	var_df = data.frame(VarCorr(x))
	varcomps = c(calc_fixedvar(x), var_df$vcov)
	names(varcomps) = c('Fixed',var_df$grp)
	varcomps
}

# A function to calculate total variance from a glmm
# Note that function r.squaredGLMM(x) in MuMIn now does this R2 calculation
# If using glmm with poisson, make sure that dummy variable with 1 level for each observation is included as random intercept
#	> This is how residual var is incorporated into sigma.random
calc_totvar = function(mod){


	if(class(mod)=='glmerMod'){
		# Drop fixed effects to get null model
		mod_null = mod
		drop.terms = attr(terms(mod), 'term.labels')
		if(length(drop.terms)>0){
			drop.form = paste('.~.-', paste(drop.terms, collapse='-'), sep='')
			mod_null = update(mod, drop.form, data=mod$model)
		}

		family = mod@resp$family$family
		link = mod@resp$family$link

		if(family=='poisson'&link=='log'){
			sigma.mod = log(1 + (1/exp(as.numeric(fixef(mod_null)))))
		}
		
		if(family=='binomial'&link=='logit'){
			sigma.mod = (pi^2)/3
		}
	}

	if(class(mod)=='lmerMod'){
		sigma.mod = 0
	}

	sigma.fixed = calc_fixedvar(mod)
	sigma.random = sum(data.frame(VarCorr(mod))$vcov) + ifelse(class(mod)=='lmerMod', attr(VarCorr(mod), 'sc')^2,0)
	
	sigma.total = sigma.fixed + sigma.random + sigma.mod # Error is incorporated into sigma.random term
	sigma.total
}


# This function creates a response variable for logistic regression P(occurance | # morphos)
# i is the trait
# lichdata is a dataframe with the trait value for each morphotype in each sample
calc_bin_pres = function(i, lichdata){
	tab = xtabs(~lichdata$SampID+lichdata[,i])
	tot = rowSums(tab)
	data.frame(Trues = tab[,2], Tot = tot)
}

# This function creates a response variable for logistic regression P(occurance | abundance)
# i is the trait
# lichdata is a dataframe with the trait value for each morphotype in each sample
# lichdata must include morphotype abundance as AbunCount
calc_bin_abun = function(i, lichdata){
	tab = xtabs(AbunCount~SampID+lichdata[,i], data=lichdata)
	tot = rowSums(tab)
	data.frame(Trues = tab[,2], Tot = tot)
}


# This function creates a response variable for multinomial regression P(occurance | # morphos)
# i is the trait
# lichdata is a dataframe with the trait value for each morphotype in each sample
calc_cat_pres = function(i, lichdata){
	tab = xtabs(~lichdata$SampID+lichdata[,i])
	tab
}

# This function creates a response variable for multinomial regression P(occurance | abundance)
# i is the trait
# lichdata is a dataframe with the trait value for each morphotype in each sample
# lichdata must include morphotype abundance as AbunCount
calc_cat_abun = function(i, lichdata){
	tab = xtabs(AbunCount~SampID+lichdata[,i], data=lichdata)
	tab
}

## A function that partitions variation between two models
## Must be supplied with R2 for [A, B, AB]
partvar2 = function(R2s){
	
	a = R2s[3]-R2s[2]
	c = R2s[3]-R2s[1]
	b = R2s[1] - a
	d = 1-R2s[3]

	this_partition = c(a,b,c,d)
	names(this_partition) = c(names(R2s)[1], 'Both', names(R2s)[2], 'Unexplained')

	this_partition
}


## A function to plot a color bar
## from: http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut))/(max-min)

    par(mar=c(0.5,6,0,0))
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='')
    axis(2, las=1)
    for (i in 1:(length(lut))) {
     y = (i-1)/scale + min
     rect(0,y,3,y+1/scale, col=lut[i], border=NA)
    }
}

