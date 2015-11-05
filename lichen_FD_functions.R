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

# A function that calculates trait diversity for a single trait
# Function will drop all missing observations.
# x : vector of trait values
# N : a vector of morphospecies frequencies, or missing if diversity not abundance-weighted
# Group : factor vector indicating which groups each observation belongs to
# vartype : string indicating whether trait has numeric or categorical values
calc_trait_div = function(x, N=NA, Group, vartype){
	# Drop missing observations
	Group = Group[!is.na(x)]
	N = N[!is.na(x)]
	x = x[!is.na(x)]

	# Frequency of each character state
	if(!is.na(N[1])){
		abun_mat = xtabs(N~Group+x)
	} else {
		abun_mat = as.matrix(table(Group, x))
	}
	
	# Relative frequency
	prop_mat = abun_mat/rowSums(abun_mat)

	# Euclidean distance between character states
	if(vartype %in% c('numeric','ordered')){
		d_mat = as.matrix(dist(as.numeric(colnames(abun_mat))))
	}
	if(vartype %in% c('factor', 'symm', 'asymm')){
		nfact = length(levels(factor(x)))
		d_mat = matrix(1, nrow=nfact, ncol=nfact)
		diag(dmat) = 0
	}

	# Calculat mpd
	mpd = calc_rao(prop_mat, d_mat)

	names(mpd) = rownames(abun_mat)
	mpd
}

# A function that calculates a null distribution for trait diversity
# reps :  number of randomizations
# see other parameters above
calc_trait_div_null = function(reps, x, N, Group, vartype){
	null_dist = sapply(1:reps, function(i){
		neworder = sample(x, length(x), replace=F)
		newvals = calc_trait_div(neworder, N, Group, vartype)
		return(newvals)
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


# A function that calculated the fixed effects variance of a lmm or lm model
calc_fixedvar = function(x){
	if(class(x)=='glmerMod' | class(x)=='lmerMod'){
		fixed.var = var(as.numeric(as.vector(fixef(x)) %*% t(model.matrix(x))))
	}
	if('lm' %in% class(x)){
		fixed.var = var(as.numeric(as.vector(coef(x)) %*% t(model.matrix(x))))
	}
	fixed.var	
}

# A function that extracts all variance components from a lmm model
get_allvar = function(x){
	var_df = data.frame(VarCorr(x))
	varcomps = c(calc_fixedvar(x), var_df$vcov)
	names(varcomps) = c('Fixed',var_df$grp)
	varcomps
}

# A function to calculates repeatability from a mixed effects model
# > based on equations in Nakagawa & Schielzeth 2010, Biological Reviews 85: 936-956
# Assumes models have been fit with additive overdispersion
# Can also calculate adjusted repeatability
# mod : fitted model object
# focalID : character name of random effect to test
# obsID : character name of random effect whose variance gives additive error variance
# adjID  : character vector of random and fixed effects to be controlled for before calculating R
calc_rpt = function(mod, focalID, adjID=NA, obsID){

	if(class(mod)=='lmerMod'){
		sigma.mod = 0
	} else {

		# Drop fixed effects to get null model for calculating distribution specific variance of log link
		mod_null = mod
		drop.terms = attr(terms(mod), 'term.labels')
		if(length(drop.terms)>0){
			drop.form = paste('.~.-', paste(drop.terms, collapse='-'), sep='')
			mod_null = update(mod, drop.form) #, data=mod@frame)
		}

		# If model is not Gaussian, get family and link function for calculating distribution specific variance
		if(class(mod)=='glmerMod'){
			family = mod@resp$family$family
			link = mod@resp$family$link
		}
		if('glm' %in% class(mod)){
			family = mod$family$family
			link = mod$family$link
		}

		if(family=='poisson'&link=='log'){
			if(class(mod)=='glmerMod') sigma.mod = log(1 + (1/exp(as.numeric(fixef(mod_null)))))
			if(class(mod)=='glm') sigma.mod = log(1 + (1/exp(as.numeric(coef(mod_null)))))
		}
		if(family=='poisson'&link=='sqrt') sigma.mod = 0.25
		if(family=='binomial'&link=='logit') sigma.mod = (pi^2)/3
		if(family=='gaussian'&link=='identity') sigma.mod = 0
	}

	all.var = get_allvar(mod)
	R.focal = all.var[focalID]/(sum(all.var) + sigma.mod)

	if(!is.na(adjID)){
		R.focal.adj = all.var[focalID]/(sum(all.var[c(focalID,obsID)]) + sigma.mod)
		R.adj = sum(all.var[adjID])/(sum(all.var) + sigma.mod)
		R.adj.focal = sum(all.var[adjID])/(sum(all.var[c(adjID,obsID)]) + sigma.mod)

		data.frame(R.focal, R.focal.adj, R.adj, R.adj.focal)
	} else {
		R.focal
	}
	
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

