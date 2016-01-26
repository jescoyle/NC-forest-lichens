## This script assesses variation in single traits for the Lichen FD Project

git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/NC-forest-lichens/'

# Load data & make data frames for analysis
source(paste(git_dir, 'load_data.R', sep=''))
source(paste(git_dir,'lichen_FD_functions.R', sep=''))

# Colors
repmode_col = c('#ff6600','#5f8dd3')

# Subset to plots only in Piedmont and Mountains
lichens_pm = subset(lichens, SiteID!='Bladen')
lichens_pm = droplevels(lichens_pm)
plots_pm = subset(plot_data, SiteID!='Bladen')
samples_pm = subset(samples, SiteID!='Bladen')
samples_pm = droplevels(samples_pm)

# Color ramp
blue2red = apply(read.csv('../../blue2red_10colramp.txt'), 1, function(x) rgb(x[1],x[2],x[3],maxColorValue=255))
blue2red = rev(blue2red)


################################################################
### Check co-variance between environmental predictors
library(lattice)

sampvars = c('Height','Angle','Bryophytes','Shedding','FurrowDepth','Aspect','pH','Density','WaterCapacity')
treevars = c('DBH','Trans_tot_cor','TreeTaxonID')

## Make histograms of each variable to check values
pdf('./Figures/Sample variable histograms.pdf', height=6, width=6)
par(mar=c(5,5,1,1))
for(x in sampvars){
	if(is.numeric(samples[,x])){
		hist(samples[,x], xlab=x, las=1, main='')
	} else {
		plot(table(samples[,x]), xlab=x, ylab='# Samples', las=1) 
	}
}
dev.off()

tree_data$PlotID = factor(tree_data$PlotID, levels=plot_data[order(plot_data$PairID, plot_data$TopoPos),'PlotID'])

pdf('./Figures/Tree variable distributions.pdf', height=6, width=9)
for(x in treevars){
	if(is.numeric(tree_data[,x])){
		par(mar=c(5,5,1,1))
		hist(tree_data[,x], xlab=x, las=1, main='')
	} else {
		use_tab = table(tree_data[,x])
		use_tab = use_tab[order(use_tab)]
		par(mar=c(10,5,1,1))
		plot(use_tab, ylab='# Trees', xlab='', las=2, axes=F)
		axis(1, at=1:length(use_tab), labels=names(use_tab), las=2)
		axis(2) 
	}
}

par(mar=c(5,5,1,1))
histogram(~DBH|PlotID, data=tree_data)
boxplot(DBH~PlotID, data=tree_data, las=2)

histogram(~Trans_tot_cor|PlotID, data=tree_data)
boxplot(Trans_tot_cor~PlotID, data=tree_data, las=2)

dev.off()

## Covariance between sample and tree level variables
pdf('./Figures/Sample-tree variable covariance.pdf', height=6, width=10)
for(x in sampvars){
	if(is.numeric(samples[,x])){
		par(mfrow=c(1,2))
		par(mar=c(5,5,1,1))
		plot(samples[,x]~DBH, ylab=x, las=1, main='', data=samples)
		plot(samples[,x]~Trans_tot_cor, ylab=x, las=1, main='', data=samples)
		
		use_mns = aggregate(samples[,x], by=list(TreeTaxon=samples$TreeTaxonID), FUN=function(x) mean(x, na.rm=T))
		use_levels = use_mns[order(use_mns$x),'TreeTaxon']
		par(mfrow=c(1,1))
		par(mar=c(10,5,1,1))
		boxplot(samples[,x]~factor(samples$TreeTaxonID, levels=use_levels), las=2, ylab=x)
	} else {
		par(mfrow=c(1,2))
		par(mar=c(5,5,1,1))
		boxplot(DBH~samples[,x], data=samples, las=2, ylab='DBH', xlab=x)	
		boxplot(Trans_tot_cor~samples[,x], data=samples, las=2, ylab='Trans_tot_cor', xlab=x)	
		
		use_mns = aggregate(as.numeric(samples[,x]), by=list(TreeTaxon=samples$TreeTaxonID), FUN=function(x) mean(x, na.rm=T))
		use_levels = use_mns[order(use_mns$x),'TreeTaxon']
		par(mfrow=c(1,1))
		par(mar=c(10,5,1,1))
		boxplot(as.numeric(samples[,x])~factor(samples$TreeTaxonID, levels=use_levels), las=2, ylab=x)
	}
}
dev.off()

## Covariance between plot level variables and sample/tree level variables
plotvars = c('RH','Soil_pH','CloudFreq_mean','CloudFreq_sd','AP','T_max','T_mean','VPD_max','Elevation','OpenPos')

pdf('./Figures/Sample-tree-plot variable covariances.pdf', height=6, width=8)
for(x in plotvars){
for(y in sampvars){
	par(mar=c(5,5,1,1))
	if(is.numeric(samples[,y])){
		plot(samples[,y]~samples[,x], xlab=x, ylab=y, las=1)
	} else {
		boxplot(samples[,x]~samples[,y], ylab=x, xlab=y, las=1)
	}
}

	plot(DBH~samples[,x], data=samples, xlab=x, las=1)
	plot(Trans_tot_cor~samples[,x], data=samples, xlab=x, las=1)

	use_mns = aggregate(as.numeric(samples[,x]), by=list(TreeTaxon=samples$TreeTaxonID), FUN=function(x) mean(x, na.rm=T))
	use_levels = use_mns[order(use_mns$x),'TreeTaxon']
	par(mar=c(10,5,1,1))
	boxplot(samples[,x]~factor(samples$TreeTaxonID, levels=use_levels), las=2, ylab=x)
}
for(y in sampvars){
	par(mar=c(5,5,1,1))
	plot(samples[,y]~Ecoregion, data=samples, ylab=y)
	plot(samples[,y]~TopoPos, data=samples, ylab=y)
}
dev.off()

## Covariance among plot-level variables

pdf('./Figures/Plot variable correlations.pdf', height=8, width=8)
pairs(plot_data[,plotvars], lower.panel=panel.smooth, upper.panel=panel.cor)

for(x in plotvars){
	boxplot(plot_data[,x]~plot_data$Ecoregion, las=1, ylab=x)
}
dev.off()

cor(plot_data[,climvars])

climvars = c('RH','CloudFreq_mean','CloudFreq_sd','AP','T_max','T_mean','VPD_max','Elevation')
plotpca = prcomp(plot_data[,climvars], scale=T)
summary(plotpca)
cor(plot_data[,climvars])


# Makes sense just to use Elevation as a proxy for temp and humidity - both are nearly linear functions
plot(CloudFreq_mean~Elevation, data=plot_data)
plot(T_mean~Elevation, data=plot_data)


################################################################
### Maximum Likelihood

library(lme4)
library(R.utils) # capitalize
library(reshape2)
library(lattice)
library(MuMIn)

###########################################
#### Water traits ~ water availability ####
# Focal traits: P(crustose), P(fruticose), P(cilia|foliose), Photobiont, Attachment, LobeDissect
# Focal predictors: Sample: WaterCapacity, Bryophytes  Plot: CloudFreq_mean, VPD_max, OpenPos 

## Binary Traits

# Create response variable for probability of crutose and fruticose forms
sum(is.na(lichens_pm$Form)) # should be 0
form_abun = xtabs(AbunCount ~ SampID + Form, data=lichens_pm)
form_pres = xtabs(~ SampID + Form, data=lichens_pm)

crust_abun = cbind(form_abun[,'crustose'], rowSums(form_abun))
crust_pres = cbind(form_pres[,'crustose'], rowSums(form_pres))
frut_abun = cbind(form_abun[,'fruticose'], rowSums(form_abun))
frut_pres = cbind(form_pres[,'fruticose'], rowSums(form_pres))

# Create response variables for P(cilia|foliose)
cilia_pres = xtabs(~SampID + Cilia, data=lichens_pm)
cilia_pres = cbind(cilia_pres[,1], rowSums(cilia_pres))
cilia_pres = cilia_pres[cilia_pres[,2]>0,] # Drop samples without foliose lichens
cilia_abun = xtabs(AbunCount~SampID + Cilia, data=lichens_pm)
cilia_abun = cbind(cilia_abun[,1], rowSums(cilia_abun))
cilia_abun = cilia_abun[cilia_abun[,2]>0,]

# Create response variable for P(cyanobacteria)
cyano_abun = xtabs(AbunCount ~ SampID + Photobiont, data=lichens_pm)
cyano_abun = cbind(cyano_abun[,'cyano'], rowSums(cyano_abun))

# Calculate sample-mean traits for numeric traits
# Create relative presence/abundance community data matrices
num_traits = c('Attachment','LobeDissect')

sampXmorph_pm = sampXmorph[samples[rownames(sampXmorph),'SiteID']!='Bladen',]

sampXmorph_pres = sampXmorph_pm>0
sampXmorph_abun = sampXmorph_pm

# Calculate sample mean traits based on pres-abs vs. abun
sampXtrait_pres = sapply(num_traits, function(i){
	y = unclass(morphos[,i]) # Note: this causes Attachment to be +1 from values in lichens
	comm = sampXmorph_pres[,rownames(morphos)]

	# Ignore lichens where the trait is missing
	comm[,is.na(y)] = 0
	y[is.na(y)] = 0

	# Calculate relative abundance
	comm = as.matrix(comm/rowSums(comm)	)
	
	# Calculate community-weighted mean
	comm%*%as.numeric(y)
})
rownames(sampXtrait_pres) = rownames(sampXmorph_pm)

sampXtrait_abun = sapply(num_traits, function(i){
	y = unclass(morphos[,i]) # Note: this causes Attachment to be +1 from values in lichens
	comm = sampXmorph_abun[,rownames(morphos)]

	# Ignore lichens where the trait is missing
	comm[,is.na(y)] = 0
	y[is.na(y)] = 0

	# Calculate relative abundance
	comm = as.matrix(comm/rowSums(comm)	)
	
	# Calculate community-weighted mean
	comm%*%as.numeric(y)
})
rownames(sampXtrait_abun) = rownames(sampXmorph_pm)

# Previous analysis showed that these traits are best modeled after log-transforming to account for strict positivity and increasing variance with mean
sampXtrait_pres = log(sampXtrait_pres)
sampXtrait_abun = log(sampXtrait_abun)

# Make dataframe of predictors
scale_vars = c('SampID','PlotID','SiteID','Ecoregion')
pred_vars1 = c('WaterCapacity','Bryophytes')
pred_vars2 = c('CloudFreq_mean','VPD_max','OpenPos')
Xdata = samples[,c(scale_vars, pred_vars1, pred_vars2)]

# Three WHC appear to be outliers, so set to NA
Xdata[with(Xdata, WaterCapacity<0.1&!is.na(WaterCapacity)),'WaterCapacity'] = NA

# Log-transform WaterCapacity because it is a ratio of water held / unit bark by mass
Xdata$WaterCapacity = log(Xdata$WaterCapacity)

# Convert Bryophytes to integer
Xdata$Bryophytes = unclass(Xdata$Bryophytes) - 1

# Center climate predictors so that marginal effects can be more easily interpreted
climMeans = colMeans(plot_data[rownames(plot_data)!='Bladen1',pred_vars2])
Xdata[,pred_vars2] = Xdata[,pred_vars2] - rep(climMeans, each=nrow(Xdata))

# Examine predictors: no major skew or outliers
par(mfrow=c(1,5))
for(i in c(pred_vars1, pred_vars2)) hist(unclass(Xdata[,i]), main=i)


# Set response variable
Y = cyano_abun # crust_abun, frut_abun, cilia_abun, cyano_abun

mods = vector('list', length(pred_vars1)*length(pred_vars2))
dim(mods) = c(length(pred_vars1),length(pred_vars2))
dimnames(mods) = list(pred_vars1, pred_vars2)

ests = array(NA, dim=c(length(pred_vars1), length(pred_vars2), 8, 6), 
	dimnames = list(pred_vars1, pred_vars2, c('b1','b2','b12','b21','b1X2','b_eco','b1Xeco','nobs'), c('est','low95','up95','P','R2','Dev_resid')))

## Binary variables
for(i in pred_vars1){
for(j in pred_vars2){
	env1 = Xdata[rownames(Y),i]
	env2 = Xdata[rownames(Y),j]
	Ecoregion = Xdata[rownames(Y),'Ecoregion']

	# Remove data with missing observations
	missing=is.na(rowSums(Y))|is.na(env1)|is.na(env2)
	use_Y = Y[!missing,]
	env1 = env1[!missing]
	env2 = env2[!missing]
	Ecoregion = Ecoregion[!missing]

	# Make models
	mod0 = glm(use_Y ~ 1, family='binomial')
	mod1 = glm(use_Y ~ env1, family='binomial')
	mod2 = glm(use_Y ~ env2, family='binomial')
	mod = glm(use_Y ~ env1 + env2, family='binomial')
	mod_X = glm(use_Y ~ env1*env2, family='binomial')
	mod_eco = glm(use_Y ~ env2 + Ecoregion, family='binomial')
	mod_ecoX = glm(use_Y ~ env1*Ecoregion, family='binomial')

	mods[i,j][[1]] = mod_X
	
	# Test terms
	p_inter = anova(mod, mod_X, test='Chisq')
	p_1 = anova(mod1, mod0, test='Chisq')	
	p_2 = anova(mod2, mod0, test='Chisq')
	p_12 = anova(mod, mod2, test='Chisq')
	p_21 = anova(mod, mod1, test='Chisq')
	p_eco = anova(mod_eco, mod2, test='Chisq')
	p_1Xeco = anova(mod_ecoX, update(mod_ecoX, .~.-env1:Ecoregion), test='Chisq')	
	test_list = list(p_1,p_2,p_12,p_21,p_inter,p_eco,p_1Xeco)

	# Save estimates
	ests[i,j,c('b1','b2','b1X2','b_eco','b1Xeco'),'est'] = c(coef(mod1)[2], coef(mod2)[2], 
		coef(mod_X)[4], coef(mod_eco)[3], coef(mod_ecoX)[4])
	#if(p_inter[2,'Pr(>Chi)'] < 0.05){
	#	ests[i,j,c('b12','b21'),'est'] = coef(mod_X)[2:3]
	#	ests[i,j,c('b12','b21'),c('low95','up95')] = confint(mod_X)[2:3,]
	#} else {
		ests[i,j,c('b12','b21'),'est'] = coef(mod)[2:3]
		ests[i,j,c('b12','b21'),c('low95','up95')] = confint(mod)[2:3,]
	#}
	
	ests[i,j,c('b1','b2'),c('low95','up95')] = t(sapply(list(mod1, mod2), function(x) confint(x)[2,]))
	ests[i,j,'b1X2',c('low95','up95')] = confint(mod_X)[4,]
	ests[i,j,'b_eco',c('low95','up95')] = confint(mod_eco)[3,]
	ests[i,j,'b1Xeco',c('low95','up95')] = confint(mod_ecoX)[4,]
	ests[i,j,'nobs','est'] = nrow(use_Y)
	ests[i,j,,'P'] = c(sapply(test_list, function(x) x[2,'Pr(>Chi)']), NA)
	ests[i,j,,'Dev_resid'] = c(sapply(list(mod1, mod2, mod, mod, mod_X, mod_eco, mod_ecoX), deviance), NA)
	ests[i,j,,'R2'] = c(sapply(list(mod1, mod2, mod, mod, mod_X, mod_eco, mod_ecoX), function(x) r.squaredLR(x)), NA)
	
}}

# Save results
#crust_ests = ests
#crust_mods = mods
#frut_ests = ests
#frut_mods = mods
#cilia_ests = ests
#cilia_mods = mods
#cyano_ests = ests
#cyano_mods = mods

## Numeric variables
Y = sampXtrait_abun[,'Attachment']

mods = vector('list', length(pred_vars1)*length(pred_vars2))
dim(mods) = c(length(pred_vars1),length(pred_vars2))
dimnames(mods) = list(pred_vars1, pred_vars2)

ests = array(NA, dim=c(length(pred_vars1), length(pred_vars2), 8, 6), 
	dimnames = list(pred_vars1, pred_vars2, c('b1','b2','b12','b21','b1X2','b_eco','b1Xeco','nobs'), c('est','low95','up95','P','R2','Dev_resid')))

for(i in pred_vars1){
for(j in pred_vars2){
	env1 = Xdata[names(Y),i]
	env2 = Xdata[names(Y),j]
	Ecoregion = droplevels(Xdata[names(Y),'Ecoregion'])

	# Remove data with missing observations
	missing=is.na(Y)|is.na(env1)|is.na(env2)
	use_Y = Y[!missing]
	env1 = env1[!missing]
	env2 = env2[!missing]
	Ecoregion = Ecoregion[!missing]
	
	# Make models
	mod0 = glm(use_Y ~ 1, family='gaussian')
	mod1 = glm(use_Y ~ env1, family='gaussian')
	mod2 = glm(use_Y ~ env2, family='gaussian')
	mod = glm(use_Y ~ env1 + env2, family='gaussian')
	mod_X = glm(use_Y ~ env1*env2, family='gaussian')
	mod_eco = glm(use_Y ~ env2 + Ecoregion, family='gaussian')
	mod_ecoX = glm(use_Y ~ env1*Ecoregion, family='gaussian')

	mods[i,j][[1]] = mod_X
	
	# Test terms
	p_inter = anova(mod, mod_X, test='Chisq')
	p_1 = anova(mod1, mod0, test='Chisq')	
	p_2 = anova(mod2, mod0, test='Chisq')
	p_12 = anova(mod, mod2, test='Chisq')
	p_21 = anova(mod, mod1, test='Chisq')
	p_eco = anova(mod_eco, mod2, test='Chisq')
	p_1Xeco = anova(mod_ecoX, update(mod_ecoX, .~.-env1:Ecoregion), test='Chisq')	
	test_list = list(p_1,p_2,p_12,p_21,p_inter,p_eco,p_1Xeco)

	# Save estimates
	ests[i,j,c('b1','b2','b1X2','b_eco','b1Xeco'),'est'] = c(coef(mod1)[2], coef(mod2)[2], 
		coef(mod_X)[4], coef(mod_eco)[3], coef(mod_ecoX)[4])
	#if(p_inter[2,'Pr(>Chi)'] < 0.05){
	#	ests[i,j,c('b12','b21'),'est'] = coef(mod_X)[2:3]
	#	ests[i,j,c('b12','b21'),c('low95','up95')] = confint(mod_X)[2:3,]
	#} else {
		ests[i,j,c('b12','b21'),'est'] = coef(mod)[2:3]
		ests[i,j,c('b12','b21'),c('low95','up95')] = confint(mod)[2:3,]
	#}
	
	ests[i,j,c('b1','b2'),c('low95','up95')] = t(sapply(list(mod1, mod2), function(x) confint(x)[2,]))
	ests[i,j,'b1X2',c('low95','up95')] = confint(mod_X)[4,]
	ests[i,j,'b_eco',c('low95','up95')] = confint(mod_eco)[3,]
	ests[i,j,'b1Xeco',c('low95','up95')] = confint(mod_ecoX)[4,]
	ests[i,j,'nobs','est'] = length(use_Y)
	ests[i,j,,'P'] = c(sapply(test_list, function(x) x[2,'Pr(>Chi)']), NA)
	ests[i,j,,'Dev_resid'] = c(sapply(list(mod1, mod2, mod, mod, mod_X, mod_eco, mod_ecoX), deviance), NA)
	ests[i,j,,'R2'] = c(sapply(list(mod1, mod2, mod, mod, mod_X, mod_eco, mod_ecoX), function(x) r.squaredLR(x)), NA) # note that in scale models I did not use adj R2
	
}}

#attach_mods = mods
#attach_ests = ests
#dissect_mods = mods
#dissect_ests = ests

save(crust_ests, crust_mods, frut_ests, frut_mods, cilia_ests, cilia_mods, cyano_ests, cyano_mods,
	attach_ests, attach_mods, dissect_ests, dissect_mods, file='water_trait_abun_models.RData')

load('water_trait_abun_models.RData')


## Make summary tables
parm_list = list(cyano_ests, crust_ests, frut_ests, cilia_ests, attach_ests, dissect_ests)
names(parm_list) = c('ProbCyano','ProbCrustose','ProbFruticose','Cilia','Attachment','LobeDissect')

# Combine all estimates into a big array
library(abind)
library(plyr)

all_ests = abind(cyano_ests, crust_ests, frut_ests, cilia_ests, attach_ests, dissect_ests, along=0.5, new.names=names(parm_list))
names(dimnames(all_ests)) = c('trait','pred.samp','pred.plot','parm','quant')

## Use corrplot to plot coef estimates
#library(corrplot)

use_col = c('red3','blue3')
use_pch = 21

molten = melt(all_ests)

# Interactions between sample- and larger-scale env predictors
inter = acast(molten, trait~pred.samp+pred.plot~quant, subset=.(parm=='b1X2' & quant%in%c('est','P')))

# Marginal effects of sample controlling for climate
tab = acast(molten, trait~pred.samp+pred.plot~quant, subset=.(parm=='b12' & quant%in%c('est','P','low95','up95')))

svg('./Figures/water_traits_marginal_effects_bark.svg', height=7.5, width=5)
par(mar=c(.5,0,.5,0))
par(oma=c(0,2,0,0))
par(mfrow=dim(tab)[1:2])

for(i in 1:dim(tab)[1]){
for(j in 1:dim(tab)[2]){

	this_est = tab[i,j,'est']
	if(j %in% 1:3) ylim = max(abs(range(tab[i,1:3,c('low95','up95')])))
	if(j %in% 4:6) ylim = max(abs(range(tab[i,4:6,c('low95','up95')])))
	ylim = ylim*c(-1.05,1.05)
	plot(0, this_est, axes=F, ylim=ylim, type='n')
	segments(par('usr')[1]*ifelse(j%%3==1, 0.7, 1), 0, par('usr')[2]*ifelse(j%%3==0, 0.7, 1), 0, lty=2)

	segments(0,tab[i,j,'low95'],0,tab[i,j,'up95'], 
		col=ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]))
	points(0, this_est, pch=use_pch, cex=2, 
		bg=ifelse(tab[i,j,'P']<0.05, ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]),'white'),
		col=ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]))
	if(j %in% c(1,4)) axis(2, line=-2, las=1)
	if(j==1) mtext(dimnames(tab)[[1]][i], 2, 1)
}}
dev.off()

# Marginal effects of climate controlling for sample
tab = acast(molten, trait~pred.plot+pred.samp~quant, subset=.(parm=='b21' & quant%in%c('est','P','low95','up95')))
inter = acast(molten, trait~pred.plot+pred.samp~quant, subset=.(parm=='b1X2' & quant%in%c('est','P')))

svg('./Figures/water_traits_marginal_effects_climate.svg', height=7.5, width=5)
par(mar=c(.5,0,.5,0))
par(oma=c(0,2,0,0))
par(mfrow=dim(tab)[1:2])

for(i in 1:dim(tab)[1]){
for(j in 1:dim(tab)[2]){
	
	this_est = tab[i,j,'est']
	if(j %in% 1:2) ylim = max(abs(range(tab[i,1:2,c('low95','up95')])))
	if(j %in% 3:4) ylim = max(abs(range(tab[i,3:4,c('low95','up95')])))	
	if(j %in% 5:6) ylim = max(abs(range(tab[i,5:6,c('low95','up95')])))
	ylim = ylim*c(-1.05,1.05)
	plot(0, this_est, axes=F, ylim=ylim, type='n')
	segments(par('usr')[1]*ifelse(j%%2==1, 0.7, 1), 0, par('usr')[2]*ifelse(j%%2==0, 0.7, 1), 0, lty=2)

	segments(0,tab[i,j,'low95'],0,tab[i,j,'up95'], 
		col=ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]))
	points(0, this_est, pch=use_pch, cex=2, 
		bg=ifelse(tab[i,j,'P']<0.05, ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]),'white'),
		col=ifelse(inter[i,j,'P']>=0.05, 'black', use_col[(inter[i,j,'est']>0)+1]))

	if(j %in% c(1,3,5)) axis(2, line=-2, las=1)
	if(j==1) mtext(dimnames(tab)[[1]][i], 2, 1)
	
}}
dev.off()


# Explanatory ability of models
bark_tab = dcast(molten, trait+pred.samp~quant, subset=.(parm=='b1' & pred.plot=='VPD_max')) # Doesn't matter which pred.plot specified
subset(bark_tab, P < 0.05)
clim_tab = dcast(molten, trait+pred.plot~quant, subset=.(parm=='b2' & pred.samp=='Bryophytes')) # Doesn't matter which pred.samp specified
subset(clim_tab, P < 0.05)

# Are there significantly different effect of bark scale variables across Ecoregions?
all_ests[,,'CloudFreq_mean','b1Xeco',c('est','P')] 
# for the effect of Bryophytes on Attachment: more positive effect of bryophytes on thallus height in Piedmont (drier)
# for the effects of WHC and Bryophytes on prob. cyano: more positive effects in Piedmont (drier)

# Explore interactions
# Do any models have a significant interaction between bark and plot scale variables?
inter[,,'P']<0.05

# Crustose: Bryophytes & OpenPos
i = 'Bryophytes'
j = 'OpenPos'

mod = crust_mods[i,j][[1]]
mod_func = function(x, level){	
	modlevel = level - climMeans[j]
	predict(mod, data.frame(env1=x, env2=modlevel), type='response')
}

yvals = crust_abun[,1] / crust_abun[,2]
climrange = range(plot_data[rownames(plot_data)!='Bladen1', j])


par(mfrow=c(1,2))
plot(yvals~Xdata[names(yvals), i], las=1, ylab='P(Crustose)', xlab=i, ylim=c(0,1), pch=16, col='#00000050')
curve(mod_func(x, climrange[1]), from=0, to=5, add=T, lwd=2, col=use_col[1])
curve(mod_func(x, climrange[2]), from=0, to=5, add=T, lwd=2, col=use_col[2])

plot(yvals~I(Xdata[names(yvals), j]+climMeans[j]), las=1, ylab='P(Crustose)', xlab=j, ylim=c(0,1), pch=16, col='#00000050')
curve(mod_func(x=0, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[1])
curve(mod_func(x=5, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[2])


## Figure 1: Attachment vs. Bryophytes with interactions
use_col = c('#0000CD','#CD0000') # 'blue3', 'red3'
use_col_trans = paste(use_col, '88', sep='')
j='VPD_max'
yvals = exp(sampXtrait_abun[,'Attachment'])
xvals = seq(0,5,length.out=200)
mod_func = function(x, level){
	pred = predict(mod, data.frame(env1=x, env2=level), type='response', se.fit=T)
	vals = pred$fit + t(qnorm(c(0.025, 0.5, 0.975))%*%t(pred$se.fit))
	return(exp(vals))
}

mod = attach_mods['Bryophytes',j][[1]]
pred_high = mod_func(xvals, max(mod$model$env2))
pred_low = mod_func(xvals, min(mod$model$env2))

svg('./Figures/attachment_bryophyte_VPD.svg', height=4, width=4.5)
par(mar=c(4,4,1,1))
par(lend=1)

plot(yvals~Xdata[names(yvals), 'Bryophytes'], las=1, ylab='Attachment', xlab='Bryophyte Cover', pch=16, col='#00000050', axes=F)
polygon(c(xvals, rev(xvals)), c(pred_low[,1], rev(pred_low[,3])), col=use_col_trans[1], border=NA)
polygon(c(xvals, rev(xvals)), c(pred_high[,1], rev(pred_high[,3])), col=use_col_trans[2], border=NA)
lines(xvals, pred_low[,2], lwd=3, col=use_col[1])
lines(xvals, pred_high[,2], lwd=3, col=use_col[2])
axis(1, at=0:5, labels=c('none','minute','few','several','many', 'covered'))
axis(2, las=1)
box()
dev.off()

## Plot single effects of env vars on traits~quant
names(bark_tab)[2] = 'pred'
names(clim_tab)[2] = 'pred'
tab = rbind(bark_tab, clim_tab)
tab = tab[order(tab$trait),]
tab[order(tab$R2, decreasing=T),]

trait_order = c('ProbCyano','ProbCrustose','ProbFruticose','Attachment','LobeDissect','Cilia')
trait_names = c('Cyanolichen','Crustose','Fruticose','Attachment','Lobe dissection','Cilia')
names(trait_names) = trait_order
pred_order = c('WaterCapacity','Bryophytes','CloudFreq_mean','VPD_max','OpenPos')
pred_names = c('WHC','Bryo','Cloud','VPD','Open')
pred_names_long = c('Bark WHC','Bryophytes','Cloud Freq.','Max. VPD','Openness')
names(pred_names_long) = pred_order

tab$pred = factor(tab$pred, levels=pred_order)
tab$trait = factor(tab$trait, levels=rev(trait_order))

# Write out table
eff_df = tab[order(tab$trait, tab$R2, decreasing=T),]
eff_df$pred = pred_names_long[as.character(eff_df$pred)]
eff_df$CI = with(eff_df, paste('(',round(low95, 3),', ', round(up95, 3),')', sep=''))
eff_df[,c('P','R2','est')] = round(eff_df[,c('P','R2','est')], 3)
write.table(eff_df[,c('trait','pred','est','CI','P','R2')], './Figures/water traits vs env ests.txt', sep='\t', quote=F, row.names=F)

# Write out table with significant interactions
inter = acast(molten, trait~pred.plot+pred.samp~quant, subset=.(parm=='b1X2' & quant%in%c('est','low95','up95','P')))
inter_molten = melt(inter)
inter_tab = dcast(inter_molten, Var1 + Var2 ~ Var3)
inter_df = data.frame(trait=inter_tab$Var1)
inter_df$pred = sapply(inter_tab$Var2, function(x){
	getvars = names(which(sapply(pred_order, function(y) grep(y, x))==1))
	paste(pred_names_long[getvars], collapse=' x ')	
})
inter_df$est = round(inter_tab$est, 3)
inter_df$CI = with(inter_tab, paste('(',round(low95, 3),', ', round(up95, 3),')', sep=''))
inter_df$P = round(inter_tab$P, 3)
write.table(inter_df, './Figures/water traits vs env interactions.txt', sep='\t', quote=F, row.names=F)
write.table(subset(inter_df, P < 0.05), './Figures/water traits vs env sig interactions.txt', sep='\t', quote=F, row.names=F)



# Grouped by Trait
svg('./Figures/water traits vs env grp by trait.svg', height=6, width=2.5)
par(mfrow=c(5,1))
par(lend=1)
k=0
for(i in trait_order){
	par(mar=c(k*2/4, 4, (4-k)*2/4, 1))
	use_data = subset(tab, trait==i)
	xvals = as.numeric(use_data$pred)
	plot(est~as.numeric(pred), data=use_data, axes=F, pch=16, 
		ylim = range(use_data[,c('up95','low95')]), xlab='', ylab=trait_names[i])
	abline(h=0, lty=2,lwd=2)
	segments(xvals, use_data$low95, xvals, use_data$up95, lwd=2)
	axis(2, las=1)	
	box()
	if(k==0) axis(3, at=1:5, labels=pred_names, tick=F, line=-.5)
	if(k==4) axis(1, at=1:5, labels=pred_names, tick=F, line=-.5)
	k = k+1
}
dev.off()

# Grouped by env var



## Plot all relationships
mod_func = function(mod, x1, x2, logT = F){
	pred = predict(mod, data.frame(env1=x1, env2=x2), type='response', se.fit=T)
	fam = mod$family$family
	if(fam=='gaussian'){
		vals = pred$fit + t(qnorm(c(0.025, 0.5, 0.975))%*%t(pred$se.fit))
	}
	if(fam=='binomial'){
		vals = pred$fit%*%t(c(1,1,1))
	}
	if(logT){ 
		return(exp(vals))
	} else {
		return(vals)
	}
}

mod_list = list(crust_mods, frut_mods, attach_mods, dissect_mods, cilia_mods)
Ydata = list(crust_abun, frut_abun, sampXtrait_abun[,1], sampXtrait_abun[,2], cilia_abun)
pred_names = c('log(WHC)', 'Bryophyte Cover', 'Cloud Freq.', 'VPD', 'Openness')

svg('./Figures/water traits vs env.svg', height=6.5, width=6.5)
layout(matrix(1:25, byrow=T, nrow=5))
par(mar=rep(.3, 4))
par(oma=c(3.5,4.5,0,0))
par(lend=1)	
for(i in 1:5){
	trait = trait_order[i]
	use_mods = mod_list[[i]]
	yvals = Ydata[[i]]
	if(i %in% 3:4) yvals = exp(yvals)
	if(length(yvals)>nrow(sampXtrait_abun)) yvals = yvals[,1]/yvals[,2]
	
	for(j in 1:5){
		env = pred_order[j]
		xvals = Xdata[names(yvals),env] + ifelse(j %in% 3:5, climMeans[env], 0)

		xrange = range(Xdata[,env], na.rm=T) + ifelse(j %in% 3:5, climMeans[env], 0)
		xpredvals = seq(xrange[1], xrange[2], length.out=200)
	
		if(j %in% 1:2){
			mod = use_mods[env,1][[1]]
			ypred = mod_func(mod, xpredvals, mean(mod$model$env2), logT=i%in%3:4)
			P = all_ests[trait, env, 1, 'b1', 'P']
			r2 = all_ests[trait, env, 1, 'b1', 'R2']
		}
		if(j %in% 3:5){
			mod = use_mods[1,env][[1]]
			ypred = mod_func(mod, mean(mod$model$env1), xpredvals-climMeans[env], logT=i%in%3:4)
			P = all_ests[trait, 1, env, 'b2', 'P']
			r2 = all_ests[trait, 1, env, 'b2', 'R2'] 
		}
		

		#par(mar=c((i-1)*2/4, (5-j)*4/4, (5-i)*2/4, (j-1)*4/4))
		
		plot(yvals~xvals, las=1, axes=F, pch=16, col='#00000020', ylab='', xlab='')
		if(P < 0.05){
			polygon(c(xpredvals, rev(xpredvals)), c(ypred[,1], rev(ypred[,3])), col='#00000050', border=NA)
			lines(xpredvals, ypred[,2], lwd=1, col='black')
		}
		box()
		if(j==1){
			axis(2, las=1)
			mtext(trait_names[i], 2, 2.5, cex=.7)
		}
		if(i==5){
			axis(1)
			mtext(pred_names[j], 1, 2, cex=.7)
		}
	}
}
dev.off()


# Attachment: Bryophytes & All Clim
yvals = exp(sampXtrait_abun[,'Attachment'])
i = 'Bryophytes'

pdf('./Figures/attachment_bryophyte_interactions.pdf', height=8, width=6)
par(mfrow=c(3,2))
par(mar=c(4,4,1,1))
for(j in pred_vars2){
	mod = attach_mods[i,j][[1]]
	mod_func = function(x, level){	
		modlevel = level - climMeans[j]
		exp(predict(mod, data.frame(env1=x, env2=modlevel), type='response'))
	}
	
	climrange = range(plot_data[rownames(plot_data)!='Bladen1', j])
	plot(yvals~Xdata[names(yvals), i], las=1, ylab='Thallus Height', xlab=i, pch=16, col='#00000050')
	curve(mod_func(x, climrange[1]), from=0, to=5, add=T, lwd=2, col=use_col[1])
	curve(mod_func(x, climrange[2]), from=0, to=5, add=T, lwd=2, col=use_col[2])

	plot(yvals~I(Xdata[names(yvals), j]+climMeans[j]), las=1, ylab='Thallus Height', xlab=j, pch=16, col='#00000050')
	curve(mod_func(x=0, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[1])
	curve(mod_func(x=5, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[2])
}
dev.off()

# LobeDissect: Bryophytes & CloudFreq
yvals = exp(sampXtrait_abun[,'LobeDissect'])
i = 'Bryophytes'
j = 'CloudFreq_mean'

pdf('./Figures/lobedissect_bryophyte_interactions.pdf', height=4, width=6)
par(mar=c(4,4,1,1))
	mod = dissect_mods[i,j][[1]]
	mod_func = function(x, level){	
		modlevel = level - climMeans[j]
		exp(predict(mod, data.frame(env1=x, env2=modlevel), type='response'))
	}
	
	climrange = range(plot_data[rownames(plot_data)!='Bladen1', j])
	par(mfrow=c(1,2))
	plot(yvals~Xdata[names(yvals), i], las=1, ylab='Lobe Length:Width', xlab=i, pch=16, col='#00000050')
	curve(mod_func(x, climrange[1]), from=0, to=5, add=T, lwd=2, col=use_col[1])
	curve(mod_func(x, climrange[2]), from=0, to=5, add=T, lwd=2, col=use_col[2])

	plot(yvals~I(Xdata[names(yvals), j]+climMeans[j]), las=1, ylab='Lobe Length:Width', xlab=j, pch=16, col='#00000050')
	curve(mod_func(x=0, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[1])
	curve(mod_func(x=5, level=x), from=climrange[1], to=climrange[2], add=T, lwd=2, col=use_col[2])

dev.off()


## Within-sample trait diversity for attachment and growth form

# Attachment
attach_dmat = as.matrix(dist(0:6))
rownames(attach_dmat) = 0:6; colnames(attach_dmat) = 0:6
attach_abun = xtabs(AbunCount~SampID + Attachment, data=lichens_pm)
attach_abun = attach_abun/rowSums(attach_abun)

attach_rao = calc_rao(attach_abun, attach_dmat)
cbind(attach_abun, attach_rao) # rao is 0 when only one attachment height

# Growth Form
form_dmat = matrix(1, nrow=4, ncol=4)
diag(form_dmat) = 0
rownames(form_dmat) = levels(factor(lichens_pm$Form))
colnames(form_dmat) = levels(factor(lichens_pm$Form))
form_abun # from above
form_abun_std = form_abun/rowSums(form_abun)
form_rao = calc_rao(form_abun_std, form_dmat)

div_df = data.frame(SampID=rownames(form_abun), attach_rao, form_rao)
write.csv(div_df, 'trait_diversity.csv', row.names=F)

div_df = merge(div_df, samples_pm)

# Plot trait diversity versus N_tot, Bryophytes
plot(attach_rao~Bryophytes, data=div_df)
summary(lm(attach_rao~Bryophytes, data=div_df))

plot(attach_rao~N_tot, data=div_df)

plot(attach_rao~ sampXtrait_abun[div_df$SampID,'Attachment'], data=div_df)

plot(form_rao~Bryophytes, data=div_df)
summary(lm(form_rao~Bryophytes, data=div_df))

plot(form_rao~N_tot, data=div_df)


# Compare to null model.

################################################
#### Reprodutive mode ~ substrate stability ####
use_lichens = subset(lichens, !(is.na(Asco)|is.na(Asexual))) # drops one lichen from Bladen1-8-S4

## Binomial regression - P(reproductive structure | N species or thalli sampled)
## Samples with no lichens not in model

## Define a variable for reproductive mode
use_lichens$RepMode = factor(colSums(t(use_lichens[,c('Asco','Asexual')])*c(1,10)), levels=c(0,1,10,11))
levels(use_lichens$RepMode) = c('none','sexual','asexual','both')

## Calculate abundance and presence of each mode across samples
rep_pres = xtabs(~SampID+RepMode, data=use_lichens)
rep_abun = xtabs(AbunCount~SampID+RepMode, data=use_lichens)
tot_pres = rowSums(rep_pres)
tot_abun = rowSums(rep_abun)

## Calculate response variables
asco_pres = rowSums(rep_pres[,c('sexual','both')])
asco_abun = rowSums(rep_abun[,c('sexual','both')])
asex_pres = rowSums(rep_pres[,c('asexual','both')])
asex_abun = rowSums(rep_abun[,c('asexual','both')])
norep_pres = rep_pres[,'none']
norep_abun = rep_abun[,'none']

# Check for missing data
missing_lich = subset(samples, !(SampID %in% names(asco_pres))&NoLichen=='')$SampID
subset(use_lichens, SampID %in% missing_lich)

# Make a dataframe of responses
sum(rownames(rep_pres)!=rownames(rep_abun)) # Should be 0

reproduction = data.frame(SampID = rownames(rep_pres), Asco_pres=asco_pres, 
	Asco_abun=asco_abun, Asex_pres=asex_pres, Asex_abun=asex_abun, 
	Norep_pres=norep_pres, Norep_abun=norep_abun, S=tot_pres, N = tot_abun)

# Save data
#write.csv(reproduction, './Data/Derived Tables/reproduction.csv', row.names=F)
reproduction = read.csv('./Data/Derived Tables/reproduction.csv')
rownames(reproduction) = reproduction$SampID

# Subset to data just in top2 samples
reproduction = subset(reproduction, SampID %in% top2)

# Merge with predictors
reproduction = merge(reproduction, samples, by='SampID', all.x=T, all.y=F)

# Only do models without Bladen1 because Trans_tot_cor not (yet) calculated
use_data = subset(reproduction, PlotID!='Bladen1')

# Convert factors to numbers
use_data$Shedding = as.numeric(use_data$Shedding)
use_data$Bryophytes = as.numeric(use_data$Bryophytes)

# Define covariates to control for
bark_covars = c('FurrowDepth','Bryophytes','pH','Shedding')
other_covars = c('Trans_tot_cor','OpenPos','Elevation')
focal_var='DBH'

# Subset data to remove missing observations
use_data_nona = subset(use_data, !is.na(rowSums(use_data[,c('DBH','Shedding',bark_covars, other_covars)])))
cor(use_data_nona[,c('DBH',bark_covars)])

# Models (non-hierarchical)
#covars = c(bark_covars, other_covars)

# Calculate models
repmod_table = calc_repmods(use_data)
repmod_table2 = calc_repmods(reproduction)

# Save models
#DBH_mods = repmod_table
#Shed_mods = repmod_table
save(DBH_mods, Shed_mods, file='./Data/reproduction model tables.RData')
load('./Data/reproduction model tables.RData')

## Do models by hand to look at summary tables
x = cbind(use_data$Asco_pres, use_data$S-use_data$Asco_pres)
asco_pres_mod = glm(x ~ ., data=use_data[,c(other_covars, 'Shedding','DBH')], family=binomial(link='logit'))

x = cbind(use_data$Asex_pres, use_data$S-use_data$Asex_pres)
asex_pres_mod = glm(x ~ ., data=use_data[,c(other_covars, 'Shedding','DBH')], family=binomial(link='logit'))

x = cbind(use_data$Norep_pres, use_data$S-use_data$Norep_pres)
norep_pres_mod = glm(x ~ ., data=use_data[,c(other_covars, 'Shedding','DBH')], family=binomial(link='logit'))

x = cbind(use_data$Asco_abun, use_data$N-use_data$Asco_abun)
asco_abun_mod = glm(x ~ ., data=use_data[,c(other_covars, 'Shedding', focal_var)], family=binomial(link='logit'))

x = cbind(use_data$Asex_abun, use_data$N-use_data$Asex_abun)
asex_abun_mod = glm(x ~ ., data=use_data[,c(other_covars, 'Shedding', focal_var)], family=binomial(link='logit'))

x = cbind(use_data$Norep_abun, use_data$N-use_data$Norep_abun)
norep_abun_mod = glm(x ~ ., data=use_data[,c(other_covars, bark_covars, focal_var)], family=binomial(link='logit'))

summary(asco_pres_mod)
summary(asco_abun_mod)
summary(asex_pres_mod)
summary(asex_abun_mod)
summary(norep_pres_mod)
summary(norep_abun_mod)

cbind(coef(asco_abun_mod), coef(asex_abun_mod))


## Plot summaries of all model results
model_levels=c('asco_pres','asco_abun','asex_pres','asex_abun','norep_pres','norep_abun')[6:1]
use_mods = DBH_mods
use_ests = melt(use_mods[,,'est',])
use_ests$model = factor(paste(use_ests$response, use_ests$count, sep='_'), levels=model_levels)
use_CI = melt(use_mods[,,c('low95','up95'),])
use_CI$model = factor(paste(use_ests$response, use_ests$count, sep='_'), levels=model_levels)

mypch=c(1,16)
mycol='black'

pdf('./Figures/DBH models.pdf', height=4, width=4)
par(mar=c(3,4,1,1))
xyplot(model~value, groups=control, data=use_ests, panel=function(x,y,groups,...){
	panel.points(x, y, col=mycol, pch=mypch[groups])
	panel.arrows(use_CI[use_CI$term=='low95','value'], use_CI[use_CI$term=='low95','model'],
		use_CI[use_CI$term=='up95','value'], use_CI[use_CI$term=='low95','model'],
		code=3, length=0.05, angle=90) 
	panel.abline(v=0, col='grey', lty=2)
}, xlim=c(-.05, 0.05), xlab='DBH Effect', ylab='')
dev.off()


## Manuscript Figure: estimated effects of shedding/DBH on reproductive mode (abun) w/ and w/o covariates
rep_mods = abind(DBH_mods, Shed_mods, along=.5)
dimnames(rep_mods)[[1]] = c('DBH','Shedding')
names(dimnames(rep_mods)) = c('Predictor','Mode','Count','Statistic','Control')
rep_abun_df = dcast(melt(rep_mods[,c('asco','asex'),'abun',,]), Predictor+Mode+Control~Statistic)

yvals = rep(1:2,2) + rep(c(.1,-.1), each=2)

svg('./Figures/rep mode abun vs Shed and DBH.svg', height=2, width=6.5)
#layout(matrix(1:3, nrow=1), widths=c(.4,.4,.2))
par(mfrow=c(1,2))
use_col = c('grey50','white')
par(mar=c(4,0,1,1))
par(oma=c(0,9,0,0))

# Shedding
use_df = subset(rep_abun_df, Predictor=='Shedding')
plot(1,1, type='n', xlim=range(use_df[,c('up95','low95')]), ylim=c(0.5,2.5), axes=F, 
	xlab='Effect of Bark Instability', ylab='')
abline(v=0, lty=2, col='grey')
segments(use_df$low95, yvals, use_df$up95, yvals)
points(use_df$est, yvals, pch=21, bg=use_col[use_df$Mode])
axis(1)
axis(2, las=1, at=1:2, labels=c('With bark covariates','No bark covariates'), tick=F)
usr=par('usr')
segments(usr[1], usr[3], usr[2], usr[3])
mtext('A', 3, 0, adj=0)

# DBH
use_df = subset(rep_abun_df, Predictor=='DBH')
plot(1,1, type='n', xlim=range(use_df[,c('up95','low95')]), ylim=c(0.5,2.5), axes=F, 
	xlab='Effect of Tree DBH', ylab='')
abline(v=0, lty=2, col='grey')
segments(use_df$low95, yvals, use_df$up95, yvals)
points(use_df$est, yvals, pch=21, bg=use_col[use_df$Mode])
axis(1)
usr=par('usr')
segments(usr[1], usr[3], usr[2], usr[3])
mtext('B', 3, 0, adj=0)

# Legend
#legend('bottomright', c('Sexual','Asexual'), pch=21, pt.bg=use_col)
dev.off()


## Another summary table of model results used in SEEC presentation
use_mods = Shed_mods
use_ests = melt(use_mods[c('asco','asex'),,'est',])
use_ests$response = factor(use_ests$response, levels=c('asco','asex'))
levels(use_ests$response) = c('Sexual','Asexual')
levels(use_ests$control) = c('No bark covariates', 'With bark covariates')
use_ests$count = factor(use_ests$count, levels=c('abun','pres'))
levels(use_ests$count) = c('Prop. of Tot. Abun.', 'Prop. of Species')

use_CI = melt(use_mods[c('asco','asex'),,c('up95','low95'),])
use_CI$model = factor(use_CI$response, levels=c('asco','asex'))
levels(use_CI$response) = c('Sexual','Asexual')

mypch=c(21,21)
mycol=c('#ff6600','#5f8dd3')
linecol = c('black','black')

jit=c(-0.1,0.1)

pdf('./Figures/Shed models by count.pdf', height=2.5, width=6.5)
xyplot(control~value|count, groups=response, data=use_ests, panel=function(x,y,groups,subscripts,...){
	panel.abline(v=0, col='grey', lty=2)
	panel.arrows(use_CI[use_CI$term=='low95','value'][subscripts], as.numeric(y)+jit[groups[subscripts]],
		use_CI[use_CI$term=='up95','value'][subscripts], as.numeric(y)+jit[groups[subscripts]],
		code=3, length=0.04, angle=90, col=linecol[groups[subscripts]], lwd=1.2) 
	panel.points(x, as.numeric(y)+jit[groups[subscripts]], col=linecol[groups[subscripts]], fill=mycol[groups[subscripts]], pch=mypch[groups[subscripts]], cex=1)
	
}, strip=strip.custom(bg='transparent'), xlim=c(-0.32,0.32), xlab='', ylab='')#,
#scales=list(alternating=0))#, 
#key=list(space='right', title='Reproductive mode', cex.title=1.2, points=list(pch=mypch, col=mycol), text=list(levels(use_ests$response))))

dev.off()




## Make plots for single models
mycols = colorRampPalette(c('black','darkorange'))(5)
mypch = 18

# Pres-abs model
asco_probs = use_data$Asco_pres/use_data$S
asex_probs = use_data$Asex_pres/use_data$S

pdf('./Figures/Reproduction vs DBH and Shedding pres-abs model with Elev OpenPos.pdf', height=4, width=10)
layout(matrix(1:3, nrow=1), widths=c(.4,.4,.2))
par(mar=c(4,4,2,0))
plot(asco_probs~DBH, data=use_data, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='Proportion of Morphospecies', main='Sexual Reproduction')
for(i in 1:5){
	curve(predict(asco_pres_mod, data.frame(DBH=x, Shedding=i, Elevation=572, OpenPos=84),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot(asex_probs~DBH, data=use_data, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='', main='Asexual Reproduction')
for(i in 1:5){
	curve(predict(asex_pres_mod, data.frame(DBH=x, Shedding=i, Elevation=572, OpenPos=84),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot.new()
legend('center', levels(reproduction$Shedding), col=mycols, pch=mypch, lwd=2, bty='n', xjust=0)
dev.off()

use_elevs = c(150,500,1300)
elev_breaks=c(0,use_elevs, 1400)
mycols = colorRampPalette(c('black','darkorange'))(5)
pdf('./Figures/Reproduction vs Elevation and OpenPos model with DBH Shedding.pdf', height=4, width=10)
layout(matrix(1:3, nrow=1), widths=c(.4,.4,.2))
par(mar=c(4,4,2,0))
plot(asco_probs~OpenPos, data=use_data, col=mycols[as.numeric(cut(Elevation, breaks=elev_breaks))], pch=mypch, las=1, xlab='Topographic Openness', ylab='Proportion of Morphospecies', main='Sexual Reproduction')
for(i in 1:length(use_elevs)){
	curve(predict(asco_pres_mod, data.frame(DBH=20, Shedding=2, Elevation=use_elevs[i], OpenPos=x),type='response'),  add=T, lwd=2, col=mycols[i+1])
}
plot(asex_probs~OpenPos, data=use_data, col=mycols[as.numeric(cut(Elevation, breaks=elev_breaks))], pch=mypch, las=1, xlab='Topographic Openness', ylab='', main='Asexual Reproduction')
for(i in 1:length(use_elevs)){
	curve(predict(asex_pres_mod, data.frame(DBH=20, Shedding=2, Elevation=use_elevs[i], OpenPos=x),type='response'),  add=T, lwd=2, col=mycols[i+1])
}
plot.new()
legend('center', paste(use_elevs, 'm'), col=mycols[2:4], pch=mypch, lwd=2, bty='n', xjust=0)
dev.off()


## This plot was used when there was an interaction between shedding and DBH
pdf('./Figures/Reproduction vs DBH and Shedding pres-abs model.pdf', height=4, width=10)
layout(matrix(1:3, nrow=1), widths=c(.4,.4,.2))
par(mar=c(4,4,2,0))
plot(asco_probs~DBH, data=reproduction, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='Proportion of Morphospecies', main='Sexual Reproduction')
for(i in 1:5){
	curve(predict(asco_pres_mod, data.frame(DBH=x, Shedding=i),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot(asex_probs~DBH, data=reproduction, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='', main='Asexual Reproduction')
for(i in 1:5){
	curve(predict(asex_pres_mod, data.frame(DBH=x, Shedding=i),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot.new()
legend('center', levels(reproduction$Shedding), col=mycols, pch=mypch, lwd=2, bty='n', xjust=0)
dev.off()

# Abun model
asco_probs = reproduction$Asco_abun/reproduction$N
asex_probs = reproduction$Asex_abun/reproduction$N

pdf('./Figures/Reproduction vs DBH and Shedding abun model.pdf', height=4, width=10)
layout(matrix(1:3, nrow=1), widths=c(.4,.4,.2))
par(mar=c(4,4,2,0))
plot(asco_probs~DBH, data=reproduction, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='Proportion of Thalli', main='Sexual Reproduction')
for(i in 1:5){
	curve(predict(asco_abun_mod, data.frame(DBH=x, Shedding=i),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot(asex_probs~DBH, data=reproduction, col=mycols[as.numeric(Shedding)], pch=mypch, las=1, xlab='DBH', ylab='', main='Asexual Reproduction')
for(i in 1:5){
	curve(predict(asex_abun_mod, data.frame(DBH=x, Shedding=i),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot.new()
legend('center', levels(reproduction$Shedding), col=mycols, pch=mypch, lwd=2, bty='n', xjust=0)
dev.off()

## Plot for SEEC presentation

dbh_est_df = melt(DBH_mods[c('asco','asex'),'abun','est',])
dbh_CI_df = melt(DBH_mods[c('asco','asex'),'abun',c('up95','low95'),])

shed_est_df = melt(DBH_mods[c('asco','asex'),'abun','est',])
she_CI_df = melt(DBH_mods[c('asco','asex'),'abun',c('up95','low95'),])


## Plot Prop reproductive mode vs shedding
# use_data should be reproduction w/o Bladen
boxplot(I(Asex_abun/N)~Shedding, data=use_data)
boxplot(I(Asco_abun/N)~Shedding, data=use_data)

# Another predictor: Model with light
asex_abun_mod = glm(x ~ DBH*as.numeric(Shedding)+Trans_tot_cor, data=asex_abun, family=binomial(link='logit'))
asco_abun_mod = glm(x ~ DBH*as.numeric(Shedding)+Trans_tot_cor, data=asco_abun, family=binomial(link='logit'))

par(mfrow=c(1,2))
plot(asex_probs~asex_abun$DBH, col=mycols[as.numeric(asex_pres$Shedding)], pch=mypch, las=1, xlab='DBH', ylab='', main='Asexual Reproduction')
for(i in 1:5){
	curve(predict(asex_abun_mod, data.frame(DBH=x, Shedding=i, Trans_tot_cor = mean(asex_abun$Trans_tot_cor, na.rm=T)),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot(asex_probs~asex_abun$Trans_tot_cor, col=mycols[as.numeric(asex_pres$Shedding)], pch=mypch, las=1, xlab='Light Availability', ylab='', main='Asexual Reproduction')
for(i in 1:5){
	curve(predict(asex_abun_mod, data.frame(DBH=mean(asex_abun$DBH, na.rm=T), Shedding=i, Trans_tot_cor = x),type='response'),  add=T, lwd=2, col=mycols[i])
}
par(mfrow=c(1,2))
plot(asco_probs~asco_abun$DBH, col=mycols[as.numeric(asco_pres$Shedding)], pch=mypch, las=1, xlab='DBH', ylab='', main='Sexual Reproduction')
for(i in 1:5){
	curve(predict(asco_abun_mod, data.frame(DBH=x, Shedding=i, Trans_tot_cor = mean(asco_abun$Trans_tot_cor, na.rm=T)),type='response'),  add=T, lwd=2, col=mycols[i])
}
plot(asco_probs~asex_abun$Trans_tot_cor, col=mycols[as.numeric(asco_pres$Shedding)], pch=mypch, las=1, xlab='Light Availability', ylab='', main='Sexual Reproduction')
for(i in 1:5){
	curve(predict(asco_abun_mod, data.frame(DBH=mean(asco_abun$DBH, na.rm=T), Shedding=i, Trans_tot_cor = x),type='response'),  add=T, lwd=2, col=mycols[i])
}

### Hierarchical models for assessing variance components
use_data = subset(reproduction, PlotID!='Bladen1')

x = cbind(use_data$Asco_pres, use_data$S-use_data$Asco_pres)
asco_pres_mod0 = glmer(x ~ (1|SiteID/PlotID/TreeID), data=use_data, family=binomial(link='logit'))

x = cbind(use_data$Asco_abun, use_data$N-use_data$Asco_abun)
asco_abun_mod0 = glmer(x ~ (1|SiteID/PlotID/TreeID), data=use_data, family=binomial(link='logit'))

x = cbind(use_data$Asex_pres, use_data$S-use_data$Asex_pres)
asex_pres_mod0 = glmer(x ~ (1|SiteID/PlotID/TreeID), data=use_data, family=binomial(link='logit'))

x = cbind(use_data$Asex_abun, use_data$N-use_data$Asex_abun)
asex_abun_mod0 = glmer(x ~ (1|SiteID/PlotID/TreeID), data=use_data, family=binomial(link='logit'))

summary(asco_pres_mod0)
summary(asex_pres_mod0)
summary(asco_abun_mod0)
summary(asex_abun_mod0)

# Plot variation among scales
mycol = colorRampPalette(c('blue','darkorange'))(9)
use_cols = c(mycol[1], rep(mycol[2:6],each=2), rep(mycol[7],4), rep(mycol[8:9], each=2))

asco_probs = reproduction$Asco_pres/reproduction$S
asex_probs = reproduction$Asex_pres/reproduction$S

pdf('./Figures/Reproductive variance across plots pres-abs.pdf', height=7, width=7)
par(mfrow=c(2,1))
par(mar=c(1,4,3,1))
boxplot(asco_probs~PlotID, data=reproduction, las=2, xaxt='n', ylab='Proportion Sexual', col=use_cols)
par(mar=c(4,4,0,1))
boxplot(asex_probs~PlotID, data=reproduction, las=2, ylab='Proportion Asexual',col=use_cols)
dev.off()

asco_probs = reproduction$Asco_abun/reproduction$N
asex_probs = reproduction$Asex_abun/reproduction$N

pdf('./Figures/Reproductive variance across plots abun.pdf', height=7, width=7)
par(mfrow=c(2,1))
par(mar=c(1,4,3,1))
boxplot(asco_probs~PlotID, data=reproduction, las=2, xaxt='n', ylab='Proportion Sexual', col=use_cols)
par(mar=c(4,4,0,1))
boxplot(asex_probs~PlotID, data=reproduction, las=2, ylab='Proportion Asexual',col=use_cols)
dev.off()

##############################################
#### Competition models: form ~ abundance ####
reproduction = read.csv('./Data/Derived Tables/reproduction.csv')
rownames(reproduction) = reproduction$SampID
reproduction_pm = subset(reproduction, SampID %in% samples[samples$SiteID!='Bladen','SampID'])


# Focal traits: crustose, attachment, asexual
# Focal predictors: total abundance and bryophyte abundance

# Create sample-level response variables for P(Crustose), Attachment, and P(Asexual) in above sections
crust_abun
sampXtrait_abun[,'Attachment']
asex_abun = reproduction_pm[,c('Asex_abun','N')]

# Make table of predictors
scale_vars = c('SampID','PlotID','SiteID','Ecoregion')
Xdata = merge(reproduction_pm[,c('SampID','N','S')], samples[,c('DBH','Bryophytes',scale_vars,'TreeTaxonID')], all.x=T)
rownames(Xdata) = Xdata$SampID

# Convert Bryophytes to integer
Xdata$Bryophytes = unclass(Xdata$Bryophytes) - 1

# Calculate derived abundance estimates
Xdata$N_cor = Xdata$N - Xdata$S
Xdata$N_avg = Xdata$N / Xdata$S

## Plot relationships across plots

Y = as.matrix(crust_abun) # asex_abun
yvals = Y[,1]/Y[,2]

Xdata = Xdata[rownames(Y),]
Xdata = droplevels(Xdata)

# rescale N so that models fit more easily
Xdata$N = Xdata$N/9

# Examine separate means models
plot(yvals~N, data=Xdata, pch=16, col='#00000050')
for(i in unique(lichens_pm$PlotID)){
	use_obs = Xdata$PlotID==i
	mod = glm(Y~N, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(N=x), type='response')

	curve(mod_func(x), from=min(Xdata$N[use_obs]), to = max(Xdata$N[use_obs]), add=T, lwd=2)
}
for(i in unique(lichens_pm$Ecoregion)){
	use_obs = Xdata$Ecoregion==i
	mod = glm(Y~N, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(N=x), type='response')

	curve(mod_func(x), from=min(Xdata$N[use_obs]), to = max(Xdata$N[use_obs]), add=T, lwd=2, col=2)
}

plot(yvals~N_avg, data=Xdata, pch=16, col='#00000050')
for(i in unique(lichens_pm$PlotID)){
	use_obs = Xdata$PlotID==i
	mod = glm(crust_abun~N_avg, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(N_avg=x), type='response')

	curve(mod_func(x), from=min(Xdata$N_avg[use_obs]), to = max(Xdata$N_avg[use_obs]), add=T, lwd=2)
}
for(i in unique(lichens_pm$Ecoregion)){
	use_obs = Xdata$Ecoregion==i
	mod = glm(crust_abun~N_avg, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(N_avg=x), type='response')

	curve(mod_func(x), from=min(Xdata$N_avg[use_obs]), to = max(Xdata$N_avg[use_obs]), add=T, lwd=2, col=2)
}

plot(yvals~Bryophytes, data=Xdata, pch=16, col='#00000050')
for(i in unique(lichens_pm$PlotID)){
	use_obs = Xdata$PlotID==i
	mod = glm(Y~Bryophytes, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(Bryophytes=x), type='response')

	curve(mod_func(x), from=0, to = 5, add=T, lwd=2)
}
for(i in unique(lichens_pm$Ecoregion)){
	use_obs = Xdata$Ecoregion==i
	mod = glm(Y~Bryophytes, data=Xdata, family='binomial', subset=use_obs)
	mod_func = function(x) predict(mod, data.frame(Bryophytes=x), type='response')

	curve(mod_func(x), from=0, to = 5, add=T, lwd=2, col=2)
}


## Models with random effect of plot
mod1 = glmer(Y ~ N + (N|PlotID), data=Xdata, family='binomial')
mod2 = glmer(Y ~ N + (1|PlotID), data=Xdata, family='binomial')
anova(mod1, mod2) # decide whether random effect of plot indicated:

mod1 = glmer(Y ~ Bryophytes + (Bryophytes|PlotID), data=Xdata, family='binomial')
mod2 = glmer(Y ~ Bryophytes + (1|PlotID), data=Xdata, family='binomial')
anova(mod1, mod2)

# Define models
crust_mod = glmer(Y ~ N + (1|PlotID), data=Xdata, family='binomial')
asex_mod = glmer(Y ~ N + (N|PlotID), data=Xdata, family='binomial')

crust_mod_bryo = glmer(Y ~ Bryophytes + (1|PlotID), data=Xdata, family='binomial')

# Generate null model
mod_null = glmer(Y ~ 1 + (1|PlotID), data=Xdata, family='binomial')
anova(crust_mod_bryo, mod_null)
anova(crust_mod, mod_null)

## Models with random effect of tree species
# Drop samples on tree species with < 10 samples observed across the data set
treefreq = table(Xdata$TreeTaxonID)
keep_trees = names(treefreq[treefreq>=10])
keep_obs = Xdata$TreeTaxonID %in% keep_trees

mod_tree1 = glmer(Y ~ N_avg + (N_avg|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)
mod_tree2 = glmer(Y ~ N_avg + (1|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)
anova(mod_tree1, mod_tree2)

crust_mod_tree = glmer(Y ~ N_avg + (1|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)
asex_mod_tree = glmer(Y ~ N_avg + (N_avg|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)

# Generate null model
mod_null_tree = glmer(Y ~ 1+ (1|TreeTaxonID), data=Xdata, family='binomial', subset=keep_obs)

# Plot effect
yvals = Y[,1]/Y[,2]

pdf('./Figures/asexual vs avg N.pdf', height=4, width=7)
par(mfrow=c(1,2))
plot(yvals~N_avg, data=Xdata, pch=16, col='#00000050', las=1, ylab='Asexual Prob.', 
	xlab='Avg. Abundance (# Squares / Species)', xlim=c(1,9), main='Tree Species Random Effects')

for(i in keep_trees){
	Nrange = range(subset(Xdata, TreeTaxonID==i)$N_avg)
	mod_func = function(x) predict(asex_mod_tree, data.frame(N_avg=x, TreeTaxonID=i), type='response')
	curve(mod_func(x), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}

mod_func = function(x) predict(asex_mod_tree, data.frame(N_avg=x), type='response', re.form=NA)
curve(mod_func(x), add=T, lwd=3)

r2 = calc_fixedvar(asex_mod_tree) / calc_totvar(asex_mod_tree)
P = anova(asex_mod_tree, mod_null_tree, test='Chisq')[2,8]
mylabel1 = bquote(italic(R)^2 == .(format(r2, digits = 1)))
mylabel2 = bquote(P == .(format(P, digits = 1)))
mtext(mylabel1, 3, 0, adj=1)
mtext(mylabel2, 3, 0, adj=0)

plot(yvals~N_avg, data=Xdata, pch=16, col='#00000050', las=1, ylab='Asexual Prob.', 
	xlab='Avg. Abundance (# Squares / Species)', xlim=c(1,9), main='Plot Random Effects')

for(i in unique(Xdata$PlotID)){
	Nrange = range(subset(Xdata, PlotID==i)$N_avg)
	mod_func = function(x) predict(asex_mod, data.frame(N_avg=x, PlotID=i), type='response')
	curve(mod_func(x), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}

mod_func = function(x) predict(asex_mod, data.frame(N_avg=x), type='response', re.form=NA)
curve(mod_func(x), add=T, lwd=3)

r2 = calc_fixedvar(asex_mod) / calc_totvar(asex_mod)
P = anova(asex_mod, mod_null, test='Chisq')[2,8]
mylabel1 = bquote(italic(R)^2 == .(format(r2, digits = 1)))
mylabel2 = bquote(P == .(format(P, digits = 1)))
mtext(mylabel1, 3, 0, adj=1)
mtext(mylabel2, 3, 0, adj=0)
dev.off()

# Change with tree size?



### Attachment

Y = sampXtrait_abun[,'Attachment'] # Make sure this is log-transformed from above
Xdata = Xdata[names(Y),]

## Models with random effect of plot
mod1 = lmer(Y ~ N + (N|PlotID), data=Xdata)
mod2 = lmer(Y ~ N + (1|PlotID), data=Xdata)
anova(mod1, mod2)

mod1 = lmer(Y ~ Bryophytes + (Bryophytes|PlotID), data=Xdata)
mod2 = lmer(Y ~ Bryophytes + (1|PlotID), data=Xdata)
anova(mod1, mod2)

attach_mod = lmer(Y ~ N + (N|PlotID), data=Xdata)
attach_mod_bryo = lmer(Y ~ Bryophytes + (Bryophytes|PlotID), data=Xdata)

# Generate null model
mod_null = lmer(Y ~ 1 + (1|PlotID), data=Xdata)

anova(attach_mod, mod_null)
anova(attach_mod_bryo, mod_null)

## Models with random effect of tree species
# Drop samples on tree species with < 10 samples observed across the data set
treefreq = table(Xdata$TreeTaxonID)
keep_trees = names(treefreq[treefreq>=10])
keep_obs = Xdata$TreeTaxonID %in% keep_trees

mod_tree1 = lmer(Y ~ N_avg + (N_avg|TreeTaxonID), data=Xdata, subset=keep_obs)
mod_tree2 = lmer(Y ~ N_avg + (1|TreeTaxonID), data=Xdata, subset=keep_obs)
anova(mod_tree1, mod_tree2)

attach_mod_tree = lmer(Y ~ N_avg + (N_avg|TreeTaxonID), data=Xdata, subset=keep_obs)

# Generate null model
mod_null_tree = lmer(Y ~ 1+ (1|TreeTaxonID), data=Xdata, subset=keep_obs)

# Plot effect
yvals = exp(Y)

pdf('./Figures/attachment vs avg N.pdf', height=4, width=7)
par(mfrow=c(1,2))
plot(yvals~N_avg, data=Xdata, pch=16, col='#00000050', las=1, ylab='', 
	xlab='Avg. Abundance (# Squares / Species)', xlim=c(1,9), main='Tree Species Random Effects')

for(i in keep_trees){
	Nrange = range(subset(Xdata, TreeTaxonID==i)$N_avg)
	mod_func = function(x) exp(predict(attach_mod_tree, data.frame(N_avg=x, TreeTaxonID=i), type='response'))
	curve(mod_func(x), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}

mod_func = function(x) exp(predict(attach_mod_tree, data.frame(N_avg=x), type='response', re.form=NA))
curve(mod_func(x), add=T, lwd=3)

r2 = calc_fixedvar(attach_mod_tree) / calc_totvar(attach_mod_tree)
P = anova(attach_mod_tree, mod_null_tree, test='Chisq')[2,8]
mylabel1 = bquote(italic(R)^2 == .(format(r2, digits = 1)))
mylabel2 = bquote(P == .(format(P, digits = 1)))
mtext(mylabel1, 3, 0, adj=1)
mtext(mylabel2, 3, 0, adj=0)

plot(yvals~N_avg, data=Xdata, pch=16, col='#00000050', las=1, ylab='Asexual Prob.', 
	xlab='Avg. Abundance (# Squares / Species)', xlim=c(1,9), main='Plot Random Effects')

for(i in unique(Xdata$PlotID)){
	Nrange = range(subset(Xdata, PlotID==i)$N_avg)
	mod_func = function(x) exp(predict(attach_mod, data.frame(N_avg=x, PlotID=i), type='response'))
	curve(mod_func(x), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}

mod_func = function(x) exp(predict(attach_mod, data.frame(N_avg=x), type='response', re.form=NA))
curve(mod_func(x), add=T, lwd=3)

r2 = calc_fixedvar(attach_mod) / calc_totvar(attach_mod)
P = anova(attach_mod, mod_null, test='Chisq')[2,8]
mylabel1 = bquote(italic(R)^2 == .(format(r2, digits = 1)))
mylabel2 = bquote(P == .(format(P, digits = 1)))
mtext(mylabel1, 3, 0, adj=1)
mtext(mylabel2, 3, 0, adj=0)
dev.off()



## Manuscript Figure: predicted attachment height and crutose lichens with N and Bryophytes
mod_func = function(mod, x, PID=NA, logT=F){
	
	if(is.na(PID)){
		use_df = data.frame(Pred=x)
		names(use_df) = names(mod@frame)[2]
		ypred = predict(mod, use_df, type='response', re.form=NA)
	} else {
		use_df = data.frame(Pred=x, PlotID=PID)
		names(use_df)[1] = names(mod@frame)[2]
		ypred = predict(mod, use_df, type='response')
	}
	if(logT){
		ypred = exp(ypred)
	} else {
		ypred = ypred*1
	}
}

# Remake data tables
Xdata = merge(reproduction_pm[,c('SampID','N','S')], samples[,c('DBH', 'Bryophytes', scale_vars,'TreeTaxonID')], all.x=T)
Xdata$N = Xdata$N/9
rownames(Xdata) = Xdata$SampID
Xdata$Bryophytes = unclass(Xdata$Bryophytes) - 1

# Color by trait diversity z-scores:
load('trait_diversity.RData')

pdf('./Figures/total abun and bryo vs attachment and crustose color by z-score.pdf', height=5, width=7.5)
layout(matrix(1:6, nrow=2, byrow=T), widths=c(0.45, 0.45, 0.1))
par(mar=c(.2, 2.5, 1.5, .5))
par(oma=c(3.5,3,0,1))

# Attachment ~ N
yvals = exp(sampXtrait_abun[,'Attachment'])

#range(trait_div[,c('Attachment'),'Z'], na.rm=T)
use_breaks = c(-4.1, seq(-4, 4, .5), 4.1)
use_col = colorRampPalette(blue2red)(length(use_breaks)-1)
use_col = paste(use_col, 'A0', sep='')
colvals = use_col[cut(trait_div[names(yvals),'Attachment','Z'],use_breaks)]
colvals[is.na(colvals)] = '#000000A0'

xvals = Xdata[names(yvals),'N']
plot(xvals, yvals, pch=16, col=colvals, las=1, ylab='', 
	xlab='', axes=F)
axis(2, las=1)
box()

for(i in unique(Xdata$PlotID)){
	Nrange = range(subset(Xdata, PlotID==i)$N)
	curve(mod_func(attach_mod, x, i, logT=T), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}
curve(mod_func(attach_mod, x, logT=T), add=T, lwd=3)
mtext('A', 3, 0, adj=0)
mtext('Attachment', 2, 3)

# Attachment ~ Bryo
xvals = Xdata[names(yvals),'Bryophytes']
plot(xvals, yvals, pch=16, col=colvals, las=1, ylab='', 
	xlab='', axes=F)
axis(2, las=1)
box()

for(i in unique(Xdata$PlotID)){
	Xrange = range(subset(Xdata, PlotID==i)$Bryophytes)
	curve(mod_func(attach_mod_bryo, x, i, logT=T), from=Xrange[1], to=Xrange[2], add=T, lwd=1, col='grey50')
}
curve(mod_func(attach_mod_bryo, x, logT=T), add=T, lwd=3)
mtext('B', 3, 0, adj=0)

# Legend
par(mar=c(.2, .5, 1.5, 2.5))
z=matrix(1:(length(use_breaks)-1),nrow=1)
x=1
y=use_breaks 
image(x,y,z,col=use_col,axes=FALSE,xlab="",ylab="")
axis(4, las=1)

# Prop. crustose ~ N
par(mar=c(.2, 2.5, 1.5, .5))
yvals = crust_abun[,1] / crust_abun[,2]

#range(trait_div[,c('Form'),'Z'], na.rm=T)
use_breaks = seq(-3, 3, .5)
use_col = colorRampPalette(blue2red)(length(use_breaks)-1)
use_col = paste(use_col, 'A0', sep='')
colvals = use_col[cut(trait_div[names(yvals),'Form','Z'],use_breaks)]
colvals[is.na(colvals)] = '#000000A0'

xvals = Xdata[names(yvals),'N']
plot(xvals, yvals, pch=16, col=colvals, las=1, ylab='', 
	xlab='', axes=F)
axis(2, las=1)
box()

for(i in unique(Xdata$PlotID)){
	Nrange = range(subset(Xdata, PlotID==i)$N)
	curve(mod_func(crust_mod, x, i), from=Nrange[1], to=Nrange[2], add=T, lwd=1, col='grey50')
}
curve(mod_func(crust_mod, x), add=T, lwd=3)
mtext('C', 3, 0, adj=0)
mtext('Prop. Crustose', 2, 3)

axis(1, at=0:5, labels=9*(0:5))
mtext(expression(N[TOT]), 1, 2.5)

# Prop. crustose ~ Bryophytes
xvals = Xdata[names(yvals),'Bryophytes']
plot(xvals, yvals, pch=16, col=colvals, las=1, ylab='', 
	xlab='', axes=F)
axis(2, las=1)
box()

for(i in unique(Xdata$PlotID)){
	Xrange = range(subset(Xdata, PlotID==i)$Bryophytes)
	curve(mod_func(crust_mod_bryo, x, i), from=Xrange[1], to=Xrange[2], add=T, lwd=1, col='grey50')
}
curve(mod_func(crust_mod_bryo, x), add=T, lwd=3)
mtext('D', 3, 0, adj=0)

axis(1, at=c(0:5), labels=levels(Xdata$Bryophytes))
mtext('Bryophyte Cover', 1, 2.5)

# Legend
par(mar=c(.2, .5, 1.5, 2.5))
z=matrix(1:(length(use_breaks)-1),nrow=1)
x=1
y=use_breaks
image(x,y,z,col=use_col,axes=FALSE,xlab="",ylab="")
axis(4, las=1)

dev.off()


### Succession?

Xdata$DBH = Xdata$DBH/10

# PropCrustose
Y = as.matrix(crust_abun) # asex_abun
yvals = Y[,1]/Y[,2]
mod_dbh = glmer(Y ~ DBH + (1|PlotID), data=Xdata[rownames(Y),], family='binomial')
mod_crust = mod_dbh
plot(yvals ~ I(Xdata$DBH*10))

# PropFoliose - decreases in most plots
Y = cbind(form_abun[,'foliose'], rowSums(form_abun))
yvals = Y[,1]/Y[,2]
mod_dbh1 = glmer(Y ~ DBH + (1|PlotID), data=Xdata[rownames(Y),], family='binomial')
mod_dbh2 = glmer(Y ~ DBH + (DBH|PlotID), data=Xdata[rownames(Y),], family='binomial')
anova(mod_dbh1, mod_dbh2)
effs = exp(coef(mod_dbh2)$PlotID)$DBH - 1
effs[order(effs)]
coef(mod_dbh2)$PlotID[order(effs),]
mod_fol = mod_dbh2

# PropFruticose - increases by 19% every 10cm
Y = cbind(form_abun[,'fruticose'], rowSums(form_abun))
yvals = Y[,1]/Y[,2]
mod_dbh1 = glmer(Y ~ DBH + (1|PlotID), data=Xdata[rownames(Y),], family='binomial')
mod_dbh2 = glmer(Y ~ DBH + (DBH|PlotID), data=Xdata[rownames(Y),], family='binomial')
anova(mod_dbh1, mod_dbh2)
summary(mod_dbh1)
exp(coef(mod_dbh1)$PlotID)
mod_frut = mod_dbh1

# Attachment - no change
Y = sampXtrait_abun[,'Attachment'] # Make sure this is log-transformed from above
mod_dbh = lmer(Y ~ DBH + (1|PlotID), data=Xdata[names(Y),])

# Bryophytes - usually increase
Y = cbind(Xdata$Bryophytes, 5)
mod_dbh1 = glmer(Y ~ DBH + (1|PlotID), data=Xdata, family='binomial')
mod_dbh2 = glmer(Y ~ DBH + (DBH|PlotID), data=Xdata, family='binomial')
anova(mod_dbh1, mod_dbh2)
mod_bryo = mod_dbh2
effs = exp(coef(mod_dbh2)$PlotID)$DBH - 1
effs[order(effs)]
coef(mod_dbh2)$PlotID[order(effs),]


## Plot how growth forms and Bryophyte abundance predicted to change with DBH
use_col = c('slateblue','firebrick1','darkorange','grey30')
dbh_vals = sapply(unique(samples_pm$PlotID), function(x){
	use_data = subset(Xdata, PlotID==x)
	range(use_data$DBH)
})
dbh_vals = t(dbh_vals); rownames(dbh_vals)=unique(samples_pm$PlotID)


pdf('./Figures/growth form and bryos vs DBH.pdf', height=6, width=4)
par(mfrow=c(2,1))
par(mar=c(.5,4,1.5,.5))
par(oma=c(4,0,0,0))

plot(0,0, xlim=c(0, max(samples_pm$DBH)), ylim=c(0,1), type='n', xlab='', ylab='Proportion',axes=F)
axis(2, las=1)
axis(1, at=seq(0, 80, 10), labels=rep('', 9))
box()

# Crustose
for(i in unique(samples_pm$PlotID)){
	xvals = seq(dbh_vals[i,1], dbh_vals[i,2], length.out=100)
	yvals = predict(mod_crust, data.frame(DBH=xvals, PlotID=i), type='response')
	lines(xvals*10, yvals, col=use_col[1], lwd=1)
}

# Foliose
for(i in unique(samples_pm$PlotID)){
	xvals = seq(dbh_vals[i,1], dbh_vals[i,2], length.out=100)
	yvals = predict(mod_fol, data.frame(DBH=xvals, PlotID=i), type='response')
	lines(xvals*10, yvals, col=use_col[2], lwd=1)
}

# Fruticose
for(i in unique(samples_pm$PlotID)){
	xvals = seq(dbh_vals[i,1], dbh_vals[i,2], length.out=100)
	yvals = predict(mod_frut, data.frame(DBH=xvals, PlotID=i), type='response')
	lines(xvals*10, yvals, col=use_col[3], lwd=1)
}
legend('topright', c('Crustose','Foliose','Fruticose'), col=use_col, lty=1,bty='n')
mtext('Proportion', 2, 3)
mtext('A. Growth form', 3,0, adj=0)

# Bryophytes
plot(0,0, xlim=c(0, max(samples_pm$DBH)), ylim=c(0,1), type='n', xlab='Tree DBH (cm)', ylab='', axes=F)
axis(1, at=seq(0, 80, 10))
axis(2, at=seq(0,1, length.out=6), labels=c('none','minute','few','several','many','covered'), las=1)
box()
for(i in unique(samples_pm$PlotID)){
	xvals = seq(dbh_vals[i,1], dbh_vals[i,2], length.out=100)
	yvals = predict(mod_bryo, data.frame(DBH=xvals, PlotID=i), type='response')
	lines(xvals*10, yvals, col=use_col[4], lwd=1)
}

mtext('B. Bryophyte cover', 3, 0, adj=0)
mtext('Tree DBH (cm)', 1, 2.5)

dev.off()

plot(Y, rank(Y))
plot(yvals, rank(yvals))
################################################################################
### Trait Diversity ###

library(abind)

## Calculate observed trait diversity within each sample for:
# Attachment, LobeDissect, Reproductive mode, Growth Form
groupby = lichens_pm$SampID
Nby = lichens_pm$AbunCount

attach_div = calc_trait_div(lichens_pm$Attachment + 1, Nby, groupby,  'numeric')
dissect_div = calc_trait_div(lichens_pm$LobeDissect, Nby, groupby, 'numeric')
asco_div = calc_trait_div(lichens_pm$Asco, Nby, groupby, 'symm')
asex_div = calc_trait_div(lichens_pm$Asexual, Nby, groupby, 'symm')
form_div = calc_trait_div(lichens_pm$Form, Nby, groupby, 'factor')

## Calculate expected trait diversity null distributions
reps = 1000
attach_nulls = calc_trait_div_null(reps, lichens_pm$Attachment + 1, Nby, groupby,  'numeric')
dissect_nulls = calc_trait_div_null(reps, lichens_pm$LobeDissect, Nby, groupby, 'numeric')
asco_nulls = calc_trait_div_null(reps, lichens_pm$Asco, Nby, groupby, 'symm')
asex_nulls = calc_trait_div_null(reps, lichens_pm$Asexual, Nby, groupby, 'symm')
form_nulls = calc_trait_div_null(reps, lichens_pm$Form, Nby, groupby, 'factor')

## Calculate expected null crustose and fruticose proportions and reproductive mode proportions

form_prop = form_abun / rowSums(form_abun)
form_prop = form_prop[names(form_div),]
form_prop_nulls = sapply(1:reps, function(i){
	neworder = sample(lichens_pm$Form, nrow(lichens_pm), replace=F)
	newtab = xtabs(Nby ~ groupby+neworder)
	newtab = newtab/rowSums(newtab)
	newtab
}, simplify='array')

rep_prop = rep_abun / rowSums(rep_abun)
rep_prop = rep_prop[names(asco_div),]
rep_prop_nulls = sapply(1:reps, function(i){
	neworder = sample(nrow(lichens_pm), replace=F)
	newrep = lichens_pm[neworder, c('Asco','Asexual')]
	newasco = xtabs(Nby ~ groupby+newrep[,'Asco'])
	newasco = newasco[,'TRUE']/rowSums(newasco)
	newasex = xtabs(Nby ~ groupby+newrep[,'Asexual'])
	newasex = newasex[,'TRUE']/rowSums(newasex)
	newtab = cbind(newasco, newasex)
	colnames(newtab) = c('Asco','Asexual')
	newtab
}, simplify='array')

save(attach_nulls, dissect_nulls, asco_nulls, asex_nulls, form_nulls, form_prop_nulls, rep_prop_nulls, file='trait_diversity_nulls.RData')

## Calculate P-values and z-scores of observed trait diversity based on nulls

# Assumes nulls and obs are all in same order (by SampID)
trait_div = array(NA, dim=c(length(attach_div), 9, 3), 
	dimnames=list(SampID=names(attach_div), Trait=c('Attachment','LobeDissect','Asco','Asexual','Form','PropCrust','PropFrut','PropAsco','PropAsex'),
		Statistic=c('Obs','Z','P')))
calc_pz = function(val, null) c(calc_p(val, null), calc_z(val, null))

for(i in dimnames(trait_div)[[1]]){
	trait_div[i,'Attachment',c('P','Z')] = calc_pz(attach_div[i], attach_nulls[i,])
	trait_div[i,'LobeDissect',c('P','Z')] = calc_pz(dissect_div[i], dissect_nulls[i,])
	trait_div[i,'Asco',c('P','Z')] = calc_pz(asco_div[i], asco_nulls[i,])
	trait_div[i,'Asexual',c('P','Z')] = calc_pz(asex_div[i], asex_nulls[i,])
	trait_div[i,'Form',c('P','Z')] = calc_pz(form_div[i], form_nulls[i,])	
	trait_div[i,'PropCrust',c('P','Z')] = calc_pz(form_prop[i,'crustose'], form_prop_nulls[i,'crustose',])
	trait_div[i,'PropFrut',c('P','Z')] = calc_pz(form_prop[i,'fruticose'], form_prop_nulls[i,'fruticose',])
	trait_div[i,'PropAsco',c('P','Z')] = calc_pz(rep_prop[i,'sexual'], rep_prop_nulls[i,'Asco',])
	trait_div[i,'PropAsex',c('P','Z')] = calc_pz(rep_prop[i,'asexual'], rep_prop_nulls[i,'Asexual',])
}

trait_div[names(attach_div),'Attachment','Obs'] = attach_div
trait_div[names(dissect_div),'LobeDissect','Obs'] = dissect_div
trait_div[names(asco_div),'Asco','Obs'] = asco_div
trait_div[names(asex_div),'Asexual','Obs'] = asex_div
trait_div[names(form_div),'Form','Obs'] = form_div
trait_div[,'PropCrust','Obs'] = form_prop[dimnames(trait_div)[[1]],'crustose']
trait_div[,'PropFrut','Obs'] = form_prop[dimnames(trait_div)[[1]],'fruticose']
trait_div[,'PropAsco','Obs'] = rep_prop[dimnames(trait_div)[[1]],'sexual']
trait_div[,'PropAsex','Obs'] = rep_prop[dimnames(trait_div)[[1]],'asexual']

save(trait_div, file='trait_diversity.RData')

load('trait_diversity.RData')

## Convergence in water/competitive traits?
use_traits = c('Attachment','LobeDissect','Form','PropCrust','PropFrut')
use_env = c('WaterCapacity','Bryophytes','CloudFreq_mean','VPD_max','OpenPos')
use_data = data.frame(samples_pm[dimnames(trait_div)[[1]], use_env], N=reproduction[dimnames(trait_div)[[1]],'N'])

use_col = c('#00000055','red')

par(mfrow=c(length(use_traits),length(use_env)+1 ))
par(mar=c(2,2,.5,.5))
par(oma=c(4,4,0,0))
for(i in use_traits){
for(j in c('N',use_env)){
	colorby = (trait_div[,i,'P'] < 0.05) + 1
	plot(trait_div[,i,'Z'] ~ as.numeric(use_data[,j]), xlab=j, ylab=i, pch=16, col=use_col[colorby])
	abline(h=c(0,2,-2), lty=2, col='grey')
	mod = lm(trait_div[,i,'Z'] ~ as.numeric(use_data[,j]))
	if(coef(summary(mod))[2,4] <0.05){
		abline(mod, lwd=2, col=2)
	}
}}

# Slight evidence that trait diversity increases and prop crustose with lichen abundance
# Slight evidence for convergence in attachment with fewer Bryophytes and drier sites
# Slight evidence for higher that expected PropCrustose with fewer bryophytes and in drier climates
# Slight evidence that growith form diversity increases with bryophyte abundance and in wetter climates
# No evidence for convergence in lob dissection
# If anything, water limitation imposes constraint, no competition



## Convergence in reproductive modes with substrate stability?
plot(trait_div[,'PropAsco','Z'], trait_div[,'PropAsex','Z'])

samp_order = dimnames(trait_div)[[1]]
colorby = (trait_div[,'Asexual','P'] < 0.05) + 1
plot(trait_div[,'Asexual','Z'] ~ as.numeric(samples_pm[samp_order,'Shedding']), pch=16, col=use_col[colorby])
plot(trait_div[,'Asexual','Z'] ~ reproduction[samp_order, 'Asex_abun'], pch=16, col=use_col[colorby])

colorby = (trait_div[,'PropAsex','P'] < 0.05) + 1
plot(trait_div[,'PropAsex','Z'] ~ as.numeric(samples_pm[samp_order,'Shedding']), pch=16, col=use_col[colorby])
abline(h=0)

colorby = (trait_div[,'Asco','P'] < 0.05) + 1
plot(trait_div[,'Asco','Z'] ~ as.numeric(samples_pm[samp_order,'Shedding']), pch=16, col=use_col[colorby])
plot(trait_div[,'Asco','Z'] ~ reproduction[samp_order, 'Asco_abun'], pch=16, col=use_col[colorby])

# no evidence that reproductive prevalence is constrained


## Figures for SI:

samp_order = dimnames(trait_div)[[1]]

# Asex, PropAsex ~ Shedding
pdf('./Figures/asex diversity vs bark stability.pdf', height=5, width=5)
par(mfrow=c(2,1))
par(mar=c(.5,5,1,1))
par(oma=c(3,0,0,0))
colorby = (trait_div[,'Asexual','P'] < 0.05) + 1
plot(trait_div[,'Asexual','Z'] ~ as.numeric(samples_pm[samp_order,'Shedding']), 
	pch=16, col=use_col[colorby], xlab='', ylab='z-score', axes=F)
axis(2, las=1); box()
mtext('A',3,0,adj=0)
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
colorby = (trait_div[,'PropAsex','P'] < 0.05) + 1
plot(trait_div[,'PropAsex','Z'] ~ as.numeric(samples_pm[samp_order,'Shedding']), 
	pch=16, col=use_col[colorby], xlab='', ylab='z-score', las=1)
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
mtext('B',3,0,adj=0)
mtext('Bark instability', 1, 2.5)
dev.off()

# Attachment ~ N, Bryophytes, VPD
pdf('./Figures/attachment diversity vs env.pdf', height=3, width=7.8)
par(mfrow=c(1,3))
par(mar=c(4.1,.5,1.5,.5))
par(oma=c(0,3,0,0))
colorby = (trait_div[,'Attachment','P'] < 0.05) + 1
yvals = trait_div[,'Attachment','Z']
plot(yvals ~ reproduction[samp_order,'N'], 
	pch=16, col=use_col[colorby], xlab=expression(N[TOT]), ylab='', las=1)
#mtext('A',3,0,adj=0)
mtext('z-score', 2, 2, cex=.8)
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')

plot(yvals ~ as.numeric(samples_pm[samp_order,'Bryophytes']), 
	pch=16, col=use_col[colorby], xlab='Bryophyte cover', ylab='', las=1, axes=F)
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
axis(1, at=1:6, labels=levels(samples$Bryophytes))
box()
#mtext('B',3,0,adj=0)

plot(yvals ~ samples_pm[samp_order,'VPD_max'], 
	pch=16, col=use_col[colorby], xlab='VPD (Pa)', ylab='', las=1, axes=F)
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
axis(1, at=seq(8,17,2), labels=seq(8,17,2)*100)
box()
#mtext('C',3,0,adj=0)

dev.off()

# LobeDissect ~ VPD
pdf('./Figures/lobe dissect diversity vs env.pdf', height=3.5, width=3.5)
par(mar=c(4.1,4.1,.5,.5))
colorby = (trait_div[,'LobeDissect','P'] < 0.05) + 1
plot(trait_div[,'LobeDissect','Z']~ samples_pm[samp_order,'VPD_max'], 
	pch=16, col=use_col[colorby], xlab='VPD (Pa)', ylab='z-score', axes=F)
axis(1, at=seq(8,17,2), labels=seq(8,17,2)*100)
axis(2, las=1)
box()
abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
dev.off()

# Form, PropCrust, PropFrut ~ N, Bryophytes, VPD
# Make sure to load use_data from above
use_traits = c('Form','PropCrust','PropFrut')
use_ynames = c('Growth form diversity','Prop. crustose','Prop. fruticose')
use_env = c('N','Bryophytes','VPD_max')
use_xnames=c('N[TOT]','Bryophyte~~cover','VPD~~(Pa)')

pdf('./Figures/growth form diversity vs env.pdf', height=7.8, width=7.8)
par(mfrow=c(3,3))
par(mar=c(.5,.5,.5,.5))
par(oma=c(4,4,0,0))
for(i in 1:3){
	yvals = trait_div[,use_traits[i],'Z']
	colorby = (trait_div[,use_traits[i],'P']<0.05) + 1
	for(j in 1:3){
		plot(yvals ~ as.numeric(use_data[,use_env[j]]), pch=16, col=use_col[colorby], xlab='', ylab='', axes=F)
		abline(h=c(-1.96,0,1.96), lty=2, col='grey30')
		box()
		if(j==1){
			axis(2, las=1)
			mtext(use_ynames[i], 2, 3)
		}
		if(i==3){
			mtext(parse(text=use_xnames[j]), 1, 3)
			if(j==1) axis(1)
			if(j==2) axis(1, at=1:6, labels=levels(use_data$Bryophytes))
			if(j==3) axis(1, at=seq(8,17,2), labels=seq(8,17,2)*100)
		}
}}
dev.off()


##

# Make Table for SI



################################################################################
### Spatial Variation in Reproduction ###
library(sp)
library(spdep)

# Read in distance matrix is already made
dmat = read.csv(paste(derive_dir,'sample_distance_matrix.csv', sep=''), row.names=1, check.names=F)

# Create distance matrix between samples, based on implicit scales of data

# A function that calculates the scale that separates two samples
# Define scaling factors (in m)
btw_samps = 0.5
btw_trees = 5
neighbor_mods = 10
btw_mods = matrix(0, nrow=10, ncol=10)
btw_mods[c(1,2,9,10),c(1,2,9,10)] = 1
btw_mods[c(2,3,8,9),c(2,3,8,9)] = 1
btw_mods[c(3,4,7,8),c(3,4,7,8)] = 1
btw_mods[c(4,5,6,7),c(4,5,6,7)] = 1
btw_mods[btw_mods==0] = 2
d3s = rbind(c(1,1,10,10,2,2,5,5), c(4,7,4,7,5,6,5,6))
for(i in 1:ncol(d3s)) btw_mods[d3s[,i],d3s[,i]] = 3 
d4s = rbind(c(1,1,10,10), c(5,6,5,6))
for(i in 1:ncol(d4s)) btw_mods[d4s[,i],d4s[,i]] = 4 
diag(btw_mods) = 0
btw_mods = btw_mods*neighbor_mods
btw_plots = spDists(as.matrix(plot_data[,c('Lon','Lat')]), longlat=T)*1000
rownames(btw_plots) = plot_data$PlotID
colnames(btw_plots) = plot_data$PlotID

use_scales = list(btw_samps = btw_samps, btw_trees = btw_trees, btw_mods=btw_mods, btw_plots=btw_plots)

calc_scale_dist = function(x1, x2, plot_info, samp_info, tree_info, scales){
	these_trees = samp_info[c(x1,x2), 'TreeID']
	these_plots = tree_info[these_trees,'PlotID']
	these_modules = tree_info[these_trees,'Module']

	this_dist = scales$btw_plots[these_plots[1],these_plots[2]]
	if(these_plots[1]==these_plots[2]){
		if(these_modules[1]==these_modules[2]){
			if(these_trees[1]==these_trees[2]){
				this_dist = this_dist + ifelse(x1==x2, 0, scales$btw_samps)
			} else {
				this_dist = this_dist + scales$btw_trees
			}
		} else {
			this_dist = this_dist + scales$btw_mods[these_modules[1],these_modules[2]]
		}
	}

	this_dist
}


# Calculate distance matrix between all samples in analysis
dmat = matrix(NA, nrow=nrow(samples), ncol=nrow(samples))
rownames(dmat) = samples$SampID
colnames(dmat) = samples$SampID

for(i in 1:nrow(dmat)){
for(j in i:ncol(dmat)){
	if(is.na(dmat[i,j])){
		s1 = rownames(dmat)[i]
		s2 = rownames(dmat)[j]
		dmat[i,j] = calc_scale_dist(s1,s2,plot_data, samp_data, tree_data, use_scales)
	}
}
print(paste('Done with',i))
}

# Make symmetric
dmat[lower.tri(dmat)] = t(dmat)[lower.tri(dmat)]

# Check structure of matrix by looking at plot
image(dmat, col=heat.colors(6), breaks=c(0,0.5,5,10,40,300,500000))
image(dmat[1:80,1:80], col=heat.colors(6), breaks=c(0,0.5,5,10,40,300,500000))

# Save distance matrix
write.csv(dmat, paste(deriv_dir, 'sample_distance_matrix.csv', sep=''), row.names=T)

## Calculate Moran's I for each distance class on p
## NOTE: MAY NOT BE CORRECT WAY TO DEAL WITH BINOMIAL DATA
reproduction = merge(reproduction, samples, by='SampID', all.x=T, all.y=F)
rownames(reproduction) = reproduction$SampID
use_data= subset(reproduction, PlotID!='Bladen1')
use_data$p = use_data$Asex_pres/use_data$S

# Subset dmat to samples of interest
use_dmat = dmat[use_data$SampID,use_data$SampID]

# Define spatial lags
intramod_breaks = c(0.0,0.5,5,10,20,30,40,300,500000.)
intramod_lags = cbind(intramod_breaks[1:(length(intramod_breaks)-1)],intramod_breaks[2:length(intramod_breaks)])

distbased_breaks = c(0,.5,5,40,500,5000,50000,500000)
distbased_lags = cbind(distbased_breaks[1:(length(distbased_breaks)-1)],distbased_breaks[2:length(distbased_breaks)])


# Calculate Moran's I at each lag
# Make sure that zvals and use_dmat are in same order
nrand = 999

moran_lagged_intramod = calc_moran_lagged(intramod_lags, use_dmat, use_data$p, nrand)
moran_lagged_distbased = calc_moran_lagged(distbased_lags, use_dmat, use_data$p, nrand)

# Plot correlogram
scale_labs = c('Tree','10m','Plot','500m','5km','50km','500km')
plot_moran_lagged(moran_lagged_distbased, scale_labs, './Figures/MI corellogram asexual pres-abs.pdf')


# Plot correlograms for sexual and asexual together
mi_sexual = calc_moran_lagged(intramod_lags, use_dmat, use_data$Asco_abun/use_data$N, nrand)
mi_asexual = calc_moran_lagged(intramod_lags, use_dmat, use_data$Asex_abun/use_data$N, nrand)

# Plot correlograms on model residuals
x = cbind(use_data$Asco_pres, use_data$S-use_data$Asco_pres)
asco_pres_mod = glm(x ~ DBH+as.numeric(Shedding)+OpenPos*Elevation, data=use_data, family=binomial(link='logit'))
x = cbind(use_data$Asex_pres, use_data$S-use_data$Asex_pres)
asex_pres_mod = glm(x ~ DBH+as.numeric(Shedding)+OpenPos*Elevation, data=use_data, family=binomial(link='logit'))

x = cbind(use_data$Asco_abun, use_data$N-use_data$Asco_abun)
asco_abun_mod = glm(x ~ ., data=use_data[,c(other_covars,'DBH')], family=binomial(link='logit'))
x = cbind(use_data$Asex_abun, use_data$N-use_data$Asex_abun)
asex_abun_mod = glm(x ~ ., data=use_data[,c(other_covars,'DBH')], family=binomial(link='logit'))

asco_res = residuals(asco_abun_mod, type='response')
asex_res = residuals(asex_abun_mod, type='response')

# Subset dmat to samples in models
use_dmat = dmat[names(asco_res),names(asco_res)]

mi_sexual_res = calc_moran_lagged(distbased_lags, use_dmat, asco_res, nrand)
mi_asexual_res = calc_moran_lagged(distbased_lags, use_dmat, asex_res, nrand)


combined = rbind(mi_sexual_res, mi_asexual_res)

mi = apply(combined, c(1,2), function(i) i[[1]]$statistic)
se = apply(combined, c(1,2), function(i) sqrt(var(i[[1]]$res[1:(length(i[[1]]$res)-1)])))
mu = apply(combined, c(1,2), function(i) mean(i[[1]]$res[1:(length(i[[1]]$res)-1)]))
x_coords = 1:ncol(mi)

mycol = repmode_col

doscales=1:4

pdf('./Figures/MI correlogram asexual-sexual abun no bark model residuals.pdf', height=5, width=3.5)
par(mar=c(3,4,1,1))
par(lend=1)

jit = 0.14

plot(x_coords[doscales]-jit, mi[1,doscales], xaxt='n', las=1, ylab="Moran's I", xlab='', 
	cex=1.5, pch=22, lwd=1.2, bg=mycol[1], ylim=range(c(2*se, -2*se, mi)), xlim=c(1-jit, doscales[length(doscales)]+jit))

# add autocorrelaton estimates for asexual
points(x_coords[doscales]+jit, mi[2,doscales], pch=22, lwd=1.2,
	bg=mycol[2], cex=1.5)

# null distribution (no autocorrelation)
arrows(x_coords[doscales]-jit, -2*se[1,doscales], x_coords[doscales]-jit, +2*se[1,doscales], code=3, angle=90, length=jit/2, lwd=1.2)	
arrows(x_coords[doscales]+jit, -2*se[2,doscales], x_coords[doscales]+jit, 2*se[2,doscales], code=3, angle=90, length=jit/2, lwd=1.2)	
points(x_coords[doscales]-jit, mu[1,doscales], pch=20, bg=mycol[1])
points(x_coords[doscales]+jit, mu[2,doscales], pch=20, bg=mycol[2])

axis(1, at=x_coords[doscales], labels=scale_labs[doscales])
legend('topright', c('Sexual','Asexual'), fill=mycol, bty='n')

dev.off()

##########################################################################
### Misc Analyses

morphos = morphos[,colnames(morphos)!='MorphID']

# Only use top2 samples with lichens
use_top2 = top2[top2 %in% rownames(sampXmorph)]
sampXmorph = sampXmorph[use_top2,]

# Remove morphotypes not in top2 samples
keep_morphs = names(colSums(sampXmorph)!=0)
morphos = morphos[keep_morphs,]
sampXmorph = sampXmorph[,keep_morphs]


## Tree ages on which different morphotypes reside
no_rep = names(which(rowSums(sampXmorph[,which(rowSums(morphos[,c('Asexual','Asco')])==0)])>0))
asco = names(which(rowSums(sampXmorph[,which(morphos$Asco)])>0))
asex = names(which(rowSums(sampXmorph[,which(morphos$Asexual)])>0))

par(mfrow=c(1,3))
hist(samples[no_rep,'DBH'], ylim=c(0,220))
hist(samples[asco,'DBH'],ylim=c(0,220),border=2)
hist(samples[asex,'DBH'],ylim=c(0,220),border=4)

## For each tree age, calculate proportion of samples without a morphotype with each reproductive mode


## For each sample, calculate proportion of abundance in each reproductive class
## sort samples by tree age
## 

repro_cat = factor(colSums(t(morphos[,c('Asco','Asexual')])*c(1,10)), levels=c(0,1,10,11))
levels(repro_cat) = c('none','sexual','asexual','both')

repro_samp_tab = apply(sampXmorph, 1, function(x) xtabs(x~repro_cat))
repro_samp_tab_scaled = t(repro_samp_tab)/colSums(repro_samp_tab)

N_tot = colSums(repro_samp_tab)
N_asco = colSums(repro_samp_tab[c('sexual','both'),])
N_asex = colSums(repro_samp_tab[c('asexual','both'),])

pdf('./Figures/Abundance of reproductive types vs DBH.pdf', height=5, width=9)
plot(N_tot~samples[colnames(repro_samp_tab),'DBH'], log='x', pch=16, xlab='DBH', ylab='Abundance')
points(N_asex~samples[colnames(repro_samp_tab),'DBH'],col=4, pch=16)
points(N_asco~samples[colnames(repro_samp_tab),'DBH'], col=2, pch=16)
legend('topright', c('Total','Sexual','Asexual'), col=c(1,2,4), pch=16, bty='n')
dev.off()

size_cat = cut(samples[colnames(repro_samp_tab),'DBH'], breaks=c(0,3,4,5,6,7,10,15,20,30,40,50,60,80))

repro_size_tab_list = apply(repro_samp_tab, 1, function(x){
	aggregate(x, by=list(DBH=size_cat), FUN=function(y) sum(y>0)/length(y))
})
repro_size_tab = sapply(repro_size_tab_list, function(df) df$x)
rownames(repro_size_tab) = levels(size_cat)

barplot(t(repro_size_tab[,1:3]), beside=T)

par(mfrow=c(2,2))
plot(repro_samp_tab_scaled[,'none']~samples[colnames(repro_samp_tab),'DBH'])
plot(repro_samp_tab_scaled[,'asexual']~samples[colnames(repro_samp_tab),'DBH'])
plot(repro_samp_tab_scaled[,'sexual']~samples[colnames(repro_samp_tab),'DBH'])
plot(repro_samp_tab_scaled[,'both']~samples[colnames(repro_samp_tab),'DBH'])









##########################################################################
##### Models of Total Abundance #####
N_tot = sapply(samples$SampID, function(x){
	use_lichens = subset(lichens, SampID==x)
	N = 0
	if(nrow(use_lichens)>0){
		N = sum(use_lichens$AbunCount)
	}
	N
})

N_avg = sapply(samples$SampID, function(x){
	use_lichens = subset(lichens, SampID==x)
	N = 0
	if(nrow(use_lichens)>0){
		N = sum(use_lichens$AbunCount)/nrow(use_lichens)
	}
	N
})
 
samples$N_tot = N_tot
samples$N_avg = N_avg

## Convert ordered factors to numeric for models
samples$Bryophytes = as.numeric(samples$Bryophytes)-1
samples$Shedding = as.numeric(samples$Shedding)

## Add reproduction variables
reproduction = read.csv('./Data/Derived Tables/reproduction.csv')
rownames(reproduction) = reproduction$SampID
use_data = samples
use_data = merge(use_data, reproduction, all.x=T, all.y=F)
rownames(use_data) = use_data$SampID

# Subset to just Piedmont and Mountains plots
samples_pm = subset(samples, SiteID!='Bladen')
samples_pm = droplevels(samples_pm)


## Plot lichen abundance versus environmental predictors
plot(N_avg~N_tot, data=samples_pm)
plot(N_tot~DBH, data=samples_pm, log='x')
plot(N_tot~Trans_tot_cor, data=samples_pm)
plot(N_tot~ Elevation, data=samples_pm)
plot(N_tot ~ Bryophytes, data=samples_pm, pch=16, col='#00000050')
mod = glm(N_tot ~ Bryophytes, data=samples_pm, family='poisson')
mod = glm(N_tot ~ VPD_max, data=samples_pm, family='poisson')

# Covariates 
bark_covars = c('FurrowDepth','Bryophytes','pH','Shedding','Density','Angle')
other_covars = c('Trans_tot_cor','OpenPos','Elevation')

tot_abun_mod = glm(N_tot~., data=samples_pm[,c('N_tot', other_covars, bark_covars, 'DBH')], family=poisson(link='log'))
tot_abun_mod_dropDBH = glm(N_tot~., data=samples_pm[,c('N_tot', other_covars, bark_covars)], family=poisson(link='log'))

summary(tot_abun_mod)
anova(tot_abun_mod, tot_abun_mod_dropDBH, test='Chisq')
exp(coef(tot_abun_mod))

start_parms = c(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
avg_abun_mod = glm(N_avg~., data=samples_pm[,c('N_avg', other_covars, bark_covars, 'DBH')], family=gaussian(link='log'), start=start_parms)
start_parms = c(3, 0, 0, 0, 0, 0, 0, 0, 0 ,0,0)
avg_abun_mod_dropDBH = glm(N_avg~., data=samples_pm[,c('N_avg', other_covars, bark_covars, 'DBH')], family=gaussian(link='log'), start=start_parms)

summary(avg_abun_mod)
anova(avg_abun_mod, avg_abun_mod_dropDBH, test='Chisq')


# Scale variables to similar variance to improve convergence
scaled_data = samples_pm
scaled_data$Elevation = scaled_data$Elevation/1000
scaled_data$DBH = scaled_data$DBH/10
scaled_data$OpenPos = scaled_data$OpenPos/10
scaled_data$Angle = scaled_data$Angle/10
scaled_data$FurrowDepth = scaled_data$FurrowDepth/10

## Fit random slopes model - HAS TROUBLE CONVERGING WITH SAMPLES_PM
tot_abun_mod_plot = glmer(N_tot~FurrowDepth+Angle+pH+Density+Bryophytes+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH + (DBH|PlotID), data=scaled_data, family=poisson(link='log'))
summary(tot_abun_mod_plot) # Site and Plot doen't affect Abun~DBH much

# Fit random intercepts model - HAS TROUBLE CONVERGING
tot_abun_mod_plot = glmer(N_tot~FurrowDepth+Angle+pH+Density+Bryophytes+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH + (1|SiteID/PlotID), data=scaled_data, family=poisson(link='log'))
summary(tot_abun_mod_plot) # Although there are differences in abundance across sites

## Fit random slopes models for different tree species
# Drop samples on tree species with < 10 samples observed across the data set
treefreq = table(samples_pm$TreeTaxonID)
keep_trees = names(treefreq[treefreq>=10]) # 5 trees, be careful with num predictors
keep_obs = samples_pm$TreeTaxonID %in% keep_trees # works b/c scaled_data made from samples_pm

mod_tree1 = glmer(N_tot ~ N_avg + (N_avg|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)
mod_tree2 = glmer(Y ~ N_avg + (1|TreeTaxonID), family='binomial', data=Xdata, subset=keep_obs)
anova(mod_tree1, mod_tree2)

mod_tree2 = glmer(N_tot ~ Trans_tot_cor + OpenPos + Elevation + DBH + (DBH|TreeTaxonID), family='poisson', data=scaled_data, subset=keep_obs)
mod_tree1 = glmer(N_tot ~ Trans_tot_cor + OpenPos + Elevation + DBH + (1|TreeTaxonID), family='poisson', data=scaled_data, subset=keep_obs)
mod_tree0 = glm(N_tot ~ Trans_tot_cor + OpenPos + Elevation + DBH, family='poisson', data=scaled_data, subset=keep_obs)

mod_tree2 = glmer(N_tot ~ FurrowDepth + Bryophytes + pH + Shedding + Density + DBH + (DBH|TreeTaxonID), family='poisson', data=scaled_data, subset=keep_obs)
mod_tree1 = glmer(N_tot ~ FurrowDepth + Bryophytes + pH + Shedding + Density + DBH + (1|TreeTaxonID), family='poisson', data=scaled_data, subset=keep_obs)
mod_tree0 = glm(N_tot ~ FurrowDepth + Bryophytes + pH + Shedding + Density + DBH, family='poisson', data=scaled_data, subset=keep_obs)

AIC(mod_tree2, mod_tree1, mod_tree0)
anova(mod_tree2, mod_tree1) # Different tree species have different slopes, even when controlling for sample and plot scale vars

# Generate null model
mod_null_tree = glmer(N_tot ~ Trans_tot_cor + OpenPos + Elevation + (1|TreeTaxonID), data=scaled_data, family='poisson', subset=keep_obs)

# Plot effect
use_col = c('grey50','red')

yvals = samples_pm[keep_obs,'N_tot']

# Identify species for which abundance increases
tree_ests = ranef(mod_tree2)$TreeTaxonID
pos_trees = rownames(subset(tree_ests, DBH >0))
tree_order = rownames(tree_ests)[order(tree_ests$DBH)]

# Confindence intervals for random effects?
ci = confint(mod_tree2, type='wald')

svg('./Figures/N_tot vs DBH for different tree species.svg', height=5.5, width=5.5)
par(lend=1)
plot(yvals~DBH, data=scaled_data[keep_obs,], pch=16, col='#00000050', ylab=expression(N[TOT]), 
	xlab='Tree Diameter (cm)', axes=F)
axis(2, las=1)
axis(1, at=seq(0,8,2), labels=10*seq(0,8,2))
box()

#for(i in keep_trees){
#	use_obs = scaled_data$TreeTaxonID==i
#	xrange = scaled_data[use_obs, 'DBH']
#	envmeans = colMeans(mod_tree2@frame[,other_covars]) # using means of variables for all tree species
#	envmeans = scaled_data[use_obs, other_covars] # using means of variables for this tree species
#	mod_func = function(x) predict(mod_tree2, data.frame(DBH=x, TreeTaxonID=i, t(envmeans)), type='response')
#	curve(mod_func(x), from=xrange[1], to=xrange[2], add=T, lwd=2, col=use_col[1+(i %in% pos_trees)])
#	#if(i %in% pos_trees) text(xrange[2], 30, labels=nrow(subset(scaled_data, TreeTaxonID==i)), pos=4, col=2)
#}

# Fitting separate means models
# Decided to show prediction without covariates, but other options are commented out and available
for(i in tree_order){
	use_obs = scaled_data$TreeTaxonID==i
	#mod = glm(N_tot ~ Trans_tot_cor + OpenPos + Elevation + DBH, family='poisson', data=scaled_data, subset=use_obs)
	#mod0 = glm(N_tot ~ Trans_tot_cor + OpenPos + Elevation, family='poisson', data=scaled_data, subset=use_obs)
	#mod = glm(N_tot ~ FurrowDepth + Bryophytes + pH + Shedding + Density + DBH, family='poisson', data=scaled_data, subset=use_obs)
	#mod0 = glm(N_tot ~ FurrowDepth + Bryophytes + pH + Shedding + Density, family='poisson', data=scaled_data, subset=use_obs)
	mod = glm(N_tot ~ DBH, family='poisson', data=scaled_data[use_obs,])
	mod0 = glm(N_tot ~ 1, family='poisson', data=scaled_data[use_obs,])
	sig = anova(mod0, mod, test='Chisq')[2,'Pr(>Chi)'] < 0.05
	sign = coef(mod)['DBH'] < 0

	xrange = range(mod$model$DBH)
	xpred = seq(xrange[1],xrange[2], .01)
	#envmeans = colMeans(mod$model[,2:(ncol(mod$model)-1)])
	#envmeans = colMeans(mod$model[,other_covars])
	#ypred = predict(mod, data.frame(DBH=xpred, t(envmeans)), type='response')
	#envmeans = mean(mod$model[,'Trans_tot_cor'])
	#ypred = predict(mod, data.frame(DBH=xpred, Trans_tot_cor=envmeans), type='response')
	ypred = predict(mod, data.frame(DBH=xpred), type='response')
	lines(xpred, ypred, col=ifelse(sig, ifelse(sign, 'red3','blue3'), 'grey50'), lwd=2)
	
	if(sig&!sign) text(xrange[2], ypred[length(ypred)], labels=strsplit(i, '_')[[1]][1], pos=4, col='black')
	#if(sig&sign) print(i)
}

# Add prediction for random-slopes model
# Make sure mod_tree2 has the right predictors (bark-level)
mod_data = mod_tree2@frame
envmeans = colMeans(mod_data[,!(colnames(mod_data) %in% c('DBH','TreeTaxonID','N_tot'))])
mod_func = function(x) predict(mod_tree2, data.frame(DBH=x, t(envmeans)), type='response', re.form=NA)
curve(mod_func(x), add=T, lwd=4)

dev.off()

## How do these tree species differ?
use_order = with(samples_pm, tapply(Shedding, TreeTaxonID, mean))
use_order = names(use_order[order(use_order)])
sp_labs = gsub('_', ' ', use_order)

pdf('./Figures/tree species shedding.pdf', height=6, width=10)
par(mar=c(10,4,3,.5))
plot(Shedding ~ as.numeric(factor(TreeTaxonID, levels=use_order)), data=samples_pm, axes=F, 
	xlab='', ylab='Bark instability', pch=15, col='#00000050')
axis(1, at=1:length(use_order), labels=sp_labs, las=2, col=1:2)
axis(2, at=1:5, labels=c('none','attached\nflakes','loose\nflakes','slight\npeeling','peeling'))
axis(3, at=1:length(use_order), labels=treefreq[use_order], las=2)
for(i in 1:length(use_order)){
	if(treefreq[use_order[i]] < 10){ abline(v=i, col='grey', lty=2) } else {
		use_obs = scaled_data$TreeTaxonID==use_order[i]
		mod = glm(N_tot ~ DBH, family='poisson', data=scaled_data[use_obs,])
		mod0 = glm(N_tot ~ 1, family='poisson', data=scaled_data[use_obs,])
		sig = anova(mod0, mod, test='Chisq')[2,'Pr(>Chi)'] < 0.05
		sign = coef(mod)['DBH'] < 0
		
		if(!sig) abline(v=i, col='grey')
		if(sig) abline(v=i, col=c('blue','red')[sign+1])
	}
}
dev.off()

# Interaction between DBH and Shedding?
mod_inter = glm(N_tot ~ DBH*Shedding, data=samples_pm, family='poisson')

plot(Shedding ~ DBH, data=samples_pm)
use_col = colorRampPalette(blue2red[2:9])(5)

pdf('./Figures/lichen abun vs DBH with Shedding.pdf', height=4, width=4.5)
par(lend=1)
par(mar=c(4.5, 4, .5, .5))
plot(N_tot ~ DBH, data=samples_pm, pch=16, col='#00000050', xlab='DBH (cm)', ylab=expression(N[TOT]), las=1)
#mod_func = function(x,i) predict(mod_inter, data.frame(DBH=x, Shedding=i), type='response')
for(i in c(1:3,5)){
	this_data = subset(samples_pm, Shedding==i)
	DBH_range = range(this_data$DBH)
	this_mod = glm(N_tot ~ DBH, data=this_data, family='poisson')
	pval = 

	mod_func = function(x) predict(this_mod, data.frame(DBH=x), type='response')
	curve(mod_func(x), from=DBH_range[1], to=DBH_range[2], col='black', add=T, lwd=5)
	curve(mod_func(x), from=DBH_range[1], to=DBH_range[2], col=use_col[i], add=T, lwd=3)
}
legend('topright', c('none','attached flakes','loose flakes','peeling'), 
	col=use_col[c(1:3, 5)], lwd=3, bty='n')
dev.off()

## How does successively decreasing tree sizes affect coef estimate for DBH?
sm_mods = sapply(4:max(use_data$DBH), function(thresh){
	use_samps = samples$DBH<=thresh
	abun_mod_sm = glm(N_tot~FurrowDepth+Angle+pH+Density+Bryophytes+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH, family=poisson(link='log'), data=samples[use_samps,])
	
	DBH_coef = as.numeric(coef(abun_mod_sm)['DBH'])
	DBH_se = coef(summary(abun_mod_sm))['DBH','Std. Error']
	CI = confint(abun_mod_sm, 'DBH')
	N = nrow(abun_mod_sm$model)
	P = anova(abun_mod_sm, test='Chisq')$'Pr(>Chi)'[rownames(anova(abun_mod_sm))=='DBH']
	
	c(Size=thresh, DBH_coef=DBH_coef, DBH_se=DBH_se, CI_low95 = CI[1], CI_up95=CI[2], N=N, P=P)
})
sm_mods = t(sm_mods)
sm_mods = as.data.frame(sm_mods)

lg_mods = sapply(68:4, function(thresh){
	use_samps = samples$DBH>=thresh
	abun_mod_sm = glm(N_tot~FurrowDepth+Angle+pH+Density+Bryophytes+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH, family=poisson(link='log'), data=samples[use_samps,])
	
	DBH_coef = as.numeric(coef(abun_mod_sm)['DBH'])
	DBH_se = coef(summary(abun_mod_sm))['DBH','Std. Error']
	CI = confint(abun_mod_sm, 'DBH')
	N = nrow(abun_mod_sm$model)
	P = anova(abun_mod_sm, test='Chisq')$'Pr(>Chi)'[rownames(anova(abun_mod_sm))=='DBH']
	
	c(Size=thresh, DBH_coef=DBH_coef, DBH_se=DBH_se, CI_low95 = CI[1], CI_up95=CI[2], N=N, P=P)
})
lg_mods = t(lg_mods)
lg_mods = as.data.frame(lg_mods)

pdf('./Figures/Tot abun vs DBH models with different size breakpoints.pdf', height=4, width=8)
plot(DBH_coef~Size, data=sm_mods, las=1,ylim=c(-.1,.12))
arrows(sm_mods$Size, sm_mods$CI_low95, sm_mods$Size, sm_mods$CI_up95, length=0.02, angle=90, code=3)
points(DBH_coef~Size, data=sm_mods[sm_mods$P<0.05,], pch=16)
abline(h=0, lty=2)
dev.off()

plot(DBH_coef~Size, data=lg_mods, las=1, ylim=c(-.05,.05), col=2)
arrows(lg_mods$Size, lg_mods$CI_low95, lg_mods$Size, lg_mods$CI_up95, length=0.02, angle=90, code=3, col=2)
abline(h=0, lty=2)
points(DBH_coef~Size, data=lg_mods[lg_mods$P<0.05,], pch=16, col=2)


## Estimate Tot_abun~DBH for 20cm windows
window_mods = sapply(3:(max(samples$DBH)-19), function(thresh){
	use_samps = samples$DBH>=thresh & samples$DBH<(thresh+20)
	abun_mod_wind = glm(N_tot~FurrowDepth+Angle+pH+Density+as.numeric(Bryophytes)+as.numeric(Shedding)+Trans_tot_cor+DBH, family=poisson(link='log'), data=samples[use_samps,])
	
	DBH_coef = as.numeric(coef(abun_mod_wind)['DBH'])
	DBH_se = coef(summary(abun_mod_wind))['DBH','Std. Error']
	CI = confint(abun_mod_wind, 'DBH')
	N = nrow(abun_mod_wind$model)
	P = anova(abun_mod_wind, test='Chisq')$'Pr(>Chi)'[rownames(anova(abun_mod_wind))=='DBH']
	
	c(Size=thresh+10, DBH_coef=DBH_coef, DBH_se=DBH_se, CI_low95 = CI[1], CI_up95=CI[2], N=N, P=P)
})
window_mods = t(window_mods)
window_mods = as.data.frame(window_mods)

plot(DBH_coef~Size, data=window_mods, las=1,ylim=c(-.1,.1))
arrows(window_mods$Size, window_mods$CI_low95, window_mods$Size, window_mods$CI_up95, length=0.02, angle=90, code=3)
points(DBH_coef~Size, data=window_mods[window_mods$P<0.05,], pch=16)
abline(h=0, lty=2)

## Plot Tot_abun ~ DBH

plot(N_tot~DBH, data=samples)

use_obs = rownames(tot_abun_mod$model)
use_data = samples[use_obs,]
colMeans(use_data[,c('FurrowDepth','Angle','pH','Density','Bryophytes','Shedding', 'Trans_tot_cor')])

tot_abun_mod_sm = glm(N_tot~FurrowDepth+Angle+pH+Density+Bryophytes+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH, family=poisson(link='log'), data=samples[samples$DBH<=10,])
summary(tot_abun_mod_sm)

pdf('./Figures/Tot N vs DBH with full model predictions.pdf', height=5, width=5)
par(lend=1)
par(mar=c(4,4,1,1))
plot(N_tot~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col='#00000044', log='x')
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, add=T, col='white')
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='cornflowerblue', add=T)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, from=3, to=10, col='white', add=T)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, from=3, to=10, col='orange', add=T)
dev.off()


## Manuscript Figure
# Total abundance plot including predictions for reproductive mode

# Subset to just Piedmont and Mountains plots
use_data = samples
use_data = merge(use_data, reproduction, all.x=T, all.y=F)
rownames(use_data) = use_data$SampID
use_data = subset(use_data, SiteID!='Bladen')
use_data = droplevels(use_data)

tot_abun_mod = glm(N_tot~., data=use_data[,c('N_tot', other_covars, 'DBH')], family=poisson(link='log'))
rep_response = c('asco','asex')
mycol= c('#ff6600','#5f8dd3','#000000')
mycol_trans = paste(mycol, '80', sep='')
mylty=c(2,1)

envmeans = colMeans(tot_abun_mod$model[,c('Trans_tot_cor','OpenPos','Elevation')])
mod_func = function(use_mod, x){
	yvals = predict(use_mod, data.frame(DBH=x, t(envmeans)), type='response', se.fit=T)
	yvals$fit + t(qnorm(c(0.025,0.5, 0.975)) %*% t(yvals$se.fit))
}

svg('./Figures/Tot N and reproduction vs DBH models no bark covars.svg', height=4, width=4)
par(lend=1)
par(mar=c(4,4,1,1))
plot(N_tot~DBH, data=use_data, las=1, pch=16, xlab='Tree DBH (cm)', 
	ylab=expression(N[TOT]), col='#00000044') #, log='x'

xvals = seq(3,79, length.out=200)
ypred = mod_func(tot_abun_mod, xvals)
#polygon(c(xvals, rev(xvals)), c(ypred[,1], rev(ypred[,3])), col=mycol_trans[3], border=NA)
lines(xvals, ypred[,2], lwd=2, col=mycol[3])

for(i in 1:length(rep_response)){
	yvar=paste(capitalize(rep_response[i]),'abun', sep='_')	
	y = cbind(use_data[,yvar], use_data$N-use_data[,yvar])
	prop_mod = glm(y ~ ., data=use_data[,c(other_covars,'DBH')], family=binomial(link='logit'))
	prop_mod0 = glm(y ~ ., data=use_data[,c(other_covars)], family=binomial(link='logit'))

	signif = anova(prop_mod0, prop_mod, test='Chisq')[2,'Pr(>Chi)']<=0.05	

	ypred = mod_func(tot_abun_mod, xvals)*N_tot_func(prop_mod, xvals)
	#polygon(c(xvals, rev(xvals)), c(ypred[,1], rev(ypred[,3])), col=mycol_trans[i], border=NA)
	lines(xvals, ypred[,2], lwd=3, col=mycol[i], lty=mylty[signif+1])
}	
dev.off()

# Avg abun models including predictions for reproductiive mode
start_parms = c(3, rep(0, length(c(other_covars,bark_covars))), 0)
avg_abun_mod = glm(N_avg~., data=samples[,c(other_covars, bark_covars, 'DBH')], family=gaussian(link='log'), start=start_parms)
use_obs = rownames(avg_abun_mod$model)
model_data = use_data[use_obs,]

pdf('./Figures/Avg N and reproduction vs DBH models.pdf', height=5, width=5)
par(lend=1)
par(mar=c(4,4,1,1))
plot(N_avg~DBH, data=model_data, las=1, pch=16, xlab='Tree Diameter (cm)', 
	ylab='Average Lichen Abundance', col='#00000044')

N_avg_func = function(x) predict(avg_abun_mod, data.frame(DBH=x, Trans_tot_cor=0, OpenPos=84, Elevation=500,FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2), type='response') # 
curve(N_avg_func,	lwd=4, col='black', add=T, from=3, to=79)

for(i in 1:length(rep_response)){
	yabun=paste(capitalize(rep_response[i]),'abun', sep='_')
	ypres=paste(capitalize(rep_response[i]),'pres', sep='_')	
	x = ifelse(model_data[,ypres]==0, 0, model_data[,yabun]/model_data[,ypres])
	x[is.na(x)] = 0
	rep_mod = glm(x ~ ., data=model_data[,c(other_covars,bark_covars,'DBH')], family=gaussian(link='log'), start=start_parms)
	rep_mod0 = glm(x ~ ., data=model_data[,c(other_covars,bark_covars)], family=gaussian(link='log'), start=start_parms[1:(length(start_parms)-1)])

	signif = anova(rep_mod0, rep_mod, test='Chisq')[2,'Pr(>Chi)']<=0.05	

	rep_N_func = function(z) predict(rep_mod, data.frame(DBH=z, Trans_tot_cor=0, OpenPos=84, Elevation=500,FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2), type='response') #
	curve(rep_N_func, from=3, to=79, lwd=4, col=mycol[i], add=T, lty=mylty[signif+1])
}	
dev.off()







# Color points by elevation
mycol = colorRampPalette(c('darkblue','blue','lightblue'))(10)
colorby = cut(use_data$Elevation, breaks=seq(0,1400,length.out=11))

par(lend=1)
par(mar=c(4,4,1,1))
plot(N_tot~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col=mycol[colorby], log='x')
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=120), type='response'),
	lwd=5, col=1, add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=120), type='response'),
	lwd=3, col=mycol[1], add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, col=1, add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=3, col=mycol[4], add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=1300), type='response'),
	lwd=5, col=1, add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=1300), type='response'),
	lwd=3, col=mycol[10], add=T)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=120), type='response'),
	lwd=5, col=1, add=T, from=3, to=10)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=120), type='response'),
	lwd=3, lty=2, col=mycol[1], add=T, from=3, to=10)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, col=1, add=T, from=3, to=10)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=3, lty=2, col=mycol[4], add=T, from=3, to=10)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=1300), type='response'),
	lwd=5, col=1, add=T, from=3, to=10)
curve(predict(tot_abun_mod_sm, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=1300), type='response'),
	lwd=3, lty=2, col=mycol[10], add=T, from=3, to=10)

## Plot model with different bryophyte cover
plot(N_tot~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col='#00000044', log='x')
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=1, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='darkblue', add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='blue', add=T)
curve(predict(tot_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=5, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='lightblue', add=T)


## Plot random slopes model
plot_data$TopoPos = factor(plot_data$TopoPos)
use_col = repmode_col

ordered_sites=unique(plot_data[ordered_plots,'SiteID'])[-1]

pdf('./Figures/Tot N vs DBH random slopes model.pdf', height=6, width=12)
layout(matrix(1:8, 2,4, byrow=T))
par(mar=c(4,4,1,1))
for(this_site in ordered_sites){
	plot(N_tot~DBH, data=subset(use_data, SiteID==this_site), xlim=c(0,80), ylim=c(0,45),
		las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col='#00000044', main='')
	for(this_plot in subset(plot_data, SiteID==this_site)$PlotID){
		curve(predict(tot_abun_mod_plot, 
				data.frame(DBH=x/10, FurrowDepth=4/10, Angle=0/10, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=plot_data[this_plot,'OpenPos']/10, Elevation=plot_data[this_plot, 'Elevation']/1000, PlotID=this_plot, SiteID=plot_data[this_plot,'SiteID']), type='response'),
			lwd=2, col=use_col[plot_data[this_plot,'TopoPos']], add=T)
	}
	this_elev = round(mean(subset(plot_data, SiteID==this_site)$Elevation), 0)
	text(80, 45, paste(this_site, ' (', this_elev, ' m)', sep=''), cex=1.5, adj=c(1,1))
}
dev.off()


plot(N_tot~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col='#00000044')
for(this_plot in subset(plot_data, TopoPos=='sheltered')$PlotID){
curve(predict(abun_mod_plot, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, PlotID=this_plot, SiteID=plot_data[this_plot,'SiteID']), type='response'),
	lwd=2, col=use_col[plot_data[this_plot,'TopoPos']], add=T)
}


## Plot Avg_abun ~ DBH
use_obs = rownames(avg_abun_mod$model)
use_data = samples[use_obs,]

pdf('./Figures/Avg N vs DBH with full model predictions.pdf', height=5, width=5)
par(lend=1)
par(mar=c(4,4,1,1))
plot(N_avg~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Lichen Abundance', col='#00000044')
curve(predict(avg_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, add=T, col='white')
curve(predict(avg_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='black', add=T)
dev.off()


## Plot Bryophyte abundance
bryo_abun_mod = glm(Bryophytes~FurrowDepth+Angle+pH+Density+Shedding+Trans_tot_cor+OpenPos+Elevation+DBH, family=poisson(link='log'), data=samples)
summary(bryo_abun_mod)


par(lend=1)
par(mar=c(4,4,1,1))
plot(Bryophytes~DBH, data=use_data, las=1, pch=16, xlab='Tree Diameter (cm)', ylab='Bryophyte Cover', col='#00000044')
curve(predict(bryo_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=5, add=T, col='white')
curve(predict(bryo_abun_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	lwd=4, col='black', add=T)



## Frequency of no lichens and dead lichens across tree ages

# Define response variables
samples$Empty = samples$N_tot == 0
samples$Dead = samples$NoLichen=='dead'

# Convert ordered factors to numeric for models
samples$Bryophytes = as.numeric(samples$Bryophytes)
samples$Shedding = as.numeric(samples$Shedding)

# Define covariates to control for
bark_covars = c('FurrowDepth','Bryophytes','pH','DBH')
other_covars = c('Trans_tot_cor','OpenPos','Elevation')
focal_var='Shedding'

nolich_mod = glm(samples$Empty~., data=samples[,c(other_covars, bark_covars,  focal_var)], family=binomial(link='logit'))
deadlich_mod = glm(samples$Dead~., data=samples[,c(other_covars, bark_covars, focal_var)], family=binomial(link='logit'))
dead_none_mod = glm(samples$Dead[samples$Empty]~., data=samples[samples$Empty,c(other_covars, bark_covars, focal_var)], family=binomial(link='logit'))

nolich_mod = glm(samples$Empty~., data=samples[,c(other_covars,  focal_var)], family=binomial(link='logit'))
deadlich_mod = glm(samples$Dead~., data=samples[,c(other_covars, focal_var)], family=binomial(link='logit'))
dead_none_mod = glm(samples$Dead[samples$Empty]~., data=samples[samples$Empty,c(other_covars,  focal_var)], family=binomial(link='logit'))


summary(nolich_mod)
summary(deadlich_mod)
summary(dead_none_mod)

# Divide samples into roughly equal size classes
size_classes=c(0,5,10,15,20,30,50,80)
table(cut(samples$DBH, breaks=size_classes))

SizeClass = cut(samples$DBH, breaks=size_classes)

empty_tab=as.matrix(xtabs(~SizeClass+Empty, data=samples))

deadprop = xtabs(Dead~SizeClass, data=samples)
deadprop = deadprop/rowSums(empty_tab)
emptyprop = empty_tab[,'TRUE']/rowSums(empty_tab)

rel_area = rowSums(empty_tab)/sum(empty_tab)
mycol=c('grey80','black')

pdf('./Figures/Prop dead lichen across tree size classes.pdf', height=3, width=7)
par(mar=c(4,4,1,1))
barplot(emptyprop, width=rel_area, col=mycol[1], las=1, border=mycol[1], 
	xlab='Diameter Size Class (cm)', ylab='Proportion of Samples', ylim=c(0,.25))
barplot(deadprop, width=rel_area, col=mycol[2], add=T, axes=F)
legend('topleft', c('No Lichen','Dead Lichen'), fill=mycol, bty='n')
dev.off()

################################################################
### Models on specific tree species

# Remove trees in Bladen1 since Trans_tot_corr not calculated
use_data = subset(reproduction, PlotID!='Bladen1')

# Convert factors to numbers
use_data$Shedding = as.numeric(use_data$Shedding)
use_data$Bryophytes = as.numeric(use_data$Bryophytes)

# Define species and responses of interest
tree_taxa = c('Acer_rubrum','Liriodendron_tulipifera','Quercus_rubra','Quercus_alba')
response_vars = c('N_tot', 'N_avg')

# Define covariates to control for
bark_covars = c('FurrowDepth','Bryophytes','pH','Shedding')
other_covars = c('Trans_tot_cor','OpenPos','Elevation')
focal_var='DBH'
all_vars = c(other_covars, bark_covars, focal_var)

singlesp_mod_tab = array(NA, dim=c(length(response_vars), 9, 2, length(tree_taxa)),
	dimnames=list(response=response_vars, term=c('est','se','up95','low95','P','Dev_resid','nobs','df','nplot'), 
	control=c('none','bark'), species=tree_taxa))

singlesp_abun_mods = vector(mode='list', length(tree_taxa))
names(singlesp_abun_mods) = tree_taxa

for(treesp in tree_taxa){
	singlesp_abun_mods[[treesp]] = vector(mode='list', length(response_vars))
	names(singlesp_abun_mods[[treesp]]) = response_vars
for(yvar in response_vars){

	# Define samples to use
	use_data = subset(samples, TreeTaxonID==treesp)

	# Fit models
	if(yvar=='N_tot'){
		mod_bark = glm(use_data[,yvar]~., data=use_data[,c(other_covars,bark_covars,focal_var)], family=poisson(link='log'))
		mod_bark0 = glm(use_data[,yvar]~., data=use_data[,c(other_covars,bark_covars)], family=poisson(link='log'))
		mod_none = glm(use_data[,yvar]~., data=use_data[,c(other_covars,focal_var)], family=poisson(link='log'))
		mod_none0 = glm(use_data[,yvar]~., data=use_data[,c(other_covars)], family=poisson(link='log'))
	}
	if(yvar=='N_avg'){
		mod_bark = glm(use_data[,yvar]~., data=use_data[,c(other_covars,bark_covars,focal_var)], family=gaussian(link='log'), start=c(3,rep(0,length(all_vars))))
		mod_none = glm(use_data[,yvar]~., data=use_data[,c(other_covars,focal_var)], family=gaussian(link='log'), start=c(3,rep(0,length(other_covars)+1)))
		mod_bark0 = glm(use_data[,yvar]~., data=use_data[,c(other_covars,bark_covars)], family=gaussian(link='log'), start=c(3,rep(0,length(all_vars)-1)))
		mod_none0 = glm(use_data[,yvar]~., data=use_data[,c(other_covars)], family=gaussian(link='log'), start=c(3,rep(0,length(other_covars))))
	}

	ptest_none = anova(mod_none0, mod_none, test='Chisq')[2,'Pr(>Chi)']
	ptest_bark = anova(mod_bark0, mod_bark, test='Chisq')[2,'Pr(>Chi)']

	singlesp_mod_tab[yvar,'est',,treesp] =  c(coef(mod_none)[focal_var], coef(mod_bark)[focal_var])
	singlesp_mod_tab[yvar,'se',,treesp] = c(summary(mod_none)$coef[focal_var,'Std. Error'], summary(mod_bark)$coef[focal_var,'Std. Error'])
	singlesp_mod_tab[yvar,'P',,treesp] = c(ptest_none, ptest_bark)
	singlesp_mod_tab[yvar,c('low95','up95'),,treesp] = cbind(as.numeric(confint(mod_none, focal_var)), as.numeric(confint(mod_bark, focal_var)))
	singlesp_mod_tab[yvar,'Dev_resid',,treesp] = c(deviance(mod_none),deviance(mod_bark))
	singlesp_mod_tab[yvar, 'nobs',,treesp] = c(nrow(mod_none$model), nrow(mod_bark$model))
	singlesp_mod_tab[yvar,'df',,treesp] = c(summary(mod_none)$df.residual, summary(mod_bark)$df.residual)
	singlesp_mod_tab[yvar,'nplot',,treesp] = length(unique(use_data$PlotID))

	singlesp_abun_mods[[treesp]][[yvar]] = list(none=mod_none, bark=mod_bark)
}}

## Reproductive models

# Define covariates to control for - this is used by calc_repmods() function
bark_covars = c('FurrowDepth','Bryophytes','pH','DBH')
other_covars = c('Trans_tot_cor','OpenPos','Elevation')
focal_var='Shedding'

singlesp_repmod_tab = array(NA, dim=c(4, 3, 2, 8, 2), 
	dimnames=list(species=tree_taxa, response=c('asco','asex','norep'), count=c('pres','abun'), term=c('est','se','up95','low95','P','Dev_resid','nobs','df'), control=c('none','bark')))

singlesp_rep_mods = vector(mode='list', length(tree_taxa))
names(singlesp_rep_mods) = tree_taxa

for(treesp in tree_taxa){
	singlesp_rep_mods[[treesp]] = vector(mode='list', length(response_vars))
	names(singlesp_rep_mods[[treesp]]) = response_vars

	this_data = subset(use_data, TreeTaxonID==treesp)
	
	these_mods = calc_repmods(this_data)

	singlesp_repmod_tab[treesp,,,,] = these_mods
}

DBH_ssp_mods = singlesp_repmod_tab 
Shed_ssp_mods = singlesp_repmod_tab

save(DBH_ssp_mods, Shed_ssp_mods, singlesp_mod_tab, singlesp_abun_mods, file='./Data/single species reproduction model tables.RData')

## Make plots
load('./Data/single species reproduction model tables.RData')


# Reproductive models
use_mods = DBH_ssp_mods
use_ests = melt(use_mods[,,,'est','none'])
use_ests$response = factor(use_ests$response, levels=c('asco','asex','norep'))
levels(use_ests$response) = c('Sexual','Asexual','No reproduction')
use_CI = melt(use_mods[,,,c('low95','up95'),'none'])
use_CI$model = factor(use_CI$response, levels=c('asco','asex','norep'))
levels(use_CI$response) = c('Sexual','Asexual','No reproduction')

mypch=c(1,16)
mycol='black'

jit=c(-0.1,0.1)

pdf('./Figures/DBH models by species.pdf', height=3, width=9)
xyplot(species~value|response, groups=count, data=use_ests, panel=function(x,y,groups,subscripts,...){
	panel.abline(v=0, col='grey', lty=2)
	panel.points(x, as.numeric(y)+jit[groups[subscripts]], col=mycol, pch=mypch[groups[subscripts]])
	panel.arrows(use_CI[use_CI$term=='low95','value'][subscripts], as.numeric(y)+jit[groups[subscripts]],
		use_CI[use_CI$term=='up95','value'][subscripts], as.numeric(y)+jit[groups[subscripts]],
		code=3, length=0.05, angle=90) 

}, strip=strip.custom(bg='transparent'), xlim=c(-0.08,0.08), xlab='DBH Effect', ylab='',
scales=list(alternating=1), 
key=list(space='right', title='Model response variable', cex.title=1.2, points=list(pch=mypch, col=mycol), text=list(c('Proportion of species','Proportion of abundance'))))

dev.off()

# Abundance models
use_mods = singlesp_mod_tab
use_P = melt(use_mods[,'P','none',])

# Refit models to all data
bark_covars = c('FurrowDepth','Bryophytes','pH','Shedding')
tot_abun_mod = glm(use_data$N_tot~., data=use_data[,c(other_covars,'DBH')], family=poisson(link='log'))
#avg_abun_mod = glm(use_data$N_avg~., data=use_data[,c(other_covars,'DBH')], family=gaussian(link='log'), start=c(3,0,0,0,0,0,0,0,0))

# Only plot data used to fit models (e.g. no misssing covariates)
model_data = use_data[rownames(tot_abun_mod$model),]

mylty = c(2,1)
mycol = c(repmode_col,'grey40')
rep_response = c('asco','asex')

pdf('./Figures/Abundance and reproduction vs DBH models single species no bark covars.pdf', height=8, width=3.5)
par(lend=2)
par(mfrow=c(4,1))
par(mar=c(2,2,1,1))
par(oma=c(2,2.5,0,0))
for(treesp in tree_taxa){
	this_data = subset(model_data, TreeTaxonID==treesp)
	
	plot(N_tot~DBH, data=this_data, las=1, pch=16, xlab='', ylab='', 
		col='#00000044', xlim=c(3,80), ylim=c(0,42))
	
	# model of total abundance
	use_mod = singlesp_abun_mods[[treesp]]$N_tot$none

	signif = subset(use_P, response=='N_tot'&species==treesp)$value<=0.05
	
	N_tot_func = function(x) predict(use_mod, data.frame(DBH=x, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response', from=3, to=max(this_data$DBH))
	curve(N_tot_func, from=3, to=max(this_data$DBH), lwd=3, col='black', add=T, lty=mylty[signif+1])

	# add model of avg abundance
	#use_mod = singlesp_abun_mods[[treesp]]$N_avg$none
	#signif = subset(use_P, response=='N_avg'&species==treesp)$value<=0.05
	#curve(predict(use_mod, data.frame(DBH=x, FurrowDepth=4, Angle=0, pH=6.3, Density=0.5, Bryophytes=3, Shedding=2, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response'),
	#	from=3, to=max(this_data$DBH),lwd=3, col='grey30', add=T, lty=mylty[signif+1])

	# proportion sexual vs asexual
	for(i in 1:length(rep_response)){
		yvar=paste(capitalize(rep_response[i]),'abun', sep='_')	
		x = cbind(this_data[,yvar], this_data$N-this_data[,yvar])
		prop_mod = glm(x ~ ., data=this_data[,c(other_covars,'DBH')], 
			family=binomial(link='logit'))
		
		signif = DBH_ssp_mods[treesp,rep_response[i],'abun','P','none']<=0.05

		rep_tot_func = function(z) N_tot_func(z)*predict(prop_mod, data.frame(DBH=z, Trans_tot_cor=0, OpenPos=84, Elevation=500), type='response')
		curve(rep_tot_func, from=3, to=max(this_data$DBH), lwd=3, col=mycol[i], add=T, lty=mylty[signif+1])
	}	

	# add species label
	splabel=paste(unlist(strsplit(treesp,'_')), collapse=' ')
	mtext(splabel, 3, -1.5, adj=.95, font=3)
}
mtext('Tree Diameter (cm)',1,1,outer=T, adj=.6)
mtext('Lichen Abundance',2,1,outer=T)
dev.off()

###############################################################
### SEM

library(lavaan)

# Define model
attach_sem = '

# Cloud freq -> Bryophyte + Attachment
Bryophytes + Attachment + N_tot ~ CloudFreq_mean

# Bryophytes -> Attachment
Attachment ~ Bryophytes

# Bryophytes -> N_tot
#N_tot ~ Bryophytes

# N -> Attachment
Attachment ~ N_tot

'

# Define data
# Make sure sampXtrait_abun is created in code above and log-transformed
model_data = data.frame(Attachment = sampXtrait_abun[,'Attachment'], samples_pm[rownames(sampXtrait_abun),])


sem_mod = sem(attach_sem, data=model_data, estimator='ML', se='robust')
summary(sem_mod, standardized=T, rsq=TRUE, fit.measures=T)
semPaths(sem_mod, 'std', intercepts=F)







################################################################
### Bayesian

library('R2OpenBUGS')
library('coda')
library('MASS')

### Logistic Regression for binary traits



################################################################
### Functions

# A panel function for plotting correlations
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r = (cor(x, y))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex * abs(r))
}

calc_moran = function(x, wmat){
	N = length(x)
	covmat = (x - mean(x))%*%t((x - mean(x)))
	numerator = sum(covmat*wmat)
	denominator = sum(diag(covmat))
	scalar = N/sum(wmat)
	MI = numerator/denominator*scalar
	MI
}
 
calc_moran_lagged = function(use_lags, use_dmat, zvals, nrand){
	lapply(1:nrow(use_lags), function(i){

	this_lag = use_lags[i,]
	
	# Derive spatial weights matrix at this lag from dmat
	wmat = use_dmat>this_lag[1] & use_dmat<=this_lag[2]
	wmat = wmat*1

	# Subset only to observations with neighbors
	zvals = zvals[rowSums(wmat)>0]
	wmat = wmat[rowSums(wmat)>0, colSums(wmat)>0]
		

	# Use spdep functions to do moran calculations
	lw = mat2listw(wmat)
	mi_boot = moran.mc(zvals, lw, nrand, return_boot=F) 
	
	mi_boot	


	# Calculate Moran's I by hand
	#N = length(zvals)
	#MI = calc_moran(use_zvals, use_wmat)
	#MI_null = -1/(N-1)

	# Calculate CI based on randomization
	#rand_p = sapply(1:nrand, function(i) sample(zvals, N, replace=F))
	#MI_rand = apply(rand_p, 2, function(x) calc_moran(x, wmat) )

	#MI_p = 2*min(c(sum(MI_rand<=MI), sum(MI_rand>=MI)))/nrand

	#list(MI=MI, MI_null=MI_null, MI_p=MI_p, MI_rand=MI_rand)

	})
}


plot_moran_lagged = function(x, lag_labels, filename){
	
	mi = sapply(x, function(i) i$statistic)
	se = sapply(x, function(i) sqrt(var(i$res[1:(length(i$res)-1)])))
	mu = sapply(x, function(i) mean(i$res[1:(length(i$res)-1)]))
	x_coords = 1:length(mi)
	
	pdf(filename, height=5, width=5)
	par(mar=c(3,4,1,1))

	plot(x_coords, mi, xaxt='n', las=1, ylab="Moran's I", xlab='', cex=1.2, pch=16, ylim=range(c(mi+2*se, mi-2*se)))
	arrows(x_coords, mi-2*se, x_coords, mi+2*se,code=3, angle=90, length=.07, lwd=1.5)	
	points(x_coords, mu, pch=8, col=2)
	axis(1, at=x_coords, labels=lag_labels)

	dev.off()
}


#  A function that calculates reproductive mode models  on a given data subset
# covariates lists must exist already: bark_covars, other_covars, focal_var

calc_repmods = function(use_data){
	
	repmod_table = array(NA, dim=c(3, 2, 8, 2), 
		dimnames=list(response=c('asco','asex','norep'), count=c('pres','abun'), term=c('est','se','up95','low95','P','Dev_resid','nobs','df'), control=c('none','bark')))

	modcombos = expand.grid(dimnames(repmod_table)[[1]], dimnames(repmod_table)[[2]])

	for(i in 1:nrow(modcombos)){
		this_response = as.character(modcombos[i,1])
		this_count = as.character(modcombos[i,2])

		response_name = capitalize(paste(this_response, this_count, sep='_'))
	
		if(this_count=='pres') x = cbind(use_data[,response_name], use_data$S-use_data[,response_name])
		if(this_count=='abun') x = cbind(use_data[,response_name], use_data$N-use_data[,response_name])
	
		mod_none0 = glm(x ~ ., data=use_data[,other_covars], family=binomial(link='logit'))
		mod_none = glm(x ~ ., data=use_data[,c(other_covars, focal_var)], family=binomial(link='logit'))
		mod_bark0 = glm(x ~ ., data=use_data[,c(other_covars, bark_covars)], family=binomial(link='logit'))
		mod_bark = glm(x ~ ., data=use_data[,c(other_covars, bark_covars, focal_var)], family=binomial(link='logit'))

		ptest_none = anova(mod_none0, mod_none, test='Chisq')[2,'Pr(>Chi)']
		ptest_bark = anova(mod_bark0, mod_bark, test='Chisq')[2,'Pr(>Chi)']

		repmod_table[this_response, this_count,'est',] =  c(coef(mod_none)[focal_var], coef(mod_bark)[focal_var])
		repmod_table[this_response, this_count,'se',] = c(summary(mod_none)$coef[focal_var,'Std. Error'], summary(mod_bark)$coef[focal_var,'Std. Error'])
		repmod_table[this_response, this_count,'P',] = c(ptest_none, ptest_bark)
		repmod_table[this_response, this_count, c('low95','up95'),] = cbind(as.numeric(confint(mod_none, focal_var)), as.numeric(confint(mod_bark, focal_var)))
		repmod_table[this_response, this_count, 'Dev_resid',] = c(deviance(mod_none),deviance(mod_bark))
		repmod_table[this_response, this_count, 'nobs',] = c(nrow(mod_none$model), nrow(mod_bark$model))
		repmod_table[this_response, this_count, 'df',] = c(summary(mod_none)$df.residual, summary(mod_bark)$df.residual)
	}

	repmod_table
}





