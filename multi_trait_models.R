## This script assesses variation in multiple traits for the Lichen FD Project


git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/'

# Load data & make data frames for analysis
source(paste(git_dir, 'load_data.R', sep=''))

abunrich = read.csv(paste(derive_dir, 'sample_richness_abundance.csv', sep=''), row.names=1)

########################################################
#### Functional Diversity

library(cluster) # daisy
library(ape) # pcoa
#library(FD)

## Calculate distance matrix among morphotypes

# Scaling factor - amount by which 2ndary traits affect distance relative to 1st traits
sfact = 0.5

use_traits = traitdf[traitdf$TraitName %in% colnames(morphos),]
rownames(use_traits) = use_traits$TraitName

# weights are (sfact)^level-1
dist_L2 = daisy(morphos, type=list(asymm=use_traits[use_traits$type=='asymm','TraitName'],
					symm=use_traits[use_traits$type=='symm','TraitName']),
		weights=sfact^(use_traits[colnames(morphos),'Level']-1)
)

dist_L2mat = as.matrix(dist_L2)

write.csv(dist_L2mat, paste(derive_dir, 'morpho_distance_matrix_sfact0.5.csv'), row.names=T)

## Calculate Functional Diversity

# Put matrices in same order
sampXmorph = sampXmorph[,colnames(dist_L2mat)]

# Calculate relative abundance
sampXmorph_scaled = as.matrix(sampXmorph/rowSums(sampXmorph))

# Calculate rao's Q
rao_L2 = calc_rao(sampXmorph_scaled, dist_L2mat)

# Add to data frame
# note that this compiles several iterations where the distance matrix has been changed
FD = data.frame(SampID=rownames(sampXmorph_scaled), rao_L2_s.5=rao_L2)
FD$rao_L2_s.25 = rao_L2
FD$rao_L2_s1 = rao_L2

write.csv(FD, './Data/raos_q_samples.csv', row.names=F)

### Compare actual FD to null models

# Read in FD data
FD = read.csv('./Data/raos_q_samples.csv')

# Read in null FD distributions from KillDevil runs and calculate p and z-scores
null_fd_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/FD NUll Sims/'

scales = c('withinplot','withinsite')
alphas = c('0.5', '0.25')
model_methods = c('swsh_samp_r','swsh_samp_c')

runs = expand.grid(scales, alphas, model_methods)

# Make array to hold data
FD_scaled = array(NA, dim=c(nrow(FD), 2, length(scales), length(alphas), length(model_methods)), 
	dimnames=list(rownames(FD), c('p','z'), scales, alphas, model_methods))

sum(colnames(this_null)== FD$SampID)

for(i in 1:nrow(runs)){
	this_file = paste('null_rao',runs[i,1], runs[i,2],runs[i,3],sep='-')
	this_null = read.csv(paste(null_fd_dir,this_file,'.csv', sep=''), check.names=F)

	# Check to make sure samples are in the same order
	if(sum(colnames(this_null)!= FD$SampID)>0){
		print(paste("ERROR: Samples in",this_file,"not same as in FD."))
	} else {

		# Choose observed rao's q based on scaling factor = 0.5 or 0.25
		if(runs[i,2]==0.5) use_obs = FD$rao_L2_s.5
		if(runs[i,2]==0.25) use_obs = FD$rao_L2_s.25
		
		# Calculate p and z-scores
		pzs = sapply(1:length(use_obs), function(j){
			x = use_obs[j]
			null = this_null[,j]
			c(p=calc_p(x, null),z=calc_z(x, null))
		})
		
		FD_scaled[,,runs[i,1], runs[i,2], runs[i,3]] = t(pzs)
	}
}

# Save array
save(FD_scaled, file=paste(null_fd_dir,'FD_scaled.RData', sep=''))
load(paste(null_fd_dir,'FD_scaled.RData', sep=''))

## Boxplots comparing FD distributions across calculation methods
# FD_scaled[SampID, p/z, withinplot/withinsite, 0.5/0.25, swsh_samp_r/swsh_samp_c]

bp_data = data.frame(FD_scaled[,'z',,'0.5','swsh_samp_c'])
boxplot(bp_data, las=1, ylab='Z-score', xaxt='n')
axis(1, at=1:4, labels=rep(c('Plot','Site'),2), tick=F, line=-.5)
#axis(1, at=c(0,1.5, 3.5), labels=c('Abundance preserved:','across morphotypes','within samples'), tick=F, line=1)


library('colorRamps')

colorby = factor(samples[FD$SampID,'SiteID'], levels=c('Bladen','Yadkin','John','Eno','Hang','Pisgah','New','High','Jeff'))
pchby = factor(plot_data[samples[FD$SampID,'PlotID'],'TopoPos'])
mycol = matlab.like(9)
mypch=c(0,16)

par(mfrow=c(1,2))
plot(bp_data[pchby=='exposed',], col=mycol[as.numeric(colorby[pchby=='exposed'])], pch=mypch[2])
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col='grey20', lwd=2)
legend('topleft',levels(colorby), col=mycol, pch=16, bty='n')
plot(bp_data[pchby=='sheltered',], col=mycol[as.numeric(colorby[pchby=='sheltered'])], pch=mypch[2])
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col='grey20', lwd=2)



# Compare effect of 2ndary trait scaling exponent on FD inference
FD_scaled = FD_rep_scaled

pdf('./Figures/FD z-scores compare 2ndary trait scaling exponent.pdf', height=8, width=8)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1))
plot(FD_scaled[,'z','withinplot',,'swsh_samp_r'], main='Swap rows within plots')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinsite',,'swsh_samp_r'], main='Swap rows within sites')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinplot',,'swsh_samp_c'], main='Swap columns within plots')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinsite',,'swsh_samp_c'], main='Swap columns within sites')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
dev.off()


# Compare effect of null model swap methods on FD inference
pdf('./Figures/FD z-scores compare null model swap methods.pdf', height=8, width=8)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1))
plot(FD_scaled[,'z','withinplot','0.25',], main='Within plots, 0.25')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinsite','0.25',], main='Within sites, 0.25')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinplot','0.5',], main='Within plots, 0.5')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
plot(FD_scaled[,'z','withinsite','0.5',], main='Within sites, 0.5')
abline(h=c(-2,2), v=c(-2,2))
abline(0,1, col=2)
dev.off()


## Summarize number of under/over dispersed samples according to different nulls
## make sure to remove Bladen site plot using code below first.
library(reshape)

dimnames(FD_scaled)[[1]] = FD$SampID
names(dimnames(FD_scaled)) = c('SampID','Statistic','Swap_scale','Alpha','Swap_method')

FD_melt = melt(FD_scaled)

FD_tab = cast(FD_melt, SampID~Swap_scale+Statistic, subset=Alpha=='0.5'&Swap_method=='swsh_samp_r')
FD_tab = merge(FD_tab, samples, all.x=T)


sum(with(FD_tab, withinplot_p<0.05&withinplot_z<0))# 6 under
sum(with(FD_tab, withinplot_p<0.05&withinplot_z>0))# 65 over
unders = FD_tab[with(FD_tab, withinplot_p<0.05&withinplot_z<0),]
overs = FD_tab[with(FD_tab, withinplot_p<0.05&withinplot_z>0),]

sum(with(FD_tab, withinsite_p<0.05&withinsite_z<0))# 6 under
sum(with(FD_tab, withinsite_p<0.05&withinsite_z>0))# 57 over

## Examine FD=0
FD = merge(FD, samples, all.x=T)
rownames(FD) = FD$SampID
FD0_samps = subset(FD, rao_L2_s.5==0)$SampID
FD_scaled[FD0_samps,'z','withinplot','0.5',] # all are negative
summary(abunrich[FD0_samps,]) # all have 1 morphotype and max 3 species



##################################
### FD z-scores vs environment ###


dimnames(FD_scaled)[[1]] = FD$SampID
names(dimnames(FD_scaled)) = c('SampID','Statistic','Swap_scale','Alpha','Swap_method')

mycol = read.csv('../../blue2red_10colramp.txt')
mycol = apply(mycol, 1, function(x) rgb(x[1],x[2],x[3], maxColorValue=255))

# Remove Bladen plot
keep_samps = samples[FD$SampID, 'SiteID']!='Bladen'
FD_scaled = FD_scaled[keep_samps,,,,]
FD = FD[keep_samps,]

## Test underdispersion on drier samples


# Use FD null based on swapping traits within the same plot and preserving sample abundance
use_z = FD_scaled[,'z','withinplot','0.5','swsh_samp_r']
use_z = FD$rao_L2_s.5

# Define predictors
sampvars = c('Bryophytes','WaterCapacity')
climvars = c('CloudFreq_mean','VPD_max','OpenPos')
Xdata = samples[FD$SampID, c(sampvars,climvars)]

# Make plots and add lines for models
pdf('./Figures/FD_0.5_r vs water variables.pdf', height=5.5, width=5.5)
par(mfrow=c(3,3))
for(j in climvars){
for(i in sampvars){
	par(mar=c(4, 4, 1, 1))

	keep_obs = !(is.na(Xdata[,i])|is.na(Xdata[,j])|is.na(use_z))
	z = use_z[keep_obs]
	sampX = unclass(Xdata[keep_obs,i])
	climX = Xdata[keep_obs,j]

	use_col = mycol[10:1]
	if(j=='CloudFreq_mean') use_col=mycol
	
	ptcols = paste(use_col[cut(climX, 10, include.lowest=T)], '50', sep='')

	plot(z~sampX, pch=16, col=ptcols, xlab=i, ylab='FD z-score', las=1)
	abline(h=c(-2,0,2), col='grey50', lty=2)

	clim_range = range(climX)
	
	mod_inter = lm(z~sampX*climX)
	P_inter = anova(mod_inter)[3,'Pr(>F)']
	if(P_inter <0.05){
		mod_func = function(x, x_j) predict(mod_inter, data.frame(sampX=x, climX=x_j))
		curve(mod_func(x, clim_range[1]), from=min(sampX), to=max(sampX), add=T, col=use_col[1], lwd=2)		
		curve(mod_func(x, clim_range[2]), from=min(sampX), to=max(sampX), add=T, col=use_col[10], lwd=2)	
	} else {
		mod = lm(z~sampX+climX)
		P = anova(mod)[1:2,'Pr(>F)']
		mod_func = function(x, x_j) predict(mod, data.frame(sampX=x, climX=x_j))

		if(P[2] < 0.05){
			curve(mod_func(x, clim_range[1]), from=min(sampX), to=max(sampX), add=T, col=use_col[1], lwd=2)		
			curve(mod_func(x, clim_range[2]), from=min(sampX), to=max(sampX), add=T, col=use_col[10], lwd=2)	
		} else {
			curve(mod_func(x, mean(samples[,j])), from=min(sampX), to=max(sampX), add=T, lwd=2, lty=(2:1)[(P[1]<0.05)+1])	
		}
	}
	
}
# Add color legend
color.bar(use_col, min=clim_range[1], max=clim_range[2], nticks=length(use_col)+1, title=j)

}
dev.off()

## Effect of topographic openness when morphotypes swapped within sites.


climX = scale(Xdata[FD$SampID, 'OpenPos'], center=T, scale=F)
use_z = FD_scaled[,'z','withinsite','0.5','swsh_samp_r']
plot(use_z~climX, xlab='Topographic Openess', ylab='FD z-score', las=1)
abline(h=c(-2,0,2), col='grey50', lty=2)

topo_mod = lm(use_z ~ SiteID + climX - 1, data=samples[FD$SampID,])
# FD increases with topographic openness



## Test underdispersion on less stable substrates
use_z = FD_scaled[,'z','withinplot','0.5','swsh_samp_r']
sampX = unclass(samples[FD$SampID, 'Shedding'])
ptcols = mycol[10:1][cut(FD$rao_L2_s.5, 10, include.lowest=T)]

pdf('./Figures/FD_0.5_r vs bark shedding.pdf', height=4, width=6)
par(mar=c(5,4,1,1))
plot(use_z~sampX, pch=16, col=paste(ptcols,'50', sep=''), las=1, xlab='Bark Instability', ylab='FD z-score', axes=F)
abline(h=c(-2,0,2), col='grey50', lty=2)
mod = lm(use_z~sampX)
xvals=1:5
points(xvals, tapply(use_z, sampX, mean), pch=3, col='blue', cex=1.5)
points(xvals, predict(mod, data.frame(sampX=xvals)), pch=1, col=2, cex=2)
legend('topright', c('Model','Means'), pch=c(1,3), col=c('red','blue'), bty='n')

axis(2, las=1)
axis(1, at=1:5, labels=levels(sampX))
box()
dev.off()


pdf('./Figures/FD_0.5 vs bark shedding.pdf', height=4.5, width=7.5)
par(mar=c(5,4,1,1))
plot(rao_L2_s.5~Shedding, data=FD, pch=16, col='#00000020', las=1, xlab='Bark Instability', ylab='FD', axes=F)
axis(2, las=1)
axis(1, at=1:5, labels=levels(sampX))
box()
dev.off()

## Test underdispersion on younger trees
use_z = FD_scaled[,'z','withinplot','0.5','swsh_samp_c']
use_data = samples[FD$SampID,]

FD_DBHmod_sm = lm(use_z[use_data$DBH<=6] ~ DBH, data=use_data[use_data$DBH<=6,])
FD_DBHmod_lg = lm(use_z[use_data$DBH>6] ~ DBH, data=use_data[use_data$DBH>6,])

pdf('./Figures/FD z-score site-scaled vs DBH.pdf', height=4.5, width=5.5)
par(mar=c(4,4,1,1))
plot(use_z~DBH, data = use_data, log='x',
	ylab='FD Z-score', xlab='Tree Diameter (cm)', las=1)
abline(h=c(-1.96,1.96), lty=2, col='grey20')
curve(predict(FD_DBHmod_sm, data.frame(DBH=x)), from=3, to=6, add=T, lwd=4)
curve(predict(FD_DBHmod_lg, data.frame(DBH=x)), from=7, to=max(use_data$DBH), add=T, lwd=4)
dev.off()

summary(FD_DBHmod_sm)
summary(FD_DBHmod_lg)



#################################################################
### FD on subsets of traits

library(lme4)

# Read in morphos so that traits are not factors
morphos = read.csv(paste(derive_dir, 'morphotypes.csv', sep=''), row.names=1)
morphos = morphos[,colnames(morphos)!='MorphID']

## FD of reproductive traits

use_traits = traitdf[grep('Asco|Asexual',traitdf$TraitName),]
rownames(use_traits) = use_traits$TraitName

# Make matrix of reproductive morphotypes
rep_morphos = unique(morphos[,use_traits$TraitName])
rownames(rep_morphos) = paste('RM', 1:nrow(rep_morphos), sep='')

write.csv(rep_morphos, paste(derive_dir, 'reproductive_morphotypes.csv',sep=''), row.names=T)
rep_morphos = read.csv(paste(derive_dir, 'reproductive_morphotypes.csv',sep=''), row.names=1)

# Generate sample X morphotype matrix
lichen_data$RMorphID = sapply(1:nrow(lichen_data), function(i){
	 match_morpho(lichen_data[i,], rep_morphos)
})

sampXrepM = xtabs(AbunCount ~ SampID + RMorphID, data=lichen_data)

write.csv(sampXrepM, paste(derive_dir,'sampXrepM.csv', sep=''), row.names=T)

# Cover morphos to factors so that distances can be calculated
sampXrepM = read.csv(paste(derive_dir,'sampXrepM.csv', sep=''), row.names=1)

# Calculate distance between morphotypes
rep_morphos$Asco = factor(rep_morphos$Asco)
rep_morphos$Asexual = factor(rep_morphos$Asexual)
rep_morphos$AscoForm = factor(rep_morphos$AscoForm)
rep_morphos$AscoCover = factor(rep_morphos$AscoCover)
rep_morphos$AsexualForm = factor(rep_morphos$AsexualForm)
rep_morphos$AsexualAbun = factor(rep_morphos$AsexualAbun, levels=c('few','several','many','covered'), ordered=T)

# weights are (sfact)^level-1
sfact=1
dist_repM = daisy(rep_morphos, type=list(asymm=use_traits[use_traits$type=='asymm','TraitName'],
					symm=use_traits[use_traits$type=='symm','TraitName']),
		weights=sfact^(use_traits[colnames(rep_morphos),'Level']-1)
)

dist_repM = as.matrix(dist_repM)

# Calculate Functional Diversity

# Put matrices in same order
use_RM = colnames(sampXrepM)
dist_repM = dist_repM[use_RM,use_RM]

# Calculate relative abundance
sampXrepM_scaled = as.matrix(sampXrepM/rowSums(sampXrepM))

# Calculate rao's Q
rao_rep = calc_rao(sampXrepM_scaled, dist_repM)

# Add to data frame
# note that this compiles several iterations where the distance matrix has been changed
#FD_rep = data.frame(SampID=rownames(sampXrepM_scaled), rao_L2_s.5=rao_rep)
FD_rep$rao_L2_s.25 = rao_rep
FD_rep$rao_L2_s1 = rao_rep

write.csv(FD_rep, './Data/raos_q_reproductive_traits_samples.csv', row.names=F)

## Compare actual FD to null models

# Read in observed FD
FD_rep = read.csv('./Data/raos_q_reproductive_traits_samples.csv')
rownames(FD_rep) = FD_rep$SampID

# Only use top2 samples with lichens
use_top2 = top2[top2 %in% rownames(sampXrepM)]
sampXrepM = sampXrepM[use_top2,]
FD_rep = FD_rep[use_top2,]

# Read in null FD distributions from KillDevil runs and calculate p and z-scores
null_fd_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/FD NUll Sims/'

scales = c('withinplot')
alphas = c('0.5', '0.25')
model_methods = c('swsh_samp_r','swsh_samp_c')

runs = expand.grid(scales, alphas, model_methods)

# Make array to hold data
FD_rep_scaled = array(NA, dim=c(nrow(FD), 2, length(scales), length(alphas), length(model_methods)), 
	dimnames=list(rownames(FD), c('p','z'), scales, alphas, model_methods))

# Check to make sure observed and null samples are in the same order after reading in the first file by hand 
sum(colnames(this_null)== FD_rep$SampID) 

for(i in 1:nrow(runs)){
	this_file = paste('null_rep_rao',runs[i,1], runs[i,2],runs[i,3],sep='-')
	this_null = read.csv(paste(null_fd_dir,this_file,'.csv', sep=''), check.names=F)

	# Check to make sure samples are in the same order
	if(sum(colnames(this_null)!= FD$SampID)>0){
		print(paste("ERROR: Samples in",this_file,"not same as in FD."))
	} else {

		# Choose observed rao's q based on scaling factor = 0.5 or 0.25
		if(runs[i,2]==0.5) use_obs = FD_rep$rao_L2_s.5
		if(runs[i,2]==0.25) use_obs = FD_rep$rao_L2_s.25
		
		# Calculate p and z-scores
		pzs = sapply(1:length(use_obs), function(j){
			x = use_obs[j]
			null = this_null[,j]
			c(p=calc_p(x, null),z=calc_z(x, null))
		})
		
		FD_rep_scaled[,,runs[i,1], runs[i,2], runs[i,3]] = t(pzs)
	}
}

# Save array
save(FD_rep_scaled, file=paste(null_fd_dir,'FD_rep_scaled.RData', sep=''))
load(paste(null_fd_dir,'FD_rep_scaled.RData', sep=''))

## Can set FD=FD_rep and FD_scaled=FD_rep_scaled to run most plots from FD section above

## Compare FD for all traits to FD_rep
FD = read.csv('./Data/raos_q_reproductive_traits_samples.csv')
rownames(FD) = FD$SampID
load(paste(null_fd_dir,'FD_scaled.RData', sep=''))

use_z1 = FD_scaled[,'z','withinplot','0.5','swsh_samp_c']
use_z2 = FD_rep_scaled[,'z','withinplot','0.5','swsh_samp_c']
plot(use_z1,use_z2, xlab='FD All Traits', ylab='FD Reproductive Traits')
abline(0,1)
abline(h=c(-1.96,1.96), v=c(-1.96,1.96), lty=2, col='grey20')

pdf('./Figures/FD all + rep traits z-score plot-scaled vs DBH.pdf', height=4.5, width=10)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))

use_z = FD_scaled[,'z','withinplot','0.5','swsh_samp_c']
use_data = samples[FD$SampID,]
FD_DBHmod_sm = lm(use_z[use_data$DBH<=6] ~ DBH, data=use_data[use_data$DBH<=6,])
FD_DBHmod_lg = lm(use_z[use_data$DBH>6] ~ DBH, data=use_data[use_data$DBH>6,])
plot(use_z~DBH, data = use_data, log='x',
	ylab='FD All Traits Z-score', xlab='Tree Diameter (cm)', las=1)
abline(h=c(-1.96,1.96), lty=2, col='grey20')
curve(predict(FD_DBHmod_sm, data.frame(DBH=x)), from=3, to=6, add=T, lwd=4)
curve(predict(FD_DBHmod_lg, data.frame(DBH=x)), from=7, to=max(use_data$DBH), add=T, lwd=4)

use_z = FD_rep_scaled[,'z','withinplot','0.5','swsh_samp_c']
use_data = samples[FD_rep$SampID,]
FD_DBHmod_sm = lm(use_z[use_data$DBH<=6] ~ DBH, data=use_data[use_data$DBH<=6,])
FD_DBHmod_lg = lm(use_z[use_data$DBH>6] ~ DBH, data=use_data[use_data$DBH>6,])
plot(use_z~DBH, data = use_data, log='x',
	ylab='FD Reproductive Traits Z-score', xlab='Tree Diameter (cm)', las=1)
abline(h=c(-1.96,1.96), lty=2, col='grey20')
curve(predict(FD_DBHmod_sm, data.frame(DBH=x)), from=3, to=6, add=T, lwd=4)
curve(predict(FD_DBHmod_lg, data.frame(DBH=x)), from=7, to=max(use_data$DBH), add=T, lwd=4)

dev.off()


## Models of FD_rep
# make sure to check that FD_rep and FD_rep_scaled are in same order (see code above)
# also make sure that variables in samples are correct type and levels
use_z = FD_rep_scaled[,'z','withinplot','0.5','swsh_samp_c']
use_data = samples[FD_rep$SampID,]

## All trees
fdmod = lm(use_z~FurrowDepth+Angle+pH+Density+WaterCapacity+as.numeric(Bryophytes)+as.numeric(Shedding)+DBH, data=use_data)
summary(fdmod) # FD significantly decreases with DBH even after other factors controlled for
anova(fdmod)

fdmod_plot = lmer(use_z~as.numeric(Bryophytes)+as.numeric(Shedding)+WaterCapacity+DBH + (DBH|SiteID/PlotID), data=use_data)
summary(fdmod_plot)
# NO Plot-level effects

## How does successively decreasing tree sizes affect coef estimate for DBH?
breakmods = sapply(4:max(use_data$DBH), function(thresh){
	use_samps = use_data$DBH<=thresh
	fdmod_sm = lm(use_z[use_samps]~FurrowDepth+Angle+pH+Density+WaterCapacity+as.numeric(Bryophytes)+as.numeric(Shedding)+DBH, data=use_data[use_samps,])
	
	DBH_coef = as.numeric(coef(fdmod_sm)['DBH'])
	DBH_se = coef(summary(fdmod_sm))['DBH','Std. Error']
	N = nrow(fdmod_sm$model)
	P = anova(fdmod_sm)$'Pr(>F)'[rownames(anova(fdmod_sm))=='DBH']
	
	c(Size=thresh, DBH_coef=DBH_coef, DBH_se=DBH_se, N=N, P=P)
})
breakmods = t(breakmods)
breakmods = as.data.frame(breakmods)

pdf('./Figures/FD rep vs DBH models with different size breakpoints.pdf', height=4, width=8)
plot(DBH_coef~Size, data=breakmods, las=1,ylim=c(-.05,.25))
arrows(breakmods$Size, breakmods$DBH_coef-(1.96*breakmods$DBH_se), breakmods$Size, breakmods$DBH_coef+(1.96*breakmods$DBH_se), length=0.02, angle=90, code=3)
points(DBH_coef~Size, data=breakmods[breakmods$P<0.05,], pch=16)
abline(h=0)
dev.off()

# Starting with large trees and including successively smaller trees
breakmods = sapply(3:70, function(thresh){
	use_samps = use_data$DBH>=thresh
	fdmod_sm = lm(use_z[use_samps]~FurrowDepth+Angle+pH+Density+WaterCapacity+as.numeric(Bryophytes)+as.numeric(Shedding)+DBH, data=use_data[use_samps,])
	
	DBH_coef = as.numeric(coef(fdmod_sm)['DBH'])
	DBH_se = coef(summary(fdmod_sm))['DBH','Std. Error']
	N = nrow(fdmod_sm$model)
	P = anova(fdmod_sm)$'Pr(>F)'[rownames(anova(fdmod_sm))=='DBH']
	
	c(Size=thresh, DBH_coef=DBH_coef, DBH_se=DBH_se, N=N, P=P)
})
breakmods = t(breakmods)
breakmods = as.data.frame(breakmods)

plot(DBH_coef~Size, data=breakmods, las=1,ylim=c(-.05,.25))
arrows(breakmods$Size, breakmods$DBH_coef-(1.96*breakmods$DBH_se), breakmods$Size, breakmods$DBH_coef+(1.96*breakmods$DBH_se), length=0.02, angle=90, code=3)
points(DBH_coef~Size, data=breakmods[breakmods$P<0.05,], pch=16)
abline(h=0)

## Looks like decrease in FD found for trees 20cm and greater.

