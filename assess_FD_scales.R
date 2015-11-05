## This script assess the variation in traits and diversity across spatial scales for the NC Lichen FD project

git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/NC-forest-lichens/'

# Load data & make data frames for analysis
source(paste(git_dir, 'load_data.R', sep=''))

# Define traits to be used in analysis
use_traits = traitdf[traitdf$TraitName %in% colnames(morphos),]
rownames(use_traits) = use_traits$TraitName

## Subset samples to piedmont and mountains only
samples_pm = subset(samples, SiteID!='Bladen')
samples_pm = droplevels(samples_pm)
lichens_pm = subset(lichens, SiteID!='Bladen')
lichens_pm = droplevels(lichens_pm)
plots_pm = subset(plot_data, PlotID!='Bladen1')

## Define variables used for variation partitioning in each section
# Choice of regional variables based on correlations between plot level climate variables:
# Strong correlation between Elevation, humidity and temperature, elevation least correlated with AP and CloudFreq_sd
sampvars = c('Angle','Bryophytes','Shedding','FurrowDepth','pH','Density','WaterCapacity')
treevars = c('DBH','Trans_tot_cor')
locvars = c(sampvars, treevars)
regvars = c('Elevation','CloudFreq_sd','AP','OpenPos','Soil_pH')


# Load repeatability statistics arrays generated later in script
load('repeatability_statistic_arrays.RData')


#################################################################
### Variation in morphotype diversity


## Summarize morphospecies and morphotype richness

# Calculate number of morphospecies per sample and tree (unique names)
samprich_s = sapply(rownames(samples), function(x) nrow(subset(lichens, SampID==x)))
treerich_s = sapply(rownames(tree_data), function(x) length(unique(subset(lichens, TreeID==x)$Name)))

# Calculate number of morphotypes per sample and tree (unique trait combos)
samprich_m = sapply(rownames(samples), function(x) sum(sampXmorph[x,]>0))
samprich_m[is.na(samprich_m)] = 0

treerich_m = sapply(rownames(tree_data), function(x){
	these_samps = rownames(subset(samples, TreeID==x))
	sum(colSums(sampXmorph[these_samps,])>0)
})
treerich_m[is.na(treerich_m)]=0

# Morphotype and morphospecies richness are tightly correlated at the sample and tree scale
plot(samprich_m~samprich_s); abline(0,1)
plot(treerich_m~treerich_s); abline(0,1)

# Make data frame
samprich = data.frame(SampID=names(samprich_m), Rich_M = samprich_m, Rich_S = samprich_s)
treerich = data.frame(TreeID=names(treerich_m), Rich_M_tree = treerich_m, Rich_S_tree = treerich_s)

# Calculate total and average abundance at the sample and tree level
sampabun = sapply(rownames(samples), function(x) sum(subset(lichens, SampID==x)$AbunCount))
sampabun = data.frame(SampID=names(sampabun), Tot_abun = sampabun)
treeabun = sapply(rownames(tree_data), function(x) sum(subset(lichens, TreeID==x)$AbunCount))
treeabun = data.frame(TreeID=names(treeabun), Tot_abun_tree = treeabun)

# Merge richness and abundance data
samp_derived = merge(samprich, sampabun)
tree_derived = merge(treerich, treeabun)

# Correct for increase in abundance with richness
samp_derived$Tot_abun_cor = samp_derived$Tot_abun - samp_derived$Rich_S + 1
samp_derived$Tot_abun_cor[samp_derived$Tot_abun==0] = 0
samp_derived$Avg_abun = samp_derived$Tot_abun/samp_derived$Rich_S
samp_derived$Avg_abun[is.na(samp_derived$Avg_abun)] = 0
plot(Tot_abun~Tot_abun_cor, data=samp_derived)
plot(Avg_abun~Tot_abun_cor, data=samp_derived)
plot(Tot_abun~Avg_abun, data=samp_derived)
plot(Avg_abun~Rich_S, data=samp_derived)

tree_derived$Tot_abun_tree_cor = tree_derived$Tot_abun_tree - tree_derived$Rich_S_tree + 1 # not quite right since some species can be counted twice per tree
tree_derived$Tot_abun_tree_cor[tree_derived$Tot_abun_tree==0] = 0
tree_derived$Avg_abun_tree = tree_derived$Tot_abun_tree/tree_derived$Rich_S_tree
tree_derived$Avg_abun_tree[is.na(tree_derived$Avg_abun_tree)] = 0

# Save data
write.csv(samp_derived, paste(derive_dir, 'sample_richness_abundance.csv', sep=''), row.names=F)
write.csv(tree_derived, paste(derive_dir, 'tree_richness_abundance.csv', sep=''), row.names=F)

# Read in data and merge with samples
samp_derived = read.csv(paste(derive_dir, 'sample_richness_abundance.csv', sep=''))
tree_derived = read.csv(paste(derive_dir, 'tree_richness_abundance.csv', sep=''))

samples = merge(samples, samp_derived, all.x=T)
samples = merge(samples, tree_derived, all.x=T)
rownames(samples) = samples$SampID

trees = merge(tree_data, tree_derived, all.x=T)
trees = merge(trees, plot_data, all.x=T)
rownames(trees) = trees$TreeID

# Re-subset piedmont and mountains samples
samples_pm = subset(samples, SiteID!='Bladen')
trees_pm = subset(trees, SiteID!='Bladen')
samples_pm = droplevels(samples_pm)
trees_pm = droplevels(trees_pm)

pdf('./Figures/Diversity across plots.pdf', height=6, width=4.5)
par(mfrow=c(2,1))
par(mar=c(0.5,4,5,1))
boxplot(samples$Rich_S~factor(samples$PlotID, levels=ordered_plots), 
	ylab='Species / Sample', axes=F)
axis(2, las=2); box()
abline(h=mean(samples$Rich_S), col=2)
par(mar=c(5,4,0.5,1))
boxplot(samples$Rich_S_tree~factor(samples$PlotID, levels=ordered_plots), 
	las=2, ylab='Species / Tree')
abline(h=mean(samples$Rich_S_tree), col=2)

par(mar=c(0.5,4,5,1))
boxplot(samples$Rich_M~factor(samples$PlotID, levels=ordered_plots), 
	ylab='Morphotypes / Sample', axes=F)
axis(2, las=2); box()
abline(h=mean(samples$Rich_M), col=2)
par(mar=c(5,4,0.5,1))
boxplot(samples$Rich_M_tree~factor(samples$PlotID, levels=ordered_plots), 
	las=2, ylab='Morphotypes / Tree')
abline(h=mean(samples$Rich_M_tree), col=2)

dev.off()

pdf('./Figures/Abundance across plots.pdf', height=6, width=4.5)
par(mfrow=c(2,1))
par(mar=c(0.5,5,5,1))
boxplot(samples_pm$Avg_abun~factor(samples_pm$PlotID, levels=ordered_plots), 
	ylab='Avg. Sample Abun.\n(Squares / Species)', axes=F)
axis(2, las=2); box()
abline(h=mean(samples_pm$Avg_abun), col=2)
par(mar=c(5,5,0.5,1))
boxplot(trees_pm$Avg_abun~factor(trees_pm$PlotID, levels=ordered_plots), 
	las=2, ylab='Avg. Tree Abun.\n(Squares / Species)')
abline(h=mean(trees_pm$Avg_abun), col=2)

par(mar=c(0.5,4,5,1))
boxplot(samples_pm$Tot_abun~factor(samples_pm$PlotID, levels=ordered_plots), 
	ylab='Total Squares / Sample', axes=F)
axis(2, las=2); box()
abline(h=mean(samples_pm$Tot_abun), col=2)
par(mar=c(5,4,0.5,1))
boxplot(trees_pm$Tot_abun~factor(trees_pm$PlotID, levels=ordered_plots), 
	las=2, ylab='Total Squares / Tree')
abline(h=mean(trees_pm$Tot_abun), col=2)

par(mar=c(0.5,5,5,1))
boxplot(samples_pm$Tot_abun_cor~factor(samples_pm$PlotID, levels=ordered_plots), 
	ylab='Corrected Abun.\n(Total Squares / Sample)', axes=F)
axis(2, las=2); box()
abline(h=mean(samples_pm$Tot_abun_cor), col=2)
par(mar=c(5,5,0.5,1))
boxplot(trees_pm$Tot_abun_cor~factor(trees_pm$PlotID, levels=ordered_plots), 
	las=2, ylab='Corrected Abun.\n(Total Squares / Tree)')
abline(h=mean(trees_pm$Tot_abun_cor), col=2)
dev.off()


## Manuscript Figure 3: boxplots of Tot_abun_cor, Rich_M, FD across plots
# Calculate elevation means
pairids = unique(plots_pm$PairID)
elev_means = with(plots_pm, tapply(Elevation, PairID, mean))

# Define plot order based on elevation and topographi position
pair_order = names(elev_means)[order(elev_means)]
plots_pm$PairID = factor(plots_pm$PairID, levels=pair_order)
plots_pm$TopoPos = factor(plots_pm$TopoPos, levels=c('sheltered','exposed'))
use_op = plots_pm[with(plots_pm, order(PairID, TopoPos)), 'PlotID']

# Define x-axis labels for plots
xlabels = sapply(1:length(elev_means), function(i){
	paste(pair_order[i], '\n', round(elev_means[pair_order[i]],0), 'm', sep='')
})

# Define colors
use_lcol = c('black','grey50')
use_col = c('grey80','white')
use_fact = plots_pm[use_op,'TopoPos']

# Calculate means across ecoregions and topo pos.
means = with(samples_pm, aggregate(samples_pm[,c('Avg_abun','Rich_M','rao_L2_s.5')], list(Ecoregion, TopoPos),
	FUN=function(x) mean(x, na.rm=T)))
colnames(means)[1:2] = c('Ecogregion','TopoPos')

# Calculate P-values for interaction using LR test
pvals = sapply(c('Avg_abun','Rich_M','rao_L2_s.5'), function(y){
	use_fam = ifelse(y=='Rich_M', 'poisson', 'gaussian')
	use_y = samples_pm[,y] + ifelse(use_fam=='gaussian', 0, 1)
	mod_inter = glm(use_y ~ Ecoregion*TopoPos, data=samples_pm, family=use_fam)
	mod_full = glm(use_y ~ Ecoregion + TopoPos, data=samples_pm, family=use_fam)
	LRtest = anova(mod_full, mod_inter, test='Chisq')
	LRtest[2,'Pr(>Chi)']
})
# Add these by hand later to control formatting.

# Define expression with plot labels
varnames = c('Average frequency','Morphotype richness','Functional diversity')
names(varnames) = c('Avg_abun','Rich_M','rao_L2_s.5')

# Loop through response variables
svg('./Figures/Community metrics across plots.svg', height = 6, width=4)
par(mfrow=c(3,1))
par(lend=1)
i = 1
for(y in c('Avg_abun','Rich_M','rao_L2_s.5')){
	# Set margins
	par(mar=c(1+i,5,3-i,1))
	
	# Add boxplots
	boxplot(samples_pm[,y] ~ factor(samples_pm$PlotID, levels=use_op), 
		ylab=varnames[y], axes=F, border=use_lcol[use_fact], col=use_col[use_fact])
	axis(2, las=2); box()
	
	# Add panel label
	mtext(c('A','B','C')[i], 3, 0, adj=-0.07)
	
	# Add group means
	segments(rep(c(8.5, 0.5), 2),  means[,y], rep(c(18.5, 8.5), 2),  means[,y], lwd=3, col=use_lcol[rep(1:2, each=2)])
	
	# Add x-axis labels if this is the last panel
	if(i==3){
		usr = par('usr')
		axis(1, at=seq(2, 18, 2)-.5, labels = pair_order, tick=F, line=1.5)
		axis(1, at=1:18, labels=round(plots_pm[use_op,'Elevation'], 0),
			tick=F, line=-.5, las=2)
		par(xpd = NA)
		text(usr[1]+.7, usr[3]-0.025, 'Elevation (m):', pos=2)
		text(usr[1]+.7, usr[3]-0.065,'Site:', pos=2)
		par(xpd = F)
	}
	
	i = i + 1
}
dev.off()


## FD vs Rich_M
par(mfrow=c(1,2))
for(i in c('Piedmont','Mountains')){
	use_data = subset(samples_pm, Ecoregion==i)
	plot(rao_L2_s.5 ~ Rich_M, data=use_data, xlim=c(0,17), ylim=c(0, 0.25),
		pch=21, bg = use_col[use_data$TopoPos], las=1)
}


## Morphotypes shared between plots

# Piedmont vs Mountains
pied_plots = rownames(subset(plot_data, Ecoregion=='Piedmont'))
mount_plots = rownames(subset(plot_data, Ecoregion=='Mountains'))
coast_plots = rownames(subset(plot_data, Ecoregion=='Coastal Plain'))

pied_morphs = sampXmorph[subset(samples, (rownames(samples) %in% rownames(sampXmorph))&(PlotID %in% pied_plots))$SampID,]
mount_morphs = sampXmorph[subset(samples, (rownames(samples) %in% rownames(sampXmorph))&(PlotID %in% mount_plots))$SampID,]
coast_morphs = sampXmorph[subset(samples, (rownames(samples) %in% rownames(sampXmorph))&(PlotID %in% coast_plots))$SampID,]

pied_morphs = pied_morphs[,colSums(pied_morphs)>0]
mount_morphs = mount_morphs[,colSums(mount_morphs)>0]
coast_morphs = coast_morphs[,colSums(coast_morphs)>0]

dim(pied_morphs)
dim(mount_morphs)
dim(coast_morphs)

length(intersect(colnames(pied_morphs), colnames(mount_morphs))) 
length(union(colnames(pied_morphs), colnames(mount_morphs)))
length(setdiff(colnames(pied_morphs), colnames(mount_morphs)))
length(setdiff(colnames(mount_morphs), colnames(pied_morphs)))

length(setdiff(colnames(coast_morphs), colnames(pied_morphs)))
length(intersect(colnames(pied_morphs), colnames(coast_morphs)))

length(setdiff(colnames(coast_morphs), colnames(mount_morphs)))
length(intersect(colnames(mount_morphs), colnames(coast_morphs)))

length(intersect(intersect(colnames(coast_morphs), colnames(pied_morphs)), colnames(mount_morphs)))


## Variation in sample-level morphotype/morphospecies richness across scales
library(lme4)
library(MuMIn)
library(rptR) # estimation of repeatbility: Nakagawa & Schielzeth 2010, Biological Reviews 85: 936-956

use_y = samples_pm$Rich_S + 1 # Go through all community metrics

hist(use_y)
hist(log(use_y))

# Check whether Poisson distribution is appropriate: 
# looks like it is for richness and avg abun, tot abundance is a bit over dispersed (Tot_abun_cor = 3.8)
vars = with(samples_pm, tapply(use_y, PlotID, var))
means = with(samples_pm, tapply(use_y, PlotID, mean))
plot(means, vars); abline(0,1)
abline(lm(vars~means-1), col=2)
lm(vars~means-1)

## Fit scale models for each community response metric- save all together in an array
#samples_pm$TreeID = factor(samples_pm$TreeID)
#samples_pm$PairID = factor(samples_pm$PairID)
#samples_pm$SiteID = factor(samples_pm$SiteID)
#samples_pm$Ecoregion = factor(samples_pm$Ecoregion)

## Will calculate repeatability on sqrt link scale so the variance is de-coupled from mean and can compare across comm metrics
samples_pm$SampID = factor(samples_pm$SampID)

comm_mets = c('Rich_S','Rich_M','Tot_abun_cor','Avg_abun','rao_L2_s.5')
scales = c('Tree','Plot','Site','Region')
scale_rpt = array(NA, dim=c(length(comm_mets), length(scales), 5), 
	dimnames=list(Metric=comm_mets, Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'))
)

for(i in comm_mets){
	use_y = samples_pm[,i] + ifelse(i=='rao_L2_s.5', 0, 1)

	if(i %in% c('Rich_S','Rich_M','Tot_abun_cor')){
		mod_tree = glmer(use_y ~ 1 + (1|PlotID) + (1|TreeID) + (1|SampID), data=samples_pm, family=poisson(link='sqrt'))
		mod_plot = glmer(use_y ~ 1 + (1|PlotID) + (1|PairID) + (1|SampID), data=samples_pm, family=poisson(link='sqrt'))
		mod_site = glmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID) + (1|SampID), data=samples_pm, family=poisson(link='sqrt'))
		mod_eco = glmer(use_y ~ 1 + (1|Ecoregion) + (1|SampID), data=samples_pm, family=poisson(link='sqrt'))
		obsID = 'SampID'
	}
	if(i %in% c('Avg_abun','rao_L2_s.5')){
		mod_tree = lmer(use_y ~ (1|PlotID) + (1|TreeID), data=samples_pm)
		mod_plot = lmer(use_y ~ 1 + (1|PlotID) + (1|PairID), data=samples_pm)
		mod_site = lmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID), data=samples_pm)
		mod_eco = lmer(use_y ~ 1 + (1|Ecoregion), data=samples_pm)
		obsID = 'Residual'
	}
	
	# Calculate adjusted repeatability (on link scale)
	scale_rpt[i,'Tree',1:4] = as.numeric(calc_rpt(mod_tree, focalID='TreeID', adjID='PlotID', obsID=obsID))
	scale_rpt[i,'Plot',1:4] = as.numeric(calc_rpt(mod_plot, focalID='PlotID', adjID='PairID', obsID=obsID))
	scale_rpt[i,'Site',1:4] = as.numeric(calc_rpt(mod_site, focalID='PairID', adjID='Ecoregion', obsID=obsID))
	scale_rpt[i,'Ecoregion','R.focal'] = as.numeric(calc_rpt(mod_eco, focalID='Ecoregion', obsID='SampID'))

	scale_rpt[i,,'R.residual'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}

scale_rpt[,'Ecoregion',c('R.focal.adj','R.adj','R.adj.focal')] = scale_rpt[,'Ecoregion','R.focal']

## Adjust for fixed effects from env predictors
env_data = samples_pm[,c(locvars, regvars)]
use_data = cbind(env_data, samples_pm[,c('SampID','TreeID','PlotID','PairID','Ecoregion')])

# Drop observations with missing env variables
use_data = na.omit(use_data)

# Model ordered data as integers and re-scale predictors to unit variance
for(i in 1:ncol(env_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
	use_data[,i] = scale(use_data[,i], center=T, scale=T)
}

scale_rpt_res = array(NA, dim=c(length(comm_mets), length(scales), 5, 2), 
	dimnames=list(Metric=comm_mets, Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'), Control = c('Scale','Env'))
)

for(i in comm_mets){
	use_y = samples_pm[rownames(use_data),i] + ifelse(i=='rao_L2_s.5', 0, 1)

	if(i %in% c('Rich_S','Rich_M','Tot_abun_cor')){
		mod_tree = glmer(use_y ~ 1 + (1|PlotID) + (1|TreeID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_plot = glmer(use_y ~ 1 + (1|PlotID) + (1|PairID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_site = glmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_eco = glmer(use_y ~ 1 + (1|Ecoregion) + (1|SampID), data=use_data, family=poisson(link='sqrt'))

		mod_tree_env = glmer(use_y ~ Angle + Bryophytes + Shedding + FurrowDepth + pH + Density + WaterCapacity + DBH + Trans_tot_cor + 
			(1|PlotID) + (1|TreeID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_plot_env = glmer(use_y ~ OpenPos + Soil_pH +
			(1|PlotID) + (1|PairID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_site_env = glmer(use_y ~ Elevation + AP + CloudFreq_sd +
			(1|Ecoregion) + (1|PairID) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		mod_eco_env = glmer(use_y ~ Elevation + AP + CloudFreq_sd +
			(1|Ecoregion) + (1|SampID), data=use_data, family=poisson(link='sqrt'))
		obsID = 'SampID'
	}
	if(i %in% c('Avg_abun','rao_L2_s.5')){
		mod_tree = lmer(use_y ~ (1|PlotID) + (1|TreeID), data=use_data)
		mod_plot = lmer(use_y ~ 1 + (1|PlotID) + (1|PairID), data=use_data)
		mod_site = lmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID), data=use_data)
		mod_eco = lmer(use_y ~ 1 + (1|Ecoregion), data=use_data)

		mod_tree_env = lmer(use_y ~ Angle + Bryophytes + Shedding + FurrowDepth + pH + Density + WaterCapacity + DBH + Trans_tot_cor + 
			(1|PlotID) + (1|TreeID), data=use_data)
		mod_plot_env = lmer(use_y ~ OpenPos + Soil_pH +
			(1|PlotID) + (1|PairID), data=use_data)
		mod_site_env = lmer(use_y ~  Elevation + AP + CloudFreq_sd +
			(1|Ecoregion) + (1|PairID), data=use_data)
		mod_eco_env = lmer(use_y ~ Elevation + AP + CloudFreq_sd +
			(1|Ecoregion), data=use_data)
		obsID = 'Residual'
	}
	
	# Calculate adjusted repeatability (on link scale)
	scale_rpt_res[i,'Tree',1:4,] = sapply(list(mod_tree, mod_tree_env), function(mod) as.numeric(calc_rpt(mod, focalID='TreeID', adjID=c('Fixed','PlotID'), obsID=obsID)))
	scale_rpt_res[i,'Plot',1:4,] = sapply(list(mod_plot, mod_plot_env), function(mod) as.numeric(calc_rpt(mod, focalID='PlotID', adjID=c('Fixed','PairID'), obsID=obsID)))
	scale_rpt_res[i,'Site',1:4,] = sapply(list(mod_site, mod_site_env), function(mod) as.numeric(calc_rpt(mod, focalID='PairID', adjID=c('Fixed','Ecoregion'), obsID=obsID)))
	scale_rpt_res[i,'Region',1:4,] = sapply(list(mod_eco, mod_eco_env), function(mod) as.numeric(calc_rpt(mod, focalID='Ecoregion', adjID='Fixed', obsID=obsID)))

	scale_rpt_res[i,,'R.residual','Scale'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))
	scale_rpt_res[i,,'R.residual','Env'] = as.numeric(sapply(c(mod_tree_env, mod_plot_env, mod_site_env, mod_eco_env), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}


# Format into table
library(reshape)
rpt_melt = melt(scale_rpt)
rpt_df = cast(rpt_melt, Metric+Model ~ Statistic, subset=Statistic %in% c('R.focal.adj','R.adj.focal','R.residual'))

rpt_df$Metric = factor(rpt_df$Metric, levels=comm_mets)
rpt_df$Model = factor(rpt_df$Model, levels=scales)
rpt_df = rpt_df[order(rpt_df$Metric, rpt_df$Model),c(1,2,5,4,3)]

write.csv(rpt_df, './Figures/scale_mods_repeatability.csv', row.names=F)

# Plot R2 components
components = scale_rpt['Rich_S',,'R.focal.adj']

pdf('./Figures/variance across scales nested sample morphospecies richness.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.3), ylab=expression(R^2), main='Sample Morphospecies Richness')
dev.off()

# Plot components after adjusting for fixed effects of env predictors
svg('./Figures/repeatability of community metrics across scales env residuals.svg', height=4, width=5)
par(mar=c(3, 4.5, 1, 1))
barplot(t(scale_rpt[,,'R.focal.adj']), beside=T, legend.text = colnames(scale_rpt), las=1, ylim=c(0,1),
	args.legend=list(x='topright', bty='n'), ylab='Variation Explained',
	names.arg=expression(R[S],R[M],N[TOT],N[AVG],FD), space=c(0.2,0.8))
barplot(t(scale_rpt_res[,,'R.focal.adj']), beside=T, add=T, density=20, col='black',
	axes=F, names.arg=expression(R[S],R[M],N[TOT],N[AVG],FD), space=c(0.2,0.8))
dev.off()




## Do models one at a time to make plots:
use_y = samples_pm$Rich_S + 1
mod_all = glmer(use_y ~ (1|SampID) + (1|TreeID) + (1|PlotID) + (1|SiteID) + (1|Ecoregion), data=samples_pm, family=poisson(link='sqrt'))
use_y = samples_pm$Avg_abun
mod_all = lmer(use_y ~ (1|TreeID) + (1|PlotID) + (1|SiteID) + (1|Ecoregion), data=samples_pm)

# Effects of ecoregion and topographic position
mod_inter = glmer(use_y ~ 1 + Ecoregion*TopoPos + (1|SiteID) + (1|TreeID) , data = samples_pm, family='poisson')
mod_full = glmer(use_y ~ 1 + Ecoregion+TopoPos + (1|SiteID) + (1|TreeID) , data = samples_pm, family='poisson')
mod_null = glmer(use_y ~ 1 + (1|SiteID) + (1|TreeID), data=samples_pm, family='poisson')

# Avg abun models- log link not supported by AIC
mod_inter = lmer(use_y ~ 1 + Ecoregion*TopoPos + (1|SiteID) + (1|TreeID) , data = samples_pm)
mod_full = lmer(use_y ~ 1 + Ecoregion+TopoPos + (1|SiteID) + (1|TreeID) , data = samples_pm)
mod_null = lmer(use_y ~ 1 + (1|SiteID) + (1|TreeID), data=samples_pm)

# Calculate variance components from model without interaction
# These are all on the scale of the link function (log)

# Calculate variance of fixed effects
sigma.eco = var(as.numeric(as.vector(fixef(mod_full)[1:2]) %*% t(model.matrix(mod_full)[,1:2])))
sigma.plot = var(as.numeric(as.vector(fixef(mod_full)[c(1,3)]) %*% t(model.matrix(mod_full)[,c(1,3)])))
sigma.fixed = calc_fixedvar(mod_full)
	
# Extract variance of random effects
sigma.tree = data.frame(VarCorr(mod_full))$vcov[1]
sigma.site = data.frame(VarCorr(mod_full))$vcov[2]

# Calculate total variance and R2 (based on Nakagawa and Schielzeth 2013)
# NO LONGER WORKS: NEED TO RECALCULATE 9/29/2015
#sigma.tot = calc_totvar(mod_full)

# Conditional R2 (variance explained by all factors)
(sigma.fixed + sigma.tree + sigma.site) / sigma.tot

# Plot variance components
components = c(sigma.tree, sigma.plot, sigma.site, sigma.eco)/sigma.tot
names(components) = c('Tree','Plot','Site','Ecoregion')

# Change names, axis labels, and ylims for each response variable
pdf('./Figures/variance across scales topo-eco sample Tot_abun_cor.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.7), ylab=expression(paste('Marginal ', R^2)), main='Sample Total Abundance')
dev.off()

# Plot effects of ecoregion and plot in the interaction model
counts = with(samples_pm, aggregate(factor(use_y-1, levels=0:max(use_y)), list(TopoPos, Ecoregion), FUN=table))
bars = counts[,-(1:2)] / rowSums(counts[,-(1:2)])

use_col = c('grey50', 'white')

pdf('./Figures/sample morphospecies richness across ecoregions and topo.pdf', height=4, width=4)
par(mar=c(4,4,1,1))

# Set up plot
plot(c(-1, 1), range(use_y-1), type='n', las=1, axes=F, xlim=c(-2,2), xlab='', ylab='Sample Morphotype Richness')
axis(1, at=c(-1,1), labels=c('Mountains','Piedmont'))
axis(2, at=seq(0,max(use_y-1), 2), las=1)
box()

# Make rectangles
rect(-1, (0:(max(use_y)-1))-.5, -1-bars[1,], (1:(max(use_y)))-.5, col=use_col[1])
rect(-1, (0:(max(use_y)-1))-.5, -1+bars[2,], (1:(max(use_y)))-.5, col=use_col[2])
rect(1, (0:(max(use_y)-1))-.5, 1-bars[3,], (1:(max(use_y)))-.5, col=use_col[1])
rect(1, (0:(max(use_y)-1))-.5, 1+bars[4,], (1:(max(use_y)))-.5, col=use_col[2])

# Add estimated effects from model with interaction
ests =  exp(fixef(mod_inter))

modmat = expand.grid(TopoPos=c('sheltered','exposed'),Ecoregion=c('Mountains','Piedmont'))
means = predict(mod_inter, modmat, re.form=~0, type='response') - 1

points(rep(c(-1,1), each=2) + rep(c(-.3, .3), 2), means, pch=3, bg=use_col, cex=2, lwd=2, lend=1)

dev.off()

## Boxplot of sample-level abundance differences across Ecoregions
# since TopoPos non-significant
summary(mod_inter)
summary(mod_full)

anova(mod_full, update(mod_full, .~.-TopoPos), test='Chisq') # test TopoPos
anova(mod_full, update(mod_full, .~.-Ecoregion), test='Chisq') 

boxplot(Tot_abun_cor~Ecoregion, data=samples_pm)


## Variation in tree-level morphotype/morphospecies richness across scales
use_y = trees_pm$Avg_abun_tree + 1 # Rich_M_tree, Rich_S_tree, Tot_abun_tree_cor, Avg_abun_tree

hist(use_y)
hist(log(use_y))

# Check whether Poisson distribution is appropriate: 
# looks like it is a bit overdispersed for richness and fairly overdispersed for total abun (6.76)
vars = with(trees_pm, tapply(use_y, PlotID, var))
means = with(trees_pm, tapply(use_y, PlotID, mean))
plot(means, vars); abline(0,1)
abline(lm(vars~means-1), col=2)
lm(vars~means-1)

mod_inter = glmer(use_y ~ 1 + Ecoregion*TopoPos + (1|SiteID), data = trees_pm, family='poisson')
mod_null = glmer(use_y ~ 1 + (1|SiteID), data=trees_pm, family='poisson')
mod_full = glmer(use_y ~ 1 + Ecoregion+TopoPos + (1|SiteID), data = trees_pm, family='poisson')

# Avg abun models
mod_inter = lmer(use_y ~ 1 + Ecoregion*TopoPos + (1|SiteID), data = trees_pm)
mod_null = lmer(use_y ~ 1 + (1|SiteID), data=trees_pm)
mod_full = lmer(use_y ~ 1 + Ecoregion+TopoPos + (1|SiteID), data = trees_pm)

# Calculate variance components from model without interaction
# These are all on the scale of the link function (log)

# Calculate variance of fixed effects
sigma.eco = var(as.numeric(as.vector(fixef(mod_full)[1:2]) %*% t(model.matrix(mod_full)[,1:2])))
sigma.plot = var(as.numeric(as.vector(fixef(mod_full)[c(1,3)]) %*% t(model.matrix(mod_full)[,c(1,3)])))
sigma.fixed = calc_fixedvar(mod_full)
	
# Extract variance of random effects
sigma.site = data.frame(VarCorr(mod_full))$vcov[1]

# Calculate total variance and R2 (based on Nakagawa and Schielzeth 2013)
sigma.tot = calc_totvar(mod_full)

# Conditional R2 (variance explained by all factors)
(sigma.fixed + sigma.site) / sigma.tot

# Plot variance components
components = c(sigma.plot, sigma.site, sigma.eco)/sigma.tot
names(components) = c('Plot','Site','Ecoregion')

pdf('./Figures/variance across scales topo-eco tree morphotype richness.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.4), ylab=expression(paste('Marginal ', R^2)), main='Tree Morphotype Richness')
dev.off()

# Plot effects of ecoregion and plot in the interaction model
counts = with(trees_pm, aggregate(factor(use_y-1, levels=0:max(use_y)), list(TopoPos, Ecoregion), FUN=table))
bars = counts[,-(1:2)] / rowSums(counts[,-(1:2)])

use_col = c('grey50', 'white')

pdf('./Figures/tree morphotype richness across ecoregions and topo.pdf', height=4, width=4)
par(mar=c(4,4,1,1))

# Set up plot
plot(c(-1, 1), range(use_y-1), type='n', las=1, axes=F, xlim=c(-2,2), xlab='', ylab='Tree Morphotype Richness')
axis(1, at=c(-1,1), labels=c('Mountains','Piedmont'))
axis(2, at=seq(0,max(use_y-1), 2), las=1)
box()

# Make rectangles
rect(-1, (0:(max(use_y)-1))-.5, -1-bars[1,], (1:(max(use_y)))-.5, col=use_col[1])
rect(-1, (0:(max(use_y)-1))-.5, -1+bars[2,], (1:(max(use_y)))-.5, col=use_col[2])
rect(1, (0:(max(use_y)-1))-.5, 1-bars[3,], (1:(max(use_y)))-.5, col=use_col[1])
rect(1, (0:(max(use_y)-1))-.5, 1+bars[4,], (1:(max(use_y)))-.5, col=use_col[2])

# Add estimated effects from model with interaction
ests =  exp(fixef(mod_inter))

modmat = expand.grid(TopoPos=c('sheltered','exposed'),Ecoregion=c('Mountains','Piedmont'))
means = predict(mod_inter, modmat, re.form=~0, type='response') - 1

points(rep(c(-1,1), each=2) + rep(c(-.3, .3), 2), means, pch=3, bg=use_col, cex=2, lwd=2, lend=1)

dev.off()

## Effects of scale itself
mod_plot = glmer(use_y ~ 1 + (1|SiteID/PlotID), data=trees_pm, family='poisson')
mod_site = glmer(use_y ~ 1 + (1|SiteID) + Ecoregion, data=trees_pm, family='poisson')
mod_eco = glm(use_y ~ Ecoregion, data=trees_pm, family='poisson')

# Avg abun models
mod_plot = lmer(use_y ~ 1 + (1|SiteID/PlotID), data=trees_pm)
mod_site = lmer(use_y ~ 1 + (1|SiteID) + Ecoregion, data=trees_pm)
mod_eco = glm(use_y ~ Ecoregion, data=trees_pm, family='gaussian')

# Plot R2 components
components = sapply(list(mod_plot, mod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
components = c(components, attr(r.squaredLR(mod_eco), 'adj.r.squared'))
names(components) = c('Plot','Site','Ecoregion')

pdf('./Figures/variance across scales nested tree morphotype richness.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.5), ylab=expression(R^2), main='Sample Morphotype Richness')
dev.off()


#############################
### Variance Partitioning ###

# Sample-scale response
use_y = samples_pm$rao_L2_s.5# + 1# Rich_M, Rich_S, Tot_abun_cor, Avg_abun, rao_L2_s.5  # don't use + 1 on FD
locvars = c(sampvars, treevars)
use_data = samples_pm[,c(locvars, regvars)]

# Tree-scale response
#use_y = trees_pm$Avg_abun_tree #+ 1 # Rich_S_tree, Rich_M_tree, Tot_abun_tree_cor, Avg_abun_tree
#locvars = treevars
#use_data = trees_pm[,c(locvars, regvars)]

# Remove observations with missing values
missing = rowSums(is.na(use_data)|is.na(use_y))>0
use_y = use_y[!missing]
use_data = use_data[!missing,]

# Model ordered data as integers
for(i in 1:ncol(use_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
}

# For count models
locmod = glm(use_y~., data=use_data[,locvars], family='poisson')
regmod = glm(use_y~., data=use_data[,regvars], family='poisson')
fullmod = glm(use_y~., data=use_data, family='poisson')
R2s = sapply(list(local = locmod, regional = regmod, both = fullmod), function(x) attr(r.squaredLR(x),'adj.r.squared'))


# For average abundance and FD
locmod = lm(use_y~., data=use_data[,locvars])
regmod = lm(use_y~., data=use_data[,regvars])
fullmod = lm(use_y~., data=use_data)
R2s = sapply(list(local = locmod, regional = regmod, both = fullmod), function(x) summary(x)$adj.r.squared)


# Calculate R2 using Nakagawa and Schielzeth 2013 mariginal R2 (even though no random effects)
# Estimates unreliable: DECIDED TO USE LR INSTEAD 9/27/2015
#R2s = sapply(list(local = locmod, regional = regmod, both = fullmod), function(x) r.squaredGLMM(x)['R2m'])

# Save results
varpart_richS = partvar2(R2s)
varpart_richM = partvar2(R2s)
varpart_totabun = partvar2(R2s)
varpart_avgabun = partvar2(R2s)
varpart_fd = partvar2(R2s)

# No longer using
#varpart_richS_tree = partvar2(R2s)
#varpart_richM_tree = partvar2(R2s)
#varpart_totabun_tree = partvar2(R2s)
#varpart_avgabun_tree = partvar2(R2s)

## Fit full models and save coefficients
metrics = c('Rich_S','Rich_M','Tot_abun_cor','Avg_abun','rao_L2_s.5')
envvars = c(locvars, regvars)

use_data = samples_pm[,envvars]

# Remove observations with missing values
missing = rowSums(is.na(use_data))>0
use_data = use_data[!missing,]

# Model ordered data as integers
# Center and scale predictors to get standardized effects
for(i in 1:ncol(use_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
	use_data[,i] = as.numeric(use_data[,i])
}
use_data = scale(use_data, center=T, scale=T)
use_data = data.frame(use_data)

varpart_mods = array(NA, dim=c(length(metrics), length(envvars), 6), 
	dimnames=list(Response=metrics, Predictor=envvars, Statistic=c('Est','SE','Dev.drop','P.drop','Dev.add','P.add')))

for(y in metrics){
	use_y = samples_pm[rownames(use_data), y] + ifelse(y=='rao_L2_s.5', 0, 1)

	if(y %in% c('Rich_S','Rich_M','Tot_abun_cor')){
		full_mod = glm(use_y ~ ., data=use_data, family='poisson')
	}

	if(y %in% c('Avg_abun','rao_L2_s.5')){
		full_mod = glm(use_y ~ ., data=use_data, family='gaussian')
	}
	
	varpart_mods[y,,'Est'] = coef(full_mod)[envvars]
	varpart_mods[y,,'SE'] = coef(summary(full_mod))[envvars,'Std. Error']
	
	for(x in envvars){
		drop_formula = paste('.~.-', x, sep='')
		add_formula = paste('use_y~', x, sep='')
		
		drop_mod = do.call('update', list(full_mod, as.formula(drop_formula)))
		add_mod = do.call('update', list(full_mod, as.formula(add_formula)))

		lrtest_drop = anova(drop_mod, full_mod, test='Chisq')
		lrtest_add = anova(add_mod, test='Chisq')

		varpart_mods[y,x,c('Dev.drop','P.drop')] = as.numeric(lrtest_drop[2,c('Deviance','Pr(>Chi)')])
		varpart_mods[y,x,c('Dev.add','P.add')] = as.numeric(lrtest_add[2,c('Deviance','Pr(>Chi)')])
	}

}

varpart_melt = melt(varpart_mods)
varpart_df = cast(varpart_melt, Response+Predictor~Statistic)
varpart_df$Response = factor(varpart_df$Response, levels=c('Rich_S','Rich_M','Tot_abun_cor','Avg_abun','rao_L2_s.5'))
varpart_df$Scale = with(varpart_df, ifelse(Predictor %in%  sampvars, 'Sample',
	ifelse(Predictor %in% treevars, 'Tree', 
	ifelse(Predictor %in% c('OpenPos','Soil_pH'), 'Plot', 'Site'))))
varpart_df = varpart_df[with(varpart_df, order(Response, 1-abs(Est))),]
varpart_df = varpart_df[,c('Response','Predictor','Scale', 'Est','SE','Dev.drop','P.drop','Dev.add','P.add')]

varnames = c('Surface angle','Bryophyte abundance','Bark stability','Bark furrow depth',
	'Bark pH', 'Bark density','Bark WHC','Tree DBH','Light transmittance','Elevation','Cloud seasonality',
	'Annual precipitation','Topographic openness','Soil pH')
names(varnames) = envvars
varpart_df$Predictor = varnames[varpart_df$Predictor]

write.csv(varpart_df, './Figures/community metric env models.csv')





## Plot scales over which important env variables vary
# Rich_M, FD: Bryophytes, Shedding, WaterCapacity, DBH, Trans_tot_cor
# Tot_abun_cor: also pH and Density
# Avg_abun: only Bryophytes, DBH, maybe WaterCapacity
use_env = c('Bryophytes','Shedding','WaterCapacity','pH','FurrowDepth','Trans_tot_cor', 'DBH')

use_data = samples_pm[,c(use_env, 'SampID','TreeID','PlotID','PairID','Ecoregion')]

# Model ordered data as integers
for(i in 1:length(use_env)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
}

scale_rpt_env = array(NA, dim=c(length(use_env), length(scales), 5), 
	dimnames=list(Predictor=use_env, Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'))
)

for(i in use_env){
	use_y = use_data[,i]

 	if(!(i %in% c('DBH','Trans_tot_cor'))) mod_tree = lmer(use_y ~ (1|PlotID) + (1|TreeID), data=use_data) # 	
	mod_plot = lmer(use_y ~ 1 + (1|PlotID) + (1|PairID), data=use_data)
	mod_site = lmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID), data=use_data)
	mod_eco = lmer(use_y ~ 1 + (1|Ecoregion), data=use_data)
	obsID = 'Residual'
	
	# Calculate adjusted repeatability (on link scale)
	if(!(i %in% c('DBH','Trans_tot_cor'))) scale_rpt_env[i,'Tree',1:4] = as.numeric(calc_rpt(mod_tree, focalID='TreeID', adjID='PlotID', obsID=obsID)) # 
	scale_rpt_env[i,'Plot',1:4] = as.numeric(calc_rpt(mod_plot, focalID='PlotID', adjID='PairID', obsID=obsID))
	scale_rpt_env[i,'Site',1:4] = as.numeric(calc_rpt(mod_site, focalID='PairID', adjID='Ecoregion', obsID=obsID))
	scale_rpt_env[i,'Region','R.focal'] = as.numeric(calc_rpt(mod_eco, focalID='Ecoregion', obsID='SampID'))

	scale_rpt_env[i,,'R.residual'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}

scale_rpt_env[,'Region',c('R.focal.adj','R.adj','R.adj.focal')] = scale_rpt_env[,'Region','R.focal']

## Manuscript Figure
svg('./Figures/repeatability of environment across scales.svg', height=4, width=5)
par(mar=c(3, 4.5, 1, 1))
barplot(t(scale_rpt_env[1:6,,'R.focal.adj']), beside=T, legend.text = colnames(scale_rpt_env), las=1, ylim=c(0,1),
	args.legend=list(x='topright', bty='n'), ylab='Variation Explained',
	names.arg=expression(Bryophytes,Stability,WHC,pH,Furrows,Light), space=c(0.2,0.8))
dev.off()







#################################################################
### Variation in morphotypes (equivalent to variation in species composition)
library(ape)
library(vegan)

dist_L2mat = read.csv(paste(derive_dir, 'morpho_distance_matrix_sfact0.5.csv'), row.names=1)

# Subset morphotypes and community matrix to just morphotypes in Piedmont and Mountains
pm_morphs = unique(subset(lichens, SiteID!='Bladen')$MorphID)
morph_dmat = dist_L2mat[pm_morphs, pm_morphs]
morphos_pm = morphos[pm_morphs,]
use_samps = as.character(samples_pm$SampID)
comm = sampXmorph[use_samps[use_samps %in% rownames(sampXmorph)], pm_morphs]

# PCOA of trait distance matrix
mds = cmdscale(morph_dmat, k=nrow(morph_dmat)-2, eig=T, add=T) # Uses Calliez correction for neg eigenvalues

# Eigenvalues of axes
eigs = mds$eig[1:(nrow(morph_dmat)-2)]
plot(1:length(eigs), eigs/sum(eigs))

# Define morphotype vectors for biplots
pcoa_vecs = mds$points

sum(rownames(pcoa_vecs)!=colnames(comm))

# Highlight a couple traits
colorby = morphos_pm[,'Attachment']
mycols = colorRampPalette(c('black','cornflowerblue'))(length(levels(colorby)))
pch1by = morphos_pm[,'Asexual']
pch2by = morphos_pm[,'Photobiont']
mypch=matrix(c(0,15,1,16),2,2)
rownames(mypch) = levels(factor(pch1by))
colnames(mypch) = levels(factor(pch2by))
use_pch = apply(cbind(pch1by,pch2by), 1, function(i) mypch[i[1],i[2]])

pdf('./Figures/Morphos PCoA colored by attachment pch asexual photobiont.pdf', height=5, width=5)
par(mar=c(4,4,0.5,0.5))
plot(pcoa_vecs[,c(1,2)], col=mycols[colorby], pch = use_pch, las=1, xlab='PCoA1', ylab='PCoA2')
dev.off()

# Make plot just with morphotypes from two samples
subset(tree_data, TreeTaxonID=='Acer_rubrum'&DBH>20)
S1 = 'High1-35-S1' # 32cm red maple
S2 = 'Eno2-18-S1' # 31cm red mapl

pdf('./Figures/Morphos PCoA High1 Eno1 red maples.pdf', height=5, width=5)
par(mar=c(4,4,0.5,0.5))
plot(pcoa_vecs[,c(1,2)], type='n', las=1, xlab='PCoA1', ylab='PCoA2')
use_Ms = which(colSums(comm[S1,])>0)
points(pcoa_vecs[use_Ms,c(1,2)], col=mycols[colorby[use_Ms]], pch = use_pch[use_Ms])
use_Ms = which(colSums(comm[S2,])>0)
points(pcoa_vecs[use_Ms,c(1,2)], col=mycols[colorby[use_Ms]], pch = use_pch[use_Ms])
dev.off()

subset(tree_data, TreeTaxonID=='Acer_rubrum'&DBH>20)

pdf('./Figures/Morphos PCoA colored by traits.pdf', height=5, width=6.5)
layout(matrix(c(1,3,2,4,5,5),2,3), widths=c(.4,.4,.2))
use_traits = traitdf[traitdf$type!='numeric','TraitName']
use_traits = use_traits[use_traits %in% colnames(morphos)]
for(i in use_traits){
	
	colorby = morphos_pm[,i]
	mycols = colorRampPalette(c('black','cornflowerblue'))(length(levels(colorby)))
	
	par(mar=c(4,4,1,1))
	plot(pcoa_vecs[,c(1,2)], col=mycols[colorby], pch=mypch[4])
	plot(pcoa_vecs[,c(3,4)], col=mycols[colorby], pch=mypch[4])
	plot(pcoa_vecs[,c(5,6)], col=mycols[colorby], pch=mypch[4])
	plot(pcoa_vecs[,c(7,8)], col=mycols[colorby], pch=mypch[4])
	par(mar=c(0,0,0,0))
	plot.new()
	legend('center', levels(colorby), col=mycols, pch=16, bty='n', title=i)
}
dev.off()

## Partition variance across scales using RDA

# Analyze Hellinger transformed community data to reduce effect of rare species
comm_hel = sqrt(comm/rowSums(comm))

# Calculate location of samples in trait space
sampXpcoa = as.matrix(comm_hel) %*% pcoa_vecs 


# Unconstrained ordination
rda_null = rda(sampXpcoa)

# Make data for conditioning on environmental variables
env_data = samples_pm[rownames(comm),c(locvars, regvars)]
use_data = cbind(env_data, samples_pm[rownames(comm),c('SampID','TreeID','PlotID','PairID','Ecoregion')])

# Remove samples that are missing env vars (41 samples)- will be slightly different from other rda
keep_rows = rowSums(is.na(env_data))==0
use_data = use_data[keep_rows,]

# Model ordered data as integers and re-scale predictors to unit variance
for(i in 1:ncol(env_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
	use_data[,i] = scale(use_data[,i], center=T, scale=T)
}

# Compute RDA with scales as predictors
rda_full = rda(sampXpcoa[keep_rows,] ~ TreeID + TopoPos + PairID + Ecoregion, data=use_data)

rda_tree = rda(sampXpcoa[keep_rows,] ~ TreeID + Condition(PlotID), data=use_data)
rda_plot = rda(sampXpcoa[keep_rows,] ~ PlotID + Condition(PairID), data=use_data)
rda_site = rda(sampXpcoa[keep_rows,] ~ PairID + Condition(Ecoregion), data=use_data)
rda_eco = rda(sampXpcoa[keep_rows,] ~ Ecoregion, data=use_data)

# Compute RDA with scales and env variables as predictors
rda_tree_env = rda(sampXpcoa[keep_rows,] ~ TreeID + Condition(PlotID) + Condition(Angle) + Condition(Bryophytes) + Condition(Shedding) + Condition(FurrowDepth) + Condition(pH) +
	Condition(Density) + Condition(WaterCapacity) + Condition(DBH) + Condition(Trans_tot_cor), data=use_data)
rda_plot_env = rda(sampXpcoa[keep_rows,] ~ PlotID + Condition(PairID) + Condition(OpenPos) + Condition(Soil_pH), data=use_data)
rda_site_env = rda(sampXpcoa[keep_rows,] ~ PairID + Condition(Ecoregion) + Condition(Elevation) + Condition(AP) + Condition(CloudFreq_sd), data=use_data)
rda_eco_env = rda(sampXpcoa[keep_rows,] ~ Ecoregion + Condition(Elevation) + Condition(AP) + Condition(CloudFreq_sd), data=use_data)


#save(rda_tree, rda_plot, rda_site, rda_eco, rda_tree_env, rda_plot_env, rda_site_env, rda_eco_env,file='RDA_models.RData')
load('RDA_models.RData')


# Plot variance components: type into manuscript table
# semi-partial R2
#components = sapply(list(rda_tree, rda_plot, rda_site, rda_eco), function(x) RsquareAdj(x)$adj.r.squared) # Variance explained is adjusted R2 from constrained component
# partial R2
calc_pR2 = function(mod){
	SSfit = mod$CCA$tot.chi
	SSres = mod$tot.chi - ifelse(is.null(mod$pCCA), 0, mod$pCCA$tot.chi)
	SSfit/SSres
}
components = sapply(list(rda_tree, rda_plot, rda_site, rda_eco), calc_pR2)
names(components) = c('Tree','Plot','Site','Ecoregion')
components_env = sapply(list(rda_tree_env, rda_plot_env, rda_site_env, rda_eco_env), calc_pR2)
names(components_env) = c('Tree','Plot','Site','Ecoregion')

pdf('./Figures/variance across scales sample morphotype composition.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.4), ylab=expression(paste('Adjusted Semi-Partial ', R^2)), main='Sample Functional Composition')
dev.off()


## Manuscript Figure: repeatability across scales:
# Note: all models omit 62 / 1440 observation for which some environmental variable(s) was missing.
bp_dat = rbind(scale_rpt_res[,,'R.focal.adj','Scale'], components)
bp_dat_res = rbind(scale_rpt_res[,,'R.focal.adj','Env'], components_env)

bp_dat = cbind(1-bp_dat[,'Tree'], bp_dat)
colnames(bp_dat)[1] = 'Sample'
bp_dat_res = cbind(rep(0,nrow(bp_dat_res)), bp_dat_res)
colnames(bp_dat_res)[1] = 'Sample'

svg('./Figures/repeatability of community metrics across scales.svg', height=4, width=5)
par(mar=c(3, 4.5, 1, 1))
barplot(t(as.matrix(bp_dat)), beside=T, legend.text = colnames(bp_dat), las=1, ylim=c(0,1),
	args.legend=list(x='topright', bty='n'), ylab='Variation Explained',
	names.arg=expression(R[S],R[M],N[TOT],N[AVG],FD,FC), space=c(0.2,0.8))
dev.off()

use_col = c('black', colorRampPalette(c('grey65','white'))(4))

plot_mets = c(1:2,4:6)
svg('./Figures/repeatability of community metrics across scales with res.svg', height=4, width=6)
par(mar=c(3, 4.5, 1, 1))
barplot(t(as.matrix(bp_dat)[plot_mets,]), beside=T, legend.text = colnames(bp_dat), las=1, ylim=c(0,1),
	col=use_col, args.legend=list(x='topright', bty='n'), ylab='Variation Explained',
	names.arg=expression(R[S],R[M],N[AVG],FD,FC), space=c(0.2,0.8))
barplot(t(as.matrix(bp_dat_res)[plot_mets,]), beside=T, add=T, density=25, col='black',
	axes=F, names.arg=expression(R[S],R[M],N[AVG],FD,FC), space=c(0.2,0.8))
dev.off()



# Effects of topographic position and ecoregion
rda_topo = rda(sampXpcoa ~ TopoPos*Ecoregion, data=samples[rownames(comm),]) 
anova(rda_topo, by='term')
vp = varpart(sampXpcoa, ~TopoPos, ~Ecoregion, data=samples[rownames(comm),])

## Local-Regional Variance Decompostion
use_data = samples_pm[,c(locvars, regvars)]

# Remove observations with missing values
missing = rowSums(is.na(use_data))>0
use_data = use_data[!missing,]
use_pcoa = sampXpcoa[rownames(use_data)[rownames(use_data) %in% rownames(sampXpcoa)],]

# Put in same order: omits samples without lichens
use_data = use_data[rownames(use_pcoa),]

# Model ordered data as integers
for(i in 1:ncol(use_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
	use_data[,i] = as.numeric(use_data[,i])
}

vp_locreg = varpart(use_pcoa, as.matrix(use_data[,locvars]), use_data[,regvars])
varpart_CompM = vp_locreg$part$indfract[,'Adj.R.squared']

# Combine and plot
varpart = cbind(varpart_richS, varpart_richM, varpart_totabun, varpart_avgabun, varpart_fd, varpart_CompM)[1:3,]
colnames(varpart) = c('Rich_S','Rich_M','Tot_abun','Avg_abun','FD','Comp_M')

## Manuscript figure
svg('./Figures/variance partition rich abun.svg', height=4, width=4)
barplot(as.matrix(varpart[,c('Rich_S','Rich_M','Avg_abun','FD','Comp_M')]), legend.text = c('Plot/Site','Both','Sample/Tree'), las=1, ylim=c(0,1),
	args.legend=list(x='topright', bty='n'), ylab='Variation Explained',
	names.arg=expression(R[S],R[M],N[AVG],FD,FC))

dev.off()

write.csv(varpart, './Figures/local reg varpart.csv')
varpart = read.csv('./Figures/local reg varpart.csv', row.names=1)

# Write out table of R2
Local = varpart[1,]+varpart[2,]
Regional = varpart[3,]+varpart[2,]
All = colSums(varpart)
write.csv(data.frame(Local, Regional, All), './Figures/local regional R2.csv')


## NOT FINISHED - DO LATER
## Plot locations of samples in trait space from PCOA axes (Manuscript Figure)
use_axes = 1:2

# Calculate percent of variation that each axis explains
ord_sum = eigs/sum(eigs)
pcts = paste(format(ord_sum[use_axes]*100, digits=2), '%', sep='')



# Fit vector of traits
fit_traits = subset(use_traits, Level==1)$TraitName
Tdata = morphos_pm[,fit_traits]
Tdata$Attachment = as.numeric(Tdata$Attachment)
Tdata$Asco = 1*as.logical(Tdata$Asco)
Tdata$Asexual = 1*as.logical(Tdata$Asexual)
ev_trait = envfit(mds, Tdata, choices=use_axes)

# Fit vectors of env variables
Xdata = use_data[,c(locvars,regvars)]
ev_env = envfit(use_pcoa, Xdata, choice=use_axes)

# Plot
par(mar=c(5,5,1,1))
paf(mfrow=c(1,2))

# Samples with env variables
use_pch = c(21,22)
pch_fact = samples[rownames(comm),'Ecoregion']
use_col = c('grey50','white')
col_fact = samples[rownames(comm),'TopoPos']

plot.new()
plot.window(xlim=range(sampXpcoa[,use_axes[1]]), ylim=range(sampXpcoa[,use_axes[2]]))
abline(h=0,v=0, lty=2)
points(sampXpcoa[,use_axes], pch=use_pch[pch_fact], bg=use_col[col_fact])
axis(1)
axis(2, las=1)
mtext(paste('PCoA',use_axes[1],' (',pcts[1], ')', sep=''),1,2)
mtext(paste('PCoA',use_axes[2],' (',pcts[2], ')', sep=''),2,2)
box()

# environment
add_arrow = function(x, col, label){
	arrows(0,0,x[1],x[2], col=col, length=0.1, lwd=2)
	text(x[1],x[2], labels=label, adj=1*(x<0), col=col)
}

# Morphotypes with traits
use_pch = c(21,22)
pch_fact = samples[rownames(comm),'Ecoregion']
use_col = c('grey50','white')
col_fact = samples[rownames(comm),'TopoPos']

plot.new()
plot.window(xlim=range(sampXpcoa[,use_axes[1]]), ylim=range(sampXpcoa[,use_axes[2]]))
abline(h=0,v=0, lty=2)
points(pcoa_vecs[,use_axes], col=2, pch=4)
axis(1)
axis(2, las=1)
mtext(paste('PCoA',use_axes[1],' (',pcts[1], ')', sep=''),1,2)
mtext(paste('PCoA',use_axes[2],' (',pcts[2], ')', sep=''),2,2)
box()

# traits: decided to just show attachment, form and reproduction
for(i in c('Attachment','Asco','Asexual')){
	add_arrow(ev_trait$vectors$arrows[i,], col='black', label=i)
}
for(i in rownames(ev_trait$factors$centroids)[grep('Form',rownames(ev_trait$factors$centroids))]){
	label = sub('Form','', i)
	text(ev_trait$factors$centroids[i,], label=label, col='black')
}

#################################################################
### Variation in single traits
library(lme4)
library(reshape)
library(MuMIn)

## Overall means and distributions across plots

pdf('./Figures/Trait distributions across plots.pdf', height=6, width=7)
for(i in 1:nrow(traitdf)){
	use_trait = traitdf[i,]
	
	par(mar=c(5,5,5,7))
	if(use_trait$type=='numeric'){
		boxplot(lichens[,use_trait$TraitName]~factor(lichens$PlotID, levels=ordered_plots),
			las=2, main=use_trait$TraitName)
	} else {
		tabled = xtabs(AbunCount~factor(PlotID, levels=ordered_plots)+lichens[,use_trait$TraitName], data=lichens)
		tabled = tabled/rowSums(tabled)
		barplot(t(tabled), main=use_trait$TraitName, las=2, ylab='Proportional of Total Abundance',
			legend.text=T, args.legend=list(x='right', inset=-0.25, bty='n'))
	}
}
dev.off()



### Binary Traits
samples_pm$SampID = factor(samples_pm$SampID)

bin_traits = c('Photobiont','Pseudocyphellae','Maculae','Asco','Asexual','Cilia','Rhizines','AscoCover','AsexualForm')
scales = c('Tree','Plot','Site','Ecoregion')
scale_rpt_bin = array(NA, dim=c(length(bin_traits), length(scales), 5), 
	dimnames=list(Trait=bin_traits, Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'))
)

for(i in bin_traits){
	# Remove missing data
	lichdata = subset(lichens_pm, !is.na(lichdata[,i]))

	# Drop levels so that samples without lichens with this trait are not modeled
	lichdata = droplevels(lichdata)
	
	# Calculate binomial response
	use_y = as.matrix(calc_bin_abun(i, lichdata))

	# Define scale variables
	Xdata = samples_pm[rownames(use_y),c('SampID','TreeID','PlotID','PairID','Ecoregion')]
	Xdata = droplevels(Xdata)	

	# Fit models
	mod_tree = glmer(use_y ~ 1 + (1|PlotID) + (1|TreeID) + (1|SampID), data=Xdata, family=binomial)
	mod_plot = glmer(use_y ~ 1 + (1|PlotID) + (1|PairID) + (1|SampID), data=Xdata, family=binomial)
	mod_site = glmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID) + (1|SampID), data=Xdata, family=binomial)
	mod_eco = glmer(use_y ~ 1 + (1|Ecoregion) + (1|SampID), data=Xdata, family=binomial)
	obsID = 'SampID'
	
	# Calculate adjusted repeatability (on link scale)
	scale_rpt_bin[i,'Tree',1:4] = as.numeric(calc_rpt(mod_tree, focalID='TreeID', adjID='PlotID', obsID=obsID))
	scale_rpt_bin[i,'Plot',1:4] = as.numeric(calc_rpt(mod_plot, focalID='PlotID', adjID='PairID', obsID=obsID))
	scale_rpt_bin[i,'Site',1:4] = as.numeric(calc_rpt(mod_site, focalID='PairID', adjID='Ecoregion', obsID=obsID))
	scale_rpt_bin[i,'Ecoregion','R.focal'] = as.numeric(calc_rpt(mod_eco, focalID='Ecoregion', obsID='SampID'))

	scale_rpt_bin[i,,'R.residual'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}
scale_rpt_bin[,'Ecoregion',c('R.focal.adj','R.adj','R.adj.focal')] = scale_rpt_bin[,'Ecoregion','R.focal']

# Save repeatability to data frame
bin_mods_df = cast(melt(scale_rpt_bin), Trait~Model, subset=Statistic=='R.focal.adj')


## OLDER VERIONS OF MODELS
# This function models variation in a binary trait across scales for both pres/abs and abun-weighted probabilities
make_scalemods_bin = function(i, lichdata, sampdata){

	# Remove missing data
	lichdata = subset(lichdata, !is.na(lichdata[,i]))

	# Drop levels so that samples without lichens with this trait are not modeled
	lichdata = droplevels(lichdata)

	# Calculate response variables: # TRUE | # OBS
	y_pres = as.matrix(calc_bin_pres(i, lichdata))
	y_abun = as.matrix(calc_bin_abun(i, lichdata))

	# Create table of predictors
	rownames(sampdata) = sampdata$SampID
	Xdata = sampdata[rownames(y_pres),c('TreeID','PlotID','SiteID','Ecoregion','TopoPos')]
	
	var_scales = sapply(list(y_pres, y_abun), function(y){
		Xdata$y = y
		mod_tree = glmer(y ~ 1 + (1|PlotID/TreeID), data=Xdata, family='binomial')
		mod_plot = glmer(y ~ 1 + (1|SiteID/PlotID), data=Xdata, family='binomial')
		mod_site = glmer(y ~ 1 + (1|SiteID) + Ecoregion, data=Xdata, family='binomial')
		mod_eco = glm(y ~ Ecoregion, data=Xdata, family='binomial')
		
		components = sapply(list(mod_tree, mod_plot, mod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
		components = c(components, attr(r.squaredLR(mod_eco), 'adj.r.squared'))
		names(components) = c('Tree','Plot','Site','Ecoregion')
		components
	})
	colnames(var_scales) = c('Presence','Abundance')

	var_scales
}

## Calculate variance components across scales for all models

bin_traits = c('Photobiont','Pseudocyphellae','Maculae','Asco','Asexual','Cilia','Rhizines','AscoCover','AsexualForm')
bin_mods = array(NA, dim=c(length(bin_traits), 4, 2), dimnames=list(trait=bin_traits, scale=c('Tree','Plot','Site','Region'), response=c('Presence','Abundance')))
for(i in bin_traits){
	bin_mods[i,,] = make_scalemods_bin(i, lichens_pm, samples_pm)
}
# Be careful interpreting models with 0s because maybe not enough observations at certain levels to fit
# Photobiont models don't fit b/c only 19 samples had a cyanolichen: 5 were in Piedmont
# Secondary trait models shouldn't be interpreted at the Tree level because there is not enough replication
#bin_mods_df = cast(melt(bin_mods), response+trait~scale)

pdf('./Figures/variation in traits across scales.pdf', height=5, width=7)
par(mar=c(7.5,4.5,1,3))
barplot(t(bin_mods[,,'Abundance']), beside=T, las=3, ylab=expression(R^2), ylim=c(0,1))
axis(4)
dev.off()

### Examine trait differences across Ecoregions
eco_mods_bin = lapply(bin_traits, function(i){
	# Remove missing data
	lichdata = subset(lichens_pm, !is.na(lichdata[,i]))

	# Drop levels so that samples without lichens with this trait are not modeled
	lichdata = droplevels(lichdata)

	# Calculate response variables: # TRUE | # OBS
	y_abun = as.matrix(calc_bin_abun(i, lichdata))

	# Create table of predictors
	Xdata = samples_pm[rownames(y_abun),]
	Xdata$y = y_abun
	
	# Fit model
	mod_eco = glm(y ~ Ecoregion, data=Xdata, family='binomial')

	mod_eco
})
names(eco_mods_bin) = bin_traits

lapply(eco_mods_bin, summary)

# Test for significant difference between Ecoregions
sapply(eco_mods_bin, function(mod){
	mod_null = update(mod, .~.-Ecoregion, data=mod$model)
	anova(mod, mod_null, test='Chisq')[2,]
})
# Maculae, Asco, Asexual, Cilia, AscoCover, AsexualForm

# Plot differences across ecoregion
for(i in c('Asco','Asexual','AscoCover')){
	lichdata = subset(lichens_pm, !is.na(lichdata[,i]))
	lichdata = droplevels(lichdata)
	y_abun = as.matrix(calc_bin_abun(i, lichdata))
	y = y_abun[,1]/y_abun[,2]
	Ecoregion = samples_pm[rownames(y_abun),'Ecoregion']
	
	pdf(paste('./Figures/sample',i,'by ecoregion.pdf'), height=4, width=4)
	par(mar=c(4,4,1,1))
	boxplot(y~Ecoregion, las=1, ylab='Proportion')
	dev.off()
}


### Categorical Traits
#library(MCMCglmm) - could use for mixed models
#library(nnet)

# We will only analyze probability of crutose and fruticose, other traits probably not meaningful or reliable
sum(is.na(lichens_pm$Form)) # should be 0
form_abun = xtabs(AbunCount ~ SampID + Form, data=lichens_pm)
form_pres = xtabs(~ SampID + Form, data=lichens_pm)

crust_abun = cbind(form_abun[,'crustose'], rowSums(form_abun))
crust_pres = cbind(form_pres[,'crustose'], rowSums(form_pres))
frut_abun = cbind(form_abun[,'fruticose'], rowSums(form_abun))
frut_pres = cbind(form_pres[,'fruticose'], rowSums(form_pres))

# Create table of predictors
Xdata = samples_pm[rownames(form_abun),c('TreeID','PlotID','PairID','Ecoregion','SampID')]

scale_rpt_cat = array(NA, dim=c(ncol(form_abun), length(scales), 5), 
	dimnames=list(Trait=colnames(form_abun), Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'))
)

for(i in colnames(form_abun)){
	# Calculate binomial response
	use_y = cbind(form_abun[,i], rowSums(form_abun))

	# Fit models
	mod_tree = glmer(use_y ~ 1 + (1|PlotID) + (1|TreeID) + (1|SampID), data=Xdata, family=binomial)
	mod_plot = glmer(use_y ~ 1 + (1|PlotID) + (1|PairID) + (1|SampID), data=Xdata, family=binomial)
	mod_site = glmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID) + (1|SampID), data=Xdata, family=binomial)
	mod_eco = glmer(use_y ~ 1 + (1|Ecoregion) + (1|SampID), data=Xdata, family=binomial)
	obsID = 'SampID'
	
	# Calculate adjusted repeatability (on link scale)
	scale_rpt_cat[i,'Tree',1:4] = as.numeric(calc_rpt(mod_tree, focalID='TreeID', adjID='PlotID', obsID=obsID))
	scale_rpt_cat[i,'Plot',1:4] = as.numeric(calc_rpt(mod_plot, focalID='PlotID', adjID='PairID', obsID=obsID))
	scale_rpt_cat[i,'Site',1:4] = as.numeric(calc_rpt(mod_site, focalID='PairID', adjID='Ecoregion', obsID=obsID))
	scale_rpt_cat[i,'Ecoregion','R.focal'] = as.numeric(calc_rpt(mod_eco, focalID='Ecoregion', obsID='SampID'))

	scale_rpt_cat[i,,'R.residual'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}
scale_rpt_cat[,'Ecoregion',c('R.focal.adj','R.adj','R.adj.focal')] = scale_rpt_cat[,'Ecoregion','R.focal']

# Save repeatability to data frame
cat_mods_df = cast(melt(scale_rpt_cat), Trait~Model, subset=Statistic=='R.focal.adj')


# OLDER VERSIONS OF MODELS
form_mods = sapply(list(crust_pres, crust_abun, frut_pres, frut_abun), function(y){
	Xdata$y = y
	mod_tree = glmer(y ~ 1 + (1|PlotID/TreeID), data=Xdata, family='binomial')
	mod_plot = glmer(y ~ 1 + (1|SiteID/PlotID), data=Xdata, family='binomial')
	mod_site = glmer(y ~ 1 + (1|SiteID) + Ecoregion, data=Xdata, family='binomial')
	mod_eco = glm(y ~ Ecoregion, data=Xdata, family='binomial')
		
	components = sapply(list(mod_tree, mod_plot, mod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
	components = c(components, attr(r.squaredLR(mod_eco), 'adj.r.squared'))
	names(components) = c('Tree','Plot','Site','Region')
	components
})
colnames(form_mods) = c('PropCrustose.Presence','PropCrustose.Abundance','PropFruticose.Presence','PropFruticose.Abundance')

form_mods_df = t(form_mods)
levelnames = t(sapply(strsplit(rownames(form_mods_df), '\\.'), function(x) c(trait=x[1], response=x[2])))
for(i in 1:ncol(levelnames)) levelnames[,i] = unlist(levelnames[,i])
form_mods_df = data.frame(levelnames[,2:1],form_mods_df)

# Test for differences between ecoregions
crust_eco = glm(crust_abun~Ecoregion, data=Xdata, family='binomial')
frut_eco = glm(frut_abun~Ecoregion, data=Xdata, family='binomial')
sapply(list(crust_eco, frut_eco), function(mod){
	mod_null = update(mod, .~.-Ecoregion, data=mod$model)
	anova(mod, mod_null, test='Chisq')[2,]
})

Xdata = droplevels(Xdata)

pdf('./Figures/sample prop growth form by ecoregion.pdf', height=4, width=8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
boxplot(I(crust_abun[,1]/crust_abun[,2])~Ecoregion, data=Xdata, las=1, ylab='Proportion Crustose')
boxplot(I(frut_abun[,1]/frut_abun[,2])~Ecoregion, data=Xdata, las=1, ylab='Proportion Fruticose')
dev.off()

### Numeric Traits
## We'll treat attachment as a numeric trait and not bother modeling reproductive effort

num_traits = c( 'Attachment', 'LobeArea', 'LobeDissect', 'AscoAbun', 'AsexualAbun')
summary(lichens_pm[,num_traits])

## Calculate sample mean traits

# Create relative presence/abundance community data matrices
sampXmorph_pres = sampXmorph>0
sampXmorph_abun = sampXmorph

# Calculate sample mean traits based on pres-abs vs. abun
sampXtrait_pres = sapply(num_traits, function(i){
	y = unclass(morphos[,i])
	comm = sampXmorph_pres[,rownames(morphos)]

	# Ignore lichens where the trait is missing
	comm[,is.na(y)] = 0
	y[is.na(y)] = 0

	# Calculate relative abundance
	comm = as.matrix(comm/rowSums(comm)	)
	
	# Calculate community-weighted mean
	comm%*%as.numeric(y)
})
sampXtrait_abun = sapply(num_traits, function(i){
	y = unclass(morphos[,i])
	comm = sampXmorph_abun[,rownames(morphos)]

	# Ignore lichens where the trait is missing
	comm[,is.na(y)] = 0
	y[is.na(y)] = 0

	# Calculate relative abundance
	comm = as.matrix(comm/rowSums(comm)	)
	
	# Calculate community-weighted mean
	comm%*%as.numeric(y)
})

## Examine distributions
layout(matrix(1:10, nrow=2, byrow=F))
for(i in num_traits){
	hist(sampXtrait_pres[,i], main=i)
	hist(sampXtrait_abun[,i], main=i)
}

# Examine distributional assumptions
sum(rownames(sampXtrait_abun)!=rownames(samples))
layout(matrix(1:10, nrow=2, byrow=F))
for(i in num_traits){
	means = tapply(sampXtrait_abun[,i], samples[rownames(sampXmorph),'PlotID'], mean, na.rm=T)
	vars = tapply(sampXtrait_abun[,i], samples[rownames(sampXmorph),'PlotID'], var, na.rm=T)
	plot(vars~means, main=i)
	means = tapply(sampXtrait_pres[,i], samples[rownames(sampXmorph),'PlotID'], mean, na.rm=T)
	vars = tapply(sampXtrait_pres[,i], samples[rownames(sampXmorph),'PlotID'], var, na.rm=T)
	plot(vars~means, main=i)
}

# Variance increases with mean for all but AsexualAbun
# Maybe model log
Ydata = sampXtrait_abun
Ydata[,1:4] = log(Ydata[,1:4])
par(mfrow=c(1,5))
for(i in num_traits){
	means = tapply(Ydata[,i], samples[rownames(sampXmorph),'PlotID'], mean, na.rm=T)
	vars = tapply(Ydata[,i], samples[rownames(sampXmorph),'PlotID'], var, na.rm=T)
	plot(vars~means, main=i)
}

# Create data frame for modeling
Ydata_abun = data.frame(sampXtrait_abun)
rownames(Ydata_abun) = rownames(sampXmorph)
Ydata_abun$AscoAbun = Ydata_abun$AscoAbun + 1 # Because this quantity is 0 whenever there was only one ascoma
Ydata_abun[,1:4] = log(Ydata_abun[,1:4]) # Because quantities are strictly positive and mean increases with variance
Ydata_pres = data.frame(sampXtrait_pres)
rownames(Ydata_pres) = rownames(sampXmorph)
Ydata_pres$AscoAbun = Ydata_pres$AscoAbun + 1 # Because this quantity is 0 whenever there was only one ascoma
Ydata_pres[,1:4] = log(Ydata_pres[,1:4]) # Because quantities are strictly positive and mean increases with variance

# Subset to plots in Piedmont and mountains
keep_samps = rownames(Ydata_abun)[rownames(Ydata_abun) %in% as.character(samples_pm$SampID)]
Ydata_pres = Ydata_pres[keep_samps,]
Ydata_abun = Ydata_abun[keep_samps,]
rownames(samples_pm) = samples_pm$SampID

# Define predictors
Xdata = samples_pm[keep_samps, c('SampID','TreeID','PlotID','PairID','Ecoregion')]

scale_rpt_num = array(NA, dim=c(ncol(Ydata_abun), length(scales), 5), 
	dimnames=list(Trait=colnames(Ydata_abun), Model=scales, Statistic = c('R.focal','R.focal.adj','R.adj','R.adj.focal','R.residual'))
)

for(i in colnames(Ydata_abun)){
	# Calculate binomial response
	use_y = Ydata_abun[,i]

	# Fit models
	mod_tree = lmer(use_y ~ 1 + (1|PlotID) + (1|TreeID), data=Xdata)
	mod_plot = lmer(use_y ~ 1 + (1|PlotID) + (1|PairID), data=Xdata)
	mod_site = lmer(use_y ~ 1 + (1|Ecoregion) + (1|PairID), data=Xdata)
	mod_eco = lmer(use_y ~ 1 + (1|Ecoregion), data=Xdata)
	obsID = 'Residual'
	
	# Calculate adjusted repeatability (on link scale)
	scale_rpt_num[i,'Tree',1:4] = as.numeric(calc_rpt(mod_tree, focalID='TreeID', adjID='PlotID', obsID=obsID))
	scale_rpt_num[i,'Plot',1:4] = as.numeric(calc_rpt(mod_plot, focalID='PlotID', adjID='PairID', obsID=obsID))
	scale_rpt_num[i,'Site',1:4] = as.numeric(calc_rpt(mod_site, focalID='PairID', adjID='Ecoregion', obsID=obsID))
	scale_rpt_num[i,'Ecoregion','R.focal'] = as.numeric(calc_rpt(mod_eco, focalID='Ecoregion', obsID='SampID'))

	scale_rpt_num[i,,'R.residual'] = as.numeric(sapply(c(mod_tree, mod_plot, mod_site, mod_eco), function(x) get_allvar(x)[obsID]/sum(get_allvar(x))))	
}
scale_rpt_num[,'Ecoregion',c('R.focal.adj','R.adj','R.adj.focal')] = scale_rpt_num[,'Ecoregion','R.focal']

# Save repeatability to data frame
num_mods_df = cast(melt(scale_rpt_num), Trait~Model, subset=Statistic=='R.focal.adj')

# OLDER MODELS
num_mods = array(NA, dim=c(length(num_traits), 4, 2), dimnames=list(trait=num_traits, scale=c('Tree','Plot','Site','Region'), response=c('Presence','Abundance')))

make_scalemods_num = function(y_pres, y_abun, sampdata){
		
	# Create table of predictors
	Xdata = sampdata[,c('TreeID','PlotID','SiteID','Ecoregion','TopoPos')]
	
	# Remove missing data
	Xdata = subset(Xdata, !is.na(y_pres))
	y_abun = y_abun[!is.na(y_pres)]
	y_pres = y_pres[!is.na(y_pres)]
	
	# Drop levels so that samples without lichens with this trait are not modeled
	Xdata = droplevels(Xdata)

	var_scales = sapply(list(y_pres, y_abun), function(y){
		Xdata$y = y
		mod_tree = lmer(y ~ 1 + (1|PlotID/TreeID), data=Xdata)
		mod_plot = lmer(y ~ 1 + (1|SiteID/PlotID), data=Xdata)
		mod_site = lmer(y ~ 1 + (1|SiteID) + Ecoregion, data=Xdata)
		mod_eco = glm(y ~ Ecoregion, data=Xdata, family='gaussian')
		components = sapply(list(mod_tree, mod_plot, mod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
		components = c(components, r.squaredLR(mod_eco))
		names(components) = c('Tree','Plot','Site','Ecoregion')
		components
	})
	colnames(var_scales) = c('Presence','Abundance')

	var_scales
}

for(i in num_traits){
	# Calculate response variables
	y_pres = Ydata_pres[,i]
	y_abun = Ydata_abun[,i]

	# Create models
	num_mods[i,,] = make_scalemods_num(y_pres, y_abun, samples_pm)
}

num_mods_df = cast(melt(num_mods), response+trait~scale)


## Test effect of Ecoregion
eco_mods_num = lapply(num_traits, function(i){
	y_abun = Ydata_abun[,i]

	# Remove missing data
	Xdata = subset(samples_pm, !is.na(y_abun))
	y_abun = y_abun[!is.na(y_abun)]
	
	# Drop levels so that samples without lichens with this trait are not modeled
	Xdata = droplevels(Xdata)

	# Fit model
	Xdata$y = y_abun
	mod_eco = glm(y ~ Ecoregion, data=Xdata, family='gaussian')

	mod_eco
})
names(eco_mods_num) = num_traits

# Test for significant difference between Ecoregions
sapply(eco_mods_num, function(mod){
	mod_null = update(mod, .~.-Ecoregion, data=mod$model)
	anova(mod, mod_null, test='Chisq')[2,]
})
# Attachment, LobeArea, LobeDissect, AscoAbun, AsexualAbun


# Plot differences across ecoregion
for(i in c('Attachment')){
	y = Ydata_abun[,i]
	Ecoregion = samples_pm[rownames(Ydata_abun),'Ecoregion']
	
	pdf(paste('./Figures/sample',i,'by ecoregion.pdf'), height=4, width=4)
	par(mar=c(4,4,1,1))
	boxplot(y~Ecoregion, las=1, ylab=paste('Mean',i))
	dev.off()
}

# Combine all trait models
trait_mods_df = rbind(bin_mods_df, cat_mods_df, num_mods_df)
trait_mods_df = trait_mods_df[,c('Trait','Tree','Plot','Site','Ecoregion')]

# Add a column for the analysis method
trait_mods_df$Type = use_traits[as.character(trait_mods_df$Trait),'type']
trait_mods_df$Type = ifelse(trait_mods_df$Type %in% c('ordered','numeric'), 'Gaussian', 'Binomial') 

# Add a column for the trait category
trait_mods_df$Category = use_traits[as.character(trait_mods_df$Trait), 'category']
trait_mods_df[grep('ose', as.character(trait_mods_df$Trait)),'Category'] = 'form'

# Reorder rows and save
trait_mods_df = trait_mods_df[with(trait_mods_df, order(Category, Type, Trait)),]
write.csv(trait_mods_df, './Figures/variation in traits across scales.csv', row.names=F)


## Save all repeatability arrays
save(scale_rpt, scale_rpt_res, scale_rpt_bin, scale_rpt_cat, scale_rpt_num, file='repeatability_statistic_arrays.RData')

# OLD MODELS
#colorder = c('response','trait','Tree','Plot','Site','Region')
#trait_mods_df = rbind(bin_mods_df[,colorder], form_mods_df[,colorder], num_mods_df[, colorder])

# Add a column for the analysis method
#trait_mods_df$type = use_traits[as.character(trait_mods_df$trait),'type']
#trait_mods_df$type = ifelse(trait_mods_df$type %in% c('ordered','numeric'), 'Gaussian', 'Binomial') 

# Add a column for the trait category
#trait_mods_df$category = use_traits[as.character(trait_mods_df$trait), 'category']
#trait_mods_df[grep('Crustose|Fruticose', as.character(trait_mods_df$trait)),'category'] = 'form'

# Reorder rows and save
#trait_mods_df = trait_mods_df[with(trait_mods_df, order(response, category, type, trait)),]
#write.csv(trait_mods_df, './Figures/variation in traits across scales.csv', row.names=F)

#################################################################
### Variation in multi-trait FD

## FD data already loaded by load_data.R script

## Variation in FD across sites and plots
#FD = FD_rep

factorby = factor(samples[FD$SampID,'PlotID'], levels=ordered_plots)

# z-score FD across plots
par(mfrow=c(1,2))
bp_data = FD_scaled[,'z','withinplot','0.5','swsh_samp_c']
boxplot(bp_data~factorby, las=2)
bp_data = FD_scaled[,'z','withinsite','0.5','swsh_samp_c']
boxplot(bp_data~factorby, las=2)

# Unscaled FD across plots
pdf('./Figures/FD rep across plots s0.5.pdf', height=4, width=5)
par(mar=c(6,4,1,1))
boxplot(FD$rao_L2_s.5~factorby, las=2)
dev.off()

# Hierarchical model of unscaled FD across region, site, plot, tree
# omit Bladen1
library(lme4)

use_data = merge(FD, samples)
use_data = merge(use_data, plot_data)
use_data = subset(use_data, PlotID!='Bladen1')
use_data_no0 = subset(use_data, rao_L2_s.5!=0)

# Drop un-used levels
use_data = droplevels(use_data)

# Examine distributional assumptions: Gaussian errors ok
means = with(use_data, tapply(rao_L2_s.5, PlotID, mean, na.rm=T))
vars = with(use_data, tapply(rao_L2_s.5, PlotID, var, na.rm=T))
plot(means, vars)
hist(use_data$rao_L2_s.5)

# Models including FD=0
FDmod_tree = lmer(rao_L2_s.5 ~ 1 + (1|PlotID/TreeID), data=use_data)
FDmod_plot = lmer(rao_L2_s.5 ~ 1 + (1|SiteID/PlotID), data=use_data)
FDmod_site = lmer(rao_L2_s.5 ~ 1 + (1|SiteID) + Ecoregion, data=use_data)
FDmod_eco = glm(rao_L2_s.5 ~ Ecoregion, family='gaussian', data=use_data)

components = sapply(list(FDmod_tree, FDmod_plot, FDmod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
components = c(components, r.squaredLR(FDmod_eco))
names(components) = c('Tree','Plot','Site','Ecoregion')

# Models with ecoregion and topographic position
mod_full = lmer(rao_L2_s.5 ~ 1 + Ecoregion + TopoPos + (1|SiteID/TreeID), data=use_data)
mod_inter = lmer(rao_L2_s.5 ~ 1 + Ecoregion*TopoPos + (1|SiteID/TreeID), data=use_data)
anova(mod_inter, mod_full) # significant negative interaction- no effect of TopoPos in Piedmont

# Calculate variance of fixed effects
sigma.eco = var(as.numeric(as.vector(fixef(mod_full)[1:2]) %*% t(model.matrix(mod_full)[,1:2])))
sigma.plot = var(as.numeric(as.vector(fixef(mod_full)[c(1,3)]) %*% t(model.matrix(mod_full)[,c(1,3)])))
sigma.fixed = calc_fixedvar(mod_full)
	
# Extract variance of random effects
sigma.tree = data.frame(VarCorr(mod_full))$vcov[1]
sigma.site = data.frame(VarCorr(mod_full))$vcov[2]

# Calculate total variance and R2 (based on Nakagawa and Schielzeth 2013)
sigma.tot = calc_totvar(mod_full)

# Conditional R2 (variance explained by all factors)
(sigma.fixed + sigma.tree + sigma.site) / sigma.tot

# Plot variance components
components = c(sigma.tree, sigma.plot, sigma.site, sigma.eco)/sigma.tot
names(components) = c('Tree','Plot','Site','Ecoregion')

pdf('./Figures/variance across scales topo-eco sample FD.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.7), ylab=expression(paste('Marginal ', R^2)), main='Sample Total Abundance')
dev.off()

# Plot effects of ecoregion and plot in the interaction model
counts = with(use_data, aggregate(rao_L2_s.5, list(TopoPos, Ecoregion), FUN=function(x) table(cut(x, breaks=seq(0,.25,0.01), include.lowest=T))))
bars = counts[,-(1:2)] / rowSums(counts[,-(1:2)])

use_col = c('grey50', 'white')

pdf('./Figures/sample FD across ecoregions and topo.pdf', height=4, width=4)
par(mar=c(4,4,1,1))

# Set up plot
plot(c(-1, 1), c(0,.25), type='n', las=1, axes=F, xlim=c(-2,2), xlab='', ylab='Sample Functional Diversity')
axis(1, at=c(-1,1), labels=c('Mountains','Piedmont'))
axis(2, las=1)
box()

# Make rectangles
rect(-1, (0:24)/100, -1-bars[1,], (1:25)/100, col=use_col[1])
rect(-1, (0:24)/100, -1+bars[2,], (1:25)/100, col=use_col[2])
rect(1, (0:24)/100, 1-bars[3,], (1:25)/100, col=use_col[1])
rect(1, (0:24)/100, 1+bars[4,], (1:25)/100, col=use_col[2])

# Add estimated effects from model with interaction
ests =  exp(fixef(mod_inter))

modmat = expand.grid(TopoPos=c('sheltered','exposed'),Ecoregion=c('Mountains','Piedmont'))
means = predict(mod_inter, modmat, re.form=~0, type='response')

points(rep(c(-1,1), each=2) + rep(c(-.3, .3), 2), means, pch=3, bg=use_col, cex=2, lwd=2, lend=1)

dev.off()




## CODE BELOW NOT DONE YET

FDmod_no0 = lmer(rao_L2_s.5~ (1|Ecoregion/SiteID/PlotID/TreeID), data=use_data_no0)
summary(FDmod_no0)

# Plot variance components
FDmod_vars = c(as.numeric(VarCorr(FDmod)), attr(VarCorr(FDmod),'sc')^2)
FDmod_vars = FDmod_vars/sum(FDmod_vars)
FDmod_no0_vars = c(as.numeric(VarCorr(FDmod_no0)), attr(VarCorr(FDmod_no0),'sc')^2)
FDmod_no0_vars = FDmod_no0_vars/sum(FDmod_no0_vars)

barmids = matrix(FDmod_no0_vars, 5, 5)
barmids[lower.tri(barmids)] = 0
barmids = colSums(barmids)
barmids = c(0,barmids[1:4])+FDmod_no0_vars/2

pdf('./Figures/FD variance components.pdf', height=4, width=5)
par(mar=c(2,4,1,1))
barplot(cbind(FDmod_vars, FDmod_no0_vars), names.arg=c('With FD=0', 'Without FD=0'),
	las=1, ylab='Proportion variance explained', xlim=c(0,3))->bp
text(2.5, barmids, c('Regions','Sites','Plots','Trees','Residual'), adj=0)
dev.off()
