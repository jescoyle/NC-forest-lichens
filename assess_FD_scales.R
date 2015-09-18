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

## Define variables used for variation partitioning in each section
# Choice of regional variables based on correlations between plot level climate variables:
# Strong correlation between Elevation, humidity and temperature, elevation least correlated with AP and CloudFreq_sd
sampvars = c('Angle','Bryophytes','Shedding','FurrowDepth','pH','Density','WaterCapacity')
treevars = c('DBH','Trans_tot_cor')
regvars = c('Elevation','CloudFreq_sd','AP','OpenPos','Soil_pH')



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

use_y = samples_pm$Avg_abun + 1 # Rich_M, Rich_S, Tot_abun_cor, 

hist(use_y)
hist(log(use_y))

# Check whether Poisson distribution is appropriate: 
# looks like it is for richness and avg abun, tot abundance is a bit over dispersed (Tot_abun_cor = 3.8)
vars = with(samples_pm, tapply(use_y, PlotID, var))
means = with(samples_pm, tapply(use_y, PlotID, mean))
plot(means, vars); abline(0,1)
abline(lm(vars~means-1), col=2)
lm(vars~means-1)

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
sigma.tot = calc_totvar(mod_full)

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


## Effects of scale itself
mod_tree = glmer(use_y ~ 1 + (1|PlotID/TreeID), data=samples_pm, family='poisson')
mod_plot = glmer(use_y ~ 1 + (1|SiteID/PlotID), data=samples_pm, family='poisson')
mod_site = glmer(use_y ~ 1 + (1|SiteID) + Ecoregion, data=samples_pm, family='poisson')
mod_eco = glm(use_y ~ Ecoregion, data=samples_pm, family='poisson')
mod_full = glmer(use_y ~ 1 + (1|SiteID/PlotID/TreeID) + Ecoregion, data=samples_pm, family='poisson')

# Avg abun models
mod_tree = lmer(use_y ~ 1 + (1|PlotID/TreeID), data=samples_pm)
mod_plot = lmer(use_y ~ 1 + (1|SiteID/PlotID), data=samples_pm)
mod_site = lmer(use_y ~ 1 + (1|SiteID) + Ecoregion, data=samples_pm)
mod_eco = glm(use_y ~ Ecoregion, data=samples_pm, family='gaussian')
mod_full = lmer(use_y ~ 1 + (1|SiteID/PlotID/TreeID) + Ecoregion, data=samples_pm)


# Plot R2 components
components = sapply(list(mod_tree, mod_plot, mod_site), function(x) data.frame(VarCorr(x))$vcov[1]/calc_totvar(x))
components = c(components, attr(r.squaredLR(mod_eco), 'adj.r.squared'))
names(components) = c('Tree','Plot','Site','Ecoregion')

pdf('./Figures/variance across scales nested sample morphospecies richness.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.3), ylab=expression(R^2), main='Sample Morphospecies Richness')
dev.off()


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

## Average abundace models
## Can't use Poisson.



## GET CI FOR VARIANCE COMPONENTS
## NOT IMPLEMENTED YET
	# Bootstrap fixed effects variance estimate only
	varboot = bootMer(mod, calc_fixedvar, nsim=1000, use.u=T, type='parametric')
	ints_boot = boot.ci(varboot, index=1, type='perc')$percent[4:5]

	# Profile likelihood confidence intervals for random and fixed effects
	pf = profile(mod_full, 1:5)

	vpf = varianceProf(pf) # converts to variance scale (by squaring)
	ints = confint(vpf, parm=1:5, method='profile')

#############################
### Variance Partitioning ###

# Sample-scale response
use_y = samples_pm$Avg_abun #+ 1# Rich_M, Rich_S, Tot_abun_cor, Avg_abun
#use_y = FD[rownames(samples_pm), 'rao_L2_s.5'] # Read in below in FD section
locvars = c(sampvars, treevars)
use_data = samples_pm[,c(locvars, regvars)]

# Tree-scale response
use_y = trees_pm$Avg_abun_tree #+ 1 # Rich_S_tree, Rich_M_tree, Tot_abun_tree_cor, Avg_abun_tree
locvars = treevars
use_data = trees_pm[,c(locvars, regvars)]

# Remove observations with missing values
missing = rowSums(is.na(use_data)|is.na(use_y))>0
use_y = use_y[!missing]
use_data = use_data[!missing,]

# Model ordered data as integers
for(i in 1:ncol(use_data)){
	if(is.ordered(use_data[,i])) use_data[,i] = unclass(use_data[,i])
}

locmod = glm(use_y~., data=use_data[,locvars], family='poisson')
regmod = glm(use_y~., data=use_data[,regvars], family='poisson')
fullmod = glm(use_y~., data=use_data, family='poisson')

# For average abundance and FD
locmod = glm(use_y~., data=use_data[,locvars], family='gaussian')
regmod = glm(use_y~., data=use_data[,regvars], family='gaussian')
fullmod = glm(use_y~., data=use_data, family='gaussian')

R2s = sapply(list(local = locmod, regional = regmod, both = fullmod), function(x) r.squaredLR(x))

# Save results
varpart_richS = partvar2(R2s)
varpart_richM = partvar2(R2s)
varpart_totabun = partvar2(R2s)
varpart_avgabun = partvar2(R2s)
varpart_fd = partvar2(R2s)

varpart_richS_tree = partvar2(R2s)
varpart_richM_tree = partvar2(R2s)
varpart_totabun_tree = partvar2(R2s)
varpart_avgabun_tree = partvar2(R2s)



#################################################################
### Variation in morphotypes (equivalent to variation in species composition)
library(ape)
library(vegan)

dist_L2mat = read.csv(paste(derive_dir, 'morpho_distance_matrix_sfact0.5.csv'), row.names=1)

# Subset morphotypes and community matrix to just morphotypes in Piedmont and Mountains
pm_morphs = unique(subset(lichens, SiteID!='Bladen')$MorphID)
morph_dmat = dist_L2mat[pm_morphs, pm_morphs]
morphos_pm = morphos[pm_morphs,]
use_samps = samples_pm$SampID
comm = sampXmorph[use_samps[use_samps %in% rownames(sampXmorph)], pm_morphs]

# PCOA of trait distance matrix
mds = pcoa(morph_dmat, correction='lingoes')

# Define morphotype vectors for biplots
pcoa_vecs = mds$vectors

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
sampXpcoa = as.matrix(comm_hel) %*% mds$vectors 

# Compute RDA with scales as predictors
rda_null = rda(sampXpcoa)

Xdata = samples[rownames(comm),c('TreeID','TopoPos','SiteID','Ecoregion')]
Xdata = apply(Xdata, 2, function(x) as.numeric(x))

rda_full = rda(sampXpcoa ~ TreeID + TopoPos + SiteID + Ecoregion, data=samples[rownames(comm),])

rda_tree = rda(sampXpcoa ~ TreeID + Condition(PlotID), data=samples[rownames(comm),])
rda_plot = rda(sampXpcoa ~ PlotID + Condition(SiteID), data=samples[rownames(comm),])
rda_site = rda(sampXpcoa ~ SiteID + Condition(Ecoregion), data=samples[rownames(comm),])
rda_eco = rda(sampXpcoa ~ Ecoregion, data=samples[rownames(comm),])


# Plot variance components
components = sapply(list(rda_tree, rda_plot, rda_site, rda_eco), function(x) RsquareAdj(x)$adj.r.squared) # Variance explained is adjusted R2 from constrained component
names(components) = c('Tree','Plot','Site','Ecoregion')

pdf('./Figures/variance across scales sample morphotype composition.pdf', height=4, width=4)
par(mar=c(3,5,2,1))
barplot(components, las=1, ylim=c(0,.4), ylab=expression(paste('Adjusted Semi-Partial ', R^2)), main='Sample Functional Composition')
dev.off()


# Effects of topographic position and ecoregion
rda_topo = rda(sampXpcoa ~ TopoPos*Ecoregion, data=samples[rownames(comm),]) 
anova(rda_topo, by='term')
vp = varpart(sampXpcoa, ~TopoPos, ~Ecoregion, data=samples[rownames(comm),])

## Local-Regional Variance Decompostion
use_data = samples_pm[,c(locvars, regvars)]

# Put in same order: omits samples without lichens
use_data = use_data[rownames(use_pcoa),]

# Remove observations with missing values
missing = rowSums(is.na(use_data))>0
use_data = use_data[!missing,]
use_pcoa = sampXpcoa[!missing,]

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

png('./Figures/variance partition rich abun.png', height=400, width=500)
barplot(varpart, legend.text = c('Plot/Site','Both','Sample/Tree'), las=1, ylim=c(0,1),
	args.legend=list(x='topleft', bty='n'), ylab='Variation Explained')
dev.off()

write.csv(varpart, './Figures/local reg varpart.csv')

# Write out table of R2
Local = varpart[1,]+varpart[2,]
Regional = varpart[3,]+varpart[2,]
All = colSums(varpart)
write.csv(data.frame(Local, Regional, All), './Figures/local regional R2.csv')

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
bin_mods_df = cast(melt(bin_mods), response+trait~scale)

pdf('./Figures/variation in traits across scales.pdf', height=5, width=7)
par(mar=c(7.5,4.5,1,3))
barplot(t(bin_mods[,,'Abundance']), beside=T, las=3, ylab=expression(R^2), ylim=c(0,1))
axis(4)
dev.off()


# Examine trait differences across Ecoregions
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
Xdata = samples[rownames(form_pres),c('TreeID','PlotID','SiteID','Ecoregion','TopoPos')]
	
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
Ydata_pres = Ydata_pres[samples_pm$SampID,]
Ydata_abun = Ydata_abun[samples_pm$SampID,]
rownames(samples_pm) = samples_pm$SampID

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
colorder = c('response','trait','Tree','Plot','Site','Region')
trait_mods_df = rbind(bin_mods_df[,colorder], form_mods_df[,colorder], num_mods_df[, colorder])

# Add a column for the analysis method
trait_mods_df$type = use_traits[as.character(trait_mods_df$trait),'type']
trait_mods_df$type = ifelse(trait_mods_df$type %in% c('ordered','numeric'), 'Gaussian', 'Binomial') 

# Add a column for the trait category
trait_mods_df$category = use_traits[as.character(trait_mods_df$trait), 'category']
trait_mods_df[grep('Crustose|Fruticose', as.character(trait_mods_df$trait)),'category'] = 'form'

# Reorder rows and save
trait_mods_df = trait_mods_df[with(trait_mods_df, order(response, category, type, trait)),]
write.csv(trait_mods_df, './Figures/variation in traits across scales.csv', row.names=F)

#################################################################
### Variation in multi-trait FD

## Load FD data caclulated in script: multi_trait_models.R
null_fd_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/FD NUll Sims/'

FD = read.csv('./Data/raos_q_samples.csv')
rownames(FD) = FD$SampID
load(paste(null_fd_dir,'FD_scaled.RData', sep=''))

# FD calculated on reproductve traits
FD_rep = read.csv('./Data/raos_q_reproductive_traits_samples.csv')
rownames(FD_rep) = FD_rep$SampID
load(paste(null_fd_dir,'FD_rep_scaled.RData', sep=''))

# Subset to top 2 samples on each tree
FD = FD[use_top2,]
FD_rep = FD_rep[use_top2,]


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
