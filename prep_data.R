### This script prepares Lichen FD project data tables exported from the SQLite database for analysis.

options(stringsAsFactors=F)

# Load functions
source('C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/lichen_FD_functions.R')


# Define data locations
mydir = 'C://Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/'
sql_dir = './Data/SQLite Tables/'
derive_dir = './Data/Derived Tables/'

setwd(mydir)

# Read in basic data tables from SQLite database
plots = read.csv(paste(sql_dir,'plots.csv', sep=''))
canopy = read.csv(paste(sql_dir,'canopy.csv', sep=''))
trees = read.csv(paste(sql_dir,'trees.csv', sep=''))
bark_traits = read.csv(paste(sql_dir,'bark_traits.csv', sep=''))
samples = read.csv(paste(sql_dir,'samples.csv', sep=''))
lichens = read.csv(paste(sql_dir,'lichens.csv', sep=''))

#######################################################
### Tree-level data

# Assess light bias associated with increasing light levels at dawn
treedata = merge(trees, canopy)

# Subset to plots for analysis
use_treedata = subset(treedata, PlotID!='Bladen1')
use_plots = subset(plots, PlotID!='Bladen1')

# Make date into date class
use_treedata$Date_time = as.POSIXlt(use_treedata$Date_time)

yvars = c('Pct_open','LAI_75deg','Trans_tot')
syms = expand.grid(pch=c(1,16),col=rainbow(9))
rownames(syms) = use_plots$PlotID[order(use_plots$PairID)]

par(mfrow=c(1,3))
for(y in yvars){
	plot(use_treedata[,y]~EVdif, data=use_treedata, ylab=y, 
		pch=syms[use_treedata$PlotID,'pch'], col=syms[use_treedata$PlotID, 'col'])
}

plot(Pct_open~Trans_tot, data=use_treedata, pch=syms[use_treedata$PlotID,'pch'], col=syms[use_treedata$PlotID, 'col'])
plot(LAI_75deg~LAI_60deg, data=use_treedata)

pdf('./Figures/Dependence of light on camera EVdif.pdf', height=4, width=10.5)
for(p in rownames(syms)){
	use_trees = subset(use_treedata, PlotID==p)
	par(mfrow=c(1,3))
	for(y in yvars){
		plot(use_trees[,y]~EVdif, data=use_trees, ylab=y, main=p)
	}
}
dev.off()
plot(EVdif~I(Date_time$hour+Date_time$min/60), data=use_treedata)

# Detrend Trans_tot based on EVdif assumming same relationship between EVdif and Trans_tot across sites
trans_mod = lm(Trans_tot~EVdif, data=treedata)

treedata$Trans_tot_cor = NA
treedata[names(resid(trans_mod)),'Trans_tot_cor'] = resid(trans_mod)

plot(Trans_tot_cor~EVdif, data=treedata)
plot(Trans_tot_cor~ISO, data=treedata)

plot(Trans_tot_cor~factor(TreeTaxonID), data=treedata)
plot(Trans_tot_cor~factor(PlotID), data=treedata, las=3)

keep_cols = c('PlotID','TreeTaxonID','Module','DBH','UncertainID', 'Trans_tot_cor')
rownames(treedata) = treedata$TreeID

write.csv(treedata[,keep_cols], paste(derive_dir,'treedata.csv', sep=''), row.names=T)

treedata = read.csv(paste(derive_dir,'treedata.csv', sep=''), row.names=1)
###############################################################
### Sample level data
library(stringr)

# NOTE: lines not needed to generate samp_data table have been commented out

## Dependence of pH on pH meter
#boxplot(pH~pH_meter, data=bark_traits)

## Trait differences among re-measured samples
#dups = table(bark_traits$SampID)
#dups_many = names(dups[dups>2])

# Remove multiple measurements that should not be analyzed
#subset(bark_traits, SampID==dups_many[9])
remove_samps = c(724,725,726,793,1318,1319,932,895,896,1815,1818,2002)
barkdata = subset(bark_traits, !(BarkMeasID %in% remove_samps))

dups = table(barkdata$SampID)
dups = names(dups[dups>1])

# Calculate bark trait mean and range across bark measurements from the same sample
btraits = c('pH','WaterCapacity','Density','MassDry')

traitrange = aggregate(barkdata[,btraits], by=list(SampID=barkdata$SampID), FUN=function(x) abs(diff(range(x))))
names(traitrange)[names(traitrange)%in%btraits] = paste(btraits,'rng',sep='_')

traitmean = aggregate(barkdata[,btraits], by=list(SampID=barkdata$SampID), FUN=mean)
names(traitmean)[names(traitmean)%in%btraits] = paste(btraits,'mn',sep='_')

meter = aggregate(barkdata$pH_meter, by=list(SampID=barkdata$SampID), FUN=function(x){
	if(length(unique(x))==1){ unique(x) } else {'Both'}
})
names(meter)[2] = 'pH_meter'

samptraits = merge(traitmean, traitrange)
samptraits = merge(samptraits, meter)

#use_data = samptraits[rowSums(samptraits[,paste(btraits, 'rng', sep='_')])>0,]

#plot(pH_rng~pH_mn, data=use_data)
#plot(WaterCapacity_rng~WaterCapacity_mn, data=use_data)
#plot(Density_rng~Density_mn, data=use_data)

#boxplot(pH_rng~pH_meter, data=use_data)
#tapply(samptraits$pH_rng, samptraits$pH_meter, mean)

#pHmod = lm(pH~factor(SampID)+factor(pH_meter), data=subset(barkdata, SampID %in% dups))
#anova(pHmod)

# Check which measurements are bad
use_data = merge(barkdata, samples[,c('SampID','Date')])
use_data$Date = as.Date(use_data$Date, format='%m/%d/%Y')
use_data$year = format(use_data$Date, '%Y')
#aggregate(use_data$pH, by=list(use_data$year, use_data$pH_meter), mean)
#boxplot(pH~year+pH_meter, data=use_data) # only samples from 2013 may be biased
#King2013 = subset(use_data, year=='2013'&pH_meter=='Kingsolver')
#plot(pH~Date, data=King2013)
#nrow(subset(King2013, !(SampID %in%dups))) # Would need to remeasure 773 + 45 samples

# Calculate rank order of samples measured under Kier vs Kingsolver pH meter
#use_samps = subset(samptraits, pH_meter=='Both')$SampID
#kings = subset(barkdata, SampID %in% use_samps & pH_meter=='Kingsolver')
#kiers = subset(barkdata, SampID %in% use_samps & pH_meter=='Kier')
#kings$king_rank = rank(kings$pH)
#kiers$kier_rank = rank(kiers$pH)
#kks = merge(kings[,c('SampID','king_rank','pH')],kiers[,c('SampID','kier_rank','pH')], by='SampID')
#plot(king_rank~kier_rank, data=kks[order(kks$king_rank),])

# Calculate mean difference among samples measured on both pH meters
use_samps = subset(samptraits, pH_meter=='Both')$SampID
#boxplot(pH~pH_meter, data=subset(barkdata, SampID %in% use_samps))
difs = sapply(use_samps, function(x){
	use_data = subset(barkdata, SampID == x)
	a = mean(use_data[use_data$pH_meter=='Kingsolver','pH'])
	b = mean(use_data[use_data$pH_meter=='Kier','pH'])
	diff(c(a,b))
})
mean_pH_diff = round(mean(difs),3)

# Remove duplicate samples from Kingsolver pH meter when measurement from Kier meter is available
use_samps = subset(samptraits, pH_meter=='Both')$SampID
remove_samps = subset(barkdata, (SampID %in% use_samps)&(pH_meter=='Kingsolver'))$BarkMeasID 
barkdata = subset(barkdata, !(BarkMeasID %in% remove_samps))

# Correct for difference between pH meters in 2013 by adding the mean difference from repeated measures
use_samps = subset(use_data, year=='2013'&pH_meter=='Kingsolver')$BarkMeasID
barkdata[barkdata$BarkMeasID %in% use_samps,'pH'] = barkdata[barkdata$BarkMeasID %in% use_samps,'pH'] + mean_pH_diff
#boxplot(pH~pH_meter, data=barkdata)

## Estimate bark trait differences among observers
#boxplot(pH~Observer, data=barkdata)
#boxplot(WaterCapacity~Observer, data=barkdata)
#boxplot(Density~Observer, data=barkdata)
# No worrisome differences

## Filter for unlikely trait values

# Negative water holding capacity
#subset(barkdata, WaterCapacity < 0 )
#subset(barkdata, WaterCapacity < 0 & Notes!='forgot to measure wet mass' )
# modified data where it looked like WetMass and Drymass values had been swapped with neighboring samples or mis-recorded

# High Water Capacity
#subset(barkdata, WaterCapacity >5)
# modified 2 records with MassWet off by 1 decimal place

# Assign remaining negative WaterCapacity values to NA
barkdata[which(barkdata$WaterCapacity<0),'WaterCapacity'] = NA

# Density high
#plot(Density~MassDry, data=barkdata)
#subset(barkdata, Density > 1)
#subset(barkdata, Density > 3)
# corrected a couple obvious errors, but most of these may are too difficult to figure out the possible error
#subset(barkdata, Density < 0.05)

### Generate sample-level environmental data table

sampdata = samples[,c('SampID','TreeID','SampNum','Height','Angle','Bryophytes','Shedding','FurrowDepth','Aspect','NoLichen','Notes')]

# Convert Bryophytes 'several' and 'some' code to same level
sampdata[sampdata$Bryophytes=='some','Bryophytes'] = 'several'

# Convert FurrowDepth to number
sampdata[sampdata$FurrowDepth=='<0.5','FurrowDepth'] = 0.25

# Average bark traits from samples measured more than once
dups = names(which(table(barkdata$SampID)>1))
singles = names(which(table(barkdata$SampID)==1))

barktraitdata = subset(barkdata, SampID %in% singles)[,c('SampID','pH','Density','WaterCapacity','Notes')]
for(this_samp in dups){
	these_bark = subset(barkdata, SampID==this_samp)
	
	masses = these_bark$MassDry

	# If one sample is <10mg, then exclude this data
	if(sum(masses>0.01)>=1){
		these_bark = subset(these_bark, MassDry > 0.01)
	}

	# Calculate average trait
	avgs = apply(these_bark[,c('pH','Density','WaterCapacity')], 2, mean)
	notes = str_trim(paste(these_bark$Notes, collapse=' '))
	
	barktraitdata = rbind(barktraitdata, data.frame(SampID=this_samp, t(avgs), Notes=notes))
}

# Merge sample level data
sampdata[which(!(sampdata$SampID %in% barktraitdata$SampID)),] # Several samples missing bark data
barktraitdata[which(!(barktraitdata$SampID %in% sampdata$SampID)),] # Corresponds to a sample not used in analysis

colnames(sampdata)[colnames(sampdata)=='Notes'] = 'Notes_samp'
colnames(barktraitdata)[colnames(barktraitdata)=='Notes'] = 'Notes_bark'
sampdata = merge(sampdata, barktraitdata, all.x=T, all.y=F)


rownames(sampdata) = sampdata$SampID
keep_cols = c('TreeID','SampNum','Height','Angle','Bryophytes','Shedding','FurrowDepth','Aspect','pH','Density','WaterCapacity','NoLichen','Notes_samp','Notes_bark')
write.csv(sampdata[,keep_cols], paste(derive_dir, 'sampdata.csv', sep=''), row.names=T)

###########################################################################
### Plot level data

library(sp)
library(raster)
library(geosphere)

# Select columns in plots table to be in plot environmental data table
plot_data = plots[,c('PlotID','Lon','Lat','Year','Ecoregion','PairID','TopoPos','RH','Soil_pH','Soil_ECEC')]

## Cloud frequency (Wilson et al. unpub)
clouds = read.csv('./Data/cloud_cover_Wilson2015.csv')
clouds = clouds[,c('PlotID','Meanannual_CF','Intraannual_SD')]
names(clouds)[2:3] = c('CloudFreq_mean', 'CloudFreq_sd')
plot_data = merge(plot_data, clouds)

## PRISM data: precip, tmax, tmean, vpdmax
prism_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/GIS Data/PRISM/Normals81/'
ap = raster(paste(prism_dir, 'PRISM_ppt_30yr_normal_800mM2_annual_bil/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil', sep=''))
tmax = raster(paste(prism_dir, 'PRISM_tmax_30yr_normal_800mM2_annual_bil/PRISM_tmax_30yr_normal_800mM2_annual_bil.bil', sep=''))
tmean = raster(paste(prism_dir, 'PRISM_tmean_30yr_normal_800mM2_annual_bil/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil', sep=''))
vpdmax = raster(paste(prism_dir, 'PRISM_vpdmax_30yr_normal_800mM2_annual_bil/PRISM_vpdmax_30yr_normal_800mM2_annual_bil.bil', sep=''))
prism = stack(ap, tmax, tmean, vpdmax)
names(prism) = c('AP','T_max','T_mean','VPD_max')
prism_data = extract(prism, plot_data[,c('Lon','Lat')])

plot_data = cbind(plot_data, prism_data)

## NED: 1/9 arc-sec DEM

# Go through each location separately because DEM downloaded as tiles
NED_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/GIS Data/DEM grids/NED/'
filenames = list.files(paste(NED_dir))

elev_data = data.frame()
for(p in plot_data$PlotID){
	# Get data for this plot
	this_plot = subset(plots, PlotID==p)

	# Read in DEM for this location
	loc = substr(p, 1, nchar(p)-1)
	this_dir = filenames[grep(loc, filenames)]
	file_contents = list.files(paste(NED_dir,this_dir, sep=''))
	this_file = file_contents[grep('.img$',file_contents)]
	this_dem = raster(paste(NED_dir, this_dir, '/',this_file, sep=''))

	# Find the center of the plot
	start_pt = this_plot[,c('Lon','Lat')]
	mid_pt = destPoint(start_pt, this_plot$Azimuth_true, 25)
	
	# Extract elevation
	elev = extract(this_dem, mid_pt)

	# Calculate topographic openness

	# Crop DEM to .01-deg around focal point
	use_ext = extent(rep(as.numeric(mid_pt), each=2)+c(-.01,.01,-.01,.01))
	this_dem = crop(this_dem, use_ext)
	
	# Calculate openness from this point
	topo = calc_openness(mid_pt, this_dem, 250)

	# Add to data frame
	elev_data = rbind(elev_data, data.frame(PlotID=p,Elevation=elev, topo))
}

plot_data = merge(plot_data, elev_data)
rownames(plot_data) = plot_data$PlotID
plot_data = plot_data[,-1]

# Save plot data
write.csv(plot_data, paste(derive_dir, 'plot_data.csv', sep=''), row.names=T)

# Plot difference in topographic exposure between paired plots
plot_data$TopoPos = factor(plot_data$TopoPos, levels=c('sheltered','exposed'))

pdf('./Figures/topographic openness of paired plots.pdf', height=4, width=7)
par(mfrow=c(1,2))
plot(OpenPos~as.numeric(TopoPos), data=plot_data, type='n', axes=F, 
	xlab='Topographic Position', ylab='Topographic Openness Index', main='Positive Openness')
use_tab = xtabs(OpenPos~PairID+TopoPos, data=plot_data)
for(i in 1:nrow(use_tab)){
	lines(c(1.25,1.75), use_tab[i,], type='o')	
	text(1.9, use_tab[i,2], rownames(use_tab)[i])
}
axis(1, at=c(1.25,1.75), labels=levels(plot_data$TopoPos))
axis(2, las=1)
plot(OpenNeg~as.numeric(TopoPos), data=plot_data, type='n', axes=F, 
	xlab='Topographic Position', ylab='Topographic Openness Index', main='Negative Openness')
use_tab = xtabs(OpenNeg~PairID+TopoPos, data=plot_data)
for(i in 1:nrow(use_tab)){
	lines(c(1.25,1.75), use_tab[i,], type='o')	
	text(1.9, use_tab[i,2], rownames(use_tab)[i])
}
axis(1, at=c(1.25,1.75), labels=levels(plot_data$TopoPos))
axis(2, las=1)
dev.off()

apply(use_tab, 1, function(x) x[2]-x[1])

###################################################################
### Lichen-level data

## Traits to be analyzed are:
## Level 1 : Attachment, Form, ColorVis, Photobiont, Pseudocyphellae, Maculae, Asco, Asexual 
## Level 2 : Cilia, Rhizines, LobeShape, LobeArea, LobeDissect, CrustForm, AscoForm, AscoCover, AscoAbun, AsexualForm, AsexualAbun

### Recode traits
lichen_data = lichens

# Foliose size and shape
lichen_data$LobeArea = lichens$LobeWidth*lichens$LobeDepth
lichen_data$LobeDissect = lichens$LobeWidth/lichens$LobeDepth

# Sexual reproductive output - note that lirellae will have larger values than apothecia
lichen_data$AscoAbun = ifelse(is.na(lichens$AscoDist)&!is.na(lichens$AscoWidth), 0,
	lichens$AscoWidth/(lichens$AscoWidth+lichens$AscoDist)) 

# Presence of Ascomata
lichen_data$Asco = !(lichen_data$AscoForm=='none')

# AscoCover
lichen_data$AscoCover = NA
lichen_data[lichen_data$AscoForm %in% c('perithecia','recessed apothecia','warts'),'AscoCover'] = 'covered'
lichen_data[lichen_data$AscoForm %in% c('apothecia','lirellae','stalked'),'AscoCover'] = 'exposed'

# Presence of Asexual
lichen_data$Asexual = !(lichen_data$AsexualForm=='none')

## Assign NA to hierarchical trait combinations that are not possible
lichen_data[lichen_data$Form!='foliose',c('Cilia','Rhizines','LobeShape','LobeDissect','LobeArea')] = NA
lichen_data[lichen_data$Form!='crustose', 'CrustForm'] = NA
lichen_data[(lichen_data$Asco==F)|is.na(lichen_data$Asco), c('AscoForm','AscoWidth','AscoHeight','AscoDist','AscoAbun')] = NA
lichen_data[(lichen_data$Asexual==F)|is.na(lichen_data$Asexual),c('AsexualForm', 'AsexualLoc', 'AsexualAbun')] = NA

# Make table for analysis
level1_traits = c('Attachment', 'Form', 'ColorVis', 'Photobiont', 'Pseudocyphellae', 'Maculae', 'Asco', 'Asexual') 
level2_traits = c('Cilia', 'Rhizines', 'LobeShape', 'LobeArea', 'LobeDissect', 'CrustForm', 'AscoForm', 'AscoCover', 'AscoAbun', 'AsexualForm', 'AsexualAbun')

traitdf = data.frame(TraitName=c(level1_traits,level2_traits), 
	Level=c(rep(1, length(level1_traits)), rep(2, length(level2_traits))))

lichen_data = lichen_data[,c('LichenID','SampID','Name','AbunCount',traitdf$TraitName)]
rownames(lichen_data) = lichen_data$LichenID
lichen_data = lichen_data[,-1]

# Write to file- will fill in the other columns manually
#write.csv(traitdf, 'trait_levels.csv', row.names=F)
write.csv(lichen_data, paste(derive_dir, 'lichendata.csv', sep=''), row.names=T)

# Read back in data
lichen_data = read.csv(paste(derive_dir, 'lichendata.csv', sep=''), row.names=1)

# Drop lichen with missing data
lichen_data = subset(lichen_data, Form!='')

## Calculate morphotype X trait matrix

# Drop ColorVis from analysis because trait may not be consistent measure of intrinsic property
use_traits = traitdf$TraitName[traitdf$TraitName!='ColorVis']

# Find unique morphotypes
morphos = unique(lichen_data[,use_traits])
morphos = morphos[order(morphos$Form, morphos$Photobiont, morphos$Asco, morphos$Asexual, morphos$Attachment, morphos$CrustForm),]
morphos$MorphID = paste('M', 1:nrow(morphos), sep='')
rownames(morphos) = morphos$MorphID

# Save morphotypes
write.csv(morphos, paste(derive_dir, 'morphotypes.csv', sep=''), row.names=T)

# Read in morphotypes
morphos = read.csv(paste(derive_dir, 'morphotypes.csv', sep=''), row.names=1)

## Calculate sample X morphotype abundance matrix
lichen_data$MorphID = sapply(1:nrow(lichen_data), function(x){
	 match_morpho(lichen_data[x,], morphos[,names(morphos)!='MorphID'])
})


# Save lichen data
write.csv(lichen_data, paste(derive_dir, 'lichendata.csv', sep=''), row.names=T)

# Read in lichen data
lichen_data = read.csv(paste(derive_dir, 'lichendata.csv', sep=''), row.names=1)

# Convert Abundance counts to number and assign abundance of 1 to missing values
lichen_data$AbunCount = as.numeric(lichen_data$AbunCount)
lichen_data[is.na(lichen_data$AbunCount),'AbunCount'] = 1

# Distinct species with the same morphotye have their abundance added
sampXmorph = xtabs(AbunCount ~ SampID + MorphID, data=lichen_data)

# Save morphotype abundance matrix
write.csv(sampXmorph, paste(derive_dir, 'sampXmorph.csv', sep=''), row.names=T)

morph_freq = colSums(sampXmorph)
plot(morph_freq[order(morph_freq, decreasing=T)])











