## This script calculates a null distribution of rao's Q for samples in the Lichen FD project

options(stringsAsFactors=F)

########################################
### Functions, Libraries and Data

library(cluster) # daisy
library(vegan) # nullmodel

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


# Sample X Morphotype Abundance matrix
sampXmorph = read.csv('sampXmorph.csv', row.names=1)

# Morphotype X Trait matrix
morphos = read.csv('morphotypes.csv', row.names=1)

# Trait info
traitdf = read.csv('trait_levels.csv')

# Location information
tree_data = read.csv('treedata.csv', row.names=1)
samp_data = read.csv('sampdata.csv', row.names=1)
plot_data = read.csv('plot_data.csv', row.names=1)

# Add IDs
IDlist = strsplit(rownames(samp_data), '-')
samp_data$PlotID = sapply(IDlist, function(x) x[1])
samp_data$SampID = rownames(samp_data)

##########################################################
### Script

morphos = morphos[,colnames(morphos)!='MorphID']

## Subset data to only sites and samples of interest
top2 = sapply(rownames(tree_data), function(x){
	these_samps = subset(samp_data, TreeID==x)
	keep_samps = rownames(these_samps[order(these_samps$Height,decreasing=T),])[1:2]		
	keep_samps
})
top2 = as.vector(top2)

# only use top2 samples with lichens
top2 = top2[top2 %in% rownames(sampXmorph)]
sampXmorph = sampXmorph[top2,]

# remove morphotypes not in top2 samples
keep_morphs = names(colSums(sampXmorph)!=0)
morphos = morphos[keep_morphs,]
sampXmorph = sampXmorph[,keep_morphs]

## Put morphotype traits into correct format
morphos$Attachment = factor(morphos$Attachment, ordered=T)
morphos$Form = factor(morphos$Form)
morphos$Photobiont = factor(morphos$Photobiont)
morphos$Pseudocyphellae = factor(morphos$Pseudocyphellae)
morphos$Maculae = factor(morphos$Maculae)
morphos$Asco = factor(morphos$Asco)
morphos$Asexual = factor(morphos$Asexual)
morphos$Cilia = factor(morphos$Cilia)
morphos$Rhizines = factor(morphos$Rhizines)
morphos$LobeShape = factor(morphos$LobeShape, levels=c('concave','flat','convex'), ordered=T)
morphos$CrustForm = factor(morphos$CrustForm)
morphos$AscoForm = factor(morphos$AscoForm)
morphos$AscoCover = factor(morphos$AscoCover)
morphos$AsexualForm = factor(morphos$AsexualForm)
morphos$AsexualAbun = factor(morphos$AsexualAbun, levels=c('few','several','many','covered'), ordered=T)

## Calculate distance matrix among morphotypes

# Scaling factor - amount by which 2ndary traits affect distance relative to 1st traits
sfact = .5

use_traits = traitdf[traitdf$TraitName %in% colnames(morphos),]
rownames(use_traits) = use_traits$TraitName

# weights are (sfact)^level-1
dist_L2 = daisy(morphos, type=list(asymm=use_traits[use_traits$type=='asymm','TraitName'],
					symm=use_traits[use_traits$type=='symm','TraitName']),
		weights=sfact^(use_traits[colnames(morphos),'Level']-1)
)

dist_L2mat = as.matrix(dist_L2)

## Generate random communities

# Put matrices in same order
sampXmorph = sampXmorph[,colnames(dist_L2mat)]

# Choose null model swap algorithm
# Generate null matrices: Shuffle species across samples in the same site but between paired plots, keeping richness and species abundance constant
# swsh_samp preserves morphotype frequencies across and within samples (e.g richness constant)
# swsh_samp_c preserves morphotype abundance distribution
# swsh_samp_r preserves abundance across samples
use_null_alg = 'swsh_samp_c'

N=10000 # Must be divisible by 1000

null_rao_mat = matrix(NA, ncol=nrow(sampXmorph), nrow=N)
colnames(null_rao_mat) = rownames(sampXmorph)

this_run = paste('acrossall', sfact, use_null_alg, sep='-')

x = sampXmorph
	
# Drop absent morphotypes to speed computation
x = x[,colSums(x)!=0]

# Subset morphotype distance matrix
use_dmat = dist_L2mat[colnames(x),colnames(x)]

# Re-start run
start_j = 301
null_rao_mat = read.csv(paste('null_rao-', this_run, '.csv', sep=''), check.names=F)


# This reduces memory requirements of the script
for(j in start_j:N){
	nulls = simulate(nullmodel(x, use_null_alg), 1)

	# Calculate rao's Q on each matrix
	null_scaled = as.matrix(nulls[,,1]/rowSums(nulls[,,1]))
	null_raos = calc_rao(null_scaled, use_dmat)

	null_rao_mat[j,] = null_raos
	
	# Save current progress
	if(j %% 50 == 0) write.csv(null_rao_mat, paste('null_rao-', this_run, '.csv', sep=''), row.names=F)
	if(j %% 50 == 0) print(paste('Done with',j*N/1000))
		
	gc()

}

# Write null model data to file
write.csv(null_rao_mat, paste('null_rao-', this_run, '.csv', sep=''), row.names=F)


quit('no')

