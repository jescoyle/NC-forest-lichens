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
sampXmorph = read.csv('sampXrepM.csv', row.names=1)

# Morphotype X Trait matrix
morphos = read.csv('reproductive_morphotypes.csv', row.names=1)

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
morphos$Asco = factor(morphos$Asco)
morphos$Asexual = factor(morphos$Asexual)
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
# Generate null matrices: Shuffle species across samples in the same plot, keeping richness and species abundance constant
# swsh_samp preserves morphotype frequencies across and within samples (e.g richness constant)
# swsh_samp_c preserves morphotype abundance distribution
# swsh_samp_r preserves abundance across samples
use_null_alg = 'swsh_samp_r'

N=10000 # must be divisible by 100

null_rao_mat = matrix(NA, ncol=nrow(sampXmorph), nrow=N)
colnames(null_rao_mat) = rownames(sampXmorph)

this_run = paste('withinplot', sfact, use_null_alg, sep='-')

for(this_plot in rownames(plot_data)){
	# Define community matrix for this plot
	these_samps = subset(samp_data, PlotID==this_plot)$SampID
	these_samps = these_samps[these_samps %in% rownames(sampXmorph)]
	x = sampXmorph[these_samps,]
	
	# Drop absent morphotypes to speed computation
	x = x[,colSums(x)!=0]


	# Subset morphotype distance matrix
	use_dmat = dist_L2mat[colnames(x),colnames(x)]

	for(j in 1:100){
		nulls = simulate(nullmodel(x, use_null_alg), N/100)

		# Calculate rao's Q on each matrix
		null_raos = sapply(1:(N/100), function(i){
			null_scaled = as.matrix(nulls[,,i]/rowSums(nulls[,,i]))
			calc_rao(null_scaled, use_dmat)
		})

		null_rao_mat[(((j-1)*(N/100))+1):(j*N/100),these_samps] = null_raos
	
		# Save current progress
		write.csv(null_rao_mat, paste('null_rep_rao-', this_run, '.csv', sep=''), row.names=F)
		if(j %% 10 == 0) print(paste('Done with',this_plot, j*N/100))
		
		gc()

	}

}

# Write null model data to file

write.csv(null_rao_mat, paste('null_rep_rao-', this_run, '.csv', sep=''), row.names=F)


quit('no')

