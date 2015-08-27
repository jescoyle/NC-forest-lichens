
options(stringsAsFactors=F)

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/'
derive_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/Derived Tables/'
bugs_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/BUGS/'
git_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/'

setwd(working_dir)

## Read in Functions
source(paste(git_dir, 'lichen_FD_functions.R', sep=''))

## Read in Data

# Morphotype X Trait matrix
morphos = read.csv(paste(derive_dir, 'morphotypes.csv', sep=''), row.names=1)

# Sample X Morphotype Abundance matrix
sampXmorph = read.csv(paste(derive_dir, 'sampXmorph.csv', sep=''), row.names=1)

# Lichen data
lichen_data = read.csv(paste(derive_dir, 'lichendata.csv', sep=''), row.names=1)

# Location information
tree_data = read.csv(paste(derive_dir, 'treedata.csv', sep=''), row.names=1)
samp_data = read.csv(paste(derive_dir, 'sampdata.csv', sep=''), row.names=1)
plot_data = read.csv(paste(derive_dir, 'plot_data.csv', sep=''), row.names=1)

# Trait info
traitdf = read.csv('trait_levels.csv')

## Make data frame for analysis
lichens = lichen_data

# Add plot and sample level IDs to lichen data
splitlist = strsplit(lichen_data$SampID, '-')
lichens$TreeID = sapply(splitlist, function(x) paste(x[1],x[2], sep='-'))
lichens$PlotID = sapply(splitlist, function(x) x[1])
lichens$SiteID = substr(lichens$PlotID, 1, nchar(lichens$PlotID)-1)

# Merge with sample, tree and plot level data
samp_data$SampID = rownames(samp_data)
tree_data$TreeID = rownames(tree_data)
plot_data$PlotID = rownames(plot_data)
plot_data$SiteID = substr(plot_data$PlotID, 1, nchar(plot_data$PlotID)-1)

# Make sample data frame for analysis
samples = samp_data
samples = merge(samples, tree_data)
samples = merge(samples, plot_data)
rownames(samples) = samples$SampID

## Subset to to just top 2 samples on each tree
top2 = sapply(rownames(tree_data), function(x){
	these_samps = subset(samp_data, TreeID==x)
	keep_samps = rownames(these_samps[order(these_samps$Height,decreasing=T),])[1:2]		
	keep_samps
})
top2 = as.vector(top2)

samples = subset(samples, SampID %in% top2)
lichens = subset(lichens, SampID %in% top2)


# For morphotype data, only use top2 samples with lichens
morphos = morphos[,colnames(morphos)!='MorphID']
use_top2 = top2[top2 %in% rownames(sampXmorph)]
sampXmorph = sampXmorph[use_top2,]

# Remove morphotypes not in top2 samples
keep_morphs = names(colSums(sampXmorph)!=0)
morphos = morphos[keep_morphs,]
sampXmorph = sampXmorph[,keep_morphs]


# Convert categorical variables to factors
samples$Bryophytes = factor(samples$Bryophytes, levels=c('none','minute','few','several','many','covered'), ordered=T)
samples$Shedding = factor(samples$Shedding, levels=c('none','attached flakes','loose flakes','slight peeling','peeling'), ordered=T)
samples$Ecoregion = factor(samples$Ecoregion, levels=c('Mountains','Piedmont','Coastal Plain'))
samples$TopoPos = factor(samples$TopoPos, levels=c('sheltered','exposed'))


# Assign levels to plotID that will make them plot in the correct order
site_order=c('Bladen','Yadkin','John','Eno','Hang','Pisgah','New','High','Jeff')
plot_data$SiteID = factor(plot_data$SiteID, levels=site_order)
ordered_plots = plot_data[order(plot_data$SiteID, plot_data$PairID, plot_data$TopoPos),'PlotID']


# Convert group IDs to factors
lichens$SampID = factor(lichens$SampID)
lichens$TreeID = factor(lichens$TreeID)
lichens$PlotID = factor(lichens$PlotID, levels=ordered_plots)
lichens$SiteID = factor(lichens$SiteID)
samples$TreeID = factor(samples$TreeID)
samples$PlotID = factor(samples$PlotID, levels=ordered_plots)
samples$SiteID = factor(samples$SiteID)

# Merge lichen and sample data
lichens = merge(lichens, samples, all.x=T)

# Change any missing abundance to 1
lichens[is.na(lichens$AbunCount),'AbunCount'] = 1

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



