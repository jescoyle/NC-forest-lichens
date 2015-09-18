## This script combines filed output by CIMES into one summary csv for all trees
## Generates the file: canopy_summary.csv
library(stringr)

options(stringsAsFactors=F)
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/'
sql_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/SQLite Tables/'

CIMES_dir = 'Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/Data/Canopy Photos/CIMES/Data/'

setwd(working_dir)

# Load data
trees = read.csv(paste(sql_dir,'trees.csv', sep=''))
plots = read.csv(paste(sql_dir,'plots.csv',sep=''))

# Generate list of radiation files
combos = expand.grid(paste('T',1:40, sep=''), plots$PlotID)
filenames = apply(combos, 1, function(x) paste('rad_',x[2],'_',x[1],'.txt', sep=''))


# Create dataframe to append data to
rad_dat = data.frame()

# Note: Bladen1 T21 is missing photo and data so will need to re-start loop after it hits this error
for(i in 742:nrow(combos)){
	f = filenames[i]
	dat = readLines(paste(CIMES_dir,f, sep=''))

	# Get Julian Days on which radiation data were calculated
	dayline = grep('days', dat)
	days = na.omit(as.numeric(unlist(strsplit(dat[dayline], '[ \t]+'))))

	# Get total PPFD transmitted through canopy on each day (mol -m2 -d)
	totlines = grep('Total', dat)
	transline = grep('TRANSMITTED', dat)
	totline = min(totlines[(totlines-transline)>0])
	tots = na.omit(as.numeric(unlist(strsplit(dat[totline], '[ \t]+'))))
	
	# Integrate rates at each day over period of time between sampled days (i.e. monthly)
	if(days[length(days)]!=365) days = c(days, 365)
	periods = days[2:length(days)] - days[1:(length(days)-1)]
	Trans_tot = sum(tots*periods)

	# Get canopy openness
	opline =grep('CANOPY OPENNESS', dat)+1
	Canopy_open = as.numeric(str_trim(dat[opline]))
	
	# Get plot and tree IDs
	PlotID = combos[i,2]
	TreeID = paste(PlotID, substring(combos[i,1], 2), sep='-')

	rad_dat = rbind(rad_dat, data.frame(PlotID, TreeID, Trans_tot, Canopy_open))

}

# Convert back to rate 
rad_dat$Trans_tot = rad_dat$Trans_tot/365

# Save data
write.csv(rad_dat, './Data/canopy_summary_CIMES.csv', row.names=F)






