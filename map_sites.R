## This script plots a map of lichen FD sites

library(raster)
library(sp)
library(rgdal)
library(lattice)
#library(mapmisc)
library(geosphere)

# Load data
source('C:/Users/jrcoyle/Documents/UNC/Projects/Lichen Functional Diversity/Analysis/GitHub/NC-forest-lichens/load_data.R')

# Load geographic data
elev = raster('../../GIS Data/WorldClim_1km/alt_30s_bil/alt.bil')
nam_outline = readOGR('../../../GIS shape files/N Am Outline', 'na_base_Lambert_Azimuthal')
ecoregions = readOGR('../../GIS Data/NA_CEC_Eco_Level2', 'NA_CEC_Eco_Level2')

# Define projection and re-project shapefiles
plot_prj = "+proj=longlat +units=km"
nam_outline = spTransform(nam_outline, CRS(plot_prj))
ecoregions = spTransform(ecoregions, CRS(plot_prj))

# Define bounding box and crop raster
plot_bbox = matrix(c(-84.5,-75.7, 33.6,36.7), nrow=2, byrow=T)
elev_NC = crop(elev, extent(plot_bbox))
plot(elev_NC)
elev_nc_proj = projectRaster(elev_NC, res=res(elev_NC), crs=CRS(plot_prj))

# Subset polygons
ecoregions_NC = subset(ecoregions, NA_L1NAME=='EASTERN TEMPERATE FORESTS')
nam_outline_NC = subset(nam_outline, PROVINCE_S %in% c('GEORGIA','TENNESSE','SOUTH CAROLINA','NORTH CAROLINA','VIRGINIA'))

# Define site dataset
site_data = subset(plot_data, TopoPos=='exposed')
coordinates(site_data) = c('Lon','Lat')
proj4string(site_data) = CRS(plot_prj)

# Add scale bar
startpt = c(-77.5, 33.7)
midpt = destPoint(startpt, 90, 50000)
endpt = destPoint(startpt, 90, 100000)

# Define ecoregion colors
colfact = c('white','grey50')[factor(site_data$Ecoregion)]

# Make map
pdf('./Figures/site map no labels.pdf', height=1.5, width=3)
trellis.par.set(list(fontsize=list(text=8, points=6)))
spplot(elev_nc_proj, col.regions=colorRampPalette(c('white','grey20')), panel=function(x,y,z,subscripts,...){
	panel.levelplot(x,y,z,subscripts,...)
	sp.polygons(nam_outline_NC, fill='transparent', col='grey35', lwd=1)
	#sp.polygons(ecoregions_NC, fill='transparent', col='black', lwd=1, lty=3, lend=2)
	sp.points(site_data, pch=16, col='white')
	sp.points(site_data, pch=1, col='black')
	#sp.text(coordinates(site_data), site_data$PairID, pos=c(4,4,2,2,2,4,2,2,4))
	panel.rect(startpt[1], startpt[2], endpt[1], startpt[2]+.07, fill='black')
	#sp.text(midpt, '100 km', pos=3)
})
dev.off()

# Make regional inset map

png('./Figures/inset map.png', height=400, width=400)
par(mar=c(0,0,0,0))
plot(nam_outline, xlim=c(-108,-72), ylim=c(25, 70), bg='transparent', fg='grey30', lwd=2)
rect(-84.5, 33.6,-75.7,36.7, lwd=4, col="#00000066")
box()
dev.off()

# SCALE BAR NOT YET WORKING 10/5/2015

