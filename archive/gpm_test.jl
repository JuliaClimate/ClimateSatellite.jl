using ClimateSatellite
using Dates
using PyCall,PyPlot

cd("/Users/natgeo-wong")
gpmrun(Date(2012,1,1),clisatroot(),["ENSO"])

cd("/Users/natgeo-wong/Documents/research/data/GPM/ENSO/2012/01/");
prcp_ENSO = ncread("gpm_ENSO_prcp_20120101.nc","prcp");
lon_ENSO = ncread("gpm_ENSO_prcp_20120101.nc","lon");
lat_ENSO = ncread("gpm_ENSO_prcp_20120101.nc","lat");

close("all")
ccrs = pyimport("cartopy.crs");
ax = subplot(projection=ccrs.Mollweide(central_longitude=200));
#ax.set_global()

lvls = convert(Array,0:5:80);
contourf(lon_ENSO,lat_ENSO,transpose(prcp_ENSO[:,:,5]),lvls,transform=ccrs.PlateCarree())

ax.coastlines()
gcf()
