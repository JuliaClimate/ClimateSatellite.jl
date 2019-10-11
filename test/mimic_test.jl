using ClimateSatellite
using Dates
using PyCall,PyPlot

cd("/Users/natgeo-wong")
#mimicrun(Date(2019,1,1),clisatroot(),["ENSO"])

cd("/Users/natgeo-wong/Documents/Research/MIMIC/ENSO/2019/01/");
tpw_ENSO = ncread("mimic_ENSO_tpw_20190101.nc","tpw");
lon_ENSO = ncread("mimic_ENSO_tpw_20190101.nc","lon");
lat_ENSO = ncread("mimic_ENSO_tpw_20190101.nc","lat");

close("all")
ccrs = pyimport("cartopy.crs");
ax = subplot(projection=ccrs.Mollweide(central_longitude=200));
#ax.set_global()

lvls = convert(Array,0:5:80);
contourf(lon_ENSO,lat_ENSO,transpose(tpw_ENSO[:,:,10]),lvls,transform=ccrs.PlateCarree())

ax.coastlines()
gcf()
