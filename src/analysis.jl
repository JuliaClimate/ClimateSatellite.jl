"""
This file stores the analysis functions of ClimateSatellite.jl that are important for
climatological analysis of the downloaded data.

Current functionalities stored in this file include:
- Yearly, seasonal and monthly means, maximum, minimum, range and standard deviation:
    - On a general climatological basis
    - For a specified hour
    - Of the diurnal variability
- Saving the results of the abovementioned analysis into yearly NetCDF files
- Retrieval of these data from their respective NetCDF files

Functionalities that are in development include:
- Zonal, meridional and domain averages

Because there is progressively more data every year, analysis will be conducted on a yearly
basis instead of for the entire time-series.  This is to ensure generality and ease-of-use
as the datasets become larger and larger.

"""

function clisatanaspatial(
    productID::AbstractString, yr::Integer;
    varname::AbstractString,
    path::AbstractString="", region::AbstractString="GLB"
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end
    if !isdir(clisatfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    rawfol = clisatrawfol(productID,Date(yr),region,path=dataroot);
    anafol = clisatanafol(productID,Date(yr),region,path=dataroot);
    if !isdir(rawfol)
        error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);
    clisatvarinfo!(info,productID,varname=varname); cd(rawfol);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the entire $(regionfullname(region)) region ..."

    lon,lat = clisatlonlat(info); rlon,rlat,rinfo = regiongridvec(region,lon,lat);
    nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"]+1;

    davg = zeros(Int16,nlon,nlat,nt+1,13); dstd = zeros(Int16,nlon,nlat,nt+1,13);
    dmax = zeros(Int16,nlon,nlat,nt+1,13); dmin = zeros(Int16,nlon,nlat,nt+1,13);

    for mo = 1 : 12; ndy = daysinmonth(yr,mo)

        ncraw = clisatrawname(productID,Date(yr,mo),region);
        ds = Dataset(ncraw,"r"); vds = ds[varname];
        raw = reshape(vds.var,nlon,nlat,(nt-1),ndy);

        for it = 1 : nt-1, ilat = 1 : nlat, ilon = 1 : nlon

            rawii = Float64(raw[ilon,ilat,it,:] .!= -32768);
            davg[ilon,ilat,it,mo] = round(Int16,mean(rawii));
            dstd[ilon,ilat,it,mo] = round(Int16,std(rawii));
            dmax[ilon,ilat,it,mo] = round(Int16,maximum(rawii));
            dmin[ilon,ilat,it,mo] = round(Int16,minimum(rawii));

        end

        davg[:,:,nt,mo] = round(Int16,mean(davg[:,:,1:nt-1,mo]));
        dstd[:,:,nt,mo] = round(Int16,mean(dstd[:,:,1:nt-1,mo]));
        dmax[:,:,nt,mo] = round(Int16,maximum(dmax[:,:,1:nt-1,mo]));
        dmin[:,:,nt,mo] = round(Int16,minimum(dmin[:,:,1:nt-1,mo]));

        raw = permutedims(raw,(1,2,4,3)); raw[raw.==-32768] .= NaN;
        varraw = dropdims(Int32(maximum(raw,dims=4))-Int32(minimum(raw,dims=4)),dims=4);
        for ilat = 1 : nlat, ilon = 1 : nlon

            rawii = Float64(varraw[ilon,ilat,:] .!= 0);
            vavg[ilon,ilat,nt+1,mo] = round(Int16,mean(rawii));
            vstd[ilon,ilat,nt+1,mo] = round(Int16,std(rawii));
            vmax[ilon,ilat,nt+1,mo] = round(Int16,maximum(rawii));
            vmin[ilon,ilat,nt+1,mo] = round(Int16,minimum(rawii));

        end

    end

    mean!(davg[:,:,:,end],davg[:,:,:,1:12]); maximum!(dmax[:,:,:,end],dmax[:,:,:,1:12]);
    mean!(dstd[:,:,:,end],dstd[:,:,:,1:12]); maximum!(dmin[:,:,:,end],dmin[:,:,:,1:12]);

    clisatanasave()

end
