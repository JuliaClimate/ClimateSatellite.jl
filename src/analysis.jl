"""
This file stores the analysis functions of ClimateSatellite.jl that are important for
climatological analysis of the downloaded data.

Current functionalities stored in this file include:
- Yearly, seasonal and monthly means, maximum, minimum, range and standard deviation:
    - On a general climatological basis
    - For a specified hour
    - Of the diurnal variability
- Zonal and meridional averages of the above
- Saving the results of the abovementioned analysis into yearly NetCDF files
- Retrieval of these data from their respective NetCDF files

Functionalities that are in development include:
- Zonal, meridional and domain averages

Because there is progressively more data every year, analysis will be conducted on a yearly
basis instead of for the entire time-series.  This is to ensure generality and ease-of-use
as the datasets become larger and larger.

"""

function clisatanalysis(
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

        davg[:,:,nt,mo] = round(Int16,mean(davg[:,:,1:nt-1,mo],dims=3));
        dstd[:,:,nt,mo] = round(Int16,mean(dstd[:,:,1:nt-1,mo],dims=3));
        dmax[:,:,nt,mo] = round(Int16,maximum(dmax[:,:,1:nt-1,mo],dims=3));
        dmin[:,:,nt,mo] = round(Int16,minimum(dmin[:,:,1:nt-1,mo],dims=3));

        raw = permutedims(raw,(1,2,4,3)); tmp = zeros(Int16,nlon,nlat,ndy);
        for it = 1 : ndy, ilat = 1 : nlat, ilon = 1 : nlon

            tmpii = raw[ilon,ilat,it,:] .!= -32768
            tmp[ilon,ilat,it] = maximum(tmpii)/2 - minimum(tmpii)/2;

        end
        for ilat = 1 : nlat, ilon = 1 : nlon

            rawii = tmp[ilon,ilat,:];
            davg[ilon,ilat,nt+1,mo] = round(Int16,mean(rawii));
            dstd[ilon,ilat,nt+1,mo] = round(Int16,std(rawii));
            dmax[ilon,ilat,nt+1,mo] = round(Int16,maximum(rawii));
            dmin[ilon,ilat,nt+1,mo] = round(Int16,minimum(rawii));

        end

    end

    mean!(davg[:,:,:,end],davg[:,:,:,1:12]); maximum!(dmax[:,:,:,end],dmax[:,:,:,1:12]);
    mean!(dstd[:,:,:,end],dstd[:,:,:,1:12]); maximum!(dmin[:,:,:,end],dmin[:,:,:,1:12]);

    grid = [rlon,rlat]
    clisatanasave([davg,dstd,dmax,dmin],grid,productID,varname,yr,path,region,info);

end

function clisatanasave(
    data::Array{Array{Int16,4},1}, grid::Vector{Any},
    productID::AbstractString, varname::AbstractString, yr::Integer,
    path::AbstractString, region::AbstractString, info::Dict
)

    fol = clisatanafol(productID,Date(yr),region,path=path)
    fnc = joinpath(fol,clisatananame(productID,varname,yr,region));
    rlon,rlat = grid; nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"];

    ds = Dataset(fnc,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat;
    ds.dim["hour"] = nt; ds.dim["month"] = 12;

    att_var = map(x->Dict(),1:2)
    for ii = 1 : 2
        att_var[ii]["units"]         = info["units"];
        att_var[ii]["standard_name"] = info["standard"];
        att_var[ii]["long_name"]     = info["variable"];
        att_var[ii]["add_offset"]    = info["offset"];
        att_var[ii]["missing_value"] = -32768;
    end
    att_var[1]["scale_factor"]  = info["scale"];
    att_var[2]["scale_factor"]  = info["scale"]*2;

    att_lon = Dict("units"=>"degrees_east","long_name"=>"longitude");
    att_lat = Dict("units"=>"degrees_north","long_name"=>"latitude");

    defVar(ds,"longitude",lon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",lat,("latitude",),attrib=att_lat)

    defVar(ds,"domain_yearly_mean_climatology",data[1][:,:,nt+1,end],
           ("longitude","latitude"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_std_climatology",data[2][:,:,nt+1,end],
           ("longitude","latitude"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_maximum_climatology",data[3][:,:,nt+1,end],
           ("longitude","latitude"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_minimum_climatology",data[4][:,:,nt+1,end],
           ("longitude","latitude"),attrib=att_var[1]);

    defVar(ds,"domain_yearly_mean_hourly",data[1][:,:,1:nt,end],
           ("longitude","latitude","hour"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_std_hourly",data[2][:,:,1:nt,end],
           ("longitude","latitude","hour"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_maximum_hourly",data[3][:,:,1:nt,end],
           ("longitude","latitude","hour"),attrib=att_var[1]);
    defVar(ds,"domain_yearly_minimum_hourly",data[4][:,:,1:nt,end],
           ("longitude","latitude","hour"),attrib=att_var[1]);

    defVar(ds,"domain_yearly_mean_diurnalvariance",data[1][:,:,nt+2,end],
           ("longitude","latitude"),attrib=att_var[2]);
    defVar(ds,"domain_yearly_std_diurnalvariance",data[2][:,:,nt+2,end],
           ("longitude","latitude"),attrib=att_var[2]);
    defVar(ds,"domain_yearly_maximum_diurnalvariance",data[3][:,:,nt+2,end],
           ("longitude","latitude"),attrib=att_var[2]);
    defVar(ds,"domain_yearly_minimum_diurnalvariance",data[4][:,:,nt+2,end],
           ("longitude","latitude"),attrib=att_var[2]);

    defVar(ds,"domain_monthly_mean_climatology",data[1][:,:,nt+1,1:12],
           ("longitude","latitude","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_std_climatology",data[2][:,:,nt+1,1:12],
           ("longitude","latitude","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_maximum_climatology",data[3][:,:,nt+1,1:12],
           ("longitude","latitude","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_minimum_climatology",data[4][:,:,nt+1,1:12],
           ("longitude","latitude","month"),attrib=att_var[1]);

    defVar(ds,"domain_monthly_mean_hourly",data[1][:,:,1:nt,1:12],
           ("longitude","latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_std_hourly",data[2][:,:,1:nt,1:12],
           ("longitude","latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_maximum_hourly",data[3][:,:,1:nt,1:12],
           ("longitude","latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"domain_monthly_minimum_hourly",data[4][:,:,1:nt,1:12],
           ("longitude","latitude","hour","month"),attrib=att_var[1]);

    defVar(ds,"domain_monthly_mean_diurnalvariance",data[1][:,:,nt+2,1:12],
           ("longitude","latitude","month"),attrib=att_var[2]);
    defVar(ds,"domain_monthly_std_diurnalvariance",data[2][:,:,nt+2,1:12],
           ("longitude","latitude","month"),attrib=att_var[2]);
    defVar(ds,"domain_monthly_maximum_diurnalvariance",data[3][:,:,nt+2,1:12],
           ("longitude","latitude","month"),attrib=att_var[2]);
    defVar(ds,"domain_monthly_minimum_diurnalvariance",data[4][:,:,nt+2,1:12],
           ("longitude","latitude","month"),attrib=att_var[2]);

    defVar(ds,"zonalavg_yearly_mean_climatology",
           dropdims(mean(data[1][:,:,nt+1,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_std_climatology",
           dropdims(mean(data[2][:,:,nt+1,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_maximum_climatology",
           dropdims(mean(data[3][:,:,nt+1,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_minimum_climatology",
           dropdims(mean(data[4][:,:,nt+1,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[1]);

    defVar(ds,"zonalavg_yearly_mean_hourly",
           dropdims(mean(data[1][:,:,1:nt,end],dims=1),dims=1),
           ("latitude","hour"),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_std_hourly",
           dropdims(mean(data[2][:,:,1:nt,end],dims=1),dims=1),
           ("latitude","hour"),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_maximum_hourly",
           dropdims(mean(data[3][:,:,1:nt,end],dims=1),dims=1),
           ("latitude","hour"),attrib=att_var[1]);
    defVar(ds,"zonalavg_yearly_minimum_hourly",
           dropdims(mean(data[4][:,:,1:nt,end],dims=1),dims=1),
           ("latitude","hour"),attrib=att_var[1]);

    defVar(ds,"zonalavg_yearly_mean_diurnalvariance",
           dropdims(mean(data[1][:,:,nt+2,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[2]);
    defVar(ds,"zonalavg_yearly_std_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[2]);
    defVar(ds,"zonalavg_yearly_maximum_diurnalvariance",
           dropdims(mean(data[3][:,:,nt+2,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[2]);
    defVar(ds,"zonalavg_yearly_minimum_diurnalvariance",
           dropdims(mean(data[4][:,:,nt+2,end],dims=1),dims=1),
           ("latitude",),attrib=att_var[2]);

    defVar(ds,"zonalavg_monthly_mean_climatology",
           dropdims(mean(data[1][:,:,nt+1,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_std_climatology",
           dropdims(mean(data[2][:,:,nt+1,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_maximum_climatology",
           dropdims(mean(data[3][:,:,nt+1,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_minimum_climatology",
           dropdims(mean(data[4][:,:,nt+1,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[1]);

    defVar(ds,"zonalavg_monthly_mean_hourly",
           dropdims(mean(data[1][:,:,1:nt,1:12],dims=1),dims=1),
           ("latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_std_hourly",
           dropdims(mean(data[2][:,:,1:nt,1:12],dims=1),dims=1),
           ("latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_maximum_hourly",
           dropdims(mean(data[3][:,:,1:nt,1:12],dims=1),dims=1),
           ("latitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"zonalavg_monthly_minimum_hourly",
           dropdims(mean(data[4][:,:,1:nt,1:12],dims=1),dims=1),
           ("latitude","hour","month"),attrib=att_var[1]);

    defVar(ds,"zonalavg_monthly_mean_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[2]);
    defVar(ds,"zonalavg_monthly_std_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[2]);
    defVar(ds,"zonalavg_monthly_maximum_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[2]);
    defVar(ds,"zonalavg_monthly_minimum_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=1),dims=1),
           ("latitude","month"),attrib=att_var[2]);

    defVar(ds,"meridionalavg_yearly_mean_climatology",
           dropdims(mean(data[1][:,:,nt+1,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_std_climatology",
           dropdims(mean(data[2][:,:,nt+1,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_maximum_climatology",
           dropdims(mean(data[3][:,:,nt+1,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_minimum_climatology",
           dropdims(mean(data[4][:,:,nt+1,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[1]);

    defVar(ds,"meridionalavg_yearly_mean_hourly",
           dropdims(mean(data[1][:,:,1:nt,end],dims=2),dims=2),
           ("longitude","hour"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_std_hourly",
           dropdims(mean(data[2][:,:,1:nt,end],dims=2),dims=2),
           ("longitude","hour"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_maximum_hourly",
           dropdims(mean(data[3][:,:,1:nt,end],dims=2),dims=2),
           ("longitude","hour"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_yearly_minimum_hourly",
           dropdims(mean(data[4][:,:,1:nt,end],dims=2),dims=2),
           ("longitude","hour"),attrib=att_var[1]);

    defVar(ds,"meridionalavg_yearly_mean_diurnalvariance",
           dropdims(mean(data[1][:,:,nt+2,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[2]);
    defVar(ds,"meridionalavg_yearly_std_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[2]);
    defVar(ds,"meridionalavg_yearly_maximum_diurnalvariance",
           dropdims(mean(data[3][:,:,nt+2,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[2]);
    defVar(ds,"meridionalavg_yearly_minimum_diurnalvariance",
           dropdims(mean(data[4][:,:,nt+2,end],dims=2),dims=2),
           ("longitude",),attrib=att_var[2]);

    defVar(ds,"meridionalavg_monthly_mean_climatology",
           dropdims(mean(data[1][:,:,nt+1,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_std_climatology",
           dropdims(mean(data[2][:,:,nt+1,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_maximum_climatology",
           dropdims(mean(data[3][:,:,nt+1,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_minimum_climatology",
           dropdims(mean(data[4][:,:,nt+1,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[1]);

    defVar(ds,"meridionalavg_monthly_mean_hourly",
           dropdims(mean(data[1][:,:,1:nt,1:12],dims=2),dims=2),
           ("longitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_std_hourly",
           dropdims(mean(data[2][:,:,1:nt,1:12],dims=2),dims=2),
           ("longitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_maximum_hourly",
           dropdims(mean(data[3][:,:,1:nt,1:12],dims=2),dims=2),
           ("longitude","hour","month"),attrib=att_var[1]);
    defVar(ds,"meridionalavg_monthly_minimum_hourly",
           dropdims(mean(data[4][:,:,1:nt,1:12],dims=2),dims=2),
           ("longitude","hour","month"),attrib=att_var[1]);

    defVar(ds,"meridionalavg_monthly_mean_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[2]);
    defVar(ds,"meridionalavg_monthly_std_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[2]);
    defVar(ds,"meridionalavg_monthly_maximum_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[2]);
    defVar(ds,"meridionalavg_monthly_minimum_diurnalvariance",
           dropdims(mean(data[2][:,:,nt+2,1:12],dims=2),dims=2),
           ("longitude","month"),attrib=att_var[2]);

    close(ds);

end
