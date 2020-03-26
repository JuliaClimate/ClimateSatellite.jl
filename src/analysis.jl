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

function clisatopNaN(f::Function,data::Vector{Int16})
    ndata = @view data[data.!=-32768];
    if ndata != []; return f(ndata); else; return -32768; end
end

function clisatanalysis(
    productID::AbstractString, yr::Integer;
    varname::AbstractString,
    path::AbstractString="", region::AbstractString="GLB"
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end
    isgeoregion(region);

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);
    clisatvarinfo!(info,productID,varname=varname);

    regfol = clisatregfol(info,region);
    rawfol = clisatrawfol(info,Date(yr),region);
    cd(rawfol);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data in $(gregionfullname(region)) region during $yr ..."

    lon,lat = clisatlonlat(info); rlon,rlat,rinfo = gregiongridvec(region,lon,lat);
    nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"]+1; grid = [rlon,rlat];

    davg = rand(Int16,nlon,nlat,nt+1,13); dstd = rand(Int16,nlon,nlat,nt+1,13);
    dmax = rand(Int16,nlon,nlat,nt+1,13); dmin = rand(Int16,nlon,nlat,nt+1,13);

    zavg = rand(Int16,nlat,nt+1,13); zstd = rand(Int16,nlat,nt+1,13);
    zmax = rand(Int16,nlat,nt+1,13); zmin = rand(Int16,nlat,nt+1,13);

    mavg = rand(Int16,nlon,nt+1,13); mstd = rand(Int16,nlon,nt+1,13);
    mmax = rand(Int16,nlon,nt+1,13); mmin = rand(Int16,nlon,nt+1,13);

    for mo = 1 : 12; ndy = daysinmonth(yr,mo)

        @info "$(Dates.now()) - Analyzing $(info["source"]) $(info["product"]) data in $(gregionfullname(region)) during $(Dates.monthname(mo)) $yr ..."
        ncraw = clisatrawname(productID,Date(yr,mo),region);
        ds = Dataset(ncraw,"r"); vds = ds[varname];
        raw = reshape(vds.var[:],nlon,nlat,(nt-1),ndy);

        @debug "$(Dates.now()) - Extracting hourly information for each month ..."
        for it = 1 : nt-1, ilat = 1 : nlat, ilon = 1 : nlon

            rawii = @view raw[ilon,ilat,it,:]; rawii = @view rawii[rawii.!=-32768];

            if rawii != [];
                davg[ilon,ilat,it,mo] = round(Int16,mean(rawii));
                dstd[ilon,ilat,it,mo] = round(Int16,std(rawii));
                dmax[ilon,ilat,it,mo] = round(Int16,maximum(rawii));
                dmin[ilon,ilat,it,mo] = round(Int16,minimum(rawii));
            else
                davg[ilon,ilat,it,mo] = -32768; dstd[ilon,ilat,it,mo] = -32768;
                dmax[ilon,ilat,it,mo] = -32768; dmin[ilon,ilat,it,mo] = -32768;
            end

        end

        @debug "$(Dates.now()) - Permuting days and hours dimensions ..."
        raw = permutedims(raw,(1,2,4,3)); drg = zeros(nlon,nlat,ndy);
        dmn = zeros(nlon,nlat,ndy);

        @debug "$(Dates.now()) - Averaging diurnal data into daily data ..."
        for idy = 1 : ndy, ilat = 1 : nlat, ilon = 1 : nlon

            tmpii = @view raw[ilon,ilat,idy,:]; tmpii = @view tmpii[tmpii.!=-32768];

            if tmpii != [];
                drg[ilon,ilat,idy] = maximum(tmpii)/2 - minimum(tmpii)/2;
                dmn[ilon,ilat,idy] = mean(tmpii);
            else
                drg[ilon,ilat,idy] = -32768;
                dmn[ilon,ilat,idy] = -32768;
            end

        end

        @debug "$(Dates.now()) - Extracting information on monthly climatology and diurnal variability ..."
        for ilat = 1 : nlat, ilon = 1 : nlon;

            dmnii = @view dmn[ilon,ilat,:]; dmnii = @view dmnii[dmnii.!=-32768];
            drgii = @view drg[ilon,ilat,:]; drgii = @view drgii[drgii.!=-32768];

            davg[ilon,ilat,nt,mo] = round.(Int16,mean(dmnii));
            dstd[ilon,ilat,nt,mo] = round.(Int16,std(dmnii));
            dmax[ilon,ilat,nt,mo] = round.(Int16,maximum(dmnii));
            dmin[ilon,ilat,nt,mo] = round.(Int16,minimum(dmnii));


            davg[ilon,ilat,nt+1,mo] = round.(Int16,mean(drgii));
            dstd[ilon,ilat,nt+1,mo] = round.(Int16,std(drgii));
            dmax[ilon,ilat,nt+1,mo] = round.(Int16,maximum(drgii));
            dmin[ilon,ilat,nt+1,mo] = round.(Int16,minimum(drgii));

        end

    end

    @info "$(Dates.now()) - Calculating yearly climatology for $(info["source"]) $(info["product"]) in $(gregionfullname(region)) during $yr ..."
    for it = 1 : nt+1, ilat = 1 : nlat, ilon = 1 : nlon

        davg[ilon,ilat,it,end] = round.(Int16,clisatopNaN(mean,davg[ilon,ilat,it,1:12]));
        dstd[ilon,ilat,it,end] = round.(Int16,clisatopNaN(mean,dstd[ilon,ilat,it,1:12]));
        dmax[ilon,ilat,it,end] = round.(Int16,clisatopNaN(mean,dmax[ilon,ilat,it,1:12]));
        dmin[ilon,ilat,it,end] = round.(Int16,clisatopNaN(mean,dmin[ilon,ilat,it,1:12]));

    end

    @info "$(Dates.now()) - Calculating zonal-averaged climatology for $(info["source"]) $(info["product"]) in $(gregionfullname(region)) during $yr ..."
    for ilat = 1 : nlat, it = 1 : nt+1, imo = 1 : 13
        zavg[ilat,it,imo] = round.(Int16,clisatopNaN(mean,davg[:,ilat,it,imo]));
        zstd[ilat,it,imo] = round.(Int16,clisatopNaN(mean,dstd[:,ilat,it,imo]));
        zmax[ilat,it,imo] = round.(Int16,clisatopNaN(mean,dmax[:,ilat,it,imo]));
        zmin[ilat,it,imo] = round.(Int16,clisatopNaN(mean,dmin[:,ilat,it,imo]));
    end

    @info "$(Dates.now()) - Calculating meridional-averaged climatology for $(info["source"]) $(info["product"]) in $(gregionfullname(region)) during $yr ..."
    for imo = 1 : 13, it = 1 : nt+1, ilon = 1 : nlon
        mavg[ilon,it,imo] = round.(Int16,clisatopNaN(mean,davg[ilon,:,it,imo]));
        mstd[ilon,it,imo] = round.(Int16,clisatopNaN(mean,dstd[ilon,:,it,imo]));
        mmax[ilon,it,imo] = round.(Int16,clisatopNaN(mean,dmax[ilon,:,it,imo]));
        mmin[ilon,it,imo] = round.(Int16,clisatopNaN(mean,dmin[ilon,:,it,imo]));
    end

    clisatanasave([davg,dstd,dmax,dmin],[zavg,zstd,zmax,zmin],[mavg,mstd,mmax,mmin],
                  grid,productID,varname,yr,path,region,info);

end

function clisatanasave(
    data::Array{Array{Int16,4},1},
    zdata::Array{Array{Int16,3},1},
    mdata::Array{Array{Int16,3},1},
    grid::Vector{<:Any},
    productID::AbstractString, varname::AbstractString, yr::Integer,
    path::AbstractString, region::AbstractString, info::Dict
)

    @info "$(Dates.now()) - Saving analysed $(info["source"]) $(info["product"]) data in $(gregionfullname(region)) for the year $yr ..."

    fol = clisatanafol(info,region)
    fnc = joinpath(fol,clisatananame(productID,varname,Date(yr),region));
    rlon,rlat = grid; nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"];

    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating NetCDF file $(fnc) for analyzed $(info["source"]) $(info["product"]) data in $yr ..."

    ds = Dataset(fnc,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat;
    ds.dim["hour"] = nt; ds.dim["month"] = 12;

    att_var = map(x->Dict(),1:2)
    for ii = 1 : 2
        att_var[ii]["units"]         = info["units"];
        att_var[ii]["standard_name"] = info["standard"];
        att_var[ii]["long_name"]     = info["variable"];
        att_var[ii]["missing_value"] = -32768;
    end
    att_var[1]["scale_factor"] = info["scale"];   att_var[1]["add_offset"] = info["offset"];
    att_var[2]["scale_factor"] = info["scale"]*2; att_var[2]["add_offset"] = 0;

    att_lon = Dict("units"=>"degrees_east","long_name"=>"longitude");
    att_lat = Dict("units"=>"degrees_north","long_name"=>"latitude");

    defVar(ds,"longitude",rlon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",rlat,("latitude",),attrib=att_lat)

    @debug "$(Dates.now()) - Saving analyzed $(info["source"]) $(info["product"]) data for $yr to NetCDF file $(fnc) ..."

    v = defVar(ds,"domain_yearly_mean_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_std_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_maximum_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_minimum_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,nt+1,end];


    v = defVar(ds,"domain_yearly_mean_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_std_hourly",Int16,
           ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_maximum_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_minimum_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,1:nt,end];


    v = defVar(ds,"domain_yearly_mean_diurnalvariance",Int16,
           ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[1][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_std_diurnalvariance",Int16,
               ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[2][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_maximum_diurnalvariance",Int16,
               ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[3][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_minimum_diurnalvariance",Int16,
           ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[4][:,:,nt+2,end];


    v = defVar(ds,"domain_monthly_mean_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_std_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_maximum_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_minimum_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,nt+1,1:12];


    v = defVar(ds,"domain_monthly_mean_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_std_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_maximum_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_minimum_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,1:nt,1:12];


    v = defVar(ds,"domain_monthly_mean_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[1][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_std_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[2][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_maximum_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[3][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_minimum_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[4][:,:,nt+2,1:12];


    v = defVar(ds,"zonalavg_yearly_mean_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = zdata[1][:,nt+1,end];

    v = defVar(ds,"zonalavg_yearly_std_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = zdata[2][:,nt+1,end];

    v = defVar(ds,"zonalavg_yearly_maximum_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = zdata[3][:,nt+1,end];

    v = defVar(ds,"zonalavg_yearly_minimum_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = zdata[4][:,nt+1,end];


    v = defVar(ds,"zonalavg_yearly_mean_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = zdata[1][:,1:nt,end];

    v = defVar(ds,"zonalavg_yearly_std_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = zdata[2][:,1:nt,end];

    v = defVar(ds,"zonalavg_yearly_maximum_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = zdata[3][:,1:nt,end];

    v = defVar(ds,"zonalavg_yearly_minimum_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = zdata[4][:,1:nt,end];


    v = defVar(ds,"zonalavg_yearly_mean_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = zdata[1][:,nt+2,end];

    v = defVar(ds,"zonalavg_yearly_std_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = zdata[2][:,nt+2,end];

    v = defVar(ds,"zonalavg_yearly_maximum_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = zdata[3][:,nt+2,end];

    v = defVar(ds,"zonalavg_yearly_minimum_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = zdata[4][:,nt+2,end];


    v = defVar(ds,"zonalavg_monthly_mean_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = zdata[1][:,nt+1,1:12];

    v = defVar(ds,"zonalavg_monthly_std_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = zdata[2][:,nt+1,1:12];

    v = defVar(ds,"zonalavg_monthly_maximum_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = zdata[3][:,nt+1,1:12];

    v = defVar(ds,"zonalavg_monthly_minimum_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = zdata[4][:,nt+1,1:12];


    v = defVar(ds,"zonalavg_monthly_mean_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = zdata[1][:,1:nt,1:12];

    v = defVar(ds,"zonalavg_monthly_std_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = zdata[2][:,1:nt,1:12];

    v = defVar(ds,"zonalavg_monthly_maximum_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = zdata[3][:,1:nt,1:12];

    v = defVar(ds,"zonalavg_monthly_minimum_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = zdata[4][:,1:nt,1:12];


    v = defVar(ds,"zonalavg_monthly_mean_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = zdata[1][:,nt+2,1:12];

    v = defVar(ds,"zonalavg_monthly_std_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = zdata[2][:,nt+2,1:12];

    v = defVar(ds,"zonalavg_monthly_maximum_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = zdata[3][:,nt+2,1:12];

    v = defVar(ds,"zonalavg_monthly_minimum_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = zdata[4][:,nt+2,1:12];


    v = defVar(ds,"meridionalavg_yearly_mean_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = mdata[1][:,nt+1,end];

    v = defVar(ds,"meridionalavg_yearly_std_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = mdata[2][:,nt+1,end];

    v = defVar(ds,"meridionalavg_yearly_maximum_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = mdata[3][:,nt+1,end];

    v = defVar(ds,"meridionalavg_yearly_minimum_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = mdata[4][:,nt+1,end];


    v = defVar(ds,"meridionalavg_yearly_mean_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = mdata[1][:,1:nt,end];

    v = defVar(ds,"meridionalavg_yearly_std_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = mdata[2][:,1:nt,end];

    v = defVar(ds,"meridionalavg_yearly_maximum_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = mdata[3][:,1:nt,end];

    v = defVar(ds,"meridionalavg_yearly_minimum_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = mdata[4][:,1:nt,end];


    v = defVar(ds,"meridionalavg_yearly_mean_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = mdata[1][:,nt+2,end];

    v = defVar(ds,"meridionalavg_yearly_std_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = mdata[2][:,nt+2,end];

    v = defVar(ds,"meridionalavg_yearly_maximum_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = mdata[3][:,nt+2,end];

    v = defVar(ds,"meridionalavg_yearly_minimum_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = mdata[4][:,nt+2,end];


    v = defVar(ds,"meridionalavg_monthly_mean_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = mdata[1][:,nt+1,1:12];

    v = defVar(ds,"meridionalavg_monthly_std_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = mdata[2][:,nt+1,1:12];

    v = defVar(ds,"meridionalavg_monthly_maximum_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = mdata[3][:,nt+1,1:12];

    v = defVar(ds,"meridionalavg_monthly_minimum_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = mdata[4][:,nt+1,1:12];


    v = defVar(ds,"meridionalavg_monthly_mean_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = mdata[1][:,1:nt,1:12];

    v = defVar(ds,"meridionalavg_monthly_std_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = mdata[2][:,1:nt,1:12];

    v = defVar(ds,"meridionalavg_monthly_maximum_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = mdata[3][:,1:nt,1:12];

    v = defVar(ds,"meridionalavg_monthly_minimum_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = mdata[4][:,1:nt,1:12];


    v = defVar(ds,"meridionalavg_monthly_mean_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = mdata[1][:,nt+2,1:12];

    v = defVar(ds,"meridionalavg_monthly_std_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = mdata[2][:,nt+2,1:12];

    v = defVar(ds,"meridionalavg_monthly_maximum_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = mdata[3][:,nt+2,1:12];

    v = defVar(ds,"meridionalavg_monthly_minimum_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = mdata[4][:,nt+2,1:12];

    close(ds);

    @info "$(Dates.now()) - Analysed $(info["source"]) $(info["product"]) data for the year $yr in $(gregionfullname(region)) has been saved into file $(fnc) and moved to the data directory $(fol)."

end
