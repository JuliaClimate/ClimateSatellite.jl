

function clisatcompile(
    productID::AbstractString;
    varname::AbstractString,
    path::AbstractString="", region::AbstractString="GLB",
    trange::Union{Integer,Vector}
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end
    isgeoregion(region);

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);
    clisatvarinfo!(info,productID,varname=varname);
    lon,lat = gpmlonlat(); rlon,rlat,_ = gregiongridvec(regID,lon,lat)
    grid = [rlon,rlat]; nlon = length(rlon); nlat = length(rlat)
    yrmin = minimum(trange); yrmax = maximum(trange); nt = length(yrmin:yrmax)

    @info "$(Dates.now()) - Preallocating arrays ..."
    csavg = zeros(nlon,nlat,nt); csrng = zeros(nlon,nlat,nt);
    csdhr = zeros(nlon,nlat,nt); csitr = zeros(nlon,nlat,nt); cssea = zeros(nlon,nlat,nt);

    @info "$(Dates.now()) - Extracting preliminarily-analyzed $(info["source"]) $(info["product"]) data ..."
    for yr = yrmin : yrmax; it = it + 1;

        csds,csvar = clisatanaread(
            varname,"domain_yearly_mean_climatology",
            Date(yr),region,info
        )
        csavg[:,:,it] = csvar[:]
        close(csds)

        ds1,csmax = clisatanaread(
            varname,"domain_yearly_maximum_climatology",
            Date(yr),region,info
        )
        ds2,csmin = clisatanaread(
            varname,"domain_yearly_minimum_climatology",
            Date(yr),region,info
        )
        csrng[:,:,it] = csmax[:] .- csmin[:]
        close(ds1); close(ds2)

        csds,csvar = clisatanaread(
            varname,"domain_monthly_mean_climatology",
            Date(yr),region,info
        )
        cssea[:,:,it] = maximum(csvar[:],dims=3) .- minimum(csvar[:],dims=3)
        close(csds)

        ds1,csmax = clisatanaread(
            varname,"domain_monthly_maximum_climatology",
            Date(yr),region,info
        )
        ds2,csmin = clisatanaread(
            varname,"domain_monthly_minimum_climatology",
            Date(yr),region,info
        )
        csitr[:,:,it] = mean(csmax[:] .- csmin[:],dims=3)
        close(ds1); close(ds2)

        csds,csvar = clisatanaread(
            varname,"domain_yearly_mean_diurnalvariance",
            Date(yr),region,info
        )
        csdhr[:,:,it] = csvar[:]
        close(csds)

        csrng[:,:,it] = csrng[:,:,it] .- (cssea[:,:,it] .+ csitr[:,:,it])

    end

    @info "$(Dates.now()) - Calculating yearly mean, and diurnal, seasonal and interannual variability ..."
    csian = dropdims(maximum(csavg,dims=3) .- minimum(csavg,dims=3),dims=3)
    csavg = dropdims(mean(csavg,dims=3),dims=3)
    csitr = dropdims(mean(csitr,dims=3),dims=3)
    cssea = dropdims(mean(cssea,dims=3),dims=3)
    csdhr = dropdims(mean(csdhr,dims=3),dims=3)

    clisatcmpsave(csavg,csdhr,csitr,cssea,csian,varname,region,grid,info)

end

function clisatcmpsave(
    csavg::Array{<:Real,2}, csdhr::Array{<:Real,2},
    csitr::Array{<:Real,2}, cssea::Array{<:Real,2}, csian::Array{<:Real,2},
    varname::AbstractString,
    region::AbstractString, grid::Vector{<:Any}, info::AbstractDitc
)

    @info "$(Dates.now()) - Saving compiled $(info["source"]) $(info["product"]) data for the year $yr in $(gregionfullname(region)) ..."

    cfol = clisatanafolder(info,region);
    fcmp = clisatcmpname(info,varname,region);
    cfnc = joinpath(cfol,fcmp);

    if isfile(cfnc)
        @info "$(Dates.now()) - Stale NetCDF file $(cfnc) detected.  Overwriting ..."
        rm(cfnc);
    end

    @debug "$(Dates.now()) - Creating NetCDF file $(cfnc) for compiled $(info["source"]) $(info["product"]) data ..."

    ds = Dataset(cfnc,"c");
    ds.dim["longitude"] = length(grid[1])
    ds.dim["latitude"]  = length(grid[2])

    nclon = defVar(ds,"longitude",Float64,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float64,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nclon[:] = grid[1]; nclat[:] = grid[2]

    ncavg = defVar(ds,"average",Float32,("longitude","latitude"),
        attrib = Dict(
            "long_name" => info["variable"],
            "full_name" => info["standard"],
            "units"     => info["units"],
    ))

    ncian = defVar(ds,"variability_interannual",Float32,("longitude","latitude"),
        attrib = Dict(
            "long_name" => info["variable"],
            "full_name" => info["standard"],
            "units"     => info["units"],
    ))

    ncsea = defVar(ds,"variability_seasonal",Float32,("longitude","latitude"),
        attrib = Dict(
            "long_name" => info["variable"],`
            "full_name" => info["standard"],
            "units"     => info["units"],`
    ))

    ncitr = defVar(ds,"variability_intraseasonal",Float32,("longitude","latitude"),
        attrib = Dict(
            "long_name" => info["variable"],
            "full_name" => info["standard"],
            "units"     => info["units"],
    ))

    ncdhr = defVar(ds,"variability_diurnal",Float32,("longitude","latitude"),
        attrib = Dict(
            "long_name" => info["variable"],
            "full_name" => info["standard"],
            "units"     => info["units"],
    ))

    ncavg[:] = csavg; ncian[:] = csian; ncsea[:] = cssea;
    ncitr[:] = csitr; ncdhr[:] = csdhr

    @info "$(Dates.now()) - Compiled $(info["source"]) $(info["product"]) data for $(gregionfullname(region)) has been saved into file $(cfnc) and moved to the data directory $(cfol)."

end
