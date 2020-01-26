"""
This file stores all the WIP general functions for ClimateSatellite.jl

"""

# Root Functions
function clisatextract(;
    productID::AbstractString, variable::AbstractString,
    start::TimeType, finish::TimeType,
    dataroot::AbstractString="",
    region::AbstractString="GLB",
    coord::Array{<:Real,1}
)

    if dataroot == ""; dataroot = clisatroot(productID); end

    info = Dict{Any,Any}("root"=>dataroot,"email"=>replace(email,"@"=>"%40"));
    clisatinfo!(info,productID);

    if !isdir(clisatfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    if sum(info["varID"] .== variable) == 0
        error("$(Dates.now()) - There is no variable identifier $(variable) for $(info["source"]) $(info["product"])")
    end

    dates = convert(Array,start:Dates.month(1):finish)
    ndays = Dates.value((finish-start)/Dates.day(1));
    datavec = Array{<:Real,2}(undef,ndays*nt);

    jj = 0;
    for dateii in dates; jj =+ 1

        fol = clisatfol(info,dateii,region);
        if !isdir(fol)
            error(("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatncname(info,dateii,region));

        ncread!(fnc,variable,dataii)
        datavec[:,:,]


    end

end
