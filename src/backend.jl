"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * detection of home directory on machine
    * opening an FTP request to arthuhou and jsimpson (GPM / TRMM data) servers
    * creation of folders (both data folders and temporary directories)

"""

# FTP Functions
function pmmftpopen(server::AbstractString,email::AbstractString)
    @info "$(Dates.now()) - Opening FTP request to $(server).pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@$(server).pps.eosdis.nasa.gov")
end

function pmmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request."
    close(ftp)
end

# Other Functions
function isprod(info::Dict,product::AbstractString)
    occursin(product,info["short"])
end

function clisatextractdate(startdate::TimeType,finish::TimeType);

    yrs = Dates.year(start);  mos = Dates.month(start);  dys = Dates.day(start);
    yrf = Dates.year(finish); mof = Dates.month(finish); dyf = Dates.day(finish);
    ndy = Dates.value((finish-start)/Dates.day(1));
    dvecs = Date(yrs,mos); dvecf = Date(yrf,mof);

    dvec = convert(Array,dvecs:Month(1):dvecf);

    return dvec,dys,dyf,ndy

end

function clisatrmtmp(flist::Array{<:AbstractString,2},fol::AbstractString)

    for ii = 1 : length(flist)
        fii = joinpath(fol,"$(flist[ii])"); if isfile(fii); rm(fii) end
    end

end
