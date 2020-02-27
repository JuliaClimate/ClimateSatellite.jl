"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * detection of home directory on machine
    * opening an FTP request to arthuhou and jsimpson (GPM / TRMM data) servers
    * creation of folders (both data folders and temporary directories)

"""

## FTP Functions
function pmmftpopen(server::AbstractString,email::AbstractString)
    @info "$(Dates.now()) - Opening FTP request to $(server).pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@$(server).pps.eosdis.nasa.gov")
end

function pmmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request."
    close(ftp)
end

## Other Functions
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

function clisattmpfol(info::Dict)

    fol = joinpath(info["root"],"tmp");

    if !isdir(fol)
        @debug "$(Dates.now()) - Creating temporary directory $(fol)."; mkpath(fol);
    end

    return fol

end

## DateString Aliasing
function yrmo2dir(date::TimeType) = Dates.format(date,dateformat"yyyy/mm") end
function yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm") end
function yr2str(date::TimeType)   = Dates.format(date,dateformat"yyyy") end
function ymd2str(date::TimeType)  = Dates.format(date,dateformat"yyyymmdd") end

## Real2Integer Conversion for Packing
function real2int16!(
    outarray::Array{Int16}, inarray::Array{<:Real};
    offset::Real, scale::Real
)

    if size(outarray) != size(inarray)
        dout = [i for i in size(outarray)];
        din  = [i for i in size(inarray)];
        if (dout[1:end-1] != din[1:end] && dout[1:end] != din[1:end-1]) ||
            prod(dout) != prod(din)
            error("$(Dates.now()) - output array is not of the same size as the input array")
        end
    end

    for ii = 1 : length(inarray)

        inarray[ii] = (inarray[ii] - offset) / scale

        if inarray[ii] < -32767; inarray[ii] = -32768
        elseif inarray[ii] > 32767; inarray[ii] = -32768
        end

        outarray[ii] = round(Int16,inarray[ii])

    end

    return

end
