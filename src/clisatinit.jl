"""
This file initializes the ClimateSatellite module by providing the basic
functions that are applicable to all the different satellites, such as the root
directory for which data is stored.

"""

function clisatroot()

    svrdir = "/n/kuangdss01/users/nwong";
    dskdir = "/Volumes/CliNat-Sat")
    docdir = "/Users/natgeo-wong/Documents/Research/"

    if     isdir(svrdir); return svrdir
    elseif isdir(dskdir); return dskdir
    elseif isdir(docdir); return docdir
    else                  error(logger,"The predefined directories in clisatroot.jl do not exist.  They are user-dependent, so please modify/customize accordingly.")
    end

end

function clisatroot(path::AbstractString)
    return path
end
