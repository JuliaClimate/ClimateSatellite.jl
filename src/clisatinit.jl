"""
This file initializes the ClimateSatellite module by providing the basic
functions that are applicable to all the different satellites.  Current
functionalities include:
    - detection of home directory
    - ftp open for arthuhou

"""

# Root Functions
function clisatroot()

    svrdir = "/n/kuangdss01/users/nwong";
    dskdir = "/Volumes/CliNat-Sat")
    docdir = "/Users/natgeo-wong/Documents/Research/"

    if     isdir(svrdir); return svrdir;
        infl("The path $(svrdir) exists and therefore can be used as a directory for ClimateSatellite data downloads.")
    elseif isdir(dskdir); return dskdir
        infl("The path $(dskdir) exists and therefore can be used as a directory for ClimateSatellite data downloads.")
    elseif isdir(docdir); return docdir
        infl("The path $(docdir) exists and therefore can be used as a directory for ClimateSatellite data downloads.")
    else
        errl("The predefined directories in clisatroot.jl do not exist.  They are user-dependent, so please modify/customize accordingly.")
    end

end

function clisatroot(path::AbstractString)
    if isdir(path)
        infl("The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads.")
    end
    return path
end

# FTP Functions
function pmmftpopen()
    email = "natgeo.wong%40outlook.com"
    infl("Opening FTP request to arthurhou.pps.eosdis.nasa.gov.")
    return FTP("ftp://$(email):$(email)@arthurhou.pps.eosdis.nasa.gov")
end

function ppmftpclose(ftp)
    infl("Closing FTP request to arthurhou.pps.eosdis.nasa.gov.")
    close(ftp)
end
