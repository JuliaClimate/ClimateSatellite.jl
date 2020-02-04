using ClimateSatellite

global_logger(ConsoleLogger(stderr,Logging.Info))

user = "natgeo.wong@outlook.com"
clisatdownload("gpmimerg",Date(2007,1),email=user,regions=["SGP"])
data,info,grid = clisatrawall("gpmimerg","prcp_rate",
                             Date(2007,1,3),Date(2007,1,14),
                             region="SGP");
data,info = clisatrawpoint("gpmimerg","prcp_rate",
                               Date(2007,1,3),Date(2007,1,14),
                               coord=[103.4,2.2],region="SGP");
data,info,grid = clisatrawgrid("gpmimerg","prcp_rate",
                                   Date(2007,1,3),Date(2007,1,14),
                                   grid=[2.5,1.9,105,103.4],region="SGP");
