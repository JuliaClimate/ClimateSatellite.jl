using ClimateSatellite

global_logger(ConsoleLogger(stderr,Logging.Info))

user = "natgeo.wong@outlook.com"
clisatdwn("gpmimerg",Date(2007,1),email=user,regions=["SGP"])
data,info,grid = clisatextractall("gpmimerg","prcp_rate",
                             Date(2007,1,3),Date(2007,1,14),
                             region="SGP");
data,info = clisatextractpoint("gpmimerg","prcp_rate",
                               Date(2007,1,3),Date(2007,1,14),
                               coord=[103.4,2.2],region="SGP");
data,info,grid = clisatextractgrid("gpmimerg","prcp_rate",
                                   Date(2007,1,3),Date(2007,1,14),
                                   grid=[2.5,1.9,105,103.4],region="SGP");
