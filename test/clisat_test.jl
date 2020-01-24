using ClimateSatellite

global_logger(ConsoleLogger(stderr,Logging.Info))

user = "natgeo.wong@outlook.com"
info = clisatdwn(Date(2007,1),productID="pmm_imerg",email=user,regions=["SEA"])
