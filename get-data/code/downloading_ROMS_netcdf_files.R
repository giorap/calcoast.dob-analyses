##### -- Download West Coast ROMS data -- #####

#### Load libraries
require(ncdf4)
require(RNetCDF)
require(RCurl) 
library(snow)
cl <- makeCluster(10)

#### Download netcdf files from 2016 to present from THREDDS server
date <- seq(as.Date("2011/1/2"), as.Date("2019/12/31"), by = "day")
# Remove inexistent dates
date <- date[-c(which(date %in% as.Date(c("2011/04/16", "2011/04/03"))))]
var <- rep(c("angle", "h", "ubarrot", "vbarrot", "zeta", "rho", "salt", "temp", "urot", "vrot", "w"), each = length(date))
urls <- paste0("https://thredds.cencoos.org/thredds/ncss/UCSC.nc?var=", 
              "temp", 
              "&disableLLSubset=on&disableProjSubset=on&horizStride=1&time=",
              date,
              "T09%3A00%3A00Z&vertCoord=&addLatLon=true&accept=netcdf"
)
parLapply(cl, urls, function(url){
  download.file(url, destfile = paste0("~/ROMS_download/ROMS_ncdf/West-Coast-ROMS_", 
                                       #strsplit(strsplit(url, "\\?var=")[[1]][2], "\\&")[[1]][1],
                                       #"_",
                                       #substr(strsplit(strsplit(url, "\\?var=")[[1]][2], "\\&")[[1]][5], 6, 15),
                                       gsub("\\https://thredds.cencoos.org/thredds/ncss/UCSC.nc\\?", "", url),
                                       ".nc"), mode="wb")
})

# Alternatively:
# clusterMap(cl, download.file, url = urls, destfile = paste0("~/ROMS_download/ROMS_ncdf/West-Coast-ROMS_", var, "_", date,  ".nc"), mode="wb")

#### Download netcdf files from Aug 2012 to Dec 2015 from west.rsoffice.com:8080 server
date <- seq(as.Date("2015/12/7"), as.Date("2015/12/31"), by = "day")
date2 <- gsub("-", "", date)
url <- paste0("http://west.rssoffice.com:8080/thredds/fileServer/roms/CA3km-nowcast/CA/ca_subCA_das_",
              date2, "09.nc"
)
mapply(download.file, url = url, destfile = paste0("CeNCOOS-ROMS/CeNCOOS-ROMS", "_all-vars", "_", date, ".nc"), mode="wb")
