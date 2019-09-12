import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
import xarray as xr
import scipy
from glob import glob
import cartopy.crs as ccrs
from pyresample.geometry import AreaDefinition
from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
from pyresample.kd_tree import resample_nearest
from math import radians, cos, sin, asin, sqrt
from scipy import spatial
import os.path
from os import path
import warnings
warnings.simplefilter('ignore') # filter some warning messages

import sys
sys.path.append('./subroutines/')
from read_routines import read_usv, get_filelist_l2p,get_orbital_data_l2p

input_iusv_start = int(input("Enter start cruise processing number 0-10: "))
input_iusv_end = int(input("Enter stop cruise processing number 0-10: "))

#intialize grid
adir = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/'
for iusv in range(input_iusv_start,input_iusv_end):
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

    for isat in range(0,1):

        ds_usv,name_usv=read_usv(adir,iusv)

        if isat==0:
            fileout = 'F:/data/cruise_data/saildrone/sss_collocations/'+name_usv+'rssv4_filesave2.nc'
        if isat==1:
            fileout = 'F:/data/cruise_data/saildrone/sss_collocations/'+name_usv+'jplv4.2_filesave2.nc'

        if path.exists(fileout):
            continue

        #search usv data
        file_save=[]
        minday,maxday = ds_usv.time[0],ds_usv.time[-1]
        usv_day = minday
        print(minday.data,maxday.data)
        while usv_day<=maxday:
            ds_day = ds_usv.sel(time=slice(usv_day-np.timedelta64(1,'D'),usv_day+np.timedelta64(1,'D')))
            ilen = ds_day.time.size
            if ilen<1:   #don't run on days without any data
                continue
            minlon,maxlon,minlat,maxlat = ds_day.lon.min().data,ds_day.lon.max().data,ds_day.lat.min().data,ds_day.lat.max().data
            filelist = get_filelist_l2p(isat, usv_day)
            x,y,z = [],[],[]
            for file in filelist:
                xlat,xlon,sat_time,var_data,sat_qc = get_orbital_data_l2p(isat,file)
                x = xlon.fillna(-89).data
                y = xlat.fillna(-89).data
                z = var_data.data
                lons,lats,data = x,y,z
                swath_def = SwathDefinition(lons, lats)
                result1 = resample_nearest(swath_def, data, area_def, radius_of_influence=20000, fill_value=None)
                da = xr.DataArray(result1,name='sss',coords={'lat':rlat,'lon':rlon},dims=('lat','lon'))
                subset = da.sel(lat = slice(maxlat,minlat),lon=slice(minlon,maxlon))
                num_obs = np.isfinite(subset).sum()
                if num_obs>0:
                    file_save = np.append(file_save,file)
            usv_day += np.timedelta64(1,'D')
        df = xr.DataArray(file_save,name='filenames')
        df.to_netcdf(fileout)
