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

#effort to combine the finding & collocating code
#intialize grid
adir = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/'
for iusv in range(4,5):
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
                ds = xr.open_dataset(file)
                ds.close()
                if isat==0:  #change RSS data to conform with JPL definitions
                    ds = ds.isel(look=0)
                    ds = ds.rename({'iqc_flag':'quality_flag','cellon':'lon','cellat':'lat','sss_smap':'smap_sss','ydim_grid':'phony_dim_0','xdim_grid':'phony_dim_1'})
                    ds['lon']=np.mod(ds.lon+180,360)-180
                if isat==1:  #change RSS data to conform with JPL definitions
                    ds = ds.rename({'row_time':'time'})

#first do a quick check using resample to project the orbit onto a grid
#and quickly see if there is any data in the cruise area on that day
#if there is, then continue to collocation
                x = ds['lon'].fillna(-89).data
                y = ds['lat'].fillna(-89).data
                z = ds['smap_sss'].data
                lons,lats,data = x,y,z
                swath_def = SwathDefinition(lons, lats)
                result1 = resample_nearest(swath_def, data, area_def, radius_of_influence=20000, fill_value=None)
                da = xr.DataArray(result1,name='sss',coords={'lat':rlat,'lon':rlon},dims=('lat','lon'))
                subset = da.sel(lat = slice(maxlat,minlat),lon=slice(minlon,maxlon))
                num_obs = np.isfinite(subset).sum()
                if num_obs<1:  #no collocations so go to next orbit
                    continue

                #stack xarray dataset then drop lon == nan
                ds2 = ds.stack(z=('phony_dim_0', 'phony_dim_1')).reset_index('z')
                #drop nan
                ds_drop = ds2.where(np.isfinite(ds2.lon),drop=True)
                lats = ds_drop.lat.data
                lons = ds_drop.lon.data
                inputdata = list(zip(lons.ravel(), lats.ravel()))
                tree = spatial.KDTree(inputdata)
                orbit_time = ds.time.max().data-np.timedelta64(1,'D')
                orbit_time2 = ds.time.max().data+np.timedelta64(1,'D')
         #       usv_subset = ds_usv.sel(time=slice(orbit_time,orbit_time2))
                ilen = ds_usv.time.size
                for iusv_index in range(ilen):
                    if (ds_usv.time[iusv_index]<orbit_time) or (ds_usv.time[iusv_index]>orbit_time2):
                        continue
                    pts = np.array([ds_usv.lon[iusv_index], ds_usv.lat[iusv_index]])
            #        pts = np.array([ds_usv.lon[iusv]+360, ds_usv.lat[iusv]])
                    tree.query(pts,k=1)
                    i = tree.query(pts)[1]
                    rdist = tree.query(pts)[0]
                    #don't use matchups more than 25 km away
                    if rdist>.25:
                        continue
                    #use .where to find the original indices of the matched data point
                    #find by matching sss and lat, just randomly chosen variables, you could use any
                    result = np.where((ds.smap_sss == ds_drop.smap_sss[i].data) & (ds.lat == ds_drop.lat[i].data))
                    listOfCoordinates = list(zip(result[0], result[1]))
                    if len(listOfCoordinates)==0:
                        continue
                    ii, jj = listOfCoordinates[0][0],listOfCoordinates[0][1]
                    if isat==0:
                        deltaTa = ((ds_usv.time[iusv_index]-ds.time[ii,jj]).data)/ np.timedelta64(1,'m')
                    if isat==1:
                        deltaTa = ((ds_usv.time[iusv_index]-ds.time[ii]).data)/ np.timedelta64(1,'m')
                    if np.abs(deltaTa)<np.abs(ds_usv.deltaT[iusv_index].data):
                        ds_usv.deltaT[iusv_index]=deltaTa
                        ds_usv.smap_sss[iusv_index]=ds.smap_sss[ii,jj]
                        ds_usv.smap_iqc_flag[iusv_index]=ds.quality_flag[ii,jj]
                        ds_usv.smap_name[iusv_index]=file
                        ds_usv.smap_dist[iusv_index]=rdist
                        ds_usv.smap_cell[iusv_index]=ii #phony_dim0
                        ds_usv.smap_scan[iusv_index]=jj #phony_dim1
                        if isat == 0:
                            ds_usv.smap_SSS_40km[iusv_index] = ds.sss_smap_40km[ii, jj]
                            ds_usv.smap_wind_speed[iusv_index] = ds.winspd[ii, jj]
                            ds_usv.smap_sss_ref[iusv_index] = ds.sss_ref[ii, jj]
                            ds_usv.smap_wind_dir[iusv_index] = ds.windir[ii, jj]
                        if isat == 1:
                            ds_usv.smap_SSS_40km[iusv_index] = -99
                            ds_usv.smap_wind_speed[iusv_index] = ds.smap_spd[ii, jj]
                            ds_usv.smap_sss_ref[iusv_index] = ds.anc_sss[ii, jj]
                            ds_usv.smap_wind_dir[iusv_index] = ds.smap_ambiguity_dir[ii, jj]
            usv_day += np.timedelta64(1,'D')
        ds_usv.to_netcdf(fileout)


