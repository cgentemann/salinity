#!/usr/bin/env python
# coding: utf-8

# # This is the in situ and SSS collocation code. 
# 

# In[1]:


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


# # Define a function to read in insitu data
# - Read in the Saildrone USV file either from a local disc or using OpenDAP.
# - add room to write collocated data to in situ dataset
# 

# In[2]:


def read_usv(iusv):

    filename_usv_list = ['https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/insitu/L2/spurs2/saildrone/SPURS2_Saildrone1005.nc',
                         'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/insitu/L2/spurs2/saildrone/SPURS2_Saildrone1006.nc',
                         'https://podaac-opendap.jpl.nasa.gov/opendap/allData/insitu/L2/saildrone/Baja/saildrone-gen_4-baja_2018-sd1002-20180411T180000-20180611T055959-1_minutes-v1.nc',
                        'F:/data/cruise_data/access/CTD_casts_ALL_NASA_update_010819.xlsx',
                        'F:/data/cruise_data/saildrone/noaa_arctic/saildrone_PMEL_Arctic_2015_126.nc',
                        'F:/data/cruise_data/saildrone/noaa_arctic/saildrone_PMEL_Arctic_2015_128.nc',
                        'F:/data/cruise_data/saildrone/noaa_arctic/saildrone_PMEL_Arctic_2016_126.nc',
                        'F:/data/cruise_data/saildrone/noaa_arctic/saildrone_PMEL_Arctic_2016_128.nc']
    name_usv_list = ['SPURS2_1005','SPURS2_1006','baja','access',
                     'arctic2015_126',
                     'arctic2015_128',
                     'arctic2016_126',
                     'arctic2016_128']   

    filename_usv = filename_usv_list[iusv]
    if iusv==3:
        df = pd.read_excel(filename_usv, sheet_name='data')
        ds_usv = df.to_xarray()
        ds_usv = ds_usv.where(ds_usv.Depth==-2,drop=True)
        ds_usv = ds_usv.swap_dims({'index':'Date'}).rename({'Date':'time','Longitude':'lon','Latitude':'lat','Salinity':'salinity'}).sortby('time')
    elif iusv<3:
        ds_usv = xr.open_dataset(filename_usv)
        ds_usv.close()
        if iusv==2:
            ds_usv = ds_usv.isel(trajectory=0).swap_dims({'obs':'time'}).rename({'longitude':'lon','latitude':'lat','SAL_MEAN':'salinity'})
            ds_usv = ds_usv.sel(time=slice('2018-04-12T02','2018-06-10T18')) #get rid of last part and first part where USV being towed
        else:
            ds_usv = ds_usv.rename({'longitude':'lon','latitude':'lat','sss':'salinity'})
    elif iusv>3:
        ds_usv = xr.open_dataset(filename_usv)
        ds_usv.close()
        ds_usv = ds_usv.isel(trajectory=0).swap_dims({'obs':'time'}).rename({'longitude':'lon','latitude':'lat','sal_mean':'salinity'})

            #    ds_usv['lon'] = ds_usv.lon.interpolate_na(dim='time',method='linear') #there are 6 nan values
#    ds_usv['lat'] = ds_usv.lat.interpolate_na(dim='time',method='linear')

    #add room to write collocated data information
    ilen = ds_usv.time.shape[0]
    ds_usv['deltaT']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_SSS']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_name']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_dist']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_ydim']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_xdim']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['smap_iqc_flag']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))

    return ds_usv,name_usv_list[iusv]


# ## explore the in situ data and quickly plot using cartopy
#     

# In[ ]:


#intialize grid
for iusv in range(8):
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

    for isat in range(0,2):

        ds_usv,name_usv = read_usv(iusv)

        if isat==0:
            sat_directory = 'F:/data/sat_data/smap/SSS/L2/RSS/V3/40km/'
    #        sat_directory = 'Z:/SalinityDensity/smap/L2/RSS/V3/SCI/40KM/'
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+name_usv+'_rss40km_filesave2.nc'
            file_end = '/*.nc'
        if isat==1:
            sat_directory = 'F:/data/sat_data/smap/SSS/L2/JPL/V4.2/'
    #        sat_directory = 'Z:/SalinityDensity/smap/L2/JPL/V4.2/'
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+name_usv+'_jplv4.2_filesave2.nc'   
            file_end = '/*.h5'

        if path.exists(fileout):
            continue
        #init filelist
        file_save=[]

        #search usv data
        minday,maxday = ds_usv.time[0],ds_usv.time[-1]
        usv_day = minday
        print(minday.data,maxday.data)
        while usv_day<=maxday:
    #        check_day = np.datetime64(str(usv_day.dt.year.data)+'-'+str(usv_day.dt.month.data).zfill(2)+'-'+str(usv_day.dt.day.data).zfill(2))
    #        usv_day1 = usv_day + np.timedelta64(1,'D')
    #        check_day1 = np.datetime64(str(usv_day1.dt.year.data)+'-'+str(usv_day1.dt.month.data).zfill(2)+'-'+str(usv_day1.dt.day.data).zfill(2))
    #        ds_day = ds_usv.sel(time=slice(check_day,check_day1))
            ds_day = ds_usv.sel(time=slice(usv_day-np.timedelta64(1,'D'),usv_day+np.timedelta64(1,'D')))
            ilen = ds_day.time.size
            if ilen<1:   #don't run on days without any data
                continue
            minlon,maxlon,minlat,maxlat = ds_day.lon.min().data,ds_day.lon.max().data,ds_day.lat.min().data,ds_day.lat.max().data
            #caluclate filelist
            filelist = glob(sat_directory+str(usv_day.dt.year.data)+'/'+str(usv_day.dt.dayofyear.data)+file_end)   
            x,y,z = [],[],[]
            for file in filelist:
                file.replace('\\','/')
                ds = xr.open_dataset(file)
                ds.close()
                if isat==0:  #change RSS data to conform with JPL definitions
                    ds = ds.isel(look=0)
                    ds = ds.rename({'cellon':'lon','cellat':'lat','sss_smap':'smap_sss'})
                    ds['lon']=np.mod(ds.lon+180,360)-180  
                x = ds.lon.fillna(-89).data 
                y = ds.lat.fillna(-89).data 
                z = ds.smap_sss.data 
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


# ## Now, loop through only the files that we know have some data in the region of interest.  Use the fast search kdtree which is part of pyresample software, but I think maybe comes originally from sci-kit-learn.
# 
# - read in the in situ data
# - read in a single orbit of satellite data
# - kdtree can't handle it when lat/lon are set to nan.  I frankly have no idea why there is orbital data for both the JPL and RSS products that have nan for the geolocation.  That isn't normal.  But, okay, let's deal with it.  
# - stack the dataset scanline and cell positions into a new variable 'z'
# - drop all variables from the dataset when the longitude is nan
# - set up the tree
# - loop through the orbital data
# - only save a match if it is less than 0.25 deg distance AND time is less than any previous match
# - save the satellite indices & some basic data onto the USV grid
# 

# In[ ]:


for num_usv in range(8):
    for isat in range(2):
        ds_usv,usv_name = read_usv(num_usv)
        if isat==0:
            filelist = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'rss40km_filesave2.nc'
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'rss40km_usv2.nc'
        if isat==1:
            filelist = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'jplv4.2_filesave2.nc'   
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'jplv42_usv2.nc'   
        df = xr.open_dataset(filelist)
        print(isat)
        for file2 in df.filenames.data:
            file = file2
            file.replace('\\','/')
            ds = xr.open_dataset(file)
            ds.close()  
            if isat==0:  #change RSS data to conform with JPL definitions
                ds = ds.isel(look=0)
                ds = ds.rename({'iqc_flag':'quality_flag','cellon':'lon','cellat':'lat','sss_smap':'smap_sss','ydim_grid':'phony_dim_0','xdim_grid':'phony_dim_1'})
                ds['lon']=np.mod(ds.lon+180,360)-180  
            if isat==1:  #change RSS data to conform with JPL definitions
                ds = ds.rename({'row_time':'time'})
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
            usv_subset = ds_usv.sel(time=slice(orbit_time,orbit_time2))
            ilen = ds_usv.time.size
            for iusv in range(ilen):
                if (ds_usv.time[iusv]<orbit_time) or (ds_usv.time[iusv]>orbit_time2):
                    continue
                pts = np.array([ds_usv.lon[iusv], ds_usv.lat[iusv]])
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
                    deltaTa = ((ds_usv.time[iusv]-ds.time[ii,jj]).data)/ np.timedelta64(1,'m')
                if isat==1:
                    deltaTa = ((ds_usv.time[iusv]-ds.time[ii]).data)/ np.timedelta64(1,'m')
                if np.abs(deltaTa)<np.abs(ds_usv.deltaT[iusv].data):
                    ds_usv.deltaT[iusv]=deltaTa
                    ds_usv.smap_SSS[iusv]=ds.smap_sss[ii,jj]
                    ds_usv.smap_iqc_flag[iusv]=ds.quality_flag[ii,jj]
                    ds_usv.smap_name[iusv]=file2
                    ds_usv.smap_dist[iusv]=rdist
                    ds_usv.smap_ydim[iusv]=ii
                    ds_usv.smap_xdim[iusv]=jj
        ds_usv.to_netcdf(fileout)


# In[ ]:


for num_usv in range(8):
    for isat in range(2):
        ds_usv,usv_name = read_usv(num_usv)
        if isat==0:
            file = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'_rss40km_usv2.nc'
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'_rss40km_usv2_norepeats.nc'
        if isat==1:
            file = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'_jplv42_usv2.nc'   
            fileout = 'F:/data/cruise_data/saildrone/sat_collocations/'+usv_name+'_jplv42_usv2_norepeats.nc'   
        ds_usv=xr.open_dataset(file)
        ds_usv.close()
        ds_usv = ds_usv.where(ds_usv.smap_SSS<10000,np.nan)
        ilen,index = ds_usv.dims['time'],0
        ds_tem = ds_usv.copy(deep=True)
        duu, duu2, duv1, duv2, dlat, dlon, dut = [],[],[],[],[],[],np.empty((),dtype='datetime64')
        index=0
        while index <= ilen-2:
            index += 1
            if np.isnan(ds_usv.smap_SSS[index]):
                continue
            if np.isnan(ds_usv.smap_xdim[index]):
                continue
            result = np.where((ds_usv.smap_xdim == ds_tem.smap_xdim[index].data) & (ds_usv.smap_ydim == ds_tem.smap_ydim[index].data))       
            duu=np.append(duu,ds_usv.smap_SSS[result[0][0]].data)
            duu2=np.append(duu2,ds_usv.smap_iqc_flag[result[0][0]].data)
            duv1=np.append(duv1,ds_usv.SAL_MEAN[result].mean().data)
            dlat=np.append(dlat,ds_usv.lat[result].mean().data)
            dlon=np.append(dlon,ds_usv.lon[result].mean().data)
            dut=np.append(dut,ds_usv.time[result].mean().data)
            ds_usv.smap_SSS[result]=np.nan
        dut2 = dut[1:]  #remove first data point which is a repeat from what array defined       
        ds_new=xr.Dataset(data_vars={'smap_SSS': ('time',duu),'smap_iqc_flag': ('time',duu2),
                                     'SAL_MEAN':('time',duv1),
                                     'lon': ('time',dlon),
                                     'lat': ('time',dlat)},
                          coords={'time':dut2})
        ds_new.to_netcdf(fileout)

