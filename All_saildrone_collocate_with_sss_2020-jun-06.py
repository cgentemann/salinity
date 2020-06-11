import sys
from os import path
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import xarray as xr
import cartopy.crs as ccrs
from pyresample.geometry import AreaDefinition
from pyresample.geometry import GridDefinition
from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
from pyresample.kd_tree import resample_nearest
from scipy import spatial
sys.path.append('../saildrone/subroutines/')
from read_routines import read_all_usv, add_coll_vars,get_filelist_l2p,get_orbital_data_l2p
import warnings
warnings.simplefilter('ignore') # filter some warning messages
from glob import glob

dir_data = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/' #'f:/data/cruise_data/saildrone/saildrone_data/'
dir_data_pattern = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/*.nc' 

data_dict = read_all_usv(dir_data_pattern)
data_dict = add_coll_vars(data_dict)

input_iusv_start = int(input("Enter start cruise processing number 0-44: "))
input_iusv_end = int(input("Enter stop cruise processing number 0-44: "))

adir = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/'
#for name in data_dict:
for iname,name in enumerate(data_dict):
    if iname<input_iusv_start:
        continue
    if iname>input_iusv_end:
        continue
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

    for isat in range(2):

        ds_usv,name_usv=data_dict[name],name

        if isat==0:
            fileout = 'F:/data/cruise_data/saildrone/sss_collocations/'+name_usv+'rssv4_filesave3.nc'
        if isat==1:
            fileout = 'F:/data/cruise_data/saildrone/sss_collocations/'+name_usv+'jplv4.2_filesave3.nc'   

        if path.exists(fileout):
            continue

        #search usv data
        minday,maxday = ds_usv.time[0],ds_usv.time[-1]
        usv_day = minday
        print(iname,name)
        print(minday.data,maxday.data)
        while usv_day<=maxday:
            ds_day = ds_usv.sel(time=slice(usv_day-np.timedelta64(1,'D'),usv_day+np.timedelta64(1,'D')))
            ilen = ds_day.time.size
            if ilen<1:   #don't run on days without any data
                continue
            minlon,maxlon,minlat,maxlat = ds_day.lon.min().data,ds_day.lon.max().data,ds_day.lat.min().data,ds_day.lat.max().data
            filelist = get_filelist_l2p(isat, usv_day)
            x,y,z = [],[],[]
            for ifile,file in enumerate(filelist):
#                if ifile!=7:
#                    continue
#            for file in filelist:
                ds = xr.open_dataset(file)
                ds.close()  
                #print('****************')
                #print(file)
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
                
                # Resample swath to a fixed 0.01 x 0.01 grid, represented by the variable grid_def:
                # https://stackoverflow.com/questions/58065055/floor-and-ceil-with-number-of-decimals
                #changed to be just the region of the usv cruise to make grid even smaller (hopefully)
                #when working with global orbital data, work with usv BUT
                #when working with granules use ds instead of ds_usv so you just do granule region
                grid_def_lon_min = np.round(ds_day.lon.min().data - 0.5 * 10**(-2), 2)
                grid_def_lon_max = np.round(ds_day.lon.max().data + 0.5 * 10**(-2), 2)
                grid_def_lat_min = np.round(ds_day.lat.min().data - 0.5 * 10**(-2), 2)
                grid_def_lat_max = np.round(ds_day.lat.max().data + 0.5 * 10**(-2), 2)
                grid_def_lons = np.arange(grid_def_lon_min,grid_def_lon_max+0.1,0.1)
                grid_def_lats = np.arange(grid_def_lat_max,grid_def_lat_min-0.1,-0.1)
                grid_mesh_lons,grid_mesh_lats = np.meshgrid(grid_def_lons,grid_def_lats)
                # Since we have the lon and lat values for the area, we define a grid instead of an area:
                # https://pyresample.readthedocs.io/en/latest/geo_def.html#griddefinition
                grid_def = GridDefinition(lons=grid_mesh_lons,lats=grid_mesh_lats)

                result1 = resample_nearest(swath_def, data, grid_def, radius_of_influence=20000, fill_value=None)
#                result1 = resample_nearest(swath_def, data, area_def, radius_of_influence=20000, fill_value=None)
                da = xr.DataArray(result1,name='sss',coords={'lat':grid_def_lats,'lon':grid_def_lons},dims=('lat','lon'))
#                da = xr.DataArray(result1,name='sss',coords={'lat':rlat,'lon':rlon},dims=('lat','lon'))

                numdata = np.isfinite(da).sum()
                if numdata<1:
                    continue

                #chekc this!!!! for salinity data, da goes from 90 to -90, so slice max,min
                #print('lat:',maxlat,minlat)
                #print('lon:',maxlon,minlon)
#                subset = da.sel(lat = slice(maxlat,minlat),lon=slice(minlon,maxlon))
#                num_obs = np.isfinite(subset).sum()
#                if num_obs<1:  #no collocations so go to next orbit
#                    continue

                #stack xarray dataset then drop lon == nan
                ds2 = ds.stack(z=('phony_dim_0', 'phony_dim_1')).reset_index('z')
                #drop nan
                ds_drop = ds2.where(np.isfinite(ds2.lon),drop=True)
                lats = ds_drop.lat.data
                lons = ds_drop.lon.data
                inputdata = list(zip(lons.ravel(), lats.ravel()))
                tree = spatial.KDTree(inputdata)
#                ilen = ds_usv.time.size
                #find indices for ds_usv that are within 12 hours of orbit max/min time
                if isat==0:
                    orbit_time = np.datetime64(ds.attrs['time_coverage_start'])-np.timedelta64(24,'h') #changed to 24 hr for sss
                    orbit_time2 = np.datetime64(ds.attrs['time_coverage_end'])+np.timedelta64(24,'h')  
                if isat==1:
                    orbit_time = ds.time[0].data-np.timedelta64(12,'h')
                    orbit_time2 = ds.time[-1].data+np.timedelta64(12,'h')        
                cond = (ds_usv.time.data>orbit_time) & (ds_usv.time.data<orbit_time2)
                item = np.argwhere(cond)
                if item.sum()<1:  #no data within 12 hr of orbit
                    continue
                for iusv_index in range(int(item[0]),int(item[-1])):
#                    if (ds_usv.time[iusv_index]<orbit_time) or (ds_usv.time[iusv_index]>orbit_time2):
#                        continue
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
                        ds_usv.smap_SSS[iusv_index]=ds.smap_sss[ii,jj]
                        ds_usv.smap_iqc_flag[iusv_index]=ds.quality_flag[ii,jj]
                        ds_usv.smap_name[iusv_index]=file
                        ds_usv.smap_dist[iusv_index]=rdist
                        ds_usv.smap_ydim[iusv_index]=ii
                        ds_usv.smap_xdim[iusv_index]=jj
            usv_day += np.timedelta64(1,'D')
        ds_usv.to_netcdf(fileout)