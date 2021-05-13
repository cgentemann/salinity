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
from read_routines import read_all_usv,read_one_usv, add_coll_vars_ds_rss, add_coll_vars_ds_jpl,get_filelist_l2p,get_orbital_data_l2p
import warnings
warnings.simplefilter('ignore') # filter some warning messages
from glob import glob


ds = xr.open_dataset('f:/data/sat_data/distance_to_land_25km.nc').rename({'i2':'lon','j2':'lat'})
ds['lat'],ds['lon']=np.arange(-89.875,89.876,.25),np.arange(-179.875,179.876,.25)#np.arange(.125,359.876,.25)
#ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
#ds = ds.sortby(ds.lon)
ds_land=ds

dir_data = 'C:/Users/gentemann/Google Drive/public/ALL_Saildrone_Data/' 
dir_data_pattern = 'C:/Users/gentemann/Google Drive/public/ALL_Saildrone_Data/*.nc' 
#dir_data = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/' #'f:/data/cruise_data/saildrone/saildrone_data/'
#dir_data_pattern = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/*.nc' 

input_iusv_start = int(input("Enter start cruise processing number 0-44: "))
input_iusv_end = int(input("Enter stop cruise processing number 0-44: "))

#adir = 'C:/Users/gentemann/Google Drive/public/2019_saildrone/'
adir = 'C:/Users/gentemann/Google Drive/public/ALL_Saildrone_Data/'

files = glob(dir_data_pattern)
print('number of file:',len(files))
for ifile,file in enumerate(files):
    print(ifile,file)    
    
def get_time_start_end(isat,ds):
    if isat==0:
        orbit_time = np.datetime64(ds.attrs['time_coverage_start'])-np.timedelta64(24,'h') #changed to 24 hr for sss
        orbit_time2 = np.datetime64(ds.attrs['time_coverage_end'])+np.timedelta64(24,'h')  
    if isat==1:
        orbit_time = ds.time[0].data-np.timedelta64(12,'h')
        orbit_time2 = ds.time[-1].data+np.timedelta64(12,'h')        
    return orbit_time,orbit_time2

area_def = load_area('areas.cfg', 'pc_world')
rlon=np.arange(-180,180,.1)
rlat=np.arange(90,-90,-.1)

#for name in data_dict:
for iname in range(input_iusv_start,input_iusv_end): #g,name in enumerate(data_dict):
       
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

        ds_usv,name_usv = read_one_usv(files[iname])

        if isat==0:
            fileout = 'F:/data/cruise_data/saildrone/sss/sss_collocations_orbital/'+name_usv+'rssv04.0_orbital.nc'
            ds_usv = add_coll_vars_ds_rss(ds_usv)
        if isat==1:
            fileout = 'F:/data/cruise_data/saildrone/sss/sss_collocations_orbital/'+name_usv+'jplv05.0_orbital.nc'   
            ds_usv = add_coll_vars_ds_jpl(ds_usv)

        #if path.exists(fileout):
        #    continue

        #search usv data
        minday,maxday = ds_usv.time[0],ds_usv.time[-1]
        usv_day = minday
        print(iname,name_usv)
        print(minday.data,maxday.data)
        while usv_day<=maxday:
            print(usv_day.data,maxday.data)
            ds_day = ds_usv.sel(time=slice(usv_day-np.timedelta64(1,'D'),usv_day+np.timedelta64(1,'D')))
            ilen = ds_day.time.size
            if ilen<1:   #don't run on days without any data
                usv_day += np.timedelta64(1,'D')
                continue
            minlon,maxlon,minlat,maxlat = ds_day.lon.min().data,ds_day.lon.max().data,ds_day.lat.min().data,ds_day.lat.max().data
            filelist = get_filelist_l2p(isat, usv_day)
            x,y,z = [],[],[]
            for ifile,file in enumerate(filelist):
                ds = xr.open_dataset(file)
                ds.close()  
                #print('****************')
                #print(file)
                if isat==0:  #change RSS data to conform with JPL definitions
                    ds = ds.isel(look=0)
                    ds = ds.rename({'iqc_flag':'quality_flag','cellon':'lon','cellat':'lat','sss_smap':'smap_sss','sss_smap_40km':'smap_sss_40km','ydim_grid':'phony_dim_0','xdim_grid':'phony_dim_1'})
                    ds['lon']=np.mod(ds.lon+180,360)-180  
                if isat==1:  #change RSS data to conform with JPL definitions
                    ds = ds.rename({'row_time':'time','ice_concentration':'fice'})

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
                grid_def_lon_min, grid_def_lon_max = np.round(ds_day.lon.min().data - 0.5 * 10**(-2), 2), np.round(ds_day.lon.max().data + 0.5 * 10**(-2), 2)
                grid_def_lat_min, grid_def_lat_max = np.round(ds_day.lat.min().data - 0.5 * 10**(-2), 2), np.round(ds_day.lat.max().data + 0.5 * 10**(-2), 2)
                grid_def_lons, grid_def_lats = np.arange(grid_def_lon_min,grid_def_lon_max+0.1,0.1), np.arange(grid_def_lat_max,grid_def_lat_min-0.1,-0.1)
                grid_mesh_lons,grid_mesh_lats = np.meshgrid(grid_def_lons,grid_def_lats)

                # Since we have the lon and lat values for the area, we define a grid instead of an area:
                # https://pyresample.readthedocs.io/en/latest/geo_def.html#griddefinition
                grid_def = GridDefinition(lons=grid_mesh_lons,lats=grid_mesh_lats)

                result1 = resample_nearest(swath_def, data, grid_def, radius_of_influence=20000, fill_value=None)
                da = xr.DataArray(result1,name='sss',coords={'lat':grid_def_lats,'lon':grid_def_lons},dims=('lat','lon'))

                numdata = np.isfinite(da).sum()
                if numdata<1:
                    continue

                #stack xarray dataset then drop lon == nan
                ds2 = ds.stack(z=('phony_dim_0', 'phony_dim_1')).reset_index('z')
                #drop nan
                ds_drop = ds2.where(np.isfinite(ds2.lon),drop=True)
                lats = ds_drop.lat.data
                lons = ds_drop.lon.data
                inputdata = list(zip(lons.ravel(), lats.ravel()))
                tree = spatial.KDTree(inputdata)

                orbit_time, orbit_time2 = get_time_start_end(isat,ds)

                cond = (ds_usv.time.data>orbit_time) & (ds_usv.time.data<orbit_time2)
                item = np.argwhere(cond)
                if item.sum()<1:  #no data within 12 hr of orbit
                    continue
                for iusv_index in range(int(item[0]),int(item[-1])):
                    pts = np.array([ds_usv.lon[iusv_index], ds_usv.lat[iusv_index]]) #pts = np.array([ds_usv.lon[iusv]+360
                    tree.query(pts,k=1)
                    i = tree.query(pts)[1]
                    rdist = tree.query(pts)[0]                   
                    if rdist>.25:    #don't use matchups more than 25 km away
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
                        if isat==0:
                            ds_usv.smap_SSS_40km[iusv_index]=ds.smap_sss_40km[ii,jj]
                            ds_usv.smap_fland[iusv_index]=ds.fland[ii,jj]
                            ds_usv.smap_rev_number[iusv_index]=int(ds.attrs['orbit_number'])
                        else:
                            ds_usv.smap_rev_number[iusv_index]=int(ds.attrs['REVNO'])
                        ds_usv.smap_iqc_flag[iusv_index]=ds.quality_flag[ii,jj].astype('int') #test int
                        ds_usv.smap_name[iusv_index]=str(file)
                        ds_usv.smap_fice[iusv_index]=ds.fice[ii,jj]
                        ds_usv.smap_dist[iusv_index]=rdist
                        ds_usv.smap_ydim[iusv_index]=ii
                        ds_usv.smap_xdim[iusv_index]=jj
            usv_day += np.timedelta64(1,'D')

        #add distance to land
        ds_usv['dist_land']=ds_land.dist_land.interp(lat=ds_usv.lat,lon=ds_usv.lon).drop({'lat','lon'})
        lnd_att={'long_name':'distance to nearest land','units':'km'}
        ds_usv['dist_land'].attrs=lnd_att    

        ds_usv.deltaT.attrs={'long_name':'time difference between saildrone and satellite','units':'minutes'}
        ds_usv.smap_SSS.attrs=ds.smap_sss.attrs
        if isat==0:
            ds_usv.smap_SSS_40km.attrs=ds.smap_sss_40km.attrs
            ds_usv.smap_fland.attrs=ds.fland.attrs
        ds_usv.smap_iqc_flag.attrs=ds.quality_flag.attrs
        ds_usv.smap_fice.attrs=ds.fice.attrs


        ds_usv.to_netcdf(fileout)
