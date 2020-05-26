#!/usr/bin/env python
# coding: utf-8


def read_usv(adir_usv, iusv):
    import xarray as xr
    import numpy as np

    filename_usv_list = ['pmel_2015_sd126-ALL-1_min-v1.nc',
                         'pmel_2015_sd128-ALL-1_min-v1.nc',
                         'pmel_2016_sd126-ALL-1_min-v1.nc',
                         'pmel_2016_sd128-ALL-1_min-v1.nc',
                         'arctic_2019_sd1033-NRT-1_min-v1.nc',
                         'arctic_2019_sd1034-NRT-1_min-v1.nc',
                         'arctic_2019_sd1035-NRT-1_min-v1.nc',
                         'arctic_2019_sd1036-NRT-1_min-v1.nc',
                         'arctic_2019_sd1037-NRT-1_min-v1.nc',
                         'saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1440_minutes-v1.1564857794963.nc',
                        'wcoast_2018_sd1024-ALL-1_min-v1.nc',
                        'wcoast_2018_sd1025-ALL-1_min-v1.nc',
                        'wcoast_2018_sd1026-ALL-1_min-v1.nc',
                        'wcoast_2018_sd1027-ALL-1_min-v1.nc',
                        'wcoast_2018_sd1028-ALL-1_min-v1.nc']
    name_usv_list = ['pmel_2015_sd126', 'pmel_2015_sd128', 'pmel_2016_sd126', 'pmel_2016_sd128',
                     'arctic2019_1033', 'arctic2019_1034', 'arctic2019_1035', 'arctic2019_1036', 'arctic2019_1037',
                     'antarctic2019','wcoast1025','wcoast1026','wcoast1027','wcoast1028','wcoast1029']

    filename_usv = adir_usv + filename_usv_list[iusv]
    print('FILEIN:', filename_usv)
    ds_usv = xr.open_dataset(filename_usv)
    ds_usv.close()
    # NEED TO FIND OUT IF wind_speed is to/from wind_direction ?
    if (iusv == 0 or iusv == 1):  # 1033
        ds_usv = ds_usv.rename(
            {'temp_air_mean': 'TEMP_AIR_MEAN', 'rh_mean': 'RH_MEAN', 'baro_pres_mean': 'BARO_PRES_MEAN',
             'sal_mean': 'SAL_MEAN', 'temp_ctd_mean': 'TEMP_CTD_MEAN', 'temp_o2_mean': 'TEMP_O2_MEAN',
             'chlor_mean': 'CHLOR_MEAN', 'gust_wnd_mean': 'GUST_WND_MEAN', 'temp_ctd_stddev': 'TEMP_CTD_STDDEV'})
        tem_att = ds_usv.wind_speed_mean.attrs
        ds_usv['wind_speed_mean'] = ds_usv.wind_speed_mean * .51444
        ds_usv.wind_speed_mean.attrs = tem_att
        ds_usv.wind_speed_mean.attrs['units'] = 'm s-1'
        uwnd = ds_usv.wind_speed_mean * np.cos(np.deg2rad(ds_usv.wind_direction_mean))
        vwnd = ds_usv.wind_speed_mean * np.sin(np.deg2rad(ds_usv.wind_direction_mean))
        ds_usv['UWND_MEAN'] = uwnd
        ds_usv.UWND_MEAN.attrs = {'standard_name': 'eastward_wind', 'long_name': 'Eastward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
        ds_usv['VWND_MEAN'] = vwnd
        ds_usv.VWND_MEAN.attrs = {'standard_name': 'northward_wind', 'long_name': 'Northward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
        ilen = ds_usv.time.shape[0]
        ds_usv['WWND_MEAN'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time}, dims=('time'))
        ds_usv.WWND_MEAN.attrs = {'standard_name': 'upward_wind_velocity', 'long_name': 'upward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
    if (iusv == 2 or iusv == 3):  # 1033
        ds_usv = ds_usv.rename(
            {'temp_air_mean': 'TEMP_AIR_MEAN', 'rh_mean': 'RH_MEAN', 'baro_pres_mean': 'BARO_PRES_MEAN',
             'sal_mean': 'SAL_MEAN', 'temp_ctd_mean': 'TEMP_CTD_MEAN', 'temp_o2_mean': 'TEMP_O2_MEAN',
             'chlor_mean': 'CHLOR_MEAN', 'gust_wnd_mean': 'GUST_WND_MEAN', 'temp_ctd_stddev': 'TEMP_CTD_STDDEV'})
        tem_att = ds_usv.wind_speed.attrs
        ds_usv['wind_speed'] = ds_usv.wind_speed * .51444
        ds_usv.wind_speed.attrs = tem_att
        ds_usv.wind_speed.attrs['units'] = 'm s-1'
        uwnd = ds_usv.wind_speed * np.cos(np.deg2rad(ds_usv.wind_direction))
        vwnd = ds_usv.wind_speed * np.sin(np.deg2rad(ds_usv.wind_direction))
        ds_usv['UWND_MEAN'] = uwnd
        ds_usv.UWND_MEAN.attrs = {'standard_name': 'eastward_wind', 'long_name': 'Eastward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
        ds_usv['VWND_MEAN'] = vwnd
        ds_usv.VWND_MEAN.attrs = {'standard_name': 'northward_wind', 'long_name': 'Northward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
        ilen = ds_usv.time.shape[0]
        ds_usv['WWND_MEAN'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time}, dims=('time'))
        ds_usv.WWND_MEAN.attrs = {'standard_name': 'upward_wind_velocity', 'long_name': 'upward wind speed',
                                  'units': 'm s-1', 'installed_height': '5.2'}
    if iusv == 4:  # 1033
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN': 'TEMP_CTD_MEAN', 'TEMP_CTD_RBR_STDDEV': 'TEMP_CTD_STDDEV',
                                'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN', 'SAL_RBR_MEAN': 'SAL_MEAN',
                                'CHLOR_WETLABS_MEAN': 'CHLOR_MEAN'})
    if iusv == 5:  # 1034
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN': 'TEMP_CTD_MEAN', 'TEMP_CTD_RBR_STDDEV': 'TEMP_CTD_STDDEV',
                                'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN', 'SAL_RBR_MEAN': 'SAL_MEAN',
                                'CHLOR_WETLABS_MEAN': 'CHLOR_MEAN'})
    if iusv == 6:  # 1035
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN': 'TEMP_CTD_MEAN', 'TEMP_CTD_RBR_STDDEV': 'TEMP_CTD_STDDEV',
                                'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN', 'SAL_RBR_MEAN': 'SAL_MEAN',
                                'CHLOR_WETLABS_MEAN': 'CHLOR_MEAN'}) #, 'WIND_MEASUREMENT_MEAN_HEIGHT': 'WIND_MEAN_HEIGHT'})
    if iusv == 7:  # 1036
        ds_usv = ds_usv.isel(time=slice(100,
                                        -1))  # ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN': 'TEMP_CTD_MEAN', 'TEMP_CTD_RBR_STDDEV': 'TEMP_CTD_STDDEV',
                                'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN', 'SAL_RBR_MEAN': 'SAL_MEAN',
                                'CHLOR_WETLABS_MEAN': 'CHLOR_MEAN'})
    if iusv == 8:  # 1037
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN': 'TEMP_CTD_MEAN', 'TEMP_CTD_RBR_STDDEV': 'TEMP_CTD_STDDEV',
                                'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN'})
    if iusv == 9:  # 1037
        ds_usv = ds_usv.isel(trajectory=0).swap_dims({'obs': 'time'}).rename(
            {'latitude': 'lat', 'longitude': 'lon', 'TEMP_O2_RBR_MEAN': 'TEMP_O2_MEAN'})  # TEMP_CTD_RBR_MEAN':'TEMP_
    if (iusv == 9 or iusv <= 3):
        ilen = ds_usv.time.shape[0]
        ds_usv['WIND_HEIGHT_MEAN'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time}, dims=('time'))
        ds_usv.WIND_HEIGHT_MEAN.attrs = {'long_name': 'Wind measurement height', 'units': 'm',
                                         'installed_height': '5.2'}
        ds_usv['WAVE_DOMINANT_PERIOD'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time},
                                                      dims=('time'))
        ds_usv.WAVE_DOMINANT_PERIOD.attrs = {
            'standard_name': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
            'long_name': 'Dominant wave period', 'units': 's', 'installed_height': '0.34'}
        ds_usv['WAVE_SIGNIFICANT_HEIGHT'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time},
                                                         dims=('time'))
        ds_usv.WAVE_SIGNIFICANT_HEIGHT.attrs = {'standard_name': 'sea_surface_wave_significant_height',
                                                'long_name': 'Significant wave height', 'units': 'm',
                                                'installed_height': '0.34'}

    # add room to write collocated data information
    ilen = ds_usv.time.shape[0]
   
    ds_usv['deltaT'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_name'] = xr.DataArray(np.empty(ilen, dtype=str), coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_dist'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_scan'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_cell'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_iqc_flag'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_sss'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_SSS_40km'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_wind_speed'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_sss_ref'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))
    ds_usv['smap_wind_dir'] = xr.DataArray(np.ones(ilen) * 999999, coords={'time': ds_usv.time}, dims=('time'))

    return ds_usv, name_usv_list[iusv]


def get_filelist_l2p(isat, day_in):
    """
    :type asat_in: object string  'smap_jpl' or 'smap_rss'
    :type day_in: object datetime
    # the more recent data is in daily directories, so easy to search
    # the older data, pre 2018 is in monthly directories so only search files for day
    """
    import datetime as dt
    from glob import glob   
    
    if isat == 0:
        sat_directory = 'F:/data/sat_data/smap/SSS/L2/RSS/V4/SCI/'
        file_end = '*.nc'
    if isat == 1:  
        sat_directory = 'F:/data/sat_data/smap/SSS/L2/JPL/V4.2/'
        file_end = '*.h5'
        
    syr = str(day_in.dt.year.data)
    smon = str(day_in.dt.month.data).zfill(2)
    sdy = str(day_in.dt.day.data).zfill(2)
    sjdy = str(day_in.dt.dayofyear.data).zfill(3)

    adir_list = sat_directory + syr + '/' + sjdy + '/' + '/*' + file_end
    filelist = glob(adir_list)

    return filelist

def get_orbital_data_l2p(isat,file):
    import xarray as xr
    import numpy as np

    file.replace('\\', '/')
    ds = xr.open_dataset(file)
    ds.close()
    if isat==0:  #change RSS data to conform with JPL definitions
        ds = ds.isel(look=0)
        ds = ds.rename({'iqc_flag':'quality_flag','cellon':'lon','cellat':'lat','sss_smap':'smap_sss','ydim_grid':'phony_dim_0','xdim_grid':'phony_dim_1'})
        ds['lon']=np.mod(ds.lon+180,360)-180
    if isat==1:  #change JPL data to conform with RSS definitions
        ds = ds.rename({'row_time':'time'})
    xlat = ds['lat']
    xlon = ds['lon']
    var_data = ds['smap_sss']
    sat_time = ds['time']
    sat_qc = ds['quality_flag']
    return xlat,xlon,sat_time,var_data,sat_qc
