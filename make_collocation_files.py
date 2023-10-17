import datetime
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

from rss_sat_readers import read_rss_satellite_daily_xr,get_sat_def
from rss_plotting.global_map import plot_global_map
from rss_plotting.plot_2d_array import plot_2d_array
from get_era5_2d import ERA5_Hourly

def fit_plane_to_7x7_map(map,return_fit=False):
    # fit a plane to the 7x7 map
    # map is a 7x7 array
    # center_lat and center_lon are the lat and lon of the center of the map
    # the plane is fit to the 49 points in the map
    # the plane is fit in a local coordinate system where the x axis is along the
    # longitude line and the y axis is along the latitude line
    # the plane is fit to the equation z = a + b*x + c*y
    # where z is the map value, x is the longitude offset from the center longitude
    # and y is the latitude offset from the center latitude
    # the function returns a,b,c
    # the function also returns the map values predicted by the plane fit
    # the plane fit is returned as a 7x7 array
    #
    # the plane fit is done using least squares
    # the least squares solution is
    # [a,b,c] = inv(A'*A)*A'*z
    # where A is a 49x3 matrix
    # A(:,1) = 1
    # A(:,2) = x
    # A(:,3) = y
    # z is a 49x1 vector of the map values
    # the plane fit is returned as a 7x7 array
    # the predicted values are returned as a 7x7 array
    # the residuals are returned as a 7x7 array

    # make the A matrix
    A = np.zeros((49,3),dtype=np.float32)
    A[:,0] = 1.0
    A[:,1] = 0.25*np.tile(np.arange(-3,4),(7,1)).T.flatten()
    A[:,2] = 0.25*np.tile(np.arange(-3,4),(7,1)).flatten()

    # make the z vector
    z = map.flatten()

    ok = np.isfinite(z)
    A_ok = A[ok,:]
    z_ok = z[ok]

    # solve for a,b,c
    ATA = np.dot(A_ok.T,A_ok)
    ATz = np.dot(A_ok.T,z_ok)
    abc = np.linalg.solve(ATA,ATz)


    # make the plane fit
    if return_fit:
        plane_fit = abc[0] + abc[1]*A[:,1] + abc[2]*A[:,2]
        plane_fit = plane_fit.reshape((7,7))
        return abc,plane_fit
    
    else:
        return abc




max_collocs = 200000
subset = 'ocean'
sat_name = 'AMSR2'

frac_num_threshold = 10.0
#subset = 'ocean_coastal'

#this is a dictionary needed by the ERA5_Hourly class
var_dict = {
    "skt": "Skin temperature",
    "tcwv": "Total column water vapour",
    "tclw": "total_column_cloud_liquid_water",
    "u10n": "10m_u_component_of_neutral_wind",
    "v10n": "10m_v_component_of_neutral_wind",
}


#read olivier's list of gnss locations
gnss_root = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/gnss_data')
gnss_loc_file = gnss_root / f'gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.{subset}.txt'

gnss_locs = pd.read_csv(gnss_loc_file, delim_whitespace=True)

for year_to_do in range(2013,2023):
    #read in the GNSS data for this year
    gnss_data_file = gnss_root / f'NGL_vapor_{year_to_do}_{subset}.csv'
    gnss_data = pd.read_csv(gnss_data_file)

    # initalize the output arrays
    gnss_lats = np.full((max_collocs),np.nan,dtype=np.float32)
    gnss_lons = np.full((max_collocs),np.nan,dtype=np.float32)

    center_lats = np.full((max_collocs),np.nan,dtype=np.float32)
    center_lons = np.full((max_collocs),np.nan,dtype=np.float32)
    station_ids = np.zeros((max_collocs),dtype='<U4')
    days_1994 = np.full((max_collocs),np.nan,dtype=np.int32)
    date_strings = np.zeros((max_collocs),dtype='<U10')

    gnss_vapor = np.full(max_collocs,np.nan,dtype=np.float32)
    gnss_sigma_vapor = np.full(max_collocs,np.nan,dtype=np.float32)

    sat_vapor_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    sat_wspd_mf_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    sat_cloud_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    sat_rain_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    sat_sst_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    sat_time_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)

    sat_vapor_at_gnss_loc = np.full(max_collocs,np.nan,dtype=np.float32)

    era5_vapor_maps = np.full((max_collocs,7,7),np.nan,dtype=np.float32)
    era5_vapor_at_gnss_loc = np.full(max_collocs,np.nan,dtype=np.float32)
    era5_vapor_at_gnss_loc_sampled = np.full(max_collocs,np.nan,dtype=np.float32)

    sat_fit_coeffs = np.full((max_collocs,3),np.nan,dtype=np.float32)
    era5_fit_coeffs = np.full((max_collocs,3),np.nan,dtype=np.float32)
    era5_fit_coeffs_sampled = np.full((max_collocs,3),np.nan,dtype=np.float32)
    num_found = 0
    sst_found = False

    start_date = datetime.date(year_to_do, 1, 1)
    end_date = datetime.date(year_to_do, 12, 31)

    date_to_do = start_date
    while date_to_do <= end_date:
        var = 'tcwv'
        try:
            variable = (var, var_dict[var])
        except KeyError:
            raise ValueError(f'Variable {var} not found in dictionary')
            
        era5_interp = ERA5_Hourly(variable,date_to_do,Path('B:\_access_temp'),rss_grid=True)
        # read in satellite data
        try:
            d = read_rss_satellite_daily_xr(year=date_to_do.year,
                                            month=date_to_do.month,
                                            day=date_to_do.day,
                                            satellite=sat_name,
                                            verbose=True)
        except FileNotFoundError:
            print('No data for ',date_to_do)
            date_to_do += datetime.timedelta(days=1)
            continue

        # find the gnss data for this day
        num_days_1994 = (date_to_do - datetime.date(1994,1,1)).days

        gnss_data_day = gnss_data.loc[gnss_data['num_days_since_1994'] == num_days_1994]


        for index,row in gnss_data_day.iterrows():
            gnss_id = row['gnss_id']
            location = gnss_locs.loc[gnss_locs['name'] == gnss_id]
            lat = location['lat'].values[0]
            lon = location['lon'].values[0]
            ilat = int(np.floor((4.0*(lat+90.0))))
            ilon = int(np.floor((4.0*(lon))))
            
            if (ilat > 4) or (ilat < 716):
                #print(f'lat={row["lat"]},lon={row["lon"]},ilat={ilat},ilon={ilon}')
                #print(fraction_valid[ilat,ilon,:])
                ilats = np.arange(ilat-3,ilat+4)
                ilons = np.arange(ilon-3,ilon+4)
                ilons[ilons<0] += 1440
                ilons[ilons>=1440] -= 1440
                ilat_grid = np.tile(ilats,(7,1)).T
                ilon_grid = np.tile(ilons,(7,1))

                for node in [0,1]:
                    sat_vap = d['vapor'][node,:,:].values
                    sat_vap_near_loc = sat_vap[ilat_grid,ilon_grid]
                    ok = np.isfinite(sat_vap_near_loc)
                    if np.sum(ok) >= frac_num_threshold:
                        gnss_lats[num_found] = lat
                        gnss_lons[num_found] = lon
                        center_lats[num_found] = -90.0 + 0.125 + 0.25*ilat
                        center_lons[num_found] = 0.125 + 0.25*ilon
                        days_1994[num_found] = num_days_1994
                        date_strings[num_found] = date_to_do.strftime('%Y-%m-%d')
                        station_ids[num_found] = gnss_id
                        gnss_vapor[num_found] = row['iwv']
                        gnss_sigma_vapor[num_found] = row['sigma']
                        sat_vapor_maps[num_found,:,:] = sat_vap_near_loc
                        sat_wspd_mf_maps[num_found,:,:] = d['wspd_mf'][node,:,:].values[ilat_grid,ilon_grid]
                        sat_cloud_maps[num_found,:,:] = d['cloud'][node,:,:].values[ilat_grid,ilon_grid]
                        sat_rain_maps[num_found,:,:] = d['rain'][node,:,:].values[ilat_grid,ilon_grid]
                        try:
                            sat_sst_maps[num_found,:,:] = d['sst'][node,:,:].values[ilat_grid,ilon_grid]
                            sst_found = True
                        except KeyError:
                            pass
                        time_maps_temp = d['time'][node,:,:].values[ilat_grid,ilon_grid].astype(np.float64)
                        time_maps_temp[time_maps_temp < 0.0001] = np.nan
                        sat_time_maps[num_found,:,:] = time_maps_temp
                        mean_time = np.nanmean(time_maps_temp)
                        std_time = np.nanstd(time_maps_temp)
                        if std_time > 5.0:
                            print(f'WARNING: std_time={std_time} for {gnss_id} on {date_to_do}')
                            print(f'time_maps_temp={time_maps_temp}')
                            num_found += 1
                            continue
                        
                        time_map_to_pass = np.ones((7,7))*mean_time

                        #interpolate the era5 data to the satellite time
                        era5_map = era5_interp.interp_7x7(ilat,ilon,time_map_to_pass)
                        era5_vapor_maps[num_found,:,:] = era5_map
                        sat_fit = fit_plane_to_7x7_map(sat_vap_near_loc)
                        sat_fit_coeffs[num_found,:] = sat_fit

                        era5_fit = fit_plane_to_7x7_map(era5_map)
                        era5_fit_coeffs[num_found,:] = era5_fit

                        era5_map_sampled = np.copy(era5_map)
                        era5_map_sampled[~np.isfinite(sat_vap_near_loc)] = np.nan
                        era5_fit_sampled = fit_plane_to_7x7_map(era5_map_sampled)
                        era5_fit_coeffs_sampled[num_found,:] = era5_fit_sampled

                        d_lat = lat-center_lats[num_found]
                        d_lon = lon-center_lons[num_found]

                        #find vapor values at the gnss_location from the 7x7 planar fits
                        sat_vap_from_fit =  sat_fit[0] +  sat_fit[1]*d_lat +  sat_fit[2]*d_lon
                        era5_vap_from_fit = era5_fit[0] + era5_fit[1]*d_lat + era5_fit[2]*d_lon
                        era5_vap_from_fit_sampled = era5_fit_sampled[0] + era5_fit_sampled[1]*d_lat + era5_fit_sampled[2]*d_lon

                        print(f'date = {date_to_do}, lat={lat},lon={lon},center_lat={center_lats[num_found]},center_lon={center_lons[num_found]}')
                        print(sat_vap_from_fit,era5_vap_from_fit,era5_vap_from_fit_sampled,gnss_vapor[num_found])
                        xvals=np.arange(-3,4)
                        yvals=np.arange(-3,4)

                        sat_vapor_at_gnss_loc[num_found] = sat_vap_from_fit
                        era5_vapor_at_gnss_loc[num_found] = era5_vap_from_fit
                        era5_vapor_at_gnss_loc_sampled[num_found] = era5_vap_from_fit_sampled

                        num_found += 1


        date_to_do += datetime.timedelta(days=1)

    gnss_lats = gnss_lats[0:num_found]
    gnss_lons = gnss_lons[0:num_found]
    center_lats = center_lats[0:num_found]
    center_lons = center_lons[0:num_found]
    days_1994 = days_1994[0:num_found]
    date_strings = date_strings[0:num_found]
    station_ids = station_ids[0:num_found]
    gnss_vapor = gnss_vapor[0:num_found]
    gnss_sigma_vapor = gnss_sigma_vapor[0:num_found]
    sat_vapor_maps = sat_vapor_maps[0:num_found,:,:]
    sat_wspd_mf_maps = sat_wspd_mf_maps[0:num_found,:,:]
    sat_cloud_maps = sat_cloud_maps[0:num_found,:,:]
    sat_rain_maps = sat_rain_maps[0:num_found,:,:]
    sat_sst_maps = sat_sst_maps[0:num_found,:,:]
    sat_time_maps = sat_time_maps[0:num_found,:,:]

    sat_vapor_at_gnss_loc = sat_vapor_at_gnss_loc[0:num_found]
    era5_vapor_maps = era5_vapor_maps[0:num_found,:,:]
    era5_vapor_at_gnss_loc = era5_vapor_at_gnss_loc[0:num_found]
    era5_vapor_at_gnss_loc_sampled = era5_vapor_at_gnss_loc_sampled[0:num_found]

    sat_fit_coeffs = sat_fit_coeffs[0:num_found,:]
    era5_fit_coeffs = era5_fit_coeffs[0:num_found,:]
    era5_fit_coeffs_sampled = era5_fit_coeffs_sampled[0:num_found,:]


    print(f'num_found={num_found}')
    print

    ds = xr.Dataset({   
        'gnss_lats': (['colloc_index'],gnss_lats,
                    {'long_name':'latitude of GNSS station',
                    'units':'degrees_north'}),
        'gnss_lons': (['colloc_index'],gnss_lons,
                    {'long_name':'longitude of GNSS station',
                    'units':'degrees_east'}),
        'center_lats': (['colloc_index'],center_lats,
                        {'long_name':'latitude of center of 7x7 map',
                        'units':'degrees_north'}),
        'center_lons': (['colloc_index'],center_lons,
                        {'long_name':'longitude of center of 7x7 map',
                        'units':'degrees_east'}),
        'days_1994': (['colloc_index'],days_1994,
                    {'long_name':'number of days since 1994-01-01'}),
        'date_strings': (['colloc_index'],date_strings,
                        {'long_name':'date as string YYYY-MM-DD'}),
        'station_ids': (['colloc_index'],station_ids,
                        {'long_name':'station id'}),
        'gnss_vapor': (['colloc_index'],gnss_vapor,
                    {'long_name':'GNSS vapor',
                        'units':'mm'}),
        'gnss_sigma_vapor': (['colloc_index'],gnss_sigma_vapor,
                            {'long_name':'GNSS sigma vapor',
                            'units':'mm'}),
        'sat_vapor_maps': (['colloc_index','dlat','dlon'],sat_vapor_maps,
                        {'long_name':'satellite vapor',
                            'units':'mm'}),
        'sat_wspd_mf_maps': (['colloc_index','dlat','dlon'],sat_wspd_mf_maps,
                            {'long_name':'satellite wind speed',
                            'units':'m/s'}),
        'sat_cloud_maps': (['colloc_index','dlat','dlon'],sat_cloud_maps,
                        {'long_name':'satellite cloud fraction',
                            'units':'mm'}),
        'sat_rain_maps': (['colloc_index','dlat','dlon'],sat_rain_maps,
                        {'long_name':'satellite rain rate',
                        'units':'mm/hr'}),
        'sat_sst_maps': (['colloc_index','dlat','dlon'],sat_sst_maps,
                        {'long_name':'satellite sea surface temperature',
                        'units':'deg C'}),
        'sat_time_maps': (['colloc_index','dlat','dlon'],sat_time_maps,
                        {'long_name':'satellite time',
                        'units':'minutes in current day'}),
        'sat_vapor_at_gnss_loc': (['colloc_index'],sat_vapor_at_gnss_loc,
                                {'long_name':'satellite vapor at GNSS location from 7x7 fit'}),
        'sat_fit_coeffs': (['colloc_index','fit_coeffs'],sat_fit_coeffs,
                        {'long_name':'a,b,c coefficients of plane fit to 7x7 map'}),
        'era5_vapor_maps': (['colloc_index','dlat','dlon'],era5_vapor_maps,
                            {'long_name':'era5 vapor',
                            'units':'mm'}),
        'era5_vapor_at_gnss_loc': (['colloc_index'],era5_vapor_at_gnss_loc,
                                {'long_name':'era5 vapor at GNSS location from 7x7 fit'}),
        'era5_fit_coeffs': (['colloc_index','fit_coeffs'],era5_fit_coeffs,
                            {'long_name':'a,b,c coefficients of plane fit to 7x7 map'}),
        'era5_vapor_at_gnss_loc_sampled': (['colloc_index'],era5_vapor_at_gnss_loc_sampled,
                                        {'long_name':'era5 vapor at GNSS location from 7x7 fit, sampled at satellite locations'}),
        'era5_fit_coeffs_sampled': (['colloc_index','fit_coeffs'],era5_fit_coeffs_sampled,
                                    {'long_name':'a,b,c coefficients of plane fit to 7x7 map, sampled at satellite locations'})},
        coords={'colloc_index': (['colloc_index'], np.arange(num_found),
                                {'long_name':'index of collocation'}),
                'dlat': (['dlat'], 0.25*np.arange(-3,4),
                        {'long_name':'offset from center lat in degrees'}),
                'dlon': (['dlon'], 0.25*np.arange(-3,4),
                        {'long_name':'offset from center lon in degrees'}),
                'fit_coeffs': (['fit_coeffs'],np.arange(3),
                            {'long_name':'a,b,c coefficients of plane fit to 7x7 map',
                                'comment':'vapor = a + b*d_lat + c*d_lon'})})

    nc_file = gnss_root / 'collocs' / f'collocations_{year_to_do}_{subset}.nc'
    ds.to_netcdf(nc_file)
