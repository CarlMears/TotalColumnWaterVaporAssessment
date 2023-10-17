import datetime
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

from rss_sat_readers import read_rss_satellite_daily_xr,get_sat_def
from rss_plotting.global_map import plot_global_map

from get_era5_2d import ERA5_Hourly,hour_weights

frac_valid_threshold = 0.2
frac_num_threshold = 10.0
max_collocs = 10000
subset = 'ocean'
#subset = 'ocean_coastal'


    

#read olivier's list of gnss locations
gnss_root = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/gnss_data')
gnss_loc_file = gnss_root / 'gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.txt'

gnss_locs = pd.read_csv(gnss_loc_file, delim_whitespace=True)

sat_name = 'AMSR2'
sat_info = get_sat_def(sat_name=sat_name)

sat_root = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/sat_data')
nc_outfile = sat_root / 'util' / f'{sat_name}_fraction_valid.nc'
fraction_valid = xr.open_dataset(nc_outfile)
fraction_valid = fraction_valid['fraction_valid'].values

#loop over the gnss locations and keep the ones that have enough valid vapor data
#for at least one month
#enough is at least frac_num_thresholds vapor fractions > frac_valid_threshold
#in a 7x7 grid centered on the gnss location

lons_test = np.arange(0.125,360,0.25)
lats_test = np.arange(-90.0+0.125,90.0,0.25)
lons_grid = np.tile(lons_test,(720,1))
lats_grid = np.tile(lats_test,(1440,1)).T

gnss_locs.reset_index(inplace=True)
for index,row in gnss_locs.iterrows():
    ilat = int(np.floor((4.0*(row['lat']+90.0))))
    ilon = int(np.floor((4.0*(row['lon']))))
    keep_this_loc = False
    if (ilat > 4) or (ilat < 716):
        #print(f'lat={row["lat"]},lon={row["lon"]},ilat={ilat},ilon={ilon}')
        #print(fraction_valid[ilat,ilon,:])
        ilats = np.arange(ilat-3,ilat+4)
        ilons = np.arange(ilon-3,ilon+4)
        ilons[ilons<0] += 1440
        ilons[ilons>=1440] -= 1440
        ilat_grid = np.tile(ilats,(7,1)).T
        ilon_grid = np.tile(ilons,(7,1))

        lon_test = lons_grid[ilat_grid,ilon_grid]
        lat_test = lats_grid[ilat_grid,ilon_grid]
        print(f'lat={row["lat"]},lon={row["lon"]},center_lat={lat_test[3,3]},center_lon={lon_test[3,3]}')

        frac_valid_near = fraction_valid[ilat_grid,ilon_grid,:]
        if subset == 'ocean_coastal':
            if np.sum(np.isfinite(frac_valid_near)) > 0:
                for month_index in np.arange(12):
                    frac_valid_near_month = frac_valid_near[:,:,month_index]
                    ok = np.all([frac_valid_near_month > frac_valid_threshold,np.isfinite(frac_valid_near_month)],axis=0)
                    if np.sum(ok) >= frac_num_threshold:
                        #print(f'lat={row["lat"]},lon={row["lon"]},month={month_index+1},frac_valid={np.sum(ok)/49.0}')
                        keep_this_loc = True
        elif subset == 'ocean':
            if np.sum(np.isfinite(frac_valid_near)) > 0:
                for month_index in np.arange(12):
                    frac_valid_near_month = frac_valid_near[:,:,month_index]
                    side1 = frac_valid_near_month[0:4,:]
                    side2 = frac_valid_near_month[4:7,:]
                    side3 = frac_valid_near_month[:,0:4]
                    side4 = frac_valid_near_month[:,4:7]
                    keep_this_loc_this_month = True
                    for side in [side1,side2,side3,side4]:
                        ok = np.all([side > frac_valid_threshold,np.isfinite(side)],axis=0)
                        if np.sum(ok) < frac_num_threshold:
                            #print(f'lat={row["lat"]},lon={row["lon"]},month={month_index+1},frac_valid={np.sum(ok)/49.0}')
                            keep_this_loc_this_month = False
                    if keep_this_loc_this_month:
                        keep_this_loc = True
        else:
            raise ValueError('subset must be ocean or ocean_coastal')

    if keep_this_loc is False:
        gnss_locs.drop(index=index,inplace=True)

gnss_locs.reset_index(inplace=True)
gnss_locs = gnss_locs.drop(columns=['level_0', 'index'])
gnss_loc_file = gnss_root / f'gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.{subset}.txt'
gnss_locs.to_csv(gnss_loc_file,sep=' ',index=False)


zero_map = np.zeros((720,1440))
fig,ax = plot_global_map(zero_map,vmin=-1,vmax=1,cmap='BrBG',plt_colorbar=True,plt_coastlines=True)


# plot a map of the remaining locations
for index,row in gnss_locs.iterrows():
    lon = row['lon']
    if lon > 180:
        lon -= 360
    ax.plot(lon,row['lat'],'o',markersize=3,color='blue')
png_file = gnss_root / f'gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.{subset}.png'
fig.savefig(png_file,dpi=300)

