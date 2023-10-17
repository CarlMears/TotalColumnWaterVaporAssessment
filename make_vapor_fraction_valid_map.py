import datetime
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from rss_plotting.global_map import plot_global_map

from rss_sat_readers import read_rss_satellite_daily_xr,get_sat_def

sat_name = 'AMSR2'
sat_info = get_sat_def(sat_name=sat_name)

start_date_str = sat_info['start_date'].split('/')
start_date = datetime.date(int(start_date_str[0]),
                           int(start_date_str[1]),
                           int(start_date_str[2]))

if sat_info['end_date'] == 'present':
    end_date = datetime.date.today()
else:
    end_date_str = sat_info['end_date']
    end_date = datetime.date(int(end_date_str[0]),
                             int(end_date_str[1]),
                             int(end_date_str[2]))

end_date = start_date + datetime.timedelta(days=65)
date_to_do = start_date

start_date = datetime.date(2013, 1, 1)
end_date = datetime.date(2022, 12, 31)

tot_vapor_num_valid = np.zeros((720,1440,12),dtype=np.int32)
tot_time_num_valid = np.zeros((720,1440,12),dtype=np.int32)
while date_to_do <= end_date:
    # read in satellite data
    try:
        d = read_rss_satellite_daily_xr(year=date_to_do.year,
                                        month=date_to_do.month,
                                        day=date_to_do.day,
                                        satellite=sat_name,
                                        verbose=True)

        for node in [0,1]:
            time_map = d['time'][node,:,:].values
            ok = np.all([np.isfinite(time_map),time_map>0],axis=0)
            tot_time_num_valid[:,:,date_to_do.month-1][ok] += 1
            vapor_map = d['vapor'][node,:,:].values
            ok = np.isfinite(vapor_map)
            tot_vapor_num_valid[:,:,date_to_do.month-1][ok] += 1
    except FileNotFoundError:
        print('No data for ',date_to_do)

    date_to_do += datetime.timedelta(days=1)



fraction_vapor_valid = tot_vapor_num_valid / tot_time_num_valid

month_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


plot_global_map(fraction_vapor_valid[:,:,6],
                title=sat_name+' fraction valid '+month_names[6],
                cmap='jet',
                vmin=0,vmax=1,
                plt_colorbar=True,
                plt_coastlines=True)

plot_global_map(fraction_vapor_valid[:,:,7],
                title=sat_name+' fraction valid '+month_names[6],
                cmap='jet',
                vmin=0,vmax=1,
                plt_colorbar=True,
                plt_coastlines=True)
plt.show()

ds_fraction_valid = xr.Dataset({'fraction_valid': (['latitude', 'longitude','month'],  fraction_vapor_valid)},
                                 coords={'latitude': (['latitude'], d['latitude'].values),
                                        'longitude': (['longitude'], d['longitude'].values),
                                        'month': (['month'], np.arange(1,13))})
ds_fraction_valid.attrs['units'] = 'fraction'
ds_fraction_valid.attrs['satellite'] = sat_name
ds_fraction_valid.attrs['start_date'] = start_date.strftime('%Y-%m-%d')
ds_fraction_valid.attrs['end_date'] = end_date.strftime('%Y-%m-%d')
ds_fraction_valid.attrs['description'] = 'fraction of valid data for '+sat_name
ds_fraction_valid.attrs['created_on'] = datetime.datetime.now().isoformat()
ds_fraction_valid.attrs['created_by'] = 'Carl Mears, RSS'
ds_fraction_valid.attrs['contact'] = 'Carl Mears (mears@remss.com)'

sat_root = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/sat_data')
nc_outfile = sat_root / 'util' / f'{sat_name}_fraction_valid.nc'
os.makedirs(os.path.dirname(nc_outfile),exist_ok=True)
ds_fraction_valid.to_netcdf(nc_outfile)
print(f'Wrote {nc_outfile}')
print
