import datetime
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

from rss_sat_readers import read_rss_satellite_daily_xr,get_sat_def
from rss_plotting.global_map import plot_global_map
from read_NGL_vapor import read_NGL_vapor



subset = 'ocean'
for year_to_do in range(2023,2024):
    gnss_root = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/gnss_data')
    gnss_loc_file = gnss_root / f'gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.{subset}.txt'

    gnss_locs = pd.read_csv(gnss_loc_file, delim_whitespace=True)

    gnss_ids = gnss_locs['name']
    jan_1_1994 = datetime.date(1994,1,1)
    gnss_dict = {}
    first = True
    for i,gnss_id in enumerate(gnss_ids):
        print(i,gnss_id)
        df = read_NGL_vapor(gnss_id, year_to_keep = year_to_do)
        if len(df) == 0:
            continue
        year = np.round(df['yyyy'].values).astype(np.int32)
        month = np.round(df['mm'].values).astype(np.int32)
        day = np.round(df['dd'].values).astype(np.int32)
        num_days_since_1994 = np.zeros_like(year,dtype=np.int32)
        for i in range(len(year)):
            num_days_since_1994[i] = (datetime.date(year[i],month[i],day[i]) - jan_1_1994).days
        
        df['num_days_since_1994'] = num_days_since_1994
        df['gnss_id'] = gnss_id
        df.drop(columns=['yyyy','mm','dd'],inplace=True)
        if first:
            df_all = df.copy()
            first = False
        else:
            df_all = pd.concat([df_all,df],ignore_index=True)

    print(f'Year: {year_to_do} {len(df_all)} gnss observations')
    file_out = gnss_root / f'NGL_vapor_{year_to_do}_{subset}.csv'
    df_all.to_csv(file_out,index=False)
    print
        

