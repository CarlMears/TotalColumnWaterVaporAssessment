

import argparse
import datetime
import os
from contextlib import suppress
from pathlib import Path
from netCDF4 import Dataset as netcdf_dataset
import numpy as np

from Era5_requests.era5_requests import era5_hourly_single_level_request


def time_interpolate_synoptic_maps(map_array, map_times, time_map):
    """Interpolate a day of gridded data.

    The map_array gridded data is interpolated, assuming it covers an entire
    day, including both the 00:00 and 24:00 time slots.

    The map_times array is a 1-d array of the times for each map in seconds
    since midnight.
    """
    sz = map_array.shape
    num_maps = sz[0]
    num_maps2 = len(map_times)

    if num_maps != num_maps2:
        raise ValueError("# of maps in map array not equal to number of map times")

    time_step = np.rint(map_times[1] - map_times[0])

    num_steps_in_day = len(map_times) - 1
    # hr_interp = np.copy(
    #     time_map / 60.0
    # )  # convert to hours after midnight ZZ without side effects

    # output array
    map_interp = np.full_like(map_array[0, :, :], np.nan)

    for interval in range(0, num_steps_in_day):
        ok = np.all(
            [(time_map >= map_times[interval]), (time_map < map_times[interval + 1])],
            axis=(0),
        )
        num_ok = np.sum(ok)
        if num_ok > 0:
            y = np.zeros((num_ok, 2))
            y[:, 0] = map_array[interval, :, :][ok]
            y[:, 1] = map_array[interval + 1, :, :][ok]

            wt_lower = (map_times[interval + 1] - time_map[ok]) / time_step

            assert np.all(
                wt_lower >= 0.0
            ), "wt_lower < 0.0"  # check to make sure weights are in bounds
            assert np.all(wt_lower <= 1.0), "wt_lower > 1.0"

            wt_upper = 1.0 - wt_lower

            y_interp = wt_lower * y[:, 0] + wt_upper * y[:, 1]
            map_interp[ok] = y_interp

    return map_interp  # ,test_map


class ERA5_Hourly():

    def __init__(self,
                 variable,
                 date_to_do,
                 temproot,
                 rss_grid=False):

        # Download ERA5 data from ECMWF for all 24 hours of day, and the first hour
        # of the next day.
        next_day = date_to_do + datetime.timedelta(hours=24)
        os.makedirs(temproot, exist_ok=True)
        file1 = era5_hourly_single_level_request(
            date=date_to_do,
            variable=variable[1],
            target_path=temproot,
            full_day=True,
            full_month=True,
            )

        # if next day is in the same month, this second request
        # refers to the same file, so no second download will be done
        file2 = era5_hourly_single_level_request(
            date=next_day,
            variable=variable[1],
            target_path=temproot,
            full_day=True,
            full_month=True,
        )

    # open the file(s), and combine the two files into a
    # 25-map array for the day being processed

        if date_to_do.month == next_day.month:
            hour_index1 = 24 * (date_to_do.day - 1)
            hour_index2 = hour_index1 + 25
            ds1 = netcdf_dataset(file1)
            var = ds1[variable[0]][hour_index1:hour_index2, :, :]

        else:
            # This is the case when the 25th hour is in the next month
            hour_index1 = 24 * (date_to_do.day - 1)
            hour_index2 = hour_index1 + 24
            ds1 = netcdf_dataset(file1)
            var_first_day = ds1[variable[0]][hour_index1:hour_index2, :, :]

            ds2 = netcdf_dataset(file2)
            var_next_day = ds2[variable[0]][0, :, :]
            var = np.concatenate(
                (var_first_day, var_next_day[np.newaxis, :, :]), axis=0
            )

       
        if rss_grid:
            var = np.flip(var, 1)
            temp = np.zeros((25, 721, 1441))
            temp[:, :, 0:1440] = var
            temp[:, :, 1440] = var[:, :, 0]
            var = (temp[:,0:720,0:1440]+
                    temp[:,0:720,1:1441]+
                    temp[:,1:721,0:1440]+
                    temp[:,1:721,1:1441])/4
        else:
            var = np.flip(var, 1)
        self.var = var
        self.var_name = variable[1]
        self.date = date_to_do


    def subset_7x7(self,ilat,ilon):

        try:
            var = self.var[:,ilat-3:ilat+4,ilon-3:ilon+4]
        except IndexError:
            var = np.full((25,7,7),np.nan)

        return var
    
    def interp_7x7(self,ilat,ilon,minutes):

        if minutes.shape == (7,7):
            var = self.subset_7x7(ilat,ilon)
            var_interp = np.zeros((7,7))
            var_interp = time_interpolate_synoptic_maps(var,
                                                        60.0*np.arange(0,25),
                                                        minutes)
            return var_interp
        else:
            raise ValueError('minutes must be 7x7 array')
    



