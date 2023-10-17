# Collocation Methods for Vapor Comparison

## Overview
For each location, retain a 7x7 pixel regions surronding the location. Retain the all parameters from the satellite file including
* Total Column Water Vapor (of course!)
* Total Column Cloud Water
* Rain Rate
* Wind Speed
* (if available) SST
* Observation time in day (I think this is more portable than other formats)
* Date of observation

## Finding Ocean Stations
For the first set of collocations, we want to use stations that are located on islands, so that the satellite observations surround the gnss station.  These stations are found, in the the context of AMSR2 observations using the following method. We start with the RSS standard 0.25 x 0.25 degree maps of geophysical retrievals (https://www.remss.com/missions/amsr/).

1. Make monthly maps of the fraction of AMSR2 observations with valid times that also have a valid vapor retrieval.
2. For each GNSS stations, examine the 7x7 (1.75 degree x 1.75 degree) region with the GNSS station in the center cell.
3. Separately examine the 3x7 (or 7x3) regions to the East, West, North, and South of the central cell.
4. For ALL of these regions, at least (frac_num_threshold) cells need to be valid at least (frac_valid_threshold).
5. For the set used here, frac_valid_threshold = 0.2 and frac_num_threshold = 10.0. This resulted in 137 stations. The stations are listed in gnss_data\gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.ocean.txt.

This is done in 
```
find_gnss_station_subsets.py
```
Note: this is a more liberal method (more stations are chosen) that we used in Mears, C. A., J. Wang, D. Smith, and F. J. Wentz, 2015, Intercomparison of Total Precipitable Water Measurements Made by Satellite-Borne Microwave Radiometers and Ground-Based GPS Instruments, Journal of Geophysical Research: Atmospheres, 120, 2492-2504.

## Reorganizing GNSS data
In order to keep the number of files simultaneously open during the collocation process, we reorganized the GNSS data into years files that contain all of the GNSS data from ocean stations in one file.  These are named with names like "NGL_vapor_2020_ocean.csv"  
This is done in 
```
assemble_GNSS_data_by_year.py
```

## Finding Collocations
For each daily satellite map, the GNSS stations in the chosen subset are looped over. For each station:
1. The 7x7 submap around the station is examined.
2. If less than frac_num_threshold (10) valid tcwv retrievals exist in the submap, no further processing is done.
3. Save GNSS tcwv and tcwv_sigma
4. Save satellite submaps for time,vapor,cloud water, wind speed, rain rate and SST.

If the satellite observation times are close to each other, then:

5. Interpolate ERA5 tcwv data to the satellite observation time.
6. Fit bilinear (a plane in vapor space) to the satellite tcwv, era5 tcwv and era5 tcwv subsampled at the locations where there is valid data.  This allows estimation of the tcwv value at the GNSS location.  Comparison of the era5 and era5_sampled case may allow estimation of the error introduced by the fit.  If the vapor is highly variable inside the 7x7 submap, the two ERA5 values may differ substantially, especially if the satellite data is relatively sparse.
7. Save the fit parameters and the fitted value at the GNSS location for all 3 fits.

Save all the data in yearly netcdf files.
collocations_2013_ocean.nc


## Questions going forward
* In many cases, the valid vapor values might not surround the location.  An obvious case is when the station is on the coast, but could also include cases near the edges of the swath, or near regions of moderate/heavy rain. The question is should we develop some sort of statistic to describe this?
* We fit a plane to the available observations to interpolate to the station locations.  Is this valid for coastal locations?  Can we use the sampled ERA5 fits to help decide this?
