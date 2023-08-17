# TotalColumnWaterVaporAssessment
This is a repository to define methods for an assessment of TCWV products.  In this document we discuss the general approach.  A specific comparison plan is contained in ComparisonPlan.md

## Products to be included
The types of water vapor products included in the study are:
* Ground Based GPS/GNSS.  These are almost always located at fixed locations on land.
* Space-Based GPS/GNSS (e.g. COSMIC, COSMIC-2).  These are located approximately globally over land and ocean.
* Microwave Radiometer.  These are only available over ice-free oceans.
* Reanalysis.  Available everywhere with various spatial and temporal resolution
  
## Time scales and collocations
* "exact" collocations, i.e. observations within defined location and temporal bounds of each other.  The advantage of this is that sampling errors are mostly eliminated, but there are likely to be fewer observations, increasing random noise.
* Maps of time-averaged observations (e.g. monthly mean maps).  These can have substantial sampling errors, both due to sampling within the time period, and diurnal sampling.
* Large spatial scale averages (e.g. global or latitude-band averaged).  Again averages from different products can have different sampling.

## Statisical Methods to consider
* summary statistics of between-product differences (e.g. mean, std. dev., rms diff.)
* trimmed versions of the above, to remove large outliers
* binned statistics with various binning variables (TCWV itself, and any interesting confounding variables)
* discountinuity detection in time series


