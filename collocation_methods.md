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

  ### Considerations
  * how many valid vapor pixels are required to keep the data?

## "Spatial Quality"
* In many cases, the valid vapor values might not surround the location.  An obvious case is when the station is on the coast, but could also include cases near the edges of the swath, or near regions of moderate/heavy rain. The question is should we develop some sort of statistic to describe this?
* In previous work, we fit a plane to the available observations to interpolate to the station locations.  Is this valid for coastal locations?
