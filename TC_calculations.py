#Module containing various functions to calculate important tropical cyclone quantities.
#Code by Rick Russotto, adapted from Jeff Strong's Matlab code.
#Started 5 October 2020

import numpy as np
import pandas as pd
import xarray as xr

#Do I need different functions for model storms vs. observations?
#How to handle the input function?

#For compatibility of variable names, input arrays of lats, lons, wind speed, 
#pressure, and whatever else, with storm and date-time dimensions; 
#all pre-processing should be done outside.

#Some calculations/plots are so simple they don't need a separate function here:
#plotting the tracks, plotting genesis points (?)

#These could also be subset by region before the calculation is done...
#but put region definitions here? 
#Or a region subsetting function?

#But what happens when a storm crosses between different regions?
#So, maybe should adopt Jeff's approach instead and apply a "tag" for the basin 
#at each storm/time grid point.
#This could be done in the pre-processing function 
#Regional subsetting won't work when lat-lon isn't a dimension coordinate for the dataset.
#Instead make region a separate variable, and mask by it...
#or make each region a Boolean variable because they might overlap--
#e.g. central vs. eastern Pacific, South Atlantic vs. Atlantic
#The region tagging should be done in notebooks that process the model and 
#observatoin data. 

#Might need these region tags in the calculations below. 

#Jeff calculated statistics by looping through the storms & times, 
#and appending to variables if they were in the correct basin.
#Could do this as well.


###   GENESIS DENSITY   ###
###   TRACK DENSITY   ###

#Input just a list of lats and lons?
#Jeff pre-allocated a grid and then, 
#at each grid point, did a count of the number of genesis points
#within a 10-x-10 degree box centerd on that point.
#Sort of a "convolution" approach.
#Have an argument for what the final grid size should be. 
#Jeff just masked out all the points outside the box, then 
#did a nansum to count the number not masked.
#Actually can have one function for genesis or tracks!
#Just input the lats and lons of either genesis points or entire tracks.
#Call it "storm density"
#Assuming lons input defined from -180 to 180.
#Jeff used a 1-degree grid.
#Does this need to be corrected for grid cell area? Doesn't look like Jeff did this
#Don't forget to normalize by number of years.
#Jeff might've used "storm days" for the track density--how does this work? 
#Divide by 4 to account for the 6-hourly time steps? 
#But, that can be accounted for outside the function. 

#Input storm_lats, storm_lons, grid_lats, grid_lons as NumPy arrays
def calcStormDensity(storm_lats, storm_lons, grid_lats, grid_lons, num_years):
    #Preallocate grid
    density_array = np.zeros((len(grid_lats), len(grid_lons)))
    index_lat = pd.Index(grid_lats, name='lat')
    index_lon = pd.Index(grid_lons, name='lon')
    density = xr.DataArray(density_array, coords=[index_lat, index_lon])
    
    #Loop through the grid points and count number within 5 degrees, both lat and lon
    #(Note, have to look periodically at least in the lon directoin)
    
    #Subset the longitudes by latitude--within 5 points 
    #(don't worry about boundaries--no TCs poleward of 85 degrees)
    for i in grid_lats:    
        if not(np.mod(i,10)):
            print('Calculating density--latitude ' + str(i))
        storm_lons_i = storm_lons[np.abs(storm_lats-i) <= 5]
        #Subset again, this time by longitude 
        #This one is tougher because have to include storms across the date line
        #Maybe just loop through the longitudes within 5 degrees of j, 
        #use an incrementing counter. Conditions different if within 5 degrees of 180
        for j in grid_lons:
            if j < -175:
                storm_lons_j_normal = storm_lons_i[np.abs(storm_lons_i-j) <= 5]
                overshoot_min = j - 5 + 360
                storm_lons_j_aux = storm_lons_i[storm_lons_i > overshoot_min]
                num_storms = len(storm_lons_j_normal) + len(storm_lons_j_aux)
            elif j > 175:
                storm_lons_j_normal = storm_lons_i[np.abs(storm_lons_i-j) <= 5]
                overshoot_max = j + 5 - 360 
                storm_lons_j_aux = storm_lons_i[storm_lons_i < overshoot_max]
                num_storms = len(storm_lons_j_normal) + len(storm_lons_j_aux)
            else: 
                storm_lons_j = storm_lons_i[np.abs(storm_lons_i-j) <= 5]
                num_storms = len(storm_lons_j)
            density.loc[dict(lat=i, lon=j)] = num_storms
    
    density_norm = density / num_years
    return density_norm


###   INTENSITY AND WIND STATISTICS   ###

###   ACE   ###
#Ace is the sum of the squared max. sustained wind speeds over a storm's lifetime, 
#divided by 10000.
#Traditionally only winds at or above 35 knots count, but can also have 
#"modified ACE" (MACE) in which all winds count (so have an argument for min winds that count).
#Sum up ACE for all storms to get for a season.

#What to input? 
#How about just an array of all the max. wind values, already with nans taken out and flattened
#(since boundaries between storms don't matter).
#Can subset outside the function.
#Although, Jeff did a cumulative time evolution of ACE as well...
#but just seasonal in each month would make more sense. 
#Can subset winds in time as well as space outside the function.
#And min. winds

#Note: Zhao winds might be in 

#In fact having this function outside the scripts might unnecessarily make things unreadible.
#But at least combine 2 lines into one. 
#Should be 1-D NumPy array input, not XArray.
def calcACE(wind_array_knots, min_wind):
    wind_array_floored = wind_array_knots[wind_array_knots >= min_wind]
    ACE = np.sum(wind_array_floored*wind_array_floored)/1.0e4
    return ACE

###   PDI   ###

#PDI is basically the same as ACE, but cubed and divided by 10^6 instead of 10^4.
def calcPDI(wind_array_knots, min_wind):
    wind_array_floored = wind_array_knots[wind_array_knots >= min_wind]
    PDI = np.sum(wind_array_floored*wind_array_floored*wind_array_floored)/1.0e6
    return PDI
