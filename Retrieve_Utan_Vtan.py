# Script adapted from Storm_Centered_Plots.py to retrieve storm-centered
# grids for u_ref and v_ref, and calculate radial and tangential components
# (u_rad and u_tan) from them.
# Written on Flurry (to have everything in one place), 
# to be transferred to Discover and run there.
# Based on Jeff's scripts.
# Save NetCDF files that can be transferred back to Flurry for additional analysis.

# Started 2020-12-07

#Import stuff
import xarray as xr
import pandas as pd
import numpy as np

#Load the model data saved by Read_Zhao_TCs.ipynb
#(for pointing to where the storms are within the 2D model output)
ds_tracks_v1 = xr.open_dataset('zhao_tracks_v1.nc')
ds_tracks_v2 = xr.open_dataset('zhao_tracks_v2.nc')

#Add variables relevant to peak intensity
ds_tracks_v1['max_wind_kts'] = ds_tracks_v1['wind_kts'].max(dim='date_time')
ds_tracks_v2['max_wind_kts'] = ds_tracks_v2['wind_kts'].max(dim='date_time')
ds_tracks_v1['min_pressure'] = ds_tracks_v1['pressure'].min(dim='date_time')
ds_tracks_v2['min_pressure'] = ds_tracks_v2['pressure'].min(dim='date_time')
ds_tracks_v1['argmax_wind_kts'] = ds_tracks_v1['wind_kts'].argmax(dim='date_time')
ds_tracks_v2['argmax_wind_kts'] = ds_tracks_v2['wind_kts'].argmax(dim='date_time')

#Easiest way to do this would be to just do it at peak intensity,
#then have a 3D NetCDF file for all the storms, 
#with each variable.

#Number of storms
num_storms = {'v1': len(ds_tracks_v1.storm), 
              'v2': len(ds_tracks_v2.storm)}

#Preallocate arrays for each variable and storm
#One key per variable; How big is the box? 15x15 degree
#(assuming everything already regridded to 0.5 degrees)
#(This is what Jeff assumes, with -30 to 30 indices. So 61x61 grid)
dict_arrays = {'v1': {'u_ref': np.zeros((61, 61, num_storms['v1'])), 
                      'v_ref': np.zeros((61, 61, num_storms['v1']))},
               'v2': {'u_ref': np.zeros((61, 61, num_storms['v2'])), 
                      'v_ref': np.zeros((61, 61, num_storms['v2']))}}
               
#Paths to the NetCDF files
paths = {'v1': '/discover/nobackup/jdstron1/JS_C180_20YearTest/tc_data/',
         'v2': '/discover/nobackup/jdstron1/JS_C180v2_20YearTest/tc_data/'}
               
#Loop through each storm,
#read the NetCDF file for the time where it had the peak intensity, 
#and fill in that entry in a 61x61 grid in the dicts
               
ds_tracks = {'v1': ds_tracks_v1, 
             'v2': ds_tracks_v2}
               
# #What the files look like when loaded on ipython on Discover:
# <xarray.DataArray 'slp' (time: 1460, lat: 360, lon: 720)>
# [378432000 values with dtype=float32]
# Coordinates:
#   * time     (time) object 1995-01-01 06:00:00 ... 1996-01-01 00:00:00
#   * lon      (lon) float32 -179.75 -179.25 -178.75 ... 178.75 179.25 179.75
#   * lat      (lat) float32 -90.0 -89.25 -88.75 -88.25 ... 88.25 88.75 89.25 90.0
# Attributes:
#     units:      mb
#     long_name:  sea-level pressure
               
#OK so this is in -180 to 180 longitudes, just like the storm tracks.

#Stuff for figuring out the indices:
month_lengths_0 = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]
month_start_indices = np.cumsum(month_lengths_0)*4
               
for v in ['v1', 'v2']:
    print('Obtaining storm-centered grid at peak intensity for version: ' + v)
    for var in dict_arrays[v].keys():
        print('Looping through years for variable: ' + var)
        opened_year=0
        for i in np.arange(num_storms[v]): 
            print(str(i)+'th storm')
            #Find the year of peak intensity to get the right file
            #print(ds_tracks[v]['year'])
            #print(ds_tracks[v]['argmax_wind_kts'])
            year_peak = ds_tracks[v]['year'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            #Find the month, day, hour to figure out the index (Jeff did something similar)
            month_peak = ds_tracks[v]['month'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            day_peak = ds_tracks[v]['day'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            hour_peak = ds_tracks[v]['hour'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            index_peak = int(month_start_indices[int(month_peak-1)]+(day_peak-1)*4+hour_peak/6-1)
            #There is an offset by 1 on 0,6,12,18 in the track data vs. 6,12,18,0 in these files.
            #If a storm peaked in the very first observation of the year: look in the previous year's file,
            #which includes it. (Should I ask Jeff about this?)
            if index_peak == -1:
                year_peak = year_peak -1 
            str_year = str(int(year_peak))
            #Find the lat and lon where the storm currently is
            lat_peak = ds_tracks[v]['lat'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            lon_peak = ds_tracks[v]['lon'].isel(storm=i, date_time=ds_tracks[v]['argmax_wind_kts'].isel(storm=i))
            #Open the file and subset to peak time--IF not opened already
            if not year_peak == opened_year:
                print('Opening NetCDF file for variable: ' + var + ' and year: ' + str_year)
                ds_gridded = xr.open_dataset(paths[v]+'atmos.'+str_year+'010100-'+str_year+'123123.'+var+'.nc')
                print('File opened')
                opened_year = year_peak
            ds_gridded_peak_time = ds_gridded.isel(time=index_peak)
            #Reduce to a 15x15 degree (61x61 pixel) grid around the storm, following how Jeff did this
            #How to deal with the date line?
            #Just use indices--Python will take care of the boundaries using negative indices
            #(No, actually it won't have to concatenate)
            #Find the lat/lon indices on the grid closest to where the storm is
            #print(lat_peak)
            #print(ds_gridded_peak_time['lat'])
            lat_index = np.argmin(np.abs(lat_peak.data-ds_gridded_peak_time['lat'].data))
            lon_index = np.argmin(np.abs(lon_peak.data-ds_gridded_peak_time['lon'].data))
            #Extract the (assuming lat, lon dimension order--check this: yes, that's the case)
            len_lat = len(ds_gridded_peak_time['lat'])
            len_lon = len(ds_gridded_peak_time['lon'])
            #Account for lon within 15 degrees of date line (don't have to worry for lat) by concatenating
            if lon_index < 30:
                print('lon_index < 30')
                ds_subset = np.concatenate([ds_gridded_peak_time[var].data[lat_index-30:lat_index+31,lon_index-30:],
                                            ds_gridded_peak_time[var].data[lat_index-30:lat_index+31,:lon_index+31]], axis=1)
            elif (len_lon-1) - lon_index < 30:
                print('len_lon - lon_index < 30')
                ds_subset = np.concatenate([ds_gridded_peak_time[var].data[lat_index-30:lat_index+31,lon_index-30:],
                                            ds_gridded_peak_time[var].data[lat_index-30:lat_index+31,:lon_index+31-len_lon]],axis=1)
            else:
                print('middle lons. lat_index: ' + str(lat_index) + ', lon_index: ' + str(lon_index))
                ds_subset = ds_gridded_peak_time[var].data[lat_index-30:lat_index+31, 
                                                           lon_index-30:lon_index+31]
               
            #Fill in the layer in the preallocated arrays
            dict_arrays[v][var][:,:,i] = ds_subset
            
            
#Calculate tangential and radial velocity and add to the dicts
#First define grid in polar coordinates

#(Comented this out because had error--now corrected--but didn't want it to fail again 
# # after taking > 1 hour to run. Can do this on Flurry instead)
# print('Defining polar coordinates grid')

# x_cart = np.linspace(-15, 15, 61)
# y_cart = np.linspace(-15, 15, 61)

# x_grid, y_grid = np.meshgrid(x_grid, y_grid)

# r_grid = np.sqrt(x_grid*x_grid + y_grid*y_grid)
# theta_grid = np.arctan2(y_grid, x_grid)

# print('r, theta from grid')
# print(r_grid)
# print(theta_grid)

# #Now calculate tangential and radial components and add to the dicts
# for v in ['v1', 'v2']:
#     dict_arrays[v]['u_rad'] = (dict_arrays[v]['u_ref']*np.cos(theta_grid) +
#                                dict_arrays[v]['v_ref']*np.sin(theta_grid))
#     dict_arrays[v]['u_tan'] = (dict_arrays[v]['v_ref']*np.cos(theta_grid) - 
#                                dict_arrays[u]['u_ref']*np.sin(theta_grid))
    
                                            
#Print some test output
for v in ['v1','v2']:
    for var in dict_arrays[v].keys():
        print('For ' + v + ', variable ' + var +':')
        print('Array shape is: ')
        print(np.shape(dict_arrays[v][var]))
                                            
               
#Construct XArray datasets from the arrays and save
for v in ['v1', 'v2']:
    da_dict=dict()
    for var in dict_arrays[v].keys():
        da_dict[var] = xr.DataArray(dict_arrays[v][var], name=var,
                                    coords=[pd.Index(np.linspace(-15, 15, 61), name='lat'), 
                                            pd.Index(np.linspace(-15, 15, 61), name='lon'),
                                            pd.Index(np.arange(num_storms[v]), name='storm')]) 
    ds_merged = xr.merge(da_dict.values())
    ds_merged.to_netcdf('storm_centered_reference_winds_'+v+'.nc')
                                            
               
               

               

               
               
