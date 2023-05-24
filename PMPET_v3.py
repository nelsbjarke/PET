import xarray as xr
import numpy as np
import os
import sys

CP = 1.004 * 10**-3   #Specific heat of air at constant pressure (C_p)in MJ/(kg K)
LV = 2.45             # Latent Heat of Vaporization for Water.  FAO uses constant 2.45 MJ/kg

# FAO56 Constants 
CN_GRASS =  900     # define reference values for grass (short crop)
CD_GRASS =  0.34    # units: s/m
CN_ALFALFA =  1600  # define reference values for alfalfa (tall crop)  
CD_ALFALFA =  0.38  # units:   s/m

def get_parent_attributes(parent_temp_path=None):
    '''Returns a string that contains all of the parent attributes for carry over into the published netcdf files'''

    # Get the monthly temperature dataset for the GCM to grab the parent attributes
    ds_parent = xr.open_dataset(parent_temp_path)

    # Create a blank string for the storage of the parent attributes
    parent_attrs_string = '{'

    # Loop through the attributes that need to be added from the parent dataset
    for attr_id in ds_parent.attrs:
        parent_attrs_string = parent_attrs_string+"'"+attr_id+"':'"+str(ds_parent.attrs[attr_id])+"'"

        if attr_id != list(ds_parent.attrs.keys())[-1]:
            parent_attrs_string = parent_attrs_string + ','

    # Append a '}' to signify the end of the dictionary for the parent attribute string
    parent_attrs_string = parent_attrs_string+'}'

    return parent_attrs_string
    
def keep_list_of_attributes(
                            ds,
                            keep_attrs_list
                             ):
    '''Keep all attributes in the list and then delete the rest'''

    # Loop through the attributes and delete unless in the keep list
    for attr_id in list(ds.attrs.keys()):
        if attr_id not in keep_attrs_list:
            del ds.attrs[attr_id]

    return ds

    
def add_attributes(
                    ds,
                   attr_dict
                   ):
    '''Adds the attributes from a dictionary that is generated outside of this function to the netcdf file'''

    # Loop through the attributes that need to be added
    for attr_id in attr_dict:
        ds.attrs[attr_id] = attr_dict[attr_id]

    return ds
    
def add_fill_miss_values(ds,
                         varID):
    '''Adds the fill value and missing values to the esxisting encoding of the netcdf file'''

    # Add the encoding for the attributes as NaN values
    ds[varID].encoding['_FillValue'] = np.nan
    ds[varID].encoding['missing_value'] = np.nan

    return ds

def scale_wind_FAO56(u_z,z):
    '''Calculate 2-meter wind from wind at another height using FAO56 log-layer scaling'''
    import numpy as np
    print('Scaling wind from height of ' + str(z) + ' m to 2m')
    u_2 = 4.87*u_z/np.log(67.8*z - 5.42)
    return u_2

def get_rnet(path_dict,variable_dict,rnet_calculated=False,units='W/m2'):

    # If rnet_calculated=True, read in Sensible and Latent heat fluxes and calculate
    #   Rnet - G from the surface energy balance: Rnet - G = SH + LH.
    #   It is assumed that SH and LH are positive-upwards so that Rnet-G is positive downwards
    #   otherwise read rnet in from a file.  Specifying a separate Ground heat flux input file is not supported.
    #
    # Unit conversion factors (FAO56 formulas require units of MJ/m^2/day)
    #    Most climate models use W/m^2.   The ERA5 reanalysis has units of J/m^2 for a specified time interval,
    #    usually a day.  See Documentation at ECMWF for more informatino,  .  
    #
    if units == 'W/m2':
        rnet_conv = 86400.0/10**6
    elif units == 'J/m2/day':
        rnet_conv = 1/10**6

    if rnet_calculated:
        sh   = xr.open_dataset(path_dict['sh_path'])
        lh   = xr.open_dataset(path_dict['lh_path'])
        rnet_temp = rnet_conv*(sh[variable_dict['sh']] + lh[variable_dict['lh']])   #we convert from W/m2 to MJ/day
        rnet=rnet_temp.to_dataset(name=variable_dict['rnet'])
    else:
        rnet = xr.open_dataset(path_dict['rnet_path'])
        rnet[variable_dict['rnet']] = rnet_conv*rnet[variable_dict['rnet']]
    return rnet

def PenMon(path_dict = None,variable_dict=None,veg_type='grass',rnet_calculated = True,rh_from_tdew=False,wind_meas_height=10,const_wind = 1.0,rnet_units='W/m2',latname='lat',lonname='lon'):
        '''Calculate FAO56 Penman Monteith  reference (potential) evapotranspiration for the given dictionary of input files'''
#
# rnet_calculated = True(default) if rnet-G is calculated from SH + LH using the Surface Energy Balance
# veg_type=alfalfa, grass (default) 
# wind_meas_height = 2.0 (default).  wind measurement height -- if different from 2m then log-layer scaling is used
# variable_dict contains the names of the different variables used in each netCDF input file
# path_dict contains path to the input files


    # look up some variable names in the variable dictionary

    temp_varname = str(variable_dict['temp'])  
    rh_varname   = str(variable_dict['rh'])
    rnet_varname = str(variable_dict['rnet'])
    ps_varname = str(variable_dict['ps'])

    # Read in data as Xarray Datasets

    temp = xr.open_dataset(path_dict['temp_path'])
    ps   = xr.open_dataset(path_dict['ps_path'])
    rnet = get_rnet(path_dict,variable_dict,rnet_calculated=rnet_calculated,units=rnet_units)

    #Reads in windspeed from a file.  If the file is not found, sets windspeed to constant (default 1.0m/s)

    if os.path.isfile(path_dict['wind_path']):
        wind_type = 'Model Average Windspeed'
        wind = xr.open_dataset(path_dict['wind_path'])
        wind=wind.assign_coords({lonname:(temp[lonname]),
                         latname:(temp[latname])})
        windspeed = wind[variable_dict['wind']]
        #rescale wind if not at 2 meters
        if wind_meas_height != 2:
            windspeed=scale_wind_FAO56(windspeed,z=wind_meas_height)
    else:
        windspeed = const_wind
        wind_type = 'Constant: ' + str(const_wind) + ' m/s'
        print('Wind File Not Found:  Setting windspeed to ' + str(const_wind) + ' m/s')

    # copy over lats and lons because small differences in these mess up xarray computations.
    rnet=rnet.assign_coords({lonname:(temp[lonname]),
                             latname:(temp[latname])})

    ps=ps.assign_coords({lonname:(temp[lonname]),
                         latname:(temp[latname])})

    # Convert to Celsius if temperature is in Kelvin
    if temp[temp_varname].attrs['units'] == 'K':
        temp[temp_varname]=temp[temp_varname]-273.15
        print('Converting Temperature from K to C')

    # Calculate psychrometric constant from surface pressure
    psychromet = CP*ps[ps_varname]/(1000*LV*.622)   # Units: kPa/K  factor of 1000 is to convert PS in Pa to KPa. Units: kPa/K

    # Specify constants for FAO56 equation for chosen reference vegetation type
    if veg_type == 'grass':
        Cn = CN_GRASS
        Cd = CD_GRASS
    elif veg_type == 'alfalfa':
        Cn = CN_ALFALFA
        Cd = CD_ALFALFA
    else:
        sys.exit('veg_type not defined')

    # Calculate saturation vapor pressure
    temp['esat'] = 0.6108*np.exp( (17.27*temp[temp_varname] )/( temp[temp_varname]+237.3 ))    #Units: kPa

    # Calculate Vapor Pressure Deficit (VPD) either using dewpoint temperature(tdew) or from relative humidity (rh) 
    if rh_from_tdew == True:
        tdew_varname = str(variable_dict['tdew'])
        tdew = xr.open_dataset(path_dict['tdew_path'])
        tdew=tdew.assign_coords({lonname:(temp[lonname]),
                         latname:(temp[latname])})
        if tdew[tdew_varname].attrs['units'] == 'K':
            tdew[tdew_varname]=tdew[tdew_varname]-273.15
            print('Converting Temperature from K to C')
        # Calculate the saturation vapor pressure at the dewpoint temperature, which is the  vapor pressure at the actual  air temperature
        tdew['esat'] = 0.6108*np.exp( (17.27*tdew[tdew_varname] )/( tdew[tdew_varname]+237.3 ))      #Units: kPa
        # Calculate the  vapor pressure deficit from the difference of vapor pressures
        temp['VPD'] = temp['esat']-tdew['esat']
        humidity_source = 'Dewpoint Temperature'  

    else:
        rh_varname   = str(variable_dict['rh'])
        rh   = xr.open_dataset(path_dict['rh_path'])
        rh=rh.assign_coords({lonname:(temp[lonname]),
                             latname:(temp[latname])})
        # Calculate the vapor pressure deficit (VPD)  from the relative humitdity and saturation vapor pressure
        temp['VPD'] = (1.0 -  rh[rh_varname] / 100.0 ) * temp['esat']        # #Units: kPa
        humidity_source = 'Relative Humidity'  

    # Calculate delta (deriv. of esat w.r.t temperature)
    temp['delt'] = 4098*temp['esat']/((temp[temp_varname]+237.3)**2)         # Units: kPa/K
  

    # Calculate PET using the Penman-Monteith formulations from FAO Pub. 56 (FAO56)

    temp['PET'] = (0.408*temp['delt']*rnet[rnet_varname]+psychromet*Cn*windspeed*temp['VPD']/(temp[temp_varname]+273))/\
          (temp['delt']+psychromet*(1 + Cd*windspeed))   # Units: mm/day

    # Separately calculate the advection and energy components separately
    temp['PETen'] = (0.408*temp['delt']*rnet[rnet_varname])/(temp['delt']+psychromet*(1 + Cd*windspeed))
    temp['PETad'] = (psychromet*Cn*windspeed*temp['VPD']/(temp[temp_varname]+273))/(temp['delt']+psychromet*(1 + Cd*windspeed))

    # Define variable attributes
    temp['PET'].attrs['Units'] = 'mm day-1'
    temp['PET'].attrs['short_name'] = 'PET'
    temp['PET'].attrs['long_name'] = 'Penman-Monteith Potential Evapotranspriation'
    temp['PET'].attrs['Description'] = 'Penman-Monteith Potential Evapotranspiration computed using FAO56 equations (1) and (3) for vegetation type = ' + veg_type
    temp['PET'].attrs['Wind Treatment'] = wind_type

    temp['PETen'].attrs['long_name'] = 'Energy Term of Penman-Monteith Potential Evapotranspriation'
    temp['PETen'].attrs['Description'] = 'Energy Term  Penman-Monteith Potential Evapotranspiration computed using FAO56 equations (1) and (3) for vegetation type = ' + veg_type

    temp['PETad'].attrs['long_name'] = 'Advective Term of Penman-Monteith Potential Evapotranspriation'
    temp['PETad'].attrs['Description'] = 'Advective Term of Penman-Monteith Potential Evapotranspiration computed using FAO56 equations (1) and (3) for vegetation type = ' + veg_type

    temp['VPD'].attrs['Units'] = 'kPa'
    temp['VPD'].attrs['long_name'] = 'Vapor Pressure Deficit'
    temp['VPD'].attrs['short_name'] = 'VPD'
    temp['VPD'].attrs['Description'] = 'Vapor Pressure Deficit computed from monthly 2m-temperature and 2m relative humidity'
    temp['VPD'].attrs['Humidity Source'] = humidity_source

    PETout = temp['PET'].to_dataset(name = 'PET',promote_attrs=True)
    VPDout = temp['VPD'].to_dataset(name = 'VPD')

    PETenout = temp['PETen'].to_dataset(name = 'PET',promote_attrs=True)
    PETadout = temp['PETad'].to_dataset(name = 'PET',promote_attrs=True)

    # copy over attributes
    PETout.attrs=temp.attrs
    VPDout.attrs=temp.attrs

    return PETout,VPDout,PETenout,PETadout

def PriestleyTaylor(path_dict = None,variable_dict=None,rnet_calculated = True,wind_meas_height=10,const_wind = 1.0,alpha=1.26,rnet_units='W/m2',latname='lat',lonname='lon'):
    '''Calculate Priestley-Taylor potential evapotranspiration for the given dictionary of input files'''

# rnet_calculated = True(default) if rnet-G is calculated from SH + LH using the Surface Energy Balance
# veg_type=alfalfa, grass (default) 
# variable_dict contains the names of the different variables used in each netCDF input file
# path_dict contains path to the input files

# look up some variable names in the variable dictionary

    temp_varname = str(variable_dict['temp'])
    rnet_varname = str(variable_dict['rnet'])
    ps_varname = str(variable_dict['ps'])

    # Read in data as Xarray Datasets

    temp = xr.open_dataset(path_dict['temp_path'])
    ps   = xr.open_dataset(path_dict['ps_path'])

    rnet = get_rnet(path_dict,variable_dict,rnet_calculated=rnet_calculated,units=rnet_units)

    # copy over lats and lons because small differences in these mess up xarray computations
    rnet=rnet.assign_coords({lonname:(temp[lonname]),
                             latname:(temp[latname])})

    ps=ps.assign_coords({lonname:(temp[lonname]),
                         latname:(temp[latname])})

    # Convert to Celsius if temperature is in Kelvin
    if temp[temp_varname].attrs['units'] == 'K':
        temp[temp_varname]=temp[temp_varname]-273.15
        print('Converting Temperature from K to C')

    # Calculate psychrometric constant from PS
    psychromet = CP*ps[ps_varname]/(1000*LV*.622)   #Units: kPa/K  factor of 1000 is to convert PS in Pa to KPa. Units: kPa/K

    # Calculate saturation vapor pressure
    temp['esat'] = 0.6108*np.exp( (17.27*temp[temp_varname] )/( temp[temp_varname]+237.3 ))  # T (Celsius), esat (kPa)

    # Calculate delta (deriv. of esat w.r.t temperature)
    temp['delt'] = 4098*temp['esat']/((temp[temp_varname]+237.3)**2)         # kPa/K

    # Calculate PET using the Priestley-Taylor formulation with constant alpha (defalut alpha = 1.26)

    temp['PET'] = (alpha*temp['delt']*rnet[rnet_varname])/(LV*(temp['delt']+psychromet))  # Units: mm/day

    # Define variable attributes
    temp['PET'].attrs['Units'] = 'mm day-1'
    temp['PET'].attrs['short_name'] = 'PET'
    temp['PET'].attrs['long_name'] = 'Priestley-Taylor Potential Evapotranspriation'
    temp['PET'].attrs['Description'] = 'Priestley-Taylor Potential Evapotranspiration computed using formulation described in Priestley & Taylor (1972) with alpha = '+str(alpha)

    PETout = temp['PET'].to_dataset(name = 'PET',promote_attrs=True)

    # copy over attributes
    PETout.attrs=temp.attrs

    return PETout
