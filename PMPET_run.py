#  PET calculation runs for CMIP6 ensemble

import PMPET_v3 as pm
import os
import datetime
import glob

def get_monthly_gcm_path(gcm,
                     expid,
                     varID):
    '''Get the monthly data for a VarID from the directory of gcm data'''
    # Define the path form
    pathform = '/Volumes/LabShare2/CMIP6/'+gcm+'/'+expid+'/'+varID+'/**/'+varID+'_*_'+gcm+'_'+expid+'*.nc'
    # Grab the file using glob
    pathID = glob.glob(pathform)[0]
    # print(pathID)

    # Load in the data

    # Return the data
    return pathID

def kept_attributes_list():
    '''Returns a list of attributes for keeping in a dataset'''
    keep_attrs_list = ['history',
                       'authors',
                       ]
    return keep_attrs_list

def list_of_gcms():
    '''This function simply stores a list of the gcms being used in this analysis to have a clean and separate location
    for the storage of the model IDs'''
    mod_id_list= ['ACCESS-CM2',
                'ACCESS-ESM1-5',
                'CESM2',
                'CMCC-CM2-SR5',
                'CMCC-ESM2',
                'EC-Earth3',
                'GFDL-ESM4',
                'INM-CM4-8',
                'INM-CM5-0',
                'IPSL-CM6A-LR',
                'MIROC6',
                'MPI-ESM1-2-HR',
                'MPI-ESM1-2-LR',
                'MRI-ESM2-0',
                'UKESM1-0-LL',
                'HadGEM3-GC31-LL',
                ]

    return mod_id_list

def refine_attributes(ds,
                      parent_temp_path,
                      varID):
    # Keeps only the selected attributes
    ds = pm.keep_list_of_attributes(ds,
                                 kept_attributes_list())

    # Gets the parent attributes and adds them in
    ds = pm.add_attributes(ds,
                        {'parent_attrs':pm.get_parent_attributes(parent_temp_path)})

    # Adds the fill and missing values
    ds = pm.add_fill_miss_values(ds,
                                  varID)

    return ds

if __name__=='__main__':



    # variable_dict contains the names of the different variables used in each netCDF file
    variable_dict = {'temp':'tas',
                         'rnet':'SRFRAD',
                         'wind':'sfcWind',
                         'rh':'hurs',
                         'sh':'hfss',
                         'lh':'hfls',
                         'ps':'ps'}

    # Template for naming files used in calculation of PET
    directory_template = '/Volumes/LabShare2/CMIP6/{gcm}/{expid}/{varid}/{variant}/'
    filename_template = '{varid}_{freq}_{gcm}_{expid}_{variant}_{grid}_{date_range}.nc'

    # Define vegetation type
    veg_type='grass'

    # Define the frequency of data


    # Define the grid dictionary (different gcms use different grids unfortunately :( )
    grid_dict = {'ACCESS-CM2':'gn',
                'ACCESS-ESM1-5':'gn',
                'BCC-CSM2-MR':'gn',
                'CAS-ESM2-0':'gn',
                 'CESM2':'gn',
                'CESM2-WACCM':'gn',
                'CMCC-CM2-SR5':'gn',
                'CMCC-ESM2':'gn',
                'CanESM5':'gn',
                'EC-Earth3':'gr',
                'EC-Earth3-Veg':'gr',
                'EC-Earth3-Veg-LR':'gr',
                'FGOALS-f3-L':'gr',
                'FGOALS-g3':'gn',
                'GFDL-ESM4':'gr1',
                'INM-CM4-8':'gr1',
                'INM-CM5-0':'gr1',
                'IPSL-CM6A-LR':'gr',
                'KACE-1-0-G':'gr',
                'MIROC6':'gn',
                'MPI-ESM1-2-HR':'gn',
                'MPI-ESM1-2-LR':'gn',
                'MRI-ESM2-0':'gn',
                'TaiESM1':'gn',
                'CNRM-ESM2-1':'gr',
                'CNRM-CM6-1':'gr',
                'UKESM1-0-LL':'gn',
                'HadGEM3-GC31-LL':'gn'}

    # Define the date range dictionary for different experiment IDs
    date_range_dict = {'historical':'185001-201412',
                       'ssp119':'201501-210012',
                       'ssp126':'201501-210012',
                       'ssp245':'201501-210012',
                       'ssp370':'201501-210012',
                       'ssp585':'201501-210012'}

    # Define the frequency label (from CMIP6 documentation)
    freq = 'Amon'

    # Loop through the gcms in the list
    gcm_list = list_of_gcms()
    for gcm in gcm_list:
        # Loop through experiment IDs
        expid_list = ['historical',
                      'ssp126',
                      'ssp245',
                      'ssp370',
                      'ssp585']
        for expid in expid_list:

            try:
                path_id = get_monthly_gcm_path(gcm,
                                               expid,
                                               'petgrs')
            except Exception as e:
                print(e)
                print(gcm+' '+expid+' has no petgrs.')
            print(gcm,expid)

            # Define the variants based on GCM
            if gcm in ['CNRM-ESM2-1',
                        'UKESM1-0-LL',
                        'CNRM-CM6-1']:
                variant = 'r1i1p1f2'
            elif gcm in ['HadGEM3-GC31-LL']:
                variant = 'r1i1p1f3'
            elif gcm in ['CESM2']:
                variant = 'r10i1p1f1'
            else:
                variant = 'r1i1p1f1'
            # Define output directory
            output_dir_petgrs = directory_template.format(gcm=gcm,expid=expid,varid='petgrs',variant=variant)

            # Check to determine if the output path exists or not
            if os.path.exists(output_dir_petgrs):
                pass
            else:
                os.makedirs(output_dir_petgrs)

            # Define output directory
            output_dir_vpd = directory_template.format(gcm=gcm,expid=expid,varid='vpd',variant=variant)

            # Check to determine if the output path exists or not
            if os.path.exists(output_dir_vpd):
                pass
            else:
                os.makedirs(output_dir_vpd)

           # Define output directory
            output_dir_petgrs_ad = directory_template.format(gcm=gcm,expid=expid,varid='petgrsad',variant=variant)

            # Check to determine if the output path exists or not
            if os.path.exists(output_dir_petgrs_ad):
                pass
            else:
                os.makedirs(output_dir_petgrs_ad)

           # Define output directory
            output_dir_petgrs_en = directory_template.format(gcm=gcm,expid=expid,varid='petgrsen',variant=variant)

            # Check to determine if the output path exists or not
            if os.path.exists(output_dir_petgrs_en):
                pass
            else:
                os.makedirs(output_dir_petgrs_en)

           # Define output directory
            output_dir_petpt = directory_template.format(gcm=gcm,expid=expid,varid='petpt',variant=variant)

            # Check to determine if the output path exists or not
            if os.path.exists(output_dir_petpt):
                pass
            else:
                os.makedirs(output_dir_petpt)
            try:
                # Input directory path dictionary
                path_dict_in = {    'rh_path'     : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'hurs'),
                                   'wind_path'   : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'sfcWind'),
                                    'lh_path'     : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'hfls'),
                                    'sh_path'     : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'hfss'),
                                    'temp_path'   : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'tas'),
                                    'rnet_path'   : '',

                                    'ps_path'     : get_monthly_gcm_path(gcm,
                                                                         expid,
                                                                         'ps')}

                path_dict_out = {   'pet_alf_path': output_dir_petgrs+filename_template.format(varid='petalf',
                                                                                     freq=freq,
                                                                                     gcm=gcm,
                                                                                     expid=expid,
                                                                                     variant=variant,
                                                                                     grid=grid_dict[gcm],
                                                                                     date_range=date_range_dict[expid]),
                                    'pet_grass_path': output_dir_petgrs+filename_template.format(varid='petgrs',
                                                                                     freq=freq,
                                                                                     gcm=gcm,
                                                                                     expid=expid,
                                                                                     variant=variant,
                                                                                     grid=grid_dict[gcm],
                                                                                     date_range=date_range_dict[expid]),
                                    'vpd_path': output_dir_vpd+filename_template.format(varid='vpd',
                                                                             freq=freq,
                                                                             gcm=gcm,
                                                                             expid=expid,
                                                                             variant=variant,
                                                                             grid=grid_dict[gcm],
                                                                             date_range=date_range_dict[expid]),
                                    'pet_grass_en_path': output_dir_petgrs_en+filename_template.format(varid='petgrsen',
                                                                             freq=freq,
                                                                             gcm=gcm,
                                                                             expid=expid,
                                                                             variant=variant,
                                                                             grid=grid_dict[gcm],
                                                                             date_range=date_range_dict[expid]),
                                    'pet_grass_ad_path': output_dir_petgrs_ad+filename_template.format(varid='petgrsad',
                                                                             freq=freq,
                                                                             gcm=gcm,
                                                                             expid=expid,
                                                                             variant=variant,
                                                                             grid=grid_dict[gcm],
                                                                             date_range=date_range_dict[expid]),
                                    'pet_pt_path': output_dir_petpt+filename_template.format(varid='petpt',
                                                                             freq=freq,
                                                                             gcm=gcm,
                                                                             expid=expid,
                                                                             variant=variant,
                                                                             grid=grid_dict[gcm],
                                                                             date_range=date_range_dict[expid])}


                if os.path.exists(path_dict_out['pet_grass_path']) is False:     # test for existence of RH file
                    # Calculate PET and VPD
                    PET,VPD,PETen,PETad=pm.PenMon(path_dict=path_dict_in,variable_dict=variable_dict,veg_type='grass',rnet_calculated = True,wind_meas_height=10.0,const_wind=1.0)

                    # Add some global attributes for history and contacts
                    PET.attrs['history']='Created from data at CU Boulder using  CU/PSL PMPET.py on ' + datetime.date.today().strftime("%B %d, %Y")
                    PET.attrs['authors']='Nels Bjarke (CU), Joe Barsugli(CU NOAA PSL)'

                    VPD.attrs['history']='Created from data at CU Boulder using  CU/PSL PMPET.py on ' + datetime.date.today().strftime("%B %d, %Y")
                    VPD.attrs['authors']='Nels Bjarke (CU), Joe Barsugli(CU NOAA PSL)'

                    PETad.attrs['history']='Created from data at CU Boulder using  CU/PSL PMPET.py on ' + datetime.date.today().strftime("%B %d, %Y")
                    PETad.attrs['authors']='Nels Bjarke (CU), Joe Barsugli(CU NOAA PSL)'

                    PETen.attrs['history']='Created from data at CU Boulder using  CU/PSL PMPET.py on ' + datetime.date.today().strftime("%B %d, %Y")
                    PETen.attrs['authors']='Nels Bjarke (CU), Joe Barsugli(CU NOAA PSL)'

                    # Change the attributes to the parent attributes of the netcdf files for input and delete the
                    # redundant attributes of the dataset

                    # Create a dictionary for renaming the variable within the dataset for advective and energetic terms
                    var_rename_dict = {'en':['PET','PETEN'],
                                       'ad':['PET','PETAD']}



                    # Rename the PET variable ID for the components
                    PETad = PETad.rename({'PET':'PETAD'})
                    PETen = PETad.rename({'PET':'PETEN'})

                    # Refine the variables
                    PET = refine_attributes(PET,
                                          parent_temp_path = path_dict_in['temp_path'],
                                          varID = 'PET')
                    PETad = refine_attributes(PETad,
                                          parent_temp_path = path_dict_in['temp_path'],
                                          varID = 'PETAD')
                    PETen = refine_attributes(PETen,
                                          parent_temp_path = path_dict_in['temp_path'],
                                          varID = 'PETEN')
                    VPD = refine_attributes(VPD,
                                          parent_temp_path = path_dict_in['temp_path'],
                                          varID = 'VPD')

                    # Write to netCDF files
                    PET.to_netcdf(path_dict_out['pet_grass_path'])
                    PETad.to_netcdf(path_dict_out['pet_grass_ad_path'])
                    PETen.to_netcdf(path_dict_out['pet_grass_en_path'])
                    VPD.to_netcdf(path_dict_out['vpd_path'])
            except Exception as e: print(e)

            try:
                if os.path.exists(path_dict_out['pet_pt_path']) is False:     # test for existence of RH file
                    # Calculate PET and VPD
                    PET=pm.PriestleyTaylor(path_dict=path_dict_in,variable_dict=variable_dict,rnet_calculated = True,wind_meas_height=10.0,const_wind=1.0,alpha=1.26)

                    # Add some global attributes for history and contacts
                    PET.attrs['history']='Created from data at CU Boulder using  CU/PSL PMPET.py on ' + datetime.date.today().strftime("%B %d, %Y")
                    PET.attrs['authors']='Nels Bjarke (CU), Joe Barsugli(CU NOAA PSL)'


                    # Refine the variables
                    PET = refine_attributes(PET,
                                          parent_temp_path = path_dict_in['temp_path'],
                                          varID = 'PET')
                    # Write to netCDF files
                    PET.to_netcdf(path_dict_out['pet_pt_path'])
                    print('Successfully calculated Priestley-Taylor PET for '+gcm +' - '+expid)

                            # print(f"Done with PET computation for ensemble member {path_dict_out["pet_path"]}.")
            except Exception as e: print(e)

