# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
#import __init__
import lib_io as io
import lib_plot as lpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lib_initialization import ImportConst, load_json
from lib_seb_smb_model import HHsubsurf
from lib_CARRA_initialization import load_CARRA_data_opt
from os import mkdir

def run_SEB_firn():
    # Read paths for weather input and output file
    parameters = load_json()
    output_path = str(parameters['output_path'])
    weather_station = str(parameters['weather_station'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(weather_station)

    # Choose data source, CARRA or AWS
    data_source = 'CARRA'
    #data_source = 'AWS'

    # Create struct c with all constant values
    c = set_constants(weather_station)

    if data_source == 'CARRA':
        # Data with weather data is created, from CARRA data
        weather_df = load_CARRA_data_opt(weather_station)[55200:58199]   
        print(weather_df.index[0])
        print(weather_df.index[-1])

    if data_source == 'AWS':
        # DataFrame with the weather data is created, from AWS data
        weather_df = io.load_promice(weather_data_input_path)[:5999] #[3858:(3858 + 8999)]
        weather_df = weather_df.set_index("time").resample("H").mean()
        weather_df = weather_df.interpolate()
        print(weather_df.index[0])
        print(weather_df.index[-1])

    # DataFrame for the surface is created, indexed with time from df_aws
    df_surface = pd.DataFrame()
    df_surface["time"] = weather_df.index
    df_surface = df_surface.set_index("time")

    # The surface values are received 
    (
        df_surface["L"],
        df_surface["LHF"],
        df_surface["SHF"],
        df_surface["theta_2m"],
        df_surface["q_2m"],
        df_surface["ws_10m"],
        df_surface["Re"],
        df_surface["melt_mweq"],
        df_surface["sublimation_mweq"],
        df_surface["SRin"],
        df_surface["SRout"],
        df_surface["LRin"],
        df_surface["LRout_mdl"],
        snowc,
        snic,
        slwc,
        T_ice,
        zrfrz,
        rhofirn,
        zsupimp,
        dgrain,
        zrogl,
        Tsurf,
        grndc,
        grndd,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
        compaction
    ) = HHsubsurf(weather_df, c)

    thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

    # Writing output
    c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
    i = 0
    succeeded = 0
    while succeeded == 0:
        try:
            mkdir(output_path + c.RunName)
            succeeded = 1
        except:
            if i == 0:
                c.RunName = c.RunName + "_" + str(i)
            else:
                c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
            i = i + 1
    
    c.OutputFolder = output_path
    # io.write_2d_netcdf(snic, 'snic', depth_act, weather_df.index, c)
    io.write_2d_netcdf(slwc, 'slwc', depth_act, weather_df.index, c)
    #io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, weather_df.index, c)
    io.write_1d_netcdf(df_surface, c)
    io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, weather_df.index, c)
    io.write_2d_netcdf(T_ice, 'T_ice', depth_act, weather_df.index, c)
    # io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, weather_df.index, RunName)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, weather_df.index, RunName)
    # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, weather_df.index, RunName)
    # io.write_2d_netcdf(compaction, 'compaction', depth_act, weather_df.index, RunName)

    # Plot output
    plt.close("all")
    #lpl.plot_summary(weather_df, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
    lpl.plot_summary(df_surface, c, 'SEB_output')
    lpl.plot_var(c.station, c.RunName, "slwc", ylim=(10, -5), zero_surf=False)
    lpl.plot_var(c.station, c.RunName, "T_ice", ylim=(10, -5), zero_surf=False)
    lpl.plot_var(c.station, c.RunName, "density_bulk", ylim=(10, -5), zero_surf=False)
    
    melt_mweq_cum = df_surface["melt_mweq"].cumsum()
    print("Cumulative sum of melt:",melt_mweq_cum[-1], " mweq")


# Constant definition
# All constant values are defined in a set of csv file in the Input folder.
# They can be modiefied there or by giving new values in the "param" variable. 
# The values of the constant given in param will overright the ones extracted 
# from the csv files. The fieldnames in param should be the same as is c.
def set_constants(weather_station):
    c = ImportConst()
    c.station = weather_station
    c.elev = 2000
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    c.ElevGrad = 0.1
    c.z_max = 50
    c.dz_ice = 1
    NumLayer = int(c.z_max / c.dz_ice)
    #c.num_lay = NumLayer
    c.num_lay = 200
    c.verbose = 1
    c.Tdeep = 250.15
    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315
    c.snowthick_ini = 1
    c.dz_ice = 1
    c.z_ice_max = 50
    c.dt_obs = 3600   
    return c 

if __name__ == "__main__":
    run_SEB_firn()


    

