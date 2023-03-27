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
from os import mkdir
import time

def run_SEB_firn():
    # time_array = []
    # cpu_start_time = time.process_time()

    # Read paths for weather input and output file
    parameters = load_json()
    output_path = str(parameters['output_path'])
    weather_station = str(parameters['weather_station'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(weather_station)

    # Create struct c with all constant values
    c = set_constants(weather_station)

    # DataFrame with the weather data is created
    df_aws = io.load_promice(weather_data_input_path)[:5999]
    df_aws = df_aws.set_index("time").resample("H").mean()

    # DataFrame for the surface is created, indexed with time from df_aws
    df_surface = pd.DataFrame()
    df_surface["time"] = df_aws.index
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
        compaction,
    ) = HHsubsurf(df_aws, c)

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
    
        # Write results to one csv file
        # data_to_csv = pd.DataFrame(df_surface["L"])
        # data_to_csv = data_to_csv.assign(LHF = df_surface["LHF"])
        # data_to_csv = data_to_csv.assign(SHF = df_surface["SHF"])
        # data_to_csv = data_to_csv.assign(theta_2m = df_surface["theta_2m"])
        # data_to_csv = data_to_csv.assign(q_2m = df_surface["q_2m"])
        # data_to_csv = data_to_csv.assign(ws_10m = df_surface["ws_10m"])
        # data_to_csv = data_to_csv.assign(Re = df_surface["Re"])
        # data_to_csv = data_to_csv.assign(melt_mweq = df_surface["melt_mweq"])
        # data_to_csv = data_to_csv.assign(sublimation_mweq = df_surface["sublimation_mweq"])

        # data_to_csv2 = pd.DataFrame(snowc)
        # data_to_csv2 = data_to_csv2.assign(snic = snic)
        # data_to_csv2 = data_to_csv2.assign(slwc = slwc)
        # data_to_csv2 = data_to_csv2.assign(T_ice = T_ice)
        # data_to_csv2 = data_to_csv2.assign(zrfrz = zrfrz)
        # data_to_csv2 = data_to_csv2.assign(rhofirn = rhofirn)
        # data_to_csv2 = data_to_csv2.assign(zsupimp = zsupimp)
        # data_to_csv2 = data_to_csv2.assign(dgrain = dgrain)
        # data_to_csv2 = data_to_csv2.assign(zrogl = zrogl)
        # data_to_csv2 = data_to_csv2.assign(Tsurf = Tsurf)
        # data_to_csv2 = data_to_csv2.assign(grndc = grndc)
        # data_to_csv2 = data_to_csv2.assign(grndd = grndd)
        # data_to_csv2 = data_to_csv2.assign(pgrndcapc = pgrndcapc)
        # data_to_csv2 = data_to_csv2.assign(pgrndhflx = pgrndhflx)
        # data_to_csv2 = data_to_csv2.assign(dH_comp = dH_comp)
        # data_to_csv2 = data_to_csv2.assign(compaction = compaction)
        # data_to_csv2 = data_to_csv2.assign(snowbkt = snowbkt)

     #   filename = 'runsebfirn_output_opt12.2_200lay_short.csv'
        # filename= '21res__opt_manual.csv'
        # pd.DataFrame(data_to_csv).to_csv(filename)


    c.OutputFolder = output_path
    # io.write_2d_netcdf(snic, 'snic', depth_act, df_aws.index, c)
   # io.write_2d_netcdf(slwc, 'slwc', depth_act, df_aws.index, c)
    #io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_aws.index, c)
   # io.write_1d_netcdf(df_surface, c)
   # io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_aws.index, c)
   # io.write_2d_netcdf(T_ice, 'T_ice', depth_act, df_aws.index, c)
    # io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_aws.index, RunName)
    # io.write_2d_netcdf(compaction, 'compaction', depth_act, df_aws.index, RunName)

    # cpu_end_time = time.process_time()
    # cpu_time = (cpu_end_time - cpu_start_time)
    # time_array.append(cpu_time)
    # print("Time in main seb firn", time_array)

    # Plot output
    # plt.close("all")
    # #lpl.plot_summary(df_aws, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
    # lpl.plot_summary(df_surface, c, 'SEB_output')
    # lpl.plot_var(c.station, c.RunName, "slwc", ylim=(10, -5), zero_surf=False)
    # lpl.plot_var(c.station, c.RunName, "T_ice", ylim=(10, -5), zero_surf=False)
    # lpl.plot_var(c.station, c.RunName, "density_bulk", ylim=(10, -5), zero_surf=False)
    #df_L_LHF_SHF = df_surface[['L','LHF', 'SHF']]
    #lpl.plot_summary(df_L_LHF_SHF, c, 'L, LHF, SHF output', var_list = ['L','LHF', 'SHF'])



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
    c.num_lay = NumLayer
    #c.num_lay = 200
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

if __name__ == '__main__':
    
    run_SEB_firn()




    

