# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lib.io as io
import lib.plot as lpl
from lib.initialization import ImportConst
from lib.seb_smb_model import HHsubsurf
from os import mkdir
import time

# def run_SEB_firn():

start_time = time.time()
print('start processing')

# importing standard values for constants
c = ImportConst()

c.station = 'KAN_U'

c.surface_input_path = "./input/weather data/data_KAN_U_2009-2019.txt"

c.surface_input_driver = "AWS_old" 
# c.surface_input_driver = "AWS_new" 
# c.surface_input_driver = "CARRA" 

# assigning constants specific to this simulation
c.z_max = 50
c.num_lay = 100
# c.lim_new_lay = c.accum_AWS/c.new_lay_frac;

df_in = io.load_surface_input_data(c.surface_input_path, driver='AWS_old')

df_in = df_in[:17520]
# #%%
# for var in ['AirPressurehPa',
#        'AirPressurehPa_Origin', 'AirTemperature1C', 'AirTemperature1C_Origin',
#        'AirTemperature2C', 'AirTemperature2C_Origin', 'RelativeHumidity1',
#        'RelativeHumidity1_Origin', 'RelativeHumidity2',
#        'RelativeHumidity2_Origin', 'WindSpeed1ms', 'WindSpeed1ms_Origin',
#        'WindSpeed2ms', 'WindSpeed2ms_Origin', 'WindDirection1d',
#        'WindDirection2d', 'ShortwaveRadiationDownWm2',
#        'ShortwaveRadiationDownWm2_Origin', 'ShortwaveRadiationUpWm2',
#        'ShortwaveRadiationUpWm2_Origin', 'Albedo', 'LongwaveRadiationDownWm2',
#        'LongwaveRadiationDownWm2_Origin', 'LongwaveRadiationUpWm2',
#        'LongwaveRadiationUpWm2_Origin', 'HeightSensorBoomm', 'HeightStakesm']:
#     df_in[[var]].plot(marker='o')
# #%%
print('start/end of input file', df_in.index[0], df_in.index[-1])
# DataFrame for the surface is created, indexed with time from df_aws
df_out = pd.DataFrame()
df_out["time"] = df_in.index
df_out = df_out.set_index("time")

print('reading inputs took %0.03f sec'%(time.time() -start_time))
start_time = time.time()
# The surface values are received 
(
    df_out["L"],
    df_out["LHF"],
    df_out["SHF"],
    df_out["theta_2m"],
    df_out["q_2m"],
    df_out["ws_10m"],
    df_out["Re"],
    df_out["melt_mweq"],
    df_out["sublimation_mweq"],
    df_out["SRin"],
    df_out["SRout"],
    df_out["LRin"],
    df_out["LRout_mdl"],
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
) = HHsubsurf(df_in, c)

print('HHsubsurf took %0.03f sec'%(time.time() -start_time))
start_time = time.time()

thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
depth_act = np.cumsum(thickness_act, 0)
density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

# Writing output
c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
i = 0
succeeded = 0
while succeeded == 0:
    try:
        mkdir(c.output_path + c.RunName)
        succeeded = 1
    except:
        if i == 0:
            c.RunName = c.RunName + "_" + str(i)
        else:
            c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
        i = i + 1

# io.write_2d_netcdf(snic, 'snic', depth_act, df_in.index, c)
io.write_2d_netcdf(slwc, 'slwc', depth_act, df_in.index, c)
#io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, df_in.index, c)
io.write_1d_netcdf(df_out, c)
io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, df_in.index, c)
io.write_2d_netcdf(T_ice, 'T_ice', depth_act, df_in.index, c)
# io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(rfrz, 'rfrz', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(dgrain, 'dgrain', depth_act, df_in.index, RunName)
# io.write_2d_netcdf(compaction, 'compaction', depth_act, df_in.index, RunName)

print('writing output files took %0.03f sec'%(time.time() -start_time))
start_time = time.time()

# Plot output
plt.close("all")
#lpl.plot_summary(df_in, c, 'input_summary', var_list = ['RelativeHumidity1','RelativeHumidity2'])
lpl.plot_summary(df_out, c, 'SEB_output')
lpl.plot_var(c.station, c.RunName, "slwc", ylim=(10, -5), zero_surf=False)
lpl.plot_var(c.station, c.RunName, "T_ice", ylim=(10, -5), zero_surf=False)
lpl.plot_var(c.station, c.RunName, "density_bulk", ylim=(10, -5), zero_surf=False)

plt.figure()
plt.plot(df_out["LRout_mdl"], 
         df_in.LongwaveRadiationUpWm2,
         marker='.',ls='None')
# melt_mweq_cum = df_out["melt_mweq"].cumsum()
# print("Cumulative sum of melt:",melt_mweq_cum[-1], " mweq")
print('plotting took %0.03f sec'%(time.time() -start_time))
start_time = time.time()


# if __name__ == "__main__":
#     run_SEB_firn()


    

