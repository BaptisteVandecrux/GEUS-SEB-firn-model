# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import numpy as np
import lib_subsurface as sub
import lib_initialization as ini
import xarray as xr
import lib_io as io
import os
from lib_initialization import ImportConst, load_json
from progressbar import progressbar

def run_GEUS_model(site, filename):
    # Read paths for weather input and output file
    parameters = load_json()
    output_path = str(parameters['output_path'])
    weather_station = str(parameters['weather_station'])
    station_info = str(parameters['constants']['station_info']['path'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(weather_station)

    #c = set_constants()

    c = ImportConst()
    c.station = site
    c.rh2oice = c.rho_water / c.rho_ice
    c.zdtime = 3600
    NumLayer = 100
    c.num_lay = NumLayer
    c.z_max = 50
    c.dz_ice = c.z_max / NumLayer
    c.verbose = 1
    #c.OutputFolder = "C://Data_save/Output firn model"
    #c.OutputFolder = output_path
    c.OutputFolder = r'C:\Users\brink\Documents\Exjobb\GEUS-SEB-firn-model\Output main_firn'

    df_info = pd.read_csv("Input/Constants/StationInfo.csv", sep=";")
    #df_info = pd.read_csv(station_info, sep=";")

    c.Tdeep = (
        df_info.loc[
            df_info["station name"] == c.station, "deep firn temperature (degC)"
        ].values[0]
        + 273.15
    )

    # c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
    c.lim_new_lay = 0.02
    c.rho_fresh_snow = 315

    #df_aws = pd.read_csv('./Input/'+c.station+'_high-res_meteo.csv', sep=';')
   
    # Original code:
    #df_aws = pd.read_csv(filename, sep=";")
    #df_aws.time = pd.to_datetime(df_aws.time)
    #df_aws = df_aws.set_index("time").resample("H").nearest()

    # New code:
    # DataFrame with the weather data is created, copied from main seb firn
    df_aws = io.load_promice(weather_data_input_path)[:3000]
    df_aws = df_aws.set_index("time").resample("H").mean()

    print("Number of NaN: ")
    print(df_aws.isnull().sum())

    time = df_aws.index.values

    (
        rhofirn,
        snowc,
        snic,
        slwc,
        dgrain,
        tsoil,
        grndc,
        grndd,
        compaction,
        zrfrz,
        zsupimp,
        ts,
        zrogl,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
    ) = ini.IniVar(time, c)

    # Original code
    # ts = df_aws.Tsurf_K.values
    # ts_out = ts.copy() * np.nan
    # net_accum = df_aws.acc_subl_mmweq.values / 1000
    # melt = df_aws.melt_mmweq.values / 1000

    # New code
    # Read parameter values from main SEB firn
    path = r'C:\Users\brink\Documents\Exjobb\GEUS-SEB-firn-model\Input\Weather data\param_output_from_SEB_KAN_U.csv'
    df_SEB_output = pd.read_csv(path, index_col=False)
    snowfall = df_SEB_output['snowfall']
    Tsurf = df_SEB_output['Tsurf']
    sublimation_mweq = df_SEB_output['sublimation_mweq']
    melt_mweq = df_SEB_output['melt_mweq']

    Tsurf_out = Tsurf.copy() * np.nan
    net_accum = [snowfall[i] - sublimation_mweq[i] for i in range(min(len(snowfall), len(sublimation_mweq)))]

    for i in progressbar(range(0, len(time))):
        # print(i)
        (
            snowc[:, i],
            snic[:, i],
            slwc[:, i],
            tsoil[:, i],
            zrfrz[:, i],
            rhofirn[:, i],
            zsupimp[:, i],
            dgrain[:, i],
            zrogl[i],
            Tsurf_out[i],
            grndc[:, i],
            grndd[:, i],
            pgrndcapc[i],
            pgrndhflx[i],
            dH_comp[i],
            snowbkt[i],
            compaction[:, i],
        ) = sub.subsurface(
            Tsurf[i].copy(),
            grndc[:, i - 1].copy(),
            grndd[:, i - 1].copy(),
            slwc[:, i - 1].copy(),
            snic[:, i - 1].copy(),
            snowc[:, i - 1].copy(),
            rhofirn[:, i - 1].copy(),
            tsoil[:, i - 1].copy(),
            dgrain[:, i - 1].copy(),
            net_accum[i].copy(),
            0,
            melt_mweq[i].copy(),
            c.Tdeep,
            snowbkt[i - 1].copy(),
            c,
        )

    # writing output
    thickness_act = snowc * (c.rho_water / rhofirn) + snic * (c.rho_water / c.rho_ice)
    depth_act = np.cumsum(thickness_act, 0)
    density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)

    c.RunName = c.station + "_" + str(c.num_lay) + "_layers"
    i = 0
    succeeded = 0

    while succeeded == 0:
        try:
            os.mkdir(c.OutputFolder + "/" + c.RunName)
            succeeded = 1
        except:
            if i == 0:
                c.RunName = c.RunName + "_" + str(i)
            else:
                c.RunName = c.RunName[: -len(str(i - 1))] + str(i)
            i = i + 1

    io.write_2d_netcdf(snowc, "snowc", depth_act, time, c)
    io.write_2d_netcdf(snic, "snic", depth_act, time, c)
    io.write_2d_netcdf(slwc, "slwc", depth_act, time, c)
    io.write_2d_netcdf(density_bulk, "density_bulk", depth_act, time, c)
    io.write_2d_netcdf(rhofirn, "rhofirn", depth_act, time, c)
    io.write_2d_netcdf(tsoil, "T_ice", depth_act, time, c)
    # io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, time, RunName)
    # io.write_2d_netcdf(rfrz, 'rfrz', depth_act, time, RunName)
    # io.write_2d_netcdf(dgrain, 'dgrain', depth_act, time, RunName)
    io.write_2d_netcdf(compaction, "compaction", depth_act, time, c)
    print(c.RunName)
    return c.RunName


if __name__ == "__main__":
    site_list = [
        "KAN_U",
        "KAN_M"
        #"FA",
    ]  # ['EGP', 'DYE-2','CP1','Saddle', 'Summit', 'KAN_U', 'NASA-SE']
    # other sites to add:   'FA'

    # parameters = load_json()
    # weather_station = str(parameters['weather_station'])
    # weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    # weather_data_input_path = weather_data_input_path_unformatted.format(weather_station)
   
    for site in site_list:
        print(site)
        #filename = "Firn viscosity/Input files/" + site + ".csv" #Old file
        filename = "./Input/Weather data/data_" + site + "_combined_hour.txt"
        #filename = weather_data_input_path
        run_name = run_GEUS_model(site, filename)
        print(run_name)
    



# %% debbuging
# i = 0
# (pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep, psnowbkt, c) = (ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
#                             slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(),
#                             rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
#                             net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)
