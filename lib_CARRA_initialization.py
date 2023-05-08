import datetime
import netCDF4 as nc
import numpy as np
import pandas as pd

import datetime
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr


def load_CARRA_data_opt(weather_station: str):
    print("- Reading data from CARRA reanalysis set -")

    folder_path = 'C:/Users/brink/Documents/Exjobb/GEUS-SEB-firn-model/Input/KAN_U_M_CARRA_data' 
    
    # finding the index of the station in the station list
    df =  pd.read_csv(folder_path+'/AWS_KAN_U_M_location.csv')
    ind_aws = df.stid[df.stid==weather_station].index.values[0]

    aws_ds = xr.Dataset()
    
    for var in ['t2m', 'al', 'ssrd', 'strd', 'sp', 'skt', 'si10', 'r2', 'tp']:
        print('extracting',var)
        file_path = folder_path + '/CARRA_' + var + '_nearest_PROMICE.nc'
        if var != 'tp':
            # extracting data at AWS
            aws_ds[var] = xr.open_dataset(file_path)[var].isel(x=ind_aws, y=ind_aws)
        else:
            tmp = xr.open_dataset(file_path)[var].isel(x=ind_aws, y=ind_aws)
            # each time step contains t + 18 and t + 6 accumulated precipitation
            # we do tp(t+18) - tp(t+6) and shift the time index by 6h
            # so tp(t) ends up containing total precip from t to t+12
            tmp = tmp.isel(step=1) - tmp.isel(step=0)
            tmp['time'] = tmp.time + pd.Timedelta('6H')
            tmp = tmp / 12  # in kg / m2 / h
            tmp = tmp.resample(time='3H').ffill()  # this average precip is then repeated at t, t+3, t+6, t+9
            tmp = tmp*3  # eventually we multiply by 3 to get the total precip during the 3h time step
            aws_ds[var] = tmp
            # because the first tp value is at 1990-09-01T06:00:00,
            # the first two time steps (T00 and T03) were left as nan by ffill 
            # and need to be filled with zeros
            aws_ds[var] = aws_ds[var].fillna(0)
    
    # unit conversion
    aws_ds['al'] = aws_ds['al']/100
    aws_ds['sp'] = aws_ds['sp']/100

    aws_ds['ssrd'] = aws_ds['ssrd'] / (3 * 3600)


    aws_ds['strd'] = aws_ds['strd'] / (3 * 3600)

    # deriving LRout (surface thermal radiation upward = stru)
    emissivity = 0.97
    aws_ds['LongwaveRadiationUpWm2'] = ((aws_ds['t2m'])**4 * (emissivity) * 5.67e-8) + (1-emissivity)*aws_ds['strd']
    aws_ds['ShortwaveRadiationUpWm2'] = aws_ds['ssrd'] * aws_ds['al'] 

    # unit conversion
    aws_ds['t2m'] = aws_ds['t2m']-273.15

    aws_ds['Snowfallmweq'] = xr.where(aws_ds.t2m >= 0, aws_ds['tp'], 0) / 1000 # conversion to m w.eq. 
    aws_ds['Rainfallmweq'] = xr.where(aws_ds.t2m < 0, 0, aws_ds['tp']) / 1000 # conversion to m w.eq. 

    aws_ds['HeightTemperature1m'] = aws_ds.al *0 + 2
    aws_ds['HeightHumidity1m'] = aws_ds.al *0 + 2
    aws_ds['HeightWindSpeed1m'] = aws_ds.al *0 + 10

    # converting to a pandas dataframe and renaming some of the columns
    df_carra = aws_ds.to_dataframe().rename(columns={
                            't2m': 'AirTemperature1C', 
                            'r2': 'RelativeHumidity1', 
                            'si10': 'WindSpeed1ms', 
                            'sp': 'AirPressurehPa', 
                            'ssrd': 'ShortwaveRadiationDownWm2',
                            'strd': 'LongwaveRadiationDownWm2'          
                        })

    # Fill null values with 0
    df_carra['ShortwaveRadiationDownWm2'] = df_carra['ShortwaveRadiationDownWm2'].fillna(0)
    df_carra['ShortwaveRadiationUpWm2'] = df_carra['ShortwaveRadiationUpWm2'].fillna(0)

    # filename = 'new_df_carra3.csv'
    # df_carra.to_csv(filename)

    return df_carra



# OLD script
# def extract_albedo(folder_path, weather_station):
#     albedo_path = folder_path + '/CARRA_' + 'al' + '_nearest_PROMICE.nc'
#     albedo_ds = nc.Dataset(albedo_path)['al'][:]
    
#     albedo_KAN_M = albedo_ds[:,0][:,0]
#     albedo_KAN_U = albedo_ds[:,1][:,1]

#     # Make albedo range between 0 and 1
#     albedo_KAN_M = [albedo_KAN_M[i]/100 for i in range(0,len(albedo_KAN_M))]
#     albedo_KAN_U = [albedo_KAN_U[i]/100 for i in range(0,len(albedo_KAN_U))]

#     if weather_station == 'KAN_M':
#         return albedo_KAN_M
    
#     if weather_station == 'KAN_U':
#         return albedo_KAN_U

# # Relative humidity    
# def extract_rh(folder_path, weather_station):
#     RH2_path = folder_path + '/CARRA_' + 'r2' + '_nearest_PROMICE.nc'
#     RH2_ds = nc.Dataset(RH2_path)['r2'][:]
#     RH2_KAN_M = RH2_ds[:,0][:,0]
#     RH2_KAN_U = RH2_ds[:,1][:,1]

#     if weather_station == 'KAN_M':
#         return RH2_KAN_M
    
#     if weather_station == 'KAN_U':
#         return RH2_KAN_U

# # Wind speed
# def extract_wind_speed(folder_path, weather_station):
#     WS_path = folder_path + '/CARRA_' + 'si10' + '_nearest_PROMICE.nc'
#     WS_ds = nc.Dataset(WS_path)['si10'][:]
#     WS_KAN_M = WS_ds[:,0][:,0]
#     WS_KAN_U = WS_ds[:,1][:,1]

#     if weather_station == 'KAN_M':
#         return WS_KAN_M
    
#     if weather_station == 'KAN_U':
#         return WS_KAN_U

# # Skin temperature / surface temperature in Kelvin
# def extract_skin_temp(folder_path, weather_station):
#     skt_path = folder_path + '/CARRA_' + 'skt' + '_nearest_PROMICE.nc'
#     skt_ds = nc.Dataset(skt_path)['skt'][:]
#     skt_KAN_M = skt_ds[:,0][:,0]
#     skt_KAN_U = skt_ds[:,1][:,1]

#     if weather_station == 'KAN_M':
#         return skt_KAN_M
    
#     if weather_station == 'KAN_U':
#         return skt_KAN_U
    
# # Surface pressure (converting pressure unit from Pa to hPa)
# def extract_surface_pressure(folder_path, weather_station):
#     sp_path = folder_path + '/CARRA_' + 'sp' + '_nearest_PROMICE.nc'
#     sp_ds = nc.Dataset(sp_path)['sp'][:]
#     sp_KAN_M = sp_ds[:,0][:,0] / 100
#     sp_KAN_U = sp_ds[:,1][:,1] / 100   

#     if weather_station == 'KAN_M':
#         return sp_KAN_M
    
#     if weather_station == 'KAN_U':
#         return sp_KAN_U

# # Temperature, given in Kelvin (transformed to Celsius to match AWS data)
# def extract_temp(folder_path, weather_station):
#     t2m_path = folder_path + '/CARRA_' + 't2m' + '_nearest_PROMICE.nc'
#     t2m_ds = nc.Dataset(t2m_path)['t2m'][:]
#     t2m_KAN_M = t2m_ds[:,0][:,0] - 273.15
#     t2m_KAN_U = t2m_ds[:,1][:,1] - 273.15

#     if weather_station == 'KAN_M':
#         return t2m_KAN_M
    
#     if weather_station == 'KAN_U':
#         return t2m_KAN_U
    
# # Surface long-wave (thermal) radiation downwards, LRin, strd
# def extract_LRin(folder_path, weather_station):
#     LRin_path = folder_path + '/CARRA_' + 'strd' + '_nearest_PROMICE.nc'
#     LRin_ds = nc.Dataset(LRin_path)['strd'][:]

#     LRin_KAN_M = LRin_ds[:,0][:,0]
#     LRin_KAN_U = LRin_ds[:,1][:,1]

#     # Convert unit from Jm-2 to Wm-2
#     step = 3 # 3hourly data
#     for i in range(0, len(LRin_KAN_U)):
#         LRin_KAN_M[i] = LRin_KAN_M[i] / (step * 3600)
#         LRin_KAN_U[i] = LRin_KAN_U[i] / (step * 3600)

#     if weather_station == 'KAN_M':
#         return LRin_KAN_M
    
#     if weather_station == 'KAN_U':
#         return LRin_KAN_U
    
# # Calculate surface long-wave (thermal) radiation upwards, LRout
# def extract_LRout(weather_station, LRin_ds, temp_ds):
#     # Ice sheet surface emissivity:
#     emissivity = 0.97

#     if weather_station == 'KAN_M':
#         LRout_KAN_M = LRin_ds.copy()
#         for i in range(0, len(LRout_KAN_M)-1):
#             LRout_KAN_M[i] = ((temp_ds[i])**4 * (emissivity) * 5.67e-8) + (1-emissivity)*LRin_ds[i]
#         return LRout_KAN_M
    
#     if weather_station == 'KAN_U':
#         LRout_KAN_U = LRin_ds.copy()
#         for i in range(0, len(LRout_KAN_U)-1):
#             LRout_KAN_U[i] = ((temp_ds[i])**4 * (emissivity) * 5.67e-8) + (1-emissivity)*LRin_ds[i]
#         return LRout_KAN_U  

# # Surface short-wave (solar) radiation downwards
# def extract_SRin(folder_path, weather_station):
#     SRin_path = folder_path + '/CARRA_' + 'ssrd' + '_nearest_PROMICE.nc'
#     SRin_ds = nc.Dataset(SRin_path)['ssrd'][:]
#     SRin_KAN_M = SRin_ds[:,0][:,0]
#     SRin_KAN_U = SRin_ds[:,1][:,1]

#     step = 3 # 3hourly data
#     # Convert unit from Jm-2 to Wm-2
#     for i in range(0, len(SRin_KAN_U)-1):
#         SRin_KAN_M[i] = SRin_KAN_M[i] / (step * 3600)
#         SRin_KAN_U[i] = SRin_KAN_U[i] / (step * 3600)

#     if weather_station == 'KAN_M':
#         return SRin_KAN_M
    
#     if weather_station == 'KAN_U':
#         return SRin_KAN_U
    
# # Surface short-wave (solar) radiation upwards
# def extract_SRout(weather_station, SRin_ds, albedo_ds):
#     if weather_station == 'KAN_M':
#         SRout_KAN_M = SRin_ds.copy()
#         for i in range(0, len(SRin_ds)):
#             SRout_KAN_M[i] = SRin_ds[i] * albedo_ds[i]
#         return SRout_KAN_M
    
#     if weather_station == 'KAN_U':
#         SRout_KAN_U = SRin_ds.copy()
#         for i in range(0, len(SRin_ds)):
#             SRout_KAN_U[i] = SRin_ds[i] * albedo_ds[i]
#         return SRout_KAN_U

# # Total precipitation
# def extract_precipitation(folder_path, weather_station):
#     # accumulated from the initial time of the forecast to the forecast time step.
#     # 12 hour time steps
#     tp_path = folder_path + '/CARRA_' + 'tp' + '_nearest_PROMICE.nc'
#     tp_ds = nc.Dataset(tp_path)['tp'][:] 
#     tp_time_ds = nc.Dataset(tp_path)['time'][:]
    
#     # Initializing variables for the difference in acc values (total precipitation)
#     diff_acc_KAN_M_arr = np.empty(len(tp_ds), dtype=object)
#     diff_acc_KAN_U_arr = np.empty(len(tp_ds), dtype=object)

#     for i in range(0, len(tp_ds)):
#         # Calculating the difference which I believe will give the total 
#         # accumulation at the time step
#         diff_acc_KAN_M = tp_ds[i][:,0][:,0][1]-tp_ds[i][:,0][:,0][0]
#         diff_acc_KAN_M_arr[i] = diff_acc_KAN_M
#         diff_acc_KAN_U = tp_ds[i][:,1][:,1][1]-tp_ds[i][:,1][:,1][0]
#         diff_acc_KAN_U_arr[i] = diff_acc_KAN_U

#     # Creating 3hour time step length
#     i_3hour = len(tp_ds)*4
#     tp_KAN_M_3hour = np.empty(i_3hour, dtype=object)
#     tp_KAN_U_3hour = np.empty(i_3hour, dtype=object)
#     tp_time_3hour = np.empty(i_3hour, dtype=object)
#     start_time = tp_time_ds[0]
#     index = 0
#     for i in range(0, len(tp_ds)):
#         # A linear approximation. Dividing the accumulation with 4, 
#         # to get a mean value for each third hour.
#         preci_per_3hour_KAN_M = diff_acc_KAN_M_arr[i] / 4
#         preci_per_3hour_KAN_U = diff_acc_KAN_U_arr[i] / 4

#         # Assigning the mean value to the four time steps within twelve hours 
#         for hour in range(0, 4):
#             tp_time_3hour[index] = start_time + (index * 3 * 3600)
#             # Converting Kg m-2 to m w. eq. (1 Kg m-2 = 1/1000 m w. eq.)
#             tp_KAN_M_3hour[index] = preci_per_3hour_KAN_M / 1000
#             tp_KAN_U_3hour[index] = preci_per_3hour_KAN_U / 1000
#             hour = hour + 1
#             index = index + 1
#             if hour == 4:
#                 hour = 0

#     if weather_station == 'KAN_M':
#         return tp_KAN_M_3hour
    
#     if weather_station == 'KAN_U':
#         return tp_KAN_U_3hour

# # Derive snowfall and rainfall from total precipitation
# def extract_snowfall_rainfall(tp_ds, temp_ds):
#     snowfall = tp_ds[:].copy()
#     rainfall = tp_ds[:].copy()

#     # If temperature > 273.15 K, precipitation is rain, otherwise it is snow
#     for i in range(0, len(temp_ds)-1):
#             if temp_ds[i] > 273.15:
#                 snowfall[i] = 0
#             else:
#                 rainfall[i] = 0

#     return snowfall, rainfall 
    
# #Heights for temp, humidity, wind stations
# def extract_heights(temp_ds, ws_ds):
#     # Height temperature
#     z_T = temp_ds.copy()
#     for i in range(0, len(temp_ds)):
#         z_T[i] = 2

#     # Height humidity
#     z_RH = z_T.copy()

#     # Height wind speed
#     z_WS = ws_ds.copy()
#     for i in range(0, len(ws_ds)):
#         z_WS[i] = 10

#     return z_T, z_RH, z_WS


# def load_CARRA_data(weather_station: str):
#     print("- Reading data from CARRA reanalysis set -")

#     folder_path = 'C:/Users/brink/Documents/Exjobb/GEUS-SEB-firn-model/Input/KAN_U_M_CARRA_data' 
       
#     # Time in seconds
#     time_ds = nc.Dataset(folder_path + '/CARRA_' + 'al' + '_nearest_PROMICE.nc')['time'][:]

#     # Extracting variables from CARRA data
#     albedo_ds = extract_albedo(folder_path, weather_station)
#     rh_ds = extract_rh(folder_path, weather_station)
#     ws_ds = extract_wind_speed(folder_path, weather_station)
#     skt_ds = extract_skin_temp(folder_path, weather_station)
#     sp_ds = extract_surface_pressure(folder_path, weather_station)
#     temp_ds = extract_temp(folder_path, weather_station)
#     LRin_ds = extract_LRin(folder_path, weather_station)
#     LRout_ds = extract_LRout(weather_station, LRin_ds, temp_ds)
#     SRin_ds = extract_SRin(folder_path, weather_station)  
#     SRout_ds = extract_SRout(weather_station, SRin_ds, albedo_ds)  
#     tp_ds = extract_precipitation(folder_path, weather_station)
#     snowfall_ds, rainfall_ds = extract_snowfall_rainfall(tp_ds, temp_ds)
#     z_T, z_RH, z_WS = extract_heights(temp_ds, ws_ds)

#     # Convert seconds from 1970 to a datetime object
#     startdate_time_1970 = datetime.datetime(1970, 1, 1, 0, 0, 0)
#     time_formatted = np.empty((len(time_ds),), dtype = object)

#     for i in range(0, len(time_ds)):
#         time_delta = datetime.timedelta(seconds=int(time_ds[i]))
#         datetime_obj = startdate_time_1970 + time_delta
#         # Format datetime object as string
#         formatted_time = datetime_obj.strftime('%Y-%m-%dT%H:%M:%S.%f')
#         time_formatted[i] = formatted_time
    
#     # If data is of different lengths
#     min_length = min(len(time_formatted), len(temp_ds), len(rh_ds), len(ws_ds), len(sp_ds), len(SRin_ds), len(LRin_ds), len(snowfall_ds))

#     df_carra = pd.DataFrame({
#                             'time': time_formatted[:min_length],
#                             'AirTemperature1C': temp_ds[:min_length], 
#                             'HeightTemperature1m': z_T[:min_length], 
#                             'RelativeHumidity1': rh_ds[:min_length], 
#                             'HeightHumidity1m': z_RH[:min_length],
#                             'WindSpeed1ms': ws_ds[:min_length], 
#                             'HeightWindSpeed1m':z_WS[:min_length], 
#                             'AirPressurehPa': sp_ds[:min_length], 
#                             'ShortwaveRadiationDownWm2': SRin_ds[:min_length],
#                             'ShortwaveRadiationUpWm2': SRout_ds[:min_length], 
#                             'LongwaveRadiationDownWm2': LRin_ds[:min_length], 
#                             'LongwaveRadiationUpWm2': LRout_ds[:min_length], 
#                             'Snowfallmweq': snowfall_ds[:min_length], 
#                             'Rainfallmweq': rainfall_ds[:min_length]             
#                         })

#     return df_carra
   
