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

    return df_carra


