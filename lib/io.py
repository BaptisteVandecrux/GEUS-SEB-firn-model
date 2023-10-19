# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import pytz

units = {
    "snowc": "m water equivalent",
    "snic": "m water equivalent",
    "slwc": "m water equivalent",
    "rhofirn": "kg m^-3",
    "density_bulk": "kg m^-3",
    "T_ice": "K",
    "compaction": "m per time step",
    "rfrz": "m water equivalent per time step",
    "dgrain": "mm",
}
long_name = {
    "snowc": "Layer snow content",
    "snic": "Layer ice content",
    "slwc": "Layer liquid water content",
    "rhofirn": "Density of snow only",
    "density_bulk": "Bulk density",
    "T_ice": "Subsurface temperature",
    "compaction": "Layer compaction",
    "rfrz": "Amount of water refrozen",
    "dgrain": "Snow grain diameter",
}


def load_promice_old(path_promice):
    """
    Loading PROMICE data for a given path into a DataFrame.
    + adding time index
    + calculating albedo
    + (optional) calculate RH with regard to water

    INTPUTS:
        path_promice: Path to the desired file containing PROMICE data [string]

    OUTPUTS:
        df: Dataframe containing PROMICE data for the desired settings [DataFrame]
    """

    df = pd.read_csv(path_promice, delim_whitespace=True)
    df["time"] = df.Year * np.nan

    try:
        df["time"] = [
            datetime.datetime(y, m, d, h).replace(tzinfo=pytz.UTC)
            for y, m, d, h in zip(
                df["Year"].values,
                df["MonthOfYear"].values,
                df["DayOfMonth"].values,
                df["HourOfDay(UTC)"].values,
            )
        ]
    except:
        df["time"] = pd.to_datetime(
            df["Year"] * 100000 + df["DayOfYear"] * 100 + df["HourOfDayUTC"],
            format="%Y%j%H",
        )

    df.set_index("time", inplace=True, drop=False)

    # set invalid values (-999) to nan
    df[df == -999.0] = np.nan
    # df['Albedo'] = df['ShortwaveRadiationUp(W/m2)'] / df['ShortwaveRadiationDown(W/m2)']
    df.loc[df["Albedo"] > 1, "Albedo"] = np.nan
    df.loc[df["Albedo"] < 0, "Albedo"] = np.nan
    # df['SnowHeight(m)'] = 2.6 - df['HeightSensorBoom(m)']
    # df['SurfaceHeight(m)'] = 1 - df['HeightStakes(m)']

    # df['RelativeHumidity_w'] = RH_ice2water(df['RelativeHumidity(%)'] ,
    #                                                    df['AirTemperature(C)'])
    
    df = df.loc[df.AirPressurehPa.first_valid_index():,:]
    df = df.loc[df.AirTemperature1C.first_valid_index():,:]

    df = df.set_index("time").resample("H").mean()
    df = df.interpolate()
    return df


def load_CARRA_data(weather_station: str):
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


def load_surface_input_data(surface_input_path, driver='AWS_old'):
    if driver == 'AWS_old':
        return load_promice_old(surface_input_path)
    if driver == 'AWS':
        return load_promice(surface_input_path)
    if driver == 'CARRA':
        return load_CARRA_data(surface_input_path)
    
    print('Driver', driver, 'not recognized')
    return None
    
    
def write_2d_netcdf(data, name_var, depth_act, time, c):
    levels = np.arange(data.shape[0])
    time_days_since = (
        pd.to_timedelta(time - np.datetime64("1900-01-01", "ns")).total_seconds().values
        / 3600
        / 24
    )

    foo = xr.DataArray(
        data, coords=[levels, time_days_since], dims=["level", "time"], name=name_var
    )
    foo.attrs["units"] = units[name_var]
    foo.attrs["long_name"] = long_name[name_var]
    foo.time.attrs["units"] = "days since 1900-01-01"
    foo.level.attrs["units"] = "index of layer (0=top)"

    depth = xr.DataArray(
        depth_act,
        coords=[levels, time_days_since],
        dims=["level", "time"],
        name="depth",
    )
    depth.attrs["units"] = "m below the surface"
    depth.attrs["long_name"] = "Depth of layer bottom"
    depth.time.attrs["units"] = "days since 1900-01-01"
    depth.level.attrs["units"] = "index of layer (0=top)"

    ds = xr.merge([foo, depth])
    ds.to_netcdf(
        c.output_path + "/" + c.RunName + "/" + c.station + "_" + name_var + ".nc"
    )


def write_1d_netcdf(data, c, var_list=None, time=None, name_file="surface"):
    var_info = pd.read_csv("Input/Constants/output_variables_info.csv").set_index(
        "short_name"
    )
    if not time:
        time = data.index
    if not var_list:
        var_list = data.columns
    time_days_since = (
        pd.to_timedelta(time - np.datetime64("1900-01-01", "ns")).total_seconds().values
        / 3600
        / 24
    )

    for name_var in var_list:
        foo = xr.DataArray(
            data[name_var].values,
            coords=[time_days_since],
            dims=["time"],
            name=name_var,
        )
        foo.attrs["units"] = var_info.loc[name_var, "units"]
        foo.attrs["long_name"] = var_info.loc[name_var, "long_name"]
        foo.time.attrs["units"] = "days since 1900-01-01"
        if name_var == var_list[0]:
            ds = foo
        else:
            ds = xr.merge([ds, foo])
    ds.to_netcdf(
        c.output_path + "/" + c.RunName + "/" + c.station + "_" + name_file + ".nc"
    )
