# BUILD_DLR.PY
# ----------------------------------------------------------------------------------------------------------------------------------
# script with functions to calculate dynamic line ratings using ieee standard IEEE 738-2023
# called in prepare_network
# Created by Megan Jones (mchjones@umich.edu) 10/03/2025
# ----------------------------------------------------------------------------------------------------------------------------------

import xarray as xr
import pandas as pd
import numpy as np
import time

import copy
import hashlib
import logging
import re
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import requests
import yaml
from snakemake.utils import update_config
import logging
import math

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------------------------------------------------------------------
# line info
# ----------------------------------------------------------------------------------------------------------------------------------

def get_elevation_open(lat, lon):
    url = f"https://api.open-elevation.com/api/v1/lookup?locations={lat},{lon}"
    response = requests.get(url)
    data = response.json()
    return data['results'][0]['elevation']

def get_bearing(lat0,lon0,lat1,lon1): # based on: https://www.movable-type.co.uk/scripts/latlong.html?from=49.243824,-121.887340&to=49.227648,-121.89631
    lat0 = math.radians(lat0)
    lon0 = math.radians(lon0)
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    bearing = math.atan2(math.cos(lat0) * math.sin(lat1) - math.sin(lat0) * math.cos(lat1) * math.cos(lon1 - lon0),
                         math.sin(lon1 - lon0) * math.cos(lat1))
    bearing_deg = (math.degrees(bearing) + 360) % 360 # this gives angle 0-360 from east as 0 for some reason?
    az = (90 - bearing_deg) % 360
    return az

def get_elevation_ref(lat,lon,ref):
    ref["dist"] = ((ref["lat"] - lat)**2 + (ref["lon"] - lon)**2)**0.5
    min_index = ref["dist"].idxmin()
    return ref.loc[min_index,"elevation"]

def get_line_info(network,elevation_ref):
    line_info = pd.DataFrame(index=network.lines.index.values,columns=['bus0','lat0','lon0','elevation0','bus1','lat1','lon1','elevation1','azimuth'])
    line_info['bus0'] = network.lines['bus0'] 
    line_info['bus1'] = network.lines['bus1']
    for idx, row in line_info.iterrows():
        bus0 = row['bus0']
        lat0 = network.buses.loc[bus0,'y']
        lon0 = network.buses.loc[bus0,'x']
        line_info.loc[idx,'lat0'] = lat0
        line_info.loc[idx,'lon0'] = lon0
        line_info.loc[idx,'elevation0'] =  get_elevation_ref(lat0,lon0,elevation_ref)

        bus1 = row['bus1']
        lat1 = network.buses.loc[bus1,'y']
        lon1 = network.buses.loc[bus1,'x']
        line_info.loc[idx,'lat1'] = lat1
        line_info.loc[idx,'lon1'] = lon1
        line_info.loc[idx,'elevation1'] =  get_elevation_ref(lat1,lon1,elevation_ref)

        az = get_bearing(lat0,lon0,lat1,lon1)
        line_info.loc[idx,'azimuth'] = az
    return line_info

# ----------------------------------------------------------------------------------------------------------------------------------
# dumb trig functions bc they give the standard in degrees
# ----------------------------------------------------------------------------------------------------------------------------------

def sind(x):
    return np.sin(np.deg2rad(x))
def cosd(x):
    return np.cos(np.deg2rad(x))
def tand(x):
    return np.tan(np.deg2rad(x))
def asind(x):
    return np.degrees(np.arcsin(x))
def acosd(x):
    return np.degrees(np.arccos(x))
def atand(x):
    return np.degrees(np.arctan(x))

# ----------------------------------------------------------------------------------------------------------------------------------
# heat balance equation functions
# ----------------------------------------------------------------------------------------------------------------------------------

def calc_I(q_c,q_r,q_s,R):
    return np.sqrt((q_c + q_r - q_s)/R)

def calc_R(R_high,R_low,T_high,T_low,T_avg):
    return ((R_high - R_low)/(T_high - T_low))*(T_avg - T_low) + R_low

def calc_reynolds(v_w, T_s, T_a, elevation, D):
    T_film = (T_s + T_a)/2
    mu = (1.458 * 10**(-6) * (T_film + 273)**1.5)/(T_film + 383.4)
    rho = (1.293 - 1.525* 10**(-4) * elevation + 6.379 * 10**(-9) * elevation**2)/(1 + 0.00367 * T_film)
    n = (D * rho * v_w)/mu
    return n, T_film, rho

def calc_Ks(T_film,phi):
    k = 2.424*10**(-2) + 7.477*10**(-5)*T_film - 4.407*10**(-9)*T_film**2
    K_ang = 1.194 - cosd(phi) + 0.194 * cosd(2*phi) + 0.368 * sind(2*phi)
    return K_ang, k

def calc_q_c(v_w, phi, T_s, T_a, elevation, D_0):
    n, T_film, rho = calc_reynolds(v_w, T_s, T_a, elevation, D_0)
    K_ang, k = calc_Ks(T_film,phi)
    if (T_s - T_a) >= 0:
        q_cn = 3.645 * rho**0.5 * D_0**0.75 * (T_s - T_a)**1.25
    else:
        q_cn = 0
    q_c1 = K_ang * (1.01 + 1.35*n**0.52) * k * (T_s - T_a)
    q_c2 = K_ang * 0.754 * n**0.6 * k * (T_s - T_a)
    options = [q_cn,q_c1,q_c2]
    return np.max(options)

def calc_q_r(T_s, T_a, D_0, em):
    return 17.8 * D_0 * em * (((T_s + 273)/100)**4 - ((T_a + 273)/100)**4)

def calc_solar_heating(h,data):
    Q_s = 0
    for i in range(7):
        Q_s = Q_s + data[i]*h**i
    k_solar = 1 + (1.148*10**(-4))*h + (-1.108*10**(-8))*h**2
    return Q_s*k_solar

def calc_q_s(alpha,A,lat,hour,day,Z_l,solar_coeff_data,q_se,D_0):
    omega = (hour - 12)*15
    delta = 23.45*sind(((284+day)/365)*360)
    chi = sind(omega)/(sind(lat)*cosd(omega) - cosd(lat)*tand(delta))
    if omega < 0:
        if chi >= 0:
            c = 0
        elif chi < 0:
            c = 180
    elif omega >= 0:
        if chi >= 0:
            c = 180
        elif chi < 0:
            c = 360
    Z_c = c + atand(chi)
    h = asind(cosd(lat)*cosd(delta)*cosd(omega) + sind(lat)*sind(delta))
    theta = acosd(cosd(h) * cosd(Z_c - Z_l))
    if q_se == None:
        q = calc_solar_heating(h,solar_coeff_data)
    else:
        q = q_se
    return alpha * q * sind(theta) * D_0

# ----------------------------------------------------------------------------------------------------------------------------------
# computation functions
# ----------------------------------------------------------------------------------------------------------------------------------

## functions called in-line (add_extra_components.py)

def trunk(df,cap):
    trunk = df.where(df <= cap,cap)
    return trunk

def get_bus_data(n):
    coord = pd.DataFrame({
        'lon': n.buses['x'],
        'lat': n.buses['y']
    }, index=n.buses.index.values)
    return coord

def get_baseline(static,const,line,solar_coeff):
    inputs = combine_inputs(static,const,line,solar_coeff)
    baseline_rating = calc_dlr(inputs)
    return baseline_rating

def calc_system_baseline(line_info,static,const,solar_coeff):
    baseline = {}
    for key, value in line_info.items():
        #logger.info(f"Calculating baseline rating for line {key}...")
        bus_dlrs = {}
        for bus in ['bus0','bus1']:
            line_data = value[bus]
            bus_dlrs[bus] = get_baseline(static,const,line_data,solar_coeff)
        baseline[str(key)] = min(bus_dlrs['bus0'],bus_dlrs['bus1'])
    logger.info("Baseline calculations complete :)")
    return baseline

def get_relative(dict,base):
    relative = {}
    for key, value in dict.items():
        baseline_rating = base[str(key)]
        relative[key] = list(value/baseline_rating)
    return pd.DataFrame(relative)

# pending
def unpack(input_data):
    R_high = input_data['R_high']
    R_low = input_data['R_low']
    T_high = input_data['T_high']
    T_low = input_data['T_low']
    T_avg = input_data['T_avg']
    v_w = input_data['v_w']
    phi = input_data['phi']
    T_s = input_data['T_s']
    T_a = input_data['T_a']
    elevation = input_data['elevation']
    D_0 = input_data['D_0']
    em = input_data['em']
    alpha = input_data['alpha']
    a = input_data['A']
    lat = input_data['lat']
    hour = input_data['hour']
    day = input_data['day']
    Z_l = input_data['Z_l']
    q_s = input_data['Q_s']
    solar_coeff_data = input_data['solar_coeff_data']
    return R_high, R_low, T_high, T_low, T_avg, v_w, phi, T_s, T_a,elevation, D_0, em, alpha, a, lat, hour, day, Z_l, q_s, solar_coeff_data

def calc_dlr(input_data):
    R_high,R_low,T_high,T_low,T_avg,v_w,phi,T_s,T_a,elevation,D_0,em,alpha,a,lat,hour,day,Z_l,q_s,solar_coeff_data = unpack(input_data)
    r = calc_R(R_high,R_low,T_high,T_low,T_avg)
    q_c = calc_q_c(v_w, phi, T_s, T_a, elevation, D_0)
    q_r = calc_q_r(T_s, T_a, D_0, em)
    q_s = calc_q_s(alpha,a,lat,hour,day,Z_l,solar_coeff_data,q_s,D_0)
    ampacity = calc_I(q_c,q_r,q_s,r)
    return ampacity

def combine_inputs(weather,const,line,solar):
    combined = {**const, **line, **weather}
    combined['solar_coeff_data'] = solar
    return combined


def calc_military(hour):
    hour_of_day = hour % 24
    return hour_of_day

def calc_day(hour):
    day = hour // 24
    return day

def k_to_c(temp):
    return temp - 273.15

def get_v_phi(u,v,Z_l):
    v_w = np.sqrt(u**2 + v**2)
    theta_w = np.degrees(np.arctan2(u, v))  
    phi = abs((theta_w - Z_l + 180) % 360 - 180)
    v_w = np.nan_to_num(v_w, nan=0.0, posinf=0.0, neginf=0.0)
    phi = np.nan_to_num(phi, nan=0.0, posinf=0.0, neginf=0.0)
    return v_w, phi

def get_worse(dlr0,dlr1):
    return [min(x, y) for x, y in zip(dlr0, dlr1)]

def calc_wus_inputs(wus):
    wus_inputs = {
        "T_a": wus['T_a'],
        "v_w": wus['v_w'],
        "Q_s": wus['Q_s'],
        "phi": wus['phi'],
        "hour": wus['hour'], 
        "day": wus['day']
    }
    return wus_inputs

def get_line_dlr_inputs(tamu_data):
    lines = {}
    buses = {}
    for idx, row in tamu_data.iterrows():
         lines[idx] = {}
         lines[idx]['bus0'] = { # sample line data for R, T, D_0, A stays the same - elevation, lat, Z_l changes
            "R_high": 8.688*10**(-5),
            "R_low": 7.283*10**(-5),
            "T_high": 75,
            "T_low": 25,
            "elevation": row['elevation0'], 
            "D_0": 0.02814, 
            "A": 4.028324458*10**(-7), 
            "lat": row['lat0'], 
            "Z_l": row['azimuth'],
        }
         lines[idx]['bus1'] = { # sample line data for R, T, D_0, A stays the same - elevation, lat, Z_l changes
            "R_high": 8.688*10**(-5),
            "R_low": 7.283*10**(-5),
            "T_high": 75,
            "T_low": 25,
            "elevation": row['elevation1'], 
            "D_0": 0.02814, 
            "A": 4.028324458*10**(-7), 
            "lat": row['lat1'], 
            "Z_l": row['azimuth'],
        }
         buses[idx] = {'bus0':row['bus0'],'bus1':row['bus1']}
    bus_list = pd.concat([tamu_data['bus1'], tamu_data['bus0']]).unique().tolist()
    return lines,buses,bus_list

def fixed_wind():
    mag = 0.61
    phi = 90
    return mag, phi

def get_wus_var(data,cluster,azimuth,wind_fix,phi_fix):
    #cluster_data = {
    #    't': data['t'][cluster],
    #    'u': data['u'][cluster],
    #    'v': data['v'][cluster],
    #    'q': data['q'][cluster],
    #} > new method to deal with the xarray approach to extracting climate data
    cluster_data = {
        't': data['T2'].sel(points=cluster).to_series(),
        'u': data['U10'].sel(points=cluster).to_series(),
        'v': data['V10'].sel(points=cluster).to_series(),
        'q': data['SWDNB'].sel(points=cluster).to_series(),
    }
    hours = cluster_data['t'].index
    inter_calc = pd.DataFrame(index=hours)
    
    inter_calc['hour'] = calc_military(hours)
    inter_calc['day'] = calc_day(hours)
    
    if wind_fix:
        inter_calc['v_w'],inter_calc['phi'] = fixed_wind() 
    else:
        inter_calc['v_w'],inter_calc['phi'] = get_v_phi(cluster_data['u'],cluster_data['v'],azimuth)
    
    if phi_fix:
        inter_calc['phi'] = 0
    
    inter_calc['T_a'] = k_to_c(cluster_data['t'])
    inter_calc['Q_s'] = cluster_data['q'] 

    return inter_calc

def get_dlr(inter_calc,line_data,constant,solar_coeff):
    dlrs = []
    for hour, row in inter_calc.iterrows():
        wus_inputs = calc_wus_inputs(row)
        input_data = combine_inputs(wus_inputs,constant,line_data,solar_coeff)
        dlr = calc_dlr(input_data)
        dlrs.append(dlr)
    return dlrs

def find_bus_info(bus_x,buses):
    for line, data in buses.items():
        if data['bus0'] == bus_x:
            return line, 'bus0'
        elif data['bus1'] == bus_x:
            return line, 'bus1'
    return None, None

def get_line_dlrs_eff(data,constant_inputs,line_info,solar_coeff_data,bus_info,bus_list,n,wind_fix,phi_fix): # new version that calculates for all the buses first
    # data: wus weather data for clusters
    # constant_inputs, solar_coeff_data: constants for dlr equations
    # line_info: info about lines in tamu system
    # ONLY WORKS IF PHI IS FIXED!!
    
    dlrs = {}
    logger.info("Starting DLR calculations...")
    start = time.time()
    
    bus_dlrs = {}
    # calculate all buses first and store
    logger.info("Calculating DLRs for all buses")
    i = 0
    b = len(bus_list)
    for bus in bus_list:
        bus_start = time.time()
        line, bus_label = find_bus_info(bus,bus_info)
        line_data = line_info[line][bus_label]
        inter_calc = get_wus_var(data,bus_info[line][bus_label],line_data['Z_l'],wind_fix,phi_fix)
        bus_dlrs[bus] = get_dlr(inter_calc,line_data,constant_inputs,solar_coeff_data)
        i = i + 1
        if i % (b // 20) == 0 or i == b:
            check = time.time()
            logger.info(f'{(i/b)*100:.3g}% done! (line time: {(check - bus_start):.2f}, walltime: {(check - start):.2f})')
    
    logger.info("Beginning line iterations")
    i = 0
    for key, value in line_info.items():
        #logger.info(f"Calculating DLR for line {key}...")
        key_start = time.time()
        line_bus_dlrs = {}
        for bus in ['bus0','bus1']:
            bus_name = bus_info[key][bus]
            line_bus_dlrs[bus] = bus_dlrs[bus_name]
        dlrs[key] = get_worse(line_bus_dlrs['bus0'],line_bus_dlrs['bus1'])
        i = i + 1
        if i % (n // 20) == 0 or i == n:
            check = time.time()
            logger.info(f'{(i/n)*100:.3g}% done! (line time: {(check - key_start):.2f}, walltime: {(check - start):.2f})')
    finish = time.time()
    logger.info(f'get_line_dlrs ran in {(finish - start):.2f}s')
    return dlrs

def gen_year(year):
    dt = pd.date_range(start=f'{year}-01-01 00:00:00', end=f'{year}-12-31 23:00:00', freq='h')
    return dt

def gen_year_no_leap(year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    
    if is_leap:
        # Generate Jan 1 - Feb 28
        dt1 = pd.date_range(start=f'{year}-01-01 00:00:00', 
                           end=f'{year}-02-28 23:00:00', freq='h')
        # Generate Mar 1 - Dec 31 (skipping Feb 29)
        dt2 = pd.date_range(start=f'{year}-03-01 00:00:00', 
                           end=f'{year}-12-31 23:00:00', freq='h')
        dt = dt1.union(dt2)
    else:
        dt = pd.date_range(start=f'{year}-01-01 00:00:00', 
                          end=f'{year}-12-31 23:00:00', freq='h')
    
    return dt

def get_fixed_year(
    wus_data_horizon: xr.Dataset,
    wus_data_prev: xr.Dataset,
    year: int
    ) -> xr.Dataset:
    # dealing with the fact that the WUS data starts on 09-01-year
    logger.info(f"Extracting and combining data for {year}")

    begin = wus_data_prev.isel(Times=slice(-5833, -1))
    end = wus_data_horizon.isel(Times=slice(0, 2928))
    combined = xr.concat([begin,end],dim='Times')

    return combined 