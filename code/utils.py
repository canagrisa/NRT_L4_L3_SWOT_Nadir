import numpy as np
import os


# Name of all current operational altimetric satellites via CMEMS per 25-03-2023
sats_dic = {'al': 'AltiKa',
            'c2n': 'Cryosat-2',
            's3a': 'Sentinel-3A',
            's3b': 'Sentinel-3B',
            's6a-hr': 'Sentinel-6A',
            'h2b': 'Haiyang-2B',
            'j3n': 'Jason-3'}

# Region of the Balearic Sea of interest for FaSt-SWOT campain
lon_min = 1; lon_max = 5; lat_min = 37.5; lat_max = 40.5
bal_coords = [lon_min, lon_max, lat_min, lat_max]


def dirtodict(dirPath):
    #From a given folder, create a dictionary with a given folder & file tree structure
    d = {}
    for i in [os.path.join(dirPath, i) for i in os.listdir(dirPath)
              if os.path.isdir(os.path.join(dirPath, i))]:
        d[os.path.basename(i)] = dirtodict(i)
    d['.files'] = [os.path.join(dirPath, i) for i in os.listdir(dirPath)
                   if os.path.isfile(os.path.join(dirPath, i))]
    return d


def geov(ds):

    # Calculate the geostrophic velocity from the absolute dynamic topography data along track

    lon = ds.longitude
    lat = ds.latitude
    lat_r = lat*np.pi/180
    adt = ds['adt']

    # Constants
    g = 9.81  # m/s^2 - gravity
    omega = 7.292115e-5  # rad/s - angular velocity of the earth
    f = 2 * omega * np.sin(lat_r[:-1])  # Coriolis parameter (array)

    lat_km = 110574  # km - length of a degree of latitude
    lon_km = 111320*np.cos(lat_r)  # km - length of a degree of longitude (array)

    if lat[-1] > lat[0]:
        sign_lat = 1
    else:
        sign_lat = -1

    if lon[-1] > lon[0]:
        sign_lon = 1
    else:
        sign_lon = -1

    # Calculate the geostrophic velocity
    dadt = np.diff(adt)
    dlon = np.diff(lon) * lon_km[:-1]
    dlat = np.diff(lat) * lat_km

    alfa = - np.append(np.arctan(dlat/dlon), 0) 

    dl = np.sqrt(dlon**2 + dlat**2)

    velVec = + sign_lon * sign_lat * (g/f) * dadt / dl

    velVec = np.append(velVec, 0)

    vx = velVec * np.sin(alfa) * sign_lat 
    vy = velVec * np.cos(alfa) * sign_lat

    ds = ds.assign(geov=(['time'], velVec))
    ds = ds.assign(geovx=(['time'], vx))
    ds = ds.assign(geovy=(['time'], vy))

    return ds


def deal_with_gaps(ds):

    # Deal with gaps in the data
    # If there is a gap of more than 2.5 times the minimum time step, set the geov to 0
    # This is a bit of a hack, but it works

    time = ds.time
    dt = np.diff(time)
    dt_norm = dt/min(dt)
    idx = np.where(dt_norm > 2.5)[0]
    idx = np.concatenate((idx, (idx+1)))
    bad_times = time.isel(time=idx)
    ds['geov'] = ds.geov.where(np.invert(ds.time.isin(bad_times)), 0)
    ds['geovx'] = ds.geovx.where(np.invert(ds.time.isin(bad_times)), 0)
    ds['geovy'] = ds.geovy.where(np.invert(ds.time.isin(bad_times)), 0)

    return ds


def split_tracks(sats_day):

    """
    Given along track data from a single day, split the data into tracks if there is a gap of more than 1 hour
    """

    for sat in sats_day:

        tracks = {}

        ds = sats_day[sat]

        time = ds.time
        hour = 3600
        dt = (np.diff(time)*1e-9).astype('float32')
        idxs = np.where(dt > hour)[0]

        if len(idxs) <= 0:
            ds = geov(ds)
            ds = deal_with_gaps(ds)
            sats_day[sat] = ds

        else:
            idx = idxs[0]
            ds_1 = ds.isel(time=slice(0, idx+1))
            ds_1 = geov(ds_1)
            if len(ds_1.time)>1:
                ds_1 = deal_with_gaps(ds_1)

            ds_2 = ds.isel(time=slice(idx+1, None))
            #ds_2 = ds_2.assign(geov=(['time'], geov(ds_2)))
            ds_2 = geov(ds_2)
            if len(ds_2.time)>1:
                ds_2 = deal_with_gaps(ds_2)
            ds_2 = deal_with_gaps(ds_2)

            tracks[1] = ds_1
            tracks[2] = ds_2

            sats_day[sat] = tracks

    return sats_day


def get_min_max(ds_nrt, sat_dic, swot, product):

    """
    Get the common min and max value of the:
     - NRT L4 data (if it is ADT)
     - L3 altimetric sat data
     - L3 SWOT data

     This is necessary to set common limits to the colorbar in the plots
    """

    max_vals_l4 = []
    min_vals_l4 = []

    if product == 'velocity':

        ds_vel = np.sqrt((ds_nrt['ugosa']*100)**2 +  (ds_nrt['vgosa']*100)**2)
        l4_max = np.float32(ds_vel.max())
        l4_min = np.float32(ds_vel.min())

    else:
        l4_max = np.float32(ds_nrt.max())
        l4_min = np.float32(ds_nrt.min())

    max_vals_l4.append(l4_max)
    min_vals_l4.append(l4_min)

    max_vals_l3 = []
    min_vals_l3 = []


    for sat in sat_dic:
        if not isinstance(sat_dic[sat], dict):
            ds = sat_dic[sat]
            max_vals_l3.append(np.float32(ds['adt'].max()))
            min_vals_l3.append(np.float32(ds['adt'].min()))
        else:
            for swath in sat_dic[sat]:
                ds = sat_dic[sat][swath]
                max_vals_l3.append(np.float32(ds['adt'].max()))
                min_vals_l3.append(np.float32(ds['adt'].min()))


    max_vals_l3_swot = []
    min_vals_l3_swot = []

    if swot:
        if not isinstance(swot['SWOT'], dict):
            ds = swot['SWOT']
            max_vals_l3_swot.append(np.float32(ds['adt'].max()))
            min_vals_l3_swot.append(np.float32(ds['adt'].min()))
        else:
            for swath in swot['SWOT']:
                ds = swot['SWOT'][swath]
                max_vals_l3_swot.append(np.float32(ds['adt'].max()))
                min_vals_l3_swot.append(np.float32(ds['adt'].min()))

    if product != 'altimetry':

        vmax_l4 = max(max_vals_l4)
        vmin_l4 = min(min_vals_l4)

        if product == 'vgos' or product == 'ugos':
            abs_max = max(abs(vmax_l4), abs(vmin_l4))
            vmax_l4 = abs_max
            vmin_l4 = -abs_max

        if product == 'velocity':
            vmin_l4 = 0

        if sat_dic or swot:
            vmax_l3 = max(max_vals_l3 + max_vals_l3_swot)
            vmin_l3 = min(min_vals_l3 + min_vals_l3_swot)
        else:
            vmax_l3 = 1
            vmin_l3 = 0

        return vmax_l3, vmin_l3, vmax_l4, vmin_l4

    
    else:
        max_vals = max_vals_l4 + max_vals_l3 + max_vals_l3_swot
        min_vals = min_vals_l4 + min_vals_l3 + min_vals_l3_swot

        vmax = max(max_vals)
        vmin = min(min_vals)

        return vmax, vmin, vmax, vmin

    