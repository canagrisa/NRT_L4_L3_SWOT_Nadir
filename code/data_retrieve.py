import os
import ftplib
from ftplib import FTP
import xarray as xr
import utils
import numpy as np


def download_nrt(username,
                 password,
                 year,
                 month,
                 day,
                 coords,
                 product='altimetry'):
    """
    Download NRT data from the Copernicus Marine Environment Monitoring Service (CMEMS)

    Specifying the product argument, the user can choose between:

    - 'altimetry': Altimetric L4-NRT data (https://data.marine.copernicus.eu/product/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/description)
    - 'sst': SST L4-NRT data (https://data.marine.copernicus.eu/product/SST_MED_SST_L3S_NRT_OBSERVATIONS_010_012/description)
    - 'chlorophyll': Chlorophyll L4-NRT data (https://data.marine.copernicus.eu/product/OCEANCOLOUR_MED_BGC_L3_NRT_009_141/description)

    """

    if not isinstance(year, (int, str)):
        raise ValueError(
            "Year argument must be an integer or a numeric string")
    if not isinstance(month, (int, str)):
        raise ValueError(
            "Month argument must be an integer or a numeric string")
    if not isinstance(day, (int, str)):
        raise ValueError("Day argument must be an integer or a numeric string")

    if isinstance(month, int):
        month = f"{month:02d}"
    if isinstance(day, int):
        day = f"{day:02d}"

    if product == 'altimetry' or product == 'ugos' or product == 'vgos' or product == 'velocity':
        print('Retrieving Altimetric L4-NRT data for {}-{}-{}...'.format(year, month, day))
        service_id = 'SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060'
        product_id = 'dataset-duacs-nrt-europe-merged-allsat-phy-l4'
        idx_i = -14
        idx_f = -12
        lon = 'longitude'
        lat = 'latitude'
    if product == 'sst':
        print('Retrieving SST L4-NRT data for {}-{}-{}...'.format(year, month, day))
        service_id = 'SST_MED_SST_L3S_NRT_OBSERVATIONS_010_012'
        product_id = 'SST_MED_SST_L3S_NRT_OBSERVATIONS_010_012_b'
        idx_i = 6
        idx_f = 8
        lon = 'lon'
        lat = 'lat'
    if product == 'chlorophyll':
        print('Retrieving chlorophyll L4-NRT data for {}-{}-{}...'.format(year, month, day))
        service_id = 'OCEANCOLOUR_MED_BGC_L3_NRT_009_141'
        product_id = 'cmems_obs-oc_med_bgc-plankton_nrt_l3-olci-300m_P1D'
        idx_i = 6
        idx_f = 8
        lon = 'lon'
        lat = 'lat'
    # if product == 'wind':
    #     print('Retrieving wind L4-NRT data for {}-{}-{}...'.format(year, month, day))
    #     service_id = 'WIND_GLO_PHY_L4_NRT_012_004'
    #     product_id = 'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H'
    #     idx_i = -23
    #     idx_f = -21
    #     lon = 'lon'
    #     lat = 'lat'


        
    # Connect to the ftp server
    ftp = FTP('nrt.cmems-du.eu')
    ftp.login(user=username, passwd=password)

    # Navigate to the directory of interest
    ftp = FTP('nrt.cmems-du.eu', username, password)

    fold = 'Core/' + service_id + '/' + \
        product_id + '/' + str(year) + '/' + month
    ftp.cwd(fold)

    # Get a list of all the files in the directory
    files = ftp.nlst()

    # Filter the list to only include the files we want and download into a temporary folder
    for file in files:
        if file[idx_i:idx_f] == day:
            out_fold = f'../temp/{product}/'
            path = out_fold+file
            if not os.path.exists(out_fold):
                os.makedirs(out_fold)
            if not os.path.isfile(path):
                with open(path, 'wb') as f:
                    ftp.retrbinary("RETR " + file, f.write)

            # Load the data into an xarray dataset
            ds = xr.open_dataset(path)

            ds = ds.where((ds[lon] > coords[0]) &
                          (ds[lon] < coords[1]) &
                          (ds[lat] > coords[2]) &
                          (ds[lat] < coords[3]),
                          drop=True)

            if product == 'sst':
                ds['sea_surface_temperature'] = ds['sea_surface_temperature'] - 273.15
                ds = ds['sea_surface_temperature']
                ds = ds.isel(time=0)

            elif product == 'chlorophyll':
                ds = ds['CHL']
                ds = ds.isel(time=0)

            elif product == 'altimetry':
                ds = ds['adt']
                ds = ds.isel(time=0)

            elif product == 'ugos':
                ds = ds['ugos']*100
                ds = ds.isel(time=0)

            elif product == 'vgos':
                ds = ds['vgos']*100
                ds = ds.isel(time=0)

            elif product == 'velocity':
                selected_vars = ["vgos",
                                 "ugos"]
                #ds = np.sqrt((ds['vgos']*100)**2 + (ds['ugos']*100)**2)
                ds = ds.isel(time=0)

            elif product == 'wind':
                selected_vars = ["eastward_stress",
                                 "northward_stress",
                                 'stress_curl']
                ds = ds[selected_vars]
                ds = ds.mean('time')

            return ds


def download_along_track(username,
                         password,
                         year,
                         month,
                         day,
                         coords):
    """
    Download along track data from the Copernicus Marine Environment Monitoring Service (CMEMS)

    Link: https://data.marine.copernicus.eu/product/SEALEVEL_EUR_PHY_L3_NRT_OBSERVATIONS_008_059/description
    """

    if not isinstance(year, (int, str)):
        raise ValueError(
            "Year argument must be an integer or a numeric string")
    if not isinstance(month, (int, str)):
        raise ValueError(
            "Month argument must be an integer or a numeric string")
    if not isinstance(day, (int, str)):
        raise ValueError("Day argument must be an integer or a numeric string")

    if isinstance(month, int):
        month = f"{month:02d}"
    if isinstance(day, int):
        day = f"{day:02d}"

    sats = utils.sats_dic
    sats = dict(sorted(sats.items()))

    # Connect to the ftp server
    ftp = FTP('nrt.cmems-du.eu')
    ftp.login(user=username, passwd=password)

    # Navigate to the directory of interest
    service_id = 'SEALEVEL_EUR_PHY_L3_NRT_OBSERVATIONS_008_059'

    datasets = {}

    print('Loading along track data:', end=' ')

    for sat in sats:

        ftp = FTP('nrt.cmems-du.eu', username, password)

        print('{}'.format(sats[sat]), end=' ')

        if sat in ['s6a-hr', 'j3n', 'h2c']:
            product_id = 'cmems_obs-sl_eur_phy-ssh_nrt_{}-l3-duacs_PT1S'.format(
                sat)
        else:
            product_id = 'dataset-duacs-nrt-europe-{}-phy-l3'.format(sat)

        fold = 'Core/' + service_id + '/' + \
            product_id + '/' + str(year) + '/' + month
        ftp.cwd(fold)

        files = ftp.nlst()

        for file in files:
            if file[-14:-12] == day:  # Could be improved
                out_fold = f'../temp/{sat}/'
                path = out_fold+file
                if not os.path.exists(out_fold):
                    os.makedirs(out_fold)
                if not os.path.isfile(path):
                    with open(path, 'wb') as f:
                        ftp.retrbinary("RETR " + file, f.write)

        ds = xr.open_dataset(path)

        ds = ds.where((ds.longitude > coords[0]) &
                      (ds.longitude < coords[1]) &
                      (ds.latitude > coords[2]) &
                      (ds.latitude < coords[3]),
                      drop=True)

        if len(ds.time) > 1:
            ds['adt'] = ds.sla_filtered + ds.mdt
            datasets[sat] = ds
            ds.attrs['title'] = f"{sats[sat]}"

    datasets = utils.split_tracks(datasets)

    if not datasets:
        print('\n - No along track altimetric data found for the specified date and region')

    return datasets


def download_along_track_swot(month, day,
                              coords,
                              out_fold='../temp/swot_l3/'):
    
    """

    Download along track SWOT-nadir data from the SWOT-ADAC data server (by 25/03/2023)

    'These data are generated from non-validated beta Level-2 products from the SWOT nadir /
    instrument during the early phases of the mission'

    Link: https://www.swot-adac.org/resources/swot-adac-products-access/
    """

    print('Loading along track SWOT-nadir data...')

    ftp = ftplib.FTP('ftp-access.aviso.altimetry.fr')
    ftp.login(user='swot_adac', passwd='swot_adac2023')
    ftp.cwd('data/Data/ALTI/DUACS_SWOT_Nadir/')
    ftp.cwd('L3_Along_track')

    files = ftp.nlst()

    for filename in files:
        path = out_fold+filename
        if not os.path.exists(out_fold):
            os.makedirs(out_fold)
        if not os.path.isfile(path):
            with open(path, 'wb') as file:
                ftp.retrbinary('RETR ' + filename, file.write)

    ds =  xr.open_mfdataset(out_fold+'*.nc')
    ds = ds.sel(time=ds.time.dt.month.isin([int(month)]))
    ds = ds.sel(time=ds.time.dt.day.isin([int(day)]))
    ds = ds.where((ds.longitude > coords[0]) &
                    (ds.longitude < coords[1]) &
                    (ds.latitude > coords[2]) &
                    (ds.latitude < coords[3]),
                    drop=True)

    if len(ds.time) > 0:
        ds['adt'] = ds.sla_unfiltered + ds.mdt
        ds = utils.split_tracks({'SWOT': ds})
        return ds
    else:
        return {}


def get_data(username, password, year, month, day, coords, product):

    # Download L4, L3 and SWOT-nadir data all at once into a temporary file

    swot = download_along_track_swot(month, day, coords)
    ds_nrt = download_nrt(username, password, year,
                          month, day, coords, product)
    sat_dic = download_along_track(
        username, password, year, month, day, coords)

    return ds_nrt, swot, sat_dic
