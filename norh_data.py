import urllib 
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
import numpy as np 
import sunpy.map
from sunpy.time import parse_time 
from astropy.io import fits
from sunpy.coordinates import frames, get_earth, sun
import glob
from astropy import units as u
import pandas as pd 

save_dir = "/Users/lahayes/QPP/interesting_event_2014-611/norh_fuk/"
def list_files(url, type_data):

    resp = urllib.request.urlopen(url)
    soup = BeautifulSoup(resp, features="lxml")

    file_link = []
    for link in soup.find_all("a", href=True):
        if link["href"].startswith(type_data):
            file_link.append(link["href"])
    return file_link


def download_ipa():
    """
    Download NORH Fujiki data for `ipa`.

    Notes
    -----
    `ipa` - 17GHz (R+L)
    `ips` - 17GHz (R-L)
    `ipz` - 34GHz (R+L)

    """
    norh_event_url = "https://solar.nro.nao.ac.jp/norh/images/event/20140611_0534/steady_fujiki/"
    files_norh = list_files(norh_event_url, "ipa")

    for f in files_norh:
        urllib.request.urlretrieve(norh_event_url+ f, save_dir + f)


def make_map_ipa(file):
    """
    Function to make a sunpy.map.GenericMap from a NoRH fits file.
    This function fixes the header keywords so that it works within
    map.

    Parameters
    ----------
    file : `str`
        path of NoRH fits file to read into a map

    Returns
    -------
    `sunpy.map.GenericMap

    """
    hdu = fits.open(file)[0]
    header = hdu.header
    data = hdu.data

    header['CUNIT1'], header['CUNIT2'] = 'arcsec', 'arcsec'
    header['CTYPE1'], header['CTYPE2'] = 'HPLN-TAN', 'HPLT-TAN'
    header['date-obs'] = header['date-obs']+'T'+header['time-obs']
    
    # get observer location
    observer = get_earth(header['date-obs'])
    observer = observer.transform_to(frames.HeliographicStonyhurst(obstime=observer.obstime))

    header['hgln_obs'] = observer.lon.to_value(u.deg)
    header['hglt_obs'] = observer.lat.to_value(u.deg)
    header['dsun_obs'] = observer.radius.to_value(u.m)

    header['rsun_obs'] = sun.angular_radius(header['date-obs']).to('arcsec').value

    norh_map = sunpy.map.Map(data, header)
    norh_map.plot_settings['cmap'] = 'BuPu'
    #norh_map.plot_settings['norm'] = LogNorm()


    return norh_map

def make_lcs_from_map(maps):
    """
    Function ot make lightcurves from a list of sunpy.map.Map or a 
    sunpy.map.MapSequence.

    Parameters
    ----------
    maps : `list`, `sunpy.map.MapSequence`

    Returns
    -------
    pd.Series of summed data from maps
    """

    data, times = [], []
    for i in range(len(maps)):
        data.append(np.sum(maps[i].data))
        times.append(maps[i].date.datetime)
    return pd.Series(data, index=times)