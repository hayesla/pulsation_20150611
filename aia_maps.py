import sunpy.map 
from sunpy.time import parse_time
from sunpy.coordinates import frames
from astropy.time import TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u 
from scipy.io import readsav
import matplotlib.pyplot as plt 
import numpy as np 
import subprocess
import glob
import mapcube_tools

date = '2014-06-11'


"""
-----------------

Read in the data
-----------------

"""
# aia_files_1700 = glob.glob('/Users/laurahayes/QPP/interesting_event_2014-611/lmsal_ssw/*1700_.fts')
# aia_files_1700.sort()
# aia_maps_1700 = sunpy.map.Map(aia_files_1700, sequence=True)

# aia_files_131 = glob.glob('/Users/laurahayes/QPP/interesting_event_2014-611/lmsal_ssw/*131_.fts')
# aia_files_131.sort()
# aia_maps_131 = sunpy.map.Map(aia_files_131, sequence=True)

# aia_files_304 = glob.glob('/Users/laurahayes/QPP/interesting_event_2014-611/lmsal_ssw/*304_.fts')
# aia_files_304.sort()
# aia_maps_304 = sunpy.map.Map(aia_files_304, sequence=True)

aia_files_131 = glob.glob('/Users/laurahayes/QPP/interesting_event_2014-611/aia_full/*131*.fits')