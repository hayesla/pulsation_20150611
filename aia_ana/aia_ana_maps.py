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
from aiapy.calibrate import register, update_pointing

date = '2014-06-11'

savedir = '/Users/laurahayes/QPP/interesting_event_2014-611/aia_prepped/'

def prep_files(wave):
	aia_files_131 = glob.glob('/Users/laurahayes/QPP/interesting_event_2014-611/aia_full/*{:s}*.fits'.format(wave))
	aia_files_131.sort()
	maps = []
	for f in aia_files_131:
		try:
			m = sunpy.map.Map(f)
			maps.append(m)
		except:
			print(f)



	map_list = []
	for m in maps:
		m_updated_pointing = update_pointing(m)

		m_registered = register(m_updated_pointing)

		#m_normalized = sunpy.map.Map(m_registered.data/m_registered.exposure_time.to(u.s).value,
		    						 #m_registered.meta)

		map_list.append(m_registered)

	map_list = sunpy.map.Map(map_list, sequence=True)

	for m in map_list:
		m.save(savedir + 'aia_{:s}_'.format(wave)+str(m.date) + '.fits')

def do_all_prep():
	prep_files('94')
	prep_files('1600')
	prep_files('1700')
	prep_files('304')
	prep_files('171')
	prep_files('335')
	prep_files('193')
	prep_files('211')

