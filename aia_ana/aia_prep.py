import sunpy.map
import glob
import aiapy.psf
from aiapy.calibrate import register, update_pointing, normalize_exposure
import os
import time
import numpy as np 


wave = "131"

files = glob.glob("/Users/laurahayes/QPP/interesting_event_2014-611/aia_full/*{:s}a*".format(wave))
files.sort()

m = sunpy.map.Map(files)

## need to make negative values 0 before deconvolve
def zero_out_neg(m):
	return m._new_instance(np.where(m.data<0, 0, m.data), m.meta)

mm = [zero_out_neg(x) for x in m]


def do_prep(mapy):

	# calculate psf if local file doesn't exist
	if not os.path.exists("/Users/laurahayes/QPP/interesting_event_2014-611/aia_ana/psf_{:s}.npy".format(wave)):
		print("calculating psf as file not found")
		if isinstance(mm, list):
			psf = aiapy.psf.psf(mm[0].wavelength)
		else: 
			psf = aiapy.psf.psf(mm.wavelength)
		np.save("psf_{:s}".format(wave), psf)

	# load in psf array and deconvolve
	psf = np.load("psf_{:s}.npy".format(wave))
	t1 = time.time()
	m_deconvolved = aiapy.psf.deconvolve(mapy, psf=psf)
	t2 = time.time() - t1
	print("deconvolution done {:f}".format(t2))

	# aia_prep to level 1.5
	t1 = time.time()
	m_updated_pointing = update_pointing(m_deconvolved)
	m_registered = register(m_updated_pointing)
	m_normalized = normalize_exposure(m_registered)
	t3 = time.time() - t1

	def reduce_mem(m):
		return m._new_instance(m.data.astype(np.float32), m.meta)

	final_maps = reduce_mem(m_normalized) 

	basepath = "/Users/laurahayes/QPP/interesting_event_2014-611/aia_prepped_deconvolve/"
	filename = "prepped_aia_{:.0f}_{:s}.fits".format(final_maps.wavelength.value, 
													 final_maps.date.strftime("%Y%m%d_%H%M%S"))
	t1 = time.time()	
	final_maps.save(basepath+filename)
	t4 = time.time() - t1

	print("it took {:f} seconds to deconvolve".format(t2))
	print("it took {:f} seconds to prep".format(t3))
	print("it took {:f} seconds to write".format(t4))


