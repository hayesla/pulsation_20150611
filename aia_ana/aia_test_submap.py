import sunpy.map
import numpy as np 
import matplotlib.pyplot as plt 
import glob
from astropy import units as u 
from astropy.coordinates import SkyCoord

wave = "131"
files_old = glob.glob("/Users/laurahayes/QPP/interesting_event_2014-611/aia_full/*{:s}a*".format(wave))
files_old.sort()

files_new = glob.glob("/Users/laurahayes/QPP/interesting_event_2014-611/aia_prepped_deconvolve/*_{:s}_*".format(wave))
files_new.sort()

maps_old = sunpy.map.Map(files_old)
maps_new = sunpy.map.Map(files_new)

def submap(m):
	tr = SkyCoord(700*u.arcsec, -150*u.arcsec, frame=m.coordinate_frame)
	bl = SkyCoord(500*u.arcsec, -300*u.arcsec, frame=m.coordinate_frame)
	return m.submap(bl, top_right=tr)

subs_old = [submap(m) for m in maps_old]
subs_new = [submap(m) for m in maps_new]


def plot_compare(m1, m2, i=0, name="test"):
	fig = plt.figure(figsize=(11, 4))
	ax1 = fig.add_subplot(1, 2, 1, projection=m1)
	ax2 = fig.add_subplot(1, 2, 2, projection=m2)

	m1 = sunpy.map.Map(m1.data/m1.exposure_time.value, m1.meta)

	m1.plot(axes=ax1, vmin=0, vmax=3000)
	m2.plot(axes=ax2, vmin=0, vmax=3000)
	plt.savefig("./plot_tests/{:s}_{:03d}.png".format(name, i))
	plt.close()


for i in range(len(subs_old)):
	plot_compare(subs_old[i], subs_new[i], i=i, name="test3")
