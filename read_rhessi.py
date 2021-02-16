from astropy.io import fits 
import datetime
import numpy as np 
import pandas as pd 

def read_rhessi(file):

	#file = 'hsi_spectrum_20140611_052832.fits'
	a = fits.open(file)

	start_time = a[0].header['DATE_OBS']
	t_start = datetime.datetime.strptime(start_time[0:10] + ' '+start_time[11:], '%Y-%m-%d %H:%M:%S.%f')
	start_time_day = datetime.datetime.strptime(str(t_start)[0:10]+' 00:00:00', '%Y-%m-%d %H:%M:%S')


	time = a[1].data['TIME']
	time_array = []
	for i in range(len(time)):
	    time_array.append(start_time_day + datetime.timedelta(seconds = time[i]))


	e_min = a[2].data['E_MIN']
	e_max = a[2].data['E_MAX']
	channel = list(zip(list(e_min), list(e_max)))

	dif_channels = a[1].data['RATE'].T
	atten_state = a[3].data['INTERVAL_ATTEN_STATE$$STATE'] 


	return time_array, dif_channels, e_min, e_max, atten_state

def find_closest(arr, val):
	return np.argmin(np.abs(arr -val))

