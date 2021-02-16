from sunpy import timeseries as ts 
from read_rstn import read_norp
from read_rhessi import read_rhessi, find_closest
import matplotlib.pyplot as plt 
from matplotlib import dates, gridspec
from matplotlib.patches import ConnectionPatch
import numpy as np 
import datetime
from astropy.io import fits
import pandas as pd 
from sunpy.time import parse_time
import sys
sys.path.append('/Users/lahayes/afino_analysis/')
from analyse_generic_series import analyse_series
ref_time = parse_time('1979-01-01').datetime

flare_ts = '2014-06-11 05:30'
flare_te = '2014-06-11 05:40'

pul_ts = '2014-06-11 05:32:30'
pul_te = '2014-06-11 05:36:30'

goes = ts.TimeSeries('go1520140611.fits')
short_goes = goes.truncate(flare_ts, flare_te)
gl = short_goes.data['xrsb']
gs = short_goes.data['xrsa']

norp_df, norp1, norp2, norp3, norp9, norp17, norp35, norp80 = read_norp('norp20140611_0534.xdr')


rhessi_time, rhessi_arr, rhessi_emin, rhessi_emax, atten_state = read_rhessi('./rhessi_data/hsi_spectrum_20140611_052832.fits')
rhessi_atten = pd.Series(atten_state[0], index=rhessi_time)
rhessi_atten[rhessi_atten==0] = 1.1
def df_rhessi_kev(e_low, e_high):
	return pd.Series(np.sum(rhessi_arr[find_closest(rhessi_emin, e_low):find_closest(rhessi_emax, e_high)], axis=0), index=rhessi_time)
rhessi_612 = df_rhessi_kev(6, 12)
rhessi_1225 = df_rhessi_kev(12, 25)
rhessi_2550 = df_rhessi_kev(25, 50)
rhessi_2535 = df_rhessi_kev(25, 35)
rhessi_35100 = df_rhessi_kev(35, 100)


save_dir = '/Users/lahayes/QPP/interesting_event_2014-611/afino_plots'
def analyse_afino(ser, descc, title, **kwargs):

	time_ar = []
	for t in ser.index:
		time_ar.append((parse_time(t).datetime - ref_time).total_seconds())
	
	timez = np.array(time_ar)
	analyse_series(timez, np.array(ser), description=descc, title=title, savedir=save_dir, **kwargs)


from scipy import fftpack
def fourier_ana(x, dt):
	N = len(x)
	x = x*np.hanning(N)
	dt = dt
	df = 1./(N*dt)
	PSD = abs(dt*fftpack.fft(x)[:int(N/2)])**2

	f = df*np.arange(N/2)
	if len(PSD) != len(f):
		PSD = abs(dt*fftpack.fft(x)[:int(N/2) + 1])**2
	return f, PSD

f_rh25, p_r25 = fourier_ana(rhessi_2535.truncate(pul_ts, pul_te), 4)
f_rh35, p_r35 = fourier_ana(rhessi_35100.truncate(pul_ts, pul_te), 4)
f_norp17, p_norp17 = fourier_ana(norp17.resample('1s').mean().truncate(pul_ts, pul_te), 1)


fig, ax = plt.subplots(3, sharex=True, figsize=(5, 10))

ax[0].loglog(1./f_rh25, p_r25, label='RHESSI 25-35keV')
ax[1].loglog(1./f_rh35, p_r35, label='RHESSI 35-100keV')
ax[2].loglog(1./f_norp17, p_norp17, label='NoRH 17GHz')

per = 1./f_norp17
ax[2].set_xlim(per[1], per[-1])
for a in ax:
	a.axvline(14, color='grey', ls='dashed')
	a.set_ylabel('PSD')
	a.legend(loc='upper right')
plt.tight_layout()

plt.savefig('period.png', dpi=200)
plt.close()

fig, ax = plt.subplots(3, sharex=True, figsize=(5, 10))

ax[0].plot(rhessi_2535, label='RHESSI 25-35keV', drawstyle='steps-mid')
ax[1].plot(rhessi_35100, label='RHESSI 35-100keV', drawstyle='steps-mid')
ax[2].plot(norp17.resample('1s').mean(), label='NoRH 17GHz')

per = 1./f_norp17
ax[2].set_xlim(pul_ts, pul_te)
ax[2].xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
for a in ax:
	# a.axvline(14, color='grey', ls='dashed')
	# a.set_ylabel('PSD')
	a.legend()
plt.tight_layout()

plt.savefig('period_lcs.png', dpi=200)
plt.close()

