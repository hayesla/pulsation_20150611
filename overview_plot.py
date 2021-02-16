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
from scipy.io import readsav


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


euve = readsav('lyman_alpha.sav')
euve_times = [parse_time(euve['UTBASE']+euve['tarray'][i], format='utime').datetime for i in range(len(euve['tarray']))]
eve_gt = pd.Series(euve['yclean'][0], index=euve_times)

evve = eve_gt.replace(-99999.0, np.nan).truncate(pul_ts, pul_te)

def plot_tog():
	fig = plt.figure(figsize=(7, 10))

	gss = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[1, 2,2])
	ax1 = plt.subplot(gss[0])
	ax2 = plt.subplot(gss[1])
	ax3 = plt.subplot(gss[2], sharex=ax2)

	ax1.plot(gl, color='r', label='GOES 1-8$\mathrm{\AA}$')
	ax1.plot(gs, color='b', label='GOES 1-8$\mathrm{\AA}$')
	ax1.set_yscale('log')
	ax1.set_xlim(flare_ts, flare_te)
	ax1.xaxis.set_minor_locator(dates.SecondLocator(interval=10))
	ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=1))
	ax1.legend(loc='upper right')
	ax1.tick_params(labelbottom=False, which='both', direction='in')
	ax1.set_ylabel('Flux (Wm$^{-2}$)')
	ax1.axvline(pul_ts, color='k')
	ax1.axvline(pul_te, color='k')


	ax2.plot(norp17, label='NoRP 17GHz', color='darkred', drawstyle='steps-mid')
	ax2.plot(norp9, label='NoRP 9GHz', color='darkblue', drawstyle='steps-mid')
	ax2.legend(loc='upper right')
	ax2.set_ylabel('Flux (SFU)')
	ax2.set_ylim(0, 200)


	x1 = dates.date2num(datetime.datetime.strptime(pul_ts, '%Y-%m-%d %H:%M:%S')) # start time
	x2 = dates.date2num(datetime.datetime.strptime(pul_te, '%Y-%m-%d %H:%M:%S')) # end time

	xyA = (x1, -0.01)
	xyB = (0.0, 1)  # x and y in axes coordinates
	coordsA = ax1.get_xaxis_transform() # x in data coords, y in axes coords
	coordsB = "axes fraction"
	con_start = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB,
						        axesA=ax1, axesB=ax2,
	                            arrowstyle="-")
	xyC = (x2, -0.01)
	xyD = (1, 1)  # x and y in axes coordinates
	coordsC = ax1.get_xaxis_transform() # x in data coords, y in axes coords
	coordsD = "axes fraction"
	con_end = ConnectionPatch(xyA=xyC, xyB=xyD, coordsA=coordsC, coordsB=coordsD,
						     axesA=ax1, axesB=ax2,
	                         arrowstyle="-")



	ax2.add_artist(con_start)
	ax2.add_artist(con_end)


	#ax3 = ax2.twinx()
	ax3.plot(rhessi_2535, drawstyle='steps-mid', label='RHESSI 25-35keV', color='grey', lw=0.8)
	ax3.plot(rhessi_35100, drawstyle='steps-mid', label='RHESSI 35-100keV', color='k', lw=0.8)
	ax3.legend(loc = 'upper right')
	ax3.set_ylabel('Counts det$^{-1}$ s$^{-1}$')


	ax3.set_xlabel('Time (UT) 2014-06-11')
	ax3.set_xlim(pul_ts, pul_te)
	ax2.xaxis.set_minor_locator(dates.SecondLocator(interval=10))
	ax2.xaxis.set_major_locator(dates.MinuteLocator(interval=1))
	ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
	ax2.tick_params(which='both', direction='in', labelbottom=False)
	#ax3.xaxis.set_tick_params(rotation=45, which='both', direction='in')
	ax3.xaxis.set_tick_params(which='both', direction='in')


	plt.tight_layout()
	plt.subplots_adjust(hspace=0.1)
	plt.savefig('overview_test.png', dpi=200)
	plt.close()


def plot_tog2():
	fig = plt.figure(figsize=(7, 8))

	gss = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1, 2])
	ax1 = plt.subplot(gss[0])
	ax2 = plt.subplot(gss[1])

	ax1.plot(gl, color='r', label='GOES 1-8$\mathrm{\AA}$')
	ax1.plot(gs, color='b', label='GOES 1-8$\mathrm{\AA}$')
	ax1.set_yscale('log')
	ax1.set_xlim(flare_ts, flare_te)
	ax1.xaxis.set_minor_locator(dates.SecondLocator(interval=10))
	ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=1))
	ax1.legend(loc='upper right')
	ax1.tick_params(labelbottom=False, which='both', direction='in')
	ax1.set_ylabel('Flux (Wm$^{-2}$)')
	ax1.axvline(pul_ts, color='k')
	ax1.axvline(pul_te, color='k')


	ax2.plot(norp17, label='NoRP 17GHz', color='darkred', drawstyle='steps-mid')
	#ax2.plot(norp9, label='NoRP 9GHz', color='darkblue', drawstyle='steps-mid')
	ax2.plot(np.nan, label='RHESSI 25-35keV', color='k')
	ax2.plot(np.nan, label='RHESSI 35-100keV', color='darkblue')

	ax2.legend(loc='upper right')
	ax2.set_ylabel('Flux (SFU)')
	ax2.set_ylim(0, 180)


	x1 = dates.date2num(datetime.datetime.strptime(pul_ts, '%Y-%m-%d %H:%M:%S')) # start time
	x2 = dates.date2num(datetime.datetime.strptime(pul_te, '%Y-%m-%d %H:%M:%S')) # end time

	xyA = (x1, -0.01)
	xyB = (0.0, 1)  # x and y in axes coordinates
	coordsA = ax1.get_xaxis_transform() # x in data coords, y in axes coords
	coordsB = "axes fraction"
	con_start = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB,
						        axesA=ax1, axesB=ax2,
	                            arrowstyle="-")
	xyC = (x2, -0.01)
	xyD = (1, 1)  # x and y in axes coordinates
	coordsC = ax1.get_xaxis_transform() # x in data coords, y in axes coords
	coordsD = "axes fraction"
	con_end = ConnectionPatch(xyA=xyC, xyB=xyD, coordsA=coordsC, coordsB=coordsD,
						     axesA=ax1, axesB=ax2,
	                         arrowstyle="-")



	ax2.add_artist(con_start)
	ax2.add_artist(con_end)


	ax3 = ax2.twinx()
	ax3.plot(rhessi_2535, drawstyle='steps-mid', label='RHESSI 25-35keV', color='k', lw=0.8)
	ax3.plot(rhessi_35100, drawstyle='steps-mid', label='RHESSI 35-100keV', color='darkblue', lw=0.8)
	#ax3.legend(loc = 'upper right')
	ax3.set_ylabel('Counts det$^{-1}$ s$^{-1}$')


	ax3.set_xlabel('Time (UT) 2014-06-11')
	ax3.set_xlim(pul_ts, pul_te)
	ax2.xaxis.set_minor_locator(dates.SecondLocator(interval=10))
	ax2.xaxis.set_major_locator(dates.MinuteLocator(interval=1))
	ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
	ax2.tick_params(which='both', direction='in', labelbottom=False)
	#ax3.xaxis.set_tick_params(rotation=45, which='both', direction='in')
	ax3.xaxis.set_tick_params(which='both', direction='in')


	plt.tight_layout()
	plt.subplots_adjust(hspace=0.1)
	plt.savefig('overview_test2.png', dpi=200)
	plt.close()


