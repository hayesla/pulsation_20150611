import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import datetime
from matplotlib import dates
from sunpy import timeseries as ts
from sunpy.time import parse_time
from scipy.io import readsav
from scipy import fftpack
import sys
sys.path.append('/Users/lahayes/afino_analysis/')
#from analyse_generic_series import analyse_series
ref_time = parse_time('1979-01-01').datetime


flare_ts = '2014-06-11 05:30'
flare_te = '2014-06-11 05:45'
def other():
	goes = ts.TimeSeries('go1520140611.fits')

	filename = '11jun14.lis.gz'

	rstn = np.genfromtxt(filename, delimiter=2*[4]+5*[2]+8*[7], dtype=('|S10', int, int, int, int, int, int,
	                            float, float,float, float,float, float,float, float),
	                            names = ['sta','year','mon','day','hour','min','sec','f1','f2','f3','f4','f5','f6',
	                                      'f7','f8'])

	times = list(map(datetime.datetime,rstn['year'],rstn['mon'],rstn['day'],
	                                           rstn['hour'],rstn['min'],rstn['sec']))

	data = np.transpose([rstn['f1'],rstn['f2'],rstn['f3'],rstn['f4'],rstn['f5'],rstn['f6'],rstn['f7'],rstn['f8']])
	df = pd.DataFrame(data, columns=['245 MHz','410 MHz','610 MHz','1.4 GHz',
	                                                          '2.7 GHz','4.9 GHz','8.8 GHz','15.4 GHz'], index = times)


	# save data to a csv
	df.sort_index(inplace=True)
	df.to_csv('san_vito_rstn_11062014.csv', header=True, index=True)


	new_df = df.truncate(flare_ts, flare_te)

	short_goes = goes.truncate(flare_ts, flare_te)
	gl = short_goes.data['xrsb']
	gs = short_goes.data['xrsa']


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

def make_rhessi_lc(file):
    #file = 'hsi_spectrum_20130515_012024.fits'
    a = fits.open(file)
    start_time = a[0].header['DATE_OBS']
    t_start = datetime.datetime.strptime(start_time[0:10] + ' '+start_time[11:], '%Y-%m-%d %H:%M:%S.%f')
    start_time_day = datetime.datetime.strptime(str(t_start)[0:10]+' 00:00:00', '%Y-%m-%d %H:%M:%S')
    #print a[1].data.columns

    time = a[1].data['TIME']
    time_array = []
    for i in range(len(time)):
        time_array.append(start_time_day + datetime.timedelta(seconds = time[i]))


    e_min = a[2].data['E_MIN']
    e_max = a[2].data['E_MAX']
    channel = zip(list(e_min), list(e_max))
    dif_channels = a[1].data['RATE'].T

    r36 = dif_channels[0]
    r612 = dif_channels[1]
    r1225 = dif_channels[2]
    r2550 = dif_channels[3]
    r50100 = dif_channels[4]
    r100300 = dif_channels[5]


    df = pandas.DataFrame({'3-6 keV': r36, '6-12 keV': r612, '12-25 keV': r1225, '25-50 keV': r2550,'50-100 keV': r50100, '100-300 keV': r100300},index=time_array)
    
    rhessi_lc = lc.LightCurve.create(df)
    return rhessi_lc

def read_norp(filey):
	a = readsav(filey, python_dict = True)

	fv = a['fv'].T
	fvavg = a['fvavg']
	fiavg = a['fiavg']
	tim = a['tim']
	mvd = a['mvd']
	fi = a['fi'].T
	freq = a['freq']
	day = a['day']

	datee = parse_time('1979-01-01 00:00').datetime + datetime.timedelta(days = int(day.DAY)-1)
	timez = []
	for i in range(len(tim.TIME)):
		tt = float(tim.TIME[i])/1000.
		timez.append(datee + datetime.timedelta(seconds = tt))


	n1 = fi[0] - fiavg[0]
	n2 = fi[1] - fiavg[1]
	n3 = fi[2] - fiavg[2]
	n9 = fi[3] - fiavg[3]
	n17 = fi[4] - fiavg[4]
	n35 = fi[5] - fiavg[5]
	n80 = fi[6] - fiavg[6]

	norp1 = pd.Series(n1, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp2 = pd.Series(n2, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp3 = pd.Series(n3, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp9 = pd.Series(n9, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp17 = pd.Series(n17, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp35 = pd.Series(n35, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()
	norp80 = pd.Series(n80, index = timez).truncate(flare_ts, flare_te)#.resample('1s').mean()

	df = {'norp1':norp1, 'norp2':norp2,'norp3':norp3,'norp9':norp9,'norp17':norp17,'norp35':norp35,'norp80':norp80}
	norp_df = pd.DataFrame(df)
	return norp_df, norp1, norp2, norp3, norp9, norp17, norp35, norp80

norp_df, norp1, norp2, norp3, norp9, norp17, norp35, norp80 = read_norp('norp20140611_0534.xdr')
#fig, ax = plt.subplots()


# rhessi = ts.TimeSeries('/Users/lahayes/sunpy/data/hsi_obssumm_20140611_045.fits')
# rhessi = rhessi.truncate(flare_ts, flare_te)
# r2550 = rhessi.data['25 - 50 keV']


def plot_norp_rhessi():
	fig, ax = plt.subplots()

	ax.plot(norp17, label='NoRP 17GHz', color='r', drawstyle='steps-mid')
	ax.set_xlim('2014-06-11 05:32:30', '2014-06-11 05:36:30')
	ax.legend(loc='upper left')
	ax.set_ylabel('Flux (SFU)')
	ax.set_xlabel('Time (UT) 2014-06-11')

	ax.xaxis.set_minor_locator(dates.SecondLocator(interval=10))
	ax.xaxis.set_major_formatter(dates.DateFormatter('%H%M%S'))
	ax.xaxis.set_tick_params(rotation=45)
	ax.tick_params(which='both', direction='in')
	ax2 = ax.twinx()
	ax2.plot(r2550, drawstyle='steps-mid', label='RHESSI 25-50keV', color='k', lw=0.8)
	ax2.legend(loc = 'upper right')
	ax2.set_ylabel('Counts det$^{-1}$ s$^{-1}$')
	ax2.tick_params(which='both', direction='in')
	plt.tight_layout()
	plt.savefig('rhessi_norp_event.png', dpi=200)
	plt.show()


save_dir = '/Users/lahayes/QPP/interesting_event_2014-611'
def analyse_afino(ser, descc, title):

	time_ar = []
	for t in ser.index:
		time_ar.append((parse_time(t).datetime - ref_time).total_seconds())
	
	timez = np.array(time_ar)
	analyse_series(timez, np.array(ser), description = descc, title = title, savedir = save_dir)

# analyse_afino(norp17, 'norp17', 'NoRP 17GHz')
# analyse_afino(gl, 'gl', 'GOES 1-8~$\mathrm{\AA}$')
# analyse_afino(r2550, 'r2550', 'RHESSI 25-50keV')


