import urllib
from bs4 import BeautifulSoup
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

save_dir = "/Users/lahayes/QPP/interesting_event_2014-611/ssw_fits/"

timerange = TimeRange('2014-06-11 05:30:00','2014-06-11 05:36:00')

url = 'https://hesperia.gsfc.nasa.gov/sdo/aia/2014/06/11/20140611_0528-0547/'

resp = urllib.request.urlopen(url)
soup = BeautifulSoup(resp)


def find_url_waves(wave):
	file_link = []
	for link in soup.find_all('a', href=True):
		if link['href'].endswith('{:d}_.fts'.format(wave)):
			file_link.append(link['href'])
	return file_link
# Use Scraper!

def list_files(url):

	resp = urllib.request.urlopen(url)
	soup = BeautifulSoup(resp)

	file_link = []
	for link in soup.find_all('a', href=True):

		file_link.append(link['href'])
	return file_link

def get_files(wave):
	sswbrowser_pattern = ('https://hesperia.gsfc.nasa.gov/sdo/aia/2014/06/11/20140611_0528-0547/'
	                      'ssw_cutout_%Y%m%d_%H%M%S_aia_{wave}_.fts')

	ssw = Scraper(sswbrowser_pattern, wave=wave)

	timerange = TimeRange('2014-06-11 05:30:00','2014-06-11 05:36:00')

	aia_131_files = ssw.filelist(timerange)

	return aia_131_files

def save_files(wave):
	files = get_files(wave)

	for f in files:
		urllib.request.urlretrieve(f, save_dir + f.split('/')[-1])

save_files(304)
save_files(1600)
save_files(1700)
save_files(94)