import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
import numpy as np 
import glob 
import sunpy.map 
from astropy import units as u 
from astropy.coordinates import SkyCoord
import subprocess
import mapcube_tools
from norh_data import make_map_ipa, make_lcs_from_map

norh_files = "/Users/lahayes/QPP/interesting_event_2014-611/norh_fuk/"
ipa_files = glob.glob(norh_files + 'ipa*')
ipa_files.sort()


# get NORH 17GHz maps
my_maps = []
for f in ipa_files[143:350]:
    m = make_map_ipa(f)
    my_maps.append(m)
my_maps = sunpy.map.Map(my_maps, sequence=True)

# get lightcurve from maps
norh_lc = make_lcs_from_map(my_maps)
save_dir = "/Users/lahayes/QPP/interesting_event_2014-611/norh_plots/"


# bl_1 = [510, -220]
# tr_1 = [570, -170]
# width_1 = np.abs(bl_1[0] - tr_1[0])
# height_1 = np.abs(bl_1[1] - tr_1[1])

# bl_2 = [604, -225]
# tr_2 = [652, -175]
# width_2 = np.abs(bl_2[0] - tr_2[0])
# height_2 = np.abs(bl_2[1] - tr_2[1])

# submap_footpoint1 = make_submap(my_maps, bl_1[0], bl_1[1], tr_1[0], tr_1[1])
# footpoint_1_lcs = make_lcs_from_map(submap_footpoint1)
# submap_footpoint2 = make_submap(my_maps, bl_2[0], bl_2[1], tr_2[0], tr_2[1])
# footpoint_2_lcs = make_lcs_from_map(submap_footpoint2)


no_map.draw_contours(levels=np.arange(10, 100, 10)*u.percent, transfor
    m=ax.get_transform(no_map.wcs)) 


def ploty(mapy,fname='norh17', make_movie=True, **kwargs):

    for i in range(len(mapy)):
        print(i)
        fig = plt.figure(figsize=(10, 5))

        ax1 = fig.add_subplot(1,2,1, projection=mapy[i])

        mapy[i].plot(title='', **kwargs)
        mapy[i].draw_contours(levels=np.arange(1, 100, 10)*u.percent, cmap='BuPu_r')
        mapy[i].draw_rectangle(SkyCoord(bl_1[0]*u.arcsec, bl_1[1]*u.arcsec, frame=mapy[i].coordinate_frame), 
                                width_1*u.arcsec, height_1*u.arcsec, color='grey') 
        mapy[i].draw_rectangle(SkyCoord(bl_2[0]*u.arcsec, bl_2[1]*u.arcsec, frame=mapy[i].coordinate_frame), 
                                width_2*u.arcsec, height_2*u.arcsec, color='blue') 
        #mapy[i].draw_grid(grid_spacing=5*u.deg)

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(norh_lc, color='k', drawstyle='steps-mid', label='17GHz')
        ax2.plot(footpoint_1_lcs, color='grey', ls='dashed')
        ax2.plot(footpoint_2_lcs, color='blue', ls='dashed')
        ax2.set_xlim(norh_lc.index[0], norh_lc.index[-1])
        ax2.tick_params(axis='x', labelrotation=45)
        ax2.axvline(mapy[i].date.datetime)
        ax2.legend(loc='upper left')

        #plt.tight_layout()
        plt.savefig(save_dir + fname+'_{:04d}.png'.format(i))
        plt.close()

    if make_movie:
        subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
                        '1920x1080', '-i', '/Users/lahayes/QPP/interesting_event_2014-611/norh_plots/{:s}_%04d.png'.format(fname), 
                        '-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', '{:s}_mov.mp4'.format(fname)])




def make_submap(maps, bl_x, bl_y, tr_x, tr_y):
    """
    Function to make submap from bottom left, top tight
    x and y coords.

    """
    subs_rect = []
    for i in range(len(maps)):
        bll = SkyCoord(bl_x*u.arcsec, bl_y*u.arcsec, frame=maps[0].coordinate_frame)
        trr = SkyCoord(tr_x*u.arcsec, tr_y*u.arcsec, frame=maps[0].coordinate_frame)    
        subs_rect.append(maps[i].submap(bll, trr))
    return sunpy.map.Map(subs_rect, sequence=True)

