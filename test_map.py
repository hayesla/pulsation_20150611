import sunpy.map
import mapcube_tools
import glob 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize
import subprocess

aia_files = glob.glob('/Users/lahayes/QPP/interesting_event_2014-611/ssw_fits/*131_.fts')
aia_files.sort()

aia_maps = sunpy.map.Map(aia_files)
aia_maps2 = []
for m in aia_maps:
    if m.exposure_time.value > 1:
        aia_maps2.append(m)

aia_maps2 = sunpy.map.Map(aia_maps2, sequence=True)

dif_maps2 = mapcube_tools.running_difference(aia_maps2)

save_dir = '/Users/lahayes/QPP/interesting_event_2014-611/aia_diff/'
def ploty(maps, fname='dif_131', make_movie=True, **kwargs):
    for i in range(len(maps)):
        print(i)
        maps[i].plot(**kwargs)

        plt.colorbar()
        plt.savefig(save_dir + fname+'_{:04d}.png'.format(i))
        plt.close()

    if make_movie:
        subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
                        '1920x1080', '-i', '/Users/lahayes/QPP/interesting_event_2014-611/aia_diff/{:s}_%04d.png'.format(fname), 
                        '-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', '{:s}_mov.mp4'.format(fname)])
