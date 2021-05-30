from IO import dataIO,SplineResampler
from kde1d import KDE1D
from scipy.optimize import minimize
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.stats import truncnorm
from processdata import ProcessVTK
from shapely.geometry import LineString, MultiPoint
from shapely.ops import split
from shapely import wkt
plt.rcParams.update({'font.size': 22})
from scipy.spatial.transform import Rotation
import pickle
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings("error")
from AttractedGrowthConeSearch import AttractedGrowthConeSearch

def get_macaque_kde(n=100):
	import pickle
	with open('original_main_trunk_kdes.pickle','rb') as f:
		return(pickle.load(f))

def set_axes_equal(ax):
	x_limits = ax.get_xlim3d()
	y_limits = ax.get_ylim3d()
	z_limits = ax.get_zlim3d()
	
	x_range = abs(x_limits[1] - x_limits[0])
	x_middle = np.mean(x_limits)
	y_range = abs(y_limits[1] - y_limits[0])
	y_middle = np.mean(y_limits)
	z_range = abs(z_limits[1] - z_limits[0])
	z_middle = np.mean(z_limits)
	
	# The plot bounding box is a sphere in the sense of the infinity
	# norm, hence I call half the max range the plot radius.
	plot_radius = 0.5*max([x_range, y_range, z_range])
	
	ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
	ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
	ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def plot_original_and_new(origs,news):
	fig = plt.figure()
	ax = Axes3D(fig)
	for orig in origs:
		xs,ys,zs = zip(*orig)
		ax.plot(xs,ys,zs,'k')
		ax.scatter(xs,ys,zs,color='r',s=2)
	
	for new in news:
		xs,ys,zs = zip(*new)
		ax.plot(xs,ys,zs,'r')	
		ax.scatter(xs,ys,zs,color='k',s=2)
	
	plt.title('Jittered (red) vs. Original (black)')
	set_axes_equal(ax)
	plt.show()

def load_finers(hemisphere):
	import pickle
	if hemisphere == 'left':
		with open('left_finer.pickle','rb') as f:
			finer = pickle.load(f)
	else:
		with open('right_finer.pickle','rb') as f:
			finer = pickle.load(f)
	
	return([np.array(list(zip(item[0],item[1],item[2]))) for item in finer])

def dump_newlines(nlines,hemisphere):
	import pickle
	with open(hemisphere+'_newfiners.pickle','wb') as f:
		pickle.dump(nlines,f)


if __name__ == "__main__":
	io = dataIO()
	kdemaker = KDE1D()
	seg_len_kde,meander_kde = get_macaque_kde()
	sampled_seg_lens = seg_len_kde.sample(10000)
	sampled_meanders = meander_kde.sample(10000)
	histology_meanders_xd,histology_meanders_logprob,kde = kdemaker.kde_1d(sampled_meanders)
	histology_seglen_xd,histology_seglen_logprob,kde = kdemaker.kde_1d(sampled_seg_lens)	
	streamlines = load_finers('right')
	streamlines = io.streamlines_to_tree(streamlines)
	sr = SplineResampler()
	streamlines = sr.resample_streamlines([stream[0] for stream in streamlines],seg_len_kde)
	newlines = []
	agc = AttractedGrowthConeSearch()
	for line in streamlines:
		newlines.append(agc.jitter_streamline(line,meander_kde,mode='backward'))
		plot_original_and_new([line],[list(newlines[-1])])
	
	dump_newlines(newlines,'right')
