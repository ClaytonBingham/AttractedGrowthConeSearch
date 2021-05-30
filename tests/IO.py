from kde1d import KDE1D
from processdata import ProcessVTK
import roots
from roots.swcToolkit import swcToolkit
from scipy.optimize import minimize
import time
import numpy as np
from shapely.geometry import LineString, MultiPoint
from shapely.ops import split
from shapely import wkt
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
import pickle

class dataIO():
	def __init__(self):
		pass


	def get_streamlines(self,fname=None):
		if not fname:
			print('dataIO.get_streamlines() method requires a filename (fname) argument')
			return([])
		
		vtk_loader = ProcessVTK()
		lines = vtk_loader.vtk_2_streamlines(fname) #'m1_hdp_upperexCF.vtk'
		return(lines)


	def streamlines_to_tree(self,streamlines):
		cffs = []
		for streamline in streamlines:
			cffs.append({})
			cffs[-1][0] = np.array(streamline,dtype=float)*10**3
		
		return(cffs)


	def trees_to_swcs(self,trees,odir):
		swctool = swcToolkit()
		for t,tree in enumerate(trees):
			swctool.to_swc(tree,target = odir+str(t)+'.swc')

class SplineResampler():
	def __init__(self):
		pass
	
	def spline_branch(self,points,target_intersegment_length=14.7061928165059):
		number_of_points = int(len(points)*self.average_segmental_dist(points)/target_intersegment_length)+1
		line = LineString([tuple(point) for point in points])
		splitter = MultiPoint([line.interpolate((float(i)/number_of_points),normalized=True) for i in range(number_of_points+1)])
		return(self.parse_spline(splitter.wkt))

	def parse_spline(self,shapelystring):
		return(np.array([np.array(list(point)) for point in np.array(wkt.loads(shapelystring))]))

	def average_segmental_dist(self,points):
		return(np.linalg.norm(np.mean(points[1:]-points[:-1],axis=0)))

	def sample_kde(self,sampler,constraints=(0,100)):
		sample = constraints[0]-1
		while sample < constraints[0] or sample > constraints[1]:
			sample = sampler.sample(1)[0][0]
		
		return(sample)

	def calculate_path_length(self,points):
		length = 0
		for p,point in enumerate(points[:-1]):
			length+=np.linalg.norm(points[p+1]-points[p])
		
		return(length)

	def find_target_path_length_in_points(self,points,target_path_length):
		for i in range(len(points)):
			if self.calculate_path_length(points[:i]) > target_path_length:
				return(points[i-1],i)
			
		return(points[0],1)

	def resample_streamline(self,stream,seg_len_kde):
		splined = self.spline_branch(stream,1)
		newstream = [splined[0]]
		while len(splined) > 1:
			seglen = 0
			target_seglen = self.sample_kde(seg_len_kde)
			if target_seglen > 6:
				target_seglen-=5
			
			if target_seglen < 0:
				target_seglen = 1
			
			endpoint,endindex = self.find_target_path_length_in_points(splined,target_seglen)			
			newstream.append(endpoint)
			splined = splined[endindex:]
		return(newstream)

	def resample_streamlines(self,streams,seg_len_kde):
		newstreams = []
		for s,stream in enumerate(streams):
			newstreams.append(self.resample_streamline(stream,seg_len_kde))
			print('resampled '+str(s+1)+' in '+str(len(streams))+' streamlines')
		
		with open('subsampled_streamlines.pickle','wb') as f:
			pickle.dump(newstreams,f)
		
		return(newstreams)

	def load_old_streamlines(self,fname='subsampled_streamlines.pickle'):
		import pickle
		with open(fname,'rb') as f:
			dat = pickle.load(f)
		
		return([{0:d} for d in dat])	



