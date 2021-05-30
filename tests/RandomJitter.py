from kde1d import KDE1D
from processdata import ProcessVTK
import numpy as np
from shapely.geometry import LineString, MultiPoint
from shapely.ops import split
from shapely import wkt

class RandomJitter():
	def __init__(self):
		pass
	
	def random_jitter_vector(self,jitter_distance):
		vec = np.array([np.random.uniform(-1,1) for i in range(3)])
		return(vec*(jitter_distance/np.linalg.norm(vec)))
	
	def jitter_point(self,midpoint,jitter_distance):
		return(midpoint+self.random_jitter_vector(jitter_distance))						
	
	def jitter_streamline(self,points,jitter_distance):
		return(np.array([self.jitter_point(point,jitter_distance) for point in points]).reshape(-1,3))

	def jitter_streamline_by_segments(self,points,jitter_distances,target_angles,error_threshold=0.75):
		newstream = np.array([])
		for p,point in enumerate(points):
			if p <2:
				newstream = np.append(newstream,point)
				continue
						
			if p>2 or p in range(len(points))[::2]:
				newpoint = self.jitter_point(point,jitter_distances[p-1])
				try:
					angle_error = np.abs(np.abs(target_angles[p-1])-np.abs(self.calculate_angle(np.array(newstream[p-1]),np.array(newstream[p-2]),np.array(newpoint))))
					print(newstream[p-1],newstream[p-2],newpoint,'points')
				except:
					print(target_angles[p-1],self.calculate_angle(newstream[p-1],newstream[p-2],newpoint),'angles')
										
				newstream = np.append(newstream,newpoint)
				continue
				
			newstream = np.append(newstream,point)

		return(newstream.reshape(-1,3))
	
	def calc_CD(self,A,B,target_theta):
		C = (A+B)/2.0
		AC = np.linalg.norm(C-A)
		AD = AC/np.cos(np.radians(target_theta))
		CD = (AD**2-AC**2)**0.5
		return(CD)
	
	def average_segmental_dist(self,points):
		return(np.linalg.norm(np.mean(np.array(points[1:])-np.array(points[:-1]),axis=0)))
	
	def calculate_target_jitter_dist(self,points,target_theta):
		return(self.calc_CD(np.array([0,0,0]),self.average_segmental_dist(points),target_theta))
	
	def calculate_target_jitter_dists(self,points,target_thetas):
		return(np.array([self.calc_CD(points[p],points[p+1],target_thetas[p]) for p in range(len(points))[:-1]]))

	def calculate_angle(self,a,b,c):
		ba = a - b
		bc = c - b
		num = np.dot(ba,bc)
		den = np.linalg.norm(ba) * np.linalg.norm(bc)
		if den == 0.0:
			return(0.0)
		
		cosine_angle_vecs = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		try:
			angle = np.arccos(cosine_angle_vecs)
		except:
			print(cosine_angle_vecs,'couldnt arccos this angle')
			angle = 0.0
		return(np.degrees(angle))
	


class EstimateRandomJitter():
	def __init__(self):
		self.rj = RandomJitter()
		self.kdemaker = KDE1D()
	
	def score_two_kdes(self,A,kdeB):
		'''
		2d normalized root mean squared error of two kernel density estimates
		
		Args: 
			xd: numpy array
				x-values of kernel to be scored
			logprob: numpy array
				y-values of kernel to be scored
				
		Returns: float 
			2d NRMSE
		
		'''
		xd,logprob = A
		comparison = np.exp(kdeB.score_samples(xd[:,None]))
		residuals = comparison-logprob
		return(((np.mean(residuals**2))**0.5)/np.abs(np.max(comparison)-np.min(comparison)))
	
	def calculate_angle(self,a,b,c):
		ba = a - b
		bc = c - b
		num = np.dot(ba,bc)
		den = np.linalg.norm(ba) * np.linalg.norm(bc)
		if den == 0.0:
			return(0.0)
		
		cosine_angle_vecs = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		try:
			angle = np.arccos(cosine_angle_vecs)
		except:
			print(cosine_angle_vecs,'couldnt arccos this angle')
			angle = 0.0
		return(np.degrees(angle))
	
	def calculate_meander_angles(self,points):
		meanders = []
		for p in range(len(points[:-2])):
			meanders.append(self.calculate_angle(points[p+1],points[p],points[p+2]))
		
		return(meanders)
	
	def get_meanders_kernel(self,points):
		meanders = self.calculate_meander_angles(points)
		xd,logprob,kde = self.kdemaker.kde_1d(meanders)
		return(xd,logprob)
			
	def evaluate_newstream(self,newstream,target_meander_kde,tries,error_threshold=0.2):
		error = self.score_two_kdes(self.get_meanders_kernel(newstream),target_meander_kde)
		if error > error_threshold:
			return(newstream,False)
		else:
			print('With '+str(tries)+'% target jitter distance forgiveness, '+str(round(error*100.0,2))+'% meander angle error was achieved')
		return(newstream,True)
	
	def execute_jitter(self,target_intersegment_length,target_meander_kde,streamline):
		streamline = self.rj.spline_branch(streamline,target_intersegment_length)
		jitter_distance = self.rj.calculate_target_jitter_dist(streamline,target_theta=10.0)
		tries = 0
		newstream,acceptable = self.evaluate_newstream(self.rj.jitter_streamline(streamline,jitter_distance),target_meander_kde,tries,error_threshold=1.10)		
		while not acceptable and tries < 100:
			newstream,acceptable = self.evaluate_newstream(self.rj.jitter_streamline(streamline,jitter_distance),target_meander_kde,tries,error_threshold=1.10)
			tries+=1
		
		
		return(newstream,self.get_meanders_kernel(newstream),np.linalg.norm(newstream[1:]-newstream[:-1],axis=1))

	def execute_jitter_by_segments(self,target_intersegment_length,target_meander_kde,streamline):
		target_angles = []
		while len(target_angles)<len(streamline):
			target_angle = target_meander_kde.sample(1)[0]
			while target_angle>80:
				target_angle = target_meander_kde.sample(1)[0]
			target_angles.append(target_angle)
		
		
		jitter_distances = self.rj.calculate_target_jitter_dists(streamline,target_angles)
		newstream = self.rj.jitter_streamline_by_segments(streamline,jitter_distances,target_angles,error_threshold=25.75)		
		return(newstream,self.get_meanders_kernel(newstream),np.linalg.norm(newstream[1:]-newstream[:-1],axis=1))
	
	def jitter_streamlines(self,target_intersegment_length,target_meander_kde,streams):
		meanders = np.array([])
		newstreams = []
		segment_lengths = np.array([])
		for s,stream in enumerate(streams):
			newstream,meander,segment_length = self.execute_jitter(target_intersegment_length,target_meander_kde,streams[s][0])
			newstreams.append(newstream)
			meanders = np.append(meanders,meander)
			segment_lengths = np.append(segment_lengths,segment_length)
		
		return(newstreams,meanders,segment_lengths)
	
	def jitter_streamlines_by_segments(self,target_intersegment_length,target_meander_kde,streams):
		meanders = np.array([])
		newstreams = []
		segment_lengths = np.array([])
		for s,stream in enumerate(streams):
			newstream,meander,segment_length = self.execute_jitter_by_segments(target_intersegment_length,target_meander_kde,streams[s][0])
			newstreams.append(newstream)
			meanders = np.append(meanders,meander)
			segment_lengths = np.append(segment_lengths,segment_length)
			print('finished streamline '+str(s)+' of '+str(len(streams)))
		
		return(newstreams,meanders,segment_lengths)
	
	def timit(self,f,*argv):
		start_time = time.time()
		if argv == None:
			return(0)
		
		else:
			first_guess = argv[0]
			params = argv[1:]
			for i in range(5):
				res = minimize(f,first_guess,args=tuple(params),method='Powell')
		
		print((time.time()-start_time)/60.0,f)
