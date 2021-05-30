import numpy as np
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings("error")

class AttractedGrowthConeSearch():
	def __init__(self):
		pass

	def cone(self,a,target_angle=20.0,cone_height=None): #normal vs. weibull for forward looking cone vs. backward looking cone
		if cone_height is None:
			a = np.linalg.norm(a[1]-a[0])
		else:
			a = cone_height
		
		b = np.arctan(np.radians(target_angle))*a*2
		h = a * (np.random.uniform(0.0,1,1000) ** (1/3))
		r = (b / a) * h * np.sqrt(1)
		
		t = 2 * np.pi * np.random.uniform(0,1,1000)
		return(np.array(list(zip(r*np.cos(t),h,r*np.sin(t)))))

		

	def shift_and_rotate_points_from_matrix(self,origin, points, R):
		"""

		Rotate a point counterclockwise by a given angle around a given origin.
		The angle should be given in radians.

		Shift cone origin to origin of vector
			
		"""
		r = Rotation.from_matrix(R)
		rotated_points = r.apply(points)
		return(rotated_points+origin)	

	def rotation_matrix_from_vectors(self,vec2):
		""" Find the rotation matrix that aligns vec1 to vec2
		:param vec1: A 3d "source" vector
		:param vec2: A 3d "destination" vector
		:return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
		"""
		vec1 = np.array([0,np.linalg.norm(vec2[1]-vec2[0]),0])
		vec2 = vec2[1]-vec2[0]
		a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
		v = np.cross(a, b)
		c = np.dot(a, b)
		s = np.linalg.norm(v)
		kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
		rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
		return(rotation_matrix)


	def calculate_angle(self,b,a,c):
		ba = a - b
		bc = c - b
		num = np.dot(ba,bc)
		den = np.linalg.norm(ba) * np.linalg.norm(bc)
		if den == 0.0:
			print('returned zero den')
			return(0.0)
		
		cosine_angle_vecs = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		try:
			angle = np.arccos(cosine_angle_vecs)
		except:
			print(cosine_angle_vecs,'couldnt arccos this angle')
			angle = 0.0
		return(np.degrees(angle))

	def return_normalized_dists(self,a,cone):
		dists = cdist([a],cone)
		return(dists[0]/np.max(dists[0]))

	def return_by_angle_error(self,a,b,cone_pnts,target_angle):
		angles = np.array([self.calculate_angle(a,b,cone_point) for cone_point in cone_pnts])
		errors = np.abs(target_angle-angles)
		return(errors/np.max(errors))
		
	def find_lowest_ranked_cone_point(self,a,b,c,cone_pnts,target_angle,direction='forward'):
		dists = self.return_normalized_dists(c,cone_pnts)
		angles = self.return_by_angle_error(a,b,cone_pnts,target_angle)
		if direction=='forward':
			return(cone_pnts[np.argmin(angles)])
		else:
			return(cone_pnts[np.argmin(dists)])

	def apply_cone_search_forward(self,p0,p1,p2,target_angle):
		cone_pnts = self.cone([p1,p2],target_angle)
		R = self.rotation_matrix_from_vectors([p1,p2])
		rotated_points = self.shift_and_rotate_points_from_matrix(p1,cone_pnts,R)
		return(rotated_points)	
		
	def apply_cone_search_backward(self,p0,p1,p2,target_angle):
		cone_pnts = self.cone([p0,p1],target_angle,cone_height=np.linalg.norm(p2-p1))
		R = self.rotation_matrix_from_vectors([p0,p1])
		rotated_points = self.shift_and_rotate_points_from_matrix(p1,cone_pnts,R)
		return(rotated_points)

	def jitter_streamline(self,line,angle_kde,mode='backward'):
		target_angles = None
		while target_angles is None or len(target_angles)<len(line):
			target_angles = [angle for angle in angle_kde.sample(int(len(line)*1.2)) if angle < 80][:len(line)]
		
		newpoints = [line[0],line[1]]
		for p,point in enumerate(line[:-2]):
			if mode == 'forward':
				cpnts = self.apply_cone_search_forward(newpoints[p],newpoints[p+1],line[p+2],target_angles[p])
			else:
				cpnts = self.apply_cone_search_backward(newpoints[p],newpoints[p+1],line[p+2],target_angles[p])
			
			newpoints.append(self.find_lowest_ranked_cone_point(newpoints[p],newpoints[p+1],line[p+2],cpnts,target_angle=target_angles[p],direction=mode))	
		
		return(np.array([list(item) for item in newpoints]))


if __name__ == "__main__":
	pass
