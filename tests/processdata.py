import numpy as np
import subprocess
import os

class ProcessVTK():
	def __init__(self):
		pass

	def surfaceVTK_to_3d_mesh(self,fname):
		stemname = fname.split('.')[0]
		ddir = os.getcwd()+'/'
		subprocess.call('gmsh -2 '+fname+' -o '+stemname+'.stl',shell=True)		
		subprocess.call('./TetWild/build/TetWild --input '+stemname+'.stl --output '+stemname+'.mesh',shell=True)

	def vtk_2_streamlines(self,fname):
		dat,pnts,lines = self.read_vtk(fname)			
		pnt_cloud = self.vtk_pntblock_2_xyz(pnts)
		xyz_lines = self.vtk_lineblock_2_xyz(pnt_cloud,lines)
		return(xyz_lines)
	
	def vtk_2_pnt_cloud(self,fname):
		dat,pnts,lines = self.read_vtk(fname)	
		pnt_cloud = self.vtk_pntblock_2_xyz(pnts)
		return(pnt_cloud)
		

	def vtk_pnt_line_2_array(self,pntline):
		ar = [float(item) for item in pntline.strip().split(' ')]
		return(np.array(ar).reshape(-1,3))

	def vtk_pntblock_2_xyz(self,block):
		pnts = self.vtk_pnt_line_2_array(block[0])
		for line in block[1:]:
			pnts = np.append(pnts,self.vtk_pnt_line_2_array(line),axis=0)
		
		return(pnts)


	def vtk_lineblock_2_xyz(self,pnts,block):
		lines_by_index = [[int(item) for item in line.strip().split(' ')[1:]] for line in block]
		lines = [None for i in lines_by_index]
		for l,line in enumerate(lines_by_index):
			lines[l] = np.array([pnts[ind] for ind in line])
		
		return(lines)

	def read_vtk(self,fname):
		with open(fname,'r') as f:
			dat = f.readlines()
		
		for line in dat:
			if 'POINTS ' in line:
				pnt_start = dat.index(line)
						
			if 'LINES ' in line:
				ln_start = dat.index(line)
		
		pnts = dat[pnt_start+1:ln_start]
		lines = dat[ln_start+1:-1]
				

		return(dat,pnts,lines)


class ProcessMSH():
	def __init__(self):
		pass


	def float_node(self,node_string):
		return([float(item) for item in node_string])

	def msh_to_points(self,fname):
		with open(fname,'r') as f:
			dat = f.readlines()
		
		for line in dat:
			if '$Nodes' in line:
				node_start = dat.index(line)+2
			
			if '$EndNodes' in line:
				node_end = dat.index(line)
		
		nodes = np.array([self.float_node(node.strip().split(' ')[1:]) for node in dat[node_start:node_end]]).flatten().reshape(-1,3)
		
		return(dat,nodes)


class CreateCloud():
	def __init__(self):
		pass

	def load_pnts_simps(self,fname):
		with open(fname,'rb') as f:
			raw = f.readlines()
		
		for line in raw:
			if 'Vertices' in line.decode('utf-8'):
				start_pnts = raw.index(line)+2
			if 'Triangles' in line.decode('utf-8'):
				end_pnts = raw.index(line)
			if 'Tetrahedra' in line.decode('utf-8'):
				start_smps = raw.index(line)+2
			if 'End' in line.decode('utf-8'):
				end_smps = raw.index(line)
		
		pnts = [d.decode('utf-8').split(' ') for d in raw[start_pnts:end_pnts]]
		for d,da in enumerate(pnts):
			pnts[d] = [float(item.strip()) for item in da if item != ''][:-1]
		
		pnts = dict(zip(range(len(pnts)),np.array(pnts)))
		
		smps = [d.decode('utf-8').split(' ') for d in raw[start_smps:end_smps]]
		for d,da in enumerate(smps):
			smps[d] = [int(item.strip()) for item in da if item != ''][:-1]
		
		return(pnts,smps)

	def rand_in_tet(self,verts):
		s = np.random.random()
		t = np.random.random()
		u = np.random.random()
		if s+t > 1.0:
			s = 1.0-s
			t = 1.0-t
		if t+u > 1.0:
			tmp = u
			u = 1.0-s-t
			t = 1.0-tmp
		elif s+t+u > 1.0:
			tmp = u
			u = s+t+u-1.0
			s = 1.0-t-tmp
		a = 1.0-s-t-u
		return(verts[0]*a+verts[1]*s+verts[2]*t+verts[3]*u)

	def get_targets_from_msh(self,fname):	
		pnts,smps = self.load_pnts_simps(fname)
		random_pnts = []
		for smp in smps:
			if self.threshold_sample(0.25):
				verts = [pnts[s-1] for s in smp]
				random_pnts.append(list(self.rand_in_tet(verts)))
		
		return(random_pnts)

	def threshold_sample(self,threshold=0):
		val = np.random.random()
		if val>threshold:
			return(True)
		else:
			return(False)
		
	def combine_clouds(self,a,b):
		return(np.append(a,b,axis=0))
