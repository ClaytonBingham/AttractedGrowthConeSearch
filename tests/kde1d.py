import numpy as np
import scipy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
from sklearn.neighbors import KernelDensity as KDE
import os


class KDE1D():
	def __init__(self):
		pass
	
	def plot_kde(self,x_d,prob,key,xlims=None,output=os.getcwd()+'/',view=True):
		fig, ax = plt.subplots()
		ax.fill_between(x_d,prob,alpha=0.5,color='k')
		ax.set_ylim(np.min(prob)-np.min(prob)*0.025,np.max(prob)*1.2)
		if xlims is not None:
			ax.set_xlim(xlims[0],xlims[1])
		
		if key == 'Segment Lengths':
			ax.set_xlim(0.0,200.0)
		else:
			ax.set_xlim(0.0,220.0)
		# Turn on the minor TICKS, which are required for the minor GRID
		ax.minorticks_on()
		# Customize the major grid
		ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
		# Customize the minor grid
		ax.grid(which='minor', linestyle=':', linewidth='0.25', color='black')
		fig.set_size_inches(8, 6,forward=True)
		plt.title(key)
		if view:
			plt.show()
		else:
			fig.savefig(output+key+'.png',dpi=250,bbox_inches='tight')	
			plt.close('all')
	
	def plot_dual_kdes(self,x_d,prob,x_d2,prob2,key,xlims=None,output=os.getcwd()+'/',view=True):
		fig, ax = plt.subplots()
		ax.fill_between(x_d,prob,alpha=0.5,color='k',label='Histology')
		ax.set_ylim(np.min(prob)-np.min(prob)*0.025,np.max(prob)*1.2)
		ax.fill_between(x_d2,prob2,alpha=0.5,color='r',label='Model')
		if xlims is not None:
			ax.set_xlim(xlims[0],xlims[1])
		
		if key == 'Segment Lengths':
			ax.set_xlim(0.0,200.0)
		else:
			ax.set_xlim(0.0,220.0)
		# Turn on the minor TICKS, which are required for the minor GRID
		ax.minorticks_on()
		# Customize the major grid
		ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
		# Customize the minor grid
		ax.grid(which='minor', linestyle=':', linewidth='0.25', color='black')
		fig.set_size_inches(8, 6,forward=True)
		plt.legend()
		plt.title(key)
		if view:
			plt.show()
		else:
			fig.savefig(output+key+'.png',dpi=250,bbox_inches='tight')	
			plt.close('all')

	def plot_triple_kdes(self,x_d,prob,x_d2,prob2,x_d3,prob3,key,xlims=None,output=os.getcwd()+'/',view=False	):
		fig, ax = plt.subplots()
#		ax.fill_between(x_d,prob,alpha=0.5,color='k',label='Histology') 
		ax.plot(x_d,prob,color='k',linewidth=3.5,label='Histology')


		ax.set_ylim(np.min(prob)-np.min(prob)*0.025,np.max(prob)*1.2)
#		ax.fill_between(x_d2,prob2,alpha=0.5,color='#ad3e3e',label='Attracted Growth Cone Search')
		ax.plot(x_d2,prob2,color='#ad3e3e',label='Attracted Growth Cone Search',linewidth=3.5)
#		ax.fill_between(x_d3,prob3,alpha=0.5,color='#3295a8',label='Random Jitter')
		ax.plot(x_d3,prob3,color='#3295a8',label='Random Jitter',linewidth=3.5)
		if xlims is not None: 
			ax.set_xlim(xlims[0],xlims[1])
		
		if key == 'Segment Lengths':
			ax.set_xlim(0.0,200.0)
		else:
			ax.set_xlim(0.0,220.0)
		# Turn on the minor TICKS, which are required for the minor GRID
		ax.minorticks_on()
		# Customize the major grid
		ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
		# Customize the minor grid
		ax.grid(which='minor', linestyle=':', linewidth='0.25', color='black')
		fig.set_size_inches(10, 8,forward=True)
		if key == 'Segment Lengths':
			ax.set_xlabel('Segment Length ('+	u"\u03BC"+'m)')
			ax.set_ylabel('Probability Density')
			plt.legend()
		if key == 'Meanders':
			ax.set_xlabel('Meander Angles (degrees)')
			ax.set_ylabel('Probability Density')
			plt.legend()
#		plt.legend()
#		plt.title(key)
		if view:
			plt.show()
		else:
			fig = plt.gcf()
			fig.set_size_inches(10, 8,forward=True)
			fig.savefig(output+key+'.png',dpi=250,bbox_inches='tight')	
			plt.close('all')


	def cost_func(self,kde,comparison,key):
		nrmse = (((np.mean((comparison[1]-np.exp(kde.score_samples(comparison[0].reshape(-1,1))))**2))**0.5)/np.max(comparison[1]))*100.0
		print('NRMSE for '+str(key)+' = '+str(nrmse)+'%')

	def kde_1d(self,data,knl='gaussian'):
		if knl == 'epanechnikov':
			bdw = np.std(data)/1.5
		
		elif knl == 'gaussian':
	#		bdw = np.nanstd(data)/3.0
			a = np.nanmin([np.nanstd(data), np.subtract(*np.nanpercentile(data, [75, 25]))/1.34])
			bdw = 0.9*a*len(data)**(-1/5)
		
		if bdw == 0.0 or np.isnan(bdw):
			return(None,None)
		
		kde = KDE(bandwidth=bdw,kernel=knl)
		x = np.array(data)
		x = x[~np.isnan(x)]
		kde.fit(x[:,None])
		
		x_d = np.linspace(np.nanmin(data)-2*np.nanstd(data),np.nanmax(data)+2*np.nanstd(data),1000)
		logprob = kde.score_samples(x_d[:,None])
		return(x_d,np.exp(logprob),kde)
