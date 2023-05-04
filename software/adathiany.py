import numpy as np
import math
import matplotlib.pyplot as plt

pi = math.pi

class Spec:
	def __init__(self, starts, ends, limits, directions):
		# starts of pass- or stopbands in Hz
		self.starts = starts
		# ends of start- or stopbands in Hz
		self.ends = ends
		# pass- or stopband threshold values
		self.limits = limits
		# pass- or stopbands? possible values: "pass" or "stop" strings
		self.directions = directions

	def cost(self, types, values):
		"""types:
				'L': inductance
				'C': capacitance
			values in SI units [Farad/Henry]
			0 value series L or parallel C means
			skipping that element.
			Starts with parallel (shunt) element."""
		minf = (min(self.starts))
		maxf = (max(self.ends))
		# frequency axis sample points (log spacing)
		faxis = []
		# no. of frequency sample points
		n = 100
		factor = (maxf/minf)**((1/n))
		for i in range(n+1):
			faxis.append(minf*factor**(i))
		# reference impedance in Ohm
		Z0 = 50
		# array of overall S21 values at the
		# frequency sample points
		S21 = []
		for i in range(n+1):
			# S-matrix of the current element
			S = np.matrix([[1,0],[0,1]])
			ss = "shunt" # start with shunt component
			for j in range(len(values)):
				if types[j] == 'L':
					Z = faxis[i]*1j*2*pi*values[j]
				else: # 'C'
					Z = 1/(faxis[i]*1j*2*pi*values[j])
				if ss == "shunt":
					Y = 1/Z
					Y0 = 1/Z0
					m = np.matrix([[-Y, 2*Y0],
						[2*Y0, -Y]])/(Y+2*Y0)
					ss = "series"
				else: # "series"
					m = np.matrix([[-Y, 2*Y0],
						[2*Y0, -Y]])/(Y+2*Y0)
					ss = "shunt"
				S=S*m
			S21.append(abs(S.item(0,1)))
		#plt.plot(faxis, S21)
		#plt.show()
		cost=0
		for i in len(self.starts):
			for j in len(faxis):
				if faxis[j]>self.starts[i] && faxis[j]<self.ends[i]:
					if self.directions[i] == "pass"
						cost += max(0, S21[j]/self.limits[i])
					else: # "stop"

spec = Spec([0.01], [1], [1], ["pass"])
spec.cost(['L', 'L', 'C'], [1, 0, 0.1])
