import numpy as np
import math
import matplotlib.pyplot as plt

pi = math.pi
inf = math.inf

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

# calculate impedance of an L or C element on the given frequency
# if the input value is positive, then Z is calculated for a "value" valued inductance [H]
# if the input value is negative, then Z is calculated for a "-value" valued capacitance [F]
    def Impedance(freq, value)
        if(value >= 0): # value is inductance in [Henry]
            return freq*1j*2*pi*value
        else: # -1*value is capacitance in [Farad]
            return 1/freq*1j*2*pi*(-value)

# calculate parallel impedance of Z1 and Z2,
# considering possible 0 and infinite values
    def Parallel(Z1, Z2):
        if(Z1==0.0 or Z2==0.0):
            return 0.0
        elif(Z1==inf or Z2 == inf):
            if(Z1==inf):
                return Z2
            else:
                return Z1
        else:
            return Z1*Z2/(Z1+Z2)

	#def cost(self, types, values):
	def cost(self, par11=inf, par12=inf, ser11=0.0, ser12=0.0,
             par21=inf, par22=inf, ser21=0.0, ser22=0.0,
             par31=inf, par32=inf, ser31=0.0, ser32=0.0):
		"""Parameters define a ladder structure,
        where each series or parallel element consists of
        two discrete, ideal L or C components in parallel
        with each other. Default values define a perfect
        all-pass filter. In the function parameters, "parXY" means
        the value of the Y-th sub-component of the X-th parallel
        element. Similarly, "serXY" means the value of the Y-th
        sub-component of the X-th series element. The ladder starts
        with a parallel element."""
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
        #TODO: process arguments in a new way
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
		plt.plot(faxis, S21)
		plt.show()
		cost = 0
		# number of freq points in the regions where the S21 is specified
		ncost = 0
		for i in range(len(self.starts)):
			for j in range(len(faxis)):
				if faxis[j]>self.starts[i] and faxis[j]<self.ends[i]:
					ncost = ncost + 1
					if self.directions[i] == "pass":
						cost += max(0, 1-S21[j]/self.limits[i])
					else: # "stop"
						cost += max(0, S21[j]/self.limits[i]-1)
		return cost/ncost

spec = Spec([0.01], [1], [1], ["pass"])
print(spec.cost(['L', 'L', 'C'], [1, 0, 0.1]))
