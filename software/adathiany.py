#import numpy as np

class Spec:
	def __init__(self, starts, ends, limits, directions):
		# starts of pass- or stopbands
		self.starts = starts
		# ends of start- or stopbands
		self.ends = ends
		# pass- or stopband threshold values
		self.limits = limits
		# pass- or stopbands?
		self.directions = directions

	def cost(self, types, values):
		"""types:
				'L': inductance
				'C': capacitance
			values in SI units [Farad/Henry]"""
		minf = (min(self.starts))
		maxf = (max(self.ends))
		faxis = []
		n = 100
		factor = (maxf/minf)**((1/n))
		for i in range(n+1):
			faxis.append(minf*factor**(i))
		for i in range(n+1):
			for j in range()


spec = Spec([10, 20, 50], [12, 22, 52], [1,1,1], ["pass", "stop", "pass"])
spec.cost(['L', 'C', 'L'], [1, 1, 1])
