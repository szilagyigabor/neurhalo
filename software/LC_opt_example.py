from bayes_opt import BayesianOptimization
from LC_opt_defs import Spec
import math

# the starts of the specified bands in Hz
starts = [0.1, 100, 10000]
# the ends of the specified bands in Hz
ends = [7, 1000, 1000000]
# maximum or minimum values of |S21| between 0 and 1 (Smax and Smin)
limits = [0.1, 0.7, 0.3]
# specify type of bands (pass- or stopbands)
directions = ["stop", "pass", "stop"]
# declare Specification
# margin is the "height" of the margin: (Smax*(1+margin))/Smax or Smin/(Smin*(1+margin))
# n is the number of frequencies for evaluating the h(f,x) error function
spec = Spec(starts, ends, limits, directions, margin=0.47, n=1000)

# specifying which L/C elements are used and their max/min values in H/F
bounds = {'par1C': (math.log(1e-7), math.log(1e-2)), 'ser1C': (math.log(1e-7), math.log(1e-2))}

# declare optimizer, specify target function (spec.cost) and bounds
# verbosity can be set to 1 or 0 to adjust textual output
# set random state for exact reproducibility between runs
optimizer = BayesianOptimization(
    f=spec.cost,
    pbounds=bounds,
    verbose=2,
    random_state=1,
)

# run optimization with init_points random samples at the beginning
# and n_iter more iterations with active learning
optimizer.maximize(
    init_points=4,
    n_iter=10,
)

# show best parameter values
print(optimizer.max)

# print transfer function for the best results
# WARNING: parameter names should be adjusted according to what is used above for it to work
spec.cost(plot=True, par1C=optimizer.max['params']['par1C'], ser1C=optimizer.max['params']['ser1C'])