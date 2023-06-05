import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

pi = math.pi
inf = math.inf
nan = math.nan

#ipmedance of the series LC subcircuit
def Impedance(f, lnC, lnL):
    match [math.isnan(lnC), math.isnan(lnL)]:
        case [True, True]:
            return 0.0
        case [True, False]:
            L = math.exp(lnL)
            return f*1j*2*pi*L
        case [False, True]:
            C = math.exp(lnC)
            return -1j*1/(f*2*pi*C)
        case _:
            C = math.exp(lnC)
            L = math.exp(lnL)
            ZC = -1j*1/(f*2*pi*C)
            ZL = 0+f*1j*2*pi*L
            return ZL*ZC/(ZL+ZC)

def Admittance(f, lnC, lnL):
    match [math.isnan(lnC), math.isnan(lnL)]:
        case [True, True]:
            return 0.0
        case [True, False]:
            L = math.exp(lnL)
            return 1/(f*1j*2*pi*L)
        case [False, True]:
            C = math.exp(lnC)
            return f*2*pi*C*1j
        case _:
            C = math.exp(lnC)
            L = math.exp(lnL)
            ZC = -1j*1/(f*2*pi*C)
            ZL = 0+f*1j*2*pi*L
            return (ZL+ZC)/ZL*ZC

class Spec:
    def __init__(self, starts, ends, limits, directions, margin, n):
        # starts of pass- or stopbands in Hz
        self.starts = starts
        # ends of start- or stopbands in Hz
        self.ends = ends
        # pass- or stopband threshold values
        self.limits = limits
        # pass- or stopbands? possible values: "pass" or "stop" strings
        self.directions = directions
        # number of frequency points
        self.n = n
        # margin to still punish solution that only barely satisfies the specification
        # on stopbands, with limit l and margin m, the margin range is from l*(1-m) to l
        # on passbands, with limit l and margin m, the margin range is from l to l*(1+m)
        self.margin = margin

    #def cost(self, types, values):
    def cost(self, plot=False, par1C=nan, par1L=nan, ser1C=nan, ser1L=nan, par2C=nan, par2L=nan, ser2C=nan, ser2L=nan, par3C=nan, par3L=nan, ser3C=nan, ser3L=nan):
        """Parameters define a ladder structure,
        where each series or parallel element consists of
        two discrete, ideal L or C components in parallel
        with each other. Default values define a perfect
        all-pass filter. In the function parameters, "parXY" means
        the value of the Y-th sub-component of the X-th parallel
        element. Similarly, "serXY" means the value of the Y-th
        sub-component of the X-th series element. The ladder starts
        with a parallel element."""
        values = [par1C, par1L, ser1C, ser1L, par2C, par2L, ser2C, ser2L, par3C, par3L, ser3C, ser3L]
        nval = len(values)
        nladder = int(nval/4)
        minf = (min(self.starts))
        maxf = (max(self.ends))
        # frequency axis sample points (log spacing)
        faxis = []
        # no. of frequency sample points
        #n = self.n
        factor = (maxf/minf)**((1/self.n))
        for i in range(self.n+1):
            faxis.append(minf*factor**(i))
        # reference impedance in Ohm, on both ports
        Z0 = 50
        Y0 = 1/Z0
        # array of overall S21 values at the frequency sample points
        S21 = []
        #process 1 parallel and 1 series element:
        for f in faxis:
            # ABCD parameter matrix of the whole system
            ABCD = np.matrix([[1,0],[0,1]])
            for l in range(nladder):
                # admittance of the two parallel components together from the current step of the ladder
                Y = Admittance(f, values[4*l+0], values[4*l+1])
                mPar = np.matrix([[1, 0],[Y, 1]])
                Z = Impedance(f, values[4*l+2], values[4*l+3])
                mSer = np.matrix([[1, Z],[0, 1]])
                ABCD = ABCD*mPar*mSer
            # 2/(A+B/Z0+C*Z0+D)
            A = ABCD.item(0,0)
            B = ABCD.item(0,1)
            C = ABCD.item(1,0)
            D = ABCD.item(1,1)
            S21.append(abs(2/(A+B/Z0+C*Z0+D)))
        # natural log of margin+1
        lnmargin=math.log(self.margin+1)
        if(plot):
            # green: passband; red: stopband; orange: ok, but close to not ok, still punished
            fig, ax = plt.subplots()
            minS21 = min(S21)
            for i in range(len(self.starts)):
                if(self.directions[i] == "pass"):
                    # region forbidden by the original specification
                    ax.add_patch(Rectangle((self.starts[i], minS21),
                                           self.ends[i]-self.starts[i],
                                           self.limits[i]-minS21,
                                           facecolor='#00aa00'))
                    # region close to original limit, but satisfying it, additional penalty region
                    ax.add_patch(Rectangle((self.starts[i], self.limits[i]),
                                           (self.ends[i]-self.starts[i]),
                                           self.limits[i]*self.margin,
                                           facecolor='orange'))
                else: # "stop"
                    # region forbidden by the original specification
                    ax.add_patch(Rectangle((self.starts[i], self.limits[i]),
                                           self.ends[i]-self.starts[i],
                                           1.0-self.limits[i],
                                           facecolor='#aa0000'))
                    # region close to original limit, but satisfying it, additional penalty region
                    ax.add_patch(Rectangle((self.starts[i], self.limits[i]*(1/(1+self.margin))),
                                           (self.ends[i]-self.starts[i]),
                                           self.limits[i]*(1-1/(1+self.margin)),
                                           facecolor='orange'))
            ax.loglog(faxis, S21, 'k-')
            plt.ylabel(r'$|S_{21}|$')
            plt.xlabel(r'$f$')
            plt.savefig("plot.pdf", dpi=120, format='pdf', bbox_inches='tight')
            plt.show()
        cost = 0
        # number of freq points in the regions where the S21 is specified
        ncost = 0
        for i in range(len(self.starts)):
            for j in range(len(faxis)):
                if faxis[j]>self.starts[i] and faxis[j]<self.ends[i]:
                    ncost = ncost + 1
                    lnlimit = math.log(self.limits[i])
                    lnS21 = math.log(S21[j])
                    if self.directions[i] == "pass":
                        cost += max(0, min(-lnS21, lnlimit+lnmargin-lnS21))
                    else: # "stop"
                        cost += max(0, lnS21-lnlimit+lnmargin)
		# negative of average cost, for function maximizing
        return -cost/ncost