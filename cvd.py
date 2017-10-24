""" Solves the Callendar--Van Dusen (CVD) equation to approximate the resistance-temperature curve associated with a resistance-temperature device (RTD).

@author: Christian T. Jacobs
"""

# Copyright (c) 2017 Christian T. Jacobs

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from sympy import *


class CallendarVanDusen(object):

    def __init__(self, standard="DIN43760", subzero=False):
        """ Set up the coefficients and symbols required in the Callendar--Van Dusen equation.
        
        :arg str standard: The standard of the RTD curve.
        :arg bool subzero: True if subzero temperatures are being considered, False otherwise.
        """
    
        self.r = Symbol("r")  # Resistance
        self.t = Symbol("t")  # Temperature

        # Coefficients from "Measuring Temperature with RTDs - A Tutorial" by National Instruments Corporation: https://newton.ex.ac.uk/teaching/CDHW/Sensors/an046.pdf
        if(standard == "DIN43760"):
            # DIN 43760 standard.
            self.A = 3.9080e-3
            self.B = -5.8019e-7
            if(subzero):
                self.C = -4.2735e-12
            else:
                self.C = 0.0
        elif(standard == "American"):
            # American standard.
            self.A = 3.9692e-3
            self.B = -5.8495e-7
            if(subzero):
                self.C = -4.2325e-12
            else:
                self.C = 0.0
        elif(standard == "ITS-90"):
            # ITS-90 standard.
            self.A = 3.9848e-3
            self.B = -5.870e-7
            if(subzero):
                self.C = -4.0000e-12
            else:
                self.C = 0.0
        else:
            raise ValueError("Unknown standard %s" % standard)

    def temperature(self, r0, r_measured):
        """ Return the temperature at a measured resistance r_measured.
        
        :arg float r0: The resistance at temperature t = 0.
        :arg float r_measured: The measured resistance.
        :rtype: float
        :returns: The temperature from the CVD equation.
        """
        cvd = Eq(self.r, r0*(1.0 + self.A*self.t + self.B*self.t**2 + self.C*(self.t-100.0)**3))
        s = solve(cvd, self.t)
        t = s[1].subs(self.r, r_measured)
        return t

cvd = CallendarVanDusen()
r_measured=200.0
t = cvd.temperature(r0=500.0, r_measured=r_measured)
print("Resistance = %g" % (r_measured))
print("Temperature = %g" % (t))
