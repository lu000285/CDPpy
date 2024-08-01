import pandas as pd
import numpy as np
from ..Fed_batch.GetterMixin import GetterMixin

class ProductMixin(GetterMixin):
    '''Product Poly. Reg. Mixin Class
    '''
    def polynomial(self, deg=3):
        # Fitting a polynomial using the calculated cumulative concentration.
        self._polyorder = deg
        poly = self.polyfit()
        self._polyfit = poly    # will be used

        # Polynomial cumulative concentration
        c = self.cumulative_polynomial()
        name = f'cumProduct_(mg)'
        self._cumulative_poly = pd.Series(data=c, name=name, dtype=float)

        # Polynomial specific rate
        spr = self.sp_rate_polynomial()
        self._sp_rate_poly = pd.Series(data=spr, name='qProduct(mg/10^9_cell/hr)', dtype=float)


    def polyfit(self):
        '''Fit a polynomial of degree deg using the calculated cumulative concentration.'''
        idx = self._idx                     # Measurement Index
        t = self._run_time_hour.values[idx] # Run Time (hrs)
        s = self._cumulative.values[idx]    # Cumulative concentration (mM)
        deg = self._polyorder

        # Fitting data to polynomial Regression (Get slopes)
        fit = np.polyfit(x=t, y=s, deg=deg)
        # Get a polynomial for cumulative Concentration.
        return np.poly1d(fit)
    
    def cumulative_polynomial(self):
        '''Return the polynomial cumulative concentration corresponding to run time t.
        '''
        # idx = self._idx            # Measurement Index
        t = self._run_time_hour.values #[idx]    # Run Time (hrs)
        return self._polyfit(t)
    
    def sp_rate_polynomial(self):
        '''Return the polynomial specific rate corresponding to run time t.
        '''
        idx = self._idx
        t = self._run_time_hour.values #[idx]     # Run Time (hrs)
        xv = self._xv.values #[idx]               # Viable Cell Concentration (10e6 cells/mL)
        v = self._v_before_sampling.values #[idx] # Culture Volume Before Sampling (mL)

        DpDt = self._polyfit.deriv()    # the first derivetive of the polynomial fit
        y = DpDt(t)         # The derivetive corresponding to run time t.

        # Initialize
        r = np.zeros(self._sample_num)
        r.fill(np.nan)

        # Calculate SP. rate
        for i in range(1, len(t)):
            r[i] = y[i] / (xv[i] * v[i]) * 1000
        return r

