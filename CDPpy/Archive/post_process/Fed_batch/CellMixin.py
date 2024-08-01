import pandas as pd
import numpy as np
from .GetterMixin import GetterMixin

class CellMixin(GetterMixin):
    '''Cell Poly. Reg. Mixin Class.
    '''
    def polynomial(self, deg=3):
        # Fitting a polynomial using the calculated cumulative concentration.
        self._polyorder = deg
        poly = self.polyfit()
        self._polyfit = poly    # will be used

        # Polynomial cumulative concentration
        c = self.cumulative_polynomial()
        name = f'cumCell_(x10^6_cells)'
        self._cumulative_poly = pd.Series(data=c, name=name, dtype=float)

        # Polynomial specific rate
        spr = self.sp_rate_polynomial()
        self._sp_rate_poly = pd.Series(data=spr, name='qCell_(1/hr)', dtype=float)


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
        t = self._run_time_hour.values[idx]         # run time (hrs)
        xv = self._xv.values[idx]                   # vialbe cell concentration (10e6 cells/ml)
        s = self._cumulative_poly.values[idx]       # Cumulative Cell Concentraion (10e6 cells/mL)
        v1 = self._v_before_sampling.values[idx]    # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling.values[idx]     # Culture Volume After Sampling (mL)

        # Initialize
        r = np.zeros(self._sample_num)
        r.fill(np.nan)
        for i in range(1, len(t)):
            x = s[i] - s[i-1]
            y = xv[i] * v1[i] + xv[i-1] * v2[i-1]
            r[i] = x / (y * 0.5 * (t[i] - t[i-1]))
        return r

    def disp_polyreg(self):
        df = pd.concat([self._polyreg_cumulative],
                        axis=1)
        print('\n************ Post Process Data -Poly. Reg. ************')
        print(df)
