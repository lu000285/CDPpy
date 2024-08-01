import numpy as np
import pandas as pd
from ..Fed_batch.GetterMixin import GetterMixin

class MetaboliteMixin(GetterMixin):
    '''Polynomial Regression Mixin Class
    '''
    def polynomial(self, deg=3):
        '''
        '''
        # Fitting a polynomial using the calculated cumulative concentration.
        self._polyorder = deg
        poly = self.polyfit()
        self._polyfit = poly    # will be used for the plot

        # Get a polynomial for cumulative concentration
        c = self.cumulative_polynomial()
        name = f'cum{self._name.capitalize()}_(mmol)'
        self._cumulative_poly = pd.Series(data=c, name=name, dtype=float)

        # Get Polyreg SP. Rate
        r = self.sp_rate_polynomial()
        name = f'q{self._name.capitalize()}_(mmol/10^9_cell/hr)'
        self._sp_rate_poly = pd.Series(data=r, name=name, dtype=float)

    # Fitting cumulative to Polynomial Regression
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
            # If there is the measured cumulative concentration.
            if (self._measured_cumulative_flag):
                r[i] = y[i] / xv[i]
            else:
                r[i] = y[i] / (xv[i] * v[i]) * 1000
        return r

    def disp_polyreg(self):
        df = pd.concat([self._polyreg_cumulative,
                        self._polyreg_sp_rate],
                        axis=1)
        print('\n************ Post Process Data -Poly. Reg. ************')
        print(df)