import pandas as pd
import numpy as np

from .GetterMixin import GetterMixin

# Rolling Polynomial Regression Mixin Class
class MetaboliteMixin(GetterMixin):
    def rolling_window_polynomial(self, degree=3, windows=4):
        idx = self._idx
        x = self._run_time_hour.values[idx]
        y = self._cumulative.values[idx]
        vcc = self._xv.values[idx]
        v = self._v_before_sampling.values[idx]

        # x_mid = self._run_time_mid
        c_mid, x_mid = self.mid_point()
        
        self._roll_polyorder = degree
        self._roll_polywindow = windows
        self._c_mid = pd.Series(data=c_mid, name='mid_concentration_(mM)')
        self._t_mid = pd.Series(data=x_mid, name='mid_run_time_(hr)')
        
        # Calculate Specific Rate from Derivetive of Polynomial Regression on Cumularive Consumption/Production
        '''
        pre = f'Roll. Poly. Reg. Order: {degree} Window: {windows}'
        title = f'{pre} q{self._name.capitalize()} (mmol/109 cell/hr)'
        q = pd.Series(data=[pd.NA] * (len(x)-1), name=title)
        '''
        q = np.zeros(len(x)-1)
        q.fill(np.nan)

        # Calculate SP. rate
        for i in range(0, len(x_mid)):
            if (i + 1) < (windows / 2):
                x_roll = x[0: windows]
                y_roll = y[0: windows]
                # Polynomial Regression for Cumulative Consumption/Production
                fit = np.polyfit(x_roll, y_roll, degree)  # Fitting data to polynomial Regression (Get slopes)
                p = np.poly1d(fit)  # Get polynomial curve for Cumulative Consumption/Production

                dp1 = p.deriv()      # first derivetive of polynomial fit

                dy = dp1(x_mid[i])      # derivetive values corresponding to x
                if(self._measured_cumulative_flag):
                    q[i] = dy / (vcc[i] + vcc[i+1])*2
                else:
                    q[i] = dy / (vcc[i] * v[i]+vcc[i+1] * v[i+1]) * 2 * 1000

            elif (i + windows / 2) > len(x):
                x_roll = x[int(len(x)-windows/2-1):len(x)]
                y_roll = y[int(len(x)-windows/2-1):len(x)]
                # Polynomial Regression for Cumulative Consumption/Production
                fit = np.polyfit(x_roll, y_roll, degree)  # Fitting data to polynomial Regression (Get slopes)
                p = np.poly1d(fit)  # Get polynomial curve for Cumulative Consumption/Production

                dp1 = p.deriv()      # first derivetive of polynomial fit

                dy = dp1(x_mid[i])      # derivetive values corresponding to x
                if(self._measured_cumulative_flag):
                    q[i] = dy / (vcc[i] + vcc[i+1]) * 2
                else:
                    q[i] = dy / (vcc[i] * v[i]+vcc[i+1] * v[i+1]) * 2 * 1000
                
            else:
                x_roll = x[int(i-windows/2+1):int(i+windows/2+1)]
                y_roll = y[int(i-windows/2+1):int(i+windows/2+1)]
                # Polynomial Regression for Cumulative Consumption/Production
                fit = np.polyfit(x_roll, y_roll, degree)  # Fitting data to polynomial Regression (Get slopes)
                p = np.poly1d(fit)  # Get polynomial curve for Cumulative Consumption/Production

                dp1 = p.deriv()      # first derivetive of polynomial fit

                dy = dp1(x_mid[i])      # derivetive values corresponding to x
                if(self._measured_cumulative_flag):
                    q[i] = dy / (vcc[i] + vcc[i+1])*2
                else:
                    q[i] = dy / (vcc[i] * v[i] + vcc[i+1] * v[i+1]) * 2 * 1000

        # Rolling window polynomial regression
        name = f'q{self._name.capitalize()}_(mmol/10^9_cell/hr)'
        self._sp_rate_rollpoly = pd.Series(data=q, name=name)


    def mid_point(self):
        '''Calculate mid points of the run time and concentration.
        '''
        idx = self._idx                 # Measurement index
        c1 = self._c_after_feed.values[idx]      # c1: conc after feeding at t
        c2 = self._c_before_feed.values[idx]     # c2: measured conc at t + 1
        t = self._run_time_hour.values[idx]      # run time hour

        # c_mid = pd.Series(data=[na] * (len(idx)-1), name=f'{self._name} CONC. MID. (mM)', dtype='float')
        c_mid = np.zeros((len(idx)-1))
        t_mid = np.zeros((len(idx)-1))

        for i in range(len(idx)-1):
            c_mid[i] = 0.5 * (c1[i] + c2[i+1])
            t_mid[i] = 0.5 * (t[i] + t[i+1])

        return c_mid, t_mid

