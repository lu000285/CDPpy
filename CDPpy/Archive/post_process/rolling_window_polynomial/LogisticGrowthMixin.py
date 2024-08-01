import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Logistic Growth Mixin Class
class LogisticGrowthMixin:
    '''
    '''
    def LogisticGrowthFit(self):
        tcc = self._xt[0:11]
        t = self._run_time_hour[0:11]
        N0 = tcc[0]
        popt, pcov = curve_fit(lambda t,K,r:LogisticGrowth_fun(t,N0,K,r), t, tcc, p0=[24,0.05])
        return popt

    def midcalc_growth_rate_calc(self):
        df = pd.DataFrame()
        t = self.mid_point()
        self._t_mid = t
        t = pd.Series(data=t, name='mid_run_time_(hr)')
        popt = self.LogisticGrowthFit()
        N0 = self._xt[0]
        mu = pd.Series([np.nan] * len(t))
        mu = (t.apply(lambda x:sp_rate_calc(x,N0,*popt))).rename('mu')

        df['sp_growth_rate_(1/hr)'] = mu
        df['mu_max'] = pd.Series([popt[1]] * len(t))
        df['N0'] = pd.Series([N0] * len(t))
        df['K_mu'] = pd.Series([popt[0]] * len(t))
        # print(f'N0: {N0}')
        # print(f'popt: {popt}')
        # print('mu:')
        # print(mu)
        self._growth_rate = df
        # return df
    
    def mid_point(self):
        '''Calculate mid points of the run time and concentration.
        '''
        idx = self.measurement_index     # Measurement index
        t_hour = self.run_time_hour[idx]      # run time hour
        t_day = self.run_time_day[idx]      # run time hour

        t_day_mid = np.zeros((len(idx)-1))
        t_hour_mid = np.zeros((len(idx)-1))

        for i in range(len(idx)-1):
            t_day_mid[i] = 0.5 * (t_day[i] + t_day[i+1])
            t_hour_mid[i] = 0.5 * (t_hour[i] + t_hour[i+1])

        return t_day_mid, t_hour_mid

# Logistic Growth Function
def LogisticGrowth_fun(t,N0,K,r):
    return K*N0/(N0+(K-N0)*np.exp(-r*t))

# Specific rate calculation: mu = dx_t/dt/x_t
def sp_rate_calc(t,N0,K,r):
    return (K-N0)*r*np.exp(-r*t)/(N0+(K-N0)*np.exp(-r*t))



