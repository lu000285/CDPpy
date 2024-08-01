import numpy as np
import pandas as pd

from .Equation import production

class MetaboliteMixin:
    '''
    '''
    def in_process(self):
        '''
        '''
        s = self.cumulative_calc()
        self._cumulative = s
        
    def cumulative_calc(self) -> pd.DataFrame:
        '''Calculate cumulative consumption or production of a metabolite.
        Parameters
        ---------
            t : numpy array
                Time (hr)
            c : numpy array
                Concentration (mg/mL) or (mmol/L).
            d : numpy array
                Dillution rate at (hr^1).
                Flow rate (L/hr) / Culture volume (L)
            f : numpy array
                Feed concentration (mg/mL) or (mmol/L).
        Returns
            s : pandas DataFrame
                cumulative consumption or production (mg/mL) or (mmol/L).
        '''
        t = self._run_time['value'].values
        c = self._conc['value'].values
        d = self._dillution_rate['value'].values
        f = self._feed_conc['value'].values
        unit = self._conc['unit'].iat[0]

        s = np.zeros_like(t)
        for i in range(1, t.size):
            s_i = production(c[i-1], c[i], t[i-1], t[i], d[i], f[i])
            s[i] = s[i-1] + s_i

        # Check consumption or not.
        if (s<0).sum() > t.size/2:
            s *= -1
            state = 'consumed'
        else:
            state = 'produced'
        
        s = pd.DataFrame(data=s, columns=['value'])
        s['unit'] = unit
        s['state'] = state
        s.index.name = 'cumulatve concentration'
        return s
    
    @property
    def get_cumulative(self):
        ''''''
        return self._cumulative
        