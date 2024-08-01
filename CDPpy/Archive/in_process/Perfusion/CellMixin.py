import numpy as np
import pandas as pd

from .Equation import production, death_rate, growth_rate

class CellMixin:
    '''
    '''
    def in_process(self, a, c):
        '''
        '''
        s = self.cumulative_calc()
        self._cumulative = s

        g_rate, d_rate = self.sp_rate_calc(a, c)
        self._sp_rate = g_rate
        self._sp_death_rate = d_rate
        
    def cumulative_calc(self) -> pd.DataFrame:
        '''Calculate cumulative production of the cell.
        Parameters
        ---------
            t : numpy array
                Time (hr)
            c : numpy array
                Viable cell concentration (10^6 cells/mL).
            d : numpy array
                Dillution rate at (hr^1).
                Flow rate (L/hr) / Culture volume (L)
        Returns
            s : pandas DataFrame
                cumulative consumption or production (10^6 cells/mL).
        '''
        t = self._run_time['value'].values
        c = self._viable_cell_conc['value'].values
        d = self._dillution_rate['value'].values
        unit = self._viable_cell_conc['unit'].iat[0]

        s = np.zeros_like(t)
        for i in range(1, t.size):
            s_i = production(c[i-1], c[i], t[i-1], t[i], d[i], 0)
            s[i] = s[i-1] + s_i
        
        s = pd.DataFrame(data=s, columns=['value'])
        s['unit'] = unit
        s['state'] = 'produced'
        s.index.name = 'cumulatve concentration'
        return s
    
    def sp_rate_calc(self, a, c) -> tuple[pd.DataFrame, pd.DataFrame]:
        '''Calculate the specific growth rate and death rate of the cell.
        Parameters
        ---------
            a : float
                a recycling factor
            c : float
                a concentration factor
            t : numpy array
                Time (hr)
            xd : numpy array
                Dead cell concentration (10^6 cells/mL).
            xv : numpy array
                Viable cell concentration (10^6 cells/mL).
            d : numpy array
                Dillution rate at (hr^1).
                Flow rate (L/hr) / Culture volume (L)
        Returns
            g_rate, d_rate : tupple of pandas DataFrames
                Specific growth rate (hr^1) and death rate (hr^1).
        '''
        t = self._run_time['value'].values
        xd = self._dead_cell_conc['value'].values
        xv = self._viable_cell_conc['value'].values
        d = self._dillution_rate['value'].values

        dr, gr = np.zeros_like(t), np.zeros_like(t)
        dr.fill(np.nan)
        gr.fill(np.nan)
        for i in range(1, t.size):
            dr[i] = death_rate(a, c, xd[i-1], xd[i], xv[i-1], xv[i], t[i-1], t[i], d[i-1], d[i])
            gr[i] = growth_rate(a, c, xv[i-1], xv[i], t[i-1], t[i], d[i-1], d[i], dr[i])
        
        d_rate = pd.DataFrame(data=dr, columns=['value'])
        d_rate['unit'] = '(hr^-1)'
        d_rate['state'] = 'death'
        d_rate.index.name = 'Specific death rate'

        g_rate = pd.DataFrame(data=gr, columns=['value'])
        g_rate['unit'] = '(hr^-1)'
        g_rate['state'] = 'growth'
        g_rate.index.name = 'Specific growth rate'
        return g_rate, d_rate
    
    @property
    def get_cumulative(self):
        ''''''
        return self._cumulative
    
    @property
    def get_sp_rate(self):
        ''''''
        return self._sp_rate
    
    @property
    def get_sp_death_rate(self):
        ''''''
        return self._sp_death_rate