import pandas as pd
import numpy as np

class GetterMixin:
    '''Getter Mixin Class for Post-Process: Polynomial regression.
    '''    
    @property
    def get_post_process(self):
        """Get Post-Process DataFrame.
        """
        t = self._run_time_hour
        cc = self._cumulative_poly
        spr = self._sp_rate_poly
        deg = self._polyorder
        data = pd.concat([t, cc, spr], axis=1)
        data['polynomialDegree'] = deg
        return data
    
    @property
    def get_post_process_data(self):
        """
        Get Post-Process DataFrame.
        """
        # Create a cumulative concetration data for a polynomial plot.
        t = self._run_time_hour.values
        t_min, t_max = t[0], t[-1]
        t_poly = np.linspace(start=t_min, stop=t_max, num=100)
        if self._polyfit is None:
            cumulative_poly = np.zeros_like(t_poly)
            cumulative_poly.fill(np.nan)
        else:
            cumulative_poly = self._polyfit(t_poly)
        df1 = pd.DataFrame(data={'runTime': t_poly, 
                                 'value': cumulative_poly})
        df1['profile'] = 'cumulative'
        
        # specific rate data
        df2 = pd.DataFrame(data={'runTime': self._run_time_hour.values, 
                                 'value': self._sp_rate_poly.values})
        df2['profile'] = 'spRate'

        df = pd.concat([df1, df2], axis=0)
        df['kind'] = np.nan
        df['method'] = 'polynomial'
        return df