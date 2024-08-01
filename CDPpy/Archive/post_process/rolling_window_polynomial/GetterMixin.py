import pandas as pd
import numpy as np

class GetterMixin:
    '''Getter Mixin Class for Post-Process: Rolling Polynomial regression.
    '''
    @property
    def get_post_process_roll(self):
        """Get Post-Process DataFrame.
        """
        t = self._t_mid
        spr = self._sp_rate_rollpoly
        deg = self._roll_polyorder
        windows = self._roll_polywindow = windows

        data = pd.concat([t, spr], axis=1)
        data['rollingPolynomialDegree'] = deg
        data['rollingPolynomialWindow'] = windows
        return data
    
    @property
    def get_post_process_roll_data(self):
        """
        Get Post-Process DataFrame.
        """
        t = self._t_mid.values
        spr = self._sp_rate_rollpoly
        
        # specific rate data
        df = pd.DataFrame(data={'runTime': t, 'value': spr})
        df['profile'] = 'spRate'
        df['kind'] = np.nan
        df['method'] = 'rollingWindowPolynomial'
        return df