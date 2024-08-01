import pandas as pd
import numpy as np
###########################################################################
# Cell Two Point Calc Mixin Class

class CellMixnTwoPt:
    '''
    '''
    # Call methods
    def post_process_twopt(self):
        # self.sp_growth_rate()
        #self.kd()
        #self.mv()

        '''self._post_data_twopt = pd.concat([self._sp_growth_rate,
                                           self._mv,
                                           self._kd],
                                           axis=1)'''
    
    # Calculates Specific growth rate
    '''def sp_growth_rate(self):
        t = self._run_time_hour         # run time (hrs)
        s = self._cumulative            # Cumulative Cell Concentraion (10e6 cells/mL)
        v1 = self._v_before_sampling    # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling     # Culture Volume After Sampling (mL)

        # Initialize
        rate = pd.Series(data=[np.nan] * len(t), name='SP. GROWTH RATE, m (hr-1) [mv-kd]') # Initialize
        
        for i in range(1, len(t)):
            x = s.iat[i] - s.iat[i-1]
            y = self._xv.iat[i]*v1.iat[i] + self._xv.iat[i-1]*v2.iat[i-1]
            rate.iat[i] = x / (y*0.5*(t.iat[i]-t.iat[i-1]))

        # SP. Rate
        self._sp_growth_rate = rate
        self._sp_rate = rate'''

    # Calculates kd value
    '''def kd(self):
        xd = self._xd                   # Dead Cell Concentration (10e6 cells/mL)
        t = self._run_time_hour         # Run Time (hrs)
        v1 = self._v_before_sampling    # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling     # Culture Volume After Sampling (mL)

        # Initialize
        kd = pd.Series(data=[np.nan] * len(t), name='kd')

        for i in range(1, len(t)):
            x = xd.iat[i]*v1.iat[i] - xd.iat[i-1]*v2.iat[i-1]
            y = self._xv.iat[i]*v1.iat[i] + self._xv.iat[i-1]*v2.iat[i-1]
            kd.iat[i] = x / (y*0.5*(t.iat[i]-t.iat[i-1]))
        self._kd = kd

    # Calculate mv
    def mv(self):
        self._mv = (self._sp_growth_rate + self._kd).rename('mv')'''

    # getters
    '''def get_kd(self):
        return self._kd

    def get_sp_growth_rate(self):
        return self._sp_growth_rate

    def get_mv(self):
        return self._mv

    def get_post_data_twopt(self):
        return self._post_data_twopt

    def disp_post_data_twopt(self):
        print('\n************ Cell Post Process Data -Two Point Calc. ************')
        print(self._post_data_twopt)'''

###########################################################################