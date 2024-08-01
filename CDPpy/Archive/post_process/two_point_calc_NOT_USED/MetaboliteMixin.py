import pandas as pd
import numpy as np
from numpy import diff

###########################################################################
class MetaboliteMixinTwoPt:
    '''
    '''
    # Calculate Specific Rate
    def sp_rate_twopt(self):
        # Get Measurement Index
        idx = self._idx
        s = self._cumulative[idx]           # Cumulative Concentration (mM)
        t = self._run_time_hour[idx]        # Run Time (hrs)
        v1 = self._v_before_sampling[idx]   # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling[idx]    # Culture Volume Before sampling (mL)
        xv = self._xv[idx]                  # Viable Cell Concentration (10e6 cells/mL)
        
        # Initialize
        rate = pd.Series(data=[np.nan] * len(self._sample_num),
                         name='q'+self._name+' (mmol/109 cell/hr)')
        
        # IF Have Direct Calculation of Cumulative
        if (self._direct_cumulative):
            # print(f'{self._name} use direct CUM for Two-Pt Calc.')
            '''dsdt = diff(s)/diff(t)
            dsdt = np.insert(dsdt, 0, np.nan)
            rate = dsdt / xv'''
            for i in range(1, len(idx)):                
                SC = (s.iat[i] - s.iat[i-1])
                Xv = 0.5*(xv.iat[i] + xv.iat[i-1])
                rate.iat[idx[i]] = (1 / Xv) * SC / (t.iat[i]-t.iat[i-1])

        else:
            for i in range(1, len(idx)):                
                SC = (s.iat[i] - s.iat[i-1]) * 1000

                # With Concentration after Feed
                if (self._use_conc_after_feed):
                    XvV = 0.5*(xv.iat[i]*v1.iat[i] + xv.iat[i-1]*v1.iat[i]) # Average
                    
                # With Feed Concentration OR Without Both
                else:
                    XvV = 0.5*(xv.iat[i]*v1.iat[i] + xv.iat[i-1]*v2.iat[i-1]) # Average

                rate.iat[idx[i]] = (1 / XvV) * SC / (t.iat[i]-t.iat[i-1])

        # SP. Rate
        self._sp_rate = rate

    ###########################################################################