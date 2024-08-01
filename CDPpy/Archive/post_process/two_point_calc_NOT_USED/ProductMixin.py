import pandas as pd
import numpy as np
###########################################################################
# IgG Two Point Calc Mixin Class
###########################################################################
class ProductTwoptMixn:
    '''
    '''
    # Call methods
    def post_process_twopt(self):
        self.igg_sp_rate_twopt()

        self._post_data_twopt = pd.concat([self._sp_rate],
                                           axis=1)
    # Calculate Specific Rate
    '''def igg_sp_rate_twopt(self):
        # Get Measurement Index
        idx = self._product_conc[self._product_conc.notnull()].index

        s = self._cumulative[idx]           # Substrate Concentration (mM)
        t = self._run_time_hour[idx]        # Run Time (hrs)
        v1 = self._v_before_sampling[idx]   # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling[idx]    # Culture Volume Before sampling (mL)
        xv = self._xv[idx]                  # Viable Cell Concentration (10e6 cells/mL)

        title = f'Two-Pt. Calc. q{self._name} (mmol/109 cell/hr)'
        rate = pd.Series(data=[np.nan] * len(self._sample_num),
                         name=title)

        for i in range(1, len(idx)):
            x = (s.iat[i] - s.iat[i-1]) * 1000
            y = xv.iat[i]*v1.iat[i] + xv.iat[i-1]*v2.iat[i-1]

            rate.iat[idx[i]] = x / (y*0.5*(t.iat[i]-t.iat[i-1]))

        self._sp_rate = rate'''

    def disp_post_data_twopt(self):
        print('\n************ IgG Post Process Data -Two Point Calc. ************')
        print(self._post_data_twopt)