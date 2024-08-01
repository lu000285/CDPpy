import pandas as pd
import numpy as np

###########################################################################
class OxygenTwoPtMixin:
    def post_process_twopt(self):
        self.sp_our()
        self.sp_oxy_cons_rate()

        self._post_data_twopt = pd.concat([self._sp_our,
                                           self._sp_rate],
                                           axis=1)

    # SP. Oxygen Uptake Rate
    '''def sp_our(self):
        # our: OUR in Measured Data
        self._sp_our = (self._our / self._xv).rename('SP. OUR (mmol/109 cells/hr)')
    '''

    # SP. Oxygen Consumption Rate
    '''def sp_oxy_cons_rate(self):
        t = self._run_time_hour                 # t: run time (hrs)
        v1 = self._v_before_sampling            # v1: culture volume before sampling (mL)
        v2 = self._v_after_sampling             # v2: culture volume after sampling (mL)
        rate = self._oxygen_consumption_rate    # rate: oxygen consumption rate in Measured Data
        c = self._oxygen_consumed               # c: oxyge consumed in Measured Data

        r = pd.Series(data=[np.nan]*len(t), name='SP. OXYGEN CONSUMPTION RATE (mmol/109cell/hr)')
        for i in range(1, len(t)):
            if rate.iat[i] < 0:
                x = c.iat[i]*v1.iat[i] - c.iat[i-1]*v2.iat[i-1]
                y = (self._vcc.iat[i]*v1.iat[i] + self._vcc.iat[i-1]*v1.iat[i-1])*1000000*0.5*(t.iat[i] - t.iat[i-1])
                r.iat[i] = x / 1000 / y *1000000000
            else:
                r.iat[i] = rate.iat[i]
        self._sp_rate = r'''

    # getters
    def get_sp_OUR(self):
        return self._sp_our

    def get_oxy_post_data(self):
        return self._post_data_twopt

    def disp_post_data_twopt(self):
        print('\n************ Oxygen Post Process Data -Two Point Calc. ************')
        print(self._post_data_twopt)
