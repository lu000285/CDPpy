import pandas as pd
import numpy as np

from CDPpy.in_process.Fed_batch.GetterMixin import GetterMixin

###########################################################################
# Oxygen Mixin Class
###########################################################################
class OxygenMixin(GetterMixin):
    '''
    
    '''
    # Call Calculation Function
    def in_process(self):
        # Calculate Cumulative Oxygen Concentration Consumed
        cc = self.cumulative_concentration
        self._cumulative = pd.Series(data=cc, name='cumOxygen_(mmol)')

        sp_our = self.sp_our
        self._sp_our = pd.Series(data=sp_our, name='sOUR_(mmol/10^9_cells/hr)')

        sp_rate = self.sp_rate
        self._sp_rate = pd.Series(data=sp_rate, name='qOxygen_(mmol/109cell/hr)')

    # Calculate Cumulative Oxygen Concentration Consumed
    @property
    def cumulative_concentration(self):
        idx = self._idx
        oxy = self._measured_consumption.values[idx]    # Concentration of Oxygen Consumed
        v1 = self._v_before_sampling.values[idx]    # Culture Volume Before Sampling
        v2 = self._v_after_sampling.values[idx]     # Culture Volume After Sampling

        # Initialize
        c = np.zeros(self._sample_num)
        c.fill(np.nan)
        c[0] = 0
        for i in range(1, len(idx)):
            c[i] = c[i-1] + (oxy[i] * v1[i] - oxy[i-1] * v2[i-1]) / 1000
        return c
    
    # Calculate SP. Oxygen Uptake Rate
    @property
    def sp_our(self):
        our = self._measured_our.values # Measured oxygen Up Take Rate
        xv = self._xv.values            # vialbe cell concentration (10e6 cells/ml)
        return (our / xv)
    
    # Calculate SP. Oxygen Consumption Rate
    @property
    def sp_rate(self):
        idx = self._idx
        t = self._run_time_hour.values[idx]                 # t: run time (hrs)
        v1 = self._v_before_sampling.values[idx]            # v1: culture volume before sampling (mL)
        v2 = self._v_after_sampling.values[idx]             # v2: culture volume after sampling (mL)
        rate = self._measured_consumption_rate.values    # rate: oxygen consumption rate in Measured Data
        c = self._measured_consumption.values[idx]             # c: oxygen consumption in Measured Data
        xv = self._xv.values[idx]            # vialbe cell concentration (10e6 cells/ml)

        r = np.zeros(self._sample_num)
        r.fill(np.nan)
        for i in range(1, len(idx)):
            if rate[i] < 0:
                x = c[i] * v1[i] - c[i-1] * v2[i-1]
                y = (xv[i] * v1[i] + xv[i-1] * v1[i-1]) * 1000000 * 0.5 * (t[i] - t[i-1])
                r[idx[i]] = x / 1000 / y * 1000000000
            else:
                r[i] = rate[i]
        return r
    
    def get_sp_OUR(self):
        return self._sp_our

    @property
    def get_in_process(self):
        """
        Get In-Process DataFrame.
        """
        if self._in_process_flag:
            t = self._run_time_hour
            cc = self._cumulative
            sp_our = self._sp_our
            spr = self._sp_rate
        return pd.concat([t, cc, sp_our, spr], axis=1)
    
    @property
    def get_in_process_data(self):
        '''Return In-Process data.
        '''
        if self._in_process_flag:
            measured_c = self._measured_consumption
            coc = self._cumulative
            measured_r = self._measured_consumption_rate
            spr = self._sp_rate
            measured_our = self._measured_our
            sp_our = self._sp_our
            data = [measured_c,
                    coc,
                    measured_r, spr,
                    measured_our, sp_our]
            profile = ['concentration',
                    'cumulative',
                    'spRate', 'spRate',
                    'spRate', 'spRate']
            kind = [np.nan,
                    np.nan,
                    'rate', 'rate',
                    'our', 'our',]
            method = [np.nan,
                    'twoPoint',
                    'measurement', 'twoPoint',
                    'measurement', 'twoPoint']
            return self.get_profile_data(data_list=data,
                                            profile_list=profile,
                                            kind_list=kind, 
                                            method_list=method)
    # Display
    @property
    def disp_in_process(self):
        if self._in_process_flag:
            data = self.get_in_process
            print('\n************ Oxygen In Process Data ************')
            print(data)


