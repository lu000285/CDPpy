import pandas as pd
import numpy as np

from .GetterMixin import GetterMixin
from CDPpy.Constants import MetaboliteNameSpace, ExpDataNamespace
from CDPpy.Constants import SPECIES_STATE

CONSTANTS = MetaboliteNameSpace()
KEY = ExpDataNamespace()

class MetaboliteMixin(GetterMixin):
    '''
    '''
    # Call Initialize Method
    def in_process(self, use_feed_conc, use_conc_after_feed):
        ''' Calculate In-Proccess for metabolites.
        '''
        self._use_feed_conc = use_feed_conc
        self._use_conc_after_feed = use_conc_after_feed

        if not use_conc_after_feed:
            self._conc_after_feed = self.conc_after_feeding_calc()

        # If measured data already has the calculated cumulative consumption/production
        if (self.measured_cumulative_flag):
            # self._conc_after_feed = self.conc_after_feeding_calc()
            s = self.measured_cumulative_conc
        else:
            # IF Experiments measure the feed Concentrations
            if (use_feed_conc):
                # self._conc_after_feed = self.conc_after_feeding_calc()
                # Calculate Cumulative Consumption/Production with Feed Concentraion
                s = self.cumulative_calc_by_feed()
            # IF Experiments measure the concentraions after feeding
            elif (use_conc_after_feed):
                # Calculate Cumulative Consumption/Production with Concentraion after Feeding
                s = self.cumulative_calc_by_conc_after_feed()
            else:
                # Calculate Concentration After Feeding
                # self._conc_after_feed = self.conc_after_feeding_calc()
                s = self.cumulative_calc_without_feed()

        run_time = self.run_time
            
        # Cumulative Consumption/Production
        self._cumulative_conc = pd.concat([run_time, s], axis=1)

        # Specific rate
        r = self.sp_rate_calc()
        self._sp_rate = pd.concat([run_time, r], axis=1)

        # Concat concentrations before and after feeding
        c1 = pd.concat([run_time, self.conc_before_feed], axis=1)
        c2 = pd.concat([run_time, self.conc_after_feed], axis=1)
        conc = pd.concat([c1, c2], axis=0)
        conc.sort_index(inplace=True)
        self._conc = conc.reset_index(drop=True)

    def conc_after_feeding_calc(self) -> pd.DataFrame:
        ''' Calculate concentration after feeding.
        '''
        idx = self.measurement_index
        s = self.conc_before_feed['value'].values[idx]# Substrate Concentration (mM)
        sf = self.feed_conc['value'].values[idx]# Substrate Feed Concentration (mM)
        f_media = self.feed_media_added['value'].values[idx]# Feed Flowrate (ml/hr)
        v = self.volume_after_sampling['value'].values[idx]# Culture Volume After Sampling (ml)
        separate_feed = self.separate_feed# Separate Feed DataFrame
        separate_f_added = self.separate_feed_sum # Separate Feed Sum Added (E.g. glutamine, glucose, etc.)

        # Check if species has the separate feed
        if separate_feed is None:
            f = f_media
        else:
            f = separate_feed['value'].values[idx]

        # Concentration After Feeding
        c = ((s*v + sf*f) / (v + f_media + separate_f_added))
        c = pd.DataFrame(data=c,columns=['value'])
        c['unit'] = CONSTANTS.CONCE_UNIT
        return c

    def cumulative_calc_by_feed(self) -> pd.DataFrame:
        ''' Calculate ccumulative consumption using the feed concentration.
        '''
        idx = self.measurement_index
        s = self.conc_before_feed['value'].values[idx]     # Substrate Concentration (mM)
        sf = self.feed_conc['value'].values[idx]           # Substrate Feed Concentration (mM)
        feed_media = self.feed_media_added['value'].values[idx]     # Feed Flowrate (ml/hr)
        v1 = self.volume_before_sampling['value'].values[idx]   # Culture Volume Before Sampling (ml)
        v2 = self.volume_after_sampling['value'].values[idx]    # Culture Volume After Sampling (ml)
        separate_feed = self.separate_feed        # Separate Feed data

        # Check if species has the separate feed
        if separate_feed is None:
            f = feed_media
        else:
            f = separate_feed['value'].values[idx]

        # Initialize
        c = np.zeros(self.samples)
        c.fill(np.nan)
        c[0] = 0

        # Consumed Substrate = sf*f + s[i-1]*v[i-1] - s[i]*v[i]
        for i in range(1, len(idx)):
            c_i = (sf[i] * f[i-1] - s[i] * v1[i] + s[i-1] * v2[i-1]) / 1000
            c[idx[i]] = c[idx[i-1]] + c_i
        if SPECIES_STATE[self.name] == 'Produced':
            c *= -1
        c = pd.DataFrame(data=c, columns=['value'])
        c['unit'] = CONSTANTS.CUMULATIVE_UNIT
        c['state'] = SPECIES_STATE[self.name]
        c['method'] = 'twoPoint'
        return c

    def cumulative_calc_by_conc_after_feed(self) -> pd.DataFrame:
        ''' Calculate cumulative consumption using the concentration after feeding
        '''
        idx = self.measurement_index
        s1 = self.conc_before_feed['value'].values[idx]    # Substrate Concentration Before Feeding (mM)
        s2 = self.conc_after_feed['value'].values[idx]     # Substrate Concentration After Feeding (mM)
        v = self.volume_before_sampling['value'].values[idx]    # Culture Volume After Feeding (ml)

        # Initialize
        c = np.zeros(self.samples)
        c.fill(np.nan)
        c[0] = 0

        for i in range(1, len(idx)):
            c_i = (s2[i-1] * v[i] - s1[i] * v[i]) / 1000
            c[idx[i]] = c[idx[i-1]] + c_i
        if SPECIES_STATE[self.name] == 'Produced':
            c *= -1
        c = pd.DataFrame(data=c, columns=['value'])
        c['unit'] = CONSTANTS.CUMULATIVE_UNIT
        c['state'] = SPECIES_STATE[self.name]
        c['method'] = 'twoPoint'
        return c

    def cumulative_calc_without_feed(self) -> pd.DataFrame:
        ''' Calculate cumulative consumption without the concentration after feeding.
        '''
        idx = self.measurement_index
        s = self.conc_before_feed['value'].values[idx]     # Substrate Concentration (mM)
        v1 = self.volume_before_sampling['value'].values[idx]   # Culture Volume Before Sampling (ml)
        v2 = self.volume_after_sampling['value'].values[idx]    # Culture Volume After Sampling (ml)

        # Initialize
        c = np.zeros(self.samples)
        c.fill(np.nan)
        c[0] = 0

        # Consumed Substrate = sf*f + s[i-1]*v[i-1] - s[i]*v[i]
        for i in range(1, len(idx)):
            c_i = (- s[i] * v1[i] + s[i-1] * v2[i-1]) / 1000
            c[idx[i]] = c[idx[i-1]] + c_i
        if SPECIES_STATE[self.name] == 'Produced':
            c *= -1
        c = pd.DataFrame(data=c, columns=['value'])
        c['unit'] = CONSTANTS.CUMULATIVE_UNIT
        c['state'] = SPECIES_STATE[self.name]
        c['method'] = 'twoPoint'
        return c
    
    # Calculate Specific Rate
    def sp_rate_calc(self) -> pd.DataFrame:
        ''' Calculate specific rate.
        '''
        idx = self.measurement_index
        s = self.cumulative_conc['value'].values[idx]           # Cumulative Concentration (mM)
        t = self.run_time_hour['value'].values[idx]        # Run Time (hrs)
        v1 = self.volume_before_sampling['value'].values[idx]   # Culture Volume Before Sampling (mL)
        v2 = self.volume_after_sampling['value'].values[idx]    # Culture Volume Before sampling (mL)
        xv = self.viable_cell_conc['value'].values[idx]  # Viable Cell Concentration (10e6 cells/mL)
        
        # Initialize
        r = np.zeros(self.samples)
        r.fill(np.nan)
        
        # If there is the measurements of the cumulative concentration.
        if (self.measured_cumulative_flag):
            for i in range(1, len(idx)):                
                c_diff = (s[i] - s[i-1])            # concentration difference
                xv_avg = 0.5 * (xv[i] + xv[i-1])    # average vcc
                r[idx[i]] = c_diff / (t[i] - t[i-1]) / xv_avg
        else:
            for i in range(1, len(idx)):                
                c_diff = (s[i] - s[i-1]) * 1000 # concentration difference
                # With the concentration after feed
                if (self.use_conc_after_feed and not self.use_feed_conc):
                    xv_avg = 0.5 * (xv[i] * v1[i] + xv[i-1] * v1[i]) # Average
                # With the feed concentration OR without both
                else:
                    xv_avg = 0.5 * (xv[i] * v1[i] + xv[i-1] * v2[i-1]) # Average
                r[idx[i]] = c_diff / (t[i]-t[i-1]) / xv_avg
        r = pd.DataFrame(data=r, columns=['value'])
        r['unit'] = CONSTANTS.SP_RATE_UNIT
        r['method'] = 'twoPoint'
        r.index.name = f'q{self.name}'
        return r

    @property
    def use_feed_conc(self):
        return self._use_feed_conc
    @property
    def use_conc_after_feed(self):
        return self._use_conc_after_feed
    
    def get_in_process_data(self):
        '''Return in-process data.
        '''
        if self._in_process_flag:
            c_before = self._c_before_feed
            c_after_feed = self._c_after_feed
            feed_c =  self._feed_c
            cc = self._cumulative
            rate = self._sp_rate

            data = [c_before, c_after_feed, feed_c,
                    cc, rate]
            profile = ['concentration', 'concentration', 'concentration',
                       'cumulative', 'spRate']
            kind = ['beforeFeed', 'afterFeed', 'feed',
                    'cumulative', 'rate']
            method = [np.nan, np.nan, np.nan,
                      'twoPoint', 'twoPoint']
            return self.get_profile_data(data_list=data,
                                         profile_list=profile,
                                         kind_list=kind, 
                                         method_list=method)

    # Display
    def disp_in_process(self):
        if self._in_process_flag:
            data = pd.concat([self._run_time_hour, self._cumulative], axis=1)
            print(data)
        
