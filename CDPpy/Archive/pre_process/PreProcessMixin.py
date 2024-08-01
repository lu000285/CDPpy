import numpy as np
import pandas as pd

from CDPpy.Constants import ExpDataNamespace as Constants

class PreProcessMixn:
    '''
    Mixin class for MeasuredData class to execute pre-process.
    Calculate run time (day, hour) from meaured data.
    Calculate culture volume before/after sampling and after feeding from meaured data.

    Methods
    -------
        run_time : class method
            Calculate run time (day, hour) from Meaured Data.
        culture_volume : class method
            Calculate culture volume before/after sampling and after feeding from Meaured Data.
    '''
    def run_time(self):
        '''Calculate run time (day, hour) from Meaured Data.
        '''
        key = Constants()
        df = self.exp_data_df

        date = df[key.DATE]
        time = df[key.TIME]
        run_time_day = df[key.RUN_TIME_D]
        run_time_hour = df[key.RUN_TIME_H]

        if run_time_day.isnull().all() and run_time_hour.isnull().all():
            if pd.api.types.is_datetime64_any_dtype(date):
                date = date.apply(lambda x: x.date())
            
            if pd.api.types.is_datetime64_any_dtype(time):
                time = time.apply(lambda x: x.time())

            if (date==time).all():
                date_time = pd.to_datetime(date.astype(str))
            else:
                date_time = pd.to_datetime(date.astype(str) + ' ' + time.astype(str))
            
            time_diff = date_time - date_time.iat[0]
            run_time_hour = time_diff.dt.total_seconds() / 3600.0
            run_time_day = time_diff.dt.total_seconds() / 3600.0 / 24.0
            # run_time_day = time_diff.dt.days
        
            df[key.DATE] = date
            df[key.TIME] = time
            df[key.RUN_TIME_D] = run_time_day
            df[key.RUN_TIME_H] = run_time_hour


    def culture_volume(self):
        '''
        Calculate culture volume before/after sampling and after feeding from Meaured Data.
        '''
        key = Constants()
        df = self.exp_data_df
        feed_data = self.separate_feed_df.fillna(0)
        init_volume = self.initial_volume

        sample_size = df[key.SAMPLE].size
        sample_volume = df[key.SAMPLE_VOLUME].fillna(0).values
        base_added = df[key.BASE_ADDED].fillna(0).values
        feed_media_added = df[key.FEED_MEDIA].fillna(0).values
        v_after_sampling = df.pop(key.VOLUME_AFTER_SAMPLE)

        # Initialize
        v_before_sampling = np.zeros(sample_size)
        v_after_feeding = np.zeros(sample_size)

        if v_after_sampling.isna().all():
            v_after_sampling = np.zeros(sample_size)
            v_before_sampling[0] = init_volume

            # Added Supplements Volume; base + feed media + feed
            feed_sum = feed_data.sum(axis=1).values
            supplements_added = base_added + feed_media_added + feed_sum

            for i in range(sample_size):
                # Volume After Sampling
                v_after_sampling[i] = v_before_sampling[i] - sample_volume[i]
                            
                # Volume After Feeding
                v_after_feeding[i] = v_after_sampling[i] + supplements_added[i]
                
                # Volume Before Sampling
                if (i < sample_size-1):
                    v_before_sampling[i+1] = v_after_feeding[i]
        else:
            x = v_after_sampling.values
            v_before_sampling = x
            v_after_feeding = x

        df[key.VOLUME_BEFORE_SAMPLE] = v_before_sampling
        df[key.VOLUME_AFTER_FEED] = v_after_feeding
        df[key.VOLUME_AFTER_SAMPLE] = v_after_sampling