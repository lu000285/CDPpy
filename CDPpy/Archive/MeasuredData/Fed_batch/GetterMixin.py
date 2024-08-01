import pandas as pd

class GetterMixin:
    '''
    Mixin class for MeasuredData Class.

    Methods
    -------
        get_pre_data : 
            Return pre process data.
    '''
    def get_pre_data(self):
        """Return pre process data.
        """
        df = self.param_df
        run_time_day = df['run_time_(days)']
        run_time_hour = df['run_time_(hrs)']
        v_before_sampling = df['volume_before_sampling_(mL)']
        v_after_sampling = df['volume_after_sampling_(mL)']
        
        self.pre_data = pd.concat([run_time_day, run_time_hour, v_before_sampling, v_after_sampling], axis=1)
        return self.pre_data