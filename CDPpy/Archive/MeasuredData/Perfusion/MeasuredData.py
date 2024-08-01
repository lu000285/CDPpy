import pandas as pd

from CDPpy.helper import split_df, input_path

TARGET_COLUMNS = [
    'Experiment Information',
    'Experimental Data', 
    'Concentration',
    'Feed Concentration'
    ]

class MeasuredData():
    '''
    '''
    def __init__(self, file_name) -> None:
        df = pd.read_excel(input_path(file_name=file_name))
        df_list = split_df(df, TARGET_COLUMNS)
        exp_df = df_list[0]
        param_df = df_list[1]
        conc_df = df_list[2]
        feed_df = df_list[3]

        self._exp_df = exp_df
        self._param_df = param_df
        self._conc_df = conc_df
        self._feed_df = feed_df

    @property
    def exp_df(self):
        return self._exp_df
    
    @property
    def param_df(self):
        return self._param_df
    
    @property
    def conc_df(self):
        return self._conc_df
    
    @property
    def feed_df(self):
        return self._feed_df