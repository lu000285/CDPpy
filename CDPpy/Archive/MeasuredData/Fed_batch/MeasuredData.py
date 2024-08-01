# Python Libraries
import pandas as pd
import numpy as np

# My Libraries
from .pre_process import get_culture_volume, get_run_time
from CDPpy.helper import split_df, input_path, create_value_unit_df

from CDPpy.constants import ExpInfoNamespace as ExpInfoConst
from CDPpy.constants import ExpDataNamespace as ExpDataConst
from CDPpy.constants import CellNameSpace, OxygenNameSpace, ProductNameSpace


# Constants
TARGET_COLUMNS = [
    'Experiment Information', 
    'Experimental Data', 
    'Cell',
    'Oxygen',
    'Product',
    'Separate Feed Added', 
    'Concentration Before Feeding',
    'Concentration After Feeding', 
    'Feed Concentration', 
    'Calculated Cumulative Concentration',
    ]

class MeasuredData:
    '''
    Store measured data in a bioprocess experiment.
    '''
    def __init__(self, data):
        '''
        Parameters
        ----------
        data : pandas.DataFrame
            Experimental information.
        '''
        infokey = ExpInfoConst()
        datakey = ExpDataConst()
        cellkey = CellNameSpace()
        oxygenkey = OxygenNameSpace()
        productkey = ProductNameSpace()

        # Split DataFrame into sub DataFrames.
        df_list = split_df(data, TARGET_COLUMNS)
        exp_info_df = df_list[0]
        exp_data_df = df_list[1]
        cell_df = df_list[2]
        oxygen_df = df_list[3]
        product_df = df_list[4]
        separate_feed_df = df_list[5]
        conc_before_feed_df = df_list[6]
        conc_after_feed_df = df_list[7]
        feed_conc_df = df_list[8]
        measured_cumulative_conc_df = df_list[9]

        # Rename columns
        exp_info_df = exp_info_df.rename(columns={'Experiment ID': 'ID'})

        # Calculate run time (day) and (hour)
        run_time_day = exp_data_df[datakey.RUN_TIME_D]
        run_time_hour = exp_data_df[datakey.RUN_TIME_H]
        if run_time_day.isnull().all() and run_time_hour.isnull().all():
            date = exp_data_df[datakey.DATE]
            time = exp_data_df[datakey.TIME]
            run_time = get_run_time(date_df=date, time_df=time)
            run_time_day = run_time[datakey.RUN_TIME_D]
            run_time_hour = run_time[datakey.RUN_TIME_H]

        # Calculate culture volume
        volume_afeter_sampling = exp_data_df[datakey.VOLUME_AFTER_SAMPLE]
        if volume_afeter_sampling.isna().all():
            sample_size = exp_data_df[datakey.SAMPLE].size
            initial_volume = exp_info_df[infokey.INITIAL_VOLUME].iat[0]
            culture_volume = get_culture_volume(sample_size=sample_size, 
                                                initial_volume=initial_volume,
                                                sample_volume_df=exp_data_df[datakey.SAMPLE_VOLUME],
                                                separate_feed_df=separate_feed_df,
                                                base_added_df=exp_data_df[datakey.BASE_ADDED],
                                                feed_media_added_df=exp_data_df[datakey.FEED_MEDIA])
            volume_before_sampling = culture_volume[datakey.VOLUME_BEFORE_SAMPLE]
            volume_afeter_sampling = culture_volume[datakey.VOLUME_AFTER_SAMPLE]
        else:
            volume_before_sampling = pd.DataFrame(data=volume_afeter_sampling, columns=[datakey.VOLUME_BEFORE_SAMPLE])

        # Create pre-process data
        pre_date = pd.concat([exp_info_df, 
                              run_time_day, 
                              run_time_hour, 
                              volume_before_sampling, 
                              volume_afeter_sampling], axis=1)

        # Store data
        self._pre_data = pre_date

        self._exp_info_df = exp_info_df
        self._exp_data_df = exp_data_df

        self._separate_feed_df = separate_feed_df
        self._conc_before_feed_df = conc_before_feed_df
        self._conc_after_feed_df = conc_after_feed_df
        self._feed_conc_df = feed_conc_df
        self._measured_cumulative_conc_df = measured_cumulative_conc_df

        # Store Class constnat
        cell_line = exp_info_df[infokey.CELL_LINE].iat[0]
        self._cell_line = cell_line if type(cell_line)==str else str(cell_line)
        id = exp_info_df[infokey.RUN_ID].iat[0]
        self._exp_id = id if type(id)==str else str(id)
        self._initial_vlume = exp_info_df[infokey.INITIAL_VOLUME].iat[0]
        self._sample_size = exp_data_df[datakey.SAMPLE].size
        self._run_time_day = create_value_unit_df(run_time_day)
        self._run_time_hour = create_value_unit_df(run_time_hour)
        self._volume_before_sampling = create_value_unit_df(volume_before_sampling)
        self._volume_after_sampling = create_value_unit_df(volume_afeter_sampling)
        self._feed_media_added = create_value_unit_df(exp_data_df[datakey.FEED_MEDIA])
        self._viable_cell_conc = create_value_unit_df(cell_df[cellkey.VIABLE])
        self._dead_cell_conc = create_value_unit_df(cell_df[cellkey.DEAD])
        self._total_cell_conc = create_value_unit_df(cell_df[cellkey.TOTAL])
        self._production = create_value_unit_df(product_df[productkey.PRODUCTION])

    @property
    def cell_line(self):
        return self._cell_line
    @property
    def id(self):
        return self._exp_id
    @property
    def initial_volume(self):
        return self._initial_vlume
    @property
    def sample_size(self):
        return self._sample_size
    @property
    def run_time_day(self):
        return self._run_time_day
    @property
    def run_time_hour(self):
        return self._run_time_hour
    @property
    def volume_before_sampling(self):
        return self._volume_before_sampling
    @property
    def volume_after_sampling(self):
        return self._volume_after_sampling
    @property
    def feed_media_added(self):
        return self._feed_media_added
    @property
    def viable_cell_conc(self):
        return self._viable_cell_conc
    @property
    def dead_cell_conc(self):
        return self._dead_cell_conc
    @property
    def total_cell_conc(self):
        return self._total_cell_conc
    @property
    def production(self):
        return self._production
    
    @property
    def exp_info_df(self):
        return self._exp_info_df
    @property
    def pre_data(self):
        return self._pre_data
    
    @property
    def exp_data_df(self):
        return self._exp_data_df
        
    @property
    def separate_feed_df(self):
        return self._separate_feed_df
    @property
    def conc_before_feed_df(self):
        return self._conc_before_feed_df
    @property
    def conc_after_feed_df(self):
        return self._conc_after_feed_df
    @property
    def feed_conc_df(self):
        return self._feed_conc_df
    @property
    def measured_cumulative_conc_df(self):
        return self._measured_cumulative_conc_df
    