import numpy as np
import pandas as pd

from CDPpy.constants import ExpDataNamespace as Constants

def get_run_time(date_df, time_df) -> pd.DataFrame:
    '''Calculate run time (day, hour) and return time DataFrame.
    '''
    key = Constants()
    if pd.api.types.is_datetime64_any_dtype(date_df):
        date_df = date_df.apply(lambda x: x.date())
    if pd.api.types.is_datetime64_any_dtype(time_df):
        time_df = time_df.apply(lambda x: x.time())

    if (date_df==time_df).all():
        date_time = pd.to_datetime(date_df.astype(str))
    else:
        date_time = pd.to_datetime(date_df.astype(str) + ' ' + time_df.astype(str))
    
    time_diff = date_time - date_time.iat[0]
    run_time_hour = time_diff.dt.total_seconds() / 3600.0
    run_time_day = time_diff.dt.total_seconds() / 3600.0 / 24.0
    # run_time_day = time_diff.dt.days

    df = pd.DataFrame()
    #df[key.DATE] = date_df
    #df[key.TIME] = time_df
    df[key.RUN_TIME_D] = run_time_day#.reset_index(drop=True)
    df[key.RUN_TIME_H] = run_time_hour#.reset_index(drop=True)
    return df

def get_culture_volume(sample_size, initial_volume, 
                       sample_volume_df, separate_feed_df, 
                       base_added_df, feed_media_added_df) -> pd.DataFrame:
    '''Calculate culture volume before/after sampling and after feeding from Meaured Data.
    '''
    key = Constants()

    # DataFrame to arrays
    feed_data = separate_feed_df.fillna(0)
    base_added = base_added_df.fillna(0).values
    feed_media_added = feed_media_added_df.fillna(0).values
    sample_volume = sample_volume_df.fillna(0).values

    # Initialize
    v_before_sampling = np.zeros(sample_size)
    v_before_sampling[0] = initial_volume
    v_after_feeding = np.zeros(sample_size)
    v_after_sampling = np.zeros(sample_size)
    
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
    
    df = pd.DataFrame()
    df[key.VOLUME_BEFORE_SAMPLE] = v_before_sampling
    # df[key.VOLUME_AFTER_FEED] = v_after_feeding
    df[key.VOLUME_AFTER_SAMPLE] = v_after_sampling
    return df
  