import pandas as pd

#############################################################################################
def rolling_regression(bio_process, order=3, windows=6):
    '''
    Calculate SP. rate for species using rolling polynomial regression.

    Parameters
    ----------
        bio_process : BioProcess object
        order : int, defalut=3, optional
            polynomial order for rolling polynomial regression.
        windows : int, default=6, optional
            data point size used for rolling polynomial regression.
    '''
    '''
    '''