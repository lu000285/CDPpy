from CDPpy.MeasuredData.MeasuredData import MeasuredData

def pre_process(bio_process, exp_data, measured_data, feed_data):
    '''
    Execute pre process.

    Parameters
    ----------
        bio_process : python object
            BioProcess object
        exp_data : pands.DataFrame
            experiment information.
        measured_data : pands.DataFrame
            measured data.
        feed_data : pands.DataFrame
            separeata feed infomation.
    '''
    # MeausredData class
    md = MeasuredData(experiment_info=exp_data,
                      measured_data=measured_data,
                      feed_data=feed_data)

    # Calculate run time
    md.run_time()

    # Calculate culture volume before/after sampling and after feeding
    md.culture_volume()


#*** End pre_process ***#