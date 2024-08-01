import pandas as pd


###########################################################################
def cumulative_calc(bio_process,
                    feed_concentration,
                    concentration_after_feed):
    '''
    Calculate cumulative consumptions/productions for species.

    Parameters
    ----------
        bio_process : BioProcess object
        feed_concentration : bool
        concentration_after_feed : bool
    '''

    # Cell
    '''cell = bio_process.get_cell()   # Get Cell object.
    cell.in_process()   # Calculate cumulative concentration
    cell.set_method_flag(method='cumulative', flag=True)    # Set cumulative flag true

    # Oxygen
    oxygen = bio_process.get_oxygen()   # Get Oxygen object.
    oxygen.in_process() # Calc. cumulative
    oxygen.set_method_flag(method='cumulative', flag=True)  # Set cumulative flag true

    # Product/IgG
    product = bio_process.get_product()
    product.in_process()    # Calc. cumulative
    product.set_method_flag(method='cumulative', flag=True) # Set cumulative flag true

    # Metabolite Cumulative
    spc_list = bio_process.get_spc_list()  # Get species list to analyze
    spc_dict = bio_process.get_spc_dict()
    spc_cumulative_df = pd.DataFrame()   # Initialize
    spc_conc_df = pd.DataFrame()
    spc_conc_after_feed_df = pd.DataFrame()

    for s in spc_list:
        s = s.upper()   # Name
        spc = spc_dict[s]   # Species object
        
        # Calculate cumulative consumption/production
        spc.in_process(feed_concentration=feed_concentration, 
                       concentration_after_feed=concentration_after_feed)
        spc.set_method_flag(method='cumulative', flag=True)
        unit = spc.get_cumulative_unit()    # Unit
        
        spc_cumulative_df['CUM '+s+' '+unit] = spc.get_cumulative()  # Add to DF
        spc_conc_df[s+' (mM) (before Feed)'] = spc.get_conc_before_feed()
        spc_conc_after_feed_df[s+' (mM) (After Feed)'] = spc.get_conc_after_feed()
        # spc_dict[s] = spc    # Add to Dictionary
    
    # Add to bp obj
    bio_process.set_spc_df(spc_cumulative_df)
    # bio_process.set_spc_dict(spc_dict=spc_dict)
    bio_process.set_spc_conc(spc_conc_df)
    bio_process.set_spc_conc_after_feed(conc_after_feed=spc_conc_after_feed_df)

    # Cumulative for Nitrogen, and AA Carbon
    bio_process.cumulative_others()

    # Set in process flag True
    bio_process.set_process_flag(process='inpro', flag=True)'''

# End cumulative_calc