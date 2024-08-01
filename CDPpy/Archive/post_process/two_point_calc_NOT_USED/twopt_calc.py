import pandas as pd

###########################################################################
def twopt_calc(bio_process):
    '''
    Calculate SP. rate for species using polynomial regression.

    Parameters
    ----------
        bio_process : BioProcess object
    '''
    method = 'twopt'

    # Cell
    cell = bio_process.get_cell()   # Get Cell object.
    cell.post_process_twopt()       # Calculate sp. rate with two-point calc.
    cell.set_method_flag(method=method, flag=True)  

    # Oxygen
    oxygen = bio_process.get_oxygen()   # Get Oxygen object.
    oxygen.post_process_twopt() # Calculate sp. rate with two-point calc.
    oxygen.set_method_flag(method=method, flag=True)    

    # IgG
    igg = bio_process.get_product() # Get Product object.
    igg.post_process_twopt()    # Calculate sp. rate with two-point calc.
    igg.set_method_flag(method=method, flag=True)   # Set flag.
    
    # Metabolites
    spc_dict = bio_process.get_spc_dict()
    spc_list = bio_process.get_spc_list()
    twopt_df = pd.DataFrame()   # Initialize

    for spc_name in spc_list:
        spc = spc_dict[spc_name.upper()]    # species object
        spc.sp_rate_twopt()                 # calculate SP. rate
        # set flag true
        spc.set_method_flag(method=method, flag=True)
        
        title = f'Two-Pt. Calc. q{spc_name.capitalize()} (mmol/109 cell/hr)'
        twopt_df[title] = spc_dict[spc_name.upper()].get_sp_rate(method='twopt')

    # Add to bp
    bio_process.set_process_data(process=method, data=twopt_df)

    # SP. rate for Nitrogen and AA Carbon
    bio_process.sp_rate_others(method=method)

    # Ratio Calc.
    bio_process.ratio_calc(method=method)

    # Set twopt flag true
    bio_process.set_process_flag(process=method, flag=True)

# End twopt_calc

