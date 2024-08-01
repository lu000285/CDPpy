import pandas as pd

from ...helper_func.helper_func import input_path

###########################################################################
def polyreg_calc(bio_process, polyorder_file):
    '''
    Calculate SP. rate for species using polynomial regression.

    Parameters
    ----------
        bio_process : BioProcess object
        polyorder_file : str
            name of a Excel file for polynomial regression order.
    '''
    method = 'polyreg'

    # Check polynomial file
    if (polyorder_file):
        path = input_path(file_name=polyorder_file)
        # polyorder df
        polyorder = pd.read_excel(io=path, index_col=0)
        polyorder.index = [name.upper() for name in polyorder.index]

    # Cell
    try:
        order = polyorder.loc['CELL'].iat[0]
    except:
        order = 3
    cell = bio_process.get_cell()
    cell.polyreg(polyorder=order)
    cell.set_method_flag(method=method, flag=True)

    # Oxygen
    try: 
        order = polyorder.loc['OXYGEN'].iat[0]
    except:
        order = 3
    oxygen = bio_process.get_oxygen()
    oxygen.polyreg(polyorder=order)
    oxygen.set_method_flag(method=method, flag=True)

    # IgG
    try:
        order = polyorder.loc['IGG'].iat[0]
    except:
        order = 3
    igg = bio_process.get_product()
    igg.polyreg(polyorder=order)
    igg.set_method_flag(method=method, flag=True)
    
    # Metabolites
    spc_dict = bio_process.get_spc_dict()
    spc_list = bio_process.get_spc_list()
    polyreg_df = pd.DataFrame()     # Initialize

    for spc_name in spc_list:
        spc = spc_dict[spc_name.upper()]    # species object
        try:
            order = polyorder.loc[spc_name].iat[0]
        except:
            order = 3
        
        spc.polyreg(polyorder=order)
        spc.set_method_flag(method=method, flag=True)

        title = f'Poly. Reg. Order: {order} q{spc_name.capitalize()} (mmol/109 cell/hr)'
        polyreg_df[title] = spc.get_sp_rate(method=method)

    # Add
    bio_process.set_process_data(process=method, data=polyreg_df)
    # SP. rate for Nitrogen and AA Carbon
    bio_process.sp_rate_others(method=method)

    # Ratio Calc
    bio_process.ratio_calc(method=method)

    # Set polyreg flag True
    bio_process.set_process_flag(process=method, flag=True)

# End polyreg_calc