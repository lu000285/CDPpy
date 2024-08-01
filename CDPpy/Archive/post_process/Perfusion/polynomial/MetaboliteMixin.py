import numpy as np
import pandas as pd

from CDPpy.Species.Perfusion.Variables import molar_mass_dict

class MetaboliteMixin:
    '''Metabolite Mixin Class for polynomial regression.
    Methods
    -------
        polynomial
    '''
    def polynomial(self, deg):
        '''Calculate the cumulative concentration and specific rate of a metabolite using polynomial regression.
        Parameters
        ----------
            deg : int
                a polynomial degree for polynomial regression.
        '''
        global molar_mass_dict

        self._polyorder = deg

        s = self._cumulative['value'].values
        t = self._run_time['value'].values
        xv = self._viable_cell_conc['value'].values
        unit = self._cumulative['unit'].iat[0]
        state = self._cumulative['state'].iat[0]

        # Fitting a polynomial
        poly_func = polyfit(t, s, deg)

        # Calculate cumulative concentration from the polynomial function
        s_poly = poly_func(t)
        s_poly = pd.DataFrame(data=s_poly, columns=['value'])
        s_poly['unit'] = unit
        s_poly['state'] = state
        s_poly['degree'] = deg
        s_poly.index.name = 'Cumulative concentration'

        # Get the derivetive of the polynomial, and evaluate the derivetive at the run time.
        poly_deriv = poly_func.deriv()
        y = poly_deriv(t)

        # Calculate the specific rate from the derivetive of the polynomial function
        r_poly = y / xv
        if unit==('(mg/ml)'): # Match the unit
            molar_mass = molar_mass_dict[self._name.capitalize()]
            r_poly = r_poly / molar_mass * 1000
        r_poly = pd.DataFrame(data=r_poly, columns=['value'])
        r_poly['unit'] = '(mmol/10^9_cells/hr)'
        r_poly['degree'] = deg
        r_poly.index.name = 'Specific rate'

        # Store the variables
        self._poly_func = poly_func
        self._cumulative_poly = s_poly
        self._sp_rate_poly = r_poly

    @property
    def get_cumulative_poly(self):
        ''''''
        return self._cumulative_poly
    
    @property
    def get_sp_rate_poly(self):
        ''''''
        return self._sp_rate_poly
    

def polyfit(t, s, deg):
    '''Fitting a polynomial of the degree to the cumulative concentration data.
    Parameters
    ----------
        t : numpy array
            Run time (hr).
        s : numpy array
            Cumulative concentration (mg/mL) or (mmol/L).
        deg : int
            a polynomial degree for polynomial regression.
    Returns
    -------
        np.poly1d()
            a polynomial function that fit to the cumulative concentration data.
    '''
    return np.poly1d(np.polyfit(x=t, y=s, deg=deg))


        