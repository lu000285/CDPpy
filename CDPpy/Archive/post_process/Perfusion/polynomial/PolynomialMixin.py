import pandas as pd

from CDPpy.helper import input_path

class PolynomialMixin:
    '''
    '''
    def polynomial(self, polyorder_file, deg=3):
        '''
        '''
        # Get species dictionary
        species = self._species

        if polyorder_file:
            path = input_path(file_name=polyorder_file)
            # polyorder df
            polyorder = pd.read_excel(io=path, index_col=0)
            polyorder.index = [name.lower() for name in polyorder.index]

        # Get Cell
        cell = species.pop('cell')

        # Metabolite
        for name, s in species.items():
            try:
                degree = polyorder.loc[name].iat[0]
            except:
                degree = deg
            s.polynomial(deg=degree)

        # Update
        species.update({'cell': cell})
    