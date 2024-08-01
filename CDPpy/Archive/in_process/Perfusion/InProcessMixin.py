import pandas as pd

class InProcessMixin:
    '''
    Mixin class for BioProcess class to do in-processing.
    Methods
    -------

    '''
    def in_process(self, a=0.25, c=3.0):
        '''
        Calculate cumulative consumptions/productions for species.
        Parameters
        ----------
            a : float, optional, default=0.25
                a recycling factor
            c : float, optional, default=3.0
                a concentration factor
        '''
        species = self._species

        # Get cell
        cell = species.pop('cell')
        cell.in_process(a, c)

        # Metabolite
        for s in species.values():
            s.in_process()

        # Update
        species.update({'cell': cell})
            
        # Set in process flag True
        # self.set_process_flag(process='inpro', flag=True)