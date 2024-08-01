import pandas as pd
import numpy as np

from .GetterMixin import GetterMixin

from .equation import integral_viable_cell, cell_production, growth_rate

from CDPpy.Constants import CellNameSpace

CONSTANTS = CellNameSpace()

class CellMixin(GetterMixin):
    '''
    '''
    # Call methods
    def in_process(self):
        '''Caluculate in-process data for cell.
        '''
        # Get run time dataframe
        run_time = self.run_time

        # IVCC
        ivcc = self.integral_viable_cell_calc()
        self._ivcc = pd.concat([run_time, ivcc], axis=1)

        # Cumulative Consumption/Production
        s = self.cumulative_calc()
        self._cumulative_conc = pd.concat([run_time, s], axis=1)

        # Specific rate
        r = self.sp_rate_calc()
        self._sp_rate = pd.concat([run_time, r], axis=1)

        # Viable, Total concentrations and Viability
        viable = pd.concat([run_time, self._viable_cell_conc], axis=1)
        viable['state'] = 'Viable'
        total = pd.concat([run_time, self._total_cell_conc], axis=1)
        total['state'] = 'Total'
        viab = pd.concat([run_time, self._viability], axis=1)
        viab['state'] = 'Viability'
        self._conc = pd.concat([viable, total, viab], axis=0).reset_index(drop=True)

        '''
        kd = self.kd
        self._kd = pd.Series(data=kd, name='kd')
        mv = self.mv
        self._mv = pd.Series(data=mv, name='mv')
        '''

    def integral_viable_cell_calc(self) -> pd.DataFrame:
        '''Calculate Integral of Viable Cell
        Parameters
        ---------
            idx : numpy array
                Indices of measurements.
            c : numpy array
                Viable cell concentration (10^6 cells/mL).
            t : numpy array
                Time (hr)
        Returns
            s : pandas DataFrame
                Integral of viable cell concentration (10^6 cells hr/mL).
        '''
        idx = self.measurement_index
        c = self.viable_cell_conc['value'].values[idx]
        t = self.run_time_hour['value'].values[idx]

        # Initialize
        s = np.zeros(t.size)
        for i in range(1, len(t)):
            s_i = integral_viable_cell(c[i-1], c[i], t[i-1], t[i])
            s[idx[i]] = s[idx[i-1]] + s_i

        s = pd.DataFrame(data=s, columns=['value'])
        s['unit'] = CONSTANTS.IVCC_UNIT
        s.index.name = CONSTANTS.IVCC
        return s
        
    # Calculate Cumulative Cell Produced
    def cumulative_calc(self):
        '''Calculate cumulative production of the cell.
        Parameters
        ---------
            idx : numpy array
                Indices of measurements.
            c : numpy array
                Viable cell concentration (10^6 cells/mL).
            v1 : numpy array
                Culture volume before sampling (mL).
            v2 : numpy array
                Culture volume after sampling (mL).
        Returns
            s : pandas DataFrame
                Cumulative production (10^6 cells).
        '''
        idx = self.measurement_index
        c = self.viable_cell_conc['value'].values[idx]
        v1 = self.volume_before_sampling['value'].values[idx]
        v2 = self.volume_after_sampling['value'].values[idx]

        # Initialize
        s = np.zeros(c.size)
        s.fill(np.nan)
        s[0] = 0.0

        for i in range(1, len(idx)):
            # c[i] = c[i-1] + xv[i] * v1[i] - xv[i-1] * v2[i-1]
            s_i = cell_production(c[i-1], c[i], v1[i], v2[i-1])
            s[i] = s_i + s[i-1]
        
        s = pd.DataFrame(data=s, columns=['value'])
        s['unit'] = CONSTANTS.CUMULATIVE_UNIT
        s.index.name = CONSTANTS.CUMULATIVE
        return s

    def sp_rate_calc(self):
        '''Calculates Specific growth rate.
        Parameters
        ---------
            idx : numpy array
                Indices of measurements.
            s : numpy array
                Cumulative cell concentration (10^6 cells).
            t : numpy array
                Time (hr).
            v1 : numpy array
                Culture volume before sampling (mL).
            v2 : numpy array
                Culture volume after sampling (mL).
        Returns
            r : pandas DataFrame
                Specific growth rate (hr^1).
        '''
        idx = self.measurement_index
        t = self.run_time_hour['value'].values[idx]
        xv = self.viable_cell_conc['value'].values[idx]                   # vialbe cell concentration (10e6 cells/ml)
        s = self.cumulative_conc['value'].values[idx]            # Cumulative Cell Concentraion (10e6 cells/mL)
        v1 = self.volume_before_sampling['value'].values[idx]
        v2 = self.volume_after_sampling['value'].values[idx]
        sample_size = self.samples

        # Initialize
        r = np.zeros(sample_size)
        r.fill(np.nan)
        for i in range(1, len(t)):
            r[i] = growth_rate(s[i-1], s[i], xv[i-1], xv[i], v1[i], v2[i-1], t[i-1], t[i])
        r = pd.DataFrame(data=r, columns=['value'])
        r['unit'] = CONSTANTS.SP_RATE_UNIT
        r['method'] = 'twoPoint'
        r.index.name = CONSTANTS.SP_RATE
        return r

    # Calculates kd value
    @property
    def kd(self):
        idx = self._idx
        xv = self._xv.values[idx]                   # Vialbe cell concentration (10e6 cells/ml)
        xd = self._xd.values[idx]                   # Dead Cell Concentration (10e6 cells/mL)
        t = self._run_time_hour.values[idx]         # Run Time (hrs)
        v1 = self._v_before_sampling.values[idx]    # Culture Volume Before Sampling (mL)
        v2 = self._v_after_sampling.values[idx]     # Culture Volume After Sampling (mL)

        # Initialize
        kd = np.zeros(self._sample_num)
        kd.fill(np.nan)
        for i in range(1, len(idx)):
            x = xd[i] * v1[i] - xd[i-1] * v2[i-1]
            y = xv[i] * v1[i] + xv[i-1] * v2[i-1]
            kd[i] = x / (y * 0.5 * (t[i] - t[i-1]))

        return kd

    # Calculate m
    @property
    def mv(self):
        r = self._sp_rate.values
        kd = self._kd.values
        return (r + kd)

    @property
    def integral_viable_cell_conc(self):
        return self._ivcc    
    @property
    def growth_rate(self):
        return self._sp_rate
        
    # Display
    @property
    def disp_in_process(self):
        if self._in_process_flag:
            data = self.get_in_process
            print('\n************ Cell In Process Data ************')
            print(data)