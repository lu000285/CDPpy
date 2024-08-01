import pandas as pd

class GetterMixin:
    '''Getter Mixin Class for In-Process
    '''
    @property
    def cumulative_conc(self):
        return self._cumulative_conc
    @property
    def sp_rate(self):
        return self._sp_rate
    @property
    def conc(self):
        return self._conc

    def get_in_process(self):
        """
        Get In-Process DataFrame.
        """
        if self._in_process_flag:
            t = self._run_time_hour
            cc = self._cumulative
            spr = self._sp_rate
        return pd.concat([t, cc, spr], axis=1)
