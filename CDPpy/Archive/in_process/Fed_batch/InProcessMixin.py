import pandas as pd

from ...Species.Fed_batch.Metabolite import Metabolite2

class InProcessMixin:
    '''
    Mixin class for BioProcess class to do in-processing.
    Methods
    -------

    '''
    def in_process(self):
        '''
        Calculate cumulative consumptions/productions for species.

        Parameters
        ----------
            feed_concentration : bool
                Set True if measured data has measurements of feed concentrations for species.
            concentration_after_feed : bool
                Set Treu if measured data has masurements of concentraions after feeding for species.
        '''
        species = self.get_species('all')
        species_list = list(species.keys())

        conc_df_list = []
        cumulative_conc_df_list = []
        sp_rate_df_list = []
        
        if 'cell' in species_list:
            species_list.remove('cell')
            cell = species['cell']
            cell.in_process()

        if 'product' in species_list:
            species_list.remove('product')
            prod = species['product']
            prod.in_process()

            conc_data = prod.conc
            conc_data['species'] = prod.name
            conc_df_list.append(conc_data)

            cumulative_data = prod.cumulative_conc
            cumulative_data['species'] = prod.name
            cumulative_conc_df_list.append(cumulative_data)

            sp_rate_data = prod.sp_rate
            sp_rate_data['species'] = prod.name
            sp_rate_df_list.append(sp_rate_data)

        if 'oxygen' in species_list:
            species_list.remove('oxygen')
            prod = species['oxygen']
            prod.in_process()

        for name in species_list:
            spc = species[name] # Species object
            spc.in_process(use_feed_conc=self.use_feed_conc, 
                           use_conc_after_feed=self.use_conc_after_feed)
            conc_data = spc.conc
            conc_data['species'] = name.capitalize()
            conc_df_list.append(conc_data)

            cumulative_data = spc.cumulative_conc
            cumulative_data['species'] = name.capitalize()
            cumulative_conc_df_list.append(cumulative_data)

            sp_rate_data = spc.sp_rate
            sp_rate_data['species'] = name.capitalize()
            sp_rate_df_list.append(sp_rate_data)
        
        self._conc_df = pd.concat(conc_df_list, axis=0).reset_index(drop=True)
        self._cumulative_conc_df = pd.concat(cumulative_conc_df_list, axis=0).reset_index(drop=True)
        self._sp_rate_df = pd.concat(sp_rate_df_list, axis=0).reset_index(drop=True)

        # Cumulative for Nitrogen, and AA Carbon
        # self.__cumulative_others()

        # Set in process flag True
        # self.set_process_flag(process='inpro', flag=True)

    #*** End inprocess ***#
    @property
    def conc(self):
        return self._conc_df
    @property
    def cumulative_conc(self):
        return self._cumulative_conc_df
    @property
    def sp_rate(self):
        return self._sp_rate_df


    #*** Provate Methods ***#
    def __cumulative_others(self):
            '''Calculate cumulative consumptions/productions for Nitrogen and AA carbon.
            '''
            # check default species list is a subset of species list of users
            if not (set(self._default_spc_list) <= set(self._spc_list)):
                return None

            df = self.get_process_data(method='inpro')
            
            ala = self._spc_dict['Alanine'.upper()].get_cumulative
            arg = self._spc_dict['Arginine'.upper()].get_cumulative
            asn = self._spc_dict['Asparagine'.upper()].get_cumulative
            asp = self._spc_dict['Aspartate'.upper()].get_cumulative
            cyt = self._spc_dict['Cystine'.upper()].get_cumulative
            gln = self._spc_dict['Glutamine'.upper()].get_cumulative
            glu = self._spc_dict['Glutamate'.upper()].get_cumulative
            gly = self._spc_dict['Glycine'.upper()].get_cumulative
            his = self._spc_dict['Histidine'.upper()].get_cumulative
            iso = self._spc_dict['Isoleucine'.upper()].get_cumulative
            leu = self._spc_dict['Leucine'.upper()].get_cumulative
            lys = self._spc_dict['Lysine'.upper()].get_cumulative
            met = self._spc_dict['Methionine'.upper()].get_cumulative
            nh3 = self._spc_dict['NH3'.upper()].get_cumulative
            phe = self._spc_dict['Phenylalanine'.upper()].get_cumulative
            pro = self._spc_dict['Proline'.upper()].get_cumulative
            ser = self._spc_dict['Serine'.upper()].get_cumulative
            thr = self._spc_dict['Threonine'.upper()].get_cumulative
            tryp = self._spc_dict['Tryptophan'.upper()].get_cumulative
            tyr = self._spc_dict['Tyrosine'.upper()].get_cumulative
            val = self._spc_dict['Valine'.upper()].get_cumulative

            # Check the length for data
            x = len(ala) + len(arg) + len(asn) + len(asp) + len(cyt)
            y = len(gln) + len(glu) + len(gly) + len(his) + len(iso)
            z = len(leu) + len(lys) + len(met) + len(nh3) + len(phe)
            w = len(pro) + len(ser) + len(thr) + len(tryp) + len(tyr) + len(val)
            total_len = x + y + z + w
            if len(ala) * 21 == total_len:
                # Nitrogen cumulative consumption/production
                x = ala*1 + arg*4 + asn*2 + asp*1  + cyt*2
                y = gln*2 + glu*1 + gly*1 + his*3  + iso*1
                z = leu*1 + lys*2 + met*1 - nh3*1  + phe*1
                w = pro*1 + ser*1 + thr*1 + tryp*2 + tyr*1 + val*1
                nitrogen_cum = (x + y + z + w).rename('CUM. Nitrogen (mmol)')

                # Nitrogen (w/o NH3, Ala) cumulative consumption/production
                nitrogen_w_o_NH3_Ala_cum = (nitrogen_cum + ala + nh3).rename('CUM. Nitrogen (w/o NH3, Ala) (mmol)')
                
                # Combine Nitrogen and Nitrogen (w/o NH3, Ala)
                #nitrogen_cum_df = pd.concat([nitrogen_cum, nitrogen_w_o_NH3_Ala_cum], axis=1)
                
                # Nitrogen Metabolite2 obj
                nitrogen = Metabolite2(name='Nitrogen',
                                       measured_data=self._md,                      
                                       cumulative=nitrogen_cum)
                nitrogen_wo = Metabolite2(name='Nitrogen (w/o NH3, Ala)',
                                          measured_data=self._md,                      
                                          cumulative=nitrogen_w_o_NH3_Ala_cum)
                # set flag true
                nitrogen.set_method_flag(method='cumulative', flag=True)
                nitrogen_wo.set_method_flag(method='cumulative', flag=True)
                
                # Add obj to spc dict and list
                self._spc_dict['NITROGEN'] = nitrogen
                self._spc_dict['Nitrogen (w/o NH3, Ala)'.upper()] = nitrogen_wo
                self._spc_list_2.append('Nitrogen'.upper())
                self._spc_list_2.append('Nitrogen (w/o NH3, Ala)'.upper())


                # Add cumulative consumption to DF
                df['CUM. Nitrogen (mmol)'] = nitrogen_cum
                df['CUM. Nitrogen (w/o NH3, Ala) (mmol)'] = nitrogen_w_o_NH3_Ala_cum


            # Check the length for data
            x = len(ala) + len(arg) + len(asn) + len(asp) + len(cyt)
            y = len(gln) + len(glu) + len(gly) + len(his) + len(iso)
            z = len(leu) + len(lys) + len(met) + len(phe) + len(pro)
            w = len(ser) + len(thr) + len(tryp) + len(tyr) + len(val)
            total_len = x + y + z + w
            if len(ala) * 20 == total_len:
                x = -ala*3 + arg*6 + asn*4   + asp*4 + cyt*6
                y =  gln*5 + glu*5 + gly*2   + his*6 + iso*6
                z =  leu*6 + lys*6 + met*5   + phe*9 + pro*5
                w =  ser*3 + thr*4 + tryp*11 + tyr*9 + val*5
                # AA carbon cumulative consumption/production
                aa_carbon_cum = (x+y+z+w).rename('CUM. AA Carbon (mmol)')

                # AA carbon Metabolite2 obj
                aa_carbon = Metabolite2(name='AA Carbon',
                                        measured_data=self._md,                      
                                        cumulative=aa_carbon_cum)
                aa_carbon.set_method_flag(method='cumulative', flag=True)

                # Add obj to spc dict and list
                self._spc_dict['AA CARBON'] = aa_carbon
                self._spc_list_2.append('aa carbon'.upper())

                # Add cumulative consumption to DF
                df['CUM. AA Carbon (mmol)'] = aa_carbon_cum

    #*** __cumulative_others ***#