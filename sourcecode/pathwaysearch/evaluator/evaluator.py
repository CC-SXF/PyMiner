# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 17:25:15 2021

@author: CC-SXF
"""

from PyQt5 import QtCore

from pathwaysearch.evaluator.endo_heter_inf import EvalEHI
from pathwaysearch.evaluator.atommap import AtomMap
from pathwaysearch.evaluator.flux import Flux


class Evaluator(QtCore.QObject):
    """
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    eva_prompts_signal = QtCore.pyqtSignal(str)

    def __init__(self, ):
        """ """
        super().__init__()


    def eva_throw_prompts(self, prompts):
        """ """
        self.eva_prompts_signal.emit(prompts)
        #


    def __sort(self, unifyPathwayList, ):
        """
            Sort the pathways by their implausible steps, exogenous steps, endogenous steps and string expression.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        #
        # '0-Feasibility';
        # '1-TotalLength'; '2-EndoLength'; '3-HeterLength'; '4-InfLength';
        # '5-AtomUtilization'; '6-AtomConservation'; '7-Flux';
        # '8-S-Details';
        # '9-I-Details'; '10-M-Details';
        #
        # '11-Route'                                   <----
        #
        #
        pathway_dict = dict()
        for pathway_info in unifyPathwayList:
            endolength = pathway_info[2]
            heterlength = pathway_info[3]
            inflength = pathway_info[4]
            pathway = pathway_info[8]
            #
            key = ""
            key += str(inflength).rjust(4, '0')
            key += str(heterlength).rjust(4, '0')
            key += str(endolength).rjust(4, '0')
            key += pathway
            #
            pathway_dict[key] = pathway_info
        #
        # Sorting...
        pathway_list = list()
        for key in sorted(pathway_dict.keys()):
            pathway_info = pathway_dict[key]
            pathway_list.append(pathway_info)
        #
        return pathway_list


    def _sort(self, unifyPathwayList, num_large = 1000000):
        """
            Sort the pathways by their implausible steps, exogenous steps, endogenous steps
        atom utilization, atom conservation and string expression.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        #
        # '0-Feasibility';
        # '1-TotalLength'; '2-EndoLength'; '3-HeterLength'; '4-InfLength';
        # '5-AtomUtilization'; '6-AtomConservation'; '7-Flux';
        # '8-S-Details';
        # '9-I-Details'; '10-M-Details';
        #
        # '11-Route'                                   <----
        #
        #
        pathway_dict = dict()
        for pathway_info in unifyPathwayList:
            endolength = pathway_info[2]
            heterlength = pathway_info[3]
            inflength = pathway_info[4]
            atom_utilization = pathway_info[5]
            atom_conservation = pathway_info[6]
            metabolic_flux = pathway_info[7]
            pathway = pathway_info[8]
            #
            if atom_utilization == None:
                atom_utilization = num_large
            else:
                atom_utilization = int(num_large - (100 * atom_utilization) )
            if atom_conservation == None:
                atom_conservation = num_large
            else:
                atom_conservation = int(num_large - (100 * atom_conservation) )
            if metabolic_flux == None:
                metabolic_flux = num_large
            else:
                # NaN
                try:
                    metabolic_flux = int(num_large - (100 * metabolic_flux) )
                except ValueError:
                    metabolic_flux = num_large
            #
            key = ""
            key += str(inflength).rjust(4, '0')
            key += ('_' + str(heterlength).rjust(4, '0') )
            key += ('_' + str(endolength).rjust(4, '0') )
            key += ('_' + str(atom_utilization).rjust(8, '0') )
            key += ('_' + str(atom_conservation).rjust(8, '0') )
            key += ('_' + str(metabolic_flux).rjust(8, '0') )
            key += ('_' + pathway )
            #
            pathway_dict[key] = pathway_info
        #
        # Sorting...
        pathway_list = list()
        for key in sorted(pathway_dict.keys()):
            pathway_info = pathway_dict[key]
            pathway_list.append(pathway_info)
        #
        return pathway_list


    def evaluator(self, pathwayEnumList,
                        rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                        organism = "eco", org_model = None, orgCompRxnId_dict = dict(),
                        org_compId_set = set(), org_rxnId_set = set(),
                        dataBase = "KEGG",
                        atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/",
                        isInfeasible = True, isTransfer = True, isFlux = True,
                        ):
        """
            Pathway evaluating and sorting.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        # '0-Feasibility';
        # '1-TotalLength'; '2-EndoLength'; '3-HeterLength'; '4-InfLength';
        # '5-AtomUtilization'; '6-AtomConservation'; '7-Flux';
        # '8-S-Details';
        # '9-I-Details'; '10-M-Details';
        #
        # '11-Route'                                   <----
        #
        #
        self.eva_throw_prompts("Endogenous/Heterogenous/Infeasible Length . . . ")
        pathway_list =  EvalEHI._evalEndoHeterInf(pathwayEnumList,
                                                  rxnDistInfoDict = rxnDistInfoDict,
                                                  org_compId_set = org_compId_set, org_rxnId_set = org_rxnId_set,
                                                  isInfeasible = isInfeasible,
                                                  )
        #
        if isTransfer:
            self.eva_throw_prompts("Atom Transfer . . . ")
            pathway_list = AtomMap._evalAtomTransfer(pathway_list,
                                                     dataBase = dataBase,
                                                     atom_mapping_file_location = atom_mapping_file_location,
                                                     )
        #
        if isFlux and (org_model != None):
            self.eva_throw_prompts("Metabolic Flux . . . ")
            pathway_list = Flux._evalFlux(pathway_list,
                                          rxnDistInfoDict = rxnDistInfoDict, compDistInfoDict = compDistInfoDict,
                                          org_model = org_model, organism = organism, orgCompRxnId_dict = orgCompRxnId_dict,
                                          )
        #
        # sorting...
        pathway_list = self._sort(pathway_list)
        #
        return pathway_list


    @classmethod
    def demoFunc(self, ):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





