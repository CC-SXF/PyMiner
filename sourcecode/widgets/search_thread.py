# -*- coding: utf-8 -*-
"""
Created on Mon May 17 19:11:35 2021

@author: CC-SXF
"""


from PyQt5 import QtCore

from pathwaysearch.pathwaysearch import PathwaySearch


class SearchThread(QtCore.QThread):
    """
        Pathway search thread.
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    finished_signal = QtCore.pyqtSignal(list)
    prompts_signal = QtCore.pyqtSignal(str)

    def __init__(self, input_value = dict()):
        """ """
        super().__init__()
        self.input_value = input_value
        #
        self.ps_instance = PathwaySearch()
        self.ps_instance.ps_prompts_signal.connect(self.throw_prompts)

    def throw_prompts(self, prompts):
        """ """
        self.prompts_signal.emit(prompts)

    # Overriding
    def run(self):
        """
        """
        #
        source = self.input_value["Source"]
        target = self.input_value["Target"]
        avoid_compounds = self.input_value["Avoid Compounds"]
        avoid_reactions = self.input_value["Avoid Reactions"]
        #
        host_organism = self.input_value["Host Organism"]
        database = self.input_value["Database"]
        maximum_length = self.input_value["Maximum Length"]
        maximum_time = self.input_value[ "Maximum Time(s)"]
        #
        diffThreshold = self.input_value["DiffThreshold"]
        methods = self.input_value["Methods"]
        #
        isTotal = self.input_value["IsTotal"]
        isShortest = self.input_value["IsShortest"]
        #
        isRetro = self.input_value["IsRetro"]
        isSmart = self.input_value["IsSmart"]
        #
        isInfeasible = self.input_value["IsInfeasible"]
        isTransfer = self.input_value["IsTransfer"]
        #
        isAerobic = self.input_value["IsAerobic"]
        isFlux = self.input_value["IsFlux"]
        #
        isUpdate = self.input_value["IsUpdate"]
        #
        resource_file = self.input_value["Resource File"]
        general_cofactors_file = self.input_value["General Cofactors File"]
        metabolic_network_file = self.input_value["Metabolic Network File"]
        atom_mapping_file_location = self.input_value["Atom Mapping File Location"]
        #
        org_model_location = "../sourcedata/models/"
        #
        #
        pathwayinfo_list = self.ps_instance.searchPathway(source, target,  avoidCompSet = avoid_compounds, avoidRxnSet = avoid_reactions,
                                                          organism = host_organism, dataBase = database, maxlength = maximum_length, maxtime = maximum_time,
                                                          diffThreshold = diffThreshold, methods = methods,
                                                          isTotal = isTotal, isShortest = isShortest,
                                                          isRetro = isRetro, isSmart = isSmart,
                                                          isInfeasible = isInfeasible, isTransfer = isTransfer,
                                                          isAerobic = isAerobic, isFlux = isFlux,
                                                          isUpdate = isUpdate,
                                                          resource_file = resource_file,
                                                          general_cofactors_file = general_cofactors_file,
                                                          metabolic_network_file = metabolic_network_file,
                                                          atom_mapping_file_location = atom_mapping_file_location,
                                                          org_model_location = org_model_location,
                                                          )
        #
        self.finished_signal.emit(pathwayinfo_list)
        self.throw_prompts("")
        self.throw_prompts("<font color='blue', font-weight='bold', size=6>End.</font>")
        self.throw_prompts("")

    @classmethod
    def demoFunc(cls,):
         """ """
         pass


if __name__ == "__main__":
    """ """
    pass





