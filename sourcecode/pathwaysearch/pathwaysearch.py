# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 16:23:43 2020

@author: CC-SXF
"""

import json
import shutil
from os import path, listdir, mkdir
from cobra.io.json import load_json_model
from cobra.util.solver import linear_reaction_coefficients
from PyQt5 import QtCore

from pathwaysearch.input.input import Input
from pathwaysearch.pwbfs_naive import PwNaBfs
from pathwaysearch.pwdfs_naive import PwNaDfs
from pathwaysearch.evaluator.evaluator import Evaluator

from pathwaysearch.save import Save


class PathwaySearch(QtCore.QObject):
    """
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    ps_prompts_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        """ """
        super().__init__()
        #
        # # self.preparation()
        #
        self.pwnabfs_instance = PwNaBfs()
        self.pwnabfs_instance.bfs_prompts_signal.connect(self.ps_throw_prompts)
        self.pwnadfs_instance = PwNaDfs()
        self.pwnadfs_instance.dfs_prompts_signal.connect(self.ps_throw_prompts)
        #
        self.evaluator_instance = Evaluator()
        self.evaluator_instance.eva_prompts_signal.connect(self.ps_throw_prompts)


    def preparation(self,):
        """ """
        #
        folder_list = listdir()
        if "temp" not in folder_list:
            mkdir("temp")
        if "results" not in folder_list:
            mkdir("results")
        #
        # clear....
        for sub_route in listdir("results"):
            sub_route = path.join("results", sub_route)
            shutil.rmtree(sub_route)
        mkdir('results/ps_results/')
        #


    def ps_throw_prompts(self, prompts):
        """ """
        self.ps_prompts_signal.emit(prompts)
        #


    def searchPathway(self, compStrSet, compStr, avoidCompSet = set(), avoidRxnSet = set(),
                            organism = "eco", dataBase = "KEGG", maxlength = 16, maxtime = 100, maxnumber = 1000000,
                            diffThreshold = 0.1, methods = "bfs_naive",
                            isTotal = True,  isShortest = False,
                            isRetro = True, isSmart = True,
                            isInfeasible = False, isTransfer = True,
                            isAerobic = True, isFlux = True,
                            isUpdate = False,
                            resource_file = "../sourcedata/KEGG_Ori_Cleaned_Half.xlsx",
                            general_cofactors_file = "datas/KEGG_General_Cofactors.xlsx",
                            metabolic_network_file = "../sourcedata/KEGG_Pathway_Search_Ori.xlsx",
                            atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/",
                            org_model_location = "../sourcedata/models/",
                            max_growth_percent = 0.8,
                            diffThresholdLevel = 10,
                            isDebug = False,
                            ):
        """
        """
        #
        self.preparation()
        #
        if (len(compStrSet) == 0) and (organism == ""):
            self.ps_throw_prompts("<font color='red'>Source and Host Organism can't be null simultaneously.</font>")
            return [list(), dict(), dict()]
        if  compStr == "":
            self.ps_throw_prompts("<font color='red'>Target can't be null.</font>")
            return [list(), dict(), dict()]
        #
        # models
        # "" : ""
        # bsu: Bacillus subtilis subsp. subtilis 168
        # eco: Escherichia coli K-12 MG1655
        # kpn: Klebsiella pneumoniae subsp. pneumoniae MGH 78578
        # ppu: Pseudomonas putida KT2440
        # sce: Saccharomyces cerevisiae S288c
        # syz: Synechocystis sp. PCC 6803(Cyanobacteria)
        # ...
        #
        with open(path.join("datas", "models.json"), "r") as file_obj:
            models = json.load(file_obj)
        json_model_file_dict = dict()
        for key, value in models.items():
            json_model_file_dict[key] = value[0]
        #
        if dataBase == "KEGG":
            #
            # Generating...
            # Generate Compound Pairs.
            Input.generateRxnCompPair(resource_file = resource_file,
                                      general_cofactors_file = general_cofactors_file,
                                      metabolic_network_file = metabolic_network_file,
                                      atom_mapping_file_location = atom_mapping_file_location,
                                      org_model_location = org_model_location,
                                      dataBase = dataBase,
                                      isUpdate = isUpdate,
                                      diffThresholdLevel = diffThresholdLevel,
                                      )
            #
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient List;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = Input.loadRxnCompPair(metabolic_network_file = metabolic_network_file,
                                                                        diffThreshold = diffThreshold,
                                                                        dataBase = dataBase,
                                                                        isUpdate = isUpdate,
                                                                        diffThresholdLevel = diffThresholdLevel,
                                                                        )
            #
            # Loading...
            orgCompRxnId_dict = Input.getOrgCompRxnId(org_model_location = org_model_location,
                                                      organism = organism,
                                                      json_model_file_dict = json_model_file_dict,
                                                      dataBase = dataBase,
                                                      )
            #
            # Loading...
            if organism == "":
                org_model = None
            else:
                json_model_file = json_model_file_dict[organism]
                json_model_file = path.join(org_model_location, json_model_file)
                org_model = load_json_model(json_model_file)
                #
                if organism == "syz":
                    # Mixotrophic Biomass Growth
                    org_model.objective = org_model.reactions.BIOMASS_Ec_SynMixo_1
                #
                if isAerobic:
                    org_model.reactions.EX_o2_e.bounds = (-1000, 1000)          # O2(aerobic)
                else:
                    org_model.reactions.EX_o2_e.bounds = (0, 1000)              # O2(anaerobic)
                org_model.reactions.EX_glc__D_e.bounds = (-10, 1000)            # D-Glucose
                #
                max_growth = org_model.slim_optimize()
                max_growth_threshold = max_growth_percent * max_growth
                max_growth_threshold = round(max_growth_threshold, 2)
                objective_reaction = list(linear_reaction_coefficients(org_model).keys())[0]
                objective_reaction.bounds = (max_growth_threshold, max_growth_threshold)
            #
            # Excluding...
            for compId in avoidCompSet:
                try:
                    compDistInfoDict[compId][2] = list()
                except KeyError:
                    pass
            for rxnId in avoidRxnSet:
                try:
                    # # rxnDistInfoDict[rxnId][3] = [{}, {}]
                    rxnDistInfoDict[rxnId][3] = {"forward":[dict(), dict()], 'backward':[dict(), dict()]}
                except KeyError:
                    pass
            #
            # Microbial Chassis
            # "bsu": "Bacillus subtilis subsp. subtilis 168"
            # # "cgl": "Corynebacterium glutamicum ATCC 13032 (Kyowa Hakko)"
            # "eco": "Escherichia coli K-12 MG1655"
            # "kpn": "Klebsiella pneumoniae subsp. pneumoniae MGH 78578"
            # # "ppa": "Komagataella phaffii GS115 (Pichia pastoris GS115)"
            # "ppu": "Pseudomonas putida KT2440"
            # "sce": "Saccharomyces cerevisiae S288c"
            # "syz": "Synechocystis sp. PCC 6803(Cyanobacteria)"
            #
            # ...
            org_compId_set = set()
            org_rxnId_set = set()
            try:
                comp_org = "".join(['Compound', '(', organism, ')'])
                org_compId_set = set(orgCompRxnId_dict[comp_org].keys())
                rxn_org = "".join(['Reaction', '(', organism, ')'])
                org_rxnId_set = set(orgCompRxnId_dict[rxn_org].keys())
            except KeyError:
                # organism = ""
                org_compId_set = set(compDistInfoDict.keys())
                org_rxnId_set = set(rxnDistInfoDict.keys())
            #
            # Exclude orphan metabolites in the microbial chassis we choose.
            # The reaction number of compound C00001, C00007 or C00080 is beyond the capacity of
            # one cell of one excel file.
            # And compound C00001(H20), C00007(O2) or C00080(H+) doesn't belong to orphan metabolites.
            nonOrphan_metId_set = {'C00001', 'C00007', 'C00080'}
            orphan_metId_set = set()
            for compId in (org_compId_set - nonOrphan_metId_set):
                try:
                    compRxnList = compDistInfoDict[compId][2]  # Reactions
                    compRxnSet = set(compRxnList)
                    if len(compRxnSet & org_rxnId_set) == 0:
                        orphan_metId_set.add(compId)
                    else:
                        pass
                except KeyError:
                    pass
            org_compId_set -= orphan_metId_set
            #
            # Pathway search based on the microbial chassis we choose.
            # Search from target to endogenous metabolites of the microbial chassis we choose.
            if len(compStrSet) == 0:
                compStrSet = org_compId_set
            #
        elif dataBase == "Rhea":
            pass
        elif dataBase == "MetaCyc":
            # Generating...
            # Generate Compound Pairs.
            Input.generateRxnCompPair(resource_file = resource_file,
                                      general_cofactors_file = general_cofactors_file,
                                      metabolic_network_file = metabolic_network_file,
                                      atom_mapping_file_location = atom_mapping_file_location,
                                      org_model_location = org_model_location,
                                      dataBase = dataBase,
                                      isUpdate = isUpdate,
                                      diffThresholdLevel = diffThresholdLevel,
                                      )
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient List;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = Input.loadRxnCompPair(metabolic_network_file = metabolic_network_file,
                                                                        diffThreshold = diffThreshold,
                                                                        dataBase = dataBase,
                                                                        isUpdate = isUpdate,
                                                                        diffThresholdLevel = diffThresholdLevel,
                                                                        )
            # Loading...
            orgCompRxnId_dict = Input.getOrgCompRxnId(org_model_location = org_model_location,
                                                      organism = organism,
                                                      json_model_file_dict = json_model_file_dict,
                                                      dataBase = dataBase,
                                                      )
            # Loading...
            if organism == "":
                org_model = None
            else:
                json_model_file = json_model_file_dict[organism]
                json_model_file = path.join(org_model_location, json_model_file)
                org_model = load_json_model(json_model_file)
                #
                if organism == "syz":
                    # Mixotrophic Biomass Growth
                    org_model.objective = org_model.reactions.BIOMASS_Ec_SynMixo_1
                #
                if isAerobic:
                    org_model.reactions.EX_o2_e.bounds = (-1000, 1000)          # O2(aerobic)
                else:
                    org_model.reactions.EX_o2_e.bounds = (0, 1000)              # O2(anaerobic)
                org_model.reactions.EX_glc__D_e.bounds = (-10, 1000)            # D-Glucose
                #
                max_growth = org_model.slim_optimize()
                max_growth_threshold = max_growth_percent * max_growth
                max_growth_threshold = round(max_growth_threshold, 2)
                objective_reaction = list(linear_reaction_coefficients(org_model).keys())[0]
                objective_reaction.bounds = (max_growth_threshold, max_growth_threshold)
            #
            # Excluding...
            for compId in avoidCompSet:
                try:
                    compDistInfoDict[compId][2] = list()
                except KeyError:
                    pass
            for rxnId in avoidRxnSet:
                try:
                    # # rxnDistInfoDict[rxnId][3] = [{}, {}]
                    rxnDistInfoDict[rxnId][3] = {"forward":[dict(), dict()], 'backward':[dict(), dict()]}
                except KeyError:
                    pass
            #
            # Microbial Chassis
            # "bsu": "Bacillus subtilis subsp. subtilis 168"
            # # "cgl": "Corynebacterium glutamicum ATCC 13032 (Kyowa Hakko)"
            # "eco": "Escherichia coli K-12 MG1655"
            # "kpn": "Klebsiella pneumoniae subsp. pneumoniae MGH 78578"
            # # "ppa": "Komagataella phaffii GS115 (Pichia pastoris GS115)"
            # "ppu": "Pseudomonas putida KT2440"
            # "sce": "Saccharomyces cerevisiae S288c"
            # "syz": "Synechocystis sp. PCC 6803(Cyanobacteria)"
            #
            # ...
            org_compId_set = set()
            org_rxnId_set = set()
            try:
                comp_org = "".join(['Compound', '(', organism, ')'])
                org_compId_set = set(orgCompRxnId_dict[comp_org].keys())
                rxn_org = "".join(['Reaction', '(', organism, ')'])
                org_rxnId_set = set(orgCompRxnId_dict[rxn_org].keys())
            except KeyError:
                # organism = ""
                org_compId_set = set(compDistInfoDict.keys())
                org_rxnId_set = set(rxnDistInfoDict.keys())
            #
            # Exclude orphan metabolites in the microbial chassis we choose.
            # The reaction number of compound PROTON, WATER, OXYGEN-MOLECULE, CO-A is beyond
            # the capacity of one cell of one excel file.
            # And compound PROTON(H+), WATER(H2O), OXYGEN-MOLECULE(O2), CO-A(CoA-SH) doesn't belong
            # to orphan metabolites.
            # PROTON(Met000001-m), WATER(Met000002-m), OXYGEN-MOLECULE(Met000003-m), CO-A(Met000007-m)
            nonOrphan_metId_set = {'Met000001-m', 'Met000002-m', 'Met000003-m', 'Met000007-m',}
            orphan_metId_set = set()
            for compId in (org_compId_set - nonOrphan_metId_set):
                try:
                    compRxnList = compDistInfoDict[compId][2]  # Reactions
                    compRxnSet = set(compRxnList)
                    if len(compRxnSet & org_rxnId_set) == 0:
                        orphan_metId_set.add(compId)
                    else:
                        pass
                except KeyError:
                    pass
            org_compId_set -= orphan_metId_set
            #
            # Pathway search based on the microbial chassis we choose.
            # Search from target to endogenous metabolites of the microbial chassis we choose.
            if len(compStrSet) == 0:
                compStrSet = org_compId_set
            #
        elif dataBase == "KndPad":
            # Generating...
            # Generate Compound Pairs.
            Input.generateRxnCompPair(resource_file = resource_file,
                                      general_cofactors_file = general_cofactors_file,
                                      metabolic_network_file = metabolic_network_file,
                                      atom_mapping_file_location = atom_mapping_file_location,
                                      org_model_location = org_model_location,
                                      dataBase = dataBase,
                                      isUpdate = isUpdate,
                                      diffThresholdLevel = diffThresholdLevel,
                                      )
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient List;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = Input.loadRxnCompPair(metabolic_network_file = metabolic_network_file,
                                                                        diffThreshold = diffThreshold,
                                                                        dataBase = dataBase,
                                                                        isUpdate = isUpdate,
                                                                        diffThresholdLevel = diffThresholdLevel,
                                                                        )
            # Loading...
            orgCompRxnId_dict = Input.getOrgCompRxnId(org_model_location = org_model_location,
                                                      organism = organism,
                                                      json_model_file_dict = json_model_file_dict,
                                                      dataBase = dataBase,
                                                      )
            # Loading...
            if organism == "":
                org_model = None
            else:
                json_model_file = json_model_file_dict[organism]
                json_model_file = path.join(org_model_location, json_model_file)
                org_model = load_json_model(json_model_file)
                #
                if organism == "syz":
                    # Mixotrophic Biomass Growth
                    org_model.objective = org_model.reactions.BIOMASS_Ec_SynMixo_1
                #
                if isAerobic:
                    org_model.reactions.EX_o2_e.bounds = (-1000, 1000)          # O2(aerobic)
                else:
                    org_model.reactions.EX_o2_e.bounds = (0, 1000)              # O2(anaerobic)
                org_model.reactions.EX_glc__D_e.bounds = (-10, 1000)            # D-Glucose
                #
                max_growth = org_model.slim_optimize()
                max_growth_threshold = max_growth_percent * max_growth
                max_growth_threshold = round(max_growth_threshold, 2)
                objective_reaction = list(linear_reaction_coefficients(org_model).keys())[0]
                objective_reaction.bounds = (max_growth_threshold, max_growth_threshold)
            #
            # Excluding...
            for compId in avoidCompSet:
                try:
                    compDistInfoDict[compId][2] = list()
                except KeyError:
                    pass
            for rxnId in avoidRxnSet:
                try:
                    # # rxnDistInfoDict[rxnId][3] = [{}, {}]
                    rxnDistInfoDict[rxnId][3] = {"forward":[dict(), dict()], 'backward':[dict(), dict()]}
                except KeyError:
                    pass
            #
            # Microbial Chassis
            # "bsu": "Bacillus subtilis subsp. subtilis 168"
            # # "cgl": "Corynebacterium glutamicum ATCC 13032 (Kyowa Hakko)"
            # "eco": "Escherichia coli K-12 MG1655"
            # "kpn": "Klebsiella pneumoniae subsp. pneumoniae MGH 78578"
            # # "ppa": "Komagataella phaffii GS115 (Pichia pastoris GS115)"
            # "ppu": "Pseudomonas putida KT2440"
            # "sce": "Saccharomyces cerevisiae S288c"
            # "syz": "Synechocystis sp. PCC 6803(Cyanobacteria)"
            #
            # ...
            org_compId_set = set()
            org_rxnId_set = set()
            try:
                comp_org = "".join(['Compound', '(', organism, ')'])
                org_compId_set = set(orgCompRxnId_dict[comp_org].keys())
                rxn_org = "".join(['Reaction', '(', organism, ')'])
                org_rxnId_set = set(orgCompRxnId_dict[rxn_org].keys())
            except KeyError:
                # organism = ""
                org_compId_set = set(compDistInfoDict.keys())
                org_rxnId_set = set(rxnDistInfoDict.keys())
            #
            # Exclude orphan metabolites in the microbial chassis we choose.
            # The reaction number of compound Met000001-m, Met000002-m, Met000003-m  or Met000007-m is beyond
            # the capacity of one cell of one excel file.
            # And compound Met000001-m(H+), Met000002-m(H2O), Met000003-m(O2), Met000007-m(CoA-SH) doesn't belong
            # to orphan metabolites.
            nonOrphan_metId_set = {'Met000001-m', 'Met000002-m', 'Met000003-m', 'Met000007-m',}
            orphan_metId_set = set()
            for compId in (org_compId_set - nonOrphan_metId_set):
                try:
                    compRxnList = compDistInfoDict[compId][2]  # Reactions
                    compRxnSet = set(compRxnList)
                    if len(compRxnSet & org_rxnId_set) == 0:
                        orphan_metId_set.add(compId)
                    else:
                        pass
                except KeyError:
                    pass
            org_compId_set -= orphan_metId_set
            #
            # Pathway search based on the microbial chassis we choose.
            # Search from target to endogenous metabolites of the microbial chassis we choose.
            if len(compStrSet) == 0:
                compStrSet = org_compId_set
            #
        else:
            pass
        #
        #
        #
        #
        self.ps_throw_prompts("")
        self.ps_throw_prompts("<font color='blue', font-weight='bold', size=6>Search:</font>")
        #
        # Retro-Search(In other words, search from target to mult-sources) will be used.
        if (len(compStrSet) >= 2):
            multSourceSet = compStrSet
            target = compStr
            # Shortest pathways search.
            if isShortest:
                if methods == "bfs_wl":
                    pass
                if methods == "bfs_naive":
                    pathwayEnumList = self.pwnabfs_instance.Target2MultSource(multSourceSet, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
            # Total pathways search.
            else:
                if methods == "bfs_naive":
                    pathwayEnumList = self.pwnabfs_instance.Target2MultSource(multSourceSet, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                if methods == "dfs_naive":
                    pathwayEnumList = self.pwnadfs_instance.Target2MultSource(multSourceSet, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
        elif isSmart:
            source = list(compStrSet)[0]
            target = compStr
            try:
                source_bf = compDistInfoDict[source][1][1]  # Branching Factor(out-degree);
            except KeyError:
                self.ps_throw_prompts(f"<font color='red'>The source({source}) does not exist in the Generalized Metabolic Space.</font>")
                self.ps_throw_prompts("<font color='red'>Please have a check on the source.</font>")
                return [list(), dict(), dict()]
            try:
                target_bf = compDistInfoDict[target][1][0]  # Branching Factor(in-degree);
            except KeyError:
                self.ps_throw_prompts(f"<font color='red'>The target({target}) does not exist in the Generalized Metabolic Space.</font>")
                self.ps_throw_prompts("<font color='red'>Please have a check on the the target.</font>")
                return [list(), dict(), dict()]
            # Shortest pathways search.
            if isShortest:
                if methods == "bfs_wl":
                    pass
                if methods == "bfs_naive":
                    if source_bf < target_bf:
                        pathwayEnumList = self.pwnabfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnabfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  isShortest,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
            # Total pathways search.
            else:
                if methods == "bfs_naive":
                    if source_bf < target_bf:
                        pathwayEnumList = self.pwnabfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnabfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  isShortest,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
                if methods == "dfs_naive":
                    if source_bf < target_bf:
                        pathwayEnumList = self.pwnadfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnadfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
        else:
            source = list(compStrSet)[0]
            target = compStr
            # Shortest pathways search.
            if isShortest:
                if methods == "bfs_wl":
                    pass
                if methods == "bfs_naive":
                    if not isRetro:
                        pathwayEnumList = self.pwnabfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnabfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  isShortest,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
            # Total pathways search.
            else:
                if methods == "bfs_naive":
                    if not isRetro:
                        pathwayEnumList = self.pwnabfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              isShortest,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnabfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  isShortest,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
                if methods == "dfs_naive":
                    if not isRetro:
                        pathwayEnumList = self.pwnadfs_instance.Source2Target(source, target, maxlength, maxtime, maxnumber,
                                                                              isTotal,
                                                                              rxnDistInfoDict, compDistInfoDict,
                                                                              isDebug = isDebug,
                                                                              )
                    else:
                        pathwayEnumList = self.pwnadfs_instance.Target2MultSource({source, }, target, maxlength, maxtime, maxnumber,
                                                                                  isTotal,
                                                                                  rxnDistInfoDict, compDistInfoDict,
                                                                                  isDebug = isDebug,
                                                                                  )
        #
        #
        #
        #
        self.ps_throw_prompts("")
        self.ps_throw_prompts("<font color='blue', font-weight='bold', size=6>Evaluation:</font>")
        #
        # Pathway evaluating and sorting.
        pathwayEnumList = self.evaluator_instance.evaluator(pathwayEnumList,
                                                            rxnDistInfoDict = rxnDistInfoDict, compDistInfoDict = compDistInfoDict,
                                                            organism = organism, org_model = org_model, orgCompRxnId_dict = orgCompRxnId_dict,
                                                            org_compId_set = org_compId_set, org_rxnId_set = org_rxnId_set,
                                                            dataBase = dataBase,
                                                            atom_mapping_file_location = atom_mapping_file_location,
                                                            isInfeasible = isInfeasible, isTransfer = isTransfer, isFlux = isFlux,
                                                            )
        #
        Save.save(pathwayEnumList)
        return [pathwayEnumList, rxnDistInfoDict, compDistInfoDict]
        #


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





