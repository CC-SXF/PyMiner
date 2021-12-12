# -*- coding: utf-8 -*-
"""
Created on Wed May 26 16:43:15 2021

@author: CC-SXF
"""

import re
# import sys
import pickle
from os import path, mkdir
from rdkit import Chem

import threading
import inspect
import ctypes
from datetime import datetime

# import math
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
# # from concurrent.futures import ProcessPoolExecutor
# from concurrent.futures import as_completed

from pathwaysearch.fingerprint.fingerprint import FingerPrint
from pathwaysearch.fingerprint.fingerprintinfo import FingerprintInfo
from pathwaysearch.evaluator.drawtool import DrawTool
from pathwaysearch.input.input_meta import InputMeta


class AtomMap():
    """
    """
    def __init__(self, ):
        """ """
        pass


    @classmethod
    def _subMatch(cls, mol = None, patt = None, single_match = True, atom_idx_dict = dict()):
        """
        """
        if single_match == True:
            atom_idx_tuple = mol.GetSubstructMatch(patt, useChirality=False)
            atom_idx_dict["sub-match"] = atom_idx_tuple
        else:
            atom_idx_tuple_tuple = mol.GetSubstructMatches(patt, useChirality=False)
            atom_idx_dict["sub-match"] = atom_idx_tuple_tuple
        return atom_idx_dict


    @classmethod
    def _async_raise(cls, thr_id, exc_type):
        """
            Raise the exception, perform cleanup if needed.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        thr_id = ctypes.c_long(thr_id)
        if not inspect.isclass(exc_type):
            exc_type = type(exc_type)
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thr_id, ctypes.py_object(exc_type))
        if res == 0:
            raise ValueError("Invalid thread id.")
        elif res != 1:
            # If it returns a number greater than one, you're in trouble,
            # and you should call it again with exc=NULL to revert the effect
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thr_id, None)
            raise SystemError("PyThreadState_SetAsyncExc failed.")


    @classmethod
    def _stop_thread(cls, thr):
        """
        """
        cls._async_raise(thr.ident, SystemExit)


    @classmethod
    def _getSubMatch(cls, mol = None, patt = None, single_match = True, time_threshold = 60):
        """
        """
        atom_idx_dict = dict()
        atom_idx_dict["sub-match"] = None
        thr = threading.Thread( target = cls._subMatch, args = (mol, patt, single_match, atom_idx_dict) )
        thr.daemon = True
        thr.start()
        #
        startTime = datetime.now()
        while True:
            if atom_idx_dict["sub-match"] != None:
                break
            # Calculate the running time.
            endTime = datetime.now()
            runTime = (endTime - startTime).total_seconds()
            if runTime > time_threshold:
                cls._stop_thread(thr)
                atom_idx_dict["sub-match"] = tuple()
                break
        #
        return atom_idx_dict["sub-match"]


    @classmethod
    def _getExMapNum(cls, compMolBlock, dataBase = "KEGG"):
        """
            Get the atom-atom mapping numbers of the atoms that belong to several substructures,
        including substructure Pi, PPi, and CoA-SH.
            And we will exclude these atoms in atom transfer analysis.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        exSubStrTupleList = FingerPrint.getExSubStrInfo(dataBase = dataBase)
        # exSubStrTupleList: [(subMol, rdkitSubBit, morganSubBit), ..., ]
        #
        compMol = Chem.MolFromMolBlock(compMolBlock, sanitize=True, removeHs=False, strictParsing=True)
        rdkitBit = FingerprintInfo.getRDKitFingerprintInfo(compMol)
        morganBit = FingerprintInfo.getMorganFingerprintInfo(compMol)
        #
        exclude_mapnum_set = set()
        for exSubStrTuple in exSubStrTupleList:
            subMol = exSubStrTuple[0]
            rdkitSubBit = exSubStrTuple[1]
            morganSubBit = exSubStrTuple[2]
            isRdkit = (set(rdkitBit) & set(rdkitSubBit))  == set(rdkitSubBit)
            isMorgan = (set(morganBit) & set(morganSubBit)) == set(morganSubBit)
            if isRdkit and isMorgan:
                # atomIdx_tuple_tuple = compMol.GetSubstructMatches(subMol, useChirality=False)
                atomIdx_tuple_tuple = cls._getSubMatch(mol = compMol, patt = subMol, single_match = False)
                for atomIdx_tuple in atomIdx_tuple_tuple:
                    for atomIdx in atomIdx_tuple:
                        atom = compMol.GetAtomWithIdx(atomIdx)
                        map_num = atom.GetAtomMapNum()
                        exclude_mapnum_set.add(map_num)
        #
        return exclude_mapnum_set


    @classmethod
    def _updateTransferAtom(cls, compMolBlock, middle_structure, mapnum_atomsymbol_dict,):
        """
            Update the identical atoms between the current reactant and product.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        mol = Chem.MolFromMolBlock(compMolBlock, sanitize=True, removeHs=False, strictParsing=True)
        #
        # >>> m.GetSubstructMatch(patt)
        # (i, j, k, )
        # Those are the atom indies in m, ordered as patt's atoms.
        # atom_idx_tuple  = mol.GetSubstructMatch(middle_structure, useChirality=False)
        atom_idx_tuple = cls._getSubMatch(mol = mol, patt = middle_structure, single_match = True)
        #
        new_mapnum_set = set()
        for idx, atom_idx in enumerate(atom_idx_tuple):
            last_map_num = middle_structure.GetAtomWithIdx(idx).GetAtomMapNum()
            if last_map_num != 0:
                new_map_num = mol.GetAtomWithIdx(atom_idx).GetAtomMapNum()
                new_mapnum_set.add(new_map_num)
        #
        mapnum_atomsymbol_dict_new = dict()
        for map_num in new_mapnum_set:
            try:
                atom_symbol = mapnum_atomsymbol_dict[map_num]
                mapnum_atomsymbol_dict_new[map_num] = atom_symbol
            except KeyError:
                pass
        #
        return mapnum_atomsymbol_dict_new


    @classmethod
    def generateMapFigure(cls, compId_left, compIdList_left, compMolBlockList_left,
                               compId_right, compIdList_right, compMolBlockList_right,
                               dataBase = "KEGG",
                               middle_structure = None,
                               source_target_atom_number_dict = {"source":1, "target":1},
                               sub_route = "",
                               ):
        """
            Analyse the transfer of non-hydrogen atom. And pathway of which no carbon atom transfered from source
        to target will be consider as invalid one.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        if source_target_atom_number_dict == dict():
            isSource = True
        else:
            isSource = False
        #
        middle_structure_next = None
        identical_atom_number = 0
        isSubstructMatchError = True
        #
        #
        # The Counts Line
        #                                     "number of atoms" "number of bonds"
        # pattern_molAtomBondNum = re.compile("([ \d]{3})([ \d]{3}).*V2000\n")
        # The Atom Block
        #                                                  "atom symbol"            "charge"  "stereo"                                                  "AAM Num"
        pattern_molAtomInfo = re.compile(" *\S+ +\S+ +\S+ +([a-zA-Z\*]+) +[\d]{0,3}([ \d]{3})([ \d]{3})[ \d]{3}[ \d]{3}[ \d]{3}[ \d]{3}[ \d]{3}[ \d]{3}([ \d]{3})[ \d]{3}[ \d]{3}\n")
        # The Bond Block
        #                                   "first atom number" "second atom number" "bond type" "bond stereo"
        # pattern_molBondInfo = re.compile("\n([ \d]{3})([ \d]{3})([ \d]{3})([ \d]{3})[ \d]{3}[ \d]{3}[ \d]{3} ")
        #
        #
        # statistics of non-hydrogen atom
        mapnum_atomsymbol_dict_list_left = list()
        for compId, compMolBlock in zip(compIdList_left, compMolBlockList_left):
            mapnum_atomsymbol_dict = dict()
            if compId == compId_left:
                atom_info_tuple_list = re.findall(pattern_molAtomInfo, compMolBlock)
                # * * * * * *
                # exclude atoms that belong to some substructures(e.g., Pi, PPi, CoA-SH)
                exclude_mapnum_set = cls._getExMapNum(compMolBlock, dataBase = dataBase)
                for atom_info_tuple in atom_info_tuple_list:
                    atom_symbol = atom_info_tuple[0]
                    map_num = int(atom_info_tuple[-1])
                    if (atom_symbol != 'H') and (map_num not in exclude_mapnum_set):
                        mapnum_atomsymbol_dict[map_num] = atom_symbol
                # source of this pathway
                if isSource:
                    # # source_target_atom_number_dict["source"] = len(mapnum_atomsymbol_dict)
                    pass
                else:
                    # * * * * * * * * * * * *
                    # Update the identical atoms between this reactant and product according to the conserved substructures generated in last cycle of analysis.
                    mapnum_atomsymbol_dict = cls._updateTransferAtom(compMolBlock,
                                                                     middle_structure,
                                                                     mapnum_atomsymbol_dict,
                                                                     )
            #
            mapnum_atomsymbol_dict_list_left.append(mapnum_atomsymbol_dict)
        #
        #
        # * * * * * *
        # Error occurs while updating the identical atoms...
        for mapnum_atomsymbol_dict in mapnum_atomsymbol_dict_list_left:
            if len(mapnum_atomsymbol_dict) != 0:
                isSubstructMatchError = False
                break
        if isSubstructMatchError:
            return [middle_structure_next, identical_atom_number, isSubstructMatchError]
        #
        #
        # statistics of non-hydrogen atom
        mapnum_atomsymbol_dict_list_right = list()
        for compId, compMolBlock in zip(compIdList_right, compMolBlockList_right):
            mapnum_atomsymbol_dict = dict()
            if compId == compId_right:
                atom_info_tuple_list = re.findall(pattern_molAtomInfo, compMolBlock)
                # * * * * * *
                # exclude atoms that belong to some substructures(e.g., Pi, PPi, CoA-SH)
                exclude_mapnum_set = cls._getExMapNum(compMolBlock, dataBase = dataBase)
                for atom_info_tuple in atom_info_tuple_list:
                    atom_symbol = atom_info_tuple[0]
                    map_num = int(atom_info_tuple[-1])
                    if (atom_symbol != 'H') and (map_num not in exclude_mapnum_set):
                        mapnum_atomsymbol_dict[map_num] = atom_symbol
                # target or intermediate metabolite of this pathway
                # # source_target_atom_number_dict["target"] = len(mapnum_atomsymbol_dict)
            #
            mapnum_atomsymbol_dict_list_right.append(mapnum_atomsymbol_dict)
        #
        #
        mapnum_set_list_left = [set(item.keys()) for item in mapnum_atomsymbol_dict_list_left]
        mapnum_set_list_right = [set(item.keys()) for item in mapnum_atomsymbol_dict_list_right]
        #
        # the maximum number of identical atom-atom mapping number of non-hydrogen atom
        # and we will exclude atoms that belong to some substructures(e.g., Pi, PPi, CoA-SH).
        idx_max_mapnum_list = [0, 0, 0, set()]
        idx_left = 0
        for compId_l, mapnum_set_l in zip(compIdList_left, mapnum_set_list_left):
            if compId_l != compId_left:
                idx_left += 1
                continue
            idx_right = 0
            for compId_r, mapnum_set_r in zip(compIdList_right, mapnum_set_list_right):
                if compId_r != compId_right:
                    idx_right += 1
                    continue
                identical_mapnum_set = (mapnum_set_l & mapnum_set_r)
                identical_atom_number = len(identical_mapnum_set)
                if identical_atom_number > idx_max_mapnum_list[2]:
                    idx_max_mapnum_list = [idx_left, idx_right, identical_atom_number, identical_mapnum_set]
                idx_right += 1
            idx_left += 1
        #
        idx_reactant = idx_max_mapnum_list[0]
        idx_product = idx_max_mapnum_list[1]
        identical_atom_number = idx_max_mapnum_list[2]
        identical_mapnum_set = idx_max_mapnum_list[3]
        mapnum_atomsymbol_dict_left = mapnum_atomsymbol_dict_list_left[idx_reactant]
        mapnum_atomsymbol_dict_right = mapnum_atomsymbol_dict_list_right[idx_product]
        #
        #
        # * * * * * *
        # analyse the species of identical transfered atoms between reactant and product
        atom_species_set = set()
        for map_num in identical_mapnum_set:
            atom_symbol = mapnum_atomsymbol_dict_right[map_num]
            atom_species_set.add(atom_symbol)
        # the identical transfered atoms doesn't include carbon atom
        if "C" not in atom_species_set:
            return [middle_structure_next, identical_atom_number, isSubstructMatchError]
        #
        #
        # source of this pathway
        if isSource:
            source_target_atom_number_dict["source"] = len(mapnum_atomsymbol_dict_left)
        # target or intermediate metabolite of this pathway
        source_target_atom_number_dict["target"] = len(mapnum_atomsymbol_dict_right)
        #
        #
        compMolBlock = compMolBlockList_left[idx_reactant]
        compMol_left = Chem.MolFromMolBlock(compMolBlock, sanitize=True, removeHs=False, strictParsing=True)
        compMolBlock = compMolBlockList_right[idx_product]
        compMol_right = Chem.MolFromMolBlock(compMolBlock, sanitize=True, removeHs=False, strictParsing=True)
        #
        # source of this pathway
        if isSource:
            DrawTool.drawMoleculeRCenter(compMol_left, compId_left, list(identical_mapnum_set), sub_route = sub_route)
        # target or intermediate metabolite of this pathway
        DrawTool.drawMoleculeRCenter(compMol_right, compId_right, list(identical_mapnum_set), sub_route = sub_route)
        #
        #
        # * * * * * * * * * * * *
        # the atom-atom mapping number of unconserved will be set zero.
        for idx in range(compMol_right.GetNumAtoms()):
            atom_right = compMol_right.GetAtomWithIdx(idx)
            mapnum_right = atom_right.GetAtomMapNum()
            if mapnum_right not in identical_mapnum_set:
                atom_right.SetAtomMapNum(0)
        middle_structure_next = compMol_right
        #
        #
        return [middle_structure_next, identical_atom_number, isSubstructMatchError]
    #


    @classmethod
    def _evalOneAtTr(cls, pathway_info,
                          rxnCoefCompSmiDict, rxnIsoReaProIdSet,
                          dataBase = "KEGG",
                          atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/",
                          isDebug = False,
                          ):
        """
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
        # length = pathway_info[0]
        pathway = pathway_info[8]
        sub_route = pathway_info[11]
        comp2rxn_list = pathway.split('-->')
        # # source = comp2rxn_list[0]
        target = comp2rxn_list[-1]
        #
        # initialize
        middle_structure = None
        source_target_atom_number_dict = dict()
        #
        for idx_cr, comp2rxn in enumerate(comp2rxn_list):
            # compound
            if idx_cr%2 == 0:
                continue
            # reaction
            rxn_id = comp2rxn
            substrate_id_left = comp2rxn_list[idx_cr-1]
            substrate_id_right = comp2rxn_list[idx_cr+1]
            #
            rxnFile = "".join([rxn_id, "_ECBLAST_smiles_AAM.rxn"])
            rxnFile = path.join(atom_mapping_file_location, rxnFile)
            rxnCoefCompSmi = rxnCoefCompSmiDict[rxn_id]
            try:
                #
                rxnInputInfo = InputMeta._rxnFile2MolSmiComp(rxn_id, rxnCoefCompSmi, rxnIsoReaProIdSet, rxnFile)
                rxnCompMolBlockList = rxnInputInfo[0]
                # rxnAddHsMnSmilesList = rxnInputInfo[1]
                rxnCompIdList = rxnInputInfo[2]
                #
                reaCompMolBlockList = rxnCompMolBlockList[0]
                proCompMolBlockList = rxnCompMolBlockList[1]
                # reaAddHsMnSmilesList = rxnAddHsMnSmilesList[0]
                # proAddHsMnSmilesList = rxnAddHsMnSmilesList[1]
                reaCompIdList = rxnCompIdList[0]
                proCompIdList = rxnCompIdList[1]
                #
                # left --> right
                if (substrate_id_left in reaCompIdList) and (substrate_id_right in proCompIdList):
                    substructure_info_list = cls.generateMapFigure(substrate_id_left, reaCompIdList, reaCompMolBlockList,
                                                                   substrate_id_right, proCompIdList, proCompMolBlockList,
                                                                   dataBase = dataBase,
                                                                   middle_structure = middle_structure,
                                                                   source_target_atom_number_dict = source_target_atom_number_dict,
                                                                   sub_route = sub_route,
                                                                   )
                    middle_structure_next = substructure_info_list[0]
                    identical_atom_number = substructure_info_list[1]
                    isSubstructMatchError = substructure_info_list[2]
                    # There are no valid atoms transfered from source to target.
                    if (middle_structure_next == None) and (not isSubstructMatchError):
                        idx_r =(idx_cr-1)//2
                        pathway_info[9][idx_r] = -1
                        pathway_info[10][idx_r].append(rxn_id)
                        pathway_info[0] = False
                        pathway_info[2] = pathway_info[9].count(1)          # Endogenous Pathway Length
                        pathway_info[3] = pathway_info[9].count(0)          # Heterogeneous Pathway Length
                        pathway_info[4] = pathway_info[9].count(-1)         # Infeasible Pathway Length
                        pathway_info[5] = 0                                 # Atom Utilization
                        pathway_info[6] = 0                                 # Atom Conservation
                        break
                    # Error occurs while matching one compound's substruture.
                    elif (middle_structure_next == None) and (isSubstructMatchError):
                        break
                    # target is reached.
                    elif substrate_id_right == target:
                        # Calculate atom utilization and atom conservation of this pathway.
                        source_atom_number = source_target_atom_number_dict['source']
                        target_atom_number = source_target_atom_number_dict['target']
                        pathway_info[5] = round(identical_atom_number/source_atom_number, 2)  # Atom Utilization
                        pathway_info[6] = round(identical_atom_number/target_atom_number, 2)  # Atom Conservation
                        break
                    # update...
                    else:
                        middle_structure = middle_structure_next
                #
                # right --> left
                elif (substrate_id_left in proCompIdList) and (substrate_id_right in reaCompIdList):
                    substructure_info_list = cls.generateMapFigure(substrate_id_left, proCompIdList, proCompMolBlockList,
                                                                   substrate_id_right, reaCompIdList, reaCompMolBlockList,
                                                                   dataBase = dataBase,
                                                                   middle_structure = middle_structure,
                                                                   source_target_atom_number_dict = source_target_atom_number_dict,
                                                                   sub_route = sub_route,
                                                                   )
                    middle_structure_next = substructure_info_list[0]
                    identical_atom_number = substructure_info_list[1]
                    isSubstructMatchError = substructure_info_list[2]
                    # There are no valid atoms transfered from source to target.
                    if (middle_structure_next == None) and (not isSubstructMatchError):
                        idx_r =(idx_cr-1)//2
                        pathway_info[9][idx_r] = -1
                        pathway_info[10][idx_r].append(rxn_id)
                        pathway_info[0] = False
                        pathway_info[2] = pathway_info[9].count(1)          # Endogenous Pathway Length
                        pathway_info[3] = pathway_info[9].count(0)          # Heterogeneous Pathway Length
                        pathway_info[4] = pathway_info[9].count(-1)         # Infeasible Pathway Length
                        pathway_info[5] = 0                                 # Atom Utilization
                        pathway_info[6] = 0                                 # Atom Conservation
                        break
                    # Error occurs while matching one compound's substruture.
                    elif (middle_structure_next == None) and (isSubstructMatchError):
                        break
                    # target is reached.
                    elif substrate_id_right == target:
                        # Calculate atom utilization and atom conservation of this pathway.
                        source_atom_number = source_target_atom_number_dict['source']
                        target_atom_number = source_target_atom_number_dict['target']
                        pathway_info[5] = round(identical_atom_number/source_atom_number, 2)  # Atom Utilization
                        pathway_info[6] = round(identical_atom_number/target_atom_number, 2)  # Atom Conservation
                        break
                    # update...
                    else:
                        middle_structure = middle_structure_next
                else:
                    # MetaCyc: PWY-6519
                    # MetaCyc: RXN-11474
                    # MetaCyc: Malonyl-acp-methyl-ester / MALONYL-ACP
                    print(substrate_id_left, ",", substrate_id_right)
                    print(reaCompIdList)
                    print(proCompIdList)
                    # sys.exit("Atom transfer evaluation...")
            except FileNotFoundError:
                break
            except Exception:
                # rdkit: ArgumentError
                # MetaCyc: CPD-19511
                # Chem.MolToSmiles(Chem.MolFromMolBlock(mol block, removeHs=False))
                break
        #
        return pathway_info


    @classmethod
    def _evalAtomTransfer(cls, pathway_list,
                               dataBase = "KEGG",
                               atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/",
                               isDebug = False,
                               ):
        """
            Calculate the details of atom transfer within the pathway from source to target.
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
        info_data = "".join([dataBase, "_Info.data"])
        info_data = path.join("temp", info_data)
        with open(info_data, "rb") as file_obj:
            [rxnCoefCompSmiDict, rxnIsoReaProIdSet] = pickle.load(file_obj)
        #
        for idx_p, pathway_info in enumerate(pathway_list):
            # Append the route of this pathway's highlighted compound structure figures
            sub_route = "".join(["ps_results_", str(idx_p+1).rjust(6, '0')])
            sub_route = path.join("results", sub_route)
            mkdir(sub_route)
            pathway_info[11] = sub_route
        #
        # """
        for pathway_info in pathway_list:
            cls._evalOneAtTr(pathway_info,
                             rxnCoefCompSmiDict, rxnIsoReaProIdSet,
                             dataBase = dataBase,
                             atom_mapping_file_location = atom_mapping_file_location,
                             isDebug = isDebug,
                             )
        # """
        #
        """
        # parallelism --> concurrency
        max_workers = multiprocessing.cpu_count()
        max_workers = math.ceil(max_workers * 0.8)
        # Use a "with statement" to ensure threads are cleaned up promptly
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Start the load operations and mark each future with pathway's string.
            futures = {executor.submit(cls._evalOneAtTr, pathway_info,
                                                         rxnCoefCompSmiDict, rxnIsoReaProIdSet,
                                                         dataBase,
                                                         atom_mapping_file_location,
                                                         isDebug):pathway_info[8] for pathway_info in pathway_list}
            for future in as_completed(futures):
                s_details = futures[future]
                try:
                    pathway_info = future.result()
                except Exception as excReason:
                    print("\nWarning: \nPathway: %s generated an exception: %s\n" % (s_details, excReason))
        # """
        #
        return pathway_list


    @classmethod
    def demoFunc(cls, ):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





