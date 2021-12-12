# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:14:32 2020

@author: CC-SXF
"""

import re
import json
from openpyxl import load_workbook

from rdkit import Chem


class InputMeta():
    """
    """
    def __init__(self):
        """ """
        pass


    @classmethod
    def _getSheetInfo(cls, excelFile, sheetName):
        """
        """
        workbook = load_workbook(excelFile)
        workbook.guess_types = True
        sheet = workbook[sheetName]
        sheetAllRows = list()
        for row in sheet.rows:
            sheetOneRow = list()
            for cell in row:
                sheetOneRow.append(cell.value)
            sheetAllRows.append(sheetOneRow)
        return sheetAllRows


    @classmethod
    def _getRxnInfoFromRxnExcel(cls, rxnExcelFile, sheetName = "InfoEx"):
        """
            Get info of every reaction's coef/coefSum/compID/canoSmiles/nonIsoSmiles and
        isoReaProReactionIDs from the (partially) cleaned reaction excel derived from data
        # base KEGG/Rhea/MetaCyc/KndPad.
        Input:
            the 'path&file' of the (partially) cleaned reaction excel derived from data base
            KEGG/Rhea/MetaCyc/KndPad ( The excle file has sheet 'InfoEx' with title:
            Entry, Names, Definition, Equation, Coefficient, Compound, SMILES,
            SMILES Equation, Atom/Bond Number, ... ).
        Output:
            [{reactionID: [reactionCoefList,reactionCoefSumList,reactionCompList,reactionCanoSmiList,
            reactionNonIsoSmiList], ...}, reactionIsoReaProIdSet].
        """
        # Analysis of input parameter.
        # Data Base Name: KEGG/Rhea/MetaCyc/KndPad
        # Input Example: "../sourceexcel/KEGG_Ori_Cleaned_Half.xlsx"
        # Read information from sheet "InfoEX" of excel file 'KEGG_Ori_Cleaned_Half.xlsx'.
        #
        InfoEx = cls._getSheetInfo(rxnExcelFile, sheetName = "InfoEx")
        #
        reactionCoefCompSmiDict = dict()
        # Reaction IDs set of which the Reactants or Products have Isomeric Compounds.
        reactionIsoReaProIdSet = set()
        for row in InfoEx[1:]:    # Exclude Title
            reactionID = row[0]
            # Parse Coefficient
            reactionCoefStr = row[4]
            reactionCoefList = json.loads(reactionCoefStr)
            reactionCoefSumList =[sum(reactionCoefList[0]), sum(reactionCoefList[1])]
            # Parse Compound
            reactionCompStr = row[5]
            reactionCompList = json.loads(reactionCompStr)
            # Parse Canonical Smiles
            reactionCanoSmiStr = row[6]
            reactionCanoSmiList = json.loads(reactionCanoSmiStr)
            # Non Isomeric Smiles
            reactionNonIsoSmiList = [[],[]]
            for canoSmi in reactionCanoSmiList[0]:
                reactionNonIsoSmiList[0].append( Chem.MolToSmiles(Chem.MolFromSmiles(canoSmi), isomericSmiles=False) )
            for canoSmi in reactionCanoSmiList[1]:
                reactionNonIsoSmiList[1].append( Chem.MolToSmiles(Chem.MolFromSmiles(canoSmi), isomericSmiles=False) )
            # Isomeric Compound In Reactants Or Products. Replace Non-Isomeric Smiles with Isomeric Smiles.
            # We will employ the same processing strategy in further process.
            if (len(set(reactionCanoSmiList[0])) != len(set(reactionNonIsoSmiList[0]))):
                print(reactionID,": Isomeric Compounds In Reactants.")
                reactionNonIsoSmiList[0] = reactionCanoSmiList[0][:]
                reactionNonIsoSmiList[1] = reactionCanoSmiList[1][:]
                reactionIsoReaProIdSet.add(reactionID)
            elif (len(set(reactionCanoSmiList[1])) != len(set(reactionNonIsoSmiList[1]))):
                print(reactionID,": Isomeric Compounds In Products.")
                reactionNonIsoSmiList[0] = reactionCanoSmiList[0][:]
                reactionNonIsoSmiList[1] = reactionCanoSmiList[1][:]
                reactionIsoReaProIdSet.add(reactionID)
            else:
                pass
            reactionCoefCompSmiDict[reactionID] = [reactionCoefList, reactionCoefSumList,
                                                   reactionCompList,
                                                   reactionCanoSmiList, reactionNonIsoSmiList]
        return [reactionCoefCompSmiDict, reactionIsoReaProIdSet]


    @classmethod
    def _getGeCosInfoFromGeCosExcel(cls, general_cofactors_xlsxFile, sheetName = "GeCos-New", isDebug = False):
        """
            Get generalized cofactors from manually established excle file.
        Input:
            general_cofactors_xlsxFile, excel file of generalized cofactors.
        Output:
            [sinFamilyList, douFamilyLL, quaFamilyLL],
            list of compounds in different families belongs to generalized cofactors.
        """
        # Analysis of input parameter.
        # Input Example: "data/KEGG_General_Cofactors.xlsx"
        # Read information from sheet 'GeCos-New' of excle file 'general_cofactors_xlsxFile'.
        #
        GeCos_Ps = cls._getSheetInfo(general_cofactors_xlsxFile, sheetName = "GeCos-New")
        #
        sinFamilyList = list()
        douFamilyList = list()
        quaFamilyList = list()
        for row in GeCos_Ps[1:]:
            compId = row[0]
            family = row[-1]
            if family == "Single":
                sinFamilyList.append(compId)
            elif family == "Double":
                douFamilyList.append(compId)
            elif family == "Quadruple":
                quaFamilyList.append(compId)
            else:
                pass
        if isDebug:
            print(sinFamilyList, "\n\n", douFamilyList, "\n\n", quaFamilyList, "\n\n")
        douFamilyLL = list()
        for i in range(0, len(douFamilyList), 2):
            douCompIds = douFamilyList[i:i+2]
            douFamilyLL.append(douCompIds)
        quaFamilyLL = list()
        for i in range(0, len(quaFamilyList), 4):
            quaCompIds = quaFamilyList[i:i+4]
            quaFamilyLL.append(quaCompIds)
        if isDebug:
            print(douFamilyLL, "\n\n", quaFamilyLL, "\n\n")
        sdqFamilyList = [sinFamilyList, douFamilyLL, quaFamilyLL]
        return sdqFamilyList


    """
    # Reactions that have difference in reactant number or product number comparing with original reaction.
    # In other words, reactions that have compound with "dot-disconnected fragments".
    # reactionIDsNumDiff = list()
    """
    @classmethod
    def _rxnFile2MolSmiComp(cls, reactionID, reactionCoefCompSmi, reactionIsoReaProIdSet, rxnFile):
        """
            Parse the rxnFile, return molBlock/addHsMapNumSmiles/compID of every compound in the reaction.
        Input:
            reactionID, reactionCoefCompSmi, reactionIsoReaProIdSet, rxnFile.
            And reactionCoefCompSmi: [reactionCoefList, reactionCoefSumList, reactionCompList,
                                      reactionCanoSmiList, reactionNonIsoSmiList]
        Output:
            [reaction_compMol List, reaction_addHsMapNumSmi List, reaction_compID List]
        """
        # Analysis of input parameter.
        reactionIdCoefCompSmiOri = reactionCoefCompSmi
        # coefReactantProductOri = reactionIdCoefCompSmiOri[0]
        numReactantProductOri = reactionIdCoefCompSmiOri[1]             # [numReactantOri, numProductOri]
        compReactantProductOri = reactionIdCoefCompSmiOri[2]            # [compIdReactantListOri, compIdProductListOri]
        # canoSmiReactantProductOri = reactionIdCoefCompSmiOri[3]
        smiReactantProductOri = reactionIdCoefCompSmiOri[4]             # [nonIsoSmilesReactantListOri, nonIsoSmilesProductListOri]
        reactantSmiCompDict = dict(zip(smiReactantProductOri[0], compReactantProductOri[0]))        # {nonIsoSmilesReactantOri: compIdReactantOri}
        productSmiCompDict = dict(zip(smiReactantProductOri[1], compReactantProductOri[1]))         # {nonIsoSmilesProductOri: compIdProductOri}
        #
        with open(rxnFile, "r") as rxnFile_Obj:
            rxnFileContent = rxnFile_Obj.read()
        #
        # Number of reactant and product.
        pattern_numReactantProduct = re.compile("[$]RXN\n\n  EC-BLAST     smiles\n\n\s*(\d+)\s+(\d+)\n[$]MOL\n")
        result_numReactantProduct = re.search(pattern_numReactantProduct, rxnFileContent)
        numReactantProductStrTuple = result_numReactantProduct.groups()
        numReactant = int(numReactantProductStrTuple[0])
        numProduct = int(numReactantProductStrTuple[1])
        numReactantProduct = [numReactant, numProduct]
        #
        # Reactions that have difference in reactant number or product number comparing with original reaction.
        # In other words, reactions that have compound with "dot-disconnected fragments".
        # Example(Reactions from KEGG):
        # R02480 : [1, 2] vs [1, 3]  &&  [[],[(M00002,M00002)]],
        # R04641 : [2, 2] vs [3, 3]  &&  [[(M00001,M00002)],[(M00001,M00004)]],
        # R05182 : [6, 5] vs [8, 7]  &&  [[(M00003,M00004),(M00003,M00005)],[(M00003,M00008),(M00003,M00009)]]
        """
        global reactionIDsNumDiff
        if numReactantProductOri != numReactantProduct:
            reactionIDsNumDiff.append(reactionID)
            # print(reactionID, ":", numReactantProductOri, "vs",  numReactantProduct)
        if reactionID == "R12422":
            print(reactionIDsNumDiff, "\n")
        """
        # # reactionIDsNumDiff = ["R02480", "R04641", "R05182"]
        #
        # MolBlock of reactant and product
        pattern_compMol = re.compile("[$]MOL\n(M\d{5}.*?M  END\n)", re.S)
        result_compMol = re.findall(pattern_compMol, rxnFileContent)
        reactant_compMol = result_compMol[:numReactant]
        product_compMol = result_compMol[-numProduct:]
        reaction_compMol = [reactant_compMol, product_compMol]
        #
        # addHsMapNumSmiles of reactant and product
        reactant_addHsMapNumSmi = []
        for compMol in reactant_compMol:
            reactant_addHsMapNumSmi.append(Chem.MolToSmiles(Chem.MolFromMolBlock(compMol, removeHs=False)))
        product_addHsMapNumSmi = []
        for compMol in product_compMol:
            product_addHsMapNumSmi.append(Chem.MolToSmiles(Chem.MolFromMolBlock(compMol, removeHs=False)))
        reaction_addHsMapNumSmi = [reactant_addHsMapNumSmi, product_addHsMapNumSmi]
        #
        # Remove atom-atom mapping number of every atom.
        pattern_atomMapNum = re.compile(":\d+")
        reactant_addHsNoMapNumSmi = []
        for addHsMapNumSmi in reactant_addHsMapNumSmi:
            reactant_addHsNoMapNumSmi.append(re.sub(pattern_atomMapNum, "", addHsMapNumSmi))
        product_addHsNoMapNumSmi = []
        for addHsMapNumSmi in product_addHsMapNumSmi:
            product_addHsNoMapNumSmi.append(re.sub(pattern_atomMapNum, "", addHsMapNumSmi))
        # Canocial Smiles
        reactant_canoSmi = []
        for addHsSmi in reactant_addHsNoMapNumSmi:
            reactant_canoSmi.append(Chem.MolToSmiles(Chem.MolFromSmiles(addHsSmi)))
        product_canoSmi = []
        for addHsSmi in product_addHsNoMapNumSmi:
            product_canoSmi.append(Chem.MolToSmiles(Chem.MolFromSmiles(addHsSmi)))
        # Non Isomeric Smiles
        reactant_nonIsoSmi = []
        product_nonIsoSmi = []
        # Isomeric Compound In Reactants Or Products.
        # Replace Non-Isomeric Smiles with Isomeric Smiles.
        # We have employed the same processing strategy in previous process.
        if reactionID in reactionIsoReaProIdSet:
            for canoSmi in reactant_canoSmi:
                reactant_nonIsoSmi.append( canoSmi )
            for canoSmi in product_canoSmi:
                product_nonIsoSmi.append( canoSmi )
        # Isomeric Smiles --> Non Isomeric Smiles
        else:
            for canoSmi in reactant_canoSmi:
                reactant_nonIsoSmi.append( Chem.MolToSmiles(Chem.MolFromSmiles(canoSmi), isomericSmiles=False) )
            for canoSmi in product_canoSmi:
                product_nonIsoSmi.append( Chem.MolToSmiles(Chem.MolFromSmiles(canoSmi), isomericSmiles=False) )
        # Non Isomeric Canonical Smiles
        reactant_canoNonIsoSmiles = []
        for nonIsoSmi in reactant_nonIsoSmi:
            reactant_canoNonIsoSmiles.append( Chem.MolToSmiles(Chem.MolFromSmiles(nonIsoSmi)) )
        product_canoNonIsoSmiles = []
        for nonIsoSmi in product_nonIsoSmi:
            product_canoNonIsoSmiles.append( Chem.MolToSmiles(Chem.MolFromSmiles(nonIsoSmi)) )
        #
        reactant_comp = []
        product_comp = []
        # Reactions that have the same reactant number and product number comparing with original reaction.
        if numReactantProductOri == numReactantProduct:
            for smiles in reactant_canoNonIsoSmiles:
                reactant_comp.append(reactantSmiCompDict[smiles])
            for smiles in product_canoNonIsoSmiles:
                product_comp.append(productSmiCompDict[smiles])
        # Reactions that have difference in reactant number or product number comparing with original reaction.
        # In other words, reactions that have compound with "dot-disconnected fragments".
        # Example(Reactions from KEGG): R02480,R04641,R05182
        else:
            reactantSmiCompDictKeysList = list(reactantSmiCompDict.keys())
            for smiles in reactant_canoNonIsoSmiles:
                try:
                    reactant_comp.append(reactantSmiCompDict[smiles])
                except KeyError:
                    for smiCompKey in reactantSmiCompDictKeysList:
                        if smiles in smiCompKey.split("."):
                            reactant_comp.append(reactantSmiCompDict[smiCompKey])
                            break
            productSmiCompDictKeysList = list(productSmiCompDict.keys())
            for smiles in product_canoNonIsoSmiles:
                try:
                    product_comp.append(productSmiCompDict[smiles])
                except KeyError:
                    for simCompKey in productSmiCompDictKeysList:
                        if smiles in simCompKey.split("."):
                            product_comp.append(productSmiCompDict[simCompKey])
                            break
        reaction_compID = [reactant_comp, product_comp]
        return [reaction_compMol, reaction_addHsMapNumSmi, reaction_compID]


    @classmethod
    def _getCompCands(self, rxnCompIdList, sdqFamilyList, isDebug = True):
        """
            Get compound candidates of one reaction for pathway search, excluding manually
        generalized cofactors.
        Input:
            rxnCompIdList, compounds list of one reaction;
            sdqFamilyList, compounds list in different families belongs to manually
        generalized cofactors;
        Output:
            usaRxnCompIdList, in the form of  [reaCompIdList,proCompIdList].
            Compounds abandoned in pathway search will be set value "######".
        """
        # Analysis of input parameter.
        reaCompIdList = rxnCompIdList[0][:]
        proCompIdList = rxnCompIdList[1][:]
        [sinFamilyList, douFamilyLL, quaFamilyLL] = sdqFamilyList
        reaCompIdListLength = len(reaCompIdList)
        proCompIdListLength = len(proCompIdList)
        nonCompIdSign = "######"
        #
        if isDebug:
            print("Compound Candidates:")
            print(reaCompIdList, " --> ", proCompIdList)
        #
        # Double family
        for i in range(reaCompIdListLength):
            reaCompId = reaCompIdList[i]
            # Has been analysed
            if reaCompId == nonCompIdSign:
                continue
            #
            reaCompIdInDoubleFlag = False
            for douFamilyList in douFamilyLL:
                if reaCompId in douFamilyList:
                    reaCompIdInDoubleFlag = True
                    proCompIdCandidate = douFamilyList[1-douFamilyList.index(reaCompId)]
                    break
            if reaCompIdInDoubleFlag == True:
                for j in range(proCompIdListLength):
                    proCompId = proCompIdList[j]
                    if proCompId == proCompIdCandidate:
                        reaCompIdList[i] = nonCompIdSign
                        proCompIdList[j] = nonCompIdSign
                        break
        #
        # Quadruple family
        for i in range(reaCompIdListLength):
            reaCompId = reaCompIdList[i]
            # Has been analysed
            if reaCompId == nonCompIdSign:
                continue
            #
            reaCompIdInQuadrupleFlag = False
            for quaFamilyList in quaFamilyLL:
                if reaCompId in quaFamilyList:
                    reaCompIdInQuadrupleFlag = True
                    quaFamilyListCandidate = quaFamilyList[::]
                    break
            reaQuaNum = 0
            proQuaNum = 0
            if reaCompIdInQuadrupleFlag == True:
                for k in range(reaCompIdListLength):
                    reaCompId = reaCompIdList[k]
                    if reaCompId in quaFamilyListCandidate:
                        reaQuaNum += 1
                for j in range(proCompIdListLength):
                    proCompId = proCompIdList[j]
                    if proCompId in quaFamilyListCandidate:
                        proQuaNum += 1
                if reaQuaNum == proQuaNum:
                    for k in range(reaCompIdListLength):
                        reaCompId = reaCompIdList[k]
                        if reaCompId in quaFamilyListCandidate:
                            reaCompIdList[k] = nonCompIdSign
                    for j in range(proCompIdListLength):
                        proCompId = proCompIdList[j]
                        if proCompId in quaFamilyListCandidate:
                            proCompIdList[j] = nonCompIdSign
        #
        # Single family
        for i in range(reaCompIdListLength):
            reaCompId = reaCompIdList[i]
            if reaCompId in sinFamilyList:
                reaCompIdList[i] = nonCompIdSign
        # Single family
        for j in range(proCompIdListLength):
            proCompId = proCompIdList[j]
            if proCompId in sinFamilyList:
                proCompIdList[j] = nonCompIdSign
        #
        if isDebug:
            print(reaCompIdList, " --> ", proCompIdList)
            print("")
        #
        usaRxnCompIdList = [reaCompIdList, proCompIdList]
        return usaRxnCompIdList


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





