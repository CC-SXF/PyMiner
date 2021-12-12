# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:44:59 2021

@author: CC-SXF
"""

import re
from rdkit import Chem
from pathwaysearch.fingerprint.fingerprintinfo import FingerprintInfo


class ReaProPair():
    """
        Get Reactant/Product Pairs for Pathway Search.
    """

    def __init__(self,):
        """ """
        pass


    @classmethod
    def _exMapNum(cls, exSubStrTupleList, compMolBlock, compAddHsSmiles, isDebug=False):
        """
            Set the atom-atom mapping numbers of the atoms that belong to some avoided substructures
          of this molecule to 0.

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        # Analysis of input parameter.
        # exSubStrTupleList: [(subMol, rdkitSubBit, morganSubBit), ..., ]
        #
        compMol = Chem.MolFromMolBlock(compMolBlock, removeHs=False)
        rdkitBit = FingerprintInfo.getRDKitFingerprintInfo(compMol)
        morganBit = FingerprintInfo.getMorganFingerprintInfo(compMol)
        #
        exMapNumList = list()
        for exSubStrTuple in exSubStrTupleList:
            subMol = exSubStrTuple[0]
            rdkitSubBit = exSubStrTuple[1]
            morganSubBit = exSubStrTuple[2]
            isRdkit = (set(rdkitBit) & set(rdkitSubBit))  == set(rdkitSubBit)
            isMorgan = (set(morganBit) & set(morganSubBit)) == set(morganSubBit)
            if isRdkit and isMorgan:
                atomIdx_tuple_tuple = compMol.GetSubstructMatches(subMol, useChirality=False)
                for atomIdx_tuple in atomIdx_tuple_tuple:
                    for atomIdx in atomIdx_tuple:
                        atom = compMol.GetAtomWithIdx(atomIdx)
                        mapNum = atom.GetAtomMapNum()
                        exMapNumList.append(mapNum)
        #
        mapNum_zero = ':0]'
        for mapNum in exMapNumList:
            mapNum = "".join([':', str(mapNum), ']'])
            compAddHsSmiles = compAddHsSmiles.replace(mapNum, mapNum_zero)
        #
        if isDebug:
            print(compAddHsSmiles)
        #
        return compAddHsSmiles


    @classmethod
    def _exSubStr(cls, usaRxnCompIdList, exSubStrTupleList,
                       rxnCompMolBlockList, rxnAddHsMnSmilesList):
        """
            Set the atom-atom mapping numbers of the atoms that belong to some avoided substructures
          of all molecules in this reaction to 0.
            We will exclude the information of these atoms in generating compound pairs.
            # Example: CoA-SH; PPi; Pi;

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        # Analysis of input parameter.
        usaReaCompIdList = usaRxnCompIdList[0]
        usaProCompIdList = usaRxnCompIdList[1]
        reaCompMolBlockList = rxnCompMolBlockList[0]
        proCompMolBlockList = rxnCompMolBlockList[1]
        reaAddHsMnSmilesList = rxnAddHsMnSmilesList[0]
        proAddHsMnSmilesList = rxnAddHsMnSmilesList[1]
        #
        for idx, compId in enumerate(usaReaCompIdList):
            if compId == "######":
                continue
            else:
                compMolBlock = reaCompMolBlockList[idx]
                compAddHsSmiles = reaAddHsMnSmilesList[idx]
                compAddHsSmiles = cls._exMapNum(exSubStrTupleList, compMolBlock, compAddHsSmiles)
                reaAddHsMnSmilesList[idx] = compAddHsSmiles
        #
        for idx, compId in enumerate(usaProCompIdList):
            if compId == "######":
                continue
            else:
                compMolBlock = proCompMolBlockList[idx]
                compAddHsSmiles = proAddHsMnSmilesList[idx]
                compAddHsSmiles = cls._exMapNum(exSubStrTupleList, compMolBlock, compAddHsSmiles)
                proAddHsMnSmilesList[idx] = compAddHsSmiles
        #
        rxnAddHsMnSmilesList = [reaAddHsMnSmilesList, proAddHsMnSmilesList]
        return rxnAddHsMnSmilesList


    @classmethod
    def _getRxnMnAtSyInfo(cls, usaRxnCompIdList, exSubStrTupleList,
                               rxnCompMolBlockList, rxnAddHsMnSmilesList):
        """
            Get the information of all molecules(exclude cofactors) in this reaction, including
          atom-atom mapping number, atom symbol, atom number, and compound ID.

            ARGUMENTS:
                ...

            RETURNS:
                [
                  [mapping-number set, atom-symbol dict, atom-number int, compound-id str], ..., ],  # Molecules of reactants
                  [mapping-number set, atom-symbol dict, atom-number int, compound-id str], ..., ]   # Molecules of products
                ]
        """
        # Analysis of input parameter.
        usaReaCompIdList = usaRxnCompIdList[0]
        usaProCompIdList = usaRxnCompIdList[1]
        # reaCompMolBlockList = rxnCompMolBlockList[0]
        # proCompMolBlockList = rxnCompMolBlockList[1]
        # reaAddHsMnSmilesList = rxnAddHsMnSmilesList[0]
        # proAddHsMnSmilesList = rxnAddHsMnSmilesList[1]
        #
        # Set the atom mapping number of one atom that belong to some avoided substructure of one compound to 0.
        # We will exclude the atoms from avoided substructures.
        # Example: CoA-SH; PPi; Pi
        rxnAddHsMnSmilesList = cls._exSubStr(usaRxnCompIdList, exSubStrTupleList, rxnCompMolBlockList, rxnAddHsMnSmilesList)
        reaAddHsMnSmilesList = rxnAddHsMnSmilesList[0]
        proAddHsMnSmilesList = rxnAddHsMnSmilesList[1]
        #
        # The pattern of atom-atom mapping number
        pattern_AtomMapNum = re.compile("\[\d*[a-zA-Z\*]+@*[-\+]*\d*:(\d+)\]")
        #
        # The statistics of the atoms of unused compounds in the reaction.
        # We will exclude the atoms from unused(or avoided) compounds.
        unUsedAtomMapNumSet = set()
        for idx, compId in enumerate(usaReaCompIdList):
            if compId == "######":
                addHsMnSmiles = reaAddHsMnSmilesList[idx]
                atomMapNumList = re.findall(pattern_AtomMapNum, addHsMnSmiles)
                unUsedAtomMapNumSet |= set(atomMapNumList)          # Union Set
        for idx, compId in enumerate(usaProCompIdList):
            if compId == "######":
                addHsMnSmiles = proAddHsMnSmilesList[idx]
                atomMapNumList = re.findall(pattern_AtomMapNum, addHsMnSmiles)
                unUsedAtomMapNumSet |= set(atomMapNumList)          # Union Set
        unUsedAtomMapNumSet = {int(mapNum) for mapNum in unUsedAtomMapNumSet}
        #
        # The atom mapping number of one atom that belong to some avoided substructure of one compound has been set to 0.
        # We will exclude the atoms in generating compound pairs.
        # Example: CoA-SH; PPi; Pi;
        unUsedAtomMapNumSet.add(0)
        #
        # The pattern of atom symbol and atom-atom mapping number.
        pattern_AtSyMn = re.compile("\[\d*([a-zA-Z\*]+)@*[-\+]*\d*:(\d+)\]")
        #
        reaMnAtSyInfoList = list()
        for idx, compId in enumerate(usaReaCompIdList):
            # [mapping-number set, mapping-number atom-symbol dict, atom-number int, compound-id str]
            mapNumAtomSymbolInfo = [set(), dict(), 0, ""]
            if compId == "######":
                pass
            else:
                addHsMnSmiles = reaAddHsMnSmilesList[idx]
                atSyMnTupleList = re.findall(pattern_AtSyMn, addHsMnSmiles)
                for atSyMnTuple in atSyMnTupleList:
                    atomSymbol = atSyMnTuple[0].capitalize()
                    mapNum = int(atSyMnTuple[1])
                    #
                    if mapNum in unUsedAtomMapNumSet:
                        continue
                    #
                    mapNumAtomSymbolInfo[0].add(mapNum)
                    mapNumAtomSymbolInfo[1][mapNum] = atomSymbol
                mapNumAtomSymbolInfo[2] = len(mapNumAtomSymbolInfo[0])
            mapNumAtomSymbolInfo[3] = compId
            reaMnAtSyInfoList.append(mapNumAtomSymbolInfo)
        #
        proMnAtSyInfoList = list()
        for idx, compId in enumerate(usaProCompIdList):
            # [mapping-number set, mapping-number atom-symbol dict, atom-number int, compound-id str]
            mapNumAtomSymbolInfo = [set(), dict(), 0, ""]
            if compId == "######":
                pass
            else:
                addHsMnSmiles = proAddHsMnSmilesList[idx]
                atSyMnTupleList = re.findall(pattern_AtSyMn, addHsMnSmiles)
                for atSyMnTuple in atSyMnTupleList:
                    atomSymbol = atSyMnTuple[0].capitalize()
                    mapNum = int(atSyMnTuple[1])
                    #
                    if mapNum in unUsedAtomMapNumSet:
                        continue
                    #
                    mapNumAtomSymbolInfo[0].add(mapNum)
                    mapNumAtomSymbolInfo[1][mapNum] = atomSymbol
                mapNumAtomSymbolInfo[2] = len(mapNumAtomSymbolInfo[0])
            mapNumAtomSymbolInfo[3] = compId
            proMnAtSyInfoList.append(mapNumAtomSymbolInfo)
        #
        rxnMnAtSyInfoList = [reaMnAtSyInfoList, proMnAtSyInfoList]
        return rxnMnAtSyInfoList


    @classmethod
    def _getUniformSimilarity(cls, rxnMnAtSyInfoList, isDebug=False):
        """
            Get the uniform similarity of all reactant/product pairs(exclude cofactors),
          and we use Dice index(one specific format of Tversky index) to calculate the similarity
          between one reactant and one product.
            We will use the uniform similarity co construct compound pairs for pathway search.

            ARGUMENTS:
                ...

            RETURNS:
                {reactant(i)/product(j): uniform similarity, ..., }
        """
        # Analysis of input parameter.
        reaMnAtSyInfoList = rxnMnAtSyInfoList[0]
        proMnAtSyInfoList = rxnMnAtSyInfoList[1]
        #
        # Generate compound pairs
        rxnCompPairSet = set()
        for item_rea in reaMnAtSyInfoList:
            compId_rea = item_rea[-1]
            for item_pro in proMnAtSyInfoList:
                compId_pro = item_pro[-1]
                if (compId_rea != '######') and (compId_pro != '######'):
                    compPair = ' & '.join([compId_rea, compId_pro])
                    rxnCompPairSet.add(compPair)
        rxnCompPairList = sorted(rxnCompPairSet)
        #
        # Initializing
        rxnSinSimDict = dict()
        rxnAtSpeDict = dict()
        for compPair in rxnCompPairList:
            rxnSinSimDict[compPair] = list()
            rxnAtSpeDict[compPair] = set()
        #
        # calculate single similarity
        for item_rea in reaMnAtSyInfoList:
            mapNumSet_rea = item_rea[0]
            atSyDict_rea = item_rea[1]
            atNumInt_rea = item_rea[2]
            compId_rea = item_rea[3]
            for item_pro in proMnAtSyInfoList:
                mapNumSet_pro = item_pro[0]
                # atSyDict_pro = item_pro[1]
                atNumInt_pro = item_pro[2]
                compId_pro = item_pro[3]
                if (compId_rea != '######') and (compId_pro != '######'):
                    compPair = ' & '.join([compId_rea, compId_pro])
                    mapNumSet_com = mapNumSet_rea & mapNumSet_pro
                    atNumInt_com = len(mapNumSet_com)
                    # Dice similarity index
                    try:
                        tvSim = (2 * atNumInt_com)/(atNumInt_rea + atNumInt_pro)
                    except ZeroDivisionError:
                        tvSim = 0
                    rxnSinSimDict[compPair].append(tvSim)
                    for mapNum in mapNumSet_com:
                        atSy = atSyDict_rea[mapNum]
                        # atSy = atSyDict_pro[mapNum]
                        rxnAtSpeDict[compPair].add(atSy)
        #
        # calculate uniform similarity
        rxnUniSimDict = dict()
        for compPair, tvSimList in rxnSinSimDict.items():
            numPair = 0
            for tvSim in tvSimList:
                if tvSim != 0:
                    numPair += 1
            try:
                uniSim = sum(tvSimList)/numPair
            except ZeroDivisionError:
                uniSim = 0
            rxnUniSimDict[compPair] = uniSim
        #
        if isDebug:
            print(rxnUniSimDict)
            print(rxnAtSpeDict)
        return [rxnUniSimDict, rxnAtSpeDict]


    @classmethod
    def _getDrPairs(cls, rxnMnAtSyInfoList, diffThresholdLevel = 8,
                         trAtom = "C", lowThr=0.3,  isDebug=False):
        """
            Get all compound pairs of this reaction of every similarity difference threshold, e.g.,
          0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 and 1.0.
            Here we use conserved atom species, and lowest similarity threshold value to further
          restrict the compound pairs

            ARGUMENTS:
                ...

            RETURNS:
                [
                  [{reactant:[product, product, ..., ], ..., }, {reactant:[product, product, ..., ], ..., }],  # 0.1
                  [{product:[reactant, reactant, ..., ], ..., }, {product:[reactant, reactant, ..., ], ..., }],  # 0.2
                  ...,
                ]
        """
        #
        [rxnUniSimDict, rxnAtSpeDict] = cls._getUniformSimilarity(rxnMnAtSyInfoList)
        #
        diffThresholdScope = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        diffThresholdList = diffThresholdScope[:diffThresholdLevel]
        #
        # Unifying...
        reaCompSet = set()
        proCompSet = set()
        # The maximum value of uniform similarity
        reaMaxUniSimDict = dict()
        proMaxUniSimDict = dict()
        for compPair, uniSim in rxnUniSimDict.items():
            (compId_rea, compId_pro) = compPair.split(' & ')
            reaCompSet.add(compId_rea)
            proCompSet.add(compId_pro)
            try:
                if reaMaxUniSimDict[compId_rea] < uniSim:
                    reaMaxUniSimDict[compId_rea] = uniSim
            except KeyError:
                reaMaxUniSimDict[compId_rea] = uniSim
            try:
                if proMaxUniSimDict[compId_pro] < uniSim:
                    proMaxUniSimDict[compId_pro] = uniSim
            except KeyError:
                proMaxUniSimDict[compId_pro] = uniSim
        #
        reaCompList = sorted(reaCompSet)
        proCompList = sorted(proCompSet)
        #
        rxnDrCompPairList = list()
        for diffThreshold in diffThresholdList:
            #
            # reactant --> product
            reaCompPairDict = dict()
            for compId_rea in reaCompList:
                reaCompPairDict[compId_rea] = list()
                uniSim_max_rea = reaMaxUniSimDict[compId_rea]
                if uniSim_max_rea == 0:
                    continue
                for compId_pro in proCompList:
                    compPair = ' & '.join([compId_rea, compId_pro])
                    uniSim = rxnUniSimDict[compPair]
                    """
                    if (abs(uniSim_max_rea-uniSim)<diffThreshold) and (trAtom in rxnAtSpeDict[compPair]):
                        if uniSim == uniSim_max_rea:
                            reaCompPairDict[compId_rea].append(compId_pro)
                        elif uniSim >= lowThr:
                            reaCompPairDict[compId_rea].append(compId_pro)
                        else:
                            pass
                    """
                    if (abs(uniSim_max_rea-uniSim)<=diffThreshold) and (trAtom in rxnAtSpeDict[compPair]):
                        reaCompPairDict[compId_rea].append(compId_pro)
            #
            # product --> reactant
            proCompPairDict = dict()
            for compId_pro in proCompList:
                proCompPairDict[compId_pro] = list()
                uniSim_max_pro = proMaxUniSimDict[compId_pro]
                if uniSim_max_pro == 0:
                    continue
                for compId_rea in reaCompList:
                    compPair = ' & '.join([compId_rea, compId_pro])
                    uniSim = rxnUniSimDict[compPair]
                    """
                    if (abs(uniSim_max_pro-uniSim)<diffThreshold) and (trAtom in rxnAtSpeDict[compPair]):
                        if uniSim == uniSim_max_pro:
                            proCompPairDict[compId_pro].append(compId_rea)
                        elif uniSim >= lowThr:
                            proCompPairDict[compId_pro].append(compId_rea)
                        else:
                            pass
                    """
                    if (abs(uniSim_max_pro-uniSim)<diffThreshold) and (trAtom in rxnAtSpeDict[compPair]):
                        proCompPairDict[compId_pro].append(compId_rea)
            #
            rxnCompPairDict = [reaCompPairDict, proCompPairDict]
            rxnDrCompPairList.append(rxnCompPairDict)
        #
        #
        rxnDrCompPairList_temp = list()
        for rxnCompPairDict in rxnDrCompPairList:
            reaCompPairDict = rxnCompPairDict[0]
            proCompPairDict = rxnCompPairDict[1]
            reaCompPairDict_temp = dict()
            proCompPairDict_temp = dict()
            for compId_rea, compId_pro_list in reaCompPairDict.items():
                reaCompPairDict_temp[compId_rea] = list()
                for compId_pro in compId_pro_list:
                    try:
                        if compId_rea in proCompPairDict[compId_pro]:
                            reaCompPairDict_temp[compId_rea].append(compId_pro)
                    except KeyError:
                        pass
            for compId_pro, compId_rea_list in proCompPairDict.items():
                proCompPairDict_temp[compId_pro] = list()
                for compId_rea in compId_rea_list:
                    try:
                        if compId_pro in reaCompPairDict[compId_rea]:
                            proCompPairDict_temp[compId_pro].append(compId_rea)
                    except KeyError:
                        pass
            rxnCompPairDict_temp = [reaCompPairDict_temp, proCompPairDict_temp]
            rxnDrCompPairList_temp.append(rxnCompPairDict_temp)
        #
        rxnDrCompPairList = rxnDrCompPairList_temp
        #
        #
        if isDebug:
            for diffThreshold, rxnCompPairDict in zip(diffThresholdList, rxnDrCompPairList):
                print(diffThreshold, ':', rxnCompPairDict)
            print('')
        return rxnDrCompPairList


    @classmethod
    def getRxnCompDrPairs(cls, usaRxnCompIdList, exSubStrTupleList,
                               rxnCompMolBlockList, rxnAddHsMnSmilesList,
                               diffThresholdLevel = 8, atomSpecies = {'C', 'N'}, isDebug=True):
        """
            Get the compound pairs of this reaction for pathway search.

            ARGUMENTS:
                usaRxnCompIdList: [usaReaCompIdList, usaProCompIdList]
                exSubStrTupleList: [(subMol, rdkitSubBit, morganSubBit), ..., ]
                rxnCompMolBlockList: [reaCompMolBlockList, proCompMolBlockList]
                rxnAddHsMnSmilesList: [reaAddHsMnSmilesList, proAddHsMnSmilesList]

            RETURNS:
                [
                  [{reactant:[product, product, ..., ], ..., }, {reactant:[product, product, ..., ], ..., }],  # 0.1
                  [{product:[reactant, reactant, ..., ], ..., }, {product:[reactant, reactant, ..., ], ..., }],  # 0.2
                  ...,
                ]
        """
        # Analysis of input parameter.
        # ...
        #
        rxnMnAtSyInfoList = cls._getRxnMnAtSyInfo(usaRxnCompIdList, exSubStrTupleList,
                                                  rxnCompMolBlockList, rxnAddHsMnSmilesList)
        #
        rxnDrCompPairList = cls._getDrPairs(rxnMnAtSyInfoList, diffThresholdLevel)
        #
        return rxnDrCompPairList




    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__name__":
    """ """
    # precursor
    # successor
    pass





