# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:06:04 2021

@author: CC-SXF
"""

import json
import copy


class Ltiod():
    """
        Local Total Out-Degree.
    """

    def __init__(self, ):
        """ """
        pass


    @classmethod
    def _getBf(cls, compRxnDist, rxnCompPairDict,
                    diffThresholdLevel = 10,
                    ):
        """
            Get the branching factor of every compound
        """
        # Initializing...
        bf_list = list()
        for idx in range(diffThresholdLevel):
            bf_list.append(0)
        #
        bf_dict = dict()
        for compId in compRxnDist.keys():
            bf_dict[compId] = copy.deepcopy(bf_list)
        #
        for compId, compRxnList  in compRxnDist.items():
            for idx in range(diffThresholdLevel):
                branch_factor = 0
                for rxnId in compRxnList:
                    rxnDir = rxnCompPairDict[rxnId][0]
                    compPairDl = rxnCompPairDict[rxnId][1][idx]
                    if (rxnDir == "=>") or (rxnDir == "<=>"):
                        try:
                            genCompList = compPairDl[0][compId]
                            branch_factor += len(genCompList)
                            continue
                        except KeyError:
                            pass
                        # None: without reactant-product pairs
                        except TypeError:
                            pass
                    if (rxnDir == "<=") or (rxnDir == "<=>"):
                        try:
                            genCompList = compPairDl[1][compId]
                            branch_factor += len(genCompList)
                            continue
                        except KeyError:
                            pass
                        # None: without reactant-product pairs
                        except TypeError:
                            pass
                bf_dict[compId][idx] = branch_factor
        #
        return bf_dict


    @classmethod
    def _getSuccessor(cls, compRxnDist, rxnCompPairDict,
                           diffThresholdLevel = 10,  step = 2,
                           ):
        """
            Get the direct and indirect products of every compound.
        """
        # Initializing...
        successor_list = list()
        for idx in range(diffThresholdLevel):
            successor_list.append(set())
        #
        successor_dict = dict()
        for compId in compRxnDist.keys():
            successor_dict[compId] = copy.deepcopy(successor_list)
        #
        for source in compRxnDist.keys():
            for idx in range(diffThresholdLevel):
                #
                pathAllVisitCompSet = set()
                pathAllVisitComp2RxnSet = set()
                lastVisitRxn2CompSet = set()
                lastVisitRxn2CompSet.add('-->'.join(['Start', source]))
                #
                for st in range(step):
                    #
                    if lastVisitRxn2CompSet == set():
                        break
                    #
                    # Start the next round of searches.
                    newComp2RxnSet = set()
                    for rxn2comp in lastVisitRxn2CompSet:
                        [rxnId_ex, compId] = rxn2comp.split('-->')
                        compRxnList = compRxnDist[compId]
                        # Exclude minimal loops.
                        comp2RxnList = ["-->".join([compId, rxnId]) for rxnId in compRxnList if rxnId != rxnId_ex]
                        newComp2RxnSet |= set(comp2RxnList)
                    newComp2RxnSet -= pathAllVisitComp2RxnSet
                    pathAllVisitComp2RxnSet |= newComp2RxnSet
                    #
                    lastVisitRxn2CompSet = set()
                    #
                    # Continue the pathway search through the newly reached reations.
                    newVisitCompSet = set()
                    for compIdrxnId in newComp2RxnSet:
                        [compId, rxnId] = compIdrxnId.split("-->")
                        rxnDir = rxnCompPairDict[rxnId][0]
                        compPairDl = rxnCompPairDict[rxnId][1][idx]
                        # The information of the reaction's compound pairs doesn't exist.
                        if compPairDl is None:
                            continue
                        if (rxnDir == "=>") or (rxnDir == "<=>"):
                            try:
                                genCompList = compPairDl[0][compId]
                                newVisitCompSet |= set(genCompList)
                                lastVisitRxn2CompList = ['-->'.join([rxnId, genCompId]) for genCompId in genCompList]
                                lastVisitRxn2CompSet |= set(lastVisitRxn2CompList)
                                continue
                            except KeyError:
                                pass
                        if (rxnDir == "<=") or (rxnDir == "<=>"):
                            try:
                                genCompList = compPairDl[1][compId]
                                newVisitCompSet |= set(genCompList)
                                lastVisitRxn2CompList = ['-->'.join([rxnId, genCompId]) for genCompId in genCompList]
                                lastVisitRxn2CompSet |= set(lastVisitRxn2CompList)
                                continue
                            except KeyError:
                                pass
                    pathAllVisitCompSet |= newVisitCompSet
                pathAllVisitCompSet -= {source}
                successor_dict[source][idx] = pathAllVisitCompSet
        #
        return successor_dict


    @classmethod
    def _getLtod(cls, compRxnDist, rxnCompPairDict,
                      diffThresholdLevel = 10,  step = 2,
                      ):
        """
            Get the local total out-degree of every compound.

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        # Initializing...
        ltod_list = list()
        for idx in range(diffThresholdLevel):
            ltod_list.append(0)
        #
        ltod_dict = dict()
        for compId in compRxnDist.keys():
            ltod_dict[compId] = copy.deepcopy(ltod_list)
        #
        for source in compRxnDist.keys():
            for idx in range(diffThresholdLevel):
                #
                pathAllVisitComp2RxnSet = set()  # compound --> reaction
                lastVisitRxn2CompSet = set()
                lastVisitRxn2CompSet.add('-->'.join(['Start', source]))
                #
                for st in range(step):
                    #
                    if lastVisitRxn2CompSet == set():
                        break
                    #
                    # Start the next round of searches.
                    newComp2RxnSet = set()
                    for rxn2comp in lastVisitRxn2CompSet:
                        [rxnId_ex, compId] = rxn2comp.split('-->')
                        compRxnList = compRxnDist[compId]
                        # Exclude minimal loops.
                        comp2RxnList = ["-->".join([compId, rxnId]) for rxnId in compRxnList if rxnId != rxnId_ex]
                        newComp2RxnSet |= set(comp2RxnList)
                    newComp2RxnSet -= pathAllVisitComp2RxnSet
                    pathAllVisitComp2RxnSet |= newComp2RxnSet
                    #
                    lastVisitRxn2CompSet = set()
                    #
                    # Continue the pathway search through the newly reached reations.
                    newComp2Rxn2CompSet = set()  # compound --> reaction --> compound
                    for compIdrxnId in newComp2RxnSet:
                        [compId, rxnId] = compIdrxnId.split("-->")
                        rxnDir = rxnCompPairDict[rxnId][0]
                        compPairDl = rxnCompPairDict[rxnId][1][idx]
                        # The information of the reaction's compound pairs doesn't exist.
                        if compPairDl is None:
                            continue
                        if (rxnDir == "=>") or (rxnDir == "<=>"):
                            try:
                                genCompList = compPairDl[0][compId]
                                #
                                comp2Rxn2CompList = ["-->".join([compIdrxnId, genCompId]) for genCompId in genCompList]
                                comp2Rxn2CompSet = set(comp2Rxn2CompList)
                                newComp2Rxn2CompSet |= comp2Rxn2CompSet
                                #
                                lastVisitRxn2CompList = ['-->'.join([rxnId, genCompId]) for genCompId in genCompList]
                                lastVisitRxn2CompSet |= set(lastVisitRxn2CompList)
                                continue
                            except KeyError:
                                pass
                        if (rxnDir == "<=") or (rxnDir == "<=>"):
                            try:
                                genCompList = compPairDl[1][compId]
                                #
                                comp2Rxn2CompList = ["-->".join([compIdrxnId, genCompId]) for genCompId in genCompList]
                                comp2Rxn2CompSet = set(comp2Rxn2CompList)
                                newComp2Rxn2CompSet |= comp2Rxn2CompSet
                                #
                                lastVisitRxn2CompList = ['-->'.join([rxnId, genCompId]) for genCompId in genCompList]
                                lastVisitRxn2CompSet |= set(lastVisitRxn2CompList)
                                continue
                            except KeyError:
                                pass
                ltod_dict[source][idx] = len(newComp2Rxn2CompSet)
        #
        return ltod_dict


    @classmethod
    def _f2b(cls, forward_comp_pair):
        """
            Prepare the compound pairs for pathway retro-search.
        """
        rpf_comp_pair_dict = forward_comp_pair[0]
        prf_comp_pair_dict = forward_comp_pair[1]
        #
        rpb_comp_pair_dict = dict()
        for comp_id_list in rpf_comp_pair_dict.values():
            for comp_id in comp_id_list:
                rpb_comp_pair_dict[comp_id] = list()
        prb_comp_pair_dict = dict()
        for comp_id_list in prf_comp_pair_dict.values():
            for comp_id in comp_id_list:
                prb_comp_pair_dict[comp_id] = list()
        #
        for key, comp_id_list in rpf_comp_pair_dict.items():
            for comp_id in comp_id_list:
                rpb_comp_pair_dict[comp_id].append(key)
        for key, comp_id_list in prf_comp_pair_dict.items():
            for comp_id in comp_id_list:
                prb_comp_pair_dict[comp_id].append(key)
        #
        backward_comp_pair = [rpb_comp_pair_dict, prb_comp_pair_dict]
        return backward_comp_pair


    @classmethod
    def _getLtid(cls, compRxnDist, rxnCompPairDict,
                      diffThresholdLevel = 10,  step = 2,
                      ):
        """
            Get the local total in-degree of every compound.

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        rxnCompPairDict_temp = dict()
        for rxnId, value in rxnCompPairDict.items():
            rxnDir, compPairList = value
            compPairList_temp = list()
            for forward_comp_pair in compPairList:
                if forward_comp_pair == None:
                    backward_comp_pair = None
                else:
                    backward_comp_pair = cls._f2b(forward_comp_pair)
                compPairList_temp.append(backward_comp_pair)
            rxnCompPairDict_temp[rxnId] = [rxnDir, compPairList_temp]
        #
        ltid_dict = cls._getLtod(compRxnDist, rxnCompPairDict_temp, diffThresholdLevel,  step)
        #
        return ltid_dict


    @classmethod
    def _gfLtiod(cls, compDistInfoList, rxnFuseInfoList,
                     diffThresholdLevel = 10, step = 2,
                     ):
        """
            Get and fuse the local total in/out-degree of every compound.

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        # Initializing...
        ltiod_list = list()
        for idx in range(diffThresholdLevel):
            # # ltiod_list.append(0)
            ltiod_list.append([0, 0])
        #
        # Update title
        # # compDistInfoList[0].insert(6, "Branching Factor(0.1 - 1.0)")
        compDistInfoList[0].insert(6, "Local Total In/Out-Degree(0.1 - 1.0)")
        #
        # Insert local total out-degree
        compRxnDist = dict()
        for compDistInfo in compDistInfoList[1:]:  # Exclude title
            compId = compDistInfo[0]
            compRxnList = json.loads(compDistInfo[6])           # Reactions
            compRxnDist[compId] = compRxnList
            compDistInfo.insert(6, copy.deepcopy(ltiod_list))
        #
        rxnCompPairDict = dict()
        for rxnFuseInfo in rxnFuseInfoList[1:]:  # Exclude title
            rxnId = rxnFuseInfo[0]
            rxnDir = rxnFuseInfo[4].strip('"')                  # Direaction
            compPairList = rxnFuseInfo[8:8+diffThresholdLevel]      # Compound Pair (0.1 - 1.0)
            compPairList = [json.loads(item) for item in compPairList]
            rxnCompPairDict[rxnId] = [rxnDir, compPairList]
        #
        # Get the branching factor of every compound.
        # # bf_dict = cls._getBf(compRxnDist, rxnCompPairDict, diffThresholdLevel)
        #
        # Get the local total in-degree of every compound.
        ltid_dict = cls._getLtid(compRxnDist, rxnCompPairDict, diffThresholdLevel,  step)
        # Get the local total out-degree of every compound.
        ltod_dict = cls._getLtod(compRxnDist, rxnCompPairDict, diffThresholdLevel,  step)
        #
        # Update the local total in/out-degree of every compound.
        for compDistInfo in compDistInfoList[1:]:
            source = compDistInfo[0]
            for idx in range(diffThresholdLevel):
                # # source_bf = bf_dict[source][idx]                    # branching factor
                # # compDistInfo[6][idx] = source_bf
                source_ltid = ltid_dict[source][idx]                # local total in degree
                source_ltod = ltod_dict[source][idx]                # local total out degree
                compDistInfo[6][idx] = [source_ltid, source_ltod]
        #
        for compDistInfo in compDistInfoList[1:]:
            compDistInfo[6] = json.dumps(compDistInfo[6])
        #
        return compDistInfoList


    @classmethod
    def demoFunc(cls, ):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





