# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 10:13:28 2020

@author: CC-SXF
"""

# import re
# import sys
# import copy
import psutil
from datetime import datetime
from PyQt5 import QtCore


class PwNaDfs(QtCore.QObject):
    """
        Pathway search tools using Depth-First Search algorithm.
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    dfs_prompts_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        """ """
        super().__init__()


    def dfs_throw_prompts(self, prompts):
        """ """
        self.dfs_prompts_signal.emit(prompts)


    def _lpSort(self, totalPathwayList, maxlength = 20):
        """
            Sort the pathways by their length and string.

            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        lengthpathway_dict = dict()
        for length in range(1, maxlength+1):
            lengthpathway_dict[length] = set()
        #
        for lengthpathway in totalPathwayList:
            length = lengthpathway[0]
            pathway = lengthpathway[-1]
            lengthpathway_dict[length].add(pathway)
        #
        totalPathwayList = list()
        for length in range(1, maxlength+1):
            pathway_list = sorted(lengthpathway_dict[length])
            lengthpathway_list = [(length, pathway) for pathway in pathway_list]
            totalPathwayList += lengthpathway_list
        #
        return totalPathwayList


    def _unifyPathway(self, pathwayList):
        """
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        # '0-Feasibility';
        # '1-TotalLength'; '2-EndoLength'; '3-HeterLength'; '4-InfLength';
        # '5-AtomUtilization'; '6-AtomConservation'; '7-Flux';
        # '8-Pathway';
        # '9-I-Details'; '10-C-Details';
        #
        unifyPathwayList = list()
        for pathway_info in pathwayList:
            length = pathway_info[0]
            pathway = pathway_info[1]
            details_1 = [1 for idx in range(length)]
            details_2 = [list() for idx in range(length)]
            pathway_info = [True,
                            length, length, 0, 0,
                            None, None, None,
                            pathway,
                            details_1, details_2,
                            ]
            unifyPathwayList.append(pathway_info)
        #
        return unifyPathwayList


    def _pwSource2Target(self, source, target, pathwaylength = 0,
                               currentCompRxnList = list(), visitCompSetList = list(), totalPathwayList = list(),
                               maxlength = 20, maxtime = 100, maxnumber = 1000000,
                               isTotal = True,
                               rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                               timer = {'start':0, "consume":0, "counter":0,  "threshold":10000},
                               memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16},
                               exDir = 'forward',
                               ):
        """
          Depth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        # The maximal value of the pathway length has been reached.
        if pathwaylength >= maxlength:
            return None
        #
        # The maximal searching time has been reached.
        if timer['consume'] >= maxtime:
                return None
        #
        # Deal with out of memory error.
        # # memory = psutil.virtual_memory()
        # # used_memory_percent = memory.percent
        # # used_memory = round(memory.used/1024/1024/1024, 2)
        if (memorizer["used_memory_percent"] > memorizer["mempercent"]) or(memorizer["used_memory"] > memorizer["memgb"]):
            return None
        #
        # The maximal number of the retrieved pathways has been reached.
        if len(totalPathwayList) > maxnumber:
            return None
        #
        # Update the counter and the running time.
        # Update the memory
        timer['counter'] += 1
        if (timer['counter'] % timer['threshold']) == 0:
            interTime = datetime.now()
            timer['consume'] = (interTime - timer['start']).total_seconds()
            #
            memory = psutil.virtual_memory()
            used_memory_percent = memory.percent
            used_memory = round(memory.used/1024/1024/1024, 2)
            memorizer["used_memory_percent"] = used_memory_percent
            memorizer["used_memory"] = used_memory
        #
        comprxn_list = currentCompRxnList[:(2*pathwaylength+1)]
        comprxn_set = set(comprxn_list)
        lastVisitComp = comprxn_list[-1]
        #
        # Explicit/implicit visited compounds will be excluded.
        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
        #   0       1      2       3      4       5      6       7      8
        visitCompSet = visitCompSetList[pathwaylength]
        """
        visitCompSet = set()
        for idx, comprxn in enumerate(comprxn_list):
            if idx%2 == 0:
                continue
            rxnId = comprxn
            # # rxnCompList = rxnDistInfoDict[rxnId][1]
            rxnCompSet = rxnDistInfoDict[rxnId][2]
            visitCompSet |= rxnCompSet
        """
        #
        # Continue the pathway search through the last visited compound
        compDistInfo = compDistInfoDict[lastVisitComp]
        # compName = compDistInfo[0]
        # compBf = compDistInfo[1]
        compRxnList = compDistInfo[2]
        compRxnSet = set(compRxnList)
        # Exclude the reactions have been visited by the current pathway.
        compRxnSet -= comprxn_set  # Difference set
        # There are no new reactions to extend the current pathway.
        if compRxnSet == set():
            return None
        #
        # Continue the pathway search through newly reached reactions.
        for rxnId in compRxnSet:
            rxnDistInfo = rxnDistInfoDict[rxnId]
            rxnDir = rxnDistInfo[0]
            # # rxnCompList = rxnDistInfo[1]
            rxnCompSet = rxnDistInfo[2]
            rxnCompPairList = rxnDistInfo[3][exDir]
            #
            # The information of this reaction's compound pairs doesn't exist.
            if rxnCompPairList is None:
                continue
            #
            if (rxnDir == '=>') or (rxnDir == '<=>'):
                try:
                    genCompIdList = rxnCompPairList[0][lastVisitComp]
                    genCompIdSet = set(genCompIdList)
                    #
                    # Every compound in this pathway must be different from the compounds
                    # that belong to its direct/indirect predecessor reactions.
                    # Exclude explicit/implicit visited compounds.
                    genCompIdSet -= visitCompSet  # Difference set
                    #
                    # There are no new compounds to extend the current pathway.
                    if genCompIdSet == set():
                        continue
                    #
                    for genCompId in genCompIdSet:
                        # Extend the current pathway.
                        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
                        #   0       1      2       3      4       5      6       7      8
                        currentCompRxnList[2*pathwaylength+1] = rxnId
                        currentCompRxnList[2*pathwaylength+2] = genCompId
                        newPl = pathwaylength + 1
                        #
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if genCompId == target:
                            newPathway = "-->".join(currentCompRxnList[:2*newPl+1])
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((newPl, newPathway))
                            elif newPl == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((newPl, newPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                            continue
                        #
                        # The maximal value of the pathway length has been reached.
                        if newPl >= maxlength:
                            continue
                        #
                        # Continue the pathway search process by using recursive algorithm.
                        visitCompSetList[newPl] = (visitCompSetList[newPl-1] | rxnCompSet)
                        self._pwSource2Target(source, target, newPl,
                                              currentCompRxnList, visitCompSetList, totalPathwayList,
                                              maxlength, maxtime, maxnumber,
                                              isTotal,
                                              rxnDistInfoDict, compDistInfoDict,
                                              timer,
                                              memorizer,
                                              )
                    continue
                except KeyError:
                    pass
            #
            if (rxnDir == "<=") or (rxnDir == '<=>'):
                try:
                    genCompIdList = rxnCompPairList[1][lastVisitComp]
                    genCompIdSet = set(genCompIdList)
                    #
                    # Every compound in this pathway must be different from the compounds
                    # that belong to its direct/indirect predecessor reactions.
                    # Exclude explicit/implicit visited compounds.
                    genCompIdSet -= visitCompSet  # Difference set
                    #
                    # There are no new compounds to extend the current pathway.
                    if genCompIdSet == set():
                        continue
                    #
                    for genCompId in genCompIdSet:
                        # Extend the current pathway.
                        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
                        #   0       1      2       3      4       5      6       7      8
                        currentCompRxnList[2*pathwaylength+1] = rxnId
                        currentCompRxnList[2*pathwaylength+2] = genCompId
                        newPl = pathwaylength + 1
                        #
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if genCompId == target:
                            newPathway = "-->".join(currentCompRxnList[:2*newPl+1])
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((newPl, newPathway))
                            elif newPl == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((newPl, newPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                            continue
                        #
                        # The maximal value of the pathway length has been reached.
                        if newPl >= maxlength:
                            continue
                        #
                        # Continue the pathway search process by using recursive algorithm.
                        visitCompSetList[newPl] = (visitCompSetList[newPl-1] | rxnCompSet)
                        self._pwSource2Target(source, target, newPl,
                                              currentCompRxnList, visitCompSetList, totalPathwayList,
                                              maxlength, maxtime, maxnumber,
                                              isTotal,
                                              rxnDistInfoDict, compDistInfoDict,
                                              timer,
                                              memorizer,
                                              )
                    continue
                except KeyError:
                    pass
        #
        return None


    def Source2Target(self, source, target,
                            maxlength = 20, maxtime = 100, maxnumber = 1000000,
                            isTotal = True,
                            rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                            timer = {"start":0, "consume":0, "counter":0,  "threshold":10000},
                            memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16},
                            isDebug = False,
                            ):
        """
            Search the pathways from source to target by using
          Depth-First Search algorithm.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        if isDebug:
            # # print("")
            print("Search the pathways from source to target (DFS):")
        else:
            # # self.dfs_throw_prompts("")
            self.dfs_throw_prompts("Search the pathways from source to target (DFS):")
        #
        if (not isinstance(source, str)) or (not isinstance(target, str)):
            if isDebug:
                print("Incorrect parameter type.")
            else:
                self.dfs_throw_prompts("Incorrect parameter type.")
            return list()
        #
        allCompSet = set(compDistInfoDict.keys())
        if source not in allCompSet:
            if isDebug:
                print(f"The source({source}) does not exist in the Generalized Metabolic Space.")
                print("Please have a check on the source.")
            else:
                self.dfs_throw_prompts(f"The source({source}) does not exist in the Generalized Metabolic Space.")
                self.dfs_throw_prompts("Please have a check on the source.")
            return list()
        if target not in allCompSet:
            if isDebug:
                print(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                print("Try to use de novo pathway design instead of pathway search.")
            else:
                self.dfs_throw_prompts(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                self.dfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
            return list()
        #
        sourceName = compDistInfoDict[source][0]
        targetName = compDistInfoDict[target][0]
        if source == target:
            if isDebug:
                print(f"The source({source}, {sourceName}) is identical with the target({target}, {targetName}).")
                print("Please have a check on the source and the target.")
            else:
                self.dfs_throw_prompts(f"The source({source}, {sourceName}) is identical with the target({target}, {targetName}).")
                self.dfs_throw_prompts("Please have a check on the source and the target.")
            return list()
        #
        # Initialization
        pathwaylength = 0
        currentCompRxnList = list()
        currentCompRxnList.append(source)
        for idx in range(1, 2*maxlength+1):
            currentCompRxnList.append("")
        visitCompSetList = list()
        for idx in range(maxlength+1):
            visitCompSetList.append(set())
        visitCompSetList[0].add(source)
        totalPathwayList=list()
        startTime = datetime.now()
        # Resetting...
        timer = {'start':startTime, "consume":0, "counter":0,  "threshold":10000}
        memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16}
        #
        self._pwSource2Target(source, target, pathwaylength,
                              currentCompRxnList, visitCompSetList, totalPathwayList,
                              maxlength, maxtime, maxnumber,
                              isTotal,
                              rxnDistInfoDict, compDistInfoDict,
                              timer,
                              memorizer,
                              )
        #
        # Calculate the running time.
        endTime = datetime.now()
        runTime = (endTime - startTime).total_seconds()
        if isDebug:
            print(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        else:
            self.dfs_throw_prompts(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        #
        # Deal with out of memory error.
        used_memory_percent = memorizer["used_memory_percent"]
        used_memory = memorizer["used_memory"]
        mempercent = memorizer["mempercent"]
        memgb = memorizer["memgb"]
        #
        if isDebug:
            print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        else:
            self.dfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        #
        # Sorting...
        # # totalPathwayList = self._lpSort(totalPathwayList)
        #
        # Unifying...
        totalPathwayList = self._unifyPathway(totalPathwayList)
        pathwayNum = len(totalPathwayList)
        #
        if pathwayNum == 0:
            if isDebug:
                print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
            else:
                self.dfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
            if runTime >= maxtime:
                if isDebug:
                    print("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                else:
                    self.dfs_throw_prompts("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
            elif (used_memory_percent > mempercent) or (used_memory > memgb):
                if isDebug:
                    # # print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    print("The memory of this platform restricts the pathway search process.")
                    pass
                else:
                    # # self.dfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    self.dfs_throw_prompts("The memory of this platform restricts the pathway search process.")
                    pass
            else:
                if isTotal:
                    if isDebug:
                        print("Try to increase the maximal length of the potential pathway and repeat the search process.")
                    else:
                        self.dfs_throw_prompts("Try to increase the maximal length of the potential pathway and repeat the search process.")
                else:
                    if isDebug:
                        print("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                    else:
                        self.dfs_throw_prompts("Try to decrease or increase the length of the potential pathway and repeat the search process.")
        else:
            if isDebug:
                print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
                print(f"The total number of the retrieved pathways is {pathwayNum}.\n")
            else:
                self.dfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
                self.dfs_throw_prompts(f"The total number of the retrieved pathways is {pathwayNum}.\n")
        #
        return totalPathwayList




    def _pwTarget2MultSource(self, multSourceSet, target, pathwaylength = 0,
                                   currentCompRxnList = list(), totalPathwayList = list(),
                                   maxlength = 20, maxtime = 100, maxnumber = 1000000,
                                   isTotal = True,
                                   rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                                   timer = {'start':0, "consume":0, "counter":0,  "threshold":10000},
                                   memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16},
                                   exDir='backward',
                                   ):
        """
            Depth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        # The maximal value of the pathway length has been reached.
        if pathwaylength >= maxlength:
            return None
        #
        # The maximal searching time has been reached.
        if timer['consume'] >= maxtime:
                return None
        #
        # Deal with out of memory error.
        # # memory = psutil.virtual_memory()
        # # used_memory_percent = memory.percent
        # # used_memory = round(memory.used/1024/1024/1024, 2)
        if (memorizer["used_memory_percent"] > memorizer["mempercent"]) or(memorizer["used_memory"] > memorizer["memgb"]):
            return None
        #
        # The maximal number of the retrieved pathways has been reached.
        if len(totalPathwayList) > maxnumber:
            return None
        #
        # Update the counter and the running time.
        # Update the memory
        timer['counter'] += 1
        if (timer['counter'] % timer['threshold']) == 0:
            interTime = datetime.now()
            timer['consume'] = (interTime - timer['start']).total_seconds()
            #
            memory = psutil.virtual_memory()
            used_memory_percent = memory.percent
            used_memory = round(memory.used/1024/1024/1024, 2)
            memorizer["used_memory_percent"] = used_memory_percent
            memorizer["used_memory"] = used_memory
        #
        #
        comprxn_list = currentCompRxnList[:(2*pathwaylength+1)]
        comprxn_set = set(comprxn_list[:-1])
        lastVisitComp = comprxn_list[-1]
        #
        # Continue the pathway search through the last visited compound
        compDistInfo = compDistInfoDict[lastVisitComp]
        # compName = compDistInfo[0]
        # compBf = compDistInfo[1]
        compRxnList = compDistInfo[2]
        compRxnSet = set(compRxnList)
        # Exclude the reactions have been visited by the current pathway.
        compRxnSet -= comprxn_set  # Difference set
        #
        # There are no new reactions to extend the current pathway.
        if compRxnSet == set():
            return None
        #
        # Continue the pathway search through newly reached reactions.
        for rxnId in compRxnSet:
            rxnDistInfo = rxnDistInfoDict[rxnId]
            rxnDir = rxnDistInfo[0]
            # rxnCompList = rxnDistInfo[1]
            rxnCompSet = rxnDistInfo[2]
            rxnCompPairList = rxnDistInfo[3][exDir]
            #
            # The information of this reaction's compound pairs doesn't exist.
            if rxnCompPairList is None:
                continue
            #
            # Explicit/implicit visited compounds will be excluded.
            # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
            #   0       1      2       3      4       5      6       7      8
            if len(rxnCompSet & comprxn_set) != 0:
                continue
            #
            if (rxnDir == '=>') or (rxnDir == '<=>'):
                try:
                    # left-->right
                    genCompIdList = rxnCompPairList[0][lastVisitComp]
                    #
                    # There are no new compounds to extend the current pathway.
                    if len(genCompIdList) == 0:
                        continue
                    #
                    for genCompId in genCompIdList:
                        # Extend the current pathway.
                        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
                        #   0       1      2       3      4       5      6       7      8
                        currentCompRxnList[2*pathwaylength+1] = rxnId
                        currentCompRxnList[2*pathwaylength+2] = genCompId
                        newPl = pathwaylength + 1
                        #
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if genCompId in multSourceSet:
                            newPathway = "-->".join(currentCompRxnList[:2*newPl+1][::-1])
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((newPl, newPathway))
                            elif newPl == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((newPl, newPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                            continue
                        #
                        # The maximal value of the pathway length has been reached.
                        if newPl >= maxlength:
                            continue
                        #
                        # Continue the pathway search process by using recursive algorithm.
                        self._pwTarget2MultSource(multSourceSet, target, newPl,
                                                  currentCompRxnList, totalPathwayList,
                                                  maxlength, maxtime, maxnumber,
                                                  isTotal,
                                                  rxnDistInfoDict, compDistInfoDict,
                                                  timer,
                                                  memorizer,
                                                  )
                    continue
                except KeyError:
                    pass
            #
            if (rxnDir == "<=") or (rxnDir == '<=>'):
                try:
                    # right-->left
                    genCompIdList = rxnCompPairList[1][lastVisitComp]
                    #
                    # There are no new compounds to extend the current pathway.
                    if len(genCompIdList) == 0:
                        continue
                    #
                    for genCompId in genCompIdList:
                        # Extend the current pathway.
                        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
                        #   0       1      2       3      4       5      6       7      8
                        currentCompRxnList[2*pathwaylength+1] = rxnId
                        currentCompRxnList[2*pathwaylength+2] = genCompId
                        newPl = pathwaylength + 1
                        #
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if genCompId in multSourceSet:
                            newPathway = "-->".join(currentCompRxnList[:2*newPl+1][::-1])
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((newPl, newPathway))
                            elif newPl == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((newPl, newPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                            continue
                        #
                        # The maximal value of the pathway length has been reached.
                        if newPl >= maxlength:
                            continue
                        #
                        # Continue the pathway search process by using recursive algorithm.
                        self._pwTarget2MultSource(multSourceSet, target, newPl,
                                                  currentCompRxnList, totalPathwayList,
                                                  maxlength, maxtime, maxnumber,
                                                  isTotal,
                                                  rxnDistInfoDict, compDistInfoDict,
                                                  timer,
                                                  memorizer,
                                                  )
                    continue
                except KeyError:
                    pass
        #
        return None


    def Target2MultSource(self, multSourceSet, target,
                                maxlength = 20, maxtime = 100, maxnumber = 1000000,
                                isTotal = True,
                                rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                                timer = {'start':0, "consume":0, "counter":0,  "threshold":10000},
                                memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16},
                                isDebug = False,
                                ):
        """
            Search the pathways from target to mult-sources by using
          Depth-First Search algorithm.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        if isDebug:
            # # print("")
            print("Search the pathways from target to mult-sources (DFS):")
        else:
            # # self.dfs_throw_prompts("")
            self.dfs_throw_prompts("Search the pathways from target to mult-sources (DFS):")
        #
        if (not isinstance(multSourceSet, set)) or (not isinstance(target, str)):
            # print("Incorrect parameter type.")
            self.dfs_throw_prompts("Incorrect parameter type.")
            return list()
        #
        allCompSet = set(compDistInfoDict.keys())
        if target not in allCompSet:
            if isDebug:
                print(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                print("Try to use de novo pathway design instead of pathway search.")
            else:
                self.dfs_throw_prompts(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                self.dfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
            return list()
        #
        sourceEliSet = {target}
        for source in multSourceSet:
            if source not in allCompSet:
                sourceEliSet.add(source)
                # print(f"The source({source}) does not exist in the Generalized Metabolic Space.")
        #
        multSourceSet -= sourceEliSet                  # Difference Set
        #
        if multSourceSet == set():
            if isDebug:
                print("The valid 'source set' in empty.")
                print("Please have a check on the 'source set'.")
            else:
                self.dfs_throw_prompts("The valid 'source set' in empty.")
                self.dfs_throw_prompts("Please have a check on the 'source set'.")
            return list()
        #
        multSourceIdNameSet = set()
        for source in list(multSourceSet)[:6]:
            sourceName = compDistInfoDict[source][0]
            sourceIdName = "".join([source, "(", sourceName, ")"])
            multSourceIdNameSet.add(sourceIdName)
        if len(multSourceSet) > 6:
            multSourceIdNameSet.add('...')
        targetName = compDistInfoDict[target][0]
        #
        # Initialization
        pathwaylength = 0
        currentCompRxnList = list()
        currentCompRxnList.append(target)
        for idx in range(1, 2*maxlength+1):
            currentCompRxnList.append("")
        totalPathwayList=list()
        startTime = datetime.now()
        # Resetting...
        timer = {'start':startTime, "consume":0, "counter":0,  "threshold":10000}
        memorizer = {"used_memory_percent":0, "used_memory":0, "mempercent":85, "memgb":16}
        #
        self._pwTarget2MultSource(multSourceSet, target, pathwaylength,
                                  currentCompRxnList, totalPathwayList,
                                  maxlength, maxtime, maxnumber,
                                  isTotal,
                                  rxnDistInfoDict, compDistInfoDict,
                                  timer,
                                  memorizer,
                                  )
        #
        # Calculate the running time.
        endTime = datetime.now()
        runTime = (endTime - startTime).total_seconds()
        if isDebug:
            print(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        else:
            self.dfs_throw_prompts(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        #
        # Deal with out of memory error.
        used_memory_percent = memorizer["used_memory_percent"]
        used_memory = memorizer["used_memory"]
        mempercent = memorizer["mempercent"]
        memgb = memorizer["memgb"]
        #
        if isDebug:
            print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        else:
            self.dfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        #
        # Sorting...
        # # totalPathwayList = self._lpSort(totalPathwayList)
        #
        # Unifying...
        totalPathwayList = self._unifyPathway(totalPathwayList)
        pathwayNum = len(totalPathwayList)
        #
        if pathwayNum == 0:
            if isDebug:
                print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
            else:
                self.dfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
            if runTime >= maxtime:
                if isDebug:
                    print("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                else:
                    self.dfs_throw_prompts("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
            elif (used_memory_percent > mempercent) or (used_memory > memgb):
                if isDebug:
                    # # print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    print("The memory of this platform restricts the pathway search process.")
                    pass
                else:
                    # # self.dfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    self.dfs_throw_prompts("The memory of this platform restricts the pathway search process.")
                    pass
            else:
                if isTotal:
                    if isDebug:
                        print("Try to increase the maximal length of the potential pathway and repeat the search process.")
                    else:
                        self.dfs_throw_prompts("Try to increase the maximal length of the potential pathway and repeat the search process.")
                else:
                    if isDebug:
                        print("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                    else:
                        self.dfs_throw_prompts("Try to decrease or increase the length of the potential pathway and repeat the search process.")
            return list()
        #
        multSourceSet = set()
        for pathway_info in totalPathwayList:
            # # pathway = pathway_info['8-Pathway']
            pathway = pathway_info[8]
            source = pathway.split('-->')[0]
            multSourceSet.add(source)
        multSourceIdNameSet = set()
        for source in multSourceSet:
            sourceName = compDistInfoDict[source][0]
            sourceIdName = "".join([source, "(", sourceName, ")"])
            multSourceIdNameSet.add(sourceIdName)
        targetName = compDistInfoDict[target][0]
        #
        if isDebug:
            print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
            print(f"The total number of the retrieved pathways is {pathwayNum}.\n")
        else:
            self.dfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
            self.dfs_throw_prompts(f"The total number of the retrieved pathways is {pathwayNum}.\n")
        #
        return totalPathwayList


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





