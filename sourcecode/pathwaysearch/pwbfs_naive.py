# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 09:12:12 2021

@author: CC-SXF
"""

import re
import sys
# import copy
import psutil
from datetime import datetime
from PyQt5 import QtCore


class PwNaBfs(QtCore.QObject):
    """
        Pathway search tools using Breadth-First Search algorithm.
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    bfs_prompts_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        """ """
        super().__init__()


    def bfs_throw_prompts(self, prompts):
        """ """
        self.bfs_prompts_signal.emit(prompts)


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


    def _extendForward(self, target, pathwaylength = 0,
                             maxlength = 20, isTotal = True,
                             currentPathway = '',
                             totalPathwayList = list(),
                             rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                             exDir = 'forward',
                             ):
        """
            Breadth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        extPathwaySet = set()
        #
        currentCompRxnList = currentPathway.split('-->')
        currentCompRxnSet = set(currentCompRxnList)
        lastVisitComp = currentCompRxnList[-1]
        #
        # Explicit/implicit visited compounds will be excluded.
        # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
        #   0       1      2       3      4       5      6       7      8
        visitCompSet = set()
        for idx, comprxn in enumerate(currentCompRxnList):
            if idx%2 == 0:
                continue
            rxnId = comprxn
            # # rxnCompList = rxnDistInfoDict[rxnId][1]
            rxnCompSet = rxnDistInfoDict[rxnId][2]
            visitCompSet |= rxnCompSet
        #
        # Continue the pathway search through the last visit compound of the current pathway.
        compDistInfo = compDistInfoDict[lastVisitComp]
        # compName = compDistInfo[0]
        # compBf = compDistInfo[1]
        compRxnList = compDistInfo[2]
        newReachRxnSet = set(compRxnList)
        # Exclude the reactions haved been visited by the current pathway.
        newReachRxnSet -= currentCompRxnSet  # Difference set
        # There are no new reactions to extend the current pathway.
        if newReachRxnSet == set():
            return extPathwaySet
        #
        # Continue the pathway search through the newly reached reations.
        for newReachRxn in newReachRxnSet:
            #
            rxnDistInfo = rxnDistInfoDict[newReachRxn]
            rxnDir = rxnDistInfo[0]
            # # rxnCompList = rxnDistInfo[1]
            # rxnCompSet = rxnDistInfo[2]
            rxnCompPairList = rxnDistInfo[3][exDir]
            #
            # The information of the reaction's compound pairs doesn't exist.
            if rxnCompPairList is None:
                continue
            #
            if (rxnDir == '=>') or (rxnDir == '<=>'):
                try:
                    # left-->right
                    newReachCompList = rxnCompPairList[0][lastVisitComp]
                    newReachCompSet = set(newReachCompList)
                    #
                    # Every compound in this pathway must be different from the compounds
                    # that belong to its direct/indirect predecessor reactions.
                    # Exclude explicit/implicit visited compounds.
                    newReachCompSet -= visitCompSet  # Difference set
                    #
                    # There are no new compounds to extend the current pathway.
                    if newReachCompSet == set():
                        continue
                    #
                    for newReachComp in newReachCompSet:
                        # Extend the current pathway.
                        extendPathway = '-->'.join([currentPathway, newReachRxn, newReachComp])
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if newReachComp == target:
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((pathwaylength, extendPathway))
                            elif pathwaylength == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((pathwaylength, extendPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                        else:
                            extPathwaySet.add(extendPathway)
                    continue
                except KeyError:
                    pass
            #
            if (rxnDir == '<=') or (rxnDir == '<=>'):
                try:
                    # right-->left
                    newReachCompList = rxnCompPairList[1][lastVisitComp]
                    newReachCompSet = set(newReachCompList)
                    #
                    # Every compound in this pathway must be different from the compounds
                    # that belong to its direct/indirect predecessor reactions.
                    # Exclude explicit/implicit visited compounds.
                    newReachCompSet -= visitCompSet  # Difference set
                    #
                    # There are no new compounds to extend the current pathway.
                    if newReachCompSet == set():
                        continue
                    #
                    for newReachComp in newReachCompSet:
                        # Extend the current pathway.
                        extendPathway = '-->'.join([currentPathway, newReachRxn, newReachComp])
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if newReachComp == target:
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((pathwaylength, extendPathway))
                            elif pathwaylength == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((pathwaylength, extendPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                        else:
                            extPathwaySet.add(extendPathway)
                    continue
                except KeyError:
                    pass
        #
        return extPathwaySet


    def _pwSource2Target(self, source, target, maxlength = 20, maxtime = 120, maxnumber = 1000000,
                               isTotal = True,
                               isShortest = False,
                               rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                               isDebug = False, isDisplay = False,
                               threshold = 10000,
                               mempercent = 85, memgb = 16,
                               ):
        """
            Breadth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        if (not isinstance(source, str)) or (not isinstance(target, str)):
            if isDebug:
                print("Incorrect parameter type.")
            else:
                self.bfs_throw_prompts("Incorrect parameter type.")
            return list()
        #
        allCompSet = set(compDistInfoDict.keys())
        if source not in allCompSet:
            if isDebug:
                print(f"The source({source}) does not exist in the Generalized Metabolic Space.")
                print("Please have a check on the source.")
            else:
                self.bfs_throw_prompts(f"The source({source}) does not exist in the Generalized Metabolic Space.")
                self.bfs_throw_prompts("Please have a check on the source.")
            return list()
        if target not in allCompSet:
            if isDebug:
                print(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                print("Try to use de novo pathway design instead of pathway search.")
            else:
                self.bfs_throw_prompts(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                self.bfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
            return list()
        #
        sourceName = compDistInfoDict[source][0]
        targetName = compDistInfoDict[target][0]
        if source == target:
            if isDebug:
                print(f"The source({source}, {sourceName}) is identical with the target({target}, {targetName}).")
                print("Please have a check on the source and the target.")
            else:
                self.bfs_throw_prompts(f"The source({source}, {sourceName}) is identical with the target({target}, {targetName}).")
                self.bfs_throw_prompts("Please have a check on the source and the target.")
            return list()
        if isDebug:
            print("Intermediate pathways' length/number/size:")
        #
        # Initialization
        pathwaylength = 0
        memory = psutil.virtual_memory()
        used_memory_percent = memory.percent
        used_memory = round(memory.used/1024/1024/1024, 2)
        totalPathwayList = list()
        currentPathwaySet = set()
        currentPathwaySet.add(source)
        extendPathwaySet = set()
        startTime = datetime.now()
        #
        while True:
            #
            # Analyzing...
            if isDebug:
                lCurPwy = str(pathwaylength)
                nCurPwy = str(len(currentPathwaySet))
                sCurPwy = ''.join([str(int(sys.getsizeof(currentPathwaySet)/1024/1024)), 'M'])
                print(lCurPwy + '/' + nCurPwy + '/' + sCurPwy, end="; ")
            #
            # The shortest pathways from source to target have been found.
            if isShortest and (totalPathwayList != list()):
                break
            # The pathway from source to target doesn't exist.
            # Or the whole generalized metabolic space has been searched.
            if currentPathwaySet == set():
                break
            # The maximal value of the pathway length has been reached.
            if pathwaylength >= maxlength:
                break
            # The maximal searching time has been reached.
            interTime = datetime.now()
            runTime = (interTime - startTime).total_seconds()
            if runTime >= maxtime:
                break
            # The maximal number of the retrieved pathways has been reached.
            if len(totalPathwayList) > maxnumber:
                break
            # Deal with out of memory error.
            # # memory = psutil.virtual_memory()
            # # used_memory_percent = memory.percent
            # # used_memory = round(memory.used/1024/1024/1024, 2)
            if (used_memory_percent > mempercent) or (used_memory > memgb):
                break
            #
            pathwaylength += 1
            extendPathwaySet = set()
            #
            # Start the next round of searches.
            for idx, currentPathway in enumerate(currentPathwaySet):
                extPathwaySet = self._extendForward(target, pathwaylength,
                                                    maxlength, isTotal,
                                                    currentPathway,
                                                    totalPathwayList,
                                                    rxnDistInfoDict, compDistInfoDict)
                extendPathwaySet |= extPathwaySet
                # The maximal searching time has been reached.
                if idx%threshold == 0:
                    interTime = datetime.now()
                    runTime = (interTime - startTime).total_seconds()
                    if runTime >= maxtime:
                        break
                    # Deal with out of memory error.
                    memory = psutil.virtual_memory()
                    used_memory_percent = memory.percent
                    used_memory = round(memory.used/1024/1024/1024, 2)
                    if (used_memory_percent > mempercent) or (used_memory > memgb):
                        break
            currentPathwaySet = extendPathwaySet
            if not isDebug:
                pass
                # print("...", end="", flush=True)
        # print("")
        #
        # Calculate the running time.
        endTime = datetime.now()
        runTime = (endTime - startTime).total_seconds()
        if isDebug:
            print("")
            print(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        else:
            self.bfs_throw_prompts(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        #
        # Deal with out of memory error.
        # # memory = psutil.virtual_memory()
        # # used_memory_percent = memory.percent
        # # used_memory = round(memory.used/1024/1024/1024, 2)
        #
        if isDebug:
            print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        else:
            self.bfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        #
        # Sorting...
        # # totalPathwayList = self._lpSort(totalPathwayList)
        #
        # Unifying...
        totalPathwayList = self._unifyPathway(totalPathwayList)
        pathwayNum = len(totalPathwayList)
        #
        if pathwayNum == 0:
            if extendPathwaySet == set():
                if isDebug:
                    print(f"The pathway from {source}({sourceName}) to {target}({targetName}) doesn't exist.")
                    print("Try to use de novo pathway design instead of pathway search.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) doesn't exist.")
                    self.bfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
                return list()
            if (used_memory_percent > mempercent) or (used_memory > memgb):
                if isDebug:
                    print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                    # # print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    print("The memory of this platform restricts the pathway search process.")
                    print("Try to use depth-first search instead of breadth-first search.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                    # # self.bfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    self.bfs_throw_prompts("The memory of this platform restricts the pathway search process.")
                    self.bfs_throw_prompts("Try to use depth-first search instead of breadth-first search.")
                return list()
            if pathwaylength >= maxlength:
                if isTotal:
                    if isDebug:
                        print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                        print("Try to increase the maximal length of the potential pathway and repeat the search process.")
                    else:
                        self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                        self.bfs_throw_prompts("Try to increase the maximal length of the potential pathway and repeat the search process.")
                else:
                    if isDebug:
                        print(f"There are no pathways from {source}({sourceName}) to {target}({targetName}) with exactly {maxlength} steps.")
                        print("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                    else:
                        self.bfs_throw_prompts(f"There are no pathways from {source}({sourceName}) to {target}({targetName}) with exactly {maxlength} steps.")
                        self.bfs_throw_prompts("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                return list()
            if runTime >= maxtime:
                if isDebug:
                    print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                    print("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has not yet been found.")
                    self.bfs_throw_prompts("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                return list()
        #
        #
        if not isShortest:
            if isDebug:
                print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
                print(f"The total number of the retrieved pathways is {pathwayNum}.\n")
            else:
                self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
                self.bfs_throw_prompts(f"The total number of the retrieved pathways is {pathwayNum}.\n")
            return totalPathwayList
        #
        if isDebug:
            print(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
            print(f"The length of the shortest pathways is {pathwaylength}.")
            print(f"The number of the shortest pathways is {pathwayNum}.\n")
        else:
            self.bfs_throw_prompts(f"The pathway from {source}({sourceName}) to {target}({targetName}) has been found.")
            self.bfs_throw_prompts(f"The length of the shortest pathways is {pathwaylength}.")
            self.bfs_throw_prompts(f"The number of the shortest pathways is {pathwayNum}.\n")
        #
        if not isDisplay:
            return totalPathwayList
        #
        numAbbr = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th" ,"8th", "9th", "10th"]
        # KEGG; KndPad;
        pattern_compId = re.compile("(C\d{5}|Met\d{6}-m)")
        for idx, pathway_info in enumerate(totalPathwayList[:6]):
            # # pathway = pathway_info['8-Pathway']
            pathway = pathway_info[8]
            print(f"The {numAbbr[idx]} pathway is:")
            pCompIdList = re.findall(pattern_compId, pathway)
            for compId in pCompIdList:
                compName = compDistInfoDict[compId][0]
                compIdName = "".join([compId, '(', compName, ')'])
                pathway = pathway.replace(compId, compIdName)
            print(pathway)
        if pathwayNum > 6:
            print("......")
        print("")
        #
        return totalPathwayList


    def Source2Target(self, source, target, maxlength = 20, maxtime = 100, maxnumber = 1000000,
                            isTotal = True,
                            isShortest = False,
                            rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                            isDebug = False, isDisplay = False,
                            threshold = 10000,
                            ):
        """
            Search the shortest pathways from source to target by using
          Breadth-First Search algorithm.

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
            print("Search the pathways from source to target (BFS):")
        else:
            # # self.bfs_throw_prompts("")
            self.bfs_throw_prompts("Search the pathways from source to target (BFS):")
        #
        pathwayEnumList = self._pwSource2Target(source, target, maxlength, maxtime, maxnumber,
                                                isTotal,
                                                isShortest,
                                                rxnDistInfoDict, compDistInfoDict,
                                                isDebug, isDisplay,
                                                threshold,
                                                )
        return pathwayEnumList




    def _extendBackward(self, multSourceSet, pathwaylength = 0,
                              maxlength = 20, isTotal = True,
                              currentPathway = '',
                              totalPathwayList = list(),
                              rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                              exDir = 'backward',
                              ):
        """
            Breadth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        extPathwaySet = set()
        #
        currentCompRxnList = currentPathway.split('-->')
        currentCompRxnSet = set(currentCompRxnList[1:])
        lastVisitComp = currentCompRxnList[0]
        #
        # Continue the pathway search through the last visit compound of the current pathway.
        compDistInfo = compDistInfoDict[lastVisitComp]
        # compName = compDistInfo[0]
        # compBf = compDistInfo[1]
        compRxnList = compDistInfo[2]
        newReachRxnSet = set(compRxnList)
        # Exclude the reactions haved been visited by the current pathway.
        newReachRxnSet -= currentCompRxnSet  # Difference set
        # There are no new reactions to extend the current pathway.
        if newReachRxnSet == set():
            return extPathwaySet
        #
        # Continue the pathway search through the newly reached reations.
        for newReachRxn in newReachRxnSet:
            #
            rxnDistInfo = rxnDistInfoDict[newReachRxn]
            rxnDir = rxnDistInfo[0]
            # # rxnCompList = rxnDistInfo[1]
            rxnCompSet = rxnDistInfo[2]
            rxnCompPairList = rxnDistInfo[3][exDir]
            #
            # The information of the reaction's compound pairs doesn't exist.
            if rxnCompPairList is None:
                continue
            #
            # Explicit/implicit visited compounds will be excluded.
            # comp1-->rxn1-->comp2-->rxn2-->comp3-->rxn3-->comp4-->rxn4-->comp5
            #   0       1      2       3      4       5      6       7      8
            if len(rxnCompSet & currentCompRxnSet) != 0:
                continue
            #
            if (rxnDir == '=>') or (rxnDir == '<=>'):
                try:
                    # right-->left
                    newReachCompList = rxnCompPairList[0][lastVisitComp]
                    #
                    # There are no new compounds to extend the current pathway.
                    if len(newReachCompList) == 0:
                        continue
                    #
                    for newReachComp in newReachCompList:
                        # Extend the current pathway.
                        extendPathway = '-->'.join([newReachComp, newReachRxn, currentPathway])
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if newReachComp in multSourceSet:
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((pathwaylength, extendPathway))
                            elif pathwaylength == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((pathwaylength, extendPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                        else:
                            extPathwaySet.add(extendPathway)
                    continue
                except KeyError:
                    pass
            #
            if (rxnDir == '<=') or (rxnDir == '<=>'):
                try:
                    # left-->right
                    newReachCompList = rxnCompPairList[1][lastVisitComp]
                    #
                    # There are no new compounds to extend the current pathway.
                    if len(newReachCompList) == 0:
                        continue
                    #
                    for newReachComp in newReachCompList:
                        # Extend the current pathway.
                        extendPathway = '-->'.join([newReachComp, newReachRxn, currentPathway])
                        # The target compound has been reached, and the new complete pathway will be saved.
                        if newReachComp in multSourceSet:
                            if isTotal:
                                # all pathways not longer than maximum length will be kept.
                                totalPathwayList.append((pathwaylength, extendPathway))
                            elif pathwaylength == maxlength:
                                # pathway equal to maximum length will be kept
                                totalPathwayList.append((pathwaylength, extendPathway))
                            else:
                                # pathway shorter than maximum length will be discarded.
                                pass
                        else:
                            extPathwaySet.add(extendPathway)
                    continue
                except KeyError:
                    pass
        #
        return extPathwaySet


    def _pwTarget2MultSource(self, multSourceSet, target, maxlength = 20, maxtime = 100, maxnumber = 1000000,
                                   isTotal = True,
                                   isShortest = False,
                                   rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                                   isDebug = False, isDisplay = False,
                                   threshold = 10000,
                                   mempercent = 85, memgb = 16,
                                   ):
        """
            Breadth-First Search.

            ARGUMENTS:
                ...
                rxnDistInfoDict: {rxnId: [Direction, Compound List, Compound Set, Compound Pairs], ..., }
                compDistInfoDict: {compId: [Name, Branching Factor, Reactions], ..., }
                ...

            RETURNS:
                ...
        """
        #
        if (not isinstance(multSourceSet, set)) or (not isinstance(target, str)):
            if isDebug:
                print("Incorrect parameter type.")
            else:
                self.bfs_throw_prompts("Incorrect parameter type.")
            return list()
        #
        allCompSet = set(compDistInfoDict.keys())
        if target not in allCompSet:
            if isDebug:
                print(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                print("Try to use de novo pathway design instead of pathway search.")
            else:
                self.bfs_throw_prompts(f"The target({target}) does not exist in the Generalized Metabolic Space.")
                self.bfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
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
                self.bfs_throw_prompts("The valid 'source set' in empty.")
                self.bfs_throw_prompts("Please have a check on the 'source set'.")
            return list()
        if isDebug:
            print("Intermediate pathways' length/number/size:")
        #
        # Initialization
        pathwaylength = 0
        memory = psutil.virtual_memory()
        used_memory_percent = memory.percent
        used_memory = round(memory.used/1024/1024/1024, 2)
        totalPathwayList = list()
        currentPathwaySet = set()
        currentPathwaySet.add(target)
        extendPathwaySet = set()
        startTime = datetime.now()
        #
        while True:
            #
            # Analyzing...
            if isDebug:
                lCurPwy = str(pathwaylength)
                nCurPwy = str(len(currentPathwaySet))
                sCurPwy = ''.join([str(int(sys.getsizeof(currentPathwaySet)/1024/1024)), 'M'])
                print(lCurPwy + '/' + nCurPwy + '/' + sCurPwy, end="; ")
            #
            # The shortest pathways from mult-sources to target have been found.
            if isShortest and (totalPathwayList != list()):
                break
            # The pathway from mult-sources to target doesn't exist.
            # Or the whole generalized metabolic space has been searched.
            if currentPathwaySet == set():
                break
            # The maximal value of the pathway length has been reached.
            if pathwaylength >= maxlength:
                break
            # The maximal searching time has been reached.
            interTime = datetime.now()
            runTime = (interTime - startTime).total_seconds()
            if runTime >= maxtime:
                break
            # The maximal number of the retrieved pathways has been reached.
            if len(totalPathwayList) > maxnumber:
                break
            # Deal with out of memory error.
            # # memory = psutil.virtual_memory()
            # # used_memory_percent = memory.percent
            # # used_memory = round(memory.used/1024/1024/1024, 2)
            if (used_memory_percent > mempercent) or (used_memory > memgb):
                break
            #
            pathwaylength += 1
            extendPathwaySet = set()
            #
            # Start the next round of reverse searches.
            for idx, currentPathway in enumerate(currentPathwaySet):
                extPathwaySet = self._extendBackward(multSourceSet, pathwaylength,
                                                     maxlength, isTotal,
                                                     currentPathway,
                                                     totalPathwayList,
                                                     rxnDistInfoDict, compDistInfoDict,
                                                     )
                extendPathwaySet |= extPathwaySet
                # The maximal searching time has been reached.
                if idx%threshold == 0:
                    interTime = datetime.now()
                    runTime = (interTime - startTime).total_seconds()
                    if runTime >= maxtime:
                        break
                    # Deal with out of memory error.
                    memory = psutil.virtual_memory()
                    used_memory_percent = memory.percent
                    used_memory = round(memory.used/1024/1024/1024, 2)
                    if (used_memory_percent > mempercent) or (used_memory > memgb):
                        break
            currentPathwaySet = extendPathwaySet
            if not isDebug:
                pass
                # print("...", end="", flush=True)
        # print("")
        #
        # Calculate the running time.
        endTime = datetime.now()
        runTime = (endTime - startTime).total_seconds()
        if isDebug:
            print("")
            print(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        else:
            self.bfs_throw_prompts(''.join(['Search Time: ', str(round(runTime, 6)), 's']))
        #
        # Deal with out of memory error.
        # # memory = psutil.virtual_memory()
        # # used_memory_percent = memory.percent
        # # used_memory = round(memory.used/1024/1024/1024, 2)
        #
        if isDebug:
            print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        else:
            self.bfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
        #
        # Sorting...
        # # totalPathwayList = self._lpSort(totalPathwayList)
        #
        # Unifying...
        totalPathwayList = self._unifyPathway(totalPathwayList)
        pathwayNum = len(totalPathwayList)
        #
        multSourceIdNameSet = set()
        for source in list(multSourceSet)[:6]:
            sourceName = compDistInfoDict[source][0]
            sourceIdName = "".join([source, "(", sourceName, ")"])
            multSourceIdNameSet.add(sourceIdName)
        if len(multSourceSet) > 6:
            multSourceIdNameSet.add("...")
        targetName = compDistInfoDict[target][0]
        #
        if pathwayNum == 0:
            if extendPathwaySet == set():
                if isDebug:
                    print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) doesn't exist.")
                    print("Try to use de novo pathway design instead of pathway search.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) doesn't exist.")
                    self.bfs_throw_prompts("Try to use de novo pathway design instead of pathway search.")
                return list()
            if (used_memory_percent > mempercent) or (used_memory > memgb):
                if isDebug:
                    print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                    # # print(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    print("The memory of this platform restricts the pathway search process.")
                    print("Try to use depth-first search instead of breadth-first search.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                    # # self.bfs_throw_prompts(f"Memory: {used_memory}GB({used_memory_percent}%) of RAM of this platform has been used.")
                    self.bfs_throw_prompts("The memory of this platform restricts the pathway search process.")
                    self.bfs_throw_prompts("Try to use depth-first search instead of breadth-first search.")
                return list()
            if pathwaylength >= maxlength:
                if isTotal:
                    if isDebug:
                        print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                        print("Try to increase the maximal length of the potential pathway and repeat the search process.")
                    else:
                        self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                        self.bfs_throw_prompts("Try to increase the maximal length of the potential pathway and repeat the search process.")
                else:
                    if isDebug:
                        print(f"There are no pathways from {multSourceIdNameSet} to {target}({targetName}) with exactly {maxlength} steps.")
                        print("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                    else:
                        self.bfs_throw_prompts(f"There are no pathways from {multSourceIdNameSet} to {target}({targetName}) with exactly {maxlength} steps.")
                        self.bfs_throw_prompts("Try to decrease or increase the length of the potential pathway and repeat the search process.")
                return list()
            if runTime >= maxtime:
                if isDebug:
                    print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                    print("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                else:
                    self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has not yet been found.")
                    self.bfs_throw_prompts("Try to increase the maximal searching time of the potential pathway and repeat the search process.")
                return list()
        #
        #
        multSourceSet = set()
        for pathway_info in totalPathwayList:
            # # pathway = pathway_info['8-Pathway']
            pathway = pathway_info[8]
            source = pathway.split('-->', 1)[0]
            multSourceSet.add(source)
        multSourceIdNameSet = set()
        for source in multSourceSet:
            sourceName = compDistInfoDict[source][0]
            sourceIdName = "".join([source, "(", sourceName, ")"])
            multSourceIdNameSet.add(sourceIdName)
        targetName = compDistInfoDict[target][0]
        #
        if not isShortest:
            if isDebug:
                print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
                print(f"The total number of the retrieved pathways is {pathwayNum}.\n")
            else:
                self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
                self.bfs_throw_prompts(f"The total number of the retrieved pathways is {pathwayNum}.\n")
            return totalPathwayList
        #
        if isDebug:
            print(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
            print(f"The length of the shortest pathways is {pathwaylength}.")
            print(f"The number of the shortest pathways is {pathwayNum}.\n")
        else:
            self.bfs_throw_prompts(f"The pathway from {multSourceIdNameSet} to {target}({targetName}) has been found.")
            self.bfs_throw_prompts(f"The length of the shortest pathways is {pathwaylength}.")
            self.bfs_throw_prompts(f"The number of the shortest pathways is {pathwayNum}.\n")
        #
        if not isDisplay:
            return totalPathwayList
        #
        numAbbr = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th" ,"8th", "9th", "10th"]
        # KEGG; KndPad;
        pattern_compId = re.compile("(C\d{5}|Met\d{6}-m)")
        for idx, pathway_info in enumerate(totalPathwayList[:6]):
            # # pathway = pathway_info['8-Pathway']
            pathway = pathway_info[8]
            print(f"The {numAbbr[idx]} pathway is:")
            pCompIdList = re.findall(pattern_compId, pathway)
            for compId in pCompIdList:
                compName = compDistInfoDict[compId][0]
                compIdName = "".join([compId, '(', compName, ')'])
                pathway = pathway.replace(compId, compIdName)
            print(pathway)
        if pathwayNum > 6:
            print("......")
        print("")
        #
        return totalPathwayList


    def Target2MultSource(self, multSourceSet, target, maxlength = 20, maxtime = 100, maxnumber = 1000000,
                                isTotal = True,
                                isShortest = False,
                                rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                                isDebug = False, isDisplay = False,
                                threshold = 10000,
                                ):
        """
            Search the shortest pathways from target to mult-source by using
          Breadth First Search algorithm.

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
            print("Search the pathways from target to mult-sources (BFS):")
        else:
            # # self.bfs_throw_prompts("")
            self.bfs_throw_prompts("Search the pathways from target to mult-sources (BFS):")
        #
        pathwayEnumList = self._pwTarget2MultSource(multSourceSet, target, maxlength, maxtime, maxnumber,
                                                    isTotal,
                                                    isShortest,
                                                    rxnDistInfoDict, compDistInfoDict,
                                                    isDebug, isDisplay,
                                                    threshold,
                                                    )
        return pathwayEnumList


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





