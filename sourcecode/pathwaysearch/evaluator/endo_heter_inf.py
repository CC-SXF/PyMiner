# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 21:49:27 2021

@author: CC-SXF
"""

# import math
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
# # from concurrent.futures import ProcessPoolExecutor
# from concurrent.futures import as_completed


class EvalEHI():
    """
    """
    def __init__(self, ):
        """ """
        pass


    @classmethod
    def _evalOneEhi(cls, pathway_info,
                         rxnDistInfoDict = dict(),
                         modify_compId_set = set(), org_rxnId_set = set(),
                         isInfeasible = False,
                         ):
        """
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
        pathway = pathway_info[8]
        pCompRxnList = pathway.split('-->')
        pathwaylength = len(pCompRxnList)//2
        endo_met_set = modify_compId_set.copy()
        for idx in range(pathwaylength):
            compId = pCompRxnList[2*idx]
            rxnId = pCompRxnList[2*idx+1]
            rxnCompList = rxnDistInfoDict[rxnId][1]
            leftCompSet = set(rxnCompList[0])
            rightCompSet = set(rxnCompList[1])
            if compId in leftCompSet:
                diffCompSet = (leftCompSet - endo_met_set)
                if len(diffCompSet) == 0:
                    # This is a biochemically plausible reaction, considering the microbial chassis we choose.
                    if rxnId in org_rxnId_set:
                        # This is an endogenous reaction.
                        pathway_info[9][idx] = 1
                    else:
                        # This is an exogenous reaction.
                        pathway_info[9][idx] = 0
                else:
                    # This is a biochemically implausible reaction, considering the microbial chassis we choose.
                    pathway_info[9][idx] = -1
                    pathway_info[10][idx] = list(diffCompSet)
                endo_met_set |= rightCompSet
            elif compId in rightCompSet:
                diffCompSet = (rightCompSet - endo_met_set)
                if len(diffCompSet) == 0:
                    # This is a biochemically plausible reaction, considering the microbial chassis we choose.
                    if rxnId in org_rxnId_set:
                        # This is an endogenous reaction.
                        pathway_info[9][idx] = 1
                    else:
                        # This is an exogenous reaction.
                        pathway_info[9][idx] = 0
                else:
                    # This is a biochemically implausible reaction, considering the microbial chassis we choose.
                    pathway_info[9][idx] = -1
                    pathway_info[10][idx] = list(diffCompSet)
                endo_met_set |= leftCompSet
            else:
                pass
        #
        pathway_info[2] = pathway_info[9].count(1)          # Endogenous Pathway Length
        pathway_info[3] = pathway_info[9].count(0)          # Heterogeneous Pathway Length
        pathway_info[4] = pathway_info[9].count(-1)         # Infeasible Pathway Length
        if pathway_info[4] != 0:
            pathway_info[0] = False
        #
        #
        # '11-Route'                                   <----
        # The relative route of this pathway's highlighted compound structure figures appended in atom transfer evaluating phase.
        pathway_info.append("")
        #
        #
        if isInfeasible:
            # All pathways will be kept, including biochemically implausible ones.
            return pathway_info
        elif pathway_info[0] == True:
            # Only biochemically plausible pathways will be kept.
            return pathway_info
        else:
            return None


    @classmethod
    def _evalEndoHeterInf(cls, unifyPathwayList,
                               rxnDistInfoDict = dict(),
                               org_compId_set = set(), org_rxnId_set = set(),
                               isInfeasible = False,
                               ):
        """
            Validate biochemical plausibility of the reactions in every pathway based on
        the  microbial chassis we choose.
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
        # Uncertain endogenous or heterogenous compound,  and we need further analysis manually.
        #
        endo_kegg_compId_set = {'C03024', "C03161",  'C00028', "C00030",  'C00138', "C00139",
                                'C00996', "C00999",  'C00125', "C00126",  'C15602', "C15603",
                                'C00342', "C00343",  'C00662', "C00667",  'C04253', "C04570",
                                'C02745', "C02869",  'C22150', "C22151",  'C00390', "C00399",
                                }
        #
        endo_metacyc_compId_set = {'Met000018-m', "Met000019-m",  'Met000363-m', "Met000364-m",         # 'C03024' / "C03161"
                                   'Met000017-m', "Met000020-m",  'Met000017-m', "Met000652-m",         # 'C00028' / "C00030"
                                   'Met000029-m', "Met000030-m",                                        # 'C00138' / "C00139"
                                   'Met000042-m', "Met000043-m",                                        # 'C00996' / "C00999"
                                   'Met000107-m', "Met000108-m",                                        # 'C00125' / "C00126"
                                   'Met000684-m', "Met001058-m",  'Met000076-m', "Met000071-m",         # 'C15602' / "C15603"
                                   'Met000129-m', "Met000130-m",                                        # 'C00342' / "C00343"
                                   'Met000153-m', "Met000154-m",                                        # 'C00662' / "C00667"
                                   'Met000069-m', "Met000070-m",                                        # 'C04253' / "C04570"
                                   'Met000220-m', "Met000158-m",                                        # 'C02745' / "C02869"
                                                                                                        # 'C22150' / "C22151"
                                   'Met000159-m', "Met000155-m",                                        # 'C00390' / "C00399"
                                   'Met000162-m', "Met000167-m",                                        #
                                   'Met000027-m', "Met000028-m",                                        # NAD(P)+  /  NAD(P)H

                                   }
        #
        endo_kndpad_compId_set = set()
        #
        modify_compId_set = (org_compId_set | endo_kegg_compId_set | endo_metacyc_compId_set | endo_kndpad_compId_set)
        #
        # """
        pathway_list = list()
        for pathway_info in unifyPathwayList:
            pathway_info = cls._evalOneEhi(pathway_info,
                                           rxnDistInfoDict,
                                           modify_compId_set, org_rxnId_set,
                                           isInfeasible,
                                           )
            if pathway_info != None:
                pathway_list.append(pathway_info)
        # """
        #
        """
        # parallelism --> concurrency
        pathway_list = list()
        max_workers = multiprocessing.cpu_count()
        max_workers = math.ceil(max_workers * 0.8)
        # Use a "with statement" to ensure threads are cleaned up promptly
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Start the load operations and mark each future with pathway's string.
            futures = {executor.submit(cls._evalOneEhi, pathway_info,
                                                        rxnDistInfoDict,
                                                        modify_compId_set,org_rxnId_set,
                                                        isInfeasible):pathway_info[8] for pathway_info in unifyPathwayList}
            for future in as_completed(futures):
                s_details = futures[future]
                try:
                    pathway_info = future.result()
                    if pathway_info != None:
                        pathway_list.append(pathway_info)
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





