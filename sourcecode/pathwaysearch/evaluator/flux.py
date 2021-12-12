# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 19:32:10 2021

@author: CC-SXF
"""

import re
import copy
from cobra import Metabolite, Reaction
from cobra.util.solver import linear_reaction_coefficients

# import math
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
# # from concurrent.futures import ProcessPoolExecutor
# from concurrent.futures import as_completed


class Flux():
    """
    """
    def __init__(self, ):
        """ """
        pass


    @classmethod
    def _getSupp(cls,):
        """
        """
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
        supp_compId_set = (endo_kegg_compId_set | endo_metacyc_compId_set | endo_kndpad_compId_set)
        #
        return supp_compId_set


    @classmethod
    def _deEqu(cls, equ, met_id="", rxn_bounds=(-1000.0, 1000.0), sep=" @#$% ", isDebug=False):
        """
        """
        # -->/<=>/<--
        pattern_equ = re.compile("( --> | <=> | <-- )")
        rxn_dir = re.search(pattern_equ, equ).group(1)
        #
        equ = re.sub(rxn_dir, sep, equ)
        (equ_left, equ_right) = equ.split(sep)
        #
        coef_left = list()
        coef_right = list()
        comp_left = list()
        comp_right = list()
        #
        coef2comp_left_list = equ_left.split(" + ")
        for coef2comp in coef2comp_left_list:
            try:
                (coef, comp) = coef2comp.split()
                coef_left.append( float(coef) )
                comp_left.append( comp )
            except ValueError:
                coef_left.append(1.0)
                comp_left.append( coef2comp.strip() )
        #
        if equ_right != "":
            coef2comp_right_list = equ_right.split(" + ")
            for coef2comp in coef2comp_right_list:
                try:
                    (coef, comp) = coef2comp.split()
                    coef_right.append( float(coef) )
                    comp_right.append( comp )
                except ValueError:
                    coef_right.append(1.0)
                    comp_right.append( coef2comp.strip() )
        #
        # unify reaction's equation
        if (met_id in comp_left) and (len(comp_right) != 0):
            coef_left, coef_right = coef_right, coef_left
            comp_left, comp_right = comp_right, comp_left
            rxn_dir_dict = {" --> ": " <-- ", " <-- ": " --> ", " <=> ": " <=> ", }
            rxn_dir = rxn_dir_dict[rxn_dir]
            (lower_bound, upper_bound) = rxn_bounds
            rxn_bounds = (-upper_bound, -lower_bound)
        #
        coef_list = [coef_left, coef_right]
        comp_list = [comp_left, comp_right]
        equ_coef_comp_dir_bd = [coef_list, comp_list, rxn_dir, rxn_bounds]
        #
        if isDebug:
            print(equ_coef_comp_dir_bd, end="\n\n")
        #
        return equ_coef_comp_dir_bd
        #


    @classmethod
    def _enEqu(cls, equ_coef_comp_dir):
        """
        """
        [coef_list, comp_list, rxn_dir] = equ_coef_comp_dir
        [coef_left, coef_right] = coef_list
        [comp_left, comp_right] = comp_list
        #
        coef2comp_left_list = list()
        for coef, comp in zip(coef_left, comp_left):
            coef2comp = "".join([str(coef), " ", comp])
            coef2comp_left_list.append(coef2comp)
        coef2comp_right_list = list()
        for coef, comp in zip(coef_right, comp_right):
            coef2comp = "".join([str(coef), " ", comp])
            coef2comp_right_list.append(coef2comp)
        #
        equ_left = " + ".join(coef2comp_left_list)
        equ_right = " + ".join(coef2comp_right_list)
        equ = rxn_dir.join([equ_left, equ_right])
        #
        return equ


    @classmethod
    def _unify_info(cls, pathway_info,
                         rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                         organism = "eco", org_model = None, orgCompRxnId_dict = dict(),
                         isDebug = False,
                         ):
        """
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
        comp_key = "".join(["Compound", "(", organism, ")"])
        reac_key = "".join(["Reaction", "(", organism, ")"])
        #
        comp2rxn_list = pathway_info[8].split('-->')
        i_details = pathway_info[9]
        #
        new_met_set = set()
        new_rxn_dict = dict()
        endo_rxn_set = set()
        #
        for idx_cr, comp2rxn in enumerate(comp2rxn_list):
            try:
                # reaction:
                rxn_id = comp2rxn
                # # direction = rxnDistInfoDict[rxn_id][0]                       # Direction
                coefficient_list = copy.deepcopy(rxnDistInfoDict[rxn_id][6])     # Coefficient List
                compound_list = copy.deepcopy(rxnDistInfoDict[rxn_id][1])        # Compound List
                #
                isEndogenous = i_details[(idx_cr-1)//2]
                if isEndogenous == 1:
                    # endogenous reaction
                    rxn_id = orgCompRxnId_dict[reac_key][rxn_id]
                    comp2rxn_list[idx_cr] = rxn_id
                    endo_rxn_set.add( rxn_id )
                    continue
                else:
                    # exogenous reaction
                    pass
                #
                # update the compounds' id of this reaction according to the organism we select.
                # record new compounds' id, and then add them to the organism we use.
                for idx_c, comp_id in enumerate(compound_list[0]):
                    try:
                        # endogenous metabolite
                        met_id = orgCompRxnId_dict[comp_key][comp_id]
                        compound_list[0][idx_c] = met_id
                    except KeyError:
                        # exogenous metabolite
                        new_met_set.add(comp_id)
                        pass
                for idx_c, comp_id in enumerate(compound_list[1]):
                    try:
                        # endogenous metabolite
                        met_id = orgCompRxnId_dict[comp_key][comp_id]
                        compound_list[1][idx_c] = met_id
                    except KeyError:
                        # exogenous metabolite
                        new_met_set.add(comp_id)
                #
                # equation reconstitution
                coefcompid_left_list = list()
                for coef, comp_id in zip(coefficient_list[0], compound_list[0]):
                    coefcompid_left_list.append( " ".join( [str(coef), comp_id] ) )
                coefcompid_right_list = list()
                for coef, comp_id in zip(coefficient_list[1], compound_list[1]):
                    coefcompid_right_list.append( " ".join( [str(coef), comp_id] ) )
                #
                # equation reconstitution
                equ_left = " + ".join(coefcompid_left_list)
                equ_right = " + ".join(coefcompid_right_list)
                equ = " <=> ".join([equ_left, equ_right])
                #
                new_rxn_dict[rxn_id] = equ
                #
            except KeyError:
                # compound:
                comp_id = comp2rxn
                try:
                    # endogenous metabolite
                    # KEGG/MetaCyc/KndPad Compound id --> BIGG Metabolite id
                    comp2rxn_list[idx_cr] = orgCompRxnId_dict[comp_key][comp_id]
                except KeyError:
                    # exogenous metabolite
                    pass
        #
        # target
        target_id = comp2rxn_list[-1]
        try:
            target_id = orgCompRxnId_dict[comp_key][target_id]
        except KeyError:
            pass
        #
        # irreversible reaction
        new_rxn_dict["Target_Product"] = "".join( [target_id, " --> "] )
        #
        return (new_met_set, new_rxn_dict, endo_rxn_set, comp2rxn_list)


    @classmethod
    def _evalOneFlux(cls, pathway_info,
                          rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                          organism = "eco", org_model = None, orgCompRxnId_dict = dict(),
                          supp_compId_set = set(),
                          isDebug = False,
                          ):
        """
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
        isFeasibility = pathway_info[0]
        if not isFeasibility:
            return pathway_info
        #
        isBalance = True
        pathway = pathway_info[8]
        comp2rxn_list = pathway.split('-->')
        for comp2rxn in comp2rxn_list:
            try:
                status = rxnDistInfoDict[comp2rxn][5]  # Status
                if status == "Unbalanced":
                    isBalance = False
                    break
            except KeyError:
                pass
        if not isBalance:
            return pathway_info
        #
        unify_info = cls._unify_info(pathway_info,
                                     rxnDistInfoDict, compDistInfoDict,
                                     organism, None, orgCompRxnId_dict,
                                     isDebug,
                                     )
        (new_met_set, new_rxn_dict, endo_rxn_set, comp2rxn_list) = unify_info
        #
        #
        with org_model:
            #
            # add new compounds
            new_model_met_list = list()
            for met_id in new_met_set:
                # Metabolite(id=None, formula=None, name='', charge=None, compartment=None)
                new_model_met_list.append( Metabolite(met_id, formula=None, name=compDistInfoDict[met_id][0], charge=None, compartment=None) )
            org_model.add_metabolites(new_model_met_list)
            #
            # add new reactions
            new_model_rxn_list = list()
            for rxn_id in new_rxn_dict.keys():
                # Reaction(id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None)
                new_model_rxn_list.append( Reaction(rxn_id, name='', subsystem='', lower_bound=-1000.0, upper_bound=1000.0) )
            org_model.add_reactions(new_model_rxn_list)
            #
            # update the equation of the new reactions added to the model.
            for new_model_rxn in new_model_rxn_list:
                rxn_id = new_model_rxn.id
                equ = new_rxn_dict[rxn_id]
                new_model_rxn.reaction = equ
                if rxn_id == "Target_Product":
                    new_model_rxn.bounds = (0.0, 1000.0)
                else:
                    new_model_rxn.bounds = (-1000.0, 1000.0)
            #
            # all metabolites(endogenous & exogenous)
            org_met_dict = dict()
            for org_met in org_model.metabolites:
                met_id = org_met.id
                org_met_dict[met_id] = org_met
            #
            # all reactions(endogenous & exogenous)
            org_rxn_dict = dict()
            org_rxnid_dict = dict()
            for org_rxn in org_model.reactions:
                rxn_id = org_rxn.id
                org_rxn_dict[rxn_id] = org_rxn
                org_rxnid_dict[rxn_id] = 0
            # objective reaction
            obj_rxn_id = list(linear_reaction_coefficients(org_model).keys())[0].id
            org_rxnid_dict[obj_rxn_id] = 1
            #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            #
            # endogenous reaction in this pathway  --> reversible
            # reconcile zero-flux of this pathway derived from inaccurate reaction direction of this micro-organism model.
            for rxn_id in endo_rxn_set:
                org_rxn = org_rxn_dict[rxn_id]
                org_rxn.bounds = (-1000.0, 1000.0)
            #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            #
            # supply uncerntain compound to the culture media
            pwy_total_met_set = set()
            for idx_cr, comp2rxn in enumerate(comp2rxn_list):
                # metabolite
                if (idx_cr%2) == 0:
                    continue
                # reaction
                rxn_id = comp2rxn
                org_rxn = org_rxn_dict[rxn_id]
                pwy_total_met_set |= set(org_rxn.metabolites)
            #
            supp_rxn_list = list()
            supp_rxn_dict = dict()
            for org_met in pwy_total_met_set:
                met_id = org_met.id
                # uncertain endogenous or heterogenous compound
                if met_id in supp_compId_set:
                    supp_rxn_id = "".join(["Tr_", met_id])
                    supp_rxn = Reaction(supp_rxn_id, name='', subsystem='', lower_bound=-1000.0, upper_bound=1000.0)
                    supp_rxn_list.append(supp_rxn)
                    supp_equ = "".join([met_id, " ", "<=>"])
                    supp_rxn_dict[supp_rxn_id] = supp_equ
            org_model.add_reactions(supp_rxn_list)
            #
            # relax the lower and upper bound of uncertain endogenous or heterogenous compound
            for supp_rxn in supp_rxn_list:
                supp_rxn_id = supp_rxn.id
                supp_equ = supp_rxn_dict[supp_rxn_id]
                supp_rxn.reaction = supp_equ
                supp_rxn.bounds = (-1000.0, 1000.0)
            #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            #
            # add new constraints to maintain main flux flowing through this pathway
            for idx_cr, comp2rxn in enumerate(comp2rxn_list):
                # source or reaction
                if (idx_cr == 0) or ((idx_cr%2) == 1):
                    continue
                # exogenous metabolite
                met_id = comp2rxn
                if met_id in new_met_set:
                    continue
                #
                # compartments
                # c: cytosol / cytoplasm
                # m: mitochondria / mitochondrion
                # n: nucleus
                # r: endoplasmic reticulum
                # g: golgi apparatus / golgi
                # v: vacuole
                # x: peroxisome / glyoxysome
                # p: periplasm / periplasmic space
                # e: extracellular space / extra-organism / extraorganism / extracellular
                compartment_set = {'c', 'm', 'n', 'r', 'g', 'v', 'x', 'p', 'e'}
                met_id_header = met_id[:-1]
                met_id_set = {"".join([met_id_header, c]) for c in compartment_set}
                #
                # the reaction that produces this metabolite
                pre_rxn_id = comp2rxn_list[idx_cr-1]
                org_rxn = org_rxn_dict[pre_rxn_id]
                rxn_equ = org_rxn.reaction
                rxn_met_id_set = set( re.findall("\S*_[cmnrgvxpe]", rxn_equ) )
                #
                try:
                    met_id = list(met_id_set & rxn_met_id_set)[0]
                except IndexError:
                    # pseudo-endogenous metabolite
                    # conflicting mapping of compound and reaction between database and organism model
                    continue
                #
                # endogenous metabolite
                org_met = org_met_dict[met_id]
                # reactions this metabolite participates in
                org_met2rxn_list = list(org_met.reactions)
                org_met2rxn_num = len(org_met2rxn_list)
                # This metabolite only participates in two reaction.
                # One produces this metabolite, and the other consumes it
                if org_met2rxn_num <= 2:
                    continue
                #
                # add pseudo-metabolites to the model of this organism
                pseudo_model_met_list = list()
                pseudo_model_met_id_list = list()
                for idx_pcr in range(org_met2rxn_num):
                    pseudo_met_id = "_".join(["pmet", str(idx_cr).rjust(3, "0"), str(idx_pcr+1).rjust(4, "0")])
                    # Metabolite(id=None, formula=None, name='', charge=None, compartment=None)
                    pseudo_model_met_list.append( Metabolite(pseudo_met_id, formula=None, name='', charge=None, compartment=None) )
                    pseudo_model_met_id_list.append(pseudo_met_id)
                org_model.add_metabolites(pseudo_model_met_list)
                #
                # add pseudo-reactions to the model of this organism
                pseudo_model_rxn_list = list()
                for idx_pcr in range(org_met2rxn_num):
                    pseudo_rxn_id = "_".join(["prxn", str(idx_cr).rjust(3, "0"), str(idx_pcr+1).rjust(4, "0")])
                    # Reaction(id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None)
                    pseudo_model_rxn_list.append( Reaction(pseudo_rxn_id, name='', subsystem='', lower_bound=0.0, upper_bound=1000.0) )
                org_model.add_reactions(pseudo_model_rxn_list)
                #
                # update the equation of the pseudo-reactions added to the model.
                for idx_pcr, pseudo_model_rxn in enumerate(pseudo_model_rxn_list):
                    pseudo_met_id = pseudo_model_met_id_list[idx_pcr]
                    pseudo_equ = "".join( [pseudo_met_id, " --> "] )
                    pseudo_model_rxn.reaction = pseudo_equ
                    pseudo_model_rxn.bounds = (0.0, 1000.0)
                #
                pseudo_model_met_id_list_copy = copy.deepcopy(pseudo_model_met_id_list)
                #
                # * * * * * * * * * * * * * * * * * * #
                # reactions this metabolite participates in
                for org_rxn in org_met2rxn_list:
                    rxn_id = org_rxn.id
                    #
                    # primary endogenous or exogenous reation in this pathway
                    if rxn_id in comp2rxn_list:
                        # delete the pseudo-metabolite which has been added to the complementary reaction
                        _ = pseudo_model_met_id_list_copy.pop(0)
                        #
                        rxn_bounds = org_rxn.bounds
                        rxn_equ = org_rxn.reaction
                        # decompose the equation of this reaction and move the main metabolite to the right of this equation
                        [rxn_coef, rxn_comp, rxn_dir, rxn_bounds] = cls._deEqu(rxn_equ, met_id, rxn_bounds)
                        # left, reactant; right, product; without additional constraint
                        if (rxn_id == pre_rxn_id) and (met_id in rxn_comp[1]) and (not org_rxnid_dict[rxn_id]):
                            met_index = rxn_comp[1].index(met_id)
                            met_coef = rxn_coef[1][met_index]
                            # coefficient normalization
                            for i, coef in enumerate(rxn_coef[0]):
                                rxn_coef[0][i] = round(coef/met_coef, 7)
                            for j , coef in enumerate(rxn_coef[1]):
                                rxn_coef[1][j] = round(coef/met_coef, 7)
                            # add pseudo-metabolites to the reaction
                            for pseudo_met_id in pseudo_model_met_id_list:
                                rxn_coef[1].append(1.0)
                                rxn_comp[1].append(pseudo_met_id)
                            # update this reaction's equation and recover its bounds
                            equ_coef_comp_dir = [rxn_coef, rxn_comp, rxn_dir]
                            rxn_equ_pseudo = cls._enEqu(equ_coef_comp_dir)
                            org_rxn.reaction = rxn_equ_pseudo
                            org_rxn.bounds = rxn_bounds
                            # additional constraint has been added to this reaction
                            org_rxnid_dict[rxn_id] = 1
                    # secondary endogenous reation in this organism
                    else:
                        # delete the pseudo-metabolite which has been added to the complementary reaction
                        pseudo_met_id = pseudo_model_met_id_list_copy.pop(0)
                        #
                        rxn_bounds = org_rxn.bounds
                        rxn_equ = org_rxn.reaction
                        # decompose the equation of this reaction and move the main metabolite to the right of this equation
                        [rxn_coef, rxn_comp, rxn_dir, rxn_bounds] = cls._deEqu(rxn_equ, met_id, rxn_bounds)
                        # right, product; without additional constraint
                        if (met_id in rxn_comp[1]) and (not org_rxnid_dict[rxn_id]):
                            met_index = rxn_comp[1].index(met_id)
                            met_coef = rxn_coef[1][met_index]
                            # coefficient normalization
                            for i, coef in enumerate(rxn_coef[0]):
                                rxn_coef[0][i] = round(coef/met_coef, 7)
                            for j , coef in enumerate(rxn_coef[1]):
                                rxn_coef[1][j] = round(coef/met_coef, 7)
                            # add pseudo-metabolite to this reaction
                            rxn_coef[0].append(1.0)
                            rxn_comp[0].append(pseudo_met_id)
                            # update this reaction's equation and recover its bounds
                            equ_coef_comp_dir = [rxn_coef, rxn_comp, rxn_dir]
                            rxn_equ_pseudo = cls._enEqu(equ_coef_comp_dir)
                            org_rxn.reaction = rxn_equ_pseudo
                            org_rxn.bounds = rxn_bounds
                            # additional constraint has been added to this reaction
                            org_rxnid_dict[rxn_id] = 1
                # * * * * * * * * * * * * * * * * * * #
                #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            #
            #
            # calculate the flux
            org_model.objective = {org_model.reactions.Target_Product: 1}
            target_flux = org_model.slim_optimize()
        #
        # metabolic flux
        if type(target_flux) == float:
            pathway_info[7] = round(target_flux, 2)
        #
        return pathway_info


    @classmethod
    def _evalFlux(cls, pathway_list,
                       rxnDistInfoDict = dict(), compDistInfoDict = dict(),
                       organism = "eco", org_model = None, orgCompRxnId_dict = dict(),
                       isDebug = False,
                       ):
        """
            Calculate the flux of every pathway by using cobrapy.
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
        # """
        #
        supp_compId_set = cls._getSupp()
        #
        for idx_p, pathway_info in enumerate(pathway_list):
            cls._evalOneFlux(pathway_info,
                             rxnDistInfoDict , compDistInfoDict,
                             organism, org_model, orgCompRxnId_dict,
                             supp_compId_set,
                             isDebug,
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
            futures = {executor.submit(cls._evalOneFlux, pathway_info,
                                                         rxnDistInfoDict, compDistInfoDict,
                                                         organism, org_model, orgCompRxnId_dict,
                                                         supp_compId_set,
                                                         isDebug):pathway_info[8] for pathway_info in pathway_list}
            for future in as_completed(futures):
                s_details = futures[future]
                try:
                    future.result()
                    # pathway_info = future.result()
                except Exception as excReason:
                    print("\nWarning: \nPathway: %s generated an exception: %s\n" % (s_details, excReason))
        # """
        #
        return pathway_list


    @classmethod
    def demoFunc(cls, ):
        """ """


if __name__ == "__main__":
    """ """
    pass





