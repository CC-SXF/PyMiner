# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 22:43:31 2021

@author: CC-SXF
"""

from pathwaysearch.pathwaysearch import PathwaySearch


class Input():
    """
    """
    source =  {"C00079", "C00082"}
    target = "C03582"
    avoidCompSet = set()
    avoidRxnSet = set()

    organism = "eco"            # "", "bsu", "eco", "kpn", "llm", "ppu", "sce", "syz"
    dataBase = "KEGG"           # "KEGG", "MetaCyc", "KndPad"
    maxlength = 4               # 1, 2, ..., 16
    maxtime = 120               # 1, 2, ..., 6000
    # maxnumber = 1000000

    diffThreshold = 0.1         # 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    methods = "bfs_naive"       # "bfs_naive", "dfs_naive"

    isTotal = True              # True, False
    isShortest = False          # True, False

    isRetro = True              # True, False
    isSmart = True              # True, False

    isInfeasible = False        # True, False
    isTransfer = True           # True, False
    isAerobic = True            # True, False
    isFlux = True               # True, False

    isUpdate = True            # True, False

    if isUpdate:
        while True:
            update_confirm = input('Do you want to update the distilled metabolic network used for pathway search ? <Yes/No>:');
            # update_confirm = update_confirm.title()
            if update_confirm == "Yes":
                isUpdate = True
                print("Warning: THE DISTILLED METABOLIC NETWORK USED FOR PATHWAY SEARCH WILL BE UPDATED.")
                break
            elif update_confirm == "No":
                isUpdate = False
                break
            else:
                continue

    if dataBase == "KEGG":
        resource_file = "../sourcedata/KEGG_Ori_Cleaned_Half.xlsx"
        general_cofactors_file = "./datas/KEGG_General_Cofactors.xlsx"
        metabolic_network_file = "../sourcedata/KEGG_Pathway_Search_Ori.xlsx"
        atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/"
    elif dataBase == "MetaCyc":
        resource_file = "../sourcedata/KndPad_Ori_Cleaned_Half.xlsx"
        general_cofactors_file = "./datas/MetaCyc_General_Cofactors.xlsx"
        metabolic_network_file = "../sourcedata/MetaCyc_Pathway_Search_Ori.xlsx"
        atom_mapping_file_location = "../KndPadAAMBCsRCsFiles_Ori_Half/"
    elif dataBase == "KndPad":
        resource_file = "../sourcedata/KndPad_Ori_Cleaned_Half.xlsx"
        general_cofactors_file = "./datas/KndPad_General_Cofactors.xlsx"
        metabolic_network_file = "../sourcedata/KndPad_Pathway_Search_Ori.xlsx"
        atom_mapping_file_location = "../KndPadAAMBCsRCsFiles_Ori_Half/"
    else:
        pass
    # org_model_location = "../sourcedata/models/"
    # isDebug = False


def search():
    """
    """
    ip = Input()

    ps_instance = PathwaySearch()

    (pathways,_,_) = ps_instance.searchPathway(ip.source, ip.target,  avoidCompSet = ip.avoidCompSet, avoidRxnSet = ip.avoidRxnSet,
                                               organism = ip.organism, dataBase = ip.dataBase, maxlength = ip.maxlength, maxtime = ip.maxtime,
                                               maxnumber = 1000000,
                                               diffThreshold = ip.diffThreshold, methods = ip.methods,
                                               isTotal = ip.isTotal, isShortest = ip.isShortest,
                                               isRetro = ip.isRetro, isSmart = ip.isSmart,
                                               isInfeasible = ip.isInfeasible, isTransfer = ip.isTransfer,
                                               isAerobic = ip.isAerobic, isFlux = ip.isFlux,
                                               isUpdate = ip.isUpdate,
                                               resource_file = ip.resource_file,
                                               general_cofactors_file = ip.general_cofactors_file,
                                               metabolic_network_file = ip.metabolic_network_file,
                                               atom_mapping_file_location = ip.atom_mapping_file_location,
                                               org_model_location = "../sourcedata/models/",
                                               isDebug = False,
                                               )
    return pathways


if __name__ == "__main__":
    """
    """
    pathways = search()





