# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:34:25 2020

@author: CC-SXF
"""

from os import path
from os import listdir
# from os import mkdir
import pickle
# import sys

from pathwaysearch.input.input_kegg import InputKegg
# # from pathwaysearch.input.input_rhea import InputRhea
from pathwaysearch.input.input_metacyc import InputMetaCyc
from pathwaysearch.input.input_kndpad import InputKndPad


class Input():
    """
        Tools for preparing datas for pathway search.
    """
    def __init__(self):
        """ """
        pass


    @classmethod
    def generateRxnCompPair(cls, resource_file = "../sourcedata/KEGG_Ori_Cleaned_Half.xlsx",
                                 general_cofactors_file = "datas/KEGG_General_Cofactors.xlsx",
                                 metabolic_network_file = "../sourcedata/KEGG_Pathway_Search_Ori.xlsx",
                                 atom_mapping_file_location = "../KEGGAAMBCsRCsFiles_Ori_Half/",
                                 org_model_location = "../sourcedata/models/",
                                 dataBase = "KEGG",
                                 isUpdate = True,
                                 diffThresholdLevel = 10,
                            ):
        """
            Generate the substrate-product pairs(reactant paris) for pathway search.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        [_, resource_xlsxFile] = path.split(resource_file)
        [dirDatas, general_cofactors_xlsxFile] = path.split(general_cofactors_file)
        [dirSourceData, xlsxFile] = path.split(metabolic_network_file)
        dirRxnFile = atom_mapping_file_location
        #

        if dataBase == "KEGG":
            if isUpdate or (xlsxFile not in listdir(dirSourceData)):
                #
                # '''
                # Distill the information from KEGG.
                # The fields of the information of compounds and reactions is given as follows.
                # Compounds:
                #  0      1      2        3           4           5       6             7      8        9
                # Entry; Names; Formula; Exact Mass; Mol Weight; SMILES; Reactions;    ChEBI; PubChem; Reference;
                # Reactions:
                #  0      1      2           3         4            5         6       7       8        9             10    11
                # Entry; Names; Definition; Equation; Coefficient; Compound; SMILES; Status; Comment; EC Number;    Rhea; Reference;
                [compDistInfoList, rxnDistInfoList] = InputKegg._getKeggCompRxnDistInfo(dirSourceData = dirSourceData,
                                                                                        resource_xlsxFile = resource_xlsxFile,
                                                                                        )
                #
                # Get info of every reaction's coef/coefSum/compID/canoSmiles/nonIsoSmiles and isoReaProReactionIDs from
                # the cleaned KEGG reaction excle(InfoEx).
                # rxnIsoReaProIdSet: Reaction IDs set of which the Reactants or Products have Isomeric Compounds.
                [rxnCoefCompSmiDict, rxnIsoReaProIdSet] = InputKegg._getKeggRxnInfo(dirSourceData = dirSourceData,
                                                                                    resource_xlsxFile = resource_xlsxFile,
                                                                                    )
                with open(path.join("temp", "KEGG_Info.data"), "wb") as file_obj:
                    pickle.dump([rxnCoefCompSmiDict, rxnIsoReaProIdSet], file_obj)

                # The direction of every reaction.
                rxnDirInfoDict = InputKegg._getKeggPsRxnDirInfo()
                #
                # Get the compound pairs of every reaction.
                rxnCompDrPairsInfoDict = InputKegg._getKeggPsRxnCompDrPairInfo(dirSourceData = dirSourceData,
                                                                               resource_xlsxFile = resource_xlsxFile,
                                                                               dirDatas = dirDatas,
                                                                               general_cofactors_xlsxFile = general_cofactors_xlsxFile,
                                                                               dirRxnFile = dirRxnFile,
                                                                               diffThresholdLevel = diffThresholdLevel,
                                                                               )
                #
                # Fusing....
                rxnFuseInfoList = InputKegg._fuseKeggPsRxnCompPairInfo(rxnDistInfoList, rxnDirInfoDict, rxnCompDrPairsInfoDict,
                                                                       diffThresholdLevel = diffThresholdLevel,
                                                                       )
                #
                # Get and fuse the local total in/out-degree of every compound.
                compFuseInfoList = InputKegg._gfLtiod(compDistInfoList, rxnFuseInfoList, diffThresholdLevel)
                #
                # """
                with open(path.join("temp", "KEGG_Debug.data"), "wb") as file_obj:
                    pickle.dump([compFuseInfoList, rxnFuseInfoList], file_obj)
                # """
                #
                # '''
                #
                with open(path.join("temp", "KEGG_Debug.data"), "rb") as file_obj:
                    [compFuseInfoList, rxnFuseInfoList] = pickle.load(file_obj)
                #
                # Saving...
                InputKegg._saveKeggPsRxnCompPairInfo(compFuseInfoList, rxnFuseInfoList,
                                                     dirSourceData = dirSourceData,
                                                     xlsxFile = xlsxFile,
                                                     diffThresholdLevel = diffThresholdLevel,
                                                     )
                #
                # Generate input prompt file(text)
                InputKegg._generateInputPrompt(compFuseInfoList, rxnFuseInfoList)
                #
                # Generate '.data' file to speed up data loading.
                InputKegg._generateDataFile(compFuseInfoList, rxnFuseInfoList,
                                            diffThresholdLevel = diffThresholdLevel,
                                            )
                #
            else:
                pass
        elif dataBase == "Rhea":
            pass
        elif dataBase == "MetaCyc":
            if isUpdate or (xlsxFile not in listdir(dirSourceData)):
                #
                # '''
                # Distill the information from MetaCyc.
                # The fields of the information of compounds and reactions is given as follows.
                # Compounds:
                #  0      1      2                 3                      4                          5       6
                # Entry; Names; Chemical Formula; Delta Gibbs(kcal/mol); Molecular Weight(Daltons); SMILES; Reactions;
                #  7       8     9      10    11        12        13       14
                # Source; MetaCyc; ChEBI; KEGG; BiGG; MetaNetX; PubChem;  Reference;
                # Reactions:
                #  0      1      2           3         4            5         6       7       8                      9
                # Entry; Names; Definition; Equation; Coefficient; Compound; SMILES; Status; Delta Gibbs(kcal/mol); EC Number;
                #  10      11    12    13    14        15       16       17
                # Source; MetaCyc; Rhea; KEGG; BiGG; MetaNetX; UniProt; Reference;
                #
                [compDistInfoList, rxnDistInfoList] = InputMetaCyc._getMetaCycCompRxnDistInfo(dirSourceData = dirSourceData,
                                                                                              resource_xlsxFile = resource_xlsxFile,
                                                                                              )
                #
                # Get info of every reaction's coef/coefSum/compID/canoSmiles/nonIsoSmiles and isoReaProReactionIDs from
                # the cleaned MetaCyc(KndPad) reaction excle(InfoEx).
                # rxnIsoReaProIdSet: Reaction IDs set of which the Reactants or Products have Isomeric Compounds.
                [rxnCoefCompSmiDict, rxnIsoReaProIdSet] = InputMetaCyc._getMetaCycRxnInfo(dirSourceData = dirSourceData,
                                                                                          resource_xlsxFile = resource_xlsxFile,
                                                                                          )
                with open(path.join("temp", "MetaCyc_Info.data"), "wb") as file_obj:
                    pickle.dump([rxnCoefCompSmiDict, rxnIsoReaProIdSet], file_obj)
                #
                # The direction of every reaction.
                rxnDirInfoDict = InputMetaCyc._getMetaCycPsRxnDirInfo()
                #
                # Get the compound pairs of every reaction.
                rxnCompDrPairsInfoDict = InputMetaCyc._getMetaCycPsRxnCompDrPairInfo(dirSourceData = dirSourceData,
                                                                                     resource_xlsxFile = resource_xlsxFile,
                                                                                     dirDatas = dirDatas,
                                                                                     general_cofactors_xlsxFile = general_cofactors_xlsxFile,
                                                                                     dirRxnFile = dirRxnFile,
                                                                                     diffThresholdLevel = diffThresholdLevel,
                                                                                     )
                #
                # Fusing....
                rxnFuseInfoList = InputMetaCyc._fuseMetaCycPsRxnCompPairInfo(rxnDistInfoList, rxnDirInfoDict, rxnCompDrPairsInfoDict,
                                                                             diffThresholdLevel = diffThresholdLevel,
                                                                             )
                #
                # Get and fuse the local total in/out-degree of every compound.
                compFuseInfoList = InputMetaCyc._gfLtiod(compDistInfoList, rxnFuseInfoList, diffThresholdLevel)
                #
                #
                with open(path.join("temp", "MetaCyc_Debug.data"), "wb") as file_obj:
                    pickle.dump([compFuseInfoList, rxnFuseInfoList], file_obj)
                #
                # '''
                #
                with open(path.join("temp", "MetaCyc_Debug.data"), "rb") as file_obj:
                    [compFuseInfoList, rxnFuseInfoList] = pickle.load(file_obj)
                #
                # Saving...
                InputMetaCyc._saveMetaCycPsRxnCompPairInfo(compFuseInfoList, rxnFuseInfoList,
                                                           dirSourceData = dirSourceData,
                                                           xlsxFile = xlsxFile,
                                                           diffThresholdLevel = diffThresholdLevel,
                                                           )
                #
                # Generate input prompt file(text)
                InputMetaCyc._generateInputPrompt(compFuseInfoList, rxnFuseInfoList)
                #
                # Generate '.data' file to speed up data loading.
                InputMetaCyc._generateDataFile(compFuseInfoList, rxnFuseInfoList,
                                               diffThresholdLevel = diffThresholdLevel,
                                               )
                #
            else:
                pass
        elif dataBase == "KndPad":
            if isUpdate or (xlsxFile not in listdir(dirSourceData)):
                #
                # '''
                # Distill the information from KndPad.
                # The fields of the information of compounds and reactions is given as follows.
                # Compounds:
                #  0      1      2                 3                      4                          5       6
                # Entry; Names; Chemical Formula; Delta Gibbs(kcal/mol); Molecular Weight(Daltons); SMILES; Reactions;
                #  7       8     9      10    11        12        13       14
                # Source; MetaCyc; ChEBI; KEGG; BiGG; MetaNetX; PubChem;  Reference;
                # Reactions:
                #  0      1      2           3         4            5         6       7       8                      9
                # Entry; Names; Definition; Equation; Coefficient; Compound; SMILES; Status; Delta Gibbs(kcal/mol); EC Number;
                #  10      11    12    13    14        15       16       17
                # Source; MetaCyc; Rhea; KEGG; BiGG; MetaNetX; UniProt; Reference;
                #
                [compDistInfoList, rxnDistInfoList] = InputKndPad._getKndPadCompRxnDistInfo(dirSourceData = dirSourceData,
                                                                                            resource_xlsxFile = resource_xlsxFile,
                                                                                            )
                # # return None
                # Get info of every reaction's coef/coefSum/compID/canoSmiles/nonIsoSmiles and isoReaProReactionIDs from
                # the cleaned KndPad reaction excle(InfoEx).
                # rxnIsoReaProIdSet: Reaction IDs set of which the Reactants or Products have Isomeric Compounds.
                [rxnCoefCompSmiDict, rxnIsoReaProIdSet] = InputKndPad._getKndPadRxnInfo(dirSourceData = dirSourceData,
                                                                                        resource_xlsxFile = resource_xlsxFile,
                                                                                        )
                with open(path.join("temp", "KndPad_Info.data"), "wb") as file_obj:
                    pickle.dump([rxnCoefCompSmiDict, rxnIsoReaProIdSet], file_obj)
                #
                # The direction of every reaction.
                rxnDirInfoDict = InputKndPad._getKndPadPsRxnDirInfo()
                #
                # Get the compound pairs of every reaction.
                rxnCompDrPairsInfoDict = InputKndPad._getKndPadPsRxnCompDrPairInfo(dirSourceData = dirSourceData,
                                                                                   resource_xlsxFile = resource_xlsxFile,
                                                                                   dirDatas = dirDatas,
                                                                                   general_cofactors_xlsxFile = general_cofactors_xlsxFile,
                                                                                   dirRxnFile = dirRxnFile,
                                                                                   diffThresholdLevel = diffThresholdLevel,
                                                                                   )
                #
                # Fusing....
                rxnFuseInfoList = InputKndPad._fuseKndPadPsRxnCompPairInfo(rxnDistInfoList, rxnDirInfoDict, rxnCompDrPairsInfoDict,
                                                                           diffThresholdLevel = diffThresholdLevel,
                                                                           )
                #
                # Get and fuse the local total in/out-degree of every compound.
                compFuseInfoList = InputKndPad._gfLtiod(compDistInfoList, rxnFuseInfoList, diffThresholdLevel)
                #
                #
                with open(path.join("temp", "KndPad_Debug.data"), "wb") as file_obj:
                    pickle.dump([compFuseInfoList, rxnFuseInfoList], file_obj)
                #
                # '''
                #
                with open(path.join("temp", "KndPad_Debug.data"), "rb") as file_obj:
                    [compFuseInfoList, rxnFuseInfoList] = pickle.load(file_obj)
                #
                # Saving...
                InputKndPad._saveKndPadPsRxnCompPairInfo(compFuseInfoList, rxnFuseInfoList,
                                                        dirSourceData = dirSourceData,
                                                        xlsxFile = xlsxFile,
                                                        diffThresholdLevel = diffThresholdLevel,
                                                        )
                #
                # Generate input prompt file(text)
                InputKndPad._generateInputPrompt(compFuseInfoList, rxnFuseInfoList)
                #
                # Generate '.data' file to speed up data loading.
                InputKndPad._generateDataFile(compFuseInfoList, rxnFuseInfoList,
                                              diffThresholdLevel = diffThresholdLevel,
                                              )
                #
            else:
                pass

        else:
            pass


    @classmethod
    def loadRxnCompPair(cls, metabolic_network_file = "../sourcedata/KEGG_Pathway_Search_Ori.xlsx",
                             diffThreshold = 0.1,
                             dataBase = "KEGG",
                             isUpdate = True,
                             diffThresholdLevel = 10,
                             ):
        """
            Load compound pairs of every reaction from KEGG, Rhea, MetaCyc or KndPad.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """

        [dirSourceData, xlsxFile] = path.split(metabolic_network_file)

        #
        if dataBase == "KEGG":
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = InputKegg.loadRxnCompPair(dirSourceData = dirSourceData,
                                                                            xlsxFile = xlsxFile,
                                                                            diffThreshold = diffThreshold,
                                                                            isUpdate = isUpdate,
                                                                            diffThresholdLevel = diffThresholdLevel,
                                                                            )
            return [rxnDistInfoDict, compDistInfoDict]
        elif dataBase == "Rhea":
            pass
        elif dataBase == "MetaCyc":
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = InputMetaCyc.loadRxnCompPair(dirSourceData = dirSourceData,
                                                                               xlsxFile = xlsxFile,
                                                                               diffThreshold = diffThreshold,
                                                                               isUpdate = isUpdate,
                                                                               diffThresholdLevel = diffThresholdLevel,
                                                                               )
            return [rxnDistInfoDict, compDistInfoDict]
        elif dataBase == "KndPad":
            # Loading...
            # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Status; Coefficient;
            # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
            [rxnDistInfoDict, compDistInfoDict] = InputKndPad.loadRxnCompPair(dirSourceData = dirSourceData,
                                                                              xlsxFile = xlsxFile,
                                                                              diffThreshold = diffThreshold,
                                                                              isUpdate = isUpdate,
                                                                              diffThresholdLevel = diffThresholdLevel,
                                                                              )
            return [rxnDistInfoDict, compDistInfoDict]
        else:
            pass
        #


    @classmethod
    def getOrgCompRxnId(cls, org_model_location = "../sourcedata/models/",
                             organism = "eco",
                             json_model_file_dict = dict(),
                             dataBase = "KEGG",
                             ):
        """
            Get the information of the selected organism from bigg.
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        # "bsu": "bsu_iYO844.json",
        # "eco": "eco_iJO1366.json",
        # "kpn": "kpn_iYL1228.json",
        # "ppu": "ppu_iJN1462.json",
        # "sce": "sce_iMM904.json",
        # "syz": "syz_iSynCJ816.json",
        # ...
        #
        if dataBase == "KEGG":
            orgCompRxnId_dict = InputKegg._getOrgCompRxnId(org_model_location = org_model_location,
                                                           organism = organism,
                                                           json_model_file_dict = json_model_file_dict,
                                                           )
            return orgCompRxnId_dict
        elif dataBase == "Rhea":
            pass
        elif dataBase == "MetaCyc":
            orgCompRxnId_dict = InputMetaCyc._getOrgCompRxnId(org_model_location = org_model_location,
                                                              organism = organism,
                                                              json_model_file_dict = json_model_file_dict,
                                                              )
            return orgCompRxnId_dict
        elif dataBase == "KndPad":
            orgCompRxnId_dict = InputKndPad._getOrgCompRxnId(org_model_location = org_model_location,
                                                             organism = organism,
                                                             json_model_file_dict = json_model_file_dict,
                                                             )
            return orgCompRxnId_dict
        else:
            pass
        #


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





