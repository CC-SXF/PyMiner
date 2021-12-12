# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 17:02:08 2021

@author: CC-SXF
"""

from rdkit import Chem
# # from pathwaysearch.fingerprint.fingerprintinfo import FingerprintInfo
from pathwaysearch.fingerprint.subfingerprintinfo import SubFingerprintInfo


class FingerPrint():
    """
    """

    def __init__(self, ):
        """ """
        pass

    @classmethod
    def _addMapNum(cls, mol):
        """ """
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            mapNum = idx + 1
            atom.SetAtomMapNum(mapNum)
        return mol

    @classmethod
    def _removeMapNum(cls, mol):
        """ """
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    @classmethod
    def _getExSubStrInfo(cls, compSmiles, pattSmarts):
        """
            Get the information of excluded substructure of one compound.
        """
        compMol = Chem.MolFromSmiles(compSmiles)
        compMol = Chem.AddHs(compMol)
        compMol = cls._addMapNum(compMol)
        #
        # charge
        mapnum_charge_dict = dict()
        for idx in range(compMol.GetNumAtoms()):
            atom = compMol.GetAtomWithIdx(idx)
            atom_mapnum = atom.GetAtomMapNum()
            atom_charge = atom.GetFormalCharge()
            mapnum_charge_dict[atom_mapnum] = atom_charge
        #
        patt = Chem.MolFromSmarts(pattSmarts)
        bondsToUse = SubFingerprintInfo._getBondsToUse(compMol, patt)
        subMol = Chem.PathToSubmol(compMol, bondsToUse)
        #
        # charge
        for idx in range(subMol.GetNumAtoms()):
            atom = subMol.GetAtomWithIdx(idx)
            atom_mapnum = atom.GetAtomMapNum()
            atom_charge = mapnum_charge_dict[atom_mapnum]
            atom.SetFormalCharge(atom_charge)
        #
        subMol = cls._removeMapNum(subMol)
        rdkitSubBit = SubFingerprintInfo.getRDKitSubFingerprintInfo(compMol, bondsToUse)
        morganSubBit = SubFingerprintInfo.getMorganSubFingerprintInfo(compMol, bondsToUse)
        return (subMol, rdkitSubBit, morganSubBit)

    @classmethod
    def getExSubStrInfo(cls, dataBase = "KEGG"):
        """
            Get the information of excluded substructures.
        """
        if dataBase == "KEGG":

            exSubStrTupleList = list()
            #
            # CoA-SH
            coa_sh = 'CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCS'
            pattSmarts = '[S]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(coa_sh, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # PPi
            ppi = 'O=P(O)(O)OP(=O)(O)O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(ppi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # Pi
            pi = 'O=P(O)(O)O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(pi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            return exSubStrTupleList
        elif dataBase == "Rhea":
            pass
        elif (dataBase == "MetaCyc"):
            exSubStrTupleList = list()
            #
            # CoA-SH
            coa_sh = 'CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS'
            pattSmarts = '[S]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(coa_sh, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # PPi
            ppi = 'O=P([O-])([O-])OP(=O)([O-])O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(ppi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # Pi
            pi = 'O=P([O-])([O-])O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(pi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            return exSubStrTupleList
        elif (dataBase == "KndPad"):
            exSubStrTupleList = list()
            #
            # CoA-SH(KEGG)
            coa_sh = 'CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCS'
            pattSmarts = '[S]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(coa_sh, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            # CoA-SH(MetaCyc)
            coa_sh = 'CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS'
            pattSmarts = '[S]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(coa_sh, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # PPi(KEGG)
            ppi = 'O=P(O)(O)OP(=O)(O)O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(ppi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            # PPi(MetaCyc)
            ppi = 'O=P([O-])([O-])OP(=O)([O-])O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(ppi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            # Pi(KEGG)
            pi = 'O=P(O)(O)O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(pi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            # Pi(MetaCyc)
            pi = 'O=P([O-])([O-])O'
            pattSmarts = '[P]-[O]-[H]'
            exSubStrTuple = cls._getExSubStrInfo(pi, pattSmarts)
            exSubStrTupleList.append(exSubStrTuple)
            #
            return exSubStrTupleList
        else:
            pass

    @classmethod
    def demoFunc(cls, ):
        """ """
        pass


if __name__ == '__main__':
    """ """
    #
    # KEGG
    #
    # fructose 1,6-biphosphate
    mol_fbp_kegg = Chem.MolFromSmiles('O=P(O)(O)OC[C@H]1OC(O)(COP(=O)(O)O)[C@@H](O)[C@@H]1O')
    mol_fbp_kegg = Chem.AddHs(mol_fbp_kegg)
    patt = FingerPrint.getExSubStrInfo('KEGG')[2][0]
    atomIdx_tuple_tuple_fbp_kegg = mol_fbp_kegg.GetSubstructMatches(patt)
    print("FBP, KEGG:", atomIdx_tuple_tuple_fbp_kegg)
    #
    # isopentenyl diphosphate
    mol_ipp_kegg = Chem.MolFromSmiles('C=C(C)CCOP(=O)(O)OP(=O)(O)O')
    mol_ipp_kegg = Chem.AddHs(mol_ipp_kegg)
    patt = FingerPrint.getExSubStrInfo('KEGG')[1][0]
    atomIdx_tuple_tuple_ipp_kegg = mol_ipp_kegg.GetSubstructMatches(patt)
    print("IPP, KEGG:", atomIdx_tuple_tuple_ipp_kegg)
    # MetaCyc
    #
    # fructose 1,6-biphosphate
    mol_fbp_metacyc = Chem.MolFromSmiles('O=P([O-])([O-])OC[C@H]1O[C@](O)(COP(=O)([O-])[O-])[C@@H](O)[C@@H]1O')
    mol_fbp_metacyc = Chem.AddHs(mol_fbp_metacyc)
    patt = FingerPrint.getExSubStrInfo('KndPad')[-3][0]
    atomIdx_tuple_tuple_fbp_metacyc = mol_fbp_metacyc.GetSubstructMatches(patt)
    print("FBP, MetaCyc:", atomIdx_tuple_tuple_fbp_metacyc)
    #
    # isopentenyl diphosphate
    mol_ipp_metacyc = Chem.MolFromSmiles('C=C(C)CCOP(=O)([O-])OP(=O)([O-])[O-]')
    mol_ipp_metacyc = Chem.AddHs(mol_ipp_metacyc)
    patt = FingerPrint.getExSubStrInfo('KndPad')[-1][0]
    atomIdx_tuple_tuple_ipp_metacyc = mol_ipp_metacyc.GetSubstructMatches(patt)
    print("IPP, MetaCyc:", atomIdx_tuple_tuple_ipp_metacyc)
    #





