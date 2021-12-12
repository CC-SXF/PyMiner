# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 10:46:46 2021

@author: CC-SXF
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from pathwaysearch.fingerprint.fingerprintinfo import FingerprintInfo


class SubFingerprintInfo():
    """
    """
    def __init__(self):
        pass


    @classmethod
    def getRDKitSubFingerprintInfo(cls, mol, subBondsToUse):
        """
        """
        # http://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint
        #
        # A dict with bits as keys and corresponding bond paths as values.
        rdkitBitInfo = dict()
        Chem.RDKFingerprint(mol, minPath=1, maxPath=7, fpSize=2048, nBitsPerHash=2, useHs=True,
                            tgtDensity=0.0, minSize=128, branchedPaths=True, useBondOrder=True,
                            atomInvariants=[], fromAtoms=[], atomBits=[], bitInfo=rdkitBitInfo)
        subBondsToUse = set(subBondsToUse)
        #
        rdkitSubBit = set()
        for bit in rdkitBitInfo.keys():
            bondpathList = rdkitBitInfo[bit]
            for bondpath in bondpathList:
                bondpath = set(bondpath)
                if (bondpath & subBondsToUse) == bondpath:
                    rdkitSubBit.add(bit)
                    break
                else:
                    pass
        rdkitSubBit = sorted(rdkitSubBit)
        # rdkitSubBit = json.dumps(rdkitSubBit)
        return rdkitSubBit


    @classmethod
    def getMorganSubFingerprintInfo(cls, mol, subBondsToUse):
        """
        """
        # http://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html#rdkit.Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect
        # Information is available about the atoms that contribute to particular bits in the Morgan fingerprint via the bitInfo argument.
        # The dictionary provided is populated with one entry per bit set in the fingerprint, the keys are the bit ids, the values are lists of (atom index, radius) tuples.
        #
        # A dict with bits as keys and corresponding tuples(atom index, radius) as values.
        morganBitInfo = dict()
        AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, invariants=[], fromAtoms=[],
                                              useChirality=False, useBondTypes=True, useFeatures=False,
                                              bitInfo=morganBitInfo, includeRedundantEnvironments=False)

        subBondsToUse = set(subBondsToUse)
        #
        subAtomsToUse = set()
        for bond_idx in subBondsToUse:
            bond = mol.GetBondWithIdx(bond_idx)
            subAtomsToUse.add(bond.GetBeginAtomIdx())
            subAtomsToUse.add(bond.GetEndAtomIdx())
        subAtomsToUse = sorted(subAtomsToUse)
        #
        atomIdxRadiusTupleSet = set()
        for atom_idx in subAtomsToUse:
            atomIdxRadiusTupleSet.add((atom_idx, 0))
            #
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=1, rootedAtAtom=atom_idx, useHs=True)
            env = set(env)
            if (env & subBondsToUse) == env:
                atomIdxRadiusTupleSet.add((atom_idx, 1))
            else:
                continue
            #
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=2, rootedAtAtom=atom_idx, useHs=True)
            env = set(env)
            if (env & subBondsToUse) == env:
                atomIdxRadiusTupleSet.add((atom_idx, 2))
            else:
                pass
        #
        morganSubBit = set()
        for bit in morganBitInfo.keys():
            atomIdxRadiusTupleTuple = morganBitInfo[bit]
            for atomIdxRadiusTuple in atomIdxRadiusTupleTuple:
                if atomIdxRadiusTuple in atomIdxRadiusTupleSet:
                    morganSubBit.add(bit)
                    break
        morganSubBit = sorted(morganSubBit)
        # morganSubBit = json.dumps(morganSubBit)
        return morganSubBit


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
    def _removeSubMapNum(cls, mol, subMol):
        """ """
        atomIdx_tuple_tuple = mol.GetSubstructMatches(subMol, useChirality=False)
        for atomIdx_tuple in atomIdx_tuple_tuple:
            for atomIdx in atomIdx_tuple:
                atom = mol.GetAtomWithIdx(atomIdx)
                atom.SetAtomMapNum(0)
        return mol


    @classmethod
    def _getBondsToUse(cls, mol, patt):
        """ """
        atom_idx_tuple = mol.GetSubstructMatch(patt, useChirality=False)
        #
        excluded_bond_set = set()
        for idx_1st in atom_idx_tuple:
            for idx_2nd in atom_idx_tuple:
                try:
                    bond_idx = mol.GetBondBetweenAtoms(idx_1st, idx_2nd).GetIdx()
                    excluded_bond_set.add(bond_idx)
                except AttributeError:
                    pass
        #
        bondsToUse = list()
        for bond_idx in range(mol.GetNumBonds()):
            if bond_idx not in excluded_bond_set:
                bondsToUse.append(bond_idx)
        #
        return bondsToUse


    @classmethod
    def demofunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    # '''
    # KEGG
    # CoA-SH
    coa_sh = 'CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCS'
    coa_sh = Chem.MolFromSmiles(coa_sh)
    coa_sh = Chem.AddHs(coa_sh)
    coa_sh = SubFingerprintInfo._addMapNum(coa_sh)
    print("CoA-SH(KEGG):")
    print("\t", Chem.MolToSmiles(coa_sh))
    # exAtomMapNumPairList = [(48, 84), ]  # *S-H
    patt = Chem.MolFromSmarts('[S]-[H]')
    bondsToUse = SubFingerprintInfo._getBondsToUse(coa_sh, patt)
    subMol_k = Chem.PathToSubmol(coa_sh, bondsToUse)
    subMol_k = SubFingerprintInfo._removeMapNum(subMol_k)
    coa_rdkitSubBit_k = SubFingerprintInfo.getRDKitSubFingerprintInfo(coa_sh, bondsToUse)
    coa_morganSubBit_k = SubFingerprintInfo.getMorganSubFingerprintInfo(coa_sh, bondsToUse)
    # Acetyl-CoA
    acetyl_coa = 'CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O'
    acetyl_coa = Chem.MolFromSmiles(acetyl_coa)
    acetyl_coa = Chem.AddHs(acetyl_coa)
    acetyl_coa = SubFingerprintInfo._addMapNum(acetyl_coa)
    acetyl_coa = SubFingerprintInfo._removeSubMapNum(acetyl_coa, subMol_k)
    ace_rdkitBit = FingerprintInfo.getRDKitFingerprintInfo(acetyl_coa)
    ace_morganBit = FingerprintInfo.getMorganFingerprintInfo(acetyl_coa)
    print((set(ace_rdkitBit) & set(coa_rdkitSubBit_k))  == set(coa_rdkitSubBit_k))
    print((set(ace_morganBit) & set(coa_morganSubBit_k)) == set(coa_morganSubBit_k))
    print(acetyl_coa.HasSubstructMatch(subMol_k, useChirality=False))
    print(acetyl_coa.GetSubstructMatches(subMol_k, useChirality=False))
    print()
    # '''
    #
    # '''
    # MetaCyc
    coa_sh = 'CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS'
    coa_sh = Chem.MolFromSmiles(coa_sh)
    coa_sh = Chem.AddHs(coa_sh)
    coa_sh = SubFingerprintInfo._addMapNum(coa_sh)
    print("CoA-SH(MetaCyc):")
    print("\t", Chem.MolToSmiles(coa_sh))
    # exAtomMapNumPairList = [(48, 80), ]  # *S-H
    patt = Chem.MolFromSmarts('[S]-[H]')
    bondsToUse = SubFingerprintInfo._getBondsToUse(coa_sh, patt)
    subMol_m = Chem.PathToSubmol(coa_sh, bondsToUse)
    subMol_m = SubFingerprintInfo._removeMapNum(subMol_m)
    coa_rdkitSubBit_m = SubFingerprintInfo.getRDKitSubFingerprintInfo(coa_sh, bondsToUse)
    coa_morganSubBit_m = SubFingerprintInfo.getMorganSubFingerprintInfo(coa_sh, bondsToUse)
    # Acetyl-CoA
    acetyl_coa = 'CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-]'
    acetyl_coa = Chem.MolFromSmiles(acetyl_coa)
    acetyl_coa = Chem.AddHs(acetyl_coa)
    acetyl_coa = SubFingerprintInfo._addMapNum(acetyl_coa)
    acetyl_coa = SubFingerprintInfo._removeSubMapNum(acetyl_coa, subMol_m)
    ace_rdkitBit = FingerprintInfo.getRDKitFingerprintInfo(acetyl_coa)
    ace_morganBit = FingerprintInfo.getMorganFingerprintInfo(acetyl_coa)
    print((set(ace_rdkitBit) & set(coa_rdkitSubBit_m))  == set(coa_rdkitSubBit_m))
    print((set(ace_morganBit) & set(coa_morganSubBit_m)) == set(coa_morganSubBit_m))
    print(acetyl_coa.HasSubstructMatch(subMol_m, useChirality=False))
    print(acetyl_coa.GetSubstructMatches(subMol_m, useChirality=False))
    print()
    # '''
    #
    # '''
    # KEGG & MetaCyc
    print("CoA-SH(KEGG/MetaCyc)")
    print("RDKit Fingerprint:", (set(coa_rdkitSubBit_k) & set(coa_rdkitSubBit_m)) == set(coa_rdkitSubBit_m))
    print("Morgan Fingerprint:", (set(coa_morganSubBit_k) & set(coa_morganSubBit_m)) == set(coa_morganSubBit_m))
    # '''





