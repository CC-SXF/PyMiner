# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 10:15:31 2021

@author: CC-SXF
"""

from rdkit import Chem
from rdkit.Chem import AllChem


class FingerprintInfo():
    """
    """
    def __init__(self):
        pass


    @classmethod
    def getRDKitFingerprintInfo(cls, mol):
        """
        """
        # http://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint
        #
        # A dict with bits as keys and corresponding bond paths as values.
        rdkitBitInfo = dict()
        Chem.RDKFingerprint(mol, minPath=1, maxPath=7, fpSize=2048, nBitsPerHash=2, useHs=True,
                            tgtDensity=0.0, minSize=128, branchedPaths=True, useBondOrder=True,
                            atomInvariants=[], fromAtoms=[], atomBits=[], bitInfo=rdkitBitInfo)
        rdkitBit = sorted(rdkitBitInfo.keys())
        return rdkitBit


    @classmethod
    def getMorganFingerprintInfo(cls, mol):
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

        morganBit = sorted(morganBitInfo.keys())
        return morganBit


    @classmethod
    def demofunc(cls,):
        """ """
        pass


if __name__ == '__main__':
    """ """
    pass





