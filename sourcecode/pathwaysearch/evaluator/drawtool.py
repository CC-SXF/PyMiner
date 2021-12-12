# -*- coding: utf-8 -*-
"""
Created on Thu May 27 21:46:34 2021

@author: CC-SXF
"""

from os import path
from copy import deepcopy

from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem.rdchem import KekulizeException
from rdkit.Chem import rdDepictor
# from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


RCenterColors = [
                    ['white', (1.0, 1.0, 1.0)],
                    ['palegreen', (0.596078431372549, 0.984313725490196, 0.596078431372549)],
                    ['khaki', (0.9411764705882353, 0.9019607843137255, 0.5490196078431373)],
                    ['paleturquoise', (0.6862745098039216, 0.9333333333333333, 0.9333333333333333)],
                    ['bisque', (1.0, 0.8941176470588236, 0.7686274509803922)],
                    ['lightblue', (0.6784313725490196, 0.8470588235294118, 0.9019607843137255)],
                    ['mistyrose', (1.0, 0.8941176470588236, 0.8823529411764706)],
                    ['pink', (1.0, 0.7529411764705882, 0.796078431372549)],
                    ['lightskyblue', (0.5294117647058824, 0.807843137254902, 0.9803921568627451)],
                    ['greenyellow', (0.6784313725490196, 1.0, 0.1843137254901961)],
                    ['plum', (0.8666666666666667, 0.6274509803921569, 0.8666666666666667)],
                    ['lightsalmon', (1.0, 0.6274509803921569, 0.47843137254901963)],
                    ['yellow', (1.0, 1.0, 0.0)],
                    ['darkseagreen', (0.5607843137254902, 0.7372549019607844, 0.5607843137254902)],
                    ['orchid', (0.8549019607843137, 0.4392156862745098, 0.8392156862745098)],
                    ['rosybrown', (0.7372549019607844, 0.5607843137254902, 0.5607843137254902)],
                    ['lime', (0.0, 1.0, 0.0)],
                    ['lightseagreen', (0.12549019607843137, 0.6980392156862745, 0.6666666666666666)],
                    ['peru', (0.803921568627451, 0.5215686274509804, 0.24705882352941178)],
                    ['green', (0.0, 0.5019607843137255, 0.0)],
                    ['indianred', (0.803921568627451, 0.3607843137254902, 0.3607843137254902)],
                    ['cyan', (0.0, 1.0, 1.0)],
                    ['dodgerblue', (0.11764705882352941, 0.5647058823529412, 1.0)],
                    ['orange', (1.0, 0.6470588235294118, 0.0)],
                    ['royalblue', (0.2549019607843137, 0.4117647058823529, 0.8823529411764706)]
                ]

from collections import namedtuple
MoleculeRCenterEnv = namedtuple('MoleculeRCenterEnv', ['mol', 'highlightAtoms', 'atomColors', 'highlightRadii', 'highlightBonds', 'bondColors', 'legend'])




class DrawTool():
    """
        Drawing Tools.
    """
    def __init__(self, ):
        """ """
        pass


    @classmethod
    def _getMoleculeRCenterEnv(cls, mol, molId = "",
                                    rCenterAtomMapNumList = list(), rCenterBondMapNumPairList = list(),
                                    rCenterAtomComponentLabelList = list(), rCenterBondComponentLabelList = list(),
                                    rCenterAtomColor = (0.0, 1.0, 0.0), rCenterBondColor = (0.0, 1.0, 0.0),
                                    baseRad = 0.6, diffColor = False,
                                    ):
        """
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        # Return the number of conformations on the molecule.
        if not mol.GetNumConformers():
            # Compute 2D coordinates for a molecule.
            # The resulting coordinates are stored on each atom of the molecule.
            rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True)
        legend = molId
        numAtoms = mol.GetNumAtoms()
        #
        # Convert atom-atom mapping number to atom index.
        atomMapNumIdxDict = dict()
        for atomIdx in range(numAtoms):
            # Atom indices start at 0
             atomMapNum = mol.GetAtomWithIdx(atomIdx).GetAtomMapNum()
             atomMapNumIdxDict[atomMapNum] = atomIdx
        #
        # Get the indices of highlight atoms.
        # Color all atoms in the reaction center of the compound.
        highlightAtoms = list()
        atomColors = dict()
        highlightRadii = dict()
        for i, atomMapNum in enumerate(rCenterAtomMapNumList):
            try:
                atomIdx = atomMapNumIdxDict[atomMapNum]
            except KeyError:
                continue
            highlightAtoms.append(atomIdx)
            if diffColor:
                atomComponentLabel = rCenterAtomComponentLabelList[i]
                atomColors[atomIdx] = RCenterColors[atomComponentLabel][1]
            else:
                atomColors[atomIdx] = rCenterAtomColor
            highlightRadii[atomIdx] = baseRad
        #
        # Convert highlight bonds' atom-atom mapping number pair to atom index pair.
        highlightBondAtomIdxPairList = list()
        for bondMapNumPair in rCenterBondMapNumPairList:
            bondAtomIdxPair = (atomMapNumIdxDict[bondMapNumPair[0]], atomMapNumIdxDict[bondMapNumPair[1]])
            highlightBondAtomIdxPairList.append(bondAtomIdxPair)
        # Get the indices of highlight bonds.
        # Color all bonds in the reaction center of the compound.
        highlightBonds = list()
        bondcolors = dict()
        for i, bondAtomIdxPair in enumerate(highlightBondAtomIdxPairList):
            bondIdx = mol.GetBondBetweenAtoms(bondAtomIdxPair[0], bondAtomIdxPair[1]).GetIdx()
            highlightBonds.append(bondIdx)
            if diffColor:
                # bondcolors[bondIdx] = (1.0, 1.0, 1.0)
                bondComponentLabel = rCenterBondComponentLabelList[i]
                bondcolors[bondIdx] = RCenterColors[bondComponentLabel][1]
            else:
                # bondcolors[bondIdx] = (1.0, 1.0, 1.0)
                bondcolors[bondIdx] = rCenterBondColor
        #
        # Does some cleanup operations on the molecule to prepare it to draw nicely.
        # The operations include: kekulization, addition of chiral Hs (so that we can draw wedges to them),
        # wedging of bonds at chiral centers, and
        # generation of a 2D conformation if the molecule does not already have a conformation.
        # Returns a modified copy of the molecule.
        """
        try:
            mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
        except KekulizeException:
            mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        """
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        #
        # Remove atom-atom mapping numbers of all atoms in this compound.
        for atomIdx in range(numAtoms):
            # Atom indices start at 0
             atom = mol.GetAtomWithIdx(atomIdx)
             atom.SetAtomMapNum(0)
        #
        return MoleculeRCenterEnv(mol, highlightAtoms, atomColors, highlightRadii,
                                  highlightBonds, bondcolors, legend)
        #


    @classmethod
    def drawMoleculeRCenter(cls, mol, molId = "",
                                 rCenterAtomMapNumList = list(), rCenterBondMapNumPairList = list(),
                                 rCenterAtomComponentLabelList = list(), rCenterBondComponentLabelList = list(),
                                 rCenterAtomColor = (0.6, 1.0, 0.6), rCenterBondColor = (0.6, 1.0, 0.6), molSize = (600, 600),
                                 baseRad = 0.6, useSVG = False, drawOptions = None, diffColor = False, highlightHs = False ,
                                 sub_route = "",
                                 ):
        """
            ARGUMENTS:
                ...

            RETURNS:
                ...
        """
        #
        """
        if len(rCenterAtomComponentLabelList) != 0:
            if len(rCenterAtomMapNumList) != len(rCenterAtomComponentLabelList):
                print("")
        else:
            diffColor = False
        #
        if len(rCenterBondComponentLabelList) != 0:
            if len(rCenterBondMapNumPairList) != len(rCenterBondComponentLabelList):
                print("")
        else:
            diffColor = False
        """
        #
        # Hide Hydrogen Atoms For Non-Sub Molecule
        if not highlightHs:
            # Avoid modifying the object.
            mol = deepcopy(mol)
            numAtoms = mol.GetNumAtoms()
            for atomIdx in range(numAtoms):
                atom = mol.GetAtomWithIdx(atomIdx)
                atomMapNum = atom.GetAtomMapNum()
                if atomMapNum in rCenterAtomMapNumList:
                    atomSymbol = atom.GetSymbol()
                    if atomSymbol == "H":
                        atom.SetIsotope(1)
            smiles= Chem.MolToSmiles(mol)
            mol = Chem.MolFromSmiles(smiles)
            numAtoms = mol.GetNumAtoms()
            for atomIdx in range(numAtoms):
                atom = mol.GetAtomWithIdx(atomIdx)
                atomMapNum = atom.GetAtomMapNum()
                if atomMapNum in rCenterAtomMapNumList:
                    atomSymbol = atom.GetSymbol()
                    if atomSymbol == "H":
                        atom.SetIsotope(0)
        #
        rcEnv = cls._getMoleculeRCenterEnv(mol, molId, rCenterAtomMapNumList, rCenterBondMapNumPairList,
                                           rCenterAtomComponentLabelList, rCenterBondComponentLabelList,
                                           rCenterAtomColor, rCenterBondColor, baseRad, diffColor)
        #
        if useSVG:
            drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
        else:
            drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])
        if drawOptions is None:
            drawopt = drawer.drawOptions()
            drawopt.continuousHighlight = True
        else:
            drawOptions.continuousHighlight = True
            drawer.SetDrawOptions(drawOptions)
        drawer.DrawMolecule(rcEnv.mol, highlightAtoms=rcEnv.highlightAtoms,
                            highlightAtomColors=rcEnv.atomColors, highlightAtomRadii=rcEnv.highlightRadii,
                            highlightBonds=rcEnv.highlightBonds, highlightBondColors=rcEnv.bondColors,
                            confId=-1, legend="")
        drawer.FinishDrawing()
        #
        #
        if useSVG:
            image_file = "".join([molId, ".svg"])
            image_file = path.join(sub_route, image_file)
            svgStr = drawer.GetDrawingText()
            with open(image_file, "w") as file_obj:
                file_obj.write(svgStr)
        else:
            image_file = "".join([molId, ".png"])
            image_file = path.join(sub_route, image_file)
            drawer.WriteDrawingText(image_file)
        #


    @classmethod
    def demoFunc(cls, ):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





