# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:45:25 2021

@author: CC-SXF
"""

import json
from openpyxl import load_workbook


class Reconcile():
    """
    """

    def __init__(self, ):
        """ """
        pass

    @classmethod
    def _getSheetInfo(cls, xlsxFile, sheetName):
        """
        """
        workbook = load_workbook(xlsxFile)
        workbook.guess_types = True
        sheet = workbook[sheetName]
        sheetAllRows = list()
        for row in sheet.rows:
            sheetOneRow = list()
            for cell in row:
                sheetOneRow.append(cell.value)
            sheetAllRows.append(sheetOneRow)
        return sheetAllRows


    @classmethod
    def _analyseKegg(cls, xlsxFile, sheetName="Reaction"):
        """
        """
        rxnInfoList = cls._getSheetInfo(xlsxFile, sheetName)
        #
        for rxnInfo in rxnInfoList[1:]: # Exclude title
            rxnId = rxnInfo[0]
            rxnComp = json.loads(rxnInfo[6])
            rxnCompPair = json.loads(rxnInfo[8])
            reaComp = rxnComp[0]
            proComp = rxnComp[1]
            rpCompPair = rxnCompPair[0]
            prCompPair = rxnCompPair[1]
            #
            for key, value in rpCompPair.items():
                if (key in reaComp) and ( (set(value) & set(proComp)) == set(value)):
                    pass
                else:
                    print("Please check on the compound pairs of reaction:", rxnId)
            for key, value in prCompPair.items():
                if (key in proComp) and ((set(value) & set(reaComp)) == set(value)):
                    pass
                else:
                    print("Please check on the compound pairs of reaction:", rxnId)
        return True


    @classmethod
    def reconcileKegg(cls, xlsxFile="", sheetName=""):
        """
        """
        if xlsxFile == "":
            xlsxFile = "datas/KEGG_Pathway_Search_Manual_Annotation.xlsx"
        if sheetName == "":
            sheetName = "Reaction"
        try:
            rxnInfoList = cls._getSheetInfo(xlsxFile, sheetName)
        except FileNotFoundError:
            rxnInfoList = list()
        #
        rxnDirCompPairDict = dict()
        for rxnInfo in rxnInfoList[1:]: # Exclude title
            rxnId = rxnInfo[0]
            rxnDir =rxnInfo[4]
            rxnCompPair = rxnInfo[8]
            rxnDirCompPairDict[rxnId] = [rxnDir, rxnCompPair]
        #
        return rxnDirCompPairDict


    @classmethod
    def reconcileMetaCyc(cls, xlsxFile="", sheetName=""):
        """
        """
        rxnDirCompPairDict = dict()
        return rxnDirCompPairDict


    @classmethod
    def reconcileKndPad(cls, xlsxFile="", sheetName=""):
        """
        """
        rxnDirCompPairDict = dict()
        return rxnDirCompPairDict


    @classmethod
    def demoFunc(cls, ):
        """ """
        pass



if __name__ == "__main__":
    """ """
    xlsxFile = "../datas/KEGG_Pathway_Search_Manual_Annotation.xlsx"
    Reconcile._analyseKegg(xlsxFile)
    rxnDirCompPairDict = Reconcile.reconcileKegg(xlsxFile)


