# -*- coding: utf-8 -*-
"""
Created on Tue May 18 23:00:12 2021

@author: CC-SXF
"""

from PyQt5 import QtCore, QtGui, QtWidgets


class AddTableItem():
    """
    """
    def __init__(self):
        """ """
        pass

    @classmethod
    def addPlainItem(cls, table, idx_row, idx_column, ic_details=""):
        """
        """
        cell_label = QtWidgets.QLabel()
        cell_label.setOpenExternalLinks(False)
        cell_label.setAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        cell_label.setText(ic_details)
        table.setCellWidget(idx_row, idx_column, cell_label)


    @classmethod
    def addColorItem(cls, table, idx_row, idx_column, last_s_details="", current_s_details=""):
        """
        """
        last_s_details_list = last_s_details.split('-->')
        current_s_details_list = current_s_details.split('-->')
        #
        if last_s_details == "":
            new_s_details_list = current_s_details_list
        else:
            for idx_diff, rxn2comp_last  in enumerate(last_s_details_list):
                rxn2comp_current = current_s_details_list[idx_diff]
                if rxn2comp_last != rxn2comp_current:
                    break
            rxn2comp_diff = current_s_details_list[idx_diff]
            rxn2comp_diff = "".join(["<font color='fuchsia'>", rxn2comp_diff, "</font>"])
            current_s_details_list[idx_diff] = rxn2comp_diff
            new_s_details_list =current_s_details_list
        #
        new_s_details = '-->'.join(new_s_details_list)
        cell_label = QtWidgets.QLabel()
        cell_label.setOpenExternalLinks(False)
        cell_label.setAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        cell_label.setText(new_s_details)
        table.setCellWidget(idx_row, idx_column, cell_label)


    @classmethod
    def addLinkItem(cls, table, idx_row, idx_column, link_content, link_postfix , database, datatype):
        """
        """
        database_dict = {"KEGG/Compound": "https://www.kegg.jp/entry/",
                         "KEGG/Reaction": "https://www.kegg.jp/entry/",
                         "ChEBI/Compound": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=",
                         "Rhea/Reaction": "https://www.rhea-db.org/rhea/",
                         "MetaCyc/Compound": "https://metacyc.org/compound?orgid=META&id=",
                         "MetaCyc/Reaction": "https://metacyc.org/META/NEW-IMAGE?type=REACTION&object=",
                         "UniProtKB/Enzyme": "https://uniprot.org/uniprot/?query=ec:",
                         }
        database_key = "/".join([database, datatype])
        link_prefix = database_dict[database_key]
        link_external = "".join([link_prefix, link_postfix])
        # text-decoration:line-through; text-decoration: overline; text-decoration:underline; text-decoration:none; text-decoration:blink;
        link_full = "".join(['''<a style='color: blue; text-decoration: underline;' href="''', link_external, '">', link_content, '</a>'])
        cell_label = QtWidgets.QLabel()
        cell_label.setOpenExternalLinks(True)
        cell_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        cell_label.setText(link_full)
        table.setCellWidget(idx_row, idx_column, cell_label)


    @classmethod
    def demoFunc(cls,):
        """ """
        pass





