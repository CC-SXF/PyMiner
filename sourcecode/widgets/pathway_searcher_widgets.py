# -*- coding: utf-8 -*-
"""
Created on Sat May 15 09:57:45 2021

@author: CC-SXF
"""

import re
import json
from os import path
# import shutil
from PyQt5 import QtCore, QtGui, QtWidgets

from rdkit import Chem
import rdkit.Chem.Draw as Chem_Draw

from widgets.new_widgets import NewEventLineEdit
from widgets.new_widgets import MultCompleter
from widgets.new_widgets import NewAppendTextBrowser
from widgets.new_widgets import NewEventTableWidget
from widgets.add_table_item import AddTableItem
from widgets.search_thread import SearchThread


class Pathway_Seracher_Widgets():
    """ """
    def __init__(self):
        """ """
        self.initSearchThread()
        self.initInputPrompt()
        self.initOrgModel()
        self.initWidgets()
        #


    def initOrgModel(self):
        """ """
        #
        # "" : ""
        # bsu: Bacillus subtilis subsp. subtilis 168
        # eco: Escherichia coli K-12 MG1655
        # kpn: Klebsiella pneumoniae subsp. pneumoniae MGH 78578
        # ppu: Pseudomonas putida KT2440
        # sce: Saccharomyces cerevisiae S288c
        # syz: Synechocystis sp. PCC 6803(Cyanobacteria)
        # ...
        #
        with open(path.join("datas", "models.json"), "r") as file_obj:
            models = json.load(file_obj)
        self.organism_list = sorted(models.keys())
        #


    def _initInputPrompt(self):
        """ """
        #
        inputprompt_dict = dict()
        inputprompt_dict["KEGG"] = dict()
        inputprompt_dict["MetaCyc"] = dict()
        inputprompt_dict["KndPad"] = dict()
        #
        try:
            with open(path.join('temp', 'KEGG_InputPrompt.json'), "r") as file_obj:
                inputprompt_kegg_dict = json.load(file_obj)
            inputprompt_dict['KEGG']['Compound'] = inputprompt_kegg_dict['Compound']
            inputprompt_dict['KEGG']['Reaction'] = inputprompt_kegg_dict['Reaction']
        except FileNotFoundError:
            inputprompt_dict['KEGG']['Compound'] = list()
            inputprompt_dict['KEGG']['Reaction'] = list()
        #
        try:
            with open(path.join('temp', 'MetaCyc_InputPrompt.json'), "r") as file_obj:
                inputprompt_metacyc_dict = json.load(file_obj)
            inputprompt_dict['MetaCyc']['Compound'] = inputprompt_metacyc_dict['Compound']
            inputprompt_dict['MetaCyc']['Reaction'] = inputprompt_metacyc_dict['Reaction']
        except FileNotFoundError:
            inputprompt_dict['MetaCyc']['Compound'] = list()
            inputprompt_dict['MetaCyc']['Reaction'] = list()
        #
        try:
            with open(path.join('temp', 'KndPad_InputPrompt.json'), "r") as file_obj:
                inputprompt_kndpad_dict = json.load(file_obj)
            inputprompt_dict['KndPad']['Compound'] = inputprompt_kndpad_dict['Compound']
            inputprompt_dict['KndPad']['Reaction'] = inputprompt_kndpad_dict['Reaction']
        except FileNotFoundError:
            inputprompt_dict['KndPad']['Compound'] = list()
            inputprompt_dict['KndPad']['Reaction'] = list()
        #
        return inputprompt_dict


    def initInputPrompt(self):
        """ """
        #
        self.inputprompt_dict = self._initInputPrompt()
        #
        self.comp_inputprompt_list = self.inputprompt_dict['KEGG']['Compound']
        self.rxn_inputprompt_list = self.inputprompt_dict['KEGG']['Reaction']
        #
        # source completer
        self.kegg_source_completer = MultCompleter(self.inputprompt_dict['KEGG']['Compound'])
        self.kegg_source_completer.setMaxVisibleItems(10)
        self.kegg_source_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kegg_source_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kegg_source_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # target completer
        self.kegg_target_completer = QtWidgets.QCompleter(self.inputprompt_dict['KEGG']['Compound'])
        self.kegg_target_completer.setMaxVisibleItems(10)
        self.kegg_target_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kegg_target_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kegg_target_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid compounds completer
        self.kegg_avoid_compounds_completer = MultCompleter(self.inputprompt_dict['KEGG']['Compound'])
        self.kegg_avoid_compounds_completer.setMaxVisibleItems(10)
        self.kegg_avoid_compounds_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kegg_avoid_compounds_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kegg_avoid_compounds_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid reactions completer
        self.kegg_avoid_reactions_completer = MultCompleter(self.inputprompt_dict['KEGG']['Reaction'])
        self.kegg_avoid_reactions_completer.setMaxVisibleItems(10)
        self.kegg_avoid_reactions_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kegg_avoid_reactions_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kegg_avoid_reactions_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        #
        # source completer
        self.metacyc_source_completer = MultCompleter(self.inputprompt_dict['MetaCyc']['Compound'])
        self.metacyc_source_completer.setMaxVisibleItems(10)
        self.metacyc_source_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.metacyc_source_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.metacyc_source_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # target completer
        self.metacyc_target_completer = QtWidgets.QCompleter(self.inputprompt_dict['MetaCyc']['Compound'])
        self.metacyc_target_completer.setMaxVisibleItems(10)
        self.metacyc_target_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.metacyc_target_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.metacyc_target_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid compounds completer
        self.metacyc_avoid_compounds_completer = MultCompleter(self.inputprompt_dict['MetaCyc']['Compound'])
        self.metacyc_avoid_compounds_completer.setMaxVisibleItems(10)
        self.metacyc_avoid_compounds_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.metacyc_avoid_compounds_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.metacyc_avoid_compounds_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid reactions completer
        self.metacyc_avoid_reactions_completer = MultCompleter(self.inputprompt_dict['MetaCyc']['Reaction'])
        self.metacyc_avoid_reactions_completer.setMaxVisibleItems(10)
        self.metacyc_avoid_reactions_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.metacyc_avoid_reactions_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.metacyc_avoid_reactions_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        #
        # source completer
        self.kndpad_source_completer = MultCompleter(self.inputprompt_dict['KndPad']['Compound'])
        self.kndpad_source_completer.setMaxVisibleItems(10)
        self.kndpad_source_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kndpad_source_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kndpad_source_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # target completer
        self.kndpad_target_completer = QtWidgets.QCompleter(self.inputprompt_dict['KndPad']['Compound'])
        self.kndpad_target_completer.setMaxVisibleItems(10)
        self.kndpad_target_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kndpad_target_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kndpad_target_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid compounds completer
        self.kndpad_avoid_compounds_completer = MultCompleter(self.inputprompt_dict['KndPad']['Compound'])
        self.kndpad_avoid_compounds_completer.setMaxVisibleItems(10)
        self.kndpad_avoid_compounds_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kndpad_avoid_compounds_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kndpad_avoid_compounds_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #
        # avoid reactions completer
        self.kndpad_avoid_reactions_completer = MultCompleter(self.inputprompt_dict['KndPad']['Reaction'])
        self.kndpad_avoid_reactions_completer.setMaxVisibleItems(10)
        self.kndpad_avoid_reactions_completer.setFilterMode(QtCore.Qt.MatchStartsWith)
        self.kndpad_avoid_reactions_completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.kndpad_avoid_reactions_completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        #


    def initWidgets(self):
        """ """
        #
        self.win_ps_top = QtWidgets.QWidget()
        self.win_ps_top.setObjectName('win_ps_top')
        self.win_ps_middle = QtWidgets.QWidget()
        self.win_ps_middle.setObjectName('win_ps_middle')
        self.win_ps_down = QtWidgets.QWidget()
        self.win_ps_down.setObjectName('win_ps_down')
        #
        self.win_ps_trisplitter =QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.win_ps_trisplitter.setObjectName('win_ps_trisplitter')
        self.win_ps_trisplitter.setHandleWidth(2)   # The width of the handle of the splitter
        self.win_ps_trisplitter.addWidget(self.win_ps_top)
        self.win_ps_trisplitter.addWidget(self.win_ps_middle)
        self.win_ps_trisplitter.addWidget(self.win_ps_down)
        self.win_ps_trisplitter.setSizes([300, 200, 500])
        #
        self.win_ps_hboxlayout = QtWidgets.QHBoxLayout()
        self.win_ps_hboxlayout.setObjectName('win_ps_hboxlayout')
        self.win_ps_hboxlayout.addWidget(self.win_ps_trisplitter)


        self.subwin_ps_top_input = QtWidgets.QWidget()
        self.subwin_ps_top_input.setObjectName("Input")
        screen_width = QtWidgets.QDesktopWidget().screenGeometry().width()
        self.h_zoom = round(screen_width/1920, 1)
        self.set_ps_top_input_widgets(h_zoom=self.h_zoom)
        #
        self.subwin_ps_top_tabwidget = QtWidgets.QTabWidget()
        self.subwin_ps_top_tabwidget.setObjectName('subwin_ps_top_tabwidget')
        self.subwin_ps_top_tabwidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.subwin_ps_top_tabwidget.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.subwin_ps_top_tabwidget.setElideMode(QtCore.Qt.ElideRight)
        self.subwin_ps_top_tabwidget.setCurrentIndex(0)
        self.subwin_ps_top_tabwidget.addTab(self.subwin_ps_top_input, QtGui.QIcon('images/Input.ico'), "Inputs")
        #
        self.win_ps_top_gridlayout = QtWidgets.QGridLayout(self.win_ps_top)
        self.win_ps_top_gridlayout.setObjectName('win_ps_top_gridlayout')
        self.win_ps_top_gridlayout.addWidget(self.subwin_ps_top_tabwidget, 0, 0)


        self.subwin_ps_middle_table = QtWidgets.QWidget()
        self.subwin_ps_middle_table.setObjectName("Pathways")
        # # screen_width = QtWidgets.QDesktopWidget().screenGeometry().width()
        # # self.h_zoom = round(screen_width/1920, 1)
        self.set_ps_middle_table_widgets(h_zoom=self.h_zoom)
        #
        self.subwin_ps_middle_tabwidget = QtWidgets.QTabWidget()
        self.subwin_ps_middle_tabwidget.setObjectName('subwin_ps_middle_tabwidget')
        self.subwin_ps_middle_tabwidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.subwin_ps_middle_tabwidget.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.subwin_ps_middle_tabwidget.setElideMode(QtCore.Qt.ElideRight)
        self.subwin_ps_middle_tabwidget.setCurrentIndex(0)
        self.subwin_ps_middle_tabwidget.addTab(self.subwin_ps_middle_table, QtGui.QIcon('images/Table.ico'), "Pathways")
        #
        self.win_ps_middle_gridlayout = QtWidgets.QGridLayout(self.win_ps_middle)
        self.win_ps_middle_gridlayout.setObjectName('win_ps_middle_gridlayout')
        self.win_ps_middle_gridlayout.addWidget(self.subwin_ps_middle_tabwidget, 0, 0)


        self.subwin_ps_down_text = QtWidgets.QWidget()
        self.subwin_ps_down_text.setObjectName("Tips")
        # # screen_width = QtWidgets.QDesktopWidget().screenGeometry().width()
        # # self.h_zoom = round(screen_width/1920, 1)
        self.set_ps_down_text_widgets(h_zoom=self.h_zoom)
        #
        self.subwin_ps_down_image = QtWidgets.QWidget()
        self.subwin_ps_down_image.setObjectName("Details")
        # # screen_width = QtWidgets.QDesktopWidget().screenGeometry().width()
        # # self.h_zoom = round(screen_width/1920, 1)
        self.set_ps_down_image_widgets(h_zoom=self.h_zoom)
        #
        self.subwin_ps_down_tabwidget = QtWidgets.QTabWidget()
        self.subwin_ps_down_tabwidget.setObjectName('subwin_ps_down_tabwidget')
        self.subwin_ps_down_tabwidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.subwin_ps_down_tabwidget.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.subwin_ps_down_tabwidget.setElideMode(QtCore.Qt.ElideRight)
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        self.subwin_ps_down_tabwidget.addTab(self.subwin_ps_down_text, QtGui.QIcon('images/Tips.ico'), "Tips")
        self.subwin_ps_down_tabwidget.addTab(self.subwin_ps_down_image, QtGui.QIcon('images/Details.ico'), "Details")
        #
        self.win_ps_down_gridlayout = QtWidgets.QGridLayout(self.win_ps_down)
        self.win_ps_down_gridlayout.setObjectName('win_ps_down_gridlayout')
        self.win_ps_down_gridlayout.addWidget(self.subwin_ps_down_tabwidget, 0, 0)


    def set_ps_top_input_widgets(self, h_zoom=1.0):
        """ """
        #
        # source
        self.source_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.source_label.setGeometry(QtCore.QRect(20*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.source_label.setObjectName("source_label")
        self.source_label.setText("Sources:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.source_label.setFont(font)
        # source
        self.source_lineEdit = NewEventLineEdit(self.subwin_ps_top_input)
        self.source_lineEdit.setGeometry(QtCore.QRect(20*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.source_lineEdit.setObjectName("source_lineEdit")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.source_lineEdit.setFont(font)
        # self.source_lineEdit.setText('D-Glyceraldehyde 3-phosphate(C00118); Pyruvate(C00022); ')  # 1
        # self.source_lineEdit.setText('L-Tryptophan(C00078); ')  # 2
        self.source_lineEdit.setText('L-Phenylalanine(C00079); L-Tyrosine(C00082); ')  # 3
        # source completer
        self.source_input_ori = ""
        self.source_lineEdit.setCompleter(self.kegg_source_completer)
        self.source_lineEdit.newKeyPressEvent.connect(self.source_update_input)

        # target
        self.target_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.target_label.setGeometry(QtCore.QRect(290*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.target_label.setObjectName("target_label")
        self.target_label.setText("Target:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.target_label.setFont(font)
        # target
        self.target_lineEdit = NewEventLineEdit(self.subwin_ps_top_input)
        self.target_lineEdit.setGeometry(QtCore.QRect(290*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.target_lineEdit.setObjectName("target_lineEdit")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.target_lineEdit.setFont(font)
        # self.target_lineEdit.setText('Taxa-4,11-diene(C11894)')  # 1
        # self.target_lineEdit.setText('Violacein(C21136)')  # 2
        self.target_lineEdit.setText('Resveratrol(C03582)')  # 3
        # target completer
        # self.target_input_ori = ""
        self.target_lineEdit.setCompleter(self.kegg_target_completer)
        self.target_lineEdit.newKeyPressEvent.connect(self.target_update_input)

        # avoid compounds
        self.avoid_compounds_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.avoid_compounds_label.setGeometry(QtCore.QRect(560*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.avoid_compounds_label.setObjectName("avoid_compounds_label")
        self.avoid_compounds_label.setText("Avoid Compounds:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.avoid_compounds_label.setFont(font)
        # avoid compounds
        self.avoid_compounds_lineEdit = NewEventLineEdit(self.subwin_ps_top_input)
        self.avoid_compounds_lineEdit.setGeometry(QtCore.QRect(560*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.avoid_compounds_lineEdit.setObjectName("avoid_compounds_lineEdit")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.avoid_compounds_lineEdit.setFont(font)
        self.avoid_compounds_lineEdit.setText("")
        # self.avoid_compounds_lineEdit.setText('Acetyl-CoA(C00024); ')
        # avoid compounds completer
        self.avoid_compounds_input_ori = ""
        self.avoid_compounds_lineEdit.setCompleter(self.kegg_avoid_compounds_completer)
        self.avoid_compounds_lineEdit.newKeyPressEvent.connect(self.avoid_compounds_update_input)

        # avoid reactions
        self.avoid_reactions_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.avoid_reactions_label.setGeometry(QtCore.QRect(830*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.avoid_reactions_label.setObjectName("avoid_reactions_label")
        self.avoid_reactions_label.setText("Avoid Reactions:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.avoid_reactions_label.setFont(font)
        # avoid reactions
        self.avoid_reactions_lineEdit = NewEventLineEdit(self.subwin_ps_top_input)
        self.avoid_reactions_lineEdit.setGeometry(QtCore.QRect(830*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.avoid_reactions_lineEdit.setObjectName("avoid_reactions_lineEdit")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.avoid_reactions_lineEdit.setFont(font)
        # # self.avoid_reactions_lineEdit.setText('R02950; ')
        # avoid reactions completer
        self.avoid_reactions_input_ori = ""
        self.avoid_reactions_lineEdit.setCompleter(self.kegg_avoid_reactions_completer)
        self.avoid_reactions_lineEdit.newKeyPressEvent.connect(self.avoid_reactions_update_input)


        # host organism
        self.host_organism_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.host_organism_label.setGeometry(QtCore.QRect(20*h_zoom, 85*h_zoom, 250*h_zoom, 20*h_zoom))
        self.host_organism_label.setObjectName("host_organism_label")
        self.host_organism_label.setText( "Host Organism:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.host_organism_label.setFont(font)
        # host organism
        self.host_organism_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.host_organism_comboBox.setGeometry(QtCore.QRect(20*h_zoom, 110*h_zoom, 250*h_zoom, 22*h_zoom))
        self.host_organism_comboBox.setObjectName("host_organism_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.host_organism_comboBox.setFont(font)
        self.host_organism_comboBox.setEditable(False)
        # # self.organism_list = ['', 'bsu', 'eco', 'kpn', 'ppu', 'sce', 'syz']
        self.host_organism_comboBox.addItems(self.organism_list)
        self.host_organism_comboBox.setCurrentIndex(2)  # Escherichia coli
        # "" : ""
        # bsu: Bacillus subtilis subsp. subtilis 168
        # eco: Escherichia coli K-12 MG1655
        # kpn: Klebsiella pneumoniae subsp. pneumoniae MGH 78578
        # ppu: Pseudomonas putida KT2440
        # sce: Saccharomyces cerevisiae S288c
        # syz: Synechocystis sp. PCC 6803(Cyanobacteria)
        #

        # database
        self.database_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.database_label.setGeometry(QtCore.QRect(290*h_zoom, 85*h_zoom, 250*h_zoom, 20*h_zoom))
        self.database_label.setObjectName("database_label")
        self.database_label.setText("Database:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.database_label.setFont(font)
        # datase
        self.database_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.database_comboBox.setGeometry(QtCore.QRect(290*h_zoom, 110*h_zoom, 250*h_zoom, 22*h_zoom))
        self.database_comboBox.setObjectName("database_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.database_comboBox.setFont(font)
        self.database_comboBox.setEditable(False)
        self.dataBase_list = ["KEGG", "MetaCyc", "KndPad", ]
        self.database_comboBox.addItems(self.dataBase_list)
        self.database_comboBox.setCurrentIndex(0)
        self.database_comboBox.activated.connect(self.selectDatabase)

        # maximum length
        self.maximum_length_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.maximum_length_label.setGeometry(QtCore.QRect(560*h_zoom, 85*h_zoom, 250*h_zoom, 20*h_zoom))
        self.maximum_length_label.setObjectName("maximum_length_label")
        self.maximum_length_label.setText("Maximum Length:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.maximum_length_label.setFont(font)
        # maximum length
        self.maximum_length_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.maximum_length_comboBox.setGeometry(QtCore.QRect(560*h_zoom, 110*h_zoom, 250*h_zoom, 22*h_zoom))
        self.maximum_length_comboBox.setObjectName("maximum_length_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.maximum_length_comboBox.setFont(font)
        self.maximum_length_comboBox.setEditable(False)
        self.maximum_length_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']
        self.maximum_length_comboBox.addItems(self.maximum_length_list)
        self.maximum_length_comboBox.setCurrentIndex(3)

        # maximum time(s)
        self.maximum_time_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.maximum_time_label.setGeometry(QtCore.QRect(830*h_zoom, 85*h_zoom, 250*h_zoom, 20*h_zoom))
        self.maximum_time_label.setObjectName("maximum_time_label")
        self.maximum_time_label.setText("Maximum Time(s):")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.maximum_time_label.setFont(font)
        # maximum time(s)
        self.maximum_time_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.maximum_time_comboBox.setGeometry(QtCore.QRect(830*h_zoom, 110*h_zoom, 250*h_zoom, 22*h_zoom))
        self.maximum_time_comboBox.setObjectName("maximum_time_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.maximum_time_comboBox.setFont(font)
        self.maximum_time_comboBox.setEditable(False)
        self.maximum_time_list = [str(i) for i in range(1, 6001)]
        self.maximum_time_comboBox.addItems(self.maximum_time_list)
        self.maximum_time_comboBox.setCurrentIndex(119)


        # similarity difference threshold
        self.diffThreshold_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.diffThreshold_label.setGeometry(QtCore.QRect(20*h_zoom, 150*h_zoom, 250*h_zoom, 20*h_zoom))
        self.diffThreshold_label.setObjectName("diffThreshold_label")
        self.diffThreshold_label.setText("Similarity Difference Threshold:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.diffThreshold_label.setFont(font)
        # similarity difference threshold
        self.diffThreshold_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.diffThreshold_comboBox.setGeometry(QtCore.QRect(20*h_zoom, 175*h_zoom, 250*h_zoom, 22*h_zoom))
        self.diffThreshold_comboBox.setObjectName("diffThreshold_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.diffThreshold_comboBox.setFont(font)
        self.diffThreshold_comboBox.setEditable(False)
        self.diffThreshold_list = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
        self.diffThreshold_comboBox.addItems(self.diffThreshold_list)
        self.diffThreshold_comboBox.setCurrentIndex(0)

        # methods
        self.methods_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.methods_label.setGeometry(QtCore.QRect(290*h_zoom, 150*h_zoom, 130*h_zoom, 20*h_zoom))
        self.methods_label.setObjectName("methods_label")
        self.methods_label.setText("Search Method:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.methods_label.setFont(font)
        # methods
        self.methods_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.methods_comboBox.setGeometry(QtCore.QRect(290*h_zoom, 175*h_zoom, 130*h_zoom, 22*h_zoom))
        self.methods_comboBox.setObjectName("methods_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.methods_comboBox.setFont(font)
        self.methods_comboBox.setEditable(False)
        # # self.methods_list = ['bfs_wl', "bfs_naive", "dfs_naive"]
        # # self.methods_list = ['bfs_wl', "bfs_naive"]
        self.methods_list = ["bfs_naive", "dfs_naive"]
        self.methods_comboBox.addItems(self.methods_list)
        self.methods_comboBox.setCurrentIndex(0)


        # isTotal
        self.isTotal_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isTotal_checkBox.setGeometry(QtCore.QRect(440*h_zoom, 150*h_zoom, 100*h_zoom, 20*h_zoom))
        self.isTotal_checkBox.setObjectName("isTotal_checkBox")
        self.isTotal_checkBox.setText("Total")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isTotal_checkBox.setFont(font)
        self.isTotal_checkBox.setTristate(False)
        self.isTotal_checkBox.setChecked(True)
        self.isTotal_checkBox.setEnabled(True)

        # isShortest
        self.isShortest_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isShortest_checkBox.setGeometry(QtCore.QRect(440*h_zoom, 175*h_zoom, 100*h_zoom, 22*h_zoom))
        self.isShortest_checkBox.setObjectName("isShortest_checkBox")
        self.isShortest_checkBox.setText("Shortest")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isShortest_checkBox.setFont(font)
        self.isShortest_checkBox.setTristate(False)
        self.isShortest_checkBox.setChecked(False)
        self.isShortest_checkBox.setEnabled(True)
        # signal and plot
        self.isShortest_checkBox.stateChanged.connect(self.reset_methods)


        # isRetro
        self.isRetro_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isRetro_checkBox.setGeometry(QtCore.QRect(560*h_zoom, 150*h_zoom, 120*h_zoom, 20*h_zoom))
        self.isRetro_checkBox.setObjectName("isRetro_checkBox")
        self.isRetro_checkBox.setText("Retro Search")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isRetro_checkBox.setFont(font)
        self.isRetro_checkBox.setTristate(False)
        self.isRetro_checkBox.setChecked(True)
        self.isRetro_checkBox.setEnabled(False)

        # isSmart
        self.isSmart_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isSmart_checkBox.setGeometry(QtCore.QRect(560*h_zoom, 175*h_zoom, 120*h_zoom, 22*h_zoom))
        self.isSmart_checkBox.setObjectName("isSmart_checkBox")
        self.isSmart_checkBox.setText("Smart Search")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isSmart_checkBox.setFont(font)
        self.isSmart_checkBox.setTristate(False)
        self.isSmart_checkBox.setChecked(True)
        self.isSmart_checkBox.setEnabled(True)
        # signal and plot
        self.isSmart_checkBox.stateChanged.connect(self.isValid_isRetro)


        # isInfeasible
        self.isInfeasible_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isInfeasible_checkBox.setGeometry(QtCore.QRect(690*h_zoom, 150*h_zoom, 120*h_zoom, 20*h_zoom))
        self.isInfeasible_checkBox.setObjectName("isInfeasible_checkBox")
        self.isInfeasible_checkBox.setText("Infeasibility")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isInfeasible_checkBox.setFont(font)
        self.isInfeasible_checkBox.setTristate(False)
        self.isInfeasible_checkBox.setChecked(False)
        self.isInfeasible_checkBox.setEnabled(True)

        # isTransfer
        self.isTransfer_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isTransfer_checkBox.setGeometry(QtCore.QRect(690*h_zoom, 175*h_zoom, 120*h_zoom, 22*h_zoom))
        self.isTransfer_checkBox.setObjectName("isTransfer_checkBox")
        self.isTransfer_checkBox.setText("Atom Trace")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isTransfer_checkBox.setFont(font)
        self.isTransfer_checkBox.setTristate(False)
        self.isTransfer_checkBox.setChecked(True)
        self.isTransfer_checkBox.setEnabled(True)


        # condition(aerobic/anaerobic)
        self.isAerobic_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isAerobic_checkBox.setGeometry(QtCore.QRect(830*h_zoom, 150*h_zoom, 130*h_zoom, 20*h_zoom))
        self.isAerobic_checkBox.setObjectName("isAerobic_checkBox")
        self.isAerobic_checkBox.setText("Aerobic Culture")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isAerobic_checkBox.setFont(font)
        self.isAerobic_checkBox.setTristate(False)
        self.isAerobic_checkBox.setChecked(True)
        self.isAerobic_checkBox.setEnabled(True)

        # isFlux
        self.isFlux_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isFlux_checkBox.setGeometry(QtCore.QRect(830*h_zoom, 175*h_zoom, 130*h_zoom, 22*h_zoom))
        self.isFlux_checkBox.setObjectName("isFlux_checkBox")
        self.isFlux_checkBox.setText("Flux")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isFlux_checkBox.setFont(font)
        self.isFlux_checkBox.setTristate(False)
        self.isFlux_checkBox.setChecked(True)
        self.isFlux_checkBox.setEnabled(True)
        # signal and plot
        self.isFlux_checkBox.stateChanged.connect(self.isValid_isAerobic)


        # isUpdate
        self.isUpdate_checkBox = QtWidgets.QCheckBox(self.subwin_ps_top_input)
        self.isUpdate_checkBox.setGeometry(QtCore.QRect(980*h_zoom, 150*h_zoom, 100*h_zoom, 20*h_zoom))
        self.isUpdate_checkBox.setObjectName("isUpdate_checkBox")
        self.isUpdate_checkBox.setText("Update")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.isUpdate_checkBox.setFont(font)
        self.isUpdate_checkBox.setTristate(False)
        self.isUpdate_checkBox.setChecked(False)
        self.isUpdate_checkBox.setEnabled(True)
        # signal and plot
        self.isUpdate_checkBox.stateChanged.connect(self.isPermit_update)


        # resource file
        self.resource_file_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.resource_file_label.setGeometry(QtCore.QRect(1100*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.resource_file_label.setObjectName("resource_file_label")
        self.resource_file_label.setText("Resource File:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.resource_file_label.setFont(font)
        self.resource_file_label.setEnabled(False)
        # resource file
        self.resource_file_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.resource_file_comboBox.setGeometry(QtCore.QRect(1100*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.resource_file_comboBox.setObjectName("resource_file_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.resource_file_comboBox.setFont(font)
        self.resource_file_comboBox.setEnabled(False)
        self.resource_file_comboBox.setEditable(False)
        self.resource_file_list = ["../sourcedata/KEGG_Ori_Cleaned_Half.xlsx"]
        self.resource_file_comboBox.addItems(self.resource_file_list)
        self.resource_file_comboBox.setCurrentIndex(0)

        # general cofactors file
        self.general_cofactors_file_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.general_cofactors_file_label.setGeometry(QtCore.QRect(1100*h_zoom, 85*h_zoom, 250*h_zoom, 20*h_zoom))
        self.general_cofactors_file_label.setObjectName("general_cofactors_file_label")
        self.general_cofactors_file_label.setText("General Cofactors File:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.general_cofactors_file_label.setFont(font)
        self.general_cofactors_file_label.setEnabled(False)
        # general cofactors file
        self.general_cofactors_file_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.general_cofactors_file_comboBox.setGeometry(QtCore.QRect(1100*h_zoom, 110*h_zoom, 250*h_zoom, 22*h_zoom))
        self.general_cofactors_file_comboBox.setObjectName("general_cofactors_file_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.general_cofactors_file_comboBox.setFont(font)
        self.general_cofactors_file_comboBox.setEnabled(False)
        self.general_cofactors_file_comboBox.setEditable(False)
        self.general_cofactors_file_list = ["datas/KEGG_General_Cofactors.xlsx", "datas/KEGG_General_Cofactors_Empty.xlsx"]
        self.general_cofactors_file_comboBox.addItems(self.general_cofactors_file_list)
        self.general_cofactors_file_comboBox.setCurrentIndex(0)
        self.general_cofactors_file_comboBox.activated.connect(self.cofactor_update_input)

        # metabolic network file
        self.metabolic_network_file_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.metabolic_network_file_label.setGeometry(QtCore.QRect(1100*h_zoom, 150*h_zoom, 250*h_zoom, 20*h_zoom))
        self.metabolic_network_file_label.setObjectName("metabolic_network_file_label")
        self.metabolic_network_file_label.setText( "Distilled Metabolic Network File:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.metabolic_network_file_label.setFont(font)
        self.metabolic_network_file_label.setEnabled(False)
        # metabolic network file
        self.metabolic_network_file_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.metabolic_network_file_comboBox.setGeometry(QtCore.QRect(1100*h_zoom, 175*h_zoom, 250*h_zoom, 22*h_zoom))
        self.metabolic_network_file_comboBox.setObjectName("metabolic_network_file_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.metabolic_network_file_comboBox.setFont(font)
        self.metabolic_network_file_comboBox.setEnabled(False)
        self.metabolic_network_file_comboBox.setEditable(False)
        self.metabolic_network_file_list = ["../sourcedata/KEGG_Pathway_Search_Ori.xlsx"]
        self.metabolic_network_file_comboBox.addItems(self.metabolic_network_file_list)
        self.metabolic_network_file_comboBox.setCurrentIndex(0)

        # atom mapping file location
        self.atom_mapping_file_location_label = QtWidgets.QLabel(self.subwin_ps_top_input)
        self.atom_mapping_file_location_label.setGeometry(QtCore.QRect(1370*h_zoom, 20*h_zoom, 250*h_zoom, 20*h_zoom))
        self.atom_mapping_file_location_label.setObjectName("atom_mapping_file_location_label")
        self.atom_mapping_file_location_label.setText("Atom Mapping Files:")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setWeight(75)
        font.setPointSize(9*h_zoom)
        self.atom_mapping_file_location_label.setFont(font)
        self.atom_mapping_file_location_label.setEnabled(False)
        # atom mappping file location
        self.atom_mapping_file_location_comboBox = QtWidgets.QComboBox(self.subwin_ps_top_input)
        self.atom_mapping_file_location_comboBox.setGeometry(QtCore.QRect(1370*h_zoom, 45*h_zoom, 250*h_zoom, 22*h_zoom))
        self.atom_mapping_file_location_comboBox.setObjectName("atom_mapping_file_location_comboBox")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        font.setPointSize(9*h_zoom)
        self.atom_mapping_file_location_comboBox.setFont(font)
        self.atom_mapping_file_location_comboBox.setEnabled(False)
        self.atom_mapping_file_location_comboBox.setEditable(False)
        self.atom_mapping_file_location_list = ["../KEGGAAMBCsRCsFiles_Ori_Half/"]
        self.atom_mapping_file_location_comboBox.addItems(self.atom_mapping_file_location_list)
        self.atom_mapping_file_location_comboBox.setCurrentIndex(0)


        # search(start/stop)
        self.search_pushButton = QtWidgets.QPushButton(self.subwin_ps_top_input, )
        # self.search_pushButton.setGeometry(QtCore.QRect(950, 175, 125, 22))
        self.search_pushButton.setGeometry(QtCore.QRect(1445*h_zoom, 91*h_zoom, 100*h_zoom, 100*h_zoom))
        self.search_pushButton.setObjectName("search_pushButton")
        self.search_pushButton.setText("Start")
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setPointSize(12*h_zoom)
        self.search_pushButton.setFont(font)
        radius = str(int(100*h_zoom*0.5))
        self.search_pushButton.setStyleSheet("".join(["QPushButton{border-radius:", radius, "px}"]) +
                                             "QPushButton{background-color:rgb(220, 220, 220)}" +
                                             "QPushButton:hover{background-color:rgb(230, 230, 230)}" +
                                             "QPushButton:pressed{background-color:rgb(230, 230, 230)}" +
                                             "QPushButton{border:2px solid rgb(200, 200, 200)}" +
                                             "QPushButton:hover{border:2px solid rgb(210, 210, 210)}" +
                                             "QPushButton:pressed{border-style: hidden;}"
                                             )
        # # signal and slot
        self.search_pushButton.clicked.connect(self.forbid_input)
        self.search_pushButton.clicked.connect(self.start_stop_search)
        #


    def set_ps_middle_table_widgets(self, h_zoom=1.0):
        """
        """
        #
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setPointSize(10*h_zoom)
        #
        # table of pathway.
        self.subwin_ps_middle_pathtable = QtWidgets.QTableWidget()
        self.subwin_ps_middle_pathtable.setObjectName('subwin_ps_middle_pathtable')
        self.subwin_ps_middle_pathtable.setRowCount(0)
        self.subwin_ps_middle_pathtable.setColumnCount(11)
        # # self.subwin_ps_middle_pathtable.setIconSize(QtCore.QSize(100, 100))
        self.subwin_ps_middle_pathtable.setHorizontalHeaderLabels(["Feasibility", "TotalLength", "EndoLength", "HeterLength", "InfLength",
                                                                   "AtomUtilization", "AtomConservation", "MetabolicFlux",
                                                                   "S-Details",
                                                                   "I-Details", "M-Details",
                                                                   ])
        self.subwin_ps_middle_pathtable.setFont(font)
        # self.subwin_ps_middle_pathtable.resizeRowsToContents()
        # self.subwin_ps_middle_pathtable.resizeColumnsToContents()
        self.subwin_ps_middle_pathtable.setAlternatingRowColors(True)
        # self.subwin_ps_middle_pathtable.verticalHeader().setVisible(False)
        # self.subwin_ps_middle_pathtable.horizontalHeader().setVisible(False)
        self.subwin_ps_middle_pathtable.horizontalHeader().setStretchLastSection(True)
        # # self.subwin_ps_middle_pathtable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.subwin_ps_middle_pathtable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.subwin_ps_middle_pathtable.setShowGrid(True)
        self.subwin_ps_middle_pathtable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        # # self.subwin_ps_middle_pathtable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        # # self.subwin_ps_middle_pathtable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectColumns)
        # """
        column_length_list = [100, 110, 110, 110, 110, 140, 160, 140, 300, 300, 200]
        column_length_list = [item*h_zoom for item in column_length_list]
        for idx_column, column_length in enumerate(column_length_list):
            self.subwin_ps_middle_pathtable.setColumnWidth(idx_column, column_length)
        # """
        #
        data = QtWidgets.QTableWidgetItem("")
        data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)


        # grid layout of pathway.
        self.subwin_ps_middle_path_gridlayout = QtWidgets.QGridLayout(self.subwin_ps_middle_table)
        self.subwin_ps_middle_path_gridlayout.setObjectName('subwin_ps_middle_path_gridlayout')
        self.subwin_ps_middle_path_gridlayout.addWidget(self.subwin_ps_middle_pathtable)

        # signal and slot of pathway.
        # # self.subwin_ps_middle_pathtable.itemClicked.connect(self.show_one_pathway)
        # # self.subwin_ps_middle_pathtable.itemDoubleClicked.connect(self.show_one_pathway)
        # # self.subwin_ps_middle_pathtable.cellClicked.connect(self.show_one_pathway)
        # # self.subwin_ps_middle_pathtable.cellDoubleClicked.connect(self.show_one_pathway)
        self.subwin_ps_middle_pathtable.cellPressed.connect(self.show_one_pathway)
        #


    def set_ps_down_text_widgets(self, h_zoom=1.0):
        """
        """
        #
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setPointSize(10*h_zoom)
        #
        # clear button of tips.
        self.subwin_ps_down_clearPushButton = QtWidgets.QPushButton()
        self.subwin_ps_down_clearPushButton.setObjectName('subwin_ps_down_clearPushButton')
        self.subwin_ps_down_clearPushButton.setText('Clear')
        # # self.subwin_ps_down_clearPushButton.setIcon(QtGui.QIcon('images/Clear.ico'))
        # # self.subwin_ps_down_clearPushButton.setAutoExclusive(False)
        self.subwin_ps_down_clearPushButton.setStyleSheet("QPushButton{background-color:rgb(220, 220, 220)}"
                                                          "QPushButton:hover{background-color:rgb(230, 230, 230)}"
                                                          "QPushButton:pressed{background-color:rgb(230, 230, 230)}"
                                                          "QPushButton{border:2px solid rgb(200, 200, 200)}"
                                                          "QPushButton:hover{border:2px solid rgb(210, 210, 210)}"
                                                          "QPushButton:pressed{border-style: hidden;}"
                                                          )
        self.subwin_ps_down_clearPushButton.setFixedSize(100, 30)
        # # self.subwin_ps_down_clearPushButton.setToolTip("Clear Prompts.")
        self.subwin_ps_down_clearPushButton.setFont(font)


        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setPointSize(10*h_zoom)
        #
        # text browser of tips.
        # # self.subwin_ps_down_textbrower = QtWidgets.QTextBrowser()
        self.subwin_ps_down_textbrower = NewAppendTextBrowser()
        self.subwin_ps_down_textbrower.setObjectName('subwin_ps_down_textbrower')
        self.subwin_ps_down_textbrower.setFont(font)
        self.subwin_ps_down_textbrower.setStyleSheet('border-width: 1px;border-style: solid;border-color: rgb(210, 210, 210);background-color: rgb(255, 255, 255);')
        self.subwin_ps_down_textbrower.setOpenExternalLinks(True)
        self.subwin_ps_down_textbrower.setText('PyMiner . . . ')
        self.subwin_ps_down_textbrower.append('')
        # # self.subwin_ps_down_textbrower.append('Author: Xinfang Song')
        # # self.subwin_ps_down_textbrower.append('E-mail: sxf15@mails.tsinghua.edu.cn')
        # # self.subwin_ps_down_textbrower.append('Address: Engineering Research Center of Intelligent Technologies and Equipments for Saving Energy and Increasing Benefits (Ministry of Education), Department of Automation, Tsinghua University, Beijing 100084, China')
        # # self.subwin_ps_down_textbrower.append('...')
        #
        # self.subwin_ps_down_textbrower.moveCursor(self.subwin_ps_down_textbrower.textCursor().End)


        # grid layout of tips.
        self.subwin_ps_down_text_gridlayout = QtWidgets.QGridLayout(self.subwin_ps_down_text)
        self.subwin_ps_down_text_gridlayout.setObjectName('subwin_ps_down_text_gridlayout')
        self.subwin_ps_down_text_gridlayout.addWidget(self.subwin_ps_down_clearPushButton, 0, 0, QtCore.Qt.AlignRight)
        self.subwin_ps_down_text_gridlayout.addWidget(self.subwin_ps_down_textbrower, 1, 0)

        # signal and slot of tips.
        self.subwin_ps_down_clearPushButton.clicked.connect(self.clear_tips)
        #


    def set_ps_down_image_widgets(self, h_zoom=1.0, icon_size=180):
        """
        """
        #
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setPointSize(10*h_zoom)
        #
        # table of image.
        # # self.subwin_ps_down_imagetable = QtWidgets.QTableWidget()
        self.subwin_ps_down_imagetable = NewEventTableWidget()
        self.subwin_ps_down_imagetable.setObjectName('subwin_ps_down_imagetable')
        self.subwin_ps_down_imagetable.setRowCount(6)
        self.subwin_ps_down_imagetable.setColumnCount(41)
        self.subwin_ps_down_imagetable.setIconSize(QtCore.QSize(icon_size, icon_size))
        """
        self.subwin_ps_down_imagetable.setHorizontalHeaderLabels(["", "1", "", "2", "", "3", "", "4", "", "5",
                                                                  "", "6", "", "7", "", "8", "", "9", "", "10",
                                                                  "", "11", "", "12", "", "13", "", "14", "", "15",
                                                                  "", "16", "", "17", "", "18", "", "19", "", "20", ""])
        """
        self.subwin_ps_down_imagetable.setHorizontalHeaderLabels([""]*41)
        self.subwin_ps_down_imagetable.setVerticalHeaderLabels(["KEGG", "ChEBI/Rhea", "MetaCyc", "KndPad", " ", " "])
        # # self.subwin_ps_down_imagetable.resizeRowsToContents()
        # # self.subwin_ps_down_imagetable.resizeColumnsToContents()
        self.subwin_ps_down_imagetable.setAlternatingRowColors(True)
        # self.subwin_ps_down_imagetable.verticalHeader().setVisible(False)
        self.subwin_ps_down_imagetable.horizontalHeader().setVisible(False)
        # # self.subwin_ps_down_imagetable.horizontalHeader().setStretchLastSection(True)
        # # self.subwin_ps_down_imagetable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.subwin_ps_down_imagetable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.subwin_ps_down_imagetable.setShowGrid(False)
        self.subwin_ps_down_imagetable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        # # self.subwin_ps_down_imagetable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        # # self.subwin_ps_down_imagetable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectColumns)
        # self.subwin_ps_down_imagetable.setFont(font)
        self.icon_last_size = icon_size
        self.subwin_ps_down_imagetable.newWheelEvent.connect(self.update_image_size)

        # grid layout of image.
        self.subwin_ps_down_image_gridlayout = QtWidgets.QGridLayout(self.subwin_ps_down_image)
        self.subwin_ps_down_image_gridlayout.setObjectName('subwin_ps_down_image_gridlayout')
        self.subwin_ps_down_image_gridlayout.addWidget(self.subwin_ps_down_imagetable)
        #




    def source_update_input(self, _):
        """
            source completer;
        """
        #
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        #
        new_input_text = self.source_lineEdit.text()
        # complete the new input
        if new_input_text.endswith(';'):
            new_input_text += ' '
        #
        # initialize
        if new_input_text == "":
            self.source_input_ori = ""
            self.source_lineEdit.setText(self.source_input_ori)
            self.subwin_ps_down_textbrower.append( "".join(["Source: ", self.source_lineEdit.text(), '\n']) )
            #
            # Retro-search will be used if source in null.
            self.isRetro_checkBox.setChecked(True)
            self.isRetro_checkBox.setEnabled(False)
            #
            # Source and Host Organism can't be null simultaneously.
            # "" : ""
            # bsu: Bacillus subtilis subsp. subtilis 168
            # eco: Escherichia coli K-12 MG1655
            # kpn: Klebsiella pneumoniae subsp. pneumoniae MGH 78578
            # ppu: Pseudomonas putida KT2440
            # sce: Saccharomyces cerevisiae S288c
            # syz: Synechocystis sp. PCC 6803(Cyanobacteria)
            # ...
            # # self.organism_list = ['bsu', 'eco', 'kpn', 'ppu', 'sce', 'syz']
            host = self.host_organism_comboBox.currentText()
            index_host = self.organism_list.index(host)
            self.host_organism_comboBox.clear()
            self.host_organism_comboBox.addItems(self.organism_list[1:])
            if index_host == 0:
                self.host_organism_comboBox.setCurrentIndex(1)  # Escherichia coli
            else:
                self.host_organism_comboBox.setCurrentIndex(index_host-1)
            #
            return None
        else:
            # Source and Host Organism can't be null simultaneously.
            # "" : ""
            # bsu: Bacillus subtilis subsp. subtilis 168
            # eco: Escherichia coli K-12 MG1655
            # kpn: Klebsiella pneumoniae subsp. pneumoniae MGH 78578
            # ppu: Pseudomonas putida KT2440
            # sce: Saccharomyces cerevisiae S288c
            # syz: Synechocystis sp. PCC 6803(Cyanobacteria)
            # ...
            # # self.organism_list = ['', 'bsu', 'eco', 'kpn', 'ppu', 'sce', 'syz']
            host = self.host_organism_comboBox.currentText()
            index_host = self.organism_list.index(host)
            self.host_organism_comboBox.clear()
            self.host_organism_comboBox.addItems(self.organism_list)
            self.host_organism_comboBox.setCurrentIndex(index_host)

        #
        if self.source_input_ori == "":
            self.source_input_ori = new_input_text
            self.source_lineEdit.setText(self.source_input_ori)
            self.subwin_ps_down_textbrower.append( "".join(["Source: ", self.source_lineEdit.text(), '\n']) )
            #
            self.isRetro_checkBox.setEnabled(True)
            #
            if self.isSmart_checkBox.isChecked():
                self.isRetro_checkBox.setEnabled(False)
            #
            return None
        #
        # delete(<--)
        if self.source_input_ori.startswith(new_input_text):
            self.source_input_ori = new_input_text
            self.source_lineEdit.setText(self.source_input_ori)
        #
        new_input_list = new_input_text.split('; ')
        # remove all null items("")
        new_input_list = [item for item in new_input_list  if item != ""]
        try:
            new_item = new_input_list[-1]
        except IndexError:
            new_item = ""
        ori_input_list = self.source_input_ori.split('; ')
        #
        # clear and append
        if (new_item not in ori_input_list) and (new_item in self.comp_inputprompt_list):
            # delete original invalid input
            self.source_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.comp_inputprompt_list:
                    self.source_input_ori += "".join([ori_item, "; "])
            # append new valid input
            self.source_input_ori += "".join([new_item, "; "])
        elif new_item in self.comp_inputprompt_list:
            # unify original input
            self.source_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.comp_inputprompt_list:
                    self.source_input_ori += "".join([ori_item, "; "])
        #
        self.source_lineEdit.setText(self.source_input_ori)
        self.subwin_ps_down_textbrower.append( "".join(["Source: ", self.source_lineEdit.text(), '\n']) )
        #
        # mult-sources --> retro-search method must be selected.
        if (len(self.source_lineEdit.text().split('; ')) >= 3):
            self.isRetro_checkBox.setChecked(True)
            self.isRetro_checkBox.setEnabled(False)
        else:
            # # self.isRetro_checkBox.setChecked(False)
            self.isRetro_checkBox.setEnabled(True)
        #
        if self.isSmart_checkBox.isChecked():
            self.isRetro_checkBox.setEnabled(False)
        #
        return None
        #


    def target_update_input(self, _):
        """
            target completer;
        """
        #
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        new_input_text = self.target_lineEdit.text()
        #
        # # invalid input
        if new_input_text.strip() not in self.comp_inputprompt_list:
            self.subwin_ps_down_textbrower.append( "".join(["Target: ", self.target_lineEdit.text()]) )
            self.subwin_ps_down_textbrower.append("Invalid Target input.\n")
        else:
            self.subwin_ps_down_textbrower.append( "".join(["Target: ", self.target_lineEdit.text(), "\n"]) )

        return None
        #


    def avoid_compounds_update_input(self, _):
        """
            avoid compounds completer;
        """
        #
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        #
        new_input_text = self.avoid_compounds_lineEdit.text()
        # complete the new input
        if new_input_text.endswith(';'):
            new_input_text += ' '
        #
        # initialize
        if new_input_text == "":
            self.avoid_compounds_input_ori = ""
            self.avoid_compounds_lineEdit.setText(self.avoid_compounds_input_ori)
            self.subwin_ps_down_textbrower.append( "".join(["Avoid Compounds: ", self.avoid_compounds_lineEdit.text(), '\n']) )
            return None
        #
        # delete(<--)
        if self.avoid_compounds_input_ori.startswith(new_input_text):
            self.avoid_compounds_input_ori = new_input_text
            self.avoid_compounds_lineEdit.setText(self.avoid_compounds_input_ori)
        #
        new_input_list = new_input_text.split('; ')
        # remove all null items("")
        new_input_list = [item for item in new_input_list  if item != ""]
        try:
            new_item = new_input_list[-1]
        except IndexError:
            new_item = ""
        ori_input_list = self.avoid_compounds_input_ori.split('; ')
        #
        # clear and append
        if (new_item not in ori_input_list) and (new_item in self.comp_inputprompt_list):
            # delete original invalid input
            self.avoid_compounds_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.comp_inputprompt_list:
                    self.avoid_compounds_input_ori += "".join([ori_item, "; "])
            # append new valid input
            self.avoid_compounds_input_ori += "".join([new_item, "; "])
        elif new_item in self.comp_inputprompt_list:
            # unify original input
            self.avoid_compounds_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.comp_inputprompt_list:
                    self.avoid_compounds_input_ori += "".join([ori_item, "; "])
        #
        self.avoid_compounds_lineEdit.setText(self.avoid_compounds_input_ori)
        self.subwin_ps_down_textbrower.append( "".join(["Avoid Compounds: ", self.avoid_compounds_lineEdit.text(), '\n']) )
        return None
        #


    def avoid_reactions_update_input(self, _):
        """
            avoid reactions completer;
        """
        #
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        #
        new_input_text = self.avoid_reactions_lineEdit.text()
        # complete the new input
        if new_input_text.endswith(';'):
            new_input_text += ' '
        #
        # initialize
        if new_input_text == "":
            self.avoid_reactions_input_ori = ""
            self.avoid_reactions_lineEdit.setText(self.avoid_reactions_input_ori)
            self.subwin_ps_down_textbrower.append( "".join(["Avoid Reactions: ", self.avoid_reactions_lineEdit.text(), '\n']) )
            return None
        #
        # delete(<--)
        if self.avoid_reactions_input_ori.startswith(new_input_text):
            self.avoid_reactions_input_ori = new_input_text
            self.avoid_reactions_lineEdit.setText(self.avoid_reactions_input_ori)
        #
        new_input_list = new_input_text.split('; ')
        # remove all null items("")
        new_input_list = [item for item in new_input_list  if item != ""]
        try:
            new_item = new_input_list[-1]
        except IndexError:
            new_item = ""
        ori_input_list = self.avoid_reactions_input_ori.split('; ')
        #
        # clear and append
        if (new_item not in ori_input_list) and (new_item in self.rxn_inputprompt_list):
            # delete original invalid input
            self.avoid_reactions_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.rxn_inputprompt_list:
                    self.avoid_reactions_input_ori += "".join([ori_item, "; "])
            # append new valid input
            self.avoid_reactions_input_ori += "".join([new_item, "; "])
        elif new_item in self.rxn_inputprompt_list:
            # unify original input
            self.avoid_reactions_input_ori = ""
            for ori_item in ori_input_list:
                if ori_item in self.rxn_inputprompt_list:
                    self.avoid_reactions_input_ori += "".join([ori_item, "; "])
        #
        self.avoid_reactions_lineEdit.setText(self.avoid_reactions_input_ori)
        self.subwin_ps_down_textbrower.append( "".join(["Avoid Reactions: ", self.avoid_reactions_lineEdit.text(), '\n']) )
        return None
        #


    def selectDatabase(self, _):
        """
        """
        data_base = self.database_comboBox.currentText()
        #
        if data_base == "KEGG":
            self.resource_file_list = ["../sourcedata/KEGG_Ori_Cleaned_Half.xlsx"]
            self.general_cofactors_file_list = ["datas/KEGG_General_Cofactors.xlsx", "datas/KEGG_General_Cofactors_Empty.xlsx"]
            self.metabolic_network_file_list = ["../sourcedata/KEGG_Pathway_Search_Ori.xlsx"]
            self.atom_mapping_file_location_list = ["../KEGGAAMBCsRCsFiles_Ori_Half/"]
            #
            #
            self.comp_inputprompt_list = self.inputprompt_dict['KEGG']['Compound']
            self.rxn_inputprompt_list = self.inputprompt_dict['KEGG']['Reaction']
            #
            # source completer
            self.source_lineEdit.setCompleter(self.kegg_source_completer)
            # target completer
            self.target_lineEdit.setCompleter(self.kegg_target_completer)
            # avoid compounds completer
            self.avoid_compounds_lineEdit.setCompleter(self.kegg_avoid_compounds_completer)
            # avoid reactions completer
            self.avoid_reactions_lineEdit.setCompleter(self.kegg_avoid_reactions_completer)
        #
        elif data_base == "MetaCyc":
            self.resource_file_list = ["../sourcedata/KndPad_Ori_Cleaned_Half.xlsx"]
            self.general_cofactors_file_list = ["datas/MetaCyc_General_Cofactors.xlsx", "datas/MetaCyc_General_Cofactors_Empty.xlsx"]
            self.metabolic_network_file_list = ["../sourcedata/MetaCyc_Pathway_Search_Ori.xlsx"]
            self.atom_mapping_file_location_list = ["../KndPadAAMBCsRCsFiles_Ori_Half/"]
            #
            #
            self.comp_inputprompt_list = self.inputprompt_dict['MetaCyc']['Compound']
            self.rxn_inputprompt_list = self.inputprompt_dict['MetaCyc']['Reaction']
            #
            # source completer
            self.source_lineEdit.setCompleter(self.metacyc_source_completer)
            # target completer
            self.target_lineEdit.setCompleter(self.metacyc_target_completer)
            # avoid compounds completer
            self.avoid_compounds_lineEdit.setCompleter(self.metacyc_avoid_compounds_completer)
            # avoid reactions completer
            self.avoid_reactions_lineEdit.setCompleter(self.metacyc_avoid_reactions_completer)
        #
        elif data_base == "KndPad":
            self.resource_file_list = ["../sourcedata/KndPad_Ori_Cleaned_Half.xlsx"]
            self.general_cofactors_file_list = ["datas/KndPad_General_Cofactors.xlsx", "datas/KndPad_General_Cofactors_Empty.xlsx"]
            self.metabolic_network_file_list = ["../sourcedata/KndPad_Pathway_Search_Ori.xlsx"]
            self.atom_mapping_file_location_list = ["../KndPadAAMBCsRCsFiles_Ori_Half/"]
            #
            #
            self.comp_inputprompt_list = self.inputprompt_dict['KndPad']['Compound']
            self.rxn_inputprompt_list = self.inputprompt_dict['KndPad']['Reaction']
            #
            # source completer
            self.source_lineEdit.setCompleter(self.kndpad_source_completer)
            # target completer
            self.target_lineEdit.setCompleter(self.kndpad_target_completer)
            # avoid compounds completer
            self.avoid_compounds_lineEdit.setCompleter(self.kndpad_avoid_compounds_completer)
            # avoid reactions completer
            self.avoid_reactions_lineEdit.setCompleter(self.kndpad_avoid_reactions_completer)
        #
        self.resource_file_comboBox.clear()
        self.resource_file_comboBox.addItems(self.resource_file_list)
        self.resource_file_comboBox.setCurrentIndex(0)
        #
        self.general_cofactors_file_comboBox.clear()
        self.general_cofactors_file_comboBox.addItems(self.general_cofactors_file_list)
        self.general_cofactors_file_comboBox.setCurrentIndex(0)
        #
        self.metabolic_network_file_comboBox.clear()
        self.metabolic_network_file_comboBox.addItems(self.metabolic_network_file_list)
        self.metabolic_network_file_comboBox.setCurrentIndex(0)
        #
        self.atom_mapping_file_location_comboBox.clear()
        self.atom_mapping_file_location_comboBox.addItems(self.atom_mapping_file_location_list)
        self.atom_mapping_file_location_comboBox.setCurrentIndex(0)
        #


    def cofactor_update_input(self, _):
        """
        """
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        general_cofactors_file = self.general_cofactors_file_comboBox.currentText()
        self.subwin_ps_down_textbrower.append("General Cofactors File: " + general_cofactors_file + "\n")
        #
        if general_cofactors_file.endswith("_Empty.xlsx"):
            isEmpty = True
        else:
            isEmpty = False
        #
        data_base = self.database_comboBox.currentText()
        if data_base == "KEGG":
            if isEmpty:
                self.metabolic_network_file_list = ["../sourcedata/KEGG_Pathway_Search_Ori_withCofactors.xlsx"]
            else:
                self.metabolic_network_file_list = ["../sourcedata/KEGG_Pathway_Search_Ori.xlsx"]
        #
        if data_base == "MetaCyc":
            if isEmpty:
                self.metabolic_network_file_list = ["../sourcedata/MetaCyc_Pathway_Search_Ori_withCofactors.xlsx"]
            else:
                self.metabolic_network_file_list = ["../sourcedata/MetaCyc_Pathway_Search_Ori.xlsx"]
        #
        if data_base == "KndPad":
            if isEmpty:
                self.metabolic_network_file_list = ["../sourcedata/KndPad_Pathway_Search_Ori_withCofactors.xlsx"]
            else:
                self.metabolic_network_file_list = ["../sourcedata/KndPad_Pathway_Search_Ori.xlsx"]
        #
        self.metabolic_network_file_comboBox.clear()
        self.metabolic_network_file_comboBox.addItems(self.metabolic_network_file_list)
        self.metabolic_network_file_comboBox.setCurrentIndex(0)
        #


    def reset_methods(self, _):
        """
        """
        if self.isShortest_checkBox.isChecked():
            # # self.methods_list = ['bfs_wl', "bfs_naive"]
            self.methods_list = ["bfs_naive"]
            self.isTotal_checkBox.setChecked(True)
            self.isTotal_checkBox.setEnabled(False)
        else:
            self.methods_list = ["bfs_naive", "dfs_naive"]
            self.isTotal_checkBox.setEnabled(True)
        self.methods_comboBox.clear()
        self.methods_comboBox.addItems(self.methods_list)
        self.methods_comboBox.setCurrentIndex(0)
        #


    def isValid_isRetro(self, _):
        """
        """
        if self.isSmart_checkBox.isChecked():
            self.isRetro_checkBox.setEnabled(False)
        else:
            self.isRetro_checkBox.setEnabled(True)
        #
        # mult-sources
        if len( self.source_lineEdit.text().split('; ') ) >= 3:
            self.isRetro_checkBox.setChecked(True)
            self.isRetro_checkBox.setEnabled(False)
        #

    def isValid_isAerobic(self, _):
        """
        """
        if self.isFlux_checkBox.isChecked():
            self.isAerobic_checkBox.setEnabled(True)
        else:
            self.isAerobic_checkBox.setEnabled(False)


    def isPermit_update(self, _):
        """
        """
        # QMessageBox.information; QMessageBox.question; QMessageBox.warning; QMessageBox.ctitical; QMessageBox.about
        if self.isUpdate_checkBox.isChecked():
            isConfirmUpdate = QtWidgets.QMessageBox.warning(self.subwin_ps_top_input,"Warning", "Do you want to update the distilled metabolic network used for pathway search ?",
                                                            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if isConfirmUpdate == QtWidgets.QMessageBox.Yes:
                # # self.isUpdate_checkBox.setChecked(True)
                self.isUpdate_checkBox.setCheckState(2)
                #
                self.resource_file_label.setEnabled(False)
                self.resource_file_comboBox.setEnabled(False)
                self.general_cofactors_file_label.setEnabled(True)
                self.general_cofactors_file_comboBox.setEnabled(True)
                self.metabolic_network_file_label.setEnabled(False)
                self.metabolic_network_file_comboBox.setEnabled(False)
                self.atom_mapping_file_location_label.setEnabled(False)
                self.atom_mapping_file_location_comboBox.setEnabled(False)
            else:
                self.isUpdate_checkBox.setChecked(False)
                self.isUpdate_checkBox.setCheckState(0)
                #
        else:
            self.resource_file_label.setEnabled(False)
            self.resource_file_comboBox.setEnabled(False)
            self.general_cofactors_file_label.setEnabled(False)
            self.general_cofactors_file_comboBox.setEnabled(False)
            self.metabolic_network_file_label.setEnabled(False)
            self.metabolic_network_file_comboBox.setEnabled(False)
            self.atom_mapping_file_location_label.setEnabled(False)
            self.atom_mapping_file_location_comboBox.setEnabled(False)
            #


    def forbid_input(self, _):
        """
        """
        #
        start_stop = self.search_pushButton.text()
        if start_stop == "Stop":
            return None
        #
        self.source_label.setEnabled(False)
        self.source_lineEdit.setEnabled(False)
        self.target_label.setEnabled(False)
        self.target_lineEdit.setEnabled(False)
        self.avoid_compounds_label.setEnabled(False)
        self.avoid_compounds_lineEdit.setEnabled(False)
        self.avoid_reactions_label.setEnabled(False)
        self.avoid_reactions_lineEdit.setEnabled(False)
        #
        self.host_organism_label.setEnabled(False)
        self.host_organism_comboBox.setEnabled(False)
        self.database_label.setEnabled(False)
        self.database_comboBox.setEnabled(False)
        self.maximum_length_label.setEnabled(False)
        self.maximum_length_comboBox.setEnabled(False)
        self.maximum_time_label.setEnabled(False)
        self.maximum_time_comboBox.setEnabled(False)
        #
        self.diffThreshold_label.setEnabled(False)
        self.diffThreshold_comboBox.setEnabled(False)
        self.methods_label.setEnabled(False)
        self.methods_comboBox.setEnabled(False)
        #
        self.isTotal_checkBox.setEnabled(False)
        self.isShortest_checkBox.setEnabled(False)
        self.isRetro_checkBox.setEnabled(False)
        self.isSmart_checkBox.setEnabled(False)
        self.isInfeasible_checkBox.setEnabled(False)
        self.isTransfer_checkBox.setEnabled(False)
        self.isAerobic_checkBox.setEnabled(False)
        self.isFlux_checkBox.setEnabled(False)
        self.isUpdate_checkBox.setEnabled(False)
        #
        self.resource_file_label.setEnabled(False)
        self.resource_file_comboBox.setEnabled(False)
        self.general_cofactors_file_label.setEnabled(False)
        self.general_cofactors_file_comboBox.setEnabled(False)
        self.metabolic_network_file_label.setEnabled(False)
        self.metabolic_network_file_comboBox.setEnabled(False)
        self.atom_mapping_file_location_label.setEnabled(False)
        self.atom_mapping_file_location_comboBox.setEnabled(False)
        #
        # # self.search_pushButton.setEnabled(False)
        #


    def permit_input(self, _):
        """
        """
        self.source_label.setEnabled(True)
        self.source_lineEdit.setEnabled(True)
        self.target_label.setEnabled(True)
        self.target_lineEdit.setEnabled(True)
        self.avoid_compounds_label.setEnabled(True)
        self.avoid_compounds_lineEdit.setEnabled(True)
        self.avoid_reactions_label.setEnabled(True)
        self.avoid_reactions_lineEdit.setEnabled(True)
        #
        self.host_organism_label.setEnabled(True)
        self.host_organism_comboBox.setEnabled(True)
        self.database_label.setEnabled(True)
        self.database_comboBox.setEnabled(True)
        self.maximum_length_label.setEnabled(True)
        self.maximum_length_comboBox.setEnabled(True)
        self.maximum_time_label.setEnabled(True)
        self.maximum_time_comboBox.setEnabled(True)
        #
        self.diffThreshold_label.setEnabled(True)
        self.diffThreshold_comboBox.setEnabled(True)
        self.methods_label.setEnabled(True)
        self.methods_comboBox.setEnabled(True)
        #
        self.isTotal_checkBox.setEnabled(True)
        self.isShortest_checkBox.setEnabled(True)
        self.isRetro_checkBox.setEnabled(True)
        self.isSmart_checkBox.setEnabled(True)
        self.isInfeasible_checkBox.setEnabled(True)
        self.isTransfer_checkBox.setEnabled(True)
        self.isAerobic_checkBox.setEnabled(True)
        self.isFlux_checkBox.setEnabled(True)
        self.isUpdate_checkBox.setEnabled(True)
        #
        self.resource_file_label.setEnabled(True)
        self.resource_file_comboBox.setEnabled(True)
        self.general_cofactors_file_label.setEnabled(True)
        self.general_cofactors_file_comboBox.setEnabled(True)
        self.metabolic_network_file_label.setEnabled(True)
        self.metabolic_network_file_comboBox.setEnabled(True)
        self.atom_mapping_file_location_label.setEnabled(True)
        self.atom_mapping_file_location_comboBox.setEnabled(True)
        #
        # shortest --> total
        if self.isShortest_checkBox.isChecked():
            self.isTotal_checkBox.setChecked(True)
            self.isTotal_checkBox.setEnabled(False)
        #
        # smart search --> retro search
        if self.isSmart_checkBox.isChecked():
            self.isRetro_checkBox.setEnabled(False)
        #
        # flux --> aerobic culture
        if not self.isFlux_checkBox.isChecked():
            self.isAerobic_checkBox.setEnabled(False)
        #
        # update --> ...
        if not self.isUpdate_checkBox.isChecked():
            self.resource_file_label.setEnabled(False)
            self.resource_file_comboBox.setEnabled(False)
            self.general_cofactors_file_label.setEnabled(False)
            self.general_cofactors_file_comboBox.setEnabled(False)
            self.metabolic_network_file_label.setEnabled(False)
            self.metabolic_network_file_comboBox.setEnabled(False)
            self.atom_mapping_file_location_label.setEnabled(False)
            self.atom_mapping_file_location_comboBox.setEnabled(False)
        #
        # sources --> retro search
        if len( self.source_lineEdit.text().split('; ') ) >= 3:
            self.isRetro_checkBox.setChecked(True)
            self.isRetro_checkBox.setEnabled(False)
        #
        # # self.search_pushButton.setEnabled(True)
        self.search_pushButton.setText("Start")
        #


    def getInputValue(self, isprompt=True):
        """
        """
        #
        pattern_comp = re.compile("(C\d{5}|Met\d{6}-[mck])")
        pattern_rxn = re.compile("(R\d{5}|Rxn\d{6}-[mrk])")
        #
        source = self.source_lineEdit.text()
        source = re.findall(pattern_comp, source)
        source = set(source)
        target = self.target_lineEdit.text()
        try:
            target = re.search(pattern_comp, target).group(0)
        except AttributeError:
            target = ""
        avoid_compounds = self.avoid_compounds_lineEdit.text()
        avoid_compounds = re.findall(pattern_comp, avoid_compounds)
        avoid_compounds = set(avoid_compounds)
        avoid_reactions = self.avoid_reactions_lineEdit.text()
        avoid_reactions = re.findall(pattern_rxn, avoid_reactions)
        avoid_reactions = set(avoid_reactions)
        #
        host_organism = self.host_organism_comboBox.currentText()
        database = self.database_comboBox.currentText()
        maximum_length = int(self.maximum_length_comboBox.currentText())
        maximum_time = int(self.maximum_time_comboBox.currentText())
        #
        diffThreshold = float(self.diffThreshold_comboBox.currentText())
        methods = self.methods_comboBox.currentText()
        #
        isTotal = self.isTotal_checkBox.isChecked()
        isShortest = self.isShortest_checkBox.isChecked()
        #
        isRetro = self.isRetro_checkBox.isChecked()
        isSmart = self.isSmart_checkBox.isChecked()
        #
        isInfeasible = self.isInfeasible_checkBox.isChecked()
        isTransfer = self.isTransfer_checkBox.isChecked()
        #
        isAerobic = self.isAerobic_checkBox.isChecked()
        isFlux = self.isFlux_checkBox.isChecked()
        #
        isUpdate = self.isUpdate_checkBox.isChecked()
        #
        resource_file = self.resource_file_comboBox.currentText()
        general_cofactors_file = self.general_cofactors_file_comboBox.currentText()
        metabolic_network_file = self.metabolic_network_file_comboBox.currentText()
        atom_mapping_file_location = self.atom_mapping_file_location_comboBox.currentText()


        if isprompt:
            self.subwin_ps_down_textbrower.append("<font color='blue', font-weight='bold', size=6>Input:</font>")
            self.subwin_ps_down_textbrower.append("Sources: "                       + json.dumps(sorted(source)).replace("[", "{").replace("]", "}") )
            self.subwin_ps_down_textbrower.append("Target: "                        + target)
            self.subwin_ps_down_textbrower.append("Avoid Compounds: "               + json.dumps(sorted(avoid_compounds)).replace("[", "{").replace("]", "}") )
            self.subwin_ps_down_textbrower.append("Avoid Reactions: "               + json.dumps(sorted(avoid_reactions)).replace("[", "{").replace("]", "}") )
            self.subwin_ps_down_textbrower.append("Host Organism: "                 + host_organism)
            self.subwin_ps_down_textbrower.append("Database: "                      + database)
            self.subwin_ps_down_textbrower.append("Maximum Length: "                + str(maximum_length) )
            self.subwin_ps_down_textbrower.append("Maximum Time(s): "               + str(maximum_time) )
            self.subwin_ps_down_textbrower.append("Similarity Difference Threshold: "   + str(diffThreshold) )
            self.subwin_ps_down_textbrower.append("Search Method: "                 + methods)
            self.subwin_ps_down_textbrower.append("Total: "                         + str(isTotal) )
            self.subwin_ps_down_textbrower.append("Shortest: "                      + str(isShortest) )
            self.subwin_ps_down_textbrower.append("Retro Search: "                  + str(isRetro) )
            self.subwin_ps_down_textbrower.append("Smart Search: "                  + str(isSmart) )
            self.subwin_ps_down_textbrower.append("Infeasibility: "                 + str(isInfeasible) )
            self.subwin_ps_down_textbrower.append("Atom Trace: "                    + str(isTransfer) )
            self.subwin_ps_down_textbrower.append("Aerobic Culture: "               + str(isAerobic) )
            self.subwin_ps_down_textbrower.append("Flux: "                          + str(isFlux) )
            self.subwin_ps_down_textbrower.append("Update: "                        + str(isUpdate) )
            self.subwin_ps_down_textbrower.append("Resource File: "                 + resource_file)
            self.subwin_ps_down_textbrower.append("General Cofactors File: "        + general_cofactors_file)
            self.subwin_ps_down_textbrower.append("Distilled Metabolic Network File: "        + metabolic_network_file)
            self.subwin_ps_down_textbrower.append("Atom Mapping Files: "            + atom_mapping_file_location)
            self.subwin_ps_down_textbrower.append("")

        input_value = {"Source": source, "Target": target, "Avoid Compounds": avoid_compounds, "Avoid Reactions": avoid_reactions,
                       "Host Organism": host_organism, "Database": database, "Maximum Length": maximum_length, "Maximum Time(s)": maximum_time,
                       "DiffThreshold": diffThreshold, "Methods": methods,
                       "IsTotal":isTotal, "IsShortest": isShortest,
                       "IsRetro": isRetro, "IsSmart": isSmart,
                       "IsInfeasible": isInfeasible, "IsTransfer": isTransfer,
                       "IsAerobic": isAerobic, "IsFlux": isFlux,
                       "IsUpdate": isUpdate,
                       "Resource File": resource_file,
                       "General Cofactors File": general_cofactors_file,
                       "Metabolic Network File": metabolic_network_file,
                       "Atom Mapping File Location": atom_mapping_file_location,
                       }
        return input_value


    def initSearchThread(self):
        """
        """
        self.search_thread = SearchThread()
        self.search_thread.Daemon = True
        self.search_thread.finished_signal.connect(self.show_all_pathways)
        self.search_thread.finished_signal.connect(self.permit_input)
        self.search_thread.prompts_signal.connect(self.append_prompts)
        #

    def start_stop_search(self, _=None, icon_size=180):
        """
        """
        #
        start_stop = self.search_pushButton.text()
        #
        # stop pathway search
        if start_stop == "Stop":
            self.search_thread.terminate()
            self.search_thread.wait()
            self.search_thread.quit()
            self.search_pushButton.setText("Start")
            self.permit_input(None)
            return None
        #
        # start pathway search()
        self.search_pushButton.setText("Stop")
        self.icon_last_size = icon_size
        #
        self.subwin_ps_middle_pathtable.setRowCount(0)
        self.cell_last_row = -1
        # Removes all items not in the headers from the view. This will also remove all selections. The table dimensions stay the same.
        self.subwin_ps_down_imagetable.clearContents()
        self.subwin_ps_down_tabwidget.setCurrentIndex(0)
        #
        # input parameters
        input_value = self.getInputValue()
        # Pathway search thread.
        self.search_thread.input_value = input_value
        self.search_thread.start()
        #


    def show_all_pathways(self, pathwayinfo_list):
        """
        """
        #
        # '0-Feasibility';
        # '1-TotalLength'; '2-EndoLength'; '3-HeterLength'; '4-InfLength';
        # '5-AtomUtilization'; '6-AtomConservation'; '7-Flux';
        # '8-S-Details';
        # '9-I-Details'; '10-M-Details';
        # '11-Route'
        #
        # analysis the results
        self.pathwayEnumList = pathwayinfo_list[0]
        # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number);
        self.rxnDistInfoDict = pathwayinfo_list[1]
        # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES
        self.compDistInfoDict = pathwayinfo_list[2]

        # PySide6.QtWidgets.QTableWidget.setRowCount(rows)
        # Sets the number of rows in this table’s model to rows.
        # If this is less than rowCount() , the data in the unwanted rows is discarded.
        self.subwin_ps_middle_pathtable.setRowCount(0)
        # # self.subwin_ps_middle_pathtable.clear()
        # # self.subwin_ps_middle_pathtable.clearContents
        #
        s_details_list = list()
        last_s_details = ""
        #
        for idx_row, pathway_info in enumerate(self.pathwayEnumList):
            self.subwin_ps_middle_pathtable.insertRow(idx_row)
            # Exclude the relative route of this pathway's highlighted compound structure figures appended in atom transfer evaluating phase.
            for idx_column, item_value in enumerate(pathway_info[:-1]):
                if type(item_value) == list:
                    item_value = json.dumps(item_value)
                else:
                    item_value = str(item_value)
                if idx_column == 8:
                    current_s_details = pathway_info[idx_column]
                    AddTableItem.addColorItem(self.subwin_ps_middle_pathtable, idx_row, idx_column, last_s_details, current_s_details)
                elif idx_column >= 9:
                    AddTableItem.addPlainItem(self.subwin_ps_middle_pathtable, idx_row, idx_column, item_value)
                else:
                    data = QtWidgets.QTableWidgetItem(item_value)
                    data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
                    self.subwin_ps_middle_pathtable.setItem(idx_row, idx_column, data)
            #
            # update
            last_s_details = current_s_details
            # record
            s_details_list.append(current_s_details)
        #
        total_comp2rxn_set = set()
        for s_details in s_details_list:
            total_comp2rxn_set |= set(s_details.split('-->'))
        #
        # generate the structure pictures of all compounds in all the pathways.
        for comp2rxn in total_comp2rxn_set:
            try:
                smiles_comp = self.compDistInfoDict[comp2rxn][4]
                image_comp = "".join(['results/ps_results/', comp2rxn, '.png'])
                mol_comp = Chem.MolFromSmiles(smiles_comp)
                Chem_Draw.MolToFile(mol_comp, image_comp, size=(600, 600))
            except KeyError:
                pass
        #
        # show the details of the first searched pathway.
        if self.subwin_ps_middle_pathtable.rowCount() != 0:
            self.show_one_pathway(cell_row = 0, cell_column = 0)
        #


    def show_one_pathway(self, cell_row=0, cell_column=0, column_s_details=8, icon_size=180, isLeftUp=True):
        """
        """
        #
        # # print("Row:", cell_row, "Column:", cell_column)
        #
        self.cell_last_row = cell_row
        # Reuse the size of icon updated by mousewheel while we choose one pathway by mouse.
        if icon_size == 180:
            icon_size = self.icon_last_size
        #
        # analysis the results
        # Direction; Compound List; Compound Set; Compound Pair(diffThreshold); Reference(KEGG, Rhea, MetaCyc, KndPad, EC Number); Coefficient List;
        # self.rxnDistInfoDict
        # Name; Branching Factor; Reactions; Reference(KEGG, ChEBI, MetaCyc, KndPad); SMILES;
        # self.compDistInfoDict
        #
        self.subwin_ps_down_tabwidget.setCurrentIndex(1)
        database = self.database_comboBox.currentText()
        #
        # reset...
        # # self.subwin_ps_down_imagetable.clearContents()
        self.subwin_ps_down_imagetable.setRowCount(0)
        self.subwin_ps_down_imagetable.setAlternatingRowColors(True)
        self.subwin_ps_down_imagetable.setHorizontalHeaderLabels([""]*41)
        #
        # reset...
        if isLeftUp:
            # scroll left
            self.subwin_ps_down_imagetable.horizontalScrollBar().setValue(self.subwin_ps_down_imagetable.horizontalScrollBar().minimum())
            # scroll up
            self.subwin_ps_down_imagetable.verticalScrollBar().setValue(self.subwin_ps_down_imagetable.verticalScrollBar().minimum())
        #
        # reset...
        if database == "KEGG":
            row_size = 4
            self.subwin_ps_down_imagetable.setRowCount(4)
            self.subwin_ps_down_imagetable.setVerticalHeaderLabels(["KEGG", "ChEBI/Rhea", " ", " "])
        if database in ["MetaCyc", "KndPad"]:
            row_size = 6
            self.subwin_ps_down_imagetable.setRowCount(6)
            self.subwin_ps_down_imagetable.setVerticalHeaderLabels(["KEGG", "ChEBI/Rhea", "MetaCyc", "KndPad", " ", " "])
        #
        # reset...
        self.subwin_ps_down_imagetable.setIconSize(QtCore.QSize(icon_size, icon_size));
        self.subwin_ps_down_imagetable.setRowHeight(row_size-1, icon_size)
        #
        # data...
        s_details = self.subwin_ps_middle_pathtable.cellWidget(cell_row, 8).text()
        s_details = re.sub('<.*?>', "", s_details)
        comp2rxn_list = s_details.split('-->')
        i_details = self.subwin_ps_middle_pathtable.cellWidget(cell_row, 9).text()
        i_details = json.loads(i_details)
        #
        # find the pathway corresponding to the current row.
        for pathway_info in self.pathwayEnumList:
            if pathway_info[8] == s_details:
                # With atom transfer evaluation.
                try:
                    sub_route = pathway_info[11]
                # Without atom transfer evaluation.
                except IndexError:
                    sub_route = ""
                break
        #
        self.subwin_ps_down_textbrower.append( "".join(["sub-route: ", sub_route]) )
        #
        for idx_column, comp2rxn in enumerate(comp2rxn_list):
            #
            # # comp2rxn = re.sub('<.*?>', "", comp2rxn)
            self.subwin_ps_down_imagetable.setColumnWidth(idx_column, icon_size)
            #
            # compound
            try:
                #
                compDistInfo = self.compDistInfoDict[comp2rxn]
                #
                idx_row = 0
                # KEGG Compound
                kegg_comp_list = compDistInfo[3][0]
                if len(kegg_comp_list) != 0:
                    kegg_comp = kegg_comp_list[0]
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, kegg_comp, kegg_comp, "KEGG", 'Compound')
                #
                idx_row += 1
                # ChEBI Compound
                chebi_comp_list = compDistInfo[3][1]
                if len(chebi_comp_list) != 0:
                    chebi_comp = chebi_comp_list[0]
                    chebi_comp_content = ":".join(["CHEBI", chebi_comp])
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, chebi_comp_content, chebi_comp, "ChEBI", 'Compound')
                #
                idx_row += 1
                # MetaCyc Compound
                metacyc_comp_list = compDistInfo[3][2]
                if (len(metacyc_comp_list) != 0) and ( database in ["MetaCyc", "KndPad"]):
                    metacyc_comp = metacyc_comp_list[0]
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, metacyc_comp, metacyc_comp, "MetaCyc", 'Compound')
                    idx_row += 1
                elif database in ["MetaCyc", "KndPad"]:
                    idx_row += 1
                #
                # KndPad Compound
                kndpad_comp_list = compDistInfo[3][3]
                if (len(kndpad_comp_list) != 0) and (database in ["MetaCyc", "KndPad"]):
                    kndpad_comp = kndpad_comp_list[0]
                    data = QtWidgets.QTableWidgetItem(kndpad_comp)
                    data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
                    self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, data)
                    idx_row += 1
                elif database in ["MetaCyc", "KndPad"]:
                    idx_row += 1
                #
                # Name
                name_comp = compDistInfo[0]
                data = QtWidgets.QTableWidgetItem(name_comp)
                data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
                self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, data)
                #
                idx_row += 1
                # Compound Structure
                image_comp =  "".join([comp2rxn, '.png'])
                image_comp = path.join(sub_route, image_comp)
                if not path.exists(image_comp):
                    image_comp = "".join([comp2rxn, '.png'])
                    image_comp = path.join('results', 'ps_results', image_comp)
                item = QtWidgets.QTableWidgetItem()
                item.setFlags(QtCore.Qt.ItemIsEnabled)
                item.setIcon(QtGui.QIcon(image_comp))
                self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, item)
                # # idx_row += 1
                #
            # reaction
            except KeyError:
                #
                rxnDistInfo = self.rxnDistInfoDict[comp2rxn]
                #
                idx_row = 0
                # KEGG Reaction
                kegg_rxn_list = rxnDistInfo[4][0]
                if len(kegg_rxn_list) != 0:
                    kegg_rxn = kegg_rxn_list[0]
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, kegg_rxn, kegg_rxn, "KEGG", 'Reaction')
                #
                idx_row += 1
                # Rhea Reaction
                rhea_rxn_list = rxnDistInfo[4][1]
                if len(rhea_rxn_list) != 0:
                    rhea_rxn = rhea_rxn_list[0]
                    rhea_rxn_content = ":".join(["RHEA", rhea_rxn])
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, rhea_rxn_content, rhea_rxn, "Rhea", 'Reaction')
                #
                idx_row += 1
                # MetaCyc Reaction
                metacyc_rxn_list = rxnDistInfo[4][2]
                if (len(metacyc_rxn_list) != 0) and ( database in ["MetaCyc", "KndPad"]):
                    metacyc_rxn = metacyc_rxn_list[0]
                    AddTableItem.addLinkItem(self.subwin_ps_down_imagetable, idx_row, idx_column, metacyc_rxn, metacyc_rxn, "MetaCyc", 'Reaction')
                    idx_row += 1
                elif database in ["MetaCyc", "KndPad"]:
                    idx_row += 1
                #
                # KndPad Reaction
                kndpad_rxn_list = rxnDistInfo[4][3]
                if (len(kndpad_rxn_list) != 0) and (database in ["MetaCyc", "KndPad"]):
                    kndpad_rxn = kndpad_rxn_list[0]
                    data = QtWidgets.QTableWidgetItem(kndpad_rxn)
                    data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
                    self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, data)
                    idx_row += 1
                elif database in ["MetaCyc", "KndPad"]:
                    idx_row += 1
                #
                # EC Number
                ec_number_list = rxnDistInfo[4][4]
                ec_number_list = [":".join(["EC", item]) for item in ec_number_list]
                ec_numbers = "; ".join(ec_number_list)
                data = QtWidgets.QTableWidgetItem(ec_numbers)
                data.setTextAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
                self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, data)
                #
                idx_row += 1
                # Reaction Direaction
                dir_rxn = rxnDistInfo[0]
                isFeasible = i_details[(idx_column-1)//2]
                dir_rxn_key = "/".join([dir_rxn, str(isFeasible)])
                dir_rxn_dict = {"=>/-1": "lr_red", "<=/-1": "lr_red", "<=>/-1": "bi_red", "=/-1": "bi_red",
                                "=>/0": "lr_blue", "<=/0": "lr_blue", "<=>/0": "bi_blue", "=/0": "bi_blue",
                                "=>/1": "lr_green", "<=/1": "lr_green", "<=>/1": "bi_green", "=/1": "bi_green",
                                }
                image_rxn_dir = dir_rxn_dict[dir_rxn_key]
                image_rxn_dir = "".join(["images/", image_rxn_dir, ".png"])
                item = QtWidgets.QTableWidgetItem()
                item.setFlags(QtCore.Qt.ItemIsEnabled)
                item.setIcon(QtGui.QIcon(image_rxn_dir))
                self.subwin_ps_down_imagetable.setItem(idx_row, idx_column, item)
                # # idx_row += 1
                #


    def update_image_size(self, xyz):
        """
            Update size by mousewheel
        """
        step = 6
        min_size = 60   # minimum size
        max_size = 600  # maximum size
        try:
            cell_row = self.cell_last_row
        except AttributeError:
            cell_row = -1
        if cell_row != -1:
            self.icon_last_size += ( step * xyz["Zoom_In_out"] )
            if self.icon_last_size < min_size:
                self.icon_last_size = min_size
            elif self.icon_last_size > max_size:
                self.icon_last_size = max_size
            self.show_one_pathway(cell_row = cell_row, icon_size = self.icon_last_size, isLeftUp = False)
        #


    def append_prompts(self, tips):
        """
        """
        self.subwin_ps_down_textbrower.append(tips)


    def clear_tips(self, _):
        """
        """
        self.subwin_ps_down_textbrower.clear()


    @classmethod
    def demoFunc(cls):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





