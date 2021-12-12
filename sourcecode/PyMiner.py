# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 19:05:47 2021

@author: CC-SXF
"""

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

from widgets.pathway_searcher_widgets import Pathway_Seracher_Widgets
from widgets.group_widgets import Group_StatusBar_Widgets


class Ui_MainWindow(object):

    def __init__(self):
        """ """
        pass

    def setupUi(self, MainWindow):
        #
        MainWindow.setObjectName("MainWindow")
        screen = QtWidgets.QDesktopWidget().screenGeometry()
        MainWindow.resize(screen.width()//2, screen.height()//2)
        MainWindow.setStyleSheet("#MainWindow{background-color: rgb(230, 230, 230)}")
        MainWindow.setWindowOpacity(1.0)
        MainWindow.setWindowTitle("PyMiner")
        MainWindow.setWindowIcon(QtGui.QIcon('images/main.ico'))
        MainWindow.setIconSize(QtCore.QSize(30, 30))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setBold(True)
        font.setPointSize(9)
        MainWindow.setFont(font)

        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        MainWindow.setCentralWidget(self.centralwidget)


        self.tabwidget_main = QtWidgets.QTabWidget(self.centralwidget)
        # self.tabwidget_main.setGeometry(QtCore.QRect(10, 10, 400, 600))
        self.tabwidget_main.setObjectName("tabwidget_main")
        self.tabwidget_main.setTabPosition(QtWidgets.QTabWidget.North)
        # This enum type defines the shape of the tabs:
        # Constant                  Description
        # QTabWidget.Rounded        The tabs are drawn with a rounded look. This is the default shape.
        # QTabWidget.Triangular     The tabs are drawn with a triangular look.
        self.tabwidget_main.setTabShape(QtWidgets.QTabWidget.Triangular)
        # This enum specifies where the ellipsis should appear when displaying texts that don’t fit:
        # Constant          Description
        # Qt.ElideLeft      The ellipsis should appear at the beginning of the text.
        # Qt.ElideRight     The ellipsis should appear at the end of the text.
        # Qt.ElideMiddle    The ellipsis should appear in the middle of the text.
        # Qt.ElideNone      Ellipsis should NOT appear in the text.
        self.tabwidget_main.setElideMode(QtCore.Qt.ElideRight)
        self.tabwidget_main.setCurrentIndex(0)


        # Pathway Searcher
        self.win_ps = QtWidgets.QWidget()
        # self.win_ps.setFont(QtGui.QFont("Microsoft YaHei", 6, QtGui.QFont.Bold))
        self.tabwidget_main.addTab(self.win_ps, QtGui.QIcon('images/pathway_searcher.ico'), "Pathway Searcher")
        self.pSearchUi = Pathway_Seracher_Widgets()
        self.win_ps.setLayout(self.pSearchUi.win_ps_hboxlayout)


        """
        #  Pathway Designer
        self.win_pd = QtWidgets.QWidget()
        # self.win_pd.setFont(QtGui.QFont("Microsoft YaHei", 6, QtGui.QFont.Bold))
        self.tabwidget_main.addTab(self.win_pd, QtGui.QIcon('images/Pathway Designer.ico'), "Pathway Designer")
        # """

        """
        # Enzyme Miner
        self.win_em = QtWidgets.QWidget()
        # self.win_em.setFont(QtGui.QFont("Microsoft YaHei", 6, QtGui.QFont.Bold))
        self.tabwidget_main.addTab(self.win_em, QtGui.QIcon('images/Enzyme Miner'), "Enzyme Miner")
        # """

        # """
        # statusbar
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        # self.statusBar.setStyleSheet("background-color: rgb(220,220,220)")
        self.statusBar.setStyleSheet('border-width: 0px;border-style: solid;border-color: rgb(210, 210, 210);background-color: rgb(220, 220, 220);')
        self.statusBar.setObjectName("statusBar")
        self.statusBar.setFont(QtGui.QFont("Times New Roman", 9, QtGui.QFont.Thin))
        MainWindow.setStatusBar(self.statusBar)
        self.statusBar.showMessage('Welcome to PyMiner . . . ',)
        #

        # date, time, memory in statusbar
        self.status_widget = Group_StatusBar_Widgets()
        self.statusBar.addPermanentWidget(self.status_widget.label_status_date)
        self.statusBar.addPermanentWidget(self.status_widget.label_status_time)
        self.statusBar.addPermanentWidget(self.status_widget.label_status_memory)
        #

        # timer
        self.timer = QtCore.QTimer(MainWindow)
        self.timer.timeout.connect(self.status_widget.showStatus)
        self.timer.start(1000)
        # """

        self.gridlayout_main =QtWidgets.QGridLayout()
        self.gridlayout_main.addWidget(self.tabwidget_main, 0, 0)
        self.centralwidget.setLayout(self.gridlayout_main)
        # QtCore.QMetaObject.connectSlotsByName(MainWindow)


    @classmethod
    def demoFunc(cls):
        """ """
        pass


if __name__ == "__main__":
    """ """
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    # # MainWindow.show()
    MainWindow.showMaximized()
    sys.exit(app.exec_())





