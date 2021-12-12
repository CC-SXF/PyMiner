# -*- coding: utf-8 -*-
"""
Created on Thu May 13 13:46:36 2021

@author: CC-SXF
"""

# import math
import psutil
from PyQt5 import QtCore, QtGui, QtWidgets


class Group_StatusBar_Widgets():
    """ """

    def __init__(self):
        """ """
        # self.label_status_date = QtWidgets.QLabel()
        # self.label_status_time = QtWidgets.QLabel()
        # self.label_status_memory = QtWidgets.QLabel()
        self.initWidgets()

    def initWidgets(self):
        """ """
        self.label_status_date = QtWidgets.QLabel()
        self.label_status_date.setAlignment(QtCore.Qt.AlignCenter)
        self.label_status_date.setFrameShape(QtWidgets.QFrame.Box)
        self.label_status_date.setFont(QtGui.QFont("Times New Roman", 9, QtGui.QFont.Thin))
        self.label_status_date.setStyleSheet('border-width: 5px;border-style: solid;border-color: rgb(210, 210, 210);background-color: rgb(215,215,215);')
        # Raised, Sunken, Plain
        self.label_status_date.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_status_date.setText(' Date: 0000-00-00 ')
        self.label_status_date.setObjectName("label_status_date")

        self.label_status_time = QtWidgets.QLabel()
        self.label_status_time.setAlignment(QtCore.Qt.AlignCenter)
        self.label_status_time.setFrameShape(QtWidgets.QFrame.Box)
        self.label_status_time.setFont(QtGui.QFont("Times New Roman", 9, QtGui.QFont.Thin))
        self.label_status_time.setStyleSheet('border-width: 5px;border-style: solid;border-color: rgb(210, 210, 210);background-color: rgb(215,215,215);')
        # Raised, Sunken, Plain
        self.label_status_time.setFrameShadow(QtWidgets.QFrame.Sunken)
        # # self.label_status_time.setText(' Time: 00:00:00 ')
        self.label_status_time.setText(' Time: 00:00 ')
        self.label_status_time.setObjectName("label_status_time")

        self.label_status_memory = QtWidgets.QLabel()
        self.label_status_memory.setAlignment(QtCore.Qt.AlignCenter)
        self.label_status_memory.setFrameShape(QtWidgets.QFrame.Box)
        self.label_status_memory.setFont(QtGui.QFont("Times New Roman", 9, QtGui.QFont.Thin))
        self.label_status_memory.setStyleSheet('border-width: 5px;border-style: solid;border-color: rgb(210, 210, 210);background-color: rgb(215,215,215);')
        # Raised, Sunken, Plain
        self.label_status_memory.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_status_memory.setText(' Memory: 00.0 % ')
        self.label_status_memory.setObjectName("label_status_memory")


    def showStatus(self):
        """ """
        # Date & Time
        date_time_new = QtCore.QDateTime.currentDateTime()
        date_time_new = date_time_new.toString("yyyy-MM-dd HH:mm:ss")
        (date_new, time_new) = date_time_new.split()
        date_new = "".join([" Date: ", date_new, " "])
        time_new = "".join([" Time: ", time_new, " "])
        # Memory
        memory_percent_new = psutil.virtual_memory().percent
        memory_percent_new = str(memory_percent_new)
        memory_percent_new = "".join([" Memory: ",  memory_percent_new, " % "])
        # update
        self.label_status_date.setText(date_new)
        self.label_status_time.setText(time_new)
        self.label_status_memory.setText(memory_percent_new)

    @classmethod
    def demoFunc(cls):
        """ """
        pass





