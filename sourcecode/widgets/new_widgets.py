# -*- coding: utf-8 -*-
"""
Created on Fri May 28 16:43:06 2021

@author: CC-SXF
"""


from PyQt5 import QtCore, QtGui, QtWidgets


class NewEventLineEdit(QtWidgets.QLineEdit):
    """
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    newKeyPressEvent = QtCore.pyqtSignal(str)
    #
    def __init__(self, *args, **kwargs):
        """ """
        super().__init__(*args, **kwargs)

    def keyPressEvent(self, event):
        # # QtWidgets.QLineEdit.keyPressEvent(self, event)
        super().keyPressEvent(event)
        if event.key() == QtCore.Qt.Key_Return:
            text = self.text()
            # # print(QtCore.QDateTime.currentDateTime().toString())
            # # print(event.key())
            # # print("@Event@:", text)
            self.newKeyPressEvent.emit(text)

    @classmethod
    def demoFunc(cls):
        """ """
        pass


class MultCompleter(QtWidgets.QCompleter):
    """
    """
    def __init__(self, *args, **kwargs):
        """ """
        super().__init__(*args, **kwargs)

    # Overriding
    def splitPath(self, in_string):
        """
        """
        # splitPath(self, str) -> List[str]
        in_string_list = in_string.split('; ')
        in_string = in_string_list[-1]
        if in_string == "":
            in_string = "@pyminer"
        return [in_string]

    @classmethod
    def demoFunc(cls):
        """ """
        pass


class NewAppendTextBrowser(QtWidgets.QTextBrowser):
    """
    """
    def __init__(self, *args, **kwargs):
        """ """
        super().__init__(*args, **kwargs)

    def append(self, string):
        """ """
        # scroll down
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
        # scroll left
        self.horizontalScrollBar().setValue(self.horizontalScrollBar().minimum())
        #
        string = "".join(["<span>", string, "</span>"])
        super().append(string)

    @classmethod
    def demoFunc(cls):
        """ """
        pass


class NewEventTableWidget(QtWidgets.QTableWidget):
    """
    """
    # Create one or more overloaded unbound signals as a class attribute.
    # PyQt5.QtCore.pyqtSignal(types[, name[, revision=0[, arguments=[]]]])
    newWheelEvent = QtCore.pyqtSignal(dict)
    #
    def __init__(self, *args, **kwargs):
        """ """
        super().__init__(*args, **kwargs)

    def wheelEvent(self, event):
        """ """
        # # super().wheelEvent(event)
        # print("@mousewheel:", event.angleDelta().y())
        y_delta = {"Zoom_In_out": event.angleDelta().y()//120}
        self.newWheelEvent.emit(y_delta)

    @classmethod
    def demoFunc(cls):
        """ """
        pass





