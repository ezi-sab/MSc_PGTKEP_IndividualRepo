# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui2.ui'
#
# Created: Mon Feb 13 15:12:14 2012
#      by: PyQt4 UI code generator 4.8.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore
from PyQt5.QtWidgets import QLabel, QPushButton, QApplication

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Form2(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(350, 118)
        self.label = QLabel(Form)
        self.label.setGeometry(QtCore.QRect(9, 9, 350, 31))
        self.label.setObjectName(_fromUtf8("label"))
        self.pushButton = QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(140, 50, 51, 23))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QApplication.translate("Form", "Form", None))
        self.label.setText(QApplication.translate("Form", "Click below to save the results in an excel file", None))
        self.pushButton.setText(QApplication.translate("Form", "Save", None))

