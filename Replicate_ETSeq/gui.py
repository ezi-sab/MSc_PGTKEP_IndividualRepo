# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'EXPAR_paper.ui'
#
# Created: Fri Feb 03 23:08:13 2012
#      by: PyQt4 UI code generator 4.8.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore
from PyQt5.QtWidgets import QGridLayout, QLabel, QPlainTextEdit, QLineEdit,QPushButton, QApplication

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(476, 230)
        self.gridLayout = QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label = QLabel(Form)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.plainTextEdit = QPlainTextEdit(Form)
        self.plainTextEdit.setObjectName(_fromUtf8("plainTextEdit"))
        self.gridLayout.addWidget(self.plainTextEdit, 1, 0, 1, 4)
        #self.label_2 = QLabel(Form)
        #self.label_2.setObjectName(_fromUtf8("label_2"))
        #self.gridLayout.addWidget(self.label_2, 2, 0, 1, 2)
        #self.lineEdit = QLineEdit(Form)
        #self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        #self.gridLayout.addWidget(self.lineEdit, 3, 0, 1, 3)
        #self.pushButton_3 = QPushButton(Form)
        #self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        #self.gridLayout.addWidget(self.pushButton_3, 3, 3, 1, 1)
        self.label_3 = QLabel(Form)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 4, 0, 1, 1)
        self.lineEdit_2 = QLineEdit(Form)
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))
        self.gridLayout.addWidget(self.lineEdit_2, 4, 1, 1, 1)
        self.label_4 = QLabel(Form)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 4, 2, 1, 1)
        self.lineEdit_3 = QLineEdit(Form)
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.gridLayout.addWidget(self.lineEdit_3, 4, 3, 1, 1)
        self.label_5 = QLabel(Form)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 5, 0, 1, 1)
        self.lineEdit_4 = QLineEdit(Form)
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.gridLayout.addWidget(self.lineEdit_4, 5, 1, 1, 1)
        self.label_6 = QLabel(Form)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout.addWidget(self.label_6, 5, 2, 1, 1)
        self.lineEdit_5 = QLineEdit(Form)
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.gridLayout.addWidget(self.lineEdit_5, 5, 3, 1, 1)
        self.pushButton = QPushButton(Form)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.gridLayout.addWidget(self.pushButton, 6, 0, 1, 2)
        self.pushButton_4 = QPushButton(Form)
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.gridLayout.addWidget(self.pushButton_4, 6, 3, 1, 1)
        self.pushButton_2 = QPushButton(Form)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.gridLayout.addWidget(self.pushButton_2, 6, 2, 1, 1)

    def setupUi_layout(self, Form):
        self.gridLayout = QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.gridLayout.addWidget(self.plainTextEdit, 1, 0, 1, 4)
        #self.gridLayout.addWidget(self.label_2, 2, 0, 1, 2)
        #self.gridLayout.addWidget(self.lineEdit, 3, 0, 1, 3)
        #self.gridLayout.addWidget(self.pushButton_3, 3, 3, 1, 1)
        self.gridLayout.addWidget(self.label_3, 4, 0, 1, 1)
        self.gridLayout.addWidget(self.lineEdit_2, 4, 1, 1, 1)
        self.gridLayout.addWidget(self.label_4, 4, 2, 1, 1)
        self.gridLayout.addWidget(self.lineEdit_3, 4, 3, 1, 1)
        self.gridLayout.addWidget(self.label_5, 5, 0, 1, 1)
        self.gridLayout.addWidget(self.lineEdit_4, 5, 1, 1, 1)
        self.gridLayout.addWidget(self.label_6, 5, 2, 1, 1)
        self.gridLayout.addWidget(self.lineEdit_5, 5, 3, 1, 1)
        self.gridLayout.addWidget(self.pushButton, 6, 0, 1, 2)
        self.gridLayout.addWidget(self.pushButton_4, 6, 3, 1, 1)
        self.gridLayout.addWidget(self.pushButton_2, 6, 2, 1, 1)

    def setupUi_nolayout(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(476, 230)
        self.label = QLabel(Form)
        self.label.setObjectName(_fromUtf8("label"))
        self.plainTextEdit = QPlainTextEdit(Form)
        self.plainTextEdit.setObjectName(_fromUtf8("plainTextEdit"))
        #self.label_2 = QLabel(Form)
        #self.label_2.setObjectName(_fromUtf8("label_2"))
        #self.lineEdit = QLineEdit(Form)
        #self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        #self.pushButton_3 = QPushButton(Form)
        #self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        self.label_3 = QLabel(Form)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.lineEdit_2 = QLineEdit(Form)
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))
        self.label_4 = QLabel(Form)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.lineEdit_3 = QLineEdit(Form)
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.label_5 = QLabel(Form)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.lineEdit_4 = QLineEdit(Form)
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.label_6 = QLabel(Form)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.lineEdit_5 = QLineEdit(Form)
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.pushButton = QPushButton(Form)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton_4 = QPushButton(Form)
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.pushButton_2 = QPushButton(Form)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))

        self.retranslateUi(Form)
        #QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.lineEdit.clear)
        self.pushButton_2.clicked.connect(self.plainTextEdit.clear)
        self.pushButton_2.clicked.connect(self.lineEdit_5.undo)
        self.pushButton_2.clicked.connect(self.lineEdit_4.undo)
        self.pushButton_2.clicked.connect(self.lineEdit_2.undo)
        self.pushButton_2.clicked.connect(self.lineEdit_3.undo)
        # QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.plainTextEdit.clear)
        # QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.lineEdit_5.undo)
        # QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.lineEdit_4.undo)
        # QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.lineEdit_2.undo)
        # QtCore.QObject.connect(self.pushButton_2, QtCore.SIGNAL(_fromUtf8("clicked()")), self.lineEdit_3.undo)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QApplication.translate("Form", "Form", None))
        self.label.setText(QApplication.translate("Form", "Input your template sequences in fasta format or \'name sequence\' format", None))
        #self.label_2.setText(QApplication.translate("Form", "Or upload your template file in fasta format", None))
        #self.pushButton_3.setText(QApplication.translate("Form", "Brower", None))
        self.label_3.setText(QApplication.translate("Form", "Minimum trigger-template Tm", None))
        self.lineEdit_2.setText(QApplication.translate("Form", "40", None))
        self.label_4.setText(QApplication.translate("Form", "Maximum trigger-template Tm", None))
        self.lineEdit_3.setText(QApplication.translate("Form", "60", None))
        self.label_5.setText(QApplication.translate("Form", "Maximum template-template Tm", None))
        self.lineEdit_4.setText(QApplication.translate("Form", "25", None))
        self.label_6.setText(QApplication.translate("Form", "Maximum template self bonds", None))
        self.lineEdit_5.setText(QApplication.translate("Form", "10", None))
        self.pushButton.setText(QApplication.translate("Form", "Submit", None))
        self.pushButton_4.setText(QApplication.translate("Form", "Help", None))
        self.pushButton_2.setText(QApplication.translate("Form", "Cancel", None))