#!/usr/bin/python -d

import sys
from PyQt4 import QtCore, QtGui
from gui2 import Ui_Form2

class MyForm2(QtGui.QMainWindow):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self, parent)
    self.ui = Ui_Form2()
    self.ui.setupUi_nolayout(self)
    QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.filebrower('teststring'))


  def filebrower(content):
    fName = QtGui.QFileDialog.getSaveFileName(self, "Save as text file", "Save as new file", self.tr("Text Files (*.xls)"))
    if fName.isEmpty() == False:
        fptr = open(fName, 'w')
        fptr.write(content)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = QtGui.QWidget()
    Newui = MyForm2()
    Newui.ui.setupUi_layout(Form)
    Form.show()
    sys.exit(app.exec_())
