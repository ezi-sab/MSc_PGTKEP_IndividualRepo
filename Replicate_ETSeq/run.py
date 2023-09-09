#!/usr/bin/python -d

import sys
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QMessageBox, QFileDialog
from gui import Ui_Form
from gui2 import Ui_Form2
from gui3 import Ui_Form3
import seqdep
import time
from Bio import Entrez, SeqIO
from xlwt import Workbook
from xlrd import open_workbook
import xlrd
from xlutils.copy import copy
import os


class MyForm(QMainWindow):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.ui = Ui_Form()
        self.ui.setupUi_nolayout(self)
        self.teststring='MyForm1'
        self.proc_ok=0
    #QtCore.QObject.connect(self.ui.pushButton_3, QtCore.SIGNAL("clicked()"), self.filebrower) #brower
    #         
        self.ui.pushButton.clicked.connect(self.prediction)
        self.ui.pushButton.clicked.connect(self.submitnewWindow)
        self.ui.pushButton_4.clicked.connect(self.help)
        # QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.prediction) #predict both thermodynamics and 2 methods
        # QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.submitnewWindow) #submit
        # QtCore.QObject.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.help) #submit

    def filebrower(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', 'c:/')
        self.ui.lineEdit.setText(fname)

    def submitnewWindow(self):
        if self.proc_ok==1: # if the xls file is sucessfully generated!
            self.myOtherWindow = MyForm2()
            self.myOtherWindow.show()
            self.myOtherWindow.setWindowTitle("Save")


    def help(self):
        self.myOtherWindow = MyForm3()
        self.myOtherWindow.show()
        self.myOtherWindow.setWindowTitle("Help")


    def prediction(self):
        global file1,file2
        temp_filename="~"+str(int(time.time()))+'.txt'
        file1=temp_filename
        time.sleep(1)
        data=self.ui.plainTextEdit.toPlainText()
        data=str(data).strip()
        dic_therm={}
        dic_jf={}    
        
        t=seqdep.check_input(data, temp_filename)
        # the input is in right format and stored in temp_filename sucessfully!
        if (t!=1 and t!=2):
            QMessageBox.critical(self, "Wrong", "Something is wrong with your input!\n"+str(t))
            exit()    
            # now temp_filename is the file contain all the template sequence
        else:
            if (t==2):
                QMessageBox.warning(self, "Warning", "some templates' names were showed twice!\n those showed up later were ignored!")
            trig_filename=seqdep.trigen(temp_filename)
            if trig_filename ==0:
                QMessageBox.critical(self, "Wrong", "NO GACTC in template sequences")
                return 0
            file2=trig_filename
	#print "trigger_file:", trig_filename,"temp_file",temp_filename
            print ("doing thermodynacis calculating right now!")
            dic_therm=seqdep.unafold_pred(temp_filename,trig_filename)
            print("doing sequence prediction right now!")
            for cur_record_o in SeqIO.parse(open(temp_filename), "fasta") :
                cur_seq=str(cur_record_o.seq).upper() #all capitalized
                dic_jf[cur_record_o.id]=seqdep.method_2_prediction(cur_seq)
            print("Gathered all information!Please save your file!")
            #print "outside dic_jf", dic_jf
            #print "outside dic_therm", dic_therm
            # print(dic_therm)
            print(dic_jf)
            seqdep.write_xls(self,file1+'simple.xls',dic_jf, dic_therm)
            self.proc_ok=1
            #print "step4", globals()

    

class MyForm2(QMainWindow):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.ui = Ui_Form2()
        self.ui.setupUi(self)
        self.ui.pushButton.clicked.connect(self.savefilebrower_all)
        # QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.savefilebrower_all)

    def savefilebrower_all(self):
        fName = QFileDialog.getSaveFileName(self, "Save as xls file", "Save as new file", self.tr("Text Files (*.xls)"))
        #print "step5", globals()
        if fName.isEmpty() == False:
	    #self.write_xls(fName,self, dic_jf, dic_therm)
            t=file1+'simple.xls'
            rb = xlrd.open_workbook(t)
            wb = copy(rb)
            #print(rb.sheet_by_index(0).cell(0,0).value)
            wb.save(fName)
            print(fName)
            print('is generated!')
            QMessageBox.about(self, "Saved", "%s is generated!" % (fName))
            seqdep.delete_file(file1)
            seqdep.delete_file(file2)
            self.close()


class MyForm3(QMainWindow):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.ui = Ui_Form3()
        self.ui.setupUi(self)



if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    Form = QWidget()
    if os.system('perl -v >perl.version'):
        QMessageBox.critical(Form,"Error", "you haven't properly installed perl yet!")
        exit()
                    
    global file1 # prefix name for template file
    global file2 # prefix name for trigger file
    global seq_list # template name and sequence dictionary
    Newui = MyForm()
    Newui.ui.setupUi_layout(Form)
    Form.setWindowTitle("ETSeq")
    Form.show()
    sys.exit(app.exec_())
