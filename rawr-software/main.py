# coding:utf-8
# Available under GPLv3 license
import sys
import os
import subprocess
import mainWindow
import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from src import seqs, generateSamples, calSupport, generateSamples_parallel, calSupport_parallel, MSA_support_csv_2_jalview_sequence_annotation
import ete3
import signal 

class WorkThread(QThread):
    trigger = pyqtSignal(int)

    def __init__(self, alnData=None, inputTree=None, multiprocess=-1, parameters=None):
        super(WorkThread, self).__init__()
        self.alnData = alnData
        self.inputTree = inputTree
        self.multiprocess = multiprocess
        self.parameters = parameters
        self.currProgress = 0

    def run(self):
        if self.multiprocess < 0:
            self.currProgress = generateSamples.generateSampleSeq(
                self.alnData, self.parameters, self.trigger, self.currProgress)
            self.currProgress = calSupport.estimateSampleAln(
                self.parameters, self.trigger, self.currProgress)
            if self.parameters["supportType"] == "MSA":
                support = calSupport.calculateMSASupport(
                    self.alnData, self.parameters, self.trigger, self.currProgress)
                supportFile = os.path.join(self.parameters["outputDir"],"MSA_Support.csv")
                calSupport.writeSupport(self.alnData, support, supportFile)
                MSA_support_csv_2_jalview_sequence_annotation.support_csv_2_jalview(supportFile,"pink")
            elif self.parameters["supportType"] == "Tree":
                calSupport.estimateSampleTree(self.parameters, self.trigger, self.currProgress)
                calSupport.calTreeSupport(self.inputTree, self.parameters)
            self.trigger.emit(100)
        else:
            self.currProgress = generateSamples_parallel.generateSampleSeq(self.alnData, self.parameters, self.multiprocess, self.trigger, self.currProgress)
            self.currProgress = calSupport_parallel.estimateSampleAln(self.parameters,self.multiprocess, self.trigger, self.currProgress)
            if self.parameters["supportType"] == "MSA":
                support = calSupport_parallel.calculateMSASupport(self.alnData, self.parameters, self.multiprocess, self.trigger, self.currProgress)
                supportFile = os.path.join(self.parameters["outputDir"],"MSA_Support.csv")
                calSupport_parallel.writeSupport(self.alnData, support, supportFile)   
                MSA_support_csv_2_jalview_sequence_annotation.support_csv_2_jalview(supportFile,"pink")  
            elif self.parameters["supportType"] == "Tree":
                calSupport_parallel.estimateSampleTree(self.parameters, self.multiprocess, self.trigger, self.currProgress)
                calSupport_parallel.calTreeSupport(self.inputTree, self.parameters)
            self.trigger.emit(100)


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(QMainWindow, self).__init__(parent)
        self.ui = mainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.rawrRadioButton.toggled['bool'].connect(
            self.ui.anchorLengthSpinBox.setDisabled)
        self.ui.rawrRadioButton.toggled['bool'].connect(
            self.ui.anchorNumberSpinBox.setDisabled)
        self.ui.rawrRadioButton.toggled['bool'].connect(
            self.ui.anchorLengthLabel.setDisabled)
        self.ui.rawrRadioButton.toggled['bool'].connect(
            self.ui.anchorNumberLabel.setDisabled)
        self.ui.seresRadioButton.toggled['bool'].connect(
            self.ui.anchorLengthSpinBox.setEnabled)
        self.ui.seresRadioButton.toggled['bool'].connect(
            self.ui.anchorNumberSpinBox.setEnabled)
        self.ui.seresRadioButton.toggled['bool'].connect(
            self.ui.anchorLengthLabel.setEnabled)
        self.ui.seresRadioButton.toggled['bool'].connect(
            self.ui.anchorNumberLabel.setEnabled)
        self.ui.msaRadioButton.toggled['bool'].connect(
            self.ui.treeFileLabel.setDisabled)
        self.ui.treeRadioButton.toggled['bool'].connect(
            self.ui.treeFilePushButton.setEnabled)
        self.parameters = {}
        self.ui.msaFilePushButton.clicked.connect(self.openMSAFile)
        self.ui.treeFilePushButton.clicked.connect(self.openTreeFile)
        self.ui.outputPathPushButton.clicked.connect(self.openFolder)
        self.ui.mpNumberSpinBox.setDisabled(True)
        self.ui.mpPushButton.toggled['bool'].connect(self.ui.mpNumberSpinBox.setEnabled)
        self.ui.buttonBox.button(
            QDialogButtonBox.Reset).clicked.connect(
            self.resetParams)
        self.ui.buttonBox.accepted.connect(self.execute)

    def execute(self):
        self.getParams()
        if self.validParams():
            self.work = WorkThread(self.alnData, self.inputTree, self.multiprocess, self.parameters)
            self.work.start()
            self.work.trigger.connect(self.display)
            self.work.finished.connect(self.showFinishSignal)

    def showFinishSignal(self):
        if self.parameters["supportType"] == "MSA":
            self.finishMessageBox("Support estimation finished. \nMSA support is written to " + os.path.join(self.parameters["outputDir"], "MSA_Support.csv"))
        elif self.parameters["supportType"] == "Tree":
            ts = ete3.TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            treeFile = os.path.join(self.parameters["outputDir"], "tree.support.txt")
            treeFigFile = os.path.join(self.parameters["outputDir"], "tree.support.png")
            if os.path.exists(treeFile):
                with open(treeFile, "r") as inf:
                    try:
                        tree = ete3.Tree(inf.readline())
                        tree.render(treeFigFile, tree_style=ts)
                        tree.show(tree_style=ts)
                        print("Tree rendering complete.")
                        self.finishMessageBox("Support estimation finished. \nTree support is written to " + os.path.join(self.parameters["outputDir"], "tree.support in both PNG and TXT"))                        
                    except BaseException:
                        self.errorMessageBox(
                            "No valid tree found in " + treeFile + ".")
        self.ui.progressBar.setValue(0)

    def display(self, value):
        self.ui.progressBar.setValue(value)

    def validParams(self):
        if "algorithm" not in self.parameters or self.parameters["algorithm"] not in [
                "RAWR", "SERES"]:
            self.errorMessageBox("Please select an algorithm, RAWR or SERES.")
            return False
        if "supportType" not in self.parameters or self.parameters["supportType"] not in [
                "MSA", "Tree"]:
            self.errorMessageBox(
                "Please select a support type, Multiple Sequence Alignment or Phylogenetic Tree.")
            return False
        if self.parameters["sampleNum"] < 1 or self.parameters["sampleNum"] > 500:
            self.errorMessageBox(
                "Please set the sample number between 1 and 500.")
            return False
        if self.multiprocess >= 0:
            if self.multiprocess < 2 or self.multiprocess > 50:
                self.errorMessageBox(
                    "Please set the number of processes for multiprocessing between 2 and 50.")
                return False
        if "alnFile" not in self.parameters or not os.path.exists(
                self.parameters["alnFile"]):
            self.errorMessageBox(
                "Please select a multiple sequence alignment file in *.fasta, *.fas, or *.fa format.")
            return False
        if "outputDir" not in self.parameters or not os.path.exists(
                self.parameters["outputDir"]):
            self.errorMessageBox(
                "Please select an output folder.")
            return False

        try:
            self.alnData = seqs.getAlnData(self.parameters["alnFile"])
        except BaseException:
            self.errorMessageBox(
                "Please select a valid multiple sequence alignment file in .fasta format.")
            return False

        if self.parameters["algorithm"] == "SERES" and self.parameters["anchorLen"] * \
                self.parameters["anchorNum"] >= self.alnData.shape[1]:
            self.errorMessageBox("Anchor length exceed alignment length.")
            return False

        if self.parameters["supportType"] == "Tree":
            try:
                if self.parameters["inputTree"]:
                    self.inputTree = ete3.Tree(self.parameters["inputTree"])
                    leaveNodes = [leaf.name for leaf in self.inputTree]
                if set(leaveNodes) != set(self.alnData.index.to_list()):
                    self.errorMessageBox(
                      "Input tree does not match input alignment.")
            except BaseException:
                self.errorMessageBox(
                    "Please select a valid input tree in newick format.")
                return False
        else:
            self.inputTree = None
        return True

    def errorMessageBox(self, msgText):
        QMessageBox.critical(
            self,
            "Error in Parameters",
            msgText,
            QMessageBox.Ok)

    def finishMessageBox(self, msgText):
        QMessageBox.information(
            self,
            "Support estimation finished!",
            msgText,
            QMessageBox.Ok)

    def getParams(self):
        print("Get params.")
        if self.ui.rawrRadioButton.isChecked():
            self.parameters["algorithm"] = "RAWR"
        elif self.ui.seresRadioButton.isChecked():
            self.parameters["algorithm"] = "SERES"
            self.parameters["anchorNum"] = int(
                self.ui.anchorNumberSpinBox.value())
            self.parameters["anchorLen"] = int(
                self.ui.anchorLengthSpinBox.value())
        self.parameters["reverseRate"] = float(
            self.ui.reverseRateSpinBox.value())
        self.parameters["sampleNum"] = int(self.ui.sampleNumberSpinBox.value())

        if self.ui.msaRadioButton.isChecked():
            self.parameters["supportType"] = "MSA"
        elif self.ui.treeRadioButton.isChecked():
            self.parameters["supportType"] = "Tree"
        if self.ui.mpPushButton.isChecked():
            self.multiprocess = int(self.ui.mpNumberSpinBox.value())
        else:
            self.multiprocess = -1

        print(self.parameters)

    def resetParams(self):
        print("Reset params.")
        self.parameters = {}
        self.ui.rawrRadioButton.setChecked(True)
        self.ui.treeRadioButton.setChecked(True)
        self.ui.mpPushButton.setChecked(False)
        self.ui.sampleNumberSpinBox.setValue(10)
        self.ui.reverseRateSpinBox.setValue(0.1)
        self.ui.anchorNumberSpinBox.setValue(20)
        self.ui.anchorLengthSpinBox.setValue(5)

    def openMSAFile(self):
        fileName, fileType = QFileDialog.getOpenFileName(
            None, 'Select MSA File', '', 'MSA file (*.fasta *.fa *.fas);; All Files (*.*)')
        self.parameters["alnFile"] = os.path.normpath(fileName)

    def openTreeFile(self):
        fileName, fileType = QFileDialog.getOpenFileName(
            None, 'Select Tree File', '', 'Tree file (*.tree *.newick);; All Files (*.*)')
        self.parameters["inputTree"] = os.path.normpath(fileName)

    def openFolder(self):
        outputDir = os.path.normpath(QFileDialog.getExistingDirectory(self, "Select output folder"))
        self.parameters["outputDir"] = outputDir

    def closeEvent(self, event):
        event.accept()  # Proceed with closing the main window
if __name__ == '__main__':
    myapp = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(myapp.exec_())
