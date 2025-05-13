#!/usr/bin/env python

from PyQt5.QtCore import QDateTime, Qt, QTimer
from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget)


class RAWRApp(QDialog):
    def __init__(self, parent=None):
        super(RAWRApp, self).__init__(parent)

        self.originalPalette = QApplication.palette()


        QApplication.setStyle(QStyleFactory.create("Fusion"))

        self.createTopLeftGroupBox()
        self.createTopRightGroupBox()
        self.createBottomLeftGroupBox()
        self.createBottomRightGroupBox()

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.topLeftGroupBox, 1, 0)
        mainLayout.addWidget(self.topRightGroupBox, 1, 1)
        mainLayout.addWidget(self.bottomLeftGroupBox, 2, 0)
        mainLayout.addWidget(self.bottomRightGroupBox, 2, 1)
        mainLayout.setRowStretch(1, 1)
        mainLayout.setRowStretch(2, 1)
        mainLayout.setColumnStretch(0, 1)
        mainLayout.setColumnStretch(1, 1)
        self.setLayout(mainLayout)

        self.setWindowTitle("RAWR and SERES sequence resampling")


    def createTopLeftGroupBox(self):
        self.topLeftGroupBox = QGroupBox("Resampling Algorithm")

        radioButton1 = QRadioButton("RAWR")
        radioButton2 = QRadioButton("SERES")
        radioButton1.setChecked(True)

        layout = QVBoxLayout()
        layout.addWidget(radioButton1)
        layout.addWidget(radioButton2)
        layout.addStretch(1)
        self.topLeftGroupBox.setLayout(layout)

    def createTopRightGroupBox(self):
        self.topRightGroupBox = QGroupBox("Support Type")

        radioButton1 = QRadioButton("Multiple sequence alignment")
        radioButton2 = QRadioButton("Phylogenetic tree")
        radioButton1.setChecked(True)

        layout = QVBoxLayout()
        layout.addWidget(radioButton1)
        layout.addWidget(radioButton2)
        layout.addStretch(1)
        self.topRightGroupBox.setLayout(layout)

    def createBottomLeftGroupBox(self):
        self.bottomLeftGroupBox = QGroupBox("Algorithm Parameters")

        spinBox = QSpinBox(self.bottomLeftGroupBox)
        spinBox.setValue(20)
        spinBoxLabel = QLabel("&Sample Number:")
        spinBoxLabel.setBuddy(spinBox)

        lineEdit = QLineEdit('0.1')
        lineEditLabel = QLabel("&Reverse Rate:")
        lineEditLabel.setBuddy(lineEdit)

        layout = QVBoxLayout()
        layout.addWidget(spinBoxLabel)
        layout.addWidget(spinBox)
        layout.addWidget(lineEditLabel)
        layout.addWidget(lineEdit)

        layout.addStretch(1)
        self.bottomLeftGroupBox.setLayout(layout)


    def createBottomRightGroupBox(self):
        self.bottomRightGroupBox = QGroupBox("Group 3")

        defaultPushButton = QPushButton("Select alignment file")
        defaultPushButton.setDefault(True)

        textEdit = QTextEdit()
        textEdit.setPlainText("Put the phylogenetic tree in\n"
                              "Newick format here.\n")
        textEditLabel = QLabel("&Input Tree:")
        textEditLabel.setBuddy(textEdit)

        layout = QGridLayout()
        layout.addWidget(defaultPushButton)
        layout.addWidget(textEditLabel)
        layout.addWidget(textEdit)
        self.bottomRightGroupBox.setLayout(layout)


# if __name__ == '__main__':
#
#     import sys
#
#     app = QApplication(sys.argv)
#     rawrApp = RAWRApp()
#     rawrApp.show()
#     sys.exit(app.exec_())