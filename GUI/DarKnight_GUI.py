import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
import openbabel
import darkchem
import scipy.misc
import sys
sys.path.append('../figures')
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QTextEdit, QTextBrowser, QMessageBox, QFrame,
                             QGridLayout, QToolTip, QPushButton, QApplication, QProgressBar)
from PyQt5.QtGui import QFont, QPixmap, QIcon
from PyQt5.QtCore import QCoreApplication, QSize, QThread,pyqtSignal
import time
import tensorflow as tf
tf.logging.set_verbosity(tf.logging.ERROR)


def load_model():
    """
    load the model that converts smiles strings to the latent space vectors.
    """
    model = darkchem.utils.load_model('../models/N7b_[M+H]')
    return model


# the preprocess functions


def struc2mol(sms):
    """
    A function to transform smiles strings to molecules with the module
    rdkit.Chem.MolFromSmiles, and return a DataFrame
    """
    save = pd.DataFrame(columns=['raw_smiles', 'smiles', 'mol'])
    save['raw_smiles'] = sms['smiles']
    for i in range(sms.shape[0]):
        save['mol'][i] = Chem.MolFromSmiles(sms['smiles'][i])
        if save['mol'][i] is None:
            save['smiles'][i] = 'Invalid smi str'
        else:
            save['smiles'][i] = sms['smiles'][i]
    return save


def Standardize_SMI(smi):
    """
    A function used to standardize smile strings. For optimize prediction result purpose.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smi)
    outMDL = obConversion.WriteString(mol)[:-2]
    return outMDL


def tranform(smi, model, path_vec, k):
    """
    The intermediate function used to tranform reactant smile string to product smile string
    """
    test = darkchem.utils.struct2vec(smi)
    test = np.array(test)
    test = test.reshape(-1, 100)
    t_l = model.encoder.predict(test)
    t_pre = t_l + path_vec
    t_pred = model.decoder.predict(t_pre)
    trs = darkchem.utils.beamsearch(t_pred, k=k)
    trs = trs.reshape(-1, 100)
    v2s = [darkchem.utils.vec2struct(trs[i]) for i in range(len(trs))]
    std = [Standardize_SMI(v2s[i]) for i in range(len(v2s))]
    return std


class Runthread(QThread):
    """"
    The class is used to calculate and produce the prediction result
    """
    _signal = pyqtSignal(str)
    # trigger = pyqtSignal(int)

    def __init__(self,sms,clear,pathvec):
        super(Runthread, self).__init__()
        self.sms = sms
        self.clear = clear
        self.path_vec = pathvec
    def __del__(self):
        self.wait()

    def run(self,k=15):
        #self.clear
        smi = self.sms
        model = load_model()
        # path_vec = load_path_vector()
        path_vec = self.path_vec
        a = ['Reactant', 'Product']
        b = []
        c = [smi]
        std = tranform(smi, model, path_vec, k)
        for j in range(len(std)):
            if std[j] == smi.upper():
                prd = std[j]
                break
            elif smi.replace('#', '') == std[j]:
                prd = std[j]
                break
            else:
                prd = std[14]
        b.append(prd)
        c.append(prd)
        out = pd.DataFrame(data=b, columns=['Product'])
        out.insert(0, 'Reactant', smi)
        df = struc2mol(pd.DataFrame(data=c, columns=['smiles']))
        df.insert(3, 'legend', a)
        out = (PandasTools.FrameToGridImage(df, column='mol', legendsCol='smiles', molsPerRow=2))
        a = np.array(out)
        scipy.misc.imsave('outfile.png', a)
        figname = 'outfile.png'
        #time.sleep(3)
        self._signal.emit(str(figname))
        # self.trigger.emit()


class ShowProgress(QThread):
    """"
    This class is designed for the progressbar
    """

    progressBarValue = pyqtSignal(int)

    def __init__(self):
        super(ShowProgress, self).__init__()


    def run(self):
        for i in range(101):
            time.sleep(0.07)
            self.progressBarValue.emit(i)
      
        
class Prediction(QWidget):
    def __init__(self):
        super(Prediction, self).__init__()

        self.setWindowIcon(QIcon('../figures/darknight_icon.ico'))

        self.progressBar = QProgressBar(self)
        self.progressBar.setGeometry(160, 280, 440, 20)

        self.input = QLineEdit(self)
        self.input.move(10,65)

        self.pathmd = QLineEdit(self)
        self.pathmd.move(10,135)

        self.reactant = QLabel('Reactant:',self)
        self.reactant.move(10,40)

        self.result = QLabel('Prediction Result:', self)
        self.result.move(160, 40)

        self.rct_type = QLabel('Reaction_Type:',self)
        self.rct_type.move(10,110)

        self.resize(600, 450)
        self.setWindowTitle("DarKnight")

        self.plabel = QLabel(self)
        self.plabel.setFixedSize(400, 200)
        self.plabel.move(160, 60)

        self.logolabel = QLabel(self)
        self.logolabel.setFixedSize(125, 125)
        self.logolabel.move(280, 310)


        self.plabel.setStyleSheet("QLabel{background:white;}"
                                 "QLabel{color:rgb(300,300,300,120);font-size:10px;font-weight:bold;font-family:YaHei;}"
                                 )

        self.logolabel.setStyleSheet(
                                  "QLabel{color:rgb(300,300,300,120);font-size:10px;font-weight:bold;font-family:YaHei;}"
                                  )


        btn = QPushButton(self)
        btn.setText("Predict")
        btn.move(10, 200)
        btn.clicked.connect(self.load_path)  # predict button connect load_model function
        btn.clicked.connect(self.start_login)
        btn.clicked.connect(self.progressbar)

        cbtn = QPushButton(self)
        cbtn.setText("Clear")
        cbtn.move(10,250)
        cbtn.clicked.connect(self.clear_smile)
        cbtn.clicked.connect(self.clear_output)


        qbtn = QPushButton(self)
        qbtn.setText("Quit")
        qbtn.move(10, 300)
        qbtn.clicked.connect(QCoreApplication.instance().quit)

        QToolTip.setFont(QFont('SansSerif', 10))
        btn.setToolTip('<b>Predict the Product</b> of input reactant')
        cbtn.setToolTip('<b>Clear</b> the input')
        qbtn.setToolTip('<b>Quit</b> the execution instantly')

        logo = QPixmap('../figures/uw_purple.jpeg').scaled(self.logolabel.width(), self.logolabel.height())
        self.logolabel.setPixmap(logo)


    def get_smi(self):
        """
        Accuqire the smile string
        """
        smi = self.input.text()
        return smi


    def clear_smile(self):
        """"
        Clear the input of smile string
        """
        self.input.clear()
        self.pathmd.clear()


    def clear_output(self):  # problem
        """"
        clear the output
        """
        self.plabel.clear()

    def load_path(self):
        """"
        load the path vector for the predicted type of chemical reaction
        """
        path = self.pathmd.text()
        path_vec = np.load(path)
        return path_vec
          
    def start_login(self):
        self.thread = Runthread(sms = self.get_smi(), clear = self.clear_output(),pathvec = self.load_path()) # add one parameter in Qthread class
        self.thread._signal.connect(self.load_graph) 
        self.thread.start()

    def load_graph(self,fig):
        jpg = QPixmap(fig).scaled(self.plabel.width(), self.plabel.height())
        self.plabel.setPixmap(jpg)

    def progressbar(self):
        self.thread_1 = ShowProgress()
        self.thread_1.progressBarValue.connect(self.copy_file)
        self.thread_1.start()

    def copy_file(self, i):
        self.progressBar.setValue(i)


    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to exit?", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    my = Prediction()
    my.show()
    sys.exit(app.exec_())


