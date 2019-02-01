import sklearn
import scipy.io.matlab as sio
import numpy as np
from matplotlib.pyplot import plot
from sklearn.linear_model import Lasso
import keras

fname = r'L:\Users\guru\Documents\hartmann_lab\data\VG3D\pilot\sample_data_201708D1t01.mat'

data = sio.loadmat(fname)
F = data['F_1k']
M = data['M_1k']
Y = data['spbool_1k'].astype('float32')
X = np.hstack((M,F))

clf = Lasso()
clf.fit(X,Y)