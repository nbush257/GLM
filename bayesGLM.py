from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib.pyplot as plt 
import numpy as np
from scipy.io.matlab import loadmat
from sklearn.linear_model import BayesianRidge

def loadIn():
	Tk().withdraw()
	fName = askopenfilename(initialdir = 'D:')

	fd = loadmat(fName)
	return fd
def noProx(X,prox):
	X = X.take(np.where(prox==0)[0],0)
	return X

def apply(X,y):
	mdl = BayesianRidge()
	if y.ndim==2 and y.shape[1] ==1:
		y = y.ravel
	if X.shape[1]>X.shape[0]:
		X = X.T

	mdl.fit(X,y)
	yhat = mdl.predict(X)
	return mdl,yhat