import sys
import os
import numpy as np
from pylab import *
import scipy
import scipy.spatial
import matplotlib.pyplot as plt
import math
import scipy.stats as st

#Code to compute correlation between variables using the Spearman method.
#Aisha Mahmoud-Perez


#Read data file. 
data = genfromtxt('BeLi.dat')
star = data[:,0]
tipo = data[:,1]
teff = data[:,2]
logn_be = data[:,3]
sig_be = data[:,4]
logn_li = data[:,5]

pl1=[]
pl2=[]
pl3=[]
pl4=[]
m = len(tipo)
for i in range(m):
	if tipo[i] == 1:
		pl1.append(teff[i])
		pl2.append(logn_be[i])
	if tipo[i] == 2:
		pl3.append(teff[i])
		pl4.append(logn_be[i])
#Create functions.
def average(x):
	assert len(x) > 0
	return float(sum(x)) / len(x)

def spearman_corr(x,y):
	assert len(x) == len(y)
	n = len(x)
	assert n > 0 
	xrank = st.rankdata(x)
	yrank = st.rankdata(y)
	avgx = average(xrank)
	avgy =  average(yrank)
	diffmult= 0
	xdiff2 = 0
	ydiff2 = 0
	for i in range(n):
		 xdiff = xrank[i] - avgx
		 ydiff =  yrank[i] - avgy
		 diffmult += xdiff * ydiff
		 xdiff2 += xdiff * xdiff
		 ydiff2 += ydiff * ydiff
	return diffmult/sqrt(xdiff2*ydiff2)


#Run Spearman for data set.
#Possible correlation of beryllium with whole sample
whole = spearman_corr(logn_be,teff)
print("from function: ",whole)
#For comparison... using Python's implemente function 
whole2  = st.spearmanr(logn_be,teff)
print("from python: ",whole2)

#Possible correlation of beryllium with planet-hosting stars
wplanet = spearman_corr(pl2,pl1)
wplanetpy = st.spearmanr(pl2,pl1)
print("from function: ",wplanet)
print("from python: ",wplanetpy)

#Possible correlation of beryllium with not planet-hosting stars
nowplanet = spearman_corr(pl4,pl3)
nowplanetpy = st.spearmanr(pl4,pl3)
print("from function: ",nowplanet)
print("from python: ",nowplanetpy)

#t
t = whole*sqrt(31/(1-whole**2))
tpla = (wplanet*sqrt(18))/sqrt(1-wplanet**2)
tnop =  (nowplanet*sqrt(11))/sqrt(1-nowplanet**2)
print('t-stats: ', t, tpla, tnop)
