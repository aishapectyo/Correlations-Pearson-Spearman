import sys
import os
import numpy as np
from pylab import *
import scipy
import scipy.spatial
import matplotlib.pyplot as plt
import math 
import scipy.stats as st

#Code to compute correlation between variables using the Pearson method.
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

def pearson_corr(x, y):
	assert len(x) == len(y)
	n = len(x)
	assert n > 0
	avgx = average(x)
	avgy =  average(y)
	diffmult= 0
	xdiff2 = 0
	ydiff2 = 0
	for i in range(n):
		xdiff = x[i] - avgx
		ydiff =  y[i] - avgy
		diffmult += xdiff * ydiff
		xdiff2 += xdiff * xdiff
		ydiff2 += ydiff * ydiff
	return diffmult/sqrt(xdiff2*ydiff2)

#Run Pearson for data set.
#Possible correlation of beryllium with whole sample
whole = pearson_corr(logn_be,teff)
print("from function: ",whole)

#For comparison... using Python's implemente function 
whole2  = st.pearsonr(logn_be,teff)
print("from python: ",whole2)


#Possible correlation of beryllium with planet-hosting stars
wplanet = pearson_corr(pl2,pl1)
wplanetpy = st.pearsonr(pl2,pl1)
print("from function: ",wplanet)
print("from python: ",wplanetpy)

#Possible correlation of beryllium with not planet-hosting stars
nowplanet = pearson_corr(pl4,pl3)
nowplanetpy = st.pearsonr(pl4,pl3)
print("from function: ",nowplanet)
print("from python: ",nowplanetpy)


#t-stat. Picked an alpha level of 0.05/95%
t = (whole*sqrt(31))/sqrt(1-whole**2)
tpla = (wplanet*sqrt(18))/sqrt(1-wplanet**2)
tnop =  (nowplanet*sqrt(11))/sqrt(1-nowplanet**2) 
print('t-stats: ', t, tpla, tnop)

#Plot data.
plt.figure(1)
plt.plot(teff,logn_be,marker='.',linestyle=' ',color='m')
ax=gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.xlabel('T (effective temperature, K)')
plt.ylabel('Be (abundance compared to solar abundance)')
plt.title('Whole sample, Be vs. T')

plt.figure(2)
plt.plot(pl1,pl2,marker='.',linestyle=' ',color='g')
ax=gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.xlabel('T (effective temperature, K)')
plt.ylabel('Be (abundance compared to solar abundance)')
plt.title('Planet-Hosting, Be vs. T')

plt.figure(3)
plt.plot(pl3,pl4,marker='.',linestyle=' ',color='g')
ax=gca()
ax.set_xlim(ax.get_xlim()[::-1])
plt.xlabel('T (effective temperature, K)')
plt.ylabel('Be (abundance compared to solar abundance)')
plt.title('No Planet-Hosting, Be vs. T')



plt.show()
