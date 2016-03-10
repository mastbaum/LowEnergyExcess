####################################################################################
### This module takes as input a dictionary of dataframes, perhaps from a low
### energy excess plotting script. The keys of dataframes are the names of the 
### samples, and the keys needed by this script (as of now) are 'cosmic' and 'nue'.
### This script chooses variables and computes a likelihood that an event is
### nue or cosmic. It appends a column to each of the dataframes in the dictionary
### titled "hypothesis", which has values 0 to 1 (cosmic to nue)
####################################################################################

from ROOT import TFile, TTree
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from collections import OrderedDict
from scipy.special import expit #Vectorized sigmoid function
from scipy import optimize

#Hypothesis function and cost function for logistic regression
def h(mytheta,myX): 
	return expit(np.dot(myX,mytheta))

#Cost function, default lambda (regularization) 0
def computeCost(mytheta,myX,myy,mylambda = 0.): 
	"""
	theta_start is an n- dimensional vector of initial theta guess
	X is matrix with n- columns and m- rows
	y is a matrix with m- rows and 1 column
	Note this includes regularization, if you set mylambda to nonzero
	"""
	term1 = np.dot(-np.array(myy).T,np.log(h(mytheta,myX)))
	term2 = np.dot((1-np.array(myy)).T,np.log(1-h(mytheta,myX)))
	regterm = (mylambda/2) * np.sum(np.dot(mytheta[1:].T,mytheta[1:])) #Skip theta0
	mym = myy.size
	return float( (1./mym) * ( np.sum(term1 - term2) + regterm ) )

#Optimization functions:
#An alternative to OCTAVE's 'fminunc' we'll use some scipy.optimize function, "fmin"
#Note "fmin" does not need to be told explicitly the derivative terms
#It only needs the cost function, and it minimizes with the "downhill simplex algorithm."
#http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.fmin.html
def optimizeTheta(mytheta,myX,myy,mylambda=0.):
	result = optimize.fmin(computeCost, x0=mytheta, args=(myX, myy, mylambda), maxiter=400, full_output=True)
	return result[0], result[1]

def computeHypothesis(dfs,input_theta=None):
	""" 
	Input is a(n) (ordered) dictionary of dataframes as described in the header
	This function modifies the dataframes directly, adding a column called "hypothesis"
	"""
	#First add another column to the dataframes, computed from existing columns
	for key in dfs.keys():
		dfs[key]['ptoverp'] = dfs[key]['_nu_pt']/dfs[key]['_nu_p']
                denom = len(dfs[key]['ptoverp'])
                dfs[key] = dfs[key].query('ptoverp > 0 and ptoverp <= 1')
                num = len(dfs[key]['ptoverp'])
                if denom-num:
                        print "NOTE! In sample %s, you just dropped %d out of %d events."%(key,denom-num,denom)

	#Build the feature matrix, X, from cosmic samples (with y = 0) and nue samples (with y = 1)
	cosX = dfs['cosmic'][['ptoverp', '_y_vtx']].as_matrix()
	cosX = np.insert(cosX,0,1,axis=1)
	cosy = np.zeros((cosX.shape[0],1))
	nueX = dfs['nue'][['ptoverp', '_y_vtx']].as_matrix()
	nueX = np.insert(nueX,0,1,axis=1)
	nuey = np.ones((nueX.shape[0],1))
	X = np.vstack((cosX,nueX))
	y = np.vstack((cosy,nuey))
	m = y.size # number of training samples

	#Feature normalizing the columns (subtract mean, divide by standard deviation)
	#Store the mean and std for later use
	#Note don't modify the original X matrix, use a copy
	stored_feature_means, stored_feature_stds = [], []
	Xnorm = X.copy()

	for icol in xrange(Xnorm.shape[1]):
		#Skip the first column
		if not icol:
			stored_feature_means.append(0.)
			stored_feature_stds.append(1.)
		else: 
			stored_feature_means.append(np.mean(Xnorm[:,icol]))
			stored_feature_stds.append(np.std(Xnorm[:,icol]))
			#Faster to not recompute the mean and std again, just used stored values
			Xnorm[:,icol] = (Xnorm[:,icol] - stored_feature_means[-1])/stored_feature_stds[-1]

	#Plug in a random initial theta and solve for the coefficients:
	initial_theta = np.zeros((Xnorm.shape[1],1))
        if input_theta:
                theta = np.array(input_theta)
                print "You input theta! Theta type is ",type(theta)
        else:
                theta, mincost = optimizeTheta(initial_theta,Xnorm,y)


	#Now that we have an optimized theta fit parameter vector, add hypothesis as a column to each df
	#Hacky insert 0 in stored feature mean and 1 in stored feature stds
	stored_feature_means[0]=0
	stored_feature_stds[0]=1
	for key in dfs.keys():
		testpoints = np.array([dfs[key]['ptoverp'],dfs[key]['_y_vtx']]).T
		testpoints = np.insert(testpoints,0,1,axis=1)
		testpointsscaled = (testpoints-stored_feature_means)/stored_feature_stds
		dfs[key]['hypothesis'] = h(np.array(theta),np.array([testpointsscaled])).T

	#Return the theta value for cross checks that minimizer has converged on the same value
	return theta





