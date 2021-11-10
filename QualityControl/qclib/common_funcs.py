#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:59:05 2019

@author: u0101486
"""

import numpy as np
import shelve
# Extra helper functions

def vcorrcoef(X,y):
	''' Vectorized correlation coefficient
     useful to calculate correlation across the whole brain
    
     Parameters
	   --------------------------------

	   [REQUIRED]
	   X: 2-D matrix [voxels x time]
		   Reshaped 4d functional data (typically at least).
	   y: 1-D vector
		   Signal to be correlated to the functional data.
    '''
    
	Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
	ym = np.mean(y)
	r_num = np.sum((X-Xm)*(y-ym),axis=1)
	r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
	r = r_num/r_den
	return r



def save_workspace(filename, workspaceGlobals):
    
    my_shelf = shelve.open(filename,'n') # 'n' for new

    for key in workspaceGlobals.keys():
        try:
            print(key)
            my_shelf[key] = workspaceGlobals[key]
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            print('ERROR shelving: {0}'.format(key))
            
        except KeyError:
            print('Key Error [{0}] '.format(key))
            
        except:
            print('Pickle Error [{0}]'.format(key))
    my_shelf.close()


def load_workspace(filename):
    my_shelf = shelve.open(filename)
    
    workspace = {}
    
    for key in my_shelf:

        workspace[key]=my_shelf[key]
    my_shelf.close()
    
    #To pass the variables to the actual workspace, just execute the code below in the file calling load_workspace
    # workspace = load_workspace(...)
    # for key in workspace.keys():
    #   globals()[key] = workspace[key]
    # del workspace
    return workspace

