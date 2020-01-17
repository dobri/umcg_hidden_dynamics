#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:47:18 2019

@author: Dobri
"""

import numpy as np
import pandas as pd
# import sys, getopt


def rmse(x,y):
    return np.sqrt(np.mean((x-y)**2))

def zscore(x):
    return (x-np.mean(x))/np.std(x)

def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation. 
    Shifted data filled with NaNs 
    
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else: 
        return datax.corr(datay.shift(lag))

# X = pd.read_csv(FILEPATH,skipinitialspace=True)
# from wcc_by_rmse import wcc_by_rmse
# print(wcc_by_rmse(X['NoteL'],X['NoteR'],X['Time'],True))
def wcc_by_rmse(x,y,t,lag_range_seconds=1,step_seconds=5,PLOTTINGFLAG=True):
    er = rmse(x,y)
    x = zscore(x)
    y = zscore(y)
    fps = int(1/np.mean(np.diff(t)))
    
    # Rolling window time lagged cross correlation
    t_start = 0
    step_size = int(fps*step_seconds)
    t_end = t_start + step_size
    rss=[]
    cmax=[]
    taus=[]
    lags = range(-int(lag_range_seconds*fps),int(lag_range_seconds*fps+1))
    lags_ms_secs = np.round(np.asarray(lags,dtype=float)/fps*1e3)
    while t_end < len(t):
        rs = [crosscorr(x.iloc[t_start:t_end], y.iloc[t_start:t_end], lag, wrap=False) for lag in lags]
        #rs = [crosscorr(x[t_start:t_end], y[t_start:t_end], lag, wrap=False) for lag in lags]
        cmax.append(max(rs))
        taus.append(lags_ms_secs[np.argmax(rs)])
        rss.append(rs)
        t_start = t_start + step_size
        t_end = t_end + step_size
    rss = pd.DataFrame(rss)
    cmax_ave = np.mean(cmax)
    tau = -np.mean(taus)
    
    if PLOTTINGFLAG==True:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(10,8))
        axes[0].plot(t,x,'-')
        axes[0].plot(t,y,'-')
        axes[0].set_xlabel('Time, s')
        axes[0].set_ylabel('z-score of MIDI notes')
        axes[0].legend(loc='upper right', shadow=False, fontsize='small')    
        #plt.show()
        
        d=rss.shape
        #f,ax = plt.subplots(figsize=(10,10)) #
        sns.heatmap(np.transpose(rss),cmap='RdBu_r',ax=axes[1])
        axes[1].set(title='Rolling windowed time-lagged cross-correlation', ylim=[0,d[1]], ylabel='Lag, ms',xlabel='Epochs')
        axes[1].set_yticks(np.linspace(1,d[1],11))
        axes[1].set_yticklabels(np.linspace(min(lags_ms_secs),max(lags_ms_secs),11))
        axes[1].set_xticks(range(0,d[0]+1,3))
        axes[1].set_xticklabels(range(0,d[0]+1,3))
        plt.tight_layout()
        plt.show()
    
    return (cmax_ave/er,cmax_ave,tau,er)
    
    #filename = 'trial_log-191205-153813_task05'
    #x = np.loadtxt(filename, delimiter=',', skiprows=1)
    #print (x)
    #
    #col_names = ('Time','Acc1X','Acc1Y','Acc1Z','AbsAcc','Sensor','CPGX','CPGY','CPGZ','Force','NoteL','NoteR','Acc2X','Acc2Y','Acc2Z','DiscrUpdate')
    #
    #
    #for col in (11,10):
    #    plt.plot(X[:,0],X[:,col],'-',label=col_names[col])
    #
    #plt.xlabel('Time, s')
    #plt.ylabel('MIDI Note')
    #plt.legend(loc='upper right', shadow=False, fontsize='small')    
    #plt.show()
    
    
#    argv = sys.argv[1:]
#    try:
#        opts, args = getopt.getopt(argv,"hf:p:",["help","file=","plotting="])
#    except getopt.GetoptError:
#        print '\n'
#        print 'wcc_by_rmse.py -h'
#        print '\n'
#        sys.exit(2)
#
#    file_path = []
#    for opt, arg in opts:
#        if opt in ("-h","--help"):
#            print'wcc_by_rmse.py --file=fullpathstring --plotting=<0,1>'
#            print '\n'
#            sys.exit()
#        if opt in ("-f", "--file"):
#            FILEPATH = arg
#        if opt in ("-e", "--epsilon"):
#            PLOTTINGFLAG = bool(float(arg))
    
