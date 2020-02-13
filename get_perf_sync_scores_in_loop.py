#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:19:57 2020
@author: Dobri@McMasterU
"""
import glob
from wcc_by_rmse import wcc_by_rmse
import pandas as pd
import numpy as np
import sys, getopt

if __name__ == '__main__':
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hl:s:p:c1:c2:",["help","transient=","lag=","step=","plotting=","col1name=","col2name="])
    except getopt.GetoptError:
        print '\n'
        sys.exit(2)

    lag_range_seconds=1
    step_seconds=5
    plotting=True
    col1name='NoteL'
    col2name='NoteR'
    transient_secs = 0

    for opt, arg in opts:
        if opt in ("-h","--help"):
            print'get_perf_sync_scores_in_loop.py --transient=<seconds> --lag=<seconds> --step=<seconds> --plotting=<0 or 1> --col1name=<var1 header label> --col2name=<var1 header label>'
            sys.exit()
        if opt in ("-t", "--transient"):
            transient_secs = float(arg)
        if opt in ("-l", "--lag"):
            lag_range_seconds = float(arg)
        if opt in ("-s", "--step"):
            step_seconds = float(arg)
        if opt in ("-p", "--plotting"):
            plotting = bool(float(arg))
        if opt in ("-c1", "--col1name"):
            col1name = arg
        if opt in ("-c2", "--col2name"):
            col2name = arg

    skip_raw=False
    flist = sorted(glob.glob('trial*'))
    slist = glob.glob('scores')
    if len(slist)>0:
        contents=pd.read_csv(slist[0],header=[0],na_values='     nan')
        if len(contents)==len(flist):
            print '\nExisting \'scores\' data file is same length as the number of raw data files.'
            print 'Are you sure you are not re-analyzing the same data?\n'
            skip_raw=True

    if skip_raw==False:
        if len(flist)==0:
            print 'No raw data files found.'
        else:
            if len(slist)==0:
                f = open('scores','a+')
                f.write("%s," % 'raw file name')
                f.write("%s," % 'score')
                f.write("%s," % 'cmax')
                f.write("%s," % 'tau')
                f.write("%s\n" % 'rmse')
                
                
            for log_file_name in flist:
                X = pd.read_csv(log_file_name,skipinitialspace=True)
                fps=1/np.mean(np.diff(X['Time']))
                
                x = X[col1name][int(fps*transient_secs):(len(X)-1)]
                y = X[col2name][int(fps*transient_secs):(len(X)-1)]
                if (np.var(x)>.0001) & (np.var(y)>.0001):
                    score,cmax,tau,err = wcc_by_rmse(x,y,
                                                 X['Time'][int(fps*transient_secs):(len(X)-1)],
                                                 lag_range_seconds,
                                                 step_seconds,
                                                 plotting)
                else:
                    score,cmax,tau,err = np.nan,np.nan,np.nan,np.nan
                    #print score,cmax,tau,err
                    
                f = open('scores','a+')
                f.write("%60s," % log_file_name)
                f.write("%8.3f," % score)
                f.write("%8.3f," % cmax)
                f.write("%8.3f," % (tau/1e3))
                f.write("%8.3f\n" % err)
                f.close()
                print "%8.3f" % (score*1e2), "%8.3f" % cmax, "%8.3f" % tau, "%8.3f" % err
    
    slist = glob.glob('scores')
    if len(slist)==0:
        print 'No \'scores\' file found.'
    else:
        print '\'scores\' data file found.'
        contents=pd.read_csv(slist[0],header=[0],na_values='     nan')
        if len(contents)>0:
            print contents
            if plotting:
                import matplotlib.pyplot as plt
            
                font = {'family':'DejaVu Sans','weight':'bold','size':20}
                plt.rc('font', **font)
                plt.rc('axes', labelsize=30)
            
                fig  = plt.figure(figsize=(20,15))
                axes = fig.add_axes([.15,.15,.75,.75])
                axes.plot(contents['score'],'-o',linewidth=5,markersize=15,label='C_{max}/RMSE')
                axes.plot(contents['cmax'],'-o',linewidth=2,markersize=15,label='C_{max}')
                axes.plot(contents['tau'],'-o',linewidth=2,markersize=15,label='tau')
                axes.plot(contents['rmse'],'-o',linewidth=2,markersize=15,label='RMSE')
                axes.set_xlabel('Trial')
                axes.set_xticks(range(0,len(contents['score']),1))
                axes.set_xticklabels(range(1,len(contents['score'])+1,1))
                axes.set_ylabel('Score')
                axes.set_title('Performance (Synchronization and Pitch-matching)')
                axes.legend(loc='upper left', shadow=False, fontsize='medium')
                plt.show()
        else:
            print '\'scores\' empty.'