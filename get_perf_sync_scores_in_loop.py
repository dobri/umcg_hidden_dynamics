#!/usr/bin/env python3
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
import matplotlib.pyplot as plt

# To verify which is the leading from neg tau.
#import numpy as np
#import pandas as pd
#from wcc_by_rmse import wcc_by_rmse
#
#t=np.linspace(0,20,100*20)
#x=np.sin(t*2*np.pi)
#x=pd.DataFrame(x)
#y=np.sin(t*2*np.pi-np.pi/4)
#y=pd.DataFrame(y)
#score,cmax,tau,err = wcc_by_rmse(x[0],y[0],t,lag_range_seconds=.5)
#tau

if __name__ == '__main__':
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hl:s:p:c1:c2:i:",["help","transient=","lag=","step=","plotting=","col1name=","col2name=","image_save="])
    except getopt.GetoptError:
        print('\n')
        sys.exit(2)

    lag_range_seconds = .5
    step_seconds = 5
    transient_secs = 10
    plotting = True
    image_save_flag = False
    col1name = 'NoteL'
    col2name = 'NoteR'

    for opt, arg in opts:
        if opt in ("-h","--help"):
            print('get_perf_sync_scores_in_loop.py --transient=<seconds> --lag=<seconds> --step=<seconds> --plotting=<0 or 1> --col1name=<var1 header label> --col2name=<var1 header label>')
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
        if opt in ("-i", "--image_save"):
            image_save_flag = bool(float(arg))

    skip_raw=False
    flist = sorted(glob.glob('trial*.txt'))
    slist = glob.glob('scores')
    if len(slist)>0:
        contents=pd.read_csv(slist[0],header=[0],na_values='     nan')
        if len(contents)==len(flist):
            print('\nExisting \'scores\' data file is same length as the number of raw data files.')
            print('Are you sure you are not re-analyzing the same data?\n')
            skip_raw=True

    if skip_raw==False:
        if len(flist)==0:
            print('No raw data files found.')
        else:
            if len(slist)==0:
                f = open('scores','a+')
                f.write("%s," % 'raw data file')
                f.write("%s," % 'score')
                f.write("%s," % 'cmax')
                f.write("%s," % 'tau')
                f.write("%s\n" % 'rmse')
                
                
            for log_file_name in flist:
                X = pd.read_csv(log_file_name,skipinitialspace=True)
                fps=1/np.mean(np.diff(X['Time']))
                
                # That's a hack, but anyway. Be careful when using later for different studies.
                if log_file_name.find('vis1')>0:
                    c1 = 'XPixStim'
                    c2 = 'XPixPart'
                else:
                    c1 = col1name
                    c2 = col2name
                # print c1,c2
                
                x = X[c1][int(fps*transient_secs):(len(X)-1)]
                y = X[c2][int(fps*transient_secs):(len(X)-1)]
                if (np.var(x)>.0001) & (np.var(y)>.0001):
                    if image_save_flag:
                        fnm = log_file_name+'_raw_wcc.png'
                    else:
                        fnm = False
                    score,cmax,tau,err = wcc_by_rmse(x,y,
                                                 X['Time'][int(fps*transient_secs):(len(X)-1)],
                                                 lag_range_seconds,
                                                 step_seconds,
                                                 plotting,fnm)
                else:
                    print('Bad raw data or an empty trial. What is wrong this time?')
                    score,cmax,tau,err = np.nan,np.nan,np.nan,np.nan
                    #print score,cmax,tau,err
                    #                    fig  = plt.figure(figsize=(20,15))
                    #                    axes = fig.add_axes([.15,.15,.75,.75])
                    #                    axes.plot(x,'-o',linewidth=5,markersize=15,label='Stim')
                    #                    axes.plot(y,'-o',linewidth=2,markersize=15,label='Participant')
                    #                    #axes.plot(contents['rmse'],'-o',linewidth=2,markersize=15,label='RMSE')
                    #                    axes.set_xlabel('Time')
                    #                    axes.set_title('Raw data from a bad trial. What is wrong here?')
                    #                    axes.legend(loc='upper left', shadow=False, fontsize='medium')
                    #                    plt.show()
                    
                f = open('scores','a+')
                f.write("%60s," % log_file_name)
                f.write("%8.3f," % score)
                f.write("%8.3f," % cmax)
                f.write("%8.3f," % (tau/1e3))
                f.write("%8.3f\n" % err)
                f.close()
                print("%8.3f" % (score*1e2), "%8.3f" % cmax, "%8.3f" % tau, "%8.3f" % err)
    
    slist = glob.glob('scores')
    if len(slist)==0:
        print('No \'scores\' file found.')
    else:
        print('\'scores\' data file found.')
        contents=pd.read_csv(slist[0],header=[0],na_values='     nan')
        if len(contents)>0:
            print(contents)
            # plotting = True
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
                #axes.plot(contents['rmse'],'-o',linewidth=2,markersize=15,label='RMSE')
                axes.set_xlabel('Trial')
                axes.set_xticks(range(0,len(contents['score']),1))
                axes.set_xticklabels(range(1,len(contents['score'])+1,1))
                axes.set_ylabel('Score')
                axes.set_title('Performance (Synchronization and Pitch-matching)')
                axes.legend(loc='upper left', shadow=False, fontsize='medium')
                if image_save_flag:
                    fnm = 'scores_'+flist[0][49:56]+flist[0][24:49]+flist[0][10:23]+'.png'
                    plt.savefig(fnm,dpi=100,bbox_inches='tight')
                else:
                    plt.show()
                
                plt.close()
        else:
            print('\'scores\' empty.')
